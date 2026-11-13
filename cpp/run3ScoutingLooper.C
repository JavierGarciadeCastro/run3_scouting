#include <algorithm>
#include <exception>
#include <cstdlib>
#include <filesystem>
#include <numeric>
#include <string>
#include <tuple>
namespace fs = std::filesystem;

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
using namespace fwlite;

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"

#include "tools/dorky.h"
#include "tools/goodrun.h"
#include "tools/tqdm.h"


bool doPartialUnblinding = true;
float partialUnblindingPercentage = 0.1; // 10% of each era

// SV selection
bool relaxedSVSel = false;
float sfSVsel = (relaxedSVSel) ? 5.0 : 1.0;
//float maxXerr=0.05*sfSVsel, maxYerr=0.05*sfSVsel, maxZerr=0.10*sfSVsel, maxChi2=3.0*sfSVsel;
float maxDXYerr=0.05*sfSVsel, maxD3Derr=0.10*sfSVsel; // for identification of overlapping SVs

// Trigger selection
bool applyL1 = false;   // set false to skip the L1 requirement
bool requireAtLeastOneSV = true;

// BToPhi specific
bool filterByB = true;

template<class T>
const T getObject(const Event& ev, const char* prodLabel, const char* outLabel = "") {
  Handle<T> obj;
  obj.getByLabel(ev,prodLabel,outLabel);
  return *obj.product();
}

template <typename T>
const T getObjectSafe(const Event& ev, const char* prodLabel, const char* outLabel = "") {
  Handle<T> obj;
  obj.getByLabel(ev, prodLabel, outLabel);
  if (!obj.isValid()) return T{};
  return *obj.product();
}

// Branch names in EDM files are "<type>_<label>_<instance>_<process>.", so the process name
// (PAT, HLT, ...) is not known in advance: match on the leading "<type>_<label>" part.
bool branchWithPrefix(TTree* tree, const char* prefix) {
  TObjArray* branches = tree->GetListOfBranches();
  for (int i=0; i<branches->GetEntries(); ++i) {
    if (TString(branches->At(i)->GetName()).BeginsWith(prefix))
      return true;
  }
  return false;
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
  return p;
}

template <typename T>
void apply_permutation_in_place(std::vector<T>& vec, const std::vector<std::size_t>& p) {
  std::vector<bool> done(vec.size());
  for (std::size_t i=0; i<vec.size(); ++i) {
    if (done[i])
      continue;
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j) {
      //std::swap(vec[prev_j], vec[j]); // fails for vector<bool> proxy refs in gcc12
      auto tmp = static_cast<T>(vec[prev_j]); vec[prev_j] = vec[j]; vec[j] = tmp;
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
}


float getMiniIsoDR(const TLorentzVector muVec, const float minDR=0.05, const float maxDR=0.2, const float kt=10.0) {
  return std::min(maxDR, std::max(minDR, kt/(float)muVec.Pt()));
}

std::tuple<float,float,float,float,float> getPFIsolation(const TLorentzVector muVec, const std::vector<Run3ScoutingParticle>& pfs, const float maxDR=0.3, const float vetoDR=0.01, const float maxdz=0.1,const bool doMini=false, const bool doMindrMuPF=true) {
  float maxdr, mindrMuPF=1e6;
  float chIso=0, nhIso=0, phIso=0, puIso=0;
  if (doMini)
    maxdr = getMiniIsoDR(muVec);
  for (auto pf : pfs) {
    TLorentzVector pfvec;
    pfvec.SetPtEtaPhiM(pf.pt(),pf.eta(),pf.phi(),0.0);
    float dr = muVec.DeltaR(pfvec);
    if (dr<mindrMuPF)
      mindrMuPF = dr;
    if (dr<vetoDR || dr>maxDR)
      continue;
    bool fromPV = !(abs(pf.vertex())>0 || pf.dz()>maxdz);
    bool isCh=false, isNh=false, isPh=false, isPU=false;
    if (fabs(pf.pdgId())==211) {
      if (fromPV) {
        isCh = true;
        chIso = chIso + pf.pt();
      }
      else {
        isPU = true;
        puIso = puIso + pf.pt();
      }
    }
    if (fabs(pf.pdgId())==130) {
      isNh  = true;
      nhIso = nhIso + pf.pt();
    }
    if (fabs(pf.pdgId())==22) {
      isPh  = true;
      phIso = phIso + pf.pt();
    }
  }
  if (doMindrMuPF)
    return {chIso,nhIso,phIso,puIso,mindrMuPF};
  else
    return {chIso,nhIso,phIso,puIso,1e6};
}

std::tuple<float,float,float> get_track_reference_point(const Run3ScoutingMuon& mu, const float dxyCorr) {
  const auto phi = mu.phi();
  const auto dsz = mu.trk_dsz();
  const auto dz = mu.trk_dz();
  const auto lmb = mu.trk_lambda();

  const auto sinphi = TMath::Sin(phi);
  const auto cosphi = TMath::Cos(phi);
  const auto sinlmb = TMath::Sin(lmb);
  const auto tanlmb = sinlmb/TMath::Cos(lmb);

  float refz = 1.0*dz;
  float refx = -sinphi*dxyCorr - (cosphi/sinlmb)*dsz + (cosphi/tanlmb)*refz;
  float refy =  cosphi*dxyCorr - (sinphi/sinlmb)*dsz + (sinphi/tanlmb)*refz;

  return {refx,refy,refz};
}

float recalculate_phi_at_DV(
  float refx, float refy, float refz, // initial position in cm
  float px, float py, float pz, // initial momentum in GeV
  int charge, // -1 or 1, matches Muon_charge
  float dvx, float dvy // DV coordinates to propagate to
  ) {
  // Performs a helix propagation from 
  // track reference point (initial position)
  // to cylinder containing DV x,y, and then reports
  // the updated phi coordinate the the muon at that point
  float B = 3.8; // b field in Tesla
  float mass = 0.10566; // muon mass in GeV
  float c = 0.29979; // speed of light in m/ns

  float P = pow(px*px + py*py + pz*pz, 0.5);
  float Pxy = pow(px*px + py*py, 0.5);
  float E = pow(P*P + mass*mass, 0.5);

  // relativistic velocities
  float vx = px/E * c;
  float vy = py/E * c;
  float vz = pz/E * c;
  float vxy = pow(vx*vx + vy*vy, 0.5);

  // larmor radius/angular frequency
  float R = 1e3*Pxy/(3.*charge*B);
  float w = vxy/R;

  float t = 0.;
  float curr_x = refx;
  float curr_y = refy;
  float current_rho = -1;
  float target_rho = pow(dvx*dvx + dvy*dvy, 0.5);
  float dt = 0.1;
  int nsteps = 0;
  while(current_rho < target_rho) {
      curr_x = refx + (vy/w)*( 1-TMath::Cos(w*t)) + (vx/w)*TMath::Sin(w*t);
      curr_y = refy + (vx/w)*(-1+TMath::Cos(w*t)) + (vy/w)*TMath::Sin(w*t);
      current_rho = pow(curr_x*curr_x + curr_y*curr_y, 0.5);
      t += dt;
      nsteps++;
      if (nsteps > 10000) {
          //std::cout << "Warning, >10000 steps in propagate_to_cylinder" << std::endl;
          break;
      }
  }
  float curr_vx = vy*TMath::Sin(w*t) + vx*TMath::Cos(w*t);
  float curr_vy = vy*TMath::Cos(w*t) - vx*TMath::Sin(w*t);
  float newphi = atan2(curr_vy, curr_vx);
  return newphi;
}

float getCorrectedPhi(const Run3ScoutingMuon& mu, const TLorentzVector muVec, const float dxyCorr, const float bestSVPosition_x, const float bestSVPosition_y) {
  const auto refs = get_track_reference_point(mu, dxyCorr);
  const auto refx = std::get<0>(refs);
  const auto refy = std::get<1>(refs);
  const auto refz = std::get<2>(refs);
  float phiCorr = recalculate_phi_at_DV(refx, refy, refz, muVec.Px(), muVec.Py(), muVec.Pz(), mu.charge(), bestSVPosition_x, bestSVPosition_y);
  return phiCorr;
}


struct GenPart {
  std::vector<float> pt, eta, phi, m;
  std::vector<float> vx, vy, vz, lxy;
  std::vector<int> status, index, pdgId, motherIndex, motherPdgId;
  std::vector<float> ct;

  void clear() {
    pt.clear(); eta.clear(); phi.clear(); m.clear();
    vx.clear(); vy.clear(); vz.clear(); lxy.clear();
    status.clear(); index.clear(); pdgId.clear(); motherIndex.clear(); motherPdgId.clear();
    ct.clear();
  }
};

struct SV {
  std::vector<unsigned int> index, ndof;
  std::vector<float> x, y, z;
  std::vector<float> xe, ye, ze;
  std::vector<float> chi2, prob, chi2Ndof;
  std::vector<float> lxy, l3d;
  std::vector<float> mindx, mindy, mindz, mindxy, mind3d;
  std::vector<float> maxdx, maxdy, maxdz, maxdxy, maxd3d;
  std::vector<bool> selected;
  std::vector<bool> onModule, onModuleWithinUnc;
  std::vector<float> closestDet_x, closestDet_y, closestDet_z;
  std::vector<float> minDistanceFromDet, minDistanceFromDet_x, minDistanceFromDet_y, minDistanceFromDet_z;

  void clear() {
    index.clear(); ndof.clear();
    x.clear(); y.clear(); z.clear();
    xe.clear(); ye.clear(); ze.clear();
    chi2.clear(); prob.clear(); chi2Ndof.clear();
    lxy.clear(); l3d.clear();
    mindx.clear(); mindy.clear(); mindz.clear(); mindxy.clear(); mind3d.clear();
    maxdx.clear(); maxdy.clear(); maxdz.clear(); maxdxy.clear(); maxd3d.clear();
    selected.clear();
    onModule.clear(); onModuleWithinUnc.clear();
    closestDet_x.clear(); closestDet_y.clear(); closestDet_z.clear();
    minDistanceFromDet.clear(), minDistanceFromDet_x.clear(), minDistanceFromDet_y.clear(), minDistanceFromDet_z.clear();
  }

  void sort() {
    auto comp = sort_permutation(prob, [](float const& a, float const& b){ return a > b; }); // Define permutations based on the prob vector
    // Apply the permutation to all vectors
    apply_permutation_in_place(index, comp);
    apply_permutation_in_place(ndof, comp);
    apply_permutation_in_place(x, comp);
    apply_permutation_in_place(y, comp);
    apply_permutation_in_place(z, comp);
    apply_permutation_in_place(xe, comp);
    apply_permutation_in_place(ye, comp);
    apply_permutation_in_place(ze, comp);
    apply_permutation_in_place(chi2, comp);
    apply_permutation_in_place(prob, comp);
    apply_permutation_in_place(chi2Ndof, comp);
    apply_permutation_in_place(lxy, comp);
    apply_permutation_in_place(l3d, comp);
    apply_permutation_in_place(mindx, comp);
    apply_permutation_in_place(mindy, comp);
    apply_permutation_in_place(mindz, comp);
    apply_permutation_in_place(mindxy, comp);
    apply_permutation_in_place(mind3d, comp);
    apply_permutation_in_place(maxdx, comp);
    apply_permutation_in_place(maxdy, comp);
    apply_permutation_in_place(maxdz, comp);
    apply_permutation_in_place(maxdxy, comp);
    apply_permutation_in_place(maxd3d, comp);
    apply_permutation_in_place(selected, comp);
    apply_permutation_in_place(onModule, comp);
    apply_permutation_in_place(onModuleWithinUnc, comp);
    apply_permutation_in_place(closestDet_x, comp);
    apply_permutation_in_place(closestDet_y, comp);
    apply_permutation_in_place(closestDet_z, comp);
  }
};

struct SVOverlap {
  std::vector<std::vector<unsigned int>> vtxIdxs;
  std::vector<float> x, y, z;
  std::vector<float> lxy, l3d;

  void clear() {
    vtxIdxs.clear();
    x.clear(); y.clear(); z.clear();
    lxy.clear(); l3d.clear();
  }
};

// Reconstructed secondary vertex from skimmer VertexMaker (for Vtx muon collection, 2024 only)
struct SVVtx {
  std::vector<unsigned int> ndof;
  std::vector<unsigned int> origIdx; // index in original verticesVtx collection
  std::vector<float> x, y, z;
  std::vector<float> xe, ye, ze;
  std::vector<float> chi2, prob, chi2Ndof;
  std::vector<float> lxy, l3d;
  std::vector<bool> selected;
  std::vector<bool> onModule, onModuleWithinUnc;
  std::vector<float> closestDet_x, closestDet_y, closestDet_z;
  std::vector<float> minDistanceFromDet, minDistanceFromDet_x, minDistanceFromDet_y, minDistanceFromDet_z;
  std::vector<float> mindx, mindy, mindz, mindxy, mind3d;
  std::vector<float> maxdx, maxdy, maxdz, maxdxy, maxd3d;

  void clear() {
    ndof.clear(); origIdx.clear();
    x.clear(); y.clear(); z.clear();
    xe.clear(); ye.clear(); ze.clear();
    chi2.clear(); prob.clear(); chi2Ndof.clear();
    lxy.clear(); l3d.clear();
    selected.clear();
    onModule.clear(); onModuleWithinUnc.clear();
    closestDet_x.clear(); closestDet_y.clear(); closestDet_z.clear();
    minDistanceFromDet.clear(); minDistanceFromDet_x.clear(); minDistanceFromDet_y.clear(); minDistanceFromDet_z.clear();
    mindx.clear(); mindy.clear(); mindz.clear(); mindxy.clear(); mind3d.clear();
    maxdx.clear(); maxdy.clear(); maxdz.clear(); maxdxy.clear(); maxd3d.clear();
  }

  void sort() {
    auto comp = sort_permutation(prob, [](float const& a, float const& b){ return a > b; });
    apply_permutation_in_place(ndof, comp);
    apply_permutation_in_place(origIdx, comp);
    apply_permutation_in_place(x, comp);
    apply_permutation_in_place(y, comp);
    apply_permutation_in_place(z, comp);
    apply_permutation_in_place(xe, comp);
    apply_permutation_in_place(ye, comp);
    apply_permutation_in_place(ze, comp);
    apply_permutation_in_place(chi2, comp);
    apply_permutation_in_place(prob, comp);
    apply_permutation_in_place(chi2Ndof, comp);
    apply_permutation_in_place(lxy, comp);
    apply_permutation_in_place(l3d, comp);
    apply_permutation_in_place(selected, comp);
    apply_permutation_in_place(onModule, comp);
    apply_permutation_in_place(onModuleWithinUnc, comp);
    apply_permutation_in_place(closestDet_x, comp);
    apply_permutation_in_place(closestDet_y, comp);
    apply_permutation_in_place(closestDet_z, comp);
    apply_permutation_in_place(minDistanceFromDet, comp);
    apply_permutation_in_place(minDistanceFromDet_x, comp);
    apply_permutation_in_place(minDistanceFromDet_y, comp);
    apply_permutation_in_place(minDistanceFromDet_z, comp);
    apply_permutation_in_place(mindx, comp);
    apply_permutation_in_place(mindy, comp);
    apply_permutation_in_place(mindz, comp);
    apply_permutation_in_place(mindxy, comp);
    apply_permutation_in_place(mind3d, comp);
    apply_permutation_in_place(maxdx, comp);
    apply_permutation_in_place(maxdy, comp);
    apply_permutation_in_place(maxdz, comp);
    apply_permutation_in_place(maxdxy, comp);
    apply_permutation_in_place(maxd3d, comp);
  }
};

float MUON_MASS = 0.10566;

bool isGlobalMuon(const unsigned int type) {
  return type & (1<<1);
}

bool isTrackerMuon(const unsigned int type) {
  return type & (1<<2);
}

bool isStandAloneMuon(const unsigned int type) {
  return type & (1<<3);
}

struct Muon {
  std::vector<std::vector<int>> vtxIdxs;
  std::vector<unsigned int> saHits, saMatchedStats;
  std::vector<unsigned int> muHits, muChambs, muCSCDT, muMatch, muMatchedStats, muExpMatchedStats, muMatchedRPC;
  std::vector<unsigned int> pixHits, stripHits;
  std::vector<unsigned int> pixLayers, trkLayers;
  std::vector<int> bestAssocSVIdx, bestAssocSVOverlapIdx;
  std::vector<int> ch;
  std::vector<int> isGlobal, isTracker, isStandAlone;
  std::vector<float> pt, eta, phi;
  std::vector<float> chi2Ndof;
  std::vector<float> ecalIso, hcalIso, trackIso;
  std::vector<float> ecalRelIso, hcalRelIso, trackRelIso;
  std::vector<float> dxy, dxye, dz, dze;
  std::vector<float> dxysig, dzsig;
  std::vector<float> phiCorr, dxyCorr;
  std::vector<float> PFIsoChg0p3, PFIsoChg0p4, PFIsoAll0p3, PFIsoAll0p4;
  std::vector<float> PFRelIsoChg0p3, PFRelIsoChg0p4, PFRelIsoAll0p3, PFRelIsoAll0p4;
  std::vector<float> mindrPF0p3, mindrPF0p4;
  std::vector<float> mindr, maxdr;
  std::vector<float> mindrJet, mindphiJet, mindetaJet;
  std::vector<TLorentzVector> vec;
  std::vector<bool> selected;
  std::vector<int> nhitsbeforesv, ncompatible, ncompatibletotal, nexpectedhits, nexpectedhitsmultiple, nexpectedhitsmultipletotal, nexpectedhitstotal;

  void clear() {
  vtxIdxs.clear();
  saHits.clear(); saMatchedStats.clear();
  muHits.clear(); muChambs.clear(); muCSCDT.clear(); muMatch.clear(); muMatchedStats.clear(); muExpMatchedStats.clear(); muMatchedRPC.clear();
  pixHits.clear(); stripHits.clear();
  pixLayers.clear(); trkLayers.clear();
  bestAssocSVIdx.clear(); bestAssocSVOverlapIdx.clear();
  ch.clear();
  isGlobal.clear(); isTracker.clear(); isStandAlone.clear();
  pt.clear(); eta.clear(); phi.clear();
  chi2Ndof.clear();
  ecalIso.clear(); hcalIso.clear(); trackIso.clear();
  ecalRelIso.clear(); hcalRelIso.clear(); trackRelIso.clear();
  dxy.clear(); dxye.clear(); dz.clear(); dze.clear();
  dxysig.clear(); dzsig.clear();
  phiCorr.clear(); dxyCorr.clear();
  PFIsoChg0p3.clear(); PFIsoAll0p3.clear();
  PFRelIsoChg0p3.clear(); PFRelIsoAll0p3.clear();
  PFIsoChg0p4.clear(); PFIsoAll0p4.clear();
  PFRelIsoChg0p4.clear(); PFRelIsoAll0p4.clear();
  mindrPF0p3.clear(); mindrPF0p4.clear();
  mindr.clear(); maxdr.clear();
  mindrJet.clear(); mindphiJet.clear(); mindetaJet.clear();
  vec.clear();
  selected.clear();
  nhitsbeforesv.clear(); ncompatible.clear(); ncompatibletotal.clear(); nexpectedhits.clear(); nexpectedhitsmultiple.clear(); nexpectedhitsmultipletotal.clear(); nexpectedhitstotal.clear();
  }

  void sort() {
    auto comp = sort_permutation(pt, [](float const& a, float const& b){ return a > b; }); // Define permutations based on the prob vector
    // Apply the permutation to all vectors
    apply_permutation_in_place(vtxIdxs, comp);
    apply_permutation_in_place(saHits, comp);
    apply_permutation_in_place(saMatchedStats, comp);
    apply_permutation_in_place(muHits, comp);
    apply_permutation_in_place(muChambs, comp);
    apply_permutation_in_place(muCSCDT, comp);
    apply_permutation_in_place(muMatch, comp);
    apply_permutation_in_place(muMatchedStats, comp);
    apply_permutation_in_place(muExpMatchedStats, comp);
    apply_permutation_in_place(muMatchedRPC, comp);
    apply_permutation_in_place(pixHits, comp);
    apply_permutation_in_place(stripHits, comp);
    apply_permutation_in_place(pixLayers, comp);
    apply_permutation_in_place(trkLayers, comp);
    apply_permutation_in_place(bestAssocSVIdx, comp);
    apply_permutation_in_place(bestAssocSVOverlapIdx, comp);
    apply_permutation_in_place(ch, comp);
    apply_permutation_in_place(isGlobal, comp);
    apply_permutation_in_place(isTracker, comp);
    apply_permutation_in_place(isStandAlone, comp);
    apply_permutation_in_place(pt, comp);
    apply_permutation_in_place(eta, comp);
    apply_permutation_in_place(phi, comp);
    apply_permutation_in_place(chi2Ndof, comp);
    apply_permutation_in_place(ecalIso, comp);
    apply_permutation_in_place(hcalIso, comp);
    apply_permutation_in_place(trackIso, comp);
    apply_permutation_in_place(ecalRelIso, comp);
    apply_permutation_in_place(hcalRelIso, comp);
    apply_permutation_in_place(trackRelIso, comp);
    apply_permutation_in_place(dxy, comp);
    apply_permutation_in_place(dxye, comp);
    apply_permutation_in_place(dz, comp);
    apply_permutation_in_place(dze, comp);
    apply_permutation_in_place(dxysig, comp);
    apply_permutation_in_place(dzsig, comp);
    apply_permutation_in_place(phiCorr, comp);
    apply_permutation_in_place(dxyCorr, comp);
    apply_permutation_in_place(PFIsoChg0p3, comp);
    apply_permutation_in_place(PFIsoAll0p3, comp);
    apply_permutation_in_place(PFRelIsoChg0p3, comp);
    apply_permutation_in_place(PFRelIsoAll0p3, comp);
    apply_permutation_in_place(mindrPF0p3, comp);
    apply_permutation_in_place(PFIsoChg0p4, comp);
    apply_permutation_in_place(PFIsoAll0p4, comp);
    apply_permutation_in_place(PFRelIsoChg0p4, comp);
    apply_permutation_in_place(PFRelIsoAll0p4, comp);
    apply_permutation_in_place(mindrPF0p4, comp);
    apply_permutation_in_place(mindr, comp);
    apply_permutation_in_place(maxdr, comp);
    apply_permutation_in_place(mindrJet, comp);
    apply_permutation_in_place(mindphiJet, comp);
    apply_permutation_in_place(mindetaJet, comp);
    apply_permutation_in_place(vec, comp);
    apply_permutation_in_place(selected, comp);
    apply_permutation_in_place(nhitsbeforesv, comp);
    apply_permutation_in_place(ncompatible, comp);
    apply_permutation_in_place(ncompatibletotal, comp);
    apply_permutation_in_place(nexpectedhits, comp);
    apply_permutation_in_place(nexpectedhitsmultiple, comp);
    apply_permutation_in_place(nexpectedhitsmultipletotal, comp);
    apply_permutation_in_place(nexpectedhitstotal, comp);
  }
};

void run3ScoutingLooper(std::vector<TString> inputFiles, TString year, TString process, const char* outdir="temp_data", TString label="") {
  // Output folders and files
  fs::create_directory(outdir);
  fs::permissions(outdir,fs::perms::owner_all | fs::perms::group_read | fs::perms::group_exec | fs::perms::others_read | fs::perms::others_exec);
  TFile* fout = new TFile(TString(outdir)+"/output_"+process+"_"+year+label+".root", "RECREATE");
  TTree* tout = new TTree("tout","Run3ScoutingTree");
  TH1F* counts = new TH1F("counts", "", 1, 0, 1);
  TH1F* sum2Weights = new TH1F("sum2Weights", "", 1, 0, 1);
  TH1F* cutflow = new TH1F("cutflow", "", 9, 0, 9);

  // Branch variables
  unsigned int run, lumi, evtn;
  float nPU, nPUTrue, wPU; 
  bool passL1, passHLT;
  float nPV, PV_x, PV_y, PV_z;
  float GenB_pt, GenB_eta;
  int nGenBplus;
  std::vector<float> GenBplus_ct;
  std::vector<int> GenBplus_type;
  int nMuons_NoVtx;
  int nMuons_Vtx;
  int nSVs_NoVtx;
  int nSVs_Vtx;
  //float ct1 = -1.0, ct2 = -1.0;
  GenPart GenParts;
  SV SVs;
  SVOverlap SVOverlaps;
  SVOverlap SVOverlapVtxs;
  int nMuon_Assoc;
  int nMuon_vtx_Assoc;
  Muon Muons;
  Muon MuonsVtx;
  SVVtx SVsVtx;
  TFile *file0 = nullptr;
  for (auto inputFile : inputFiles) {
    file0 = TFile::Open(inputFile);
    if (file0 && !file0->IsZombie() && file0->Get("Events")) break;
    std::cout << "Cannot open " << inputFile << " for header probe, trying next..." << std::endl;
    if (file0) { file0->Close(); file0 = nullptr; }
  }
  if (!file0) { std::cout << "No readable input file in this chunk, exiting." << std::endl; fout->Close(); std::exit(1); }
  Event ev0(file0);
  ev0.toBegin();
  auto l1Names = getObject<std::vector<std::string>>(ev0, "triggerMaker", "l1name");
  bool l1fired[l1Names.size()] = {false};

  auto hltNames = getObject<std::vector<std::string>>(ev0, "triggerMaker", "hltname");
  bool hltfired[hltNames.size()] = {false};


  // Branch definition
  tout->Branch("run", &run);
  tout->Branch("lumi", &lumi);
  tout->Branch("evtn", &evtn);
  tout->Branch("passL1", &passL1);

  std::cout << "Found the following L1 seeds:" << std::endl;
  std::cout << " - Total number of seeds " << l1Names.size() << std::endl;
  for (unsigned int iL1=0; iL1<l1Names.size(); ++iL1) {
    tout->Branch(TString(l1Names[iL1]), &l1fired[iL1]);
    std::cout << TString(l1Names[iL1]) << std::endl;
  }

  std::cout << "Found the following HLT paths:" << std::endl;
  std::cout << " - Total number of paths " << hltNames.size() << std::endl;
  for (unsigned int iHLT=0; iHLT<hltNames.size(); ++iHLT) {
    tout->Branch(TString(hltNames[iHLT]), &hltfired[iHLT]);
    std::cout << TString(hltNames[iHLT]) << std::endl;
  }


  tout->Branch("passHLT", &passHLT);

  tout->Branch("wPU", &wPU);
  tout->Branch("nPUTrue", &nPUTrue);

  tout->Branch("nPV", &nPV);
  tout->Branch("PV_x", &PV_x);
  tout->Branch("PV_y", &PV_y);
  tout->Branch("PV_z", &PV_z);

  //Gen particle branches
  tout->Branch("GenPart_pt", &GenParts.pt);
  tout->Branch("GenPart_eta", &GenParts.eta);
  tout->Branch("GenPart_phi", &GenParts.phi);
  tout->Branch("GenPart_m", &GenParts.m);
  tout->Branch("GenPart_vx", &GenParts.vx);
  tout->Branch("GenPart_vy", &GenParts.vy);
  tout->Branch("GenPart_vz", &GenParts.vz);
  tout->Branch("GenPart_lxy", &GenParts.lxy);
  tout->Branch("GenPart_status", &GenParts.status);
  tout->Branch("GenPart_index", &GenParts.index);
  tout->Branch("GenPart_pdgId", &GenParts.pdgId);
  tout->Branch("GenPart_motherIndex", &GenParts.motherIndex);
  tout->Branch("GenPart_motherPdgId", &GenParts.motherPdgId);
  tout->Branch("GenPart_ct", &GenParts.ct);
  tout->Branch("GenB_pt", &GenB_pt);
  tout->Branch("GenB_eta", &GenB_eta);
  tout->Branch("nGenBplus", &nGenBplus);
  tout->Branch("GenBplus_ct", &GenBplus_ct);
  tout->Branch("GenBplus_type", &GenBplus_type);

  //NoVtx SV branches
  tout->Branch("nSVs_NoVtx", &nSVs_NoVtx);
  tout->Branch("SV_index", &SVs.index);
  tout->Branch("SV_ndof", &SVs.ndof);
  tout->Branch("SV_x", &SVs.x);
  tout->Branch("SV_y", &SVs.y);
  tout->Branch("SV_z", &SVs.z);
  tout->Branch("SV_xe", &SVs.xe);
  tout->Branch("SV_ye", &SVs.ye);
  tout->Branch("SV_ze", &SVs.ze);
  tout->Branch("SV_chi2", &SVs.chi2);
  tout->Branch("SV_prob", &SVs.prob);
  tout->Branch("SV_chi2Ndof", &SVs.chi2Ndof);
  tout->Branch("SV_lxy", &SVs.lxy);
  tout->Branch("SV_l3d", &SVs.l3d);
  tout->Branch("SV_mindx", &SVs.mindx);
  tout->Branch("SV_mindy", &SVs.mindy);
  tout->Branch("SV_mindz", &SVs.mindz);
  tout->Branch("SV_mindxy", &SVs.mindxy);
  tout->Branch("SV_mind3d", &SVs.mind3d);
  tout->Branch("SV_maxdx", &SVs.maxdx);
  tout->Branch("SV_maxdy", &SVs.maxdy);
  tout->Branch("SV_maxdz", &SVs.maxdz);
  tout->Branch("SV_maxdxy", &SVs.maxdxy);
  tout->Branch("SV_maxd3d", &SVs.maxd3d);
  tout->Branch("SV_selected", &SVs.selected);
  tout->Branch("SV_onModule", &SVs.onModule);
  tout->Branch("SV_onModuleWithinUnc", &SVs.onModuleWithinUnc);
  tout->Branch("SV_closestDet_x", &SVs.closestDet_x);
  tout->Branch("SV_closestDet_y", &SVs.closestDet_y);
  tout->Branch("SV_closestDet_z", &SVs.closestDet_z);
  tout->Branch("SV_minDistanceFromDet", &SVs.minDistanceFromDet);
  tout->Branch("SV_minDistanceFromDet_x", &SVs.minDistanceFromDet_x);
  tout->Branch("SV_minDistanceFromDet_y", &SVs.minDistanceFromDet_y);
  tout->Branch("SV_minDistanceFromDet_z", &SVs.minDistanceFromDet_z);

  //NoVtx OSV branches
  tout->Branch("SVOverlap_vtxIdxs", &SVOverlaps.vtxIdxs);
  tout->Branch("SVOverlap_x", &SVOverlaps.x);
  tout->Branch("SVOverlap_y", &SVOverlaps.y);
  tout->Branch("SVOverlap_z", &SVOverlaps.z);
  tout->Branch("SVOverlap_lxy", &SVOverlaps.lxy);
  tout->Branch("SVOverlap_l3d", &SVOverlaps.l3d);

  //NoVtx Muon branches
  tout->Branch("nMuons_NoVtx", &nMuons_NoVtx);
  tout->Branch("nMuon_Assoc", &nMuon_Assoc);
  tout->Branch("nMuon_vtx_Assoc", &nMuon_vtx_Assoc);
  tout->Branch("Muon_vtxIdxs", &Muons.vtxIdxs);
  tout->Branch("Muon_saHits", &Muons.saHits);
  tout->Branch("Muon_saMatchedStats", &Muons.saMatchedStats);
  tout->Branch("Muon_muHits", &Muons.muHits);
  tout->Branch("Muon_muChambs", &Muons.muChambs);
  tout->Branch("Muon_muCSCDT", &Muons.muCSCDT);
  tout->Branch("Muon_muMatch", &Muons.muMatch);
  tout->Branch("Muon_muMatchedStats", &Muons.muMatchedStats);
  tout->Branch("Muon_muExpMatchedStats", &Muons.muExpMatchedStats);
  tout->Branch("Muon_muMatchedRPC", &Muons.muMatchedRPC);
  tout->Branch("Muon_pixHits", &Muons.pixHits);
  tout->Branch("Muon_stripHits", &Muons.stripHits);
  tout->Branch("Muon_pixLayers", &Muons.pixLayers);
  tout->Branch("Muon_trkLayers", &Muons.trkLayers);
  tout->Branch("Muon_pt", &Muons.pt);
  tout->Branch("Muon_eta", &Muons.eta);
  tout->Branch("Muon_phi", &Muons.phi);
  tout->Branch("Muon_ch", &Muons.ch);
  tout->Branch("Muon_bestAssocSVIdx", &Muons.bestAssocSVIdx);
  tout->Branch("Muon_bestAssocSVOverlapIdx", &Muons.bestAssocSVOverlapIdx);
  tout->Branch("Muon_isGlobal", &Muons.isGlobal);
  tout->Branch("Muon_isTracker", &Muons.isTracker);
  tout->Branch("Muon_isStandAlone", &Muons.isStandAlone);
  tout->Branch("Muon_chi2Ndof", &Muons.chi2Ndof);
  tout->Branch("Muon_ecalIso", &Muons.ecalIso);
  tout->Branch("Muon_hcalIso", &Muons.hcalIso);
  tout->Branch("Muon_trackIso", &Muons.trackIso);
  tout->Branch("Muon_ecalRelIso", &Muons.ecalRelIso);
  tout->Branch("Muon_hcalRelIso", &Muons.hcalRelIso);
  tout->Branch("Muon_trackRelIso", &Muons.trackRelIso);
  tout->Branch("Muon_dxy", &Muons.dxy);
  tout->Branch("Muon_dxye", &Muons.dxye);
  tout->Branch("Muon_dz", &Muons.dz);
  tout->Branch("Muon_dze", &Muons.dze);
  tout->Branch("Muon_dxysig", &Muons.dxysig);
  tout->Branch("Muon_dzsig", &Muons.dzsig);
  tout->Branch("Muon_phiCorr", &Muons.phiCorr);
  tout->Branch("Muon_dxyCorr", &Muons.dxyCorr);
  tout->Branch("Muon_PFIsoChg0p3", &Muons.PFIsoChg0p3);
  tout->Branch("Muon_PFIsoAll0p3", &Muons.PFIsoAll0p3);
  tout->Branch("Muon_PFRelIsoChg0p3", &Muons.PFRelIsoChg0p3);
  tout->Branch("Muon_PFRelIsoAll0p3", &Muons.PFRelIsoAll0p3);
  tout->Branch("Muon_mindrPF0p3", &Muons.mindrPF0p3);
  tout->Branch("Muon_PFIsoChg0p4", &Muons.PFIsoChg0p4);
  tout->Branch("Muon_PFIsoAll0p4", &Muons.PFIsoAll0p4);
  tout->Branch("Muon_PFRelIsoChg0p4", &Muons.PFRelIsoChg0p4);
  tout->Branch("Muon_PFRelIsoAll0p4", &Muons.PFRelIsoAll0p4);
  tout->Branch("Muon_mindrPF0p4", &Muons.mindrPF0p4);
  tout->Branch("Muon_mindr", &Muons.mindr);
  tout->Branch("Muon_maxdr", &Muons.maxdr);
  tout->Branch("Muon_mindrJet", &Muons.mindrJet);
  tout->Branch("Muon_mindphiJet", &Muons.mindphiJet);
  tout->Branch("Muon_mindetaJet", &Muons.mindetaJet);
  tout->Branch("Muon_vec", &Muons.vec);
  tout->Branch("Muon_selected", &Muons.selected);
  tout->Branch("Muon_nhitsbeforesv", &Muons.nhitsbeforesv);
  tout->Branch("Muon_ncompatible", &Muons.ncompatible);
  tout->Branch("Muon_ncompatibletotal", &Muons.ncompatibletotal);
  tout->Branch("Muon_nexpectedhits", &Muons.nexpectedhits);
  tout->Branch("Muon_nexpectedhitsmultiple", &Muons.nexpectedhitsmultiple);
  tout->Branch("Muon_nexpectedhitsmultipletotal", &Muons.nexpectedhitsmultipletotal);
  tout->Branch("Muon_nexpectedhitstotal", &Muons.nexpectedhitstotal);

  //Vtx SV branches (2024 only)
  tout->Branch("nSVs_Vtx", &nSVs_Vtx);
  tout->Branch("SV_vtx_ndof", &SVsVtx.ndof);
  tout->Branch("SV_vtx_x", &SVsVtx.x);
  tout->Branch("SV_vtx_y", &SVsVtx.y);
  tout->Branch("SV_vtx_z", &SVsVtx.z);
  tout->Branch("SV_vtx_xe", &SVsVtx.xe);
  tout->Branch("SV_vtx_ye", &SVsVtx.ye);
  tout->Branch("SV_vtx_ze", &SVsVtx.ze);
  tout->Branch("SV_vtx_chi2", &SVsVtx.chi2);
  tout->Branch("SV_vtx_prob", &SVsVtx.prob);
  tout->Branch("SV_vtx_chi2Ndof", &SVsVtx.chi2Ndof);
  tout->Branch("SV_vtx_lxy", &SVsVtx.lxy);
  tout->Branch("SV_vtx_l3d", &SVsVtx.l3d);
  tout->Branch("SV_vtx_selected", &SVsVtx.selected);
  tout->Branch("SV_vtx_onModule", &SVsVtx.onModule);
  tout->Branch("SV_vtx_onModuleWithinUnc", &SVsVtx.onModuleWithinUnc);
  tout->Branch("SV_vtx_closestDet_x", &SVsVtx.closestDet_x);
  tout->Branch("SV_vtx_closestDet_y", &SVsVtx.closestDet_y);
  tout->Branch("SV_vtx_closestDet_z", &SVsVtx.closestDet_z);
  tout->Branch("SV_vtx_minDistanceFromDet", &SVsVtx.minDistanceFromDet);
  tout->Branch("SV_vtx_minDistanceFromDet_x", &SVsVtx.minDistanceFromDet_x);
  tout->Branch("SV_vtx_minDistanceFromDet_y", &SVsVtx.minDistanceFromDet_y);
  tout->Branch("SV_vtx_minDistanceFromDet_z", &SVsVtx.minDistanceFromDet_z);
  tout->Branch("SV_vtx_index", &SVsVtx.origIdx);
  tout->Branch("SV_vtx_mindx", &SVsVtx.mindx);
  tout->Branch("SV_vtx_mindy", &SVsVtx.mindy);
  tout->Branch("SV_vtx_mindz", &SVsVtx.mindz);
  tout->Branch("SV_vtx_mindxy", &SVsVtx.mindxy);
  tout->Branch("SV_vtx_mind3d", &SVsVtx.mind3d);
  tout->Branch("SV_vtx_maxdx", &SVsVtx.maxdx);
  tout->Branch("SV_vtx_maxdy", &SVsVtx.maxdy);
  tout->Branch("SV_vtx_maxdz", &SVsVtx.maxdz);
  tout->Branch("SV_vtx_maxdxy", &SVsVtx.maxdxy);
  tout->Branch("SV_vtx_maxd3d", &SVsVtx.maxd3d);

  //Vtx OSV branches
  tout->Branch("SVOverlap_vtx_vtxIdxs", &SVOverlapVtxs.vtxIdxs);
  tout->Branch("SVOverlap_vtx_x",   &SVOverlapVtxs.x);
  tout->Branch("SVOverlap_vtx_y",   &SVOverlapVtxs.y);
  tout->Branch("SVOverlap_vtx_z",   &SVOverlapVtxs.z);
  tout->Branch("SVOverlap_vtx_lxy", &SVOverlapVtxs.lxy);
  tout->Branch("SVOverlap_vtx_l3d", &SVOverlapVtxs.l3d);

  //Vtx Muon branches
  tout->Branch("nMuons_Vtx", &nMuons_Vtx);
  tout->Branch("Muon_vtx_vtxIdxs", &MuonsVtx.vtxIdxs);
  tout->Branch("Muon_vtx_saHits", &MuonsVtx.saHits);
  tout->Branch("Muon_vtx_saMatchedStats", &MuonsVtx.saMatchedStats);
  tout->Branch("Muon_vtx_muHits", &MuonsVtx.muHits);
  tout->Branch("Muon_vtx_muChambs", &MuonsVtx.muChambs);
  tout->Branch("Muon_vtx_muCSCDT", &MuonsVtx.muCSCDT);
  tout->Branch("Muon_vtx_muMatch", &MuonsVtx.muMatch);
  tout->Branch("Muon_vtx_muMatchedStats", &MuonsVtx.muMatchedStats);
  tout->Branch("Muon_vtx_muExpMatchedStats", &MuonsVtx.muExpMatchedStats);
  tout->Branch("Muon_vtx_muMatchedRPC", &MuonsVtx.muMatchedRPC);
  tout->Branch("Muon_vtx_pixHits", &MuonsVtx.pixHits);
  tout->Branch("Muon_vtx_stripHits", &MuonsVtx.stripHits);
  tout->Branch("Muon_vtx_pixLayers", &MuonsVtx.pixLayers);
  tout->Branch("Muon_vtx_trkLayers", &MuonsVtx.trkLayers);
  tout->Branch("Muon_vtx_pt", &MuonsVtx.pt);
  tout->Branch("Muon_vtx_eta", &MuonsVtx.eta);
  tout->Branch("Muon_vtx_phi", &MuonsVtx.phi);
  tout->Branch("Muon_vtx_ch", &MuonsVtx.ch);
  tout->Branch("Muon_vtx_isGlobal", &MuonsVtx.isGlobal);
  tout->Branch("Muon_vtx_isTracker", &MuonsVtx.isTracker);
  tout->Branch("Muon_vtx_isStandAlone", &MuonsVtx.isStandAlone);
  tout->Branch("Muon_vtx_chi2Ndof", &MuonsVtx.chi2Ndof);
  tout->Branch("Muon_vtx_ecalIso", &MuonsVtx.ecalIso);
  tout->Branch("Muon_vtx_hcalIso", &MuonsVtx.hcalIso);
  tout->Branch("Muon_vtx_trackIso", &MuonsVtx.trackIso);
  tout->Branch("Muon_vtx_ecalRelIso", &MuonsVtx.ecalRelIso);
  tout->Branch("Muon_vtx_hcalRelIso", &MuonsVtx.hcalRelIso);
  tout->Branch("Muon_vtx_trackRelIso", &MuonsVtx.trackRelIso);
  tout->Branch("Muon_vtx_dxy", &MuonsVtx.dxy);
  tout->Branch("Muon_vtx_dxye", &MuonsVtx.dxye);
  tout->Branch("Muon_vtx_dz", &MuonsVtx.dz);
  tout->Branch("Muon_vtx_dze", &MuonsVtx.dze);
  tout->Branch("Muon_vtx_dxysig", &MuonsVtx.dxysig);
  tout->Branch("Muon_vtx_dzsig", &MuonsVtx.dzsig);
  tout->Branch("Muon_vtx_PFIsoChg0p3", &MuonsVtx.PFIsoChg0p3);
  tout->Branch("Muon_vtx_PFIsoAll0p3", &MuonsVtx.PFIsoAll0p3);
  tout->Branch("Muon_vtx_PFRelIsoChg0p3", &MuonsVtx.PFRelIsoChg0p3);
  tout->Branch("Muon_vtx_PFRelIsoAll0p3", &MuonsVtx.PFRelIsoAll0p3);
  tout->Branch("Muon_vtx_mindrPF0p3", &MuonsVtx.mindrPF0p3);
  tout->Branch("Muon_vtx_PFIsoChg0p4", &MuonsVtx.PFIsoChg0p4);
  tout->Branch("Muon_vtx_PFIsoAll0p4", &MuonsVtx.PFIsoAll0p4);
  tout->Branch("Muon_vtx_PFRelIsoChg0p4", &MuonsVtx.PFRelIsoChg0p4);
  tout->Branch("Muon_vtx_PFRelIsoAll0p4", &MuonsVtx.PFRelIsoAll0p4);
  tout->Branch("Muon_vtx_mindrPF0p4", &MuonsVtx.mindrPF0p4);
  tout->Branch("Muon_vtx_mindr", &MuonsVtx.mindr);
  tout->Branch("Muon_vtx_maxdr", &MuonsVtx.maxdr);
  tout->Branch("Muon_vtx_mindrJet", &MuonsVtx.mindrJet);
  tout->Branch("Muon_vtx_mindphiJet", &MuonsVtx.mindphiJet);
  tout->Branch("Muon_vtx_mindetaJet", &MuonsVtx.mindetaJet);
  tout->Branch("Muon_vtx_vec", &MuonsVtx.vec);
  tout->Branch("Muon_vtx_selected", &MuonsVtx.selected);
  tout->Branch("Muon_vtx_bestAssocSVVtxIdx", &MuonsVtx.bestAssocSVIdx);
  tout->Branch("Muon_vtx_bestAssocSVOverlapVtxIdx", &MuonsVtx.bestAssocSVOverlapIdx);
  tout->Branch("Muon_vtx_nhitsbeforesv", &MuonsVtx.nhitsbeforesv);
  tout->Branch("Muon_vtx_ncompatible", &MuonsVtx.ncompatible);
  tout->Branch("Muon_vtx_ncompatibletotal", &MuonsVtx.ncompatibletotal);
  tout->Branch("Muon_vtx_nexpectedhits", &MuonsVtx.nexpectedhits);
  tout->Branch("Muon_vtx_nexpectedhitsmultiple", &MuonsVtx.nexpectedhitsmultiple);
  tout->Branch("Muon_vtx_nexpectedhitsmultipletotal", &MuonsVtx.nexpectedhitsmultipletotal);
  tout->Branch("Muon_vtx_nexpectedhitstotal", &MuonsVtx.nexpectedhitstotal);
  tout->Branch("Muon_vtx_phiCorr", &MuonsVtx.phiCorr);
  tout->Branch("Muon_vtx_dxyCorr", &MuonsVtx.dxyCorr);


  // Event setup
  TRandom3 rndm_partialUnblinding(42);
  tqdm bar;

  bool isMC = true;
  if ( process.Contains("Data") )
    isMC = false;

  if ( !isMC ) {
    if ( year == "2022" )
      set_goodrun_file_json("../data/Cert_Collisions2022_355100_362760_Golden.json");
    else if ( year == "2023" )
      set_goodrun_file_json("../data/Cert_Collisions2023_366442_370790_Golden.json");
    else if ( year == "2024" )
      set_goodrun_file_json("../data/Cert_Collisions2024_378981_386951_Golden.json");
  }

  edm::LumiReWeighting lumi_weights;
  if (isMC) {
    if (year == "2022" && process.Contains("2022postEE")) {
      lumi_weights = edm::LumiReWeighting("pileup/MCPileupHistogram2022postEE.root", "pileup/DataPileupHistogram2022.root", "pileup", "pileup");
    } else if (year == "2022" && process.Contains("2022")) {
      lumi_weights = edm::LumiReWeighting("pileup/MCPileupHistogram2022.root", "pileup/DataPileupHistogram2022.root", "pileup", "pileup");
    } else if (year == "2023" && process.Contains("2023BPix")) {
      lumi_weights = edm::LumiReWeighting("pileup/MCPileupHistogram2023BPix.root", "pileup/DataPileupHistogram2023.root", "pileup", "pileup");
    } else if (year == "2024") {
      lumi_weights = edm::LumiReWeighting("pileup/MCPileupHistogram2024.root", "pileup/DataPileupHistogram2024.root", "pileup", "pileup");
    } else {
      lumi_weights = edm::LumiReWeighting("pileup/MCPileupHistogram2023.root", "pileup/DataPileupHistogram2023.root", "pileup", "pileup");
    }
  }

  const char* genPartsLabel = (year == "2024" || year == "2025") ? "prunedGenParticles" : "genParticles";
  const char* puInfoLabel   = (year == "2024" || year == "2025") ? "slimmedAddPileupInfo" : "addPileupInfo";

  // File loop
  unsigned int iFile = 1;
  for (auto inputFile : inputFiles) {
    std::cout << "File number: " << iFile << "\t" << inputFile <<  "\n";
    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        std::cout << "File is in zombie state, skipping..." << std::endl;
	continue;
    }
    TTree* evtree = (TTree*)file->Get("Events");
    if (!evtree) { std::cout << "No Events tree in file, skipping..." << std::endl; file->Close(); continue; }

    // MC produced by re-running the HLT privately keeps the un-pruned GEN record and the
    // un-slimmed PU summary (both written by the HLT process), while samples that went
    // through the central MiniAOD step carry the pruned/slimmed ones. Resolve per file so
    // the same looper runs on either.
    if (isMC) {
      genPartsLabel = branchWithPrefix(evtree, "recoGenParticles_prunedGenParticles") ? "prunedGenParticles"
                    : branchWithPrefix(evtree, "recoGenParticles_genParticles")       ? "genParticles"
                    : genPartsLabel;
      puInfoLabel   = branchWithPrefix(evtree, "PileupSummaryInfos_slimmedAddPileupInfo") ? "slimmedAddPileupInfo"
                    : branchWithPrefix(evtree, "PileupSummaryInfos_addPileupInfo")        ? "addPileupInfo"
                    : puInfoLabel;
      if (iFile == 1)
        std::cout << "GEN collection: " << genPartsLabel << ", PU collection: " << puInfoLabel << std::endl;
    }
    auto nEventsFile = evtree->GetEntries();
    std::cout << "Input events: " << nEventsFile <<  "\n";
    Event ev(file);

    // Event loop
    unsigned int iEv = 0;
    int nSaved = 0;
    int nGoodRun = 0;
    int nDuplicate = 0;
    int nFraction = 0;
    int nL1 = 0;
    int nHLT = 0;
    int nPreMu = 0;
    int nSV = 0;
    try {
    for (ev.toBegin(); ! ev.atEnd(); ++ev) {
      
      //Clear all variables
      GenParts.clear();
      SVs.clear();
      SVOverlaps.clear();
      Muons.clear();
      SVsVtx.clear();
      SVOverlapVtxs.clear();
      MuonsVtx.clear();
      iEv++;

      bar.progress(iEv, nEventsFile);

      auto evAux = ev.eventAuxiliary();
      auto eID = evAux.id();
      run  = eID.run();
      lumi = eID.luminosityBlock();
      evtn = eID.event();

      // Only for b-hadrons (filter and pt-reweighting)
      GenB_pt = -1.;
      GenB_eta = -1.;
      if (isMC  && inputFile.Contains("BToPhi_MPhi") && filterByB) {
        //std::cout << "BToPhi sample identified: applying 1 b-hadron filtering\n"; // Uncomment only for testing
        auto genparts_forB = getObject<std::vector<reco::GenParticle>>(ev, genPartsLabel, "");
        int nbhadron = 0;
        reco::GenParticle bhadron;
        for (unsigned int iGen=0; iGen<genparts_forB.size(); iGen++) {
          auto candidate = genparts_forB[iGen];
          if (!candidate.isLastCopy())
            continue;
          if (abs(candidate.pdgId())!=521 && 
              abs(candidate.pdgId())!=511 &&
              abs(candidate.pdgId())!=531 &&
              abs(candidate.pdgId())!=541 &&
              abs(candidate.pdgId())!=5122)
                continue;
          for (unsigned int jGen=0; jGen<genparts_forB.size(); jGen++) {
            auto dcand = genparts_forB[jGen];
            if (!dcand.isLastCopy() || abs(dcand.pdgId())!=6000211)
              continue;
            if (dcand.motherRef().index() == iGen) {
              bhadron = genparts_forB[iGen]; // last b-hadron identified
              nbhadron++;
              break;
            }
          }
        }
        //
        //std::cout << nbhadron << "\t" << bhadron.pt() << std::endl;
        if (nbhadron==1) {
          // last b-hadron is the only b-hadron
          if (bhadron.pt() < 5. || abs(bhadron.eta()) > 2.8) {
            continue;
          } else {
            // We keep the event
            //std::cout << "We keep the event, nbhadron = " << nbhadron << ", pt = " << bhadron.pt() << std::endl;
            GenB_pt = bhadron.pt();
            GenB_eta = bhadron.eta();
          }
        } else {
          continue;
        }

      }
      
      if (isMC && process.Contains("_Par_")) {
        auto genparts_forMuMu = getObject<std::vector<reco::GenParticle>>(ev, genPartsLabel, "");
        int nMuon = 0;
        int nQualifying = 0;
        for (unsigned int iGen=0; iGen<genparts_forMuMu.size(); iGen++) {
          auto candidate = genparts_forMuMu[iGen];
          if (candidate.status()!=1 || candidate.pt()<=2.0 || abs(candidate.eta())>=2.5)
            continue;
          nQualifying++;
          if (abs(candidate.pdgId())==13)
            nMuon++;
        }
        if (nMuon<1 || nQualifying<2)
          continue;
      }
      
      //
      if (isMC){
        //genWeight = (float) genEvtInfo->weight(); // EventInfo not available
        //sum2Weights->Fill(0.5, genWeight*genWeight);
        sum2Weights->Fill(0.5);
      }
        counts->Fill(0.5);

      // JSON and duplicate removal
      if ( !isMC ) {
        if ( !(goodrun(run, lumi)) )
          continue;
        nGoodRun++;
        duplicate_removal::DorkyEventIdentifier id(run, evtn, lumi);
        if ( is_duplicate(id) )
          continue;
        nDuplicate++;
        if ( doPartialUnblinding && rndm_partialUnblinding.Rndm() > partialUnblindingPercentage )
          continue;
        nFraction++;
      }

      // PU reweighting
      nPU = -1;
      nPUTrue = -1;
      wPU = 1.;
      if (isMC) {
        auto puInfoH = getObjectSafe<std::vector<PileupSummaryInfo>>(ev, puInfoLabel);
        if (puInfoH.size() > 0) {
          for(size_t i=0;i<puInfoH.size();++i) {
            if( puInfoH.at(i).getBunchCrossing() == 0) {
                 nPU = puInfoH.at(i).getPU_NumInteractions();
                 nPUTrue = puInfoH.at(i).getTrueNumInteractions();
                 continue;
            }
          }
          wPU = lumi_weights.weight(nPUTrue);
        }
      }

      // L1 selection
      auto l1s = getObjectSafe<std::vector<bool>>(ev, "triggerMaker", "l1result");
      auto l1Prescales = getObjectSafe<std::vector<double>>(ev, "triggerMaker", "l1prescale");
      passL1 = false;
      for (unsigned int iL1=0; iL1<l1s.size(); ++iL1) {
        l1fired[iL1] = l1s[iL1];
        if (l1s[iL1]==true && l1Prescales[iL1]==1) { // L1 trigger fired and is not prescaled
          passL1 = true;
        }
      }
      if (applyL1 && !passL1)
	      continue;
      nL1++;

      // HLT selection
      auto hlts = getObjectSafe<std::vector<bool>>(ev, "triggerMaker", "hltresult");
      passHLT = false;
      for (auto hlt : hlts) {
        if (hlt==true) {
          passHLT = true;
          break;
        }
      }
      if (!passHLT)
	      continue;
      nHLT++;

      // PV selection
      auto pvs = getObjectSafe<std::vector<Run3ScoutingVertex>>(ev, "hltScoutingPrimaryVertexPacker", "primaryVtx");
      nPV = 0; PV_x = 0; PV_y = 0; PV_z = 0;

      nPV = pvs.size();
      if (pvs.size() > 0 && pvs[0].isValidVtx()) {
        PV_x = pvs[0].x();
        PV_y = pvs[0].y();
        PV_z = pvs[0].z();
      }

      // GenPart
      std::vector<int> matchIndex;
      std::vector<float> cts;
      nGenBplus = 0;
      GenBplus_ct.clear();
      GenBplus_type.clear();
      if (isMC) {
        auto genparts = getObject<std::vector<reco::GenParticle>>(ev, genPartsLabel, "");
        for (unsigned int iGen=0; iGen<genparts.size(); iGen++) {
          auto genpart = genparts[iGen];
          if (abs(genpart.pdgId())!=13 && // Muon
              abs(genpart.pdgId())!=999999 && // Dark photon
              abs(genpart.pdgId())!=9900015 && // Dark photon (inconsistent with datacards)
              abs(genpart.pdgId())!=4900111 && // Dark pion
              abs(genpart.pdgId())!=4900211 && // Dark pion
              abs(genpart.pdgId())!=4900221 && // Dark eta
              abs(genpart.pdgId())!=4900113 && // Dark rho
              abs(genpart.pdgId())!=4900213 && // Dark rho
              abs(genpart.pdgId())!=4900101 && // Dark mass
              abs(genpart.pdgId())!=4900102 && // Dark mass
              abs(genpart.pdgId())!=1023 && // Dark photon from Higgs
              abs(genpart.pdgId())!=25 && // Higgs
              abs(genpart.pdgId())!=6000211 && // Scalar from b-hadrons
              abs(genpart.pdgId())!=443 && // JPsi
              abs(genpart.pdgId())!=521) // BPlus
            continue;
          if (!genpart.isLastCopy())
            continue;

          int motherIdx = -1, motherPdgId = 0;
          reco::GenParticle firstCopy=genpart; // Default value
          if (genpart.numberOfMothers() > 0) {
            motherIdx = genpart.motherRef().index();
            while (genparts[motherIdx].pdgId() == genpart.pdgId()) {
              firstCopy = genparts[motherIdx];
              motherIdx = genparts[motherIdx].motherRef().index();
            }
            motherPdgId = genparts[motherIdx].pdgId();
          }
          GenParts.pt.push_back(genpart.pt());
          GenParts.eta.push_back(genpart.eta());
          GenParts.phi.push_back(genpart.phi());
          GenParts.m.push_back(genpart.mass());
          GenParts.vx.push_back(genpart.vx());
          GenParts.vy.push_back(genpart.vy());
          GenParts.vz.push_back(genpart.vz());
          GenParts.lxy.push_back(TMath::Hypot(genpart.vx(), genpart.vy()));
          GenParts.status.push_back(genpart.status());
          GenParts.index.push_back(iGen);
          GenParts.pdgId.push_back(genpart.pdgId());
          GenParts.motherIndex.push_back(motherIdx);
          GenParts.motherPdgId.push_back(motherPdgId);
          // Generator-level information for lifetime reweighting
          if (genpart.pdgId()==13) {
            if (motherPdgId==1023 || motherPdgId==9900015 || motherPdgId==6000211) {
              float vx1 = firstCopy.vx() * 10.; // mm
              float vy1 = firstCopy.vy() * 10.; // mm
              float vx0 = genparts[motherIdx].vx() * 10.; // mm
              float vy0 = genparts[motherIdx].vy() * 10.; // mm
              matchIndex.push_back(motherIdx);
              cts.push_back ( ((vx1 - vx0)*genparts[motherIdx].px() + (vy1 - vy0)*genparts[motherIdx].py())*genparts[motherIdx].mass()/(genparts[motherIdx].pt()*genparts[motherIdx].pt()) ); // in mm
            }
          }
          // Generator-level information for lifetime reweighting in special B+ case
          if (abs(genpart.pdgId())==521 && genpart.numberOfDaughters()>0 && abs(genparts[genpart.daughterRef(0).index()].pdgId())!=521) {
            nGenBplus++;
            float vx0 = genpart.vx() * 10.; // mm
            float vy0 = genpart.vy() * 10.; // mm
            float vx1 = genparts[genpart.daughterRef(0).index()].vx() * 10.; // mm
            float vy1 = genparts[genpart.daughterRef(0).index()].vy() * 10.; // mm
            GenBplus_ct.push_back(((vx1 - vx0)*genpart.px() + (vy1 - vy0)*genpart.py())*genpart.mass()/(genpart.pt()*genpart.pt()));
            std::vector<int> daughterPdgIds;
            int BPlus_type = 3;
            for (unsigned int id=0; id<genpart.numberOfDaughters(); id++) {
              daughterPdgIds.push_back(abs(genparts[genpart.daughterRef(id).index()].pdgId()));
            }
            if ((std::find(daughterPdgIds.begin(), daughterPdgIds.end(), 11) != daughterPdgIds.end()) || (std::find(daughterPdgIds.begin(), daughterPdgIds.end(), 13) != daughterPdgIds.end()) || (std::find(daughterPdgIds.begin(), daughterPdgIds.end(), 15) != daughterPdgIds.end()))
              BPlus_type = 0;
            else if (std::find(daughterPdgIds.begin(), daughterPdgIds.end(), 443) != daughterPdgIds.end())
              BPlus_type = 1;
            else if (std::find(daughterPdgIds.begin(), daughterPdgIds.end(), 100443) != daughterPdgIds.end())
              BPlus_type = 2;
            else 
              BPlus_type = 3;
            GenBplus_type.push_back(BPlus_type);
          }
        }
        // Set the ct's values for the genParticles
        GenParts.ct = std::vector<float>(GenParts.pt.size(), 1.0);
        for (int jGen = 0; jGen < matchIndex.size(); jGen++) {
            auto it = std::find(GenParts.index.begin(), GenParts.index.end(), matchIndex.at(jGen));
            if (it != GenParts.index.end()) {
                auto position = std::distance(GenParts.index.begin(), it); // Calculate index
                GenParts.ct[position] = cts.at(jGen);
            }
        }
      }

      // NoVtx SV loop
      const char* svPackerLabel  = (year == "2024") ? "hltScoutingMuonPackerNoVtx" : "hltScoutingMuonPacker";
      const char* hitMakerLabel  = (year == "2024") ? "hitMakerNoVtx" : "hitMaker";
      auto svs = getObjectSafe<std::vector<Run3ScoutingVertex>>(ev, svPackerLabel, "displacedVtx");
      auto dvonmodule          = svs.size() > 0 ? getObject<std::vector<bool>>(ev, hitMakerLabel, "dvonmodule") : std::vector<bool>{};
      auto dvonmodulewithinunc = svs.size() > 0 ? getObject<std::vector<bool>>(ev, hitMakerLabel, "dvonmodulewithinunc") : std::vector<bool>{};
      auto dvdetxmind          = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvdetxmind") : std::vector<float>{};
      auto dvdetymind          = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvdetymind") : std::vector<float>{};
      auto dvdetzmind          = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvdetzmind") : std::vector<float>{};
      auto dvmindfromdet       = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvmindfromdet") : std::vector<float>{};
      auto dvmindfromdetx      = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvmindfromdetx") : std::vector<float>{};
      auto dvmindfromdety      = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvmindfromdety") : std::vector<float>{};
      auto dvmindfromdetz      = svs.size() > 0 ? getObject<std::vector<float>>(ev, hitMakerLabel, "dvmindfromdetz") : std::vector<float>{};

      unsigned int nSVs = svs.size();
      for (unsigned int iSV=0; iSV<nSVs; ++iSV) {
        auto sv = svs[iSV];
        if (!(sv.isValidVtx()))
          continue;

        float x=sv.x(), y=sv.y(), z=sv.z();
        float xe=sv.xError(), ye=sv.yError(), ze=sv.zError();
        float chi2=sv.chi2(), ndof=sv.ndof();

        SVs.index.push_back(iSV);
        SVs.ndof.push_back(ndof);
        SVs.x.push_back(x);
        SVs.y.push_back(y);
        SVs.z.push_back(z);
        SVs.xe.push_back(xe);
        SVs.ye.push_back(ye);
        SVs.ze.push_back(ze);
        SVs.chi2.push_back(chi2);
        SVs.prob.push_back(TMath::Prob(chi2, ndof));
        SVs.chi2Ndof.push_back(chi2/ndof);
        float lxy = TMath::Sqrt((x-PV_x)*(x-PV_x)+(y-PV_y)*(y-PV_y));
        SVs.lxy.push_back(lxy);
        SVs.l3d.push_back(TMath::Sqrt(lxy*lxy+(z-PV_z)*(z-PV_z)));
        SVs.onModule.push_back(dvonmodule.at(iSV));
        SVs.onModuleWithinUnc.push_back(dvonmodulewithinunc.at(iSV));
        SVs.closestDet_x.push_back(dvdetxmind.at(iSV));
        SVs.closestDet_y.push_back(dvdetymind.at(iSV));
        SVs.closestDet_z.push_back(dvdetzmind.at(iSV));
        SVs.minDistanceFromDet.push_back(dvmindfromdet.at(iSV));
        SVs.minDistanceFromDet_x.push_back(dvmindfromdetx.at(iSV));
        SVs.minDistanceFromDet_y.push_back(dvmindfromdety.at(iSV));
        SVs.minDistanceFromDet_z.push_back(dvmindfromdetz.at(iSV));

        float mindx=1e6, mindy=1e6, mindz=1e6, mindxy=1e6, mind3d=1e6;
        float maxdx=-1, maxdy=-1, maxdz=-1, maxdxy=-1, maxd3d=-1;
        for (unsigned int jSV=0; jSV<nSVs; ++jSV) {
          if (jSV==iSV)
            continue;
          auto svOther = svs[jSV];
          if (!(svOther.isValidVtx()))
            continue;

          float xOther=svOther.x(), yOther=svOther.y(), zOther=svOther.z();
          float dx = (x-xOther)*(x-xOther); // Squared for the moment
          float dy = (y-yOther)*(y-yOther); // Squared for the moment
          float dz = (z-zOther)*(z-zOther); // Squared for the moment
          float dxy = TMath::Sqrt(dx+dy);
          float d3d = TMath::Sqrt(dx+dy+dz);
          dx = TMath::Sqrt(dx); // Back to the proper value
          dy = TMath::Sqrt(dy); // Back to the proper value
          dz = TMath::Sqrt(dz); // Back to the proper value

          if (dx<mindx) mindx = dx;
          if (dx>maxdx) maxdx = dx;
          if (dy<mindy) mindy = dy;
          if (dy>maxdy) maxdy = dy;
          if (dz<mindz) mindz = dz;
          if (dz>maxdz) maxdz = dz;
          if (dxy<mindxy) mindxy = dxy;
          if (dxy>maxdxy) maxdxy = dxy;
          if (d3d<mind3d) mind3d = d3d;
          if (d3d>maxd3d) maxd3d = d3d;
        }
        SVs.mindx.push_back(mindx);
        SVs.mindy.push_back(mindy);
        SVs.mindz.push_back(mindz);
        SVs.mindxy.push_back(mindxy);
        SVs.mind3d.push_back(mind3d);
        SVs.maxdx.push_back(maxdx);
        SVs.maxdy.push_back(maxdy);
        SVs.maxdz.push_back(maxdz);
        SVs.maxdxy.push_back(maxdxy);
        SVs.maxd3d.push_back(maxd3d);
        SVs.selected.push_back( (chi2/ndof < 10.0) );
      }
      SVs.sort();
      nSVs_NoVtx = (int)SVs.x.size();

      //Overlapping NoVtx SVs
      for (unsigned int iSV=0; iSV<SVs.x.size(); ++iSV) {
        if (SVs.selected[iSV]==0)
          continue;

        float sumOfProb = 0.0;
        std::vector<unsigned int> vtxIdxs_temp = {};

        bool usedSV=false;
        for (unsigned int iSVOverlap=0; iSVOverlap<SVOverlaps.vtxIdxs.size(); iSVOverlap++) {
          for (auto SVOverlapVtxIdx : SVOverlaps.vtxIdxs[iSVOverlap]) {
            if (iSV==SVOverlapVtxIdx) {
              usedSV=true;
              break;
            }
          }
        }
        if (usedSV)
          continue;

        float x=SVs.x[iSV], y=SVs.y[iSV], z=SVs.z[iSV];
        float xe=SVs.xe[iSV], ye=SVs.ye[iSV], ze=SVs.ze[iSV];
        for (unsigned int jSV=iSV+1; jSV<SVs.x.size(); ++jSV) {
          if (SVs.selected[jSV]==0)
            continue;

          float xOther=SVs.x[jSV], yOther=SVs.y[jSV], zOther=SVs.z[jSV];
          float dx = x - xOther;
          float dy = y - yOther;
          float dz = z - zOther;
	        float xeOther=SVs.xe[jSV], yeOther=SVs.ye[jSV], zeOther=SVs.ze[jSV];
          float dxy = dx*dx + dy*dy;
          float d3d = dxy*dxy + dz*dz;
	        float xeOverlap = std::max(xe, xeOther);
	        float yeOverlap = std::max(ye, yeOther);
	        float zeOverlap = std::max(ze, zeOther);
	        float dxyeOverlap = TMath::Sqrt(xeOverlap*xeOverlap+yeOverlap*yeOverlap);
	        float d3deOverlap = TMath::Sqrt(dxyeOverlap*dxyeOverlap+zeOverlap*zeOverlap);
          if ( fabs(dx)<xeOverlap && fabs(dy)<yeOverlap && fabs(dz)<zeOverlap && TMath::Sqrt(dxy)<std::min(dxyeOverlap, maxDXYerr) && TMath::Sqrt(d3d)<std::min(d3deOverlap, maxD3Derr) ) {
            vtxIdxs_temp.push_back(jSV);
            sumOfProb += SVs.prob[jSV];
          }
        }
        if (vtxIdxs_temp.size() > 0) {
          // Add initial SV to the SVOverlap
          vtxIdxs_temp.insert(vtxIdxs_temp.begin(),iSV);
          sumOfProb += SVs.prob[iSV];

          // Compute and fill SVOverlap properties
          SVOverlaps.vtxIdxs.push_back(vtxIdxs_temp);
          float xOverlap=0.0, yOverlap=0.0, zOverlap=0.0;
          for (auto vtxIdx : vtxIdxs_temp) {
            xOverlap += SVs.prob[vtxIdx]/sumOfProb * SVs.x[vtxIdx];
            yOverlap += SVs.prob[vtxIdx]/sumOfProb * SVs.y[vtxIdx];
            zOverlap += SVs.prob[vtxIdx]/sumOfProb * SVs.z[vtxIdx];
          }
          SVOverlaps.x.push_back(xOverlap);
          SVOverlaps.y.push_back(yOverlap);
          SVOverlaps.z.push_back(zOverlap);
          float lxyOverlap = TMath::Sqrt((xOverlap-PV_x)*(xOverlap-PV_x)+(yOverlap-PV_y)*(yOverlap-PV_y));
          SVOverlaps.lxy.push_back(lxyOverlap);
          SVOverlaps.l3d.push_back(TMath::Sqrt(lxyOverlap*lxyOverlap+(zOverlap-PV_z)*(zOverlap-PV_z)));
        }
      }

      // NoVtx Muon selection
      auto musnoVtx = (year == "2024") ?
        getObjectSafe<std::vector<Run3ScoutingMuon>>(ev, "hltScoutingMuonPackerNoVtx") :
        getObjectSafe<std::vector<Run3ScoutingMuon>>(ev, "hltScoutingMuonPacker");
      auto jets = getObjectSafe<std::vector<Run3ScoutingPFJet>>(ev, "hltScoutingPFPacker");
      auto pfs = getObjectSafe<std::vector<Run3ScoutingParticle>>(ev, "hltScoutingPFPacker");
      auto nhitsbeforesv = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "nhitsbeforesv");
      auto ncompatible = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "ncompatible");
      auto ncompatibletotal = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "ncompatibletotal");
      auto nexpectedhits = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "nexpectedhits");
      auto nexpectedhitsmultiple = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "nexpectedhitsmultiple");
      auto nexpectedhitsmultipletotal = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "nexpectedhitsmultipletotal");
      auto nexpectedhitstotal = getObjectSafe<std::vector<std::vector<int>>>(ev, hitMakerLabel, "nexpectedhitstotal");
      std::vector<std::vector<int>> vtxIndxVtx;
      std::vector<std::vector<int>> nhitsbeforesv_vtx, ncompatible_vtx, ncompatibletotal_vtx;
      std::vector<std::vector<int>> nexpectedhits_vtx, nexpectedhitsmultiple_vtx, nexpectedhitsmultipletotal_vtx, nexpectedhitstotal_vtx;
      unsigned int nMusnoVtx = musnoVtx.size();
      nMuons_NoVtx = (int)nMusnoVtx;
      nMuon_Assoc=0;

      for (unsigned int iMu=0; iMu<nMusnoVtx; ++iMu) {
        auto mu = musnoVtx[iMu];
        std::vector<int> matchedAndSelVtxIdxs;
        int bestAssocSVIdx    = -1; // original index — for mu.vtxIndx() comparisons
        int bestAssocSVPosIdx = -1; // position index — for SVs array access
        for (auto matchedVtxIdx : mu.vtxIndx()) {
          for (unsigned int iSV=0; iSV<SVs.index.size(); ++iSV) {
            if (matchedVtxIdx==SVs.index[iSV] && SVs.selected[iSV]) {
              matchedAndSelVtxIdxs.push_back(matchedVtxIdx);
              if (bestAssocSVPosIdx==-1 || SVs.prob[iSV] > SVs.prob[bestAssocSVPosIdx]) {
                bestAssocSVPosIdx = iSV;
                bestAssocSVIdx    = matchedVtxIdx;
              }
            }
          }
        }

        if (matchedAndSelVtxIdxs.size() < 1)
          continue;
        nMuon_Assoc++;

        float pt=mu.pt(), eta=mu.eta(), phi=mu.phi();

        Muons.vtxIdxs.push_back(matchedAndSelVtxIdxs);
        Muons.bestAssocSVIdx.push_back(bestAssocSVIdx);

        int bestAssocSVOverlapIdx=-1;
        for (unsigned int iSVOverlap=0; iSVOverlap<SVOverlaps.vtxIdxs.size(); iSVOverlap++) {
          for (auto SVOverlapVtxIdx : SVOverlaps.vtxIdxs[iSVOverlap]) {
            if (bestAssocSVIdx==SVOverlapVtxIdx)
              bestAssocSVOverlapIdx = iSVOverlap;
          }
        }
        Muons.bestAssocSVOverlapIdx.push_back(bestAssocSVOverlapIdx);

        Muons.saHits.push_back(mu.nValidStandAloneMuonHits());
        Muons.saMatchedStats.push_back(mu.nStandAloneMuonMatchedStations());
        Muons.muHits.push_back(mu.nValidRecoMuonHits());
        Muons.muChambs.push_back(mu.nRecoMuonChambers());
        Muons.muCSCDT.push_back(mu.nRecoMuonChambersCSCorDT());
        Muons.muMatch.push_back(mu.nRecoMuonMatches());
        Muons.muMatchedStats.push_back(mu.nRecoMuonMatchedStations());
        Muons.muExpMatchedStats.push_back(mu.nRecoMuonExpectedMatchedStations());
        Muons.muMatchedRPC.push_back(mu.nRecoMuonMatchedRPCLayers());
        Muons.pixHits.push_back(mu.nValidPixelHits());
        Muons.stripHits.push_back(mu.nValidStripHits());
        Muons.pixLayers.push_back(mu.nPixelLayersWithMeasurement());
        Muons.trkLayers.push_back(mu.nTrackerLayersWithMeasurement());
        Muons.pt.push_back(pt);
        Muons.eta.push_back(eta);
        Muons.phi.push_back(phi);
        Muons.ch.push_back(mu.charge());
        Muons.isGlobal.push_back(isGlobalMuon(mu.type()));
        Muons.isTracker.push_back(isTrackerMuon(mu.type()));
        Muons.isStandAlone.push_back(isStandAloneMuon(mu.type()));
        Muons.chi2Ndof.push_back(mu.normalizedChi2());
        Muons.ecalIso.push_back(mu.ecalIso());
        Muons.hcalIso.push_back(mu.hcalIso());
        Muons.trackIso.push_back(mu.trackIso());
        Muons.ecalRelIso.push_back(mu.ecalIso()/pt);
        Muons.hcalRelIso.push_back(mu.hcalIso()/pt);
        Muons.trackRelIso.push_back(mu.trackIso()/pt);
        Muons.dxy.push_back(mu.trk_dxy());
        Muons.dxye.push_back(mu.trk_dxyError());
        Muons.dz.push_back(mu.trk_dz());
        Muons.dze.push_back(mu.trk_dzError());
        Muons.dxysig.push_back(mu.trk_dxy()/mu.trk_dxyError());
        Muons.dzsig.push_back(mu.trk_dz()/mu.trk_dzError());
        Muons.selected.push_back(pt>3.0 && fabs(eta)<2.4 && mu.normalizedChi2()<3.0);

        for (unsigned int iDV=0; iDV<mu.vtxIndx().size(); ++iDV) {
	        if (mu.vtxIndx().at(iDV)==bestAssocSVIdx) {
            Muons.nhitsbeforesv.push_back(nhitsbeforesv.at(iMu).at(iDV));
            Muons.ncompatible.push_back(ncompatible.at(iMu).at(iDV));
            Muons.ncompatibletotal.push_back(ncompatibletotal.at(iMu).at(iDV));
            Muons.nexpectedhits.push_back(nexpectedhits.at(iMu).at(iDV));
            Muons.nexpectedhitsmultiple.push_back(nexpectedhitsmultiple.at(iMu).at(iDV));
            Muons.nexpectedhitsmultipletotal.push_back(nexpectedhitsmultipletotal.at(iMu).at(iDV));
            Muons.nexpectedhitstotal.push_back(nexpectedhitstotal.at(iMu).at(iDV));
            break;
          }
        }

        float bestSVPosition_x, bestSVPosition_y;
        if (bestAssocSVOverlapIdx!=-1) {
          bestSVPosition_x = SVOverlaps.x[bestAssocSVOverlapIdx];
          bestSVPosition_y = SVOverlaps.y[bestAssocSVOverlapIdx];
        }
        else {
          bestSVPosition_x = SVs.x[bestAssocSVPosIdx];
          bestSVPosition_y = SVs.y[bestAssocSVPosIdx];
        }
        float dxyCorr = -(bestSVPosition_x - PV_x)*TMath::Sin(phi) + (bestSVPosition_y - PV_y)*TMath::Cos(phi);
        Muons.dxyCorr.push_back(dxyCorr);
        TLorentzVector muVec; muVec.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
        float phiCorr = getCorrectedPhi(mu, muVec, dxyCorr, bestSVPosition_x, bestSVPosition_y);
        Muons.phiCorr.push_back(phiCorr);
        muVec.SetPtEtaPhiM(pt, eta, phiCorr, MUON_MASS);
        Muons.vec.push_back(muVec);

        const auto pfIsos0p3 = getPFIsolation(muVec, pfs, 0.3);
        const auto pfIsoChg0p3 = std::get<0>(pfIsos0p3);
        const auto pfIsoAll0p3 = pfIsoChg0p3 + std::max(0.0, std::get<1>(pfIsos0p3)+std::get<2>(pfIsos0p3)-0.5*std::get<3>(pfIsos0p3));
        Muons.PFIsoChg0p3.push_back(pfIsoChg0p3);
        Muons.PFIsoAll0p3.push_back(pfIsoAll0p3);;
        Muons.PFRelIsoChg0p3.push_back(pfIsoChg0p3/pt);
        Muons.PFRelIsoAll0p3.push_back(pfIsoAll0p3/pt);
        Muons.mindrPF0p3.push_back(std::get<4>(pfIsos0p3));

        const auto pfIsos0p4 = getPFIsolation(muVec, pfs, 0.4);
        const auto pfIsoChg0p4 = std::get<0>(pfIsos0p4);
        const auto pfIsoAll0p4 = pfIsoChg0p4 + std::max(0.0, std::get<1>(pfIsos0p4)+std::get<2>(pfIsos0p4)-0.5*std::get<3>(pfIsos0p4));
        Muons.PFIsoChg0p4.push_back(pfIsoChg0p4);
        Muons.PFIsoAll0p4.push_back(pfIsoAll0p4);;
        Muons.PFRelIsoChg0p4.push_back(pfIsoChg0p4/pt);
        Muons.PFRelIsoAll0p4.push_back(pfIsoAll0p4/pt);
        Muons.mindrPF0p4.push_back(std::get<4>(pfIsos0p4));

        float mindr=1e6;
        float maxdr=-1;
        for (unsigned int jMu=0; jMu<nMusnoVtx; ++jMu) {
          if (jMu==iMu)
            continue;
          auto muOther = musnoVtx[jMu];
          bool matchedAndSelVtx = false;
          for (auto matchedVtxIdx : muOther.vtxIndx()) {
            for (auto selVtxIdx : SVs.index) {
              if (matchedVtxIdx==selVtxIdx) {
                matchedAndSelVtx = true;
                break;
              }
              if (matchedAndSelVtx)
                break;
            }
          }
          if (!matchedAndSelVtx)
            continue;
          if (!(fabs(muOther.eta())<2.4))
            continue;
          TLorentzVector muVecOther; muVecOther.SetPtEtaPhiM(muOther.pt(), muOther.eta(), muOther.phi(), MUON_MASS);
          float dr = muVec.DeltaR(muVecOther);
          if (dr<mindr) mindr = dr;
          if (dr>maxdr) maxdr = dr;
        }
        Muons.mindr.push_back(mindr);
        Muons.maxdr.push_back(maxdr);

        float mindrJet=1e6;
        float mindphiJet=1e6;
        float mindetaJet=1e6;
        for (auto jet : jets) {
          if (!(jet.pt()>20.0 && fabs(jet.eta())<3.0))
            continue;
          TLorentzVector jetVec; jetVec.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.m());
          float dr = muVec.DeltaR(jetVec);
          if (dr<mindrJet) {
	    mindrJet = dr;
	    mindphiJet = fabs(muVec.DeltaPhi(jetVec));
	    mindetaJet = fabs(muVec.Eta()-jetVec.Eta());
	        }
        }
        Muons.mindrJet.push_back(mindrJet);
        Muons.mindphiJet.push_back(mindphiJet);
        Muons.mindetaJet.push_back(mindetaJet);
      }

      nMuons_Vtx = 0;
      nSVs_Vtx   = 0;
      if (year == "2024") { //If the year is 2024, process Vtx collection

        //Vtx SVs
        auto svsVtx_raw = getObjectSafe<std::vector<Run3ScoutingVertex>>(ev, "vertexMakerVtx", "verticesVtx");
        vtxIndxVtx = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "vertexMakerVtx", "vtxIndxVtx") : std::vector<std::vector<int>>{};
        nhitsbeforesv_vtx       = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "nhitsbeforesv")             : std::vector<std::vector<int>>{};
        ncompatible_vtx         = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "ncompatible")               : std::vector<std::vector<int>>{};
        ncompatibletotal_vtx    = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "ncompatibletotal")          : std::vector<std::vector<int>>{};
        nexpectedhits_vtx       = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "nexpectedhits")             : std::vector<std::vector<int>>{};
        nexpectedhitsmultiple_vtx         = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "nexpectedhitsmultiple")         : std::vector<std::vector<int>>{};
        nexpectedhitsmultipletotal_vtx    = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "nexpectedhitsmultipletotal")    : std::vector<std::vector<int>>{};
        nexpectedhitstotal_vtx  = svsVtx_raw.size() > 0 ? getObject<std::vector<std::vector<int>>>(ev, "hitMakerVtx", "nexpectedhitstotal")        : std::vector<std::vector<int>>{};
        auto dvonmodule_vtx          = svsVtx_raw.size() > 0 ? getObject<std::vector<bool>>(ev, "hitMakerVtx", "dvonmodule")          : std::vector<bool>{};
        auto dvonmodulewithinunc_vtx = svsVtx_raw.size() > 0 ? getObject<std::vector<bool>>(ev, "hitMakerVtx", "dvonmodulewithinunc") : std::vector<bool>{};
        auto dvdetxmind_vtx          = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvdetxmind")         : std::vector<float>{};
        auto dvdetymind_vtx          = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvdetymind")         : std::vector<float>{};
        auto dvdetzmind_vtx          = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvdetzmind")         : std::vector<float>{};
        auto dvmindfromdet_vtx       = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvmindfromdet")      : std::vector<float>{};
        auto dvmindfromdetx_vtx      = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvmindfromdetx")     : std::vector<float>{};
        auto dvmindfromdety_vtx      = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvmindfromdety")     : std::vector<float>{};
        auto dvmindfromdetz_vtx      = svsVtx_raw.size() > 0 ? getObject<std::vector<float>>(ev, "hitMakerVtx", "dvmindfromdetz")     : std::vector<float>{};

        for (unsigned int iSV = 0; iSV < svsVtx_raw.size(); ++iSV) {
          const auto& sv = svsVtx_raw[iSV];
          if (!sv.isValidVtx()) continue;
          float x=sv.x(), y=sv.y(), z=sv.z();
          float xe=sv.xError(), ye=sv.yError(), ze=sv.zError();
          float chi2=sv.chi2(), ndof=sv.ndof();

          SVsVtx.origIdx.push_back(iSV);
          SVsVtx.ndof.push_back((unsigned int)ndof);
          SVsVtx.x.push_back(x); SVsVtx.y.push_back(y); SVsVtx.z.push_back(z);
          SVsVtx.xe.push_back(xe); SVsVtx.ye.push_back(ye); SVsVtx.ze.push_back(ze);
          SVsVtx.chi2.push_back(chi2);
          SVsVtx.prob.push_back(TMath::Prob(chi2, ndof));
          SVsVtx.chi2Ndof.push_back(chi2/ndof);
          float lxy = TMath::Sqrt((x-PV_x)*(x-PV_x)+(y-PV_y)*(y-PV_y));
          SVsVtx.lxy.push_back(lxy);
          SVsVtx.l3d.push_back(TMath::Sqrt(lxy*lxy+(z-PV_z)*(z-PV_z)));
          SVsVtx.selected.push_back(chi2/ndof < 10.0);
          SVsVtx.onModule.push_back(dvonmodule_vtx.at(iSV));
          SVsVtx.onModuleWithinUnc.push_back(dvonmodulewithinunc_vtx.at(iSV));
          SVsVtx.closestDet_x.push_back(dvdetxmind_vtx.at(iSV));
          SVsVtx.closestDet_y.push_back(dvdetymind_vtx.at(iSV));
          SVsVtx.closestDet_z.push_back(dvdetzmind_vtx.at(iSV));
          SVsVtx.minDistanceFromDet.push_back(dvmindfromdet_vtx.at(iSV));
          SVsVtx.minDistanceFromDet_x.push_back(dvmindfromdetx_vtx.at(iSV));
          SVsVtx.minDistanceFromDet_y.push_back(dvmindfromdety_vtx.at(iSV));
          SVsVtx.minDistanceFromDet_z.push_back(dvmindfromdetz_vtx.at(iSV));
          float mindx_v=1e6, mindy_v=1e6, mindz_v=1e6, mindxy_v=1e6, mind3d_v=1e6;
          float maxdx_v=-1,  maxdy_v=-1,  maxdz_v=-1,  maxdxy_v=-1,  maxd3d_v=-1;

          for (unsigned int jSV=0; jSV<svsVtx_raw.size(); ++jSV) {
            if (jSV==iSV) continue;
            auto svOther = svsVtx_raw[jSV];
            if (!svOther.isValidVtx()) continue;
            float xOther=svOther.x(), yOther=svOther.y(), zOther=svOther.z();
            float dx = (x-xOther)*(x-xOther);
            float dy = (y-yOther)*(y-yOther);
            float dz = (z-zOther)*(z-zOther);
            float dxy = TMath::Sqrt(dx+dy);
            float d3d = TMath::Sqrt(dx+dy+dz);
            dx = TMath::Sqrt(dx);
            dy = TMath::Sqrt(dy);
            dz = TMath::Sqrt(dz);
            if (dx<mindx_v) mindx_v=dx;
            if (dx>maxdx_v) maxdx_v=dx;
            if (dy<mindy_v) mindy_v=dy;
            if (dy>maxdy_v) maxdy_v=dy;
            if (dz<mindz_v) mindz_v=dz;
            if (dz>maxdz_v) maxdz_v=dz;
            if (dxy<mindxy_v) mindxy_v=dxy;
            if (dxy>maxdxy_v) maxdxy_v=dxy;
            if (d3d<mind3d_v) mind3d_v=d3d;
            if (d3d>maxd3d_v) maxd3d_v=d3d;
          }

          SVsVtx.mindx.push_back(mindx_v);
          SVsVtx.mindy.push_back(mindy_v);
          SVsVtx.mindz.push_back(mindz_v);
          SVsVtx.mindxy.push_back(mindxy_v);
          SVsVtx.mind3d.push_back(mind3d_v);
          SVsVtx.maxdx.push_back(maxdx_v);
          SVsVtx.maxdy.push_back(maxdy_v);
          SVsVtx.maxdz.push_back(maxdz_v);
          SVsVtx.maxdxy.push_back(maxdxy_v);
          SVsVtx.maxd3d.push_back(maxd3d_v);
        }
        SVsVtx.sort();
        nSVs_Vtx = (int)SVsVtx.x.size();

        //Overlapping Vtx SVs
        for (unsigned int iSV=0; iSV<SVsVtx.x.size(); ++iSV) {
          if (SVsVtx.selected[iSV]==0)
            continue;
          float sumOfProb = 0.0;
          std::vector<unsigned int> vtxIdxs_temp = {};
          bool usedSV=false;
          for (unsigned int iSVOverlap=0; iSVOverlap<SVOverlapVtxs.vtxIdxs.size(); iSVOverlap++) {
            for (auto SVOverlapVtxIdx : SVOverlapVtxs.vtxIdxs[iSVOverlap]) {
              if (iSV==SVOverlapVtxIdx) {
                usedSV=true;
                break;
              }
            }
          }
          if (usedSV)
            continue;
          float x=SVsVtx.x[iSV], y=SVsVtx.y[iSV], z=SVsVtx.z[iSV];
          float xe=SVsVtx.xe[iSV], ye=SVsVtx.ye[iSV], ze=SVsVtx.ze[iSV];
          for (unsigned int jSV=iSV+1; jSV<SVsVtx.x.size(); ++jSV) {
            if (SVsVtx.selected[jSV]==0)
              continue;
            float xOther=SVsVtx.x[jSV], yOther=SVsVtx.y[jSV], zOther=SVsVtx.z[jSV];
            float dx = x - xOther;
            float dy = y - yOther;
            float dz = z - zOther;
            float xeOther=SVsVtx.xe[jSV], yeOther=SVsVtx.ye[jSV], zeOther=SVsVtx.ze[jSV];
            float dxy = dx*dx + dy*dy;
            float d3d = dxy*dxy + dz*dz;
            float xeOverlap = std::max(xe, xeOther);
            float yeOverlap = std::max(ye, yeOther);
            float zeOverlap = std::max(ze, zeOther);
            float dxyeOverlap = TMath::Sqrt(xeOverlap*xeOverlap+yeOverlap*yeOverlap);
            float d3deOverlap = TMath::Sqrt(dxyeOverlap*dxyeOverlap+zeOverlap*zeOverlap);
            if ( fabs(dx)<xeOverlap && fabs(dy)<yeOverlap && fabs(dz)<zeOverlap && TMath::Sqrt(dxy)<std::min(dxyeOverlap, maxDXYerr) && TMath::Sqrt(d3d)<std::min(d3deOverlap, maxD3Derr) ) {
              vtxIdxs_temp.push_back(jSV);
              sumOfProb += SVsVtx.prob[jSV];
            }
          }

          if (vtxIdxs_temp.size() > 0) {
            vtxIdxs_temp.insert(vtxIdxs_temp.begin(), iSV);
            sumOfProb += SVsVtx.prob[iSV];
            SVOverlapVtxs.vtxIdxs.push_back(vtxIdxs_temp);
            float xOverlap=0.0, yOverlap=0.0, zOverlap=0.0;
            for (auto vtxIdx : vtxIdxs_temp) {
              xOverlap += SVsVtx.prob[vtxIdx]/sumOfProb * SVsVtx.x[vtxIdx];
              yOverlap += SVsVtx.prob[vtxIdx]/sumOfProb * SVsVtx.y[vtxIdx];
              zOverlap += SVsVtx.prob[vtxIdx]/sumOfProb * SVsVtx.z[vtxIdx];
            }
            SVOverlapVtxs.x.push_back(xOverlap);
            SVOverlapVtxs.y.push_back(yOverlap);
            SVOverlapVtxs.z.push_back(zOverlap);
            float lxyOverlap = TMath::Sqrt((xOverlap-PV_x)*(xOverlap-PV_x)+(yOverlap-PV_y)*(yOverlap-PV_y));
            SVOverlapVtxs.lxy.push_back(lxyOverlap);
            SVOverlapVtxs.l3d.push_back(TMath::Sqrt(lxyOverlap*lxyOverlap+(zOverlap-PV_z)*(zOverlap-PV_z)));
          }
        }

        // Muon Vtx collection
        auto musVtx = getObjectSafe<std::vector<Run3ScoutingMuon>>(ev, "hltScoutingMuonPackerVtx");

        nMuon_vtx_Assoc = 0;
        unsigned int nMusVtx = musVtx.size();
        nMuons_Vtx = (int)nMusVtx;
        for (unsigned int iMu=0; iMu<nMusVtx; ++iMu) {
          auto mu = musVtx[iMu];

          float pt=mu.pt(), eta=mu.eta(), phi=mu.phi();
          TLorentzVector muVec; muVec.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);

          std::vector<int> vtxIdxsVec(mu.vtxIndx().begin(), mu.vtxIndx().end());
          MuonsVtx.vtxIdxs.push_back(vtxIdxsVec);

          MuonsVtx.saHits.push_back(mu.nValidStandAloneMuonHits());
          MuonsVtx.saMatchedStats.push_back(mu.nStandAloneMuonMatchedStations());
          MuonsVtx.muHits.push_back(mu.nValidRecoMuonHits());
          MuonsVtx.muChambs.push_back(mu.nRecoMuonChambers());
          MuonsVtx.muCSCDT.push_back(mu.nRecoMuonChambersCSCorDT());
          MuonsVtx.muMatch.push_back(mu.nRecoMuonMatches());
          MuonsVtx.muMatchedStats.push_back(mu.nRecoMuonMatchedStations());
          MuonsVtx.muExpMatchedStats.push_back(mu.nRecoMuonExpectedMatchedStations());
          MuonsVtx.muMatchedRPC.push_back(mu.nRecoMuonMatchedRPCLayers());
          MuonsVtx.pixHits.push_back(mu.nValidPixelHits());
          MuonsVtx.stripHits.push_back(mu.nValidStripHits());
          MuonsVtx.pixLayers.push_back(mu.nPixelLayersWithMeasurement());
          MuonsVtx.trkLayers.push_back(mu.nTrackerLayersWithMeasurement());
          MuonsVtx.pt.push_back(pt);
          MuonsVtx.eta.push_back(eta);
          MuonsVtx.phi.push_back(phi);
          MuonsVtx.ch.push_back(mu.charge());
          MuonsVtx.isGlobal.push_back(isGlobalMuon(mu.type()));
          MuonsVtx.isTracker.push_back(isTrackerMuon(mu.type()));
          MuonsVtx.isStandAlone.push_back(isStandAloneMuon(mu.type()));
          MuonsVtx.chi2Ndof.push_back(mu.normalizedChi2());
          MuonsVtx.ecalIso.push_back(mu.ecalIso());
          MuonsVtx.hcalIso.push_back(mu.hcalIso());
          MuonsVtx.trackIso.push_back(mu.trackIso());
          MuonsVtx.ecalRelIso.push_back(mu.ecalIso()/pt);
          MuonsVtx.hcalRelIso.push_back(mu.hcalIso()/pt);
          MuonsVtx.trackRelIso.push_back(mu.trackIso()/pt);
          MuonsVtx.dxy.push_back(mu.trk_dxy());
          MuonsVtx.dxye.push_back(mu.trk_dxyError());
          MuonsVtx.dz.push_back(mu.trk_dz());
          MuonsVtx.dze.push_back(mu.trk_dzError());
          MuonsVtx.dxysig.push_back(mu.trk_dxy()/mu.trk_dxyError());
          MuonsVtx.dzsig.push_back(mu.trk_dz()/mu.trk_dzError());
          MuonsVtx.selected.push_back(pt>3.0 && fabs(eta)<2.4 && mu.normalizedChi2()<3.0);

          const auto pfIsos0p3Vtx = getPFIsolation(muVec, pfs, 0.3);
          const auto pfIsoChg0p3Vtx = std::get<0>(pfIsos0p3Vtx);
          const auto pfIsoAll0p3Vtx = pfIsoChg0p3Vtx + std::max(0.0, std::get<1>(pfIsos0p3Vtx)+std::get<2>(pfIsos0p3Vtx)-0.5*std::get<3>(pfIsos0p3Vtx));
          MuonsVtx.PFIsoChg0p3.push_back(pfIsoChg0p3Vtx);
          MuonsVtx.PFIsoAll0p3.push_back(pfIsoAll0p3Vtx);
          MuonsVtx.PFRelIsoChg0p3.push_back(pfIsoChg0p3Vtx/pt);
          MuonsVtx.PFRelIsoAll0p3.push_back(pfIsoAll0p3Vtx/pt);
          MuonsVtx.mindrPF0p3.push_back(std::get<4>(pfIsos0p3Vtx));

          const auto pfIsos0p4Vtx = getPFIsolation(muVec, pfs, 0.4);
          const auto pfIsoChg0p4Vtx = std::get<0>(pfIsos0p4Vtx);
          const auto pfIsoAll0p4Vtx = pfIsoChg0p4Vtx + std::max(0.0, std::get<1>(pfIsos0p4Vtx)+std::get<2>(pfIsos0p4Vtx)-0.5*std::get<3>(pfIsos0p4Vtx));
          MuonsVtx.PFIsoChg0p4.push_back(pfIsoChg0p4Vtx);
          MuonsVtx.PFIsoAll0p4.push_back(pfIsoAll0p4Vtx);
          MuonsVtx.PFRelIsoChg0p4.push_back(pfIsoChg0p4Vtx/pt);
          MuonsVtx.PFRelIsoAll0p4.push_back(pfIsoAll0p4Vtx/pt);
          MuonsVtx.mindrPF0p4.push_back(std::get<4>(pfIsos0p4Vtx));

          float mindr=1e6;
          float maxdr=-1;
          for (unsigned int jMu=0; jMu<nMusVtx; ++jMu) {
            if (jMu==iMu) continue;
            auto muOther = musVtx[jMu];
            if (!(fabs(muOther.eta())<2.4)) continue;
            TLorentzVector muVecOther; muVecOther.SetPtEtaPhiM(muOther.pt(), muOther.eta(), muOther.phi(), MUON_MASS);
            float dr = muVec.DeltaR(muVecOther);
            if (dr<mindr) mindr = dr;
            if (dr>maxdr) maxdr = dr;
          }
          MuonsVtx.mindr.push_back(mindr);
          MuonsVtx.maxdr.push_back(maxdr);

          float mindrJet=1e6;
          float mindphiJet=1e6;
          float mindetaJet=1e6;
          for (auto jet : jets) {
            if (!(jet.pt()>20.0 && fabs(jet.eta())<3.0)) continue;
            TLorentzVector jetVec; jetVec.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.m());
            float dr = muVec.DeltaR(jetVec);
            if (dr<mindrJet) {
              mindrJet = dr;
              mindphiJet = fabs(muVec.DeltaPhi(jetVec));
              mindetaJet = fabs(muVec.Eta()-jetVec.Eta());
            }
          }
          MuonsVtx.mindrJet.push_back(mindrJet);
          MuonsVtx.mindphiJet.push_back(mindphiJet);
          MuonsVtx.mindetaJet.push_back(mindetaJet);
          MuonsVtx.vec.push_back(muVec);

          int bestIdx = -1;
          float bestProb = -1;
          if (iMu < vtxIndxVtx.size()) {
            for (auto origVtxIdx : vtxIndxVtx[iMu]) {
              for (unsigned int iSV = 0; iSV < SVsVtx.origIdx.size(); ++iSV) {
                if ((int)SVsVtx.origIdx[iSV] == origVtxIdx && SVsVtx.selected[iSV]) {
                  if (SVsVtx.prob[iSV] > bestProb) {
                    bestProb = SVsVtx.prob[iSV];
                    bestIdx = iSV;
                  }
                }
              }
            }
          }
          MuonsVtx.bestAssocSVIdx.push_back(bestIdx);

          int bestAssocSVOverlapVtxIdx=-1;
          for (unsigned int iSVOverlap=0; iSVOverlap<SVOverlapVtxs.vtxIdxs.size(); iSVOverlap++) {
            for (auto SVOverlapVtxIdx : SVOverlapVtxs.vtxIdxs[iSVOverlap]) {
              if (bestIdx==(int)SVOverlapVtxIdx)
                bestAssocSVOverlapVtxIdx = iSVOverlap;
            }
          }
          MuonsVtx.bestAssocSVOverlapIdx.push_back(bestAssocSVOverlapVtxIdx);

          if (bestIdx >= 0) nMuon_vtx_Assoc++;
          
          if (bestIdx >= 0 && iMu < nhitsbeforesv_vtx.size()) {
            int bestOrigIdx = (int)SVsVtx.origIdx[bestIdx];
            bool foundDV = false;
            for (unsigned int iDV=0; iMu<vtxIndxVtx.size() && iDV<vtxIndxVtx[iMu].size(); ++iDV) {
              if (vtxIndxVtx[iMu].at(iDV) == bestOrigIdx) {
                MuonsVtx.nhitsbeforesv.push_back(nhitsbeforesv_vtx.at(iMu).at(iDV));
                MuonsVtx.ncompatible.push_back(ncompatible_vtx.at(iMu).at(iDV));
                MuonsVtx.ncompatibletotal.push_back(ncompatibletotal_vtx.at(iMu).at(iDV));
                MuonsVtx.nexpectedhits.push_back(nexpectedhits_vtx.at(iMu).at(iDV));
                MuonsVtx.nexpectedhitsmultiple.push_back(nexpectedhitsmultiple_vtx.at(iMu).at(iDV));
                MuonsVtx.nexpectedhitsmultipletotal.push_back(nexpectedhitsmultipletotal_vtx.at(iMu).at(iDV));
                MuonsVtx.nexpectedhitstotal.push_back(nexpectedhitstotal_vtx.at(iMu).at(iDV));
                foundDV = true;
                break;
              }
            }
            if (!foundDV) {
              MuonsVtx.nhitsbeforesv.push_back(-1);
              MuonsVtx.ncompatible.push_back(-1);
              MuonsVtx.ncompatibletotal.push_back(-1);
              MuonsVtx.nexpectedhits.push_back(-1);
              MuonsVtx.nexpectedhitsmultiple.push_back(-1);
              MuonsVtx.nexpectedhitsmultipletotal.push_back(-1);
              MuonsVtx.nexpectedhitstotal.push_back(-1);
            }
          } else {
            MuonsVtx.nhitsbeforesv.push_back(-1);
            MuonsVtx.ncompatible.push_back(-1);
            MuonsVtx.ncompatibletotal.push_back(-1);
            MuonsVtx.nexpectedhits.push_back(-1);
            MuonsVtx.nexpectedhitsmultiple.push_back(-1);
            MuonsVtx.nexpectedhitsmultipletotal.push_back(-1);
            MuonsVtx.nexpectedhitstotal.push_back(-1);
          }

          if (bestIdx >= 0) {
            float svx = SVsVtx.x[bestIdx], svy = SVsVtx.y[bestIdx];
            float dxyCorr = -(svx - PV_x)*TMath::Sin(phi) + (svy - PV_y)*TMath::Cos(phi);
            MuonsVtx.dxyCorr.push_back(dxyCorr);
            TLorentzVector muVecForCorr; muVecForCorr.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
            float phiCorr = getCorrectedPhi(mu, muVecForCorr, dxyCorr, svx, svy);
            MuonsVtx.phiCorr.push_back(phiCorr);
            MuonsVtx.vec.back().SetPtEtaPhiM(pt, eta, phiCorr, MUON_MASS);
          } 
          else {
            MuonsVtx.dxyCorr.push_back(0.);
            MuonsVtx.phiCorr.push_back(phi);
          }
        }
      }

      nPreMu++;
      if (Muons.pt.size() < 2 && MuonsVtx.pt.size() < 2)
        continue;

      if (requireAtLeastOneSV && SVs.x.empty() && SVsVtx.x.empty())
        continue;
      nSV++;

      Muons.sort();
      MuonsVtx.sort();

      tout->Fill();
      nSaved++;
    }
    } catch (const std::exception& e) {
      std::cout << "WARNING: read error in " << inputFile << " at event " << iEv
                << " (" << e.what() << ") -- skipping rest of file, keeping "
                << nSaved << " saved" << std::endl;
    }
    iFile++;
    std::cout << "Events saved: " << nSaved <<  "\n";
    std::cout << "Events good: " << nGoodRun <<  "\n";
    std::cout << "Events no-duplicate: " << nDuplicate <<  "\n";
    std::cout << "Events fraction: " << nFraction <<  "\n";
    std::cout << "Events pass L1: " << nL1 <<  "\n";
    std::cout << "Events pass HLT: " << nHLT <<  "\n";
    std::cout << "Events pre-mu: " << nPreMu <<  "\n";
    std::cout << "Events with >=1 SV: " << nSV <<  "\n";
    std::cout << "Events saved: " << nSaved <<  "\n";
    std::cout<<"\n\n";

    cutflow->SetBinContent(1, counts->GetBinContent(1));
    cutflow->SetBinContent(2, cutflow->GetBinContent(2) + nSaved);
    cutflow->SetBinContent(3, cutflow->GetBinContent(3) + nGoodRun);
    cutflow->SetBinContent(4, cutflow->GetBinContent(4) + nDuplicate);
    cutflow->SetBinContent(5, cutflow->GetBinContent(5) + nFraction);
    cutflow->SetBinContent(6, cutflow->GetBinContent(6) + nL1);
    cutflow->SetBinContent(7, cutflow->GetBinContent(7) + nHLT);
    cutflow->SetBinContent(8, cutflow->GetBinContent(8) + nPreMu);
    cutflow->SetBinContent(9, cutflow->GetBinContent(9) + nSV);

  }

  bar.finish();
  fout->Write();
  fout->Close();

  return;
}
