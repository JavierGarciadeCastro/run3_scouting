#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"

class ScoutingMuonL3PreFilterEmulator : public edm::stream::EDProducer<> {
public:
  explicit ScoutingMuonL3PreFilterEmulator(const edm::ParameterSet&);
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonToken_;
  const int minN_;
  const double maxEta_, minPt_;
  const std::string outputName_;
};

ScoutingMuonL3PreFilterEmulator::ScoutingMuonL3PreFilterEmulator(const edm::ParameterSet& p)
    : muonToken_(consumes<std::vector<Run3ScoutingMuon>>(p.getParameter<edm::InputTag>("muons"))),
      minN_(p.getParameter<int>("MinN")),
      maxEta_(p.getParameter<double>("MaxEta")),
      minPt_(p.getParameter<double>("MinPt")),
      outputName_(p.getParameter<std::string>("outputName")) {
  produces<bool>(outputName_).setBranchAlias(outputName_);
}

void ScoutingMuonL3PreFilterEmulator::produce(edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<std::vector<Run3ScoutingMuon>> muons;
  iEvent.getByToken(muonToken_, muons);

  int n = 0;
  if (muons.isValid()) {
    for (const auto& mu : *muons) {
      if (std::abs(mu.eta()) > maxEta_)
        continue;
      if (mu.pt() < minPt_)
        continue;
      ++n;
    }
  }

  iEvent.put(std::make_unique<bool>(n >= minN_), outputName_);
}

void ScoutingMuonL3PreFilterEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription d;
  d.add<edm::InputTag>("muons", edm::InputTag("hltScoutingMuonPackerNoVtx"));
  // defaults = hltDoubleMuonL3PreFilteredScoutingNoVtx
  d.add<int>("MinN", 2);
  d.add<double>("MaxEta", 2.5);
  d.add<double>("MinPt", 0.0);  
  d.add<std::string>("outputName", "DoubleMuonNoVtxv4");  // no underscores allowed
  descriptions.add("scoutingMuonL3PreFilterEmulator", d);
}

DEFINE_FWK_MODULE(ScoutingMuonL3PreFilterEmulator);
