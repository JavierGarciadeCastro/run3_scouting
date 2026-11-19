#include "Scouting/NtupleMaker/plugins/TriggerMaker.h" 
#include "TMath.h"

using namespace edm;
using namespace std;

TriggerMaker::TriggerMaker(const edm::ParameterSet& iConfig) :
  doL1_(iConfig.getParameter<bool>("doL1")),
  doTriggerObjects_(iConfig.getParameter<bool>("doTriggerObjects")),
  isMiniAOD_(iConfig.getParameter<bool>("isMiniAOD")),
  triggerCache_(triggerExpression::Data(iConfig.getParameterSet("triggerConfiguration"), consumesCollector())),
  vtriggerAlias_(iConfig.getParameter<vector<string>>("triggerAlias")),
  vtriggerSelection_(iConfig.getParameter<vector<string>>("triggerSelection")),
  l1GtUtils_(nullptr)
{
  vtriggerSelector_.reserve(vtriggerSelection_.size());
  for (auto const& vt:vtriggerSelection_) vtriggerSelector_.push_back(triggerExpression::parse(vt));

  if (doL1_){
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<InputTag>("AlgInputTag"));
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1GtUtils_ = std::make_shared<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), l1t::UseEventSetupIn::RunAndEvent);
  }

  if (doTriggerObjects_) {
    edm::InputTag trigObjTag = isMiniAOD_
        ? edm::InputTag("slimmedPatTrigger", "", "PAT")
        : edm::InputTag("patTrigger");
    triggerObjectsToken = consumes<pat::TriggerObjectStandAloneCollection>(trigObjTag);
    triggerResultsToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","", "HLT"));

    produces<std::vector<std::string> >("trigObjsfilters").setBranchAlias("trigObjs_filters");
    produces<std::vector<int> >("trigObjsid").setBranchAlias("trigObjs_id");
    produces<std::vector<float> >("trigObjspt").setBranchAlias("trigObjs_pt");
    produces<std::vector<float> >("trigObjseta").setBranchAlias("trigObjs_eta");
    produces<std::vector<float> >("trigObjsphi").setBranchAlias("trigObjs_phi");
    produces<std::vector<float> >("trigObjsmass").setBranchAlias("trigObjs_mass");
    produces<std::vector<int> >("trigObjspassLast").setBranchAlias("trigObjs_passLast");
  }

  produces<std::vector<std::string> >("l1name").setBranchAlias("l1_name");
  produces<std::vector<bool> >("l1result").setBranchAlias("l1_result");
  produces<std::vector<double> >("l1prescale").setBranchAlias("l1_prescale");

  produces<std::vector<std::string> >("hltname").setBranchAlias("hlt_name");
  produces<std::vector<bool> >("hltresult").setBranchAlias("hlt_result");
}

TriggerMaker::~TriggerMaker(){
  for (auto& vts:vtriggerSelector_) delete vts;
}

void TriggerMaker::beginJob(){}

void TriggerMaker::endJob(){}

void TriggerMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}

void TriggerMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  unique_ptr<vector<string> > l1_name(new vector<string>);
  unique_ptr<vector<bool> > l1_result(new vector<bool>);
  unique_ptr<vector<double> > l1_prescale(new vector<double>);

  unique_ptr<vector<string> > hlt_name(new vector<string>);
  unique_ptr<vector<bool> > hlt_result(new vector<bool>);

  if (triggerCache_.setEvent(iEvent, iSetup)){
    auto trigAlias=vtriggerAlias_.cbegin();
    for (auto& vts:vtriggerSelector_){
      bool result = false;
      if (vts){
        if (triggerCache_.configurationUpdated()) vts->init(triggerCache_);
        result = (*vts)(triggerCache_);
      }
      hlt_result->push_back(result);
      hlt_name->push_back(*trigAlias);
      trigAlias++;
    }
  }

  if (doL1_){
    // L1TGlobalUtil silently returns all-false decisions if the AlgBlk is missing: fail loudly instead
    edm::Handle<BXVector<GlobalAlgBlk>> algBlkH;
    iEvent.getByToken(algToken_, algBlkH);
    if (!algBlkH.isValid())
      throw cms::Exception("TriggerMaker::produce: GlobalAlgBlk product not found (check AlgInputTag / process name)");
    l1GtUtils_->retrieveL1(iEvent, iSetup, algToken_);


    for (auto const& l1seed:l1Seeds_){
      bool l1htbit = 0;
      double prescale = -1;
      l1GtUtils_->getFinalDecisionByName(l1seed, l1htbit);
      l1GtUtils_->getPrescaleByName(l1seed, prescale);
      l1_result->push_back(l1htbit);
      l1_name->push_back(l1seed);
      l1_prescale->push_back(prescale);
    }
  }

  unique_ptr<std::vector<int> > trigObjsid(new std::vector<int>);
  unique_ptr<std::vector<float > > trigObjspt(new std::vector<float>);
  unique_ptr<std::vector<float > > trigObjseta(new std::vector<float>);
  unique_ptr<std::vector<float > > trigObjsphi(new std::vector<float>);
  unique_ptr<std::vector<float > > trigObjsmass(new std::vector<float>);
  unique_ptr<std::vector<int > > trigObjspassLast(new std::vector<int>);
  unique_ptr<std::vector<std::string> > trigObjsfilters(new std::vector<std::string>);

  if (doTriggerObjects_) {

      iEvent.getByToken(triggerResultsToken, triggerResultsH_);
      if (!triggerResultsH_.isValid()) throw cms::Exception("TriggerMaker::produce: error getting TriggerResults product from Event!");

      iEvent.getByToken(triggerObjectsToken, triggerObjectStandAlonesH_);
      if (!triggerObjectStandAlonesH_.isValid()) {
          throw cms::Exception("TriggerMaker::produce: error getting TriggerObjectsStandAlone product from Event!");
      }

      triggerNames_ = iEvent.triggerNames(*triggerResultsH_);
      unsigned int nTriggers = triggerResultsH_->size();
      // find the scouting trigger index
      int itrig = -1;
      for(unsigned int i = 0; i < nTriggers; ++i){
          const string& name = triggerNames_.triggerName(i);
          //should check this part further
          if (name.find("DST_Run3_PFScoutingPixelTracking") != std::string::npos || name.find("DST_Run3_DoubleMu3_PFScoutingPixelTracking") != std::string::npos || name.find("DST_PFScouting_DoubleMuon") != std::string::npos) {
              itrig = i;
              break;
          }
      }
      if (itrig < 0) {
          throw cms::Exception("TriggerMaker::produce: Couldn't find DST_Run3_PFScoutingPixelTracking or DST_Run3_DoubleMu3_PFScoutingPixelTracking or DST_PFScouting_DoubleMuon in the list of triggers");
      }

      pat::TriggerObjectStandAlone TO;
      for ( uint iobj = 0; iobj < triggerObjectStandAlonesH_->size(); iobj++ ) {

          TO = triggerObjectStandAlonesH_->at(iobj);
          TO.unpackPathNames( triggerNames_ );
          TO.unpackFilterLabels( iEvent,*triggerResultsH_ );

              const string& name = triggerNames_.triggerName(itrig);

              // save all objects associated to path, regardless of final result. And save filter names 
              if ( TO.hasPathName(name, false, false ) ) {
                  int storeID = 0;
                  std::vector<int> IDs = TO.filterIds();
                  if (IDs.size() == 1) storeID = IDs[0];
                  else if (IDs.size() > 1) {
                      // Making some arbitrary choices
                      if ( IDs[0]==85 || IDs[1]==85 ) storeID = 85; // when in doubt call it jet (and not the bjet, 86)
                      if ( IDs[0]==92 || IDs[1]==92 ) storeID = 92; // when in doubt call it cluster (and not Photon, 81, or Electron, 82)
                  }
                  bool saveFilters = false;
                  std::string filterslist = "";
                  if ( IDs.size() > 0 ) {
                      int id = abs(IDs[0]);
                      // From: TriggerTypeDefs.h
                      //       TriggerL1Mu           = -81,
                      //       TriggerL1NoIsoEG      = -82,
                      //       TriggerL1IsoEG        = -83,
                      //       TriggerPhoton         = +81,
                      //       TriggerElectron       = +82,
                      //       TriggerMuon           = +83,
                      //       TriggerCluster        = +92,
                      if ( id == 81 || id == 82 || id == 83 || IDs[0] == 92) saveFilters = true;
                      if ( IDs.size() > 1 ) {
                          int id = abs(IDs[1]);
                          if ( id == 81 || id == 82 || id == 83 || IDs[1] == 92) saveFilters = true;
                      }
                  }
                  if (saveFilters) {
                      std::vector< std::string > filter_labels = TO.filterLabels();
                      for (uint j  = 0; j < filter_labels.size(); j++) {
                          filterslist += filter_labels[j];
                          filterslist += " ";
                      }
                  }

                  trigObjsid->push_back(storeID);
                  trigObjspt->push_back(TO.p4().pt());
                  trigObjseta->push_back(TO.p4().eta());
                  trigObjsphi->push_back(TO.p4().phi());
                  trigObjsmass->push_back(TO.p4().mass());
                  trigObjspassLast->push_back(TO.hasPathName(name, true));
                  trigObjsfilters->push_back(filterslist);
              } // hasPathName

      } // end of loop over trigger objects

  } // end if(doTriggerObjects_)

  iEvent.put(std::move(l1_name), "l1name");
  iEvent.put(std::move(l1_result), "l1result");
  iEvent.put(std::move(l1_prescale), "l1prescale");
  iEvent.put(std::move(hlt_name), "hltname");
  iEvent.put(std::move(hlt_result), "hltresult");
   
  if (doTriggerObjects_) {
      iEvent.put(std::move(trigObjsfilters),"trigObjsfilters");
      iEvent.put(std::move(trigObjsid),"trigObjsid");
      iEvent.put(std::move(trigObjspt),"trigObjspt");
      iEvent.put(std::move(trigObjseta),"trigObjseta");
      iEvent.put(std::move(trigObjsphi),"trigObjsphi");
      iEvent.put(std::move(trigObjsmass),"trigObjsmass");
      iEvent.put(std::move(trigObjspassLast),"trigObjspassLast");
  }

}

DEFINE_FWK_MODULE(TriggerMaker);
