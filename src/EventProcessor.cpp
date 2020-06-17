//  EventProcessor.cpp
//
//  Created by Jeremi Niedziela on 30/01/2019.

#include "EventProcessor.hpp"
#include "Logger.hpp"

EventProcessor eventProcessor = EventProcessor();

EventProcessor::EventProcessor()
{
  singleNamesUint = {
    "run",
    "lumi",
    "wasTagged"
  };
  
  singleNamesUlongLong = {
    "evt"
  };
  
  singleNamesFloat = {
    "vertex_x",
    "vertex_y",
    "vertex_z",
    "xsec",
    "wgtsum",
    "genWeight",
    "met_sumEt",
    "met_pt",
    "met_genPt",
    "met_mass",
    "met_phi",
    "met_eta",
    "metNoMu_pt",
    "met_jecUp_pt",
    "met_jecDown_pt",
    "metNoMu_mass",
    "metNoMu_phi",
    "metNoMu_eta",
  };
  
  arrayNamesFloat = {
    "GenChargino_pt",
    "GenChargino_eta",
    "GenChargino_phi",
    "GenChargino_mass",
    "GenChargino_charge",
  };
  
  singleNamesInt = {
    "nVert",
    "nJet30",
    "nJet30a",
    "nTauGood",
    "nGenChargino",
    "HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
    "HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
    "Flag_goodVertices",
    "Flag_BadPFMuonFilter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_eeBadScFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_ecalBadCalibFilter",
    "Flag_globalTightHalo2016Filter"
  };
  
  arrayNamesFriendInt = {
    "pion_charge",
    "chargino_charge",
    "chargino_nTrackerLayers",
    "pion_simHits_subDet",
    "chargino_simHits_subDet",
    "pixelCluster_subDet",
    "stripCluster_subDet",
    "pionCluster_subDet",
    "generalTrack_nLoops",
    "generalTrack_isLooper",
    "generalTrack_nPionHits",
  };
  
  arrayNamesFriendFloat = {
    "pion_vx",
    "pion_vy",
    "pion_vz",
    "pion_px",
    "pion_py",
    "pion_pz",
    "chargino_eta",
    "chargino_phi",
    "chargino_pt",
    "chargino_mass",
    "pion_simHits_x",
    "pion_simHits_y",
    "pion_simHits_z",
    "pion_simHits_t",
    "chargino_simHits_x",
    "chargino_simHits_y",
    "chargino_simHits_z",
    "pixelCluster_x",
    "pixelCluster_y",
    "pixelCluster_z",
    "stripCluster_x",
    "stripCluster_y",
    "stripCluster_z",
    "stripCluster_ex",
    "stripCluster_ey",
    "stripCluster_ez",
    "pionCluster_x",
    "pionCluster_y",
    "pionCluster_z",
    "pionCluster_ex",
    "pionCluster_ey",
    "pionCluster_ez",
    "pixelCluster_charge",
    "stripCluster_charge",
    "pionCluster_charge",
    "generalTrack_px",
    "generalTrack_py",
    "generalTrack_pz",
    "generalTrack_d0",
    "generalTrack_charge",
    "generalTrack_chi2",
    "generalTrack_eta",
    "generalTrack_phi",
    "generalTrack_nHits",
    "generalTrack_nMissingHits",
  };
}

EventProcessor::~EventProcessor()
{
  
}

void EventProcessor::ApplyTrackCut(shared_ptr<Event> event, const TrackCut &cut, vector<int> *trackCutReasons)
{
  
  
  for(auto track = event->tracks.begin(); track != event->tracks.end();){
    bool passesGenCuts = true;
    
    if(cut.requiresCorrectNlayers || cut.requiresCorrectCharge){
      
      // Find gen track matching this rec track
      double minDeltaR = inf;
      Track matchingGenTrack;
      
      for(auto &genTrack : event->genCharginoTrack){
        double deltaR = sqrt(pow(genTrack.GetEta() - (*track)->GetEta(), 2) +
                             pow(genTrack.GetPhi() - (*track)->GetPhi(), 2));
        
        if(deltaR < minDeltaR){
          minDeltaR = deltaR;
          matchingGenTrack = genTrack;
        }
      }
      
      // check conditions
      if(cut.requiresCorrectNlayers){
        if(matchingGenTrack.GetNtrackerLayers() != (*track)->GetNtrackerLayers()) passesGenCuts = false;
      }
      if(cut.requiresCorrectCharge){
        if(matchingGenTrack.GetCharge() != (*track)->GetCharge()) passesGenCuts = false;
      }
    }
    
    if(!trackProcessor.IsPassingCut(*track,cut, trackCutReasons) || !passesGenCuts) track = event->tracks.erase(track);
    else                                                           track++;
  }
}

void EventProcessor::ApplyJetCut(shared_ptr<Event> event, const JetCut &cut)
{
  for(auto jet = event->jets.begin(); jet != event->jets.end();){
    if(!jetProcessor.IsPassingCut(*jet,cut)) jet = event->jets.erase(jet);
    else                                     jet++;
  }
}

void EventProcessor::ApplyLeptonCut(shared_ptr<Event> event, const LeptonCut &cut)
{
  for(auto lepton = event->leptons.begin(); lepton != event->leptons.end();){
    if(!leptonProcessor.IsPassingCut(*lepton, cut))  lepton = event->leptons.erase(lepton);
    else                                             lepton++;
  }
}


bool EventProcessor::IsPassingCut(const shared_ptr<Event> event, const EventCut &cut, vector<int> *cutReasons)
{
  int cutThroughIter=0;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 0
  
  // check the trigger
  if(cut.requiresMetNoMuTrigger && !event->metNoMuTrigger) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 1
  
  // check MET filters
  if(cut.requiresPassingAllFilters){
    if(   !event->flag_goodVertices
       || !event->flag_globalTightHalo2016
       || !event->flag_HBHEnoise
       || !event->flag_HBHEnoiseIso
       || !event->flag_EcalDeadCell
       || !event->flag_badPFmuon
//     || flag_badChargedCandidate // not recommended
//     || !event->flag_eeBadSc     // not suggested
       || !event->flag_ecalBadCalib){
      return false;
    }
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 2
  
  // check number of objects
  if(cut.nGenPions.IsOutside((uint)event->genPionHelices.size())) return false;
  
  if(cut.nGenCharginos.IsOutside(event->nGenChargino)) return false;
  if(cut.nLeptons.IsOutside(event->GetNleptons()))                return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 3
  
  if(cut.nTaus.IsOutside(event->nTau))  return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 4
  
  vector<shared_ptr<Lepton>> muons;
  for(auto l : event->leptons){
    if(abs(l->GetPid()) == 13) muons.push_back(l);
  }
  
  // check number of muons
  if(cut.nMuons.IsOutside((int)muons.size())) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 5
  
  // make sure they have an opposite sign
  if(cut.requiresTwoOpositeMuons){
    if(muons.size() != 2) return false;
    if(muons[0]->GetPid() != -muons[1]->GetPid()) return false;
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 6
  
  // apply tight muon cuts (tightID flag, pt > 20 GeV, isolation < 0.15
  if(cut.requiresTightMuon){
    auto tightMuonCut = LeptonCut();
    tightMuonCut.SetRelativeIsolation(range<double>(-inf,0.15));
    tightMuonCut.SetRequireTightID(true);
    tightMuonCut.SetPt(range<double>(20.0,inf));
    
    bool atLeastOneTightMuon = false;
    for(auto muon : muons){
      if(leptonProcessor.IsPassingCut(muon, tightMuonCut)) atLeastOneTightMuon = true;
    }
    if(!atLeastOneTightMuon) return false;
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 7
  
  // check that invariant mass of muons is close to Z mass
  if(cut.requiresMuonsFromZ){
    TLorentzVector muon1vector, muon2vector, muonVectorSum;
    
    if(muons.size() != 2){
      cout<<"ERROR -- requested muons to come from Z decay, but there is "<<muons.size()<<" muons in the event!!"<<endl;
      cout<<"This event will be discarded!! Maybe you should require exactly two muons in the event?"<<endl;
      if(cutReasons) cutReasons->at(cutThroughIter++)++; // 7
      return false;
    }
    
    muon1vector.SetPtEtaPhiM(muons[0]->GetPt(),muons[0]->GetEta(),muons[0]->GetPhi(), 0.1057);
    muon2vector.SetPtEtaPhiM(muons[1]->GetPt(),muons[1]->GetEta(),muons[1]->GetPhi(), 0.1057);
    
    muonVectorSum += muon1vector;
    muonVectorSum += muon2vector;
    
    if(muonVectorSum.M() < 60. || muonVectorSum.M() > 120.) return false;
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 8
  
  if(cut.metPt.IsOutside(event->metPt)) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 9
  
  if(cut.metNoMuPt.IsOutside(event->metNoMuPt)) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 10
  
  // Check if jets are not too close to muons
  if(cut.jetMuonDeltaPhi.GetMin() > 0.0){
    vector<TLorentzVector> muonVectors;
    double muonMass = 0.1057;
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), muonMass);
      muonVectors.push_back(mv);
    }
    
    for(int iJet=0;iJet<event->GetNcentralJets();iJet++){
      shared_ptr<Jet> j = event->GetJet(iJet);
      if(!j->IsForward()){
        
        TLorentzVector jetVector;
        jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
        
        for(auto muonVector : muonVectors){
          if(cut.jetMuonDeltaPhi.IsOutside(fabs(jetVector.DeltaR(muonVector)))) return false;
        }
      }
    }
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 11
  
  // Check if iso tracks are not too close to muons
  if(cut.trackMuonDeltaPhi.GetMin() > 0.0){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      double muonMass = 0.1057;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), muonMass);
      muonVectors.push_back(mv);
    }
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      shared_ptr<Track> t = event->tracks[iTrack];
      
      TLorentzVector trackVector;
      trackVector.SetPtEtaPhiM(t->GetPt(), t->GetEta(), t->GetPhi(), t->GetMass());
      
      for(auto muonVector : muonVectors){
        if(cut.trackMuonDeltaPhi.IsOutside(fabs(trackVector.DeltaR(muonVector)))) return false;
      }
    }
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 12
  
  // check number of tracks and jets after removing those that are too close to muons
  if(cut.nJets.IsOutside(event->GetNcentralJets())) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 13
  
  if(cut.nTracks.IsOutside(event->GetNtracks())) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 14
  
  // find the jet with the highest pt meeting leading jet criteria
  shared_ptr<Jet> leadingJet = nullptr;
  double highestPt = -inf;
  
  for(auto jet : event->jets){
      if(jet->GetPt() > highestPt){
        highestPt = jet->GetPt();
        leadingJet = jet;
      }
  }
  
  if(!leadingJet) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 15
  
  if(cut.leadingJetPt.IsOutside(leadingJet->GetPt())        ||
     cut.leadingJetEta.IsOutside(leadingJet->GetEta())      ||
     cut.leadingJetChHEF.IsOutside(leadingJet->GetChHEF())  ||
     cut.leadingJetNeHEF.IsOutside(leadingJet->GetNeHEF())){
    leadingJet = nullptr;
  }
  
  // check if there is a leading jet in the event
  if(!leadingJet) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 16

  if(cut.jetMetDeltaPhi.GetMin() > 0.0){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(event->metNoMuPt, event->metNoMuEta, event->metNoMuPhi, event->metNoMuMass);
    
    for(auto j : event->jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(cut.jetMetDeltaPhi.IsOutside(fabs(metVector.DeltaPhi(jetVector)) )) return false;
    }
  }
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 17
  
  bool pionPassed = false;
  
  if(cut.genPionsPt.GetMin() < -99999) pionPassed = true;
  else{
    for(auto &pion : event->GetGenPionHelices()){
      if(cut.genPionsPt.IsInside(pion.GetMomentum().GetTransverse())){
        pionPassed=true;
        break;
      }
    }
  }
  if(!pionPassed) return false;
  if(cutReasons) cutReasons->at(cutThroughIter++)++; // 18
  
  survivingEvents.push_back(event);
  return true;
}

shared_ptr<Event> EventProcessor::GetEventFromTree(xtracks::EDataType dataType, int setIter, int year, TTree *friendTree, TTree *prefireTree, TH1D *metWeights)
{
  for(auto &name_val : singleValuesInt ){
    if(name_val.second < -999999){
      cout<<"ERROR -- branch "<<name_val.first<<" was not read correctly!"<<endl;
    }
  }
  
  auto event = make_shared<Event>();
  
  double lumi = config.params["total_luminosity_"+to_string(year)] * 1000.; // transform from fb^-1 to pb^-1
  double weight = lumi * singleValuesFloat["genWeight"] / singleValuesFloat["wgtsum"];
  
  //    static map<string,set<double>> wgts;
  //
  //    if(*_genWgt != 1.0 && wgts[fileName].find(*_genWgt) == wgts[fileName].end() ){
  //      wgts[fileName].insert(*_genWgt);
  //      cout<<*_genWgt<<"\t"<<fileName<<endl;
  //    }
  
  if(dataType == xtracks::kBackground){
    weight *= singleValuesFloat["xsec"];
  }
  if(dataType == xtracks::kSignal){
    // it's not clear how to calculate weights for the signal...
    
    // cross section for given signal (stored in fb, here transformed to pb to match background units
    weight *= 0.001 * (signalCrossSectionOneTrack.at((ESignal)setIter) +
                       signalCrossSectionTwoTracks.at((ESignal)setIter));
    //      if(*_nGenChargino == 1){
    //        weight *= 0.001 * signalCrossSectionOneTrack[iSig]; // cross section for given signal (stored in fb, here transformed to pb to match background units
    //      }
    //      else if(*_nGenChargino == 2){
    //        weight *= 0.001 * signalCrossSectionTwoTracks[iSig];
    //      }
    //      else{
    //        cout<<"WARNING -- number of generator-level charginos different than 1 or 2"<<endl;
    //      }
  }
  else if(dataType == xtracks::kData){
    weight = 1;
  }
  
  if(prefireTree){
    // uberhack to select needed entry directly...
    prefireTree->Draw("Entry$>>hist(Entries$,0,Entries$)",
                      ("lumi=="+to_string(singleValuesUint["lumi"])+
                       "&&run=="+to_string(singleValuesUint["run"])+
                       "&&evt=="+to_string(singleValuesUlonglong["evt"])).c_str(),
                      "goff");
    
    TH1I *hist = (TH1I*)gDirectory->Get("hist");
    Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
    if(hist) delete hist;
    
    if(iEntry==0){
      Log(0)<<"Event with run: "<<singleValuesUint["run"]<<"\tlumi: "<<singleValuesUint["lumi"]<<"\tevNumber: "<<to_string(singleValuesUlonglong["evt"])<<" not found in prefire tree!!\n";
    }
    else{
      prefireTree->GetEntry(iEntry);
      weight *= prefireWeight;
    }
  }
  if(metWeights){
    double metNoMu = singleValuesFloat["metNoMu_pt"];
    double scale = metWeights->GetBinContent(metWeights->GetXaxis()->FindFixBin(metNoMu));
    weight *= scale;
  }
  
  event->weight   = weight;
  event->dataType = dataType;
  event->setIter  = setIter;
  
  event->lumiSection = singleValuesUint["lumi"];
  event->runNumber   = singleValuesUint["run"];
  event->wasTagged   = singleValuesUint["wasTagged"];
  event->eventNumber = singleValuesUlonglong["evt"];
  
  event->nVertices                = singleValuesInt["nVert"];
  event->nJet30                   = singleValuesInt["nJet30"];
  event->nJet30a                  = singleValuesInt["nJet30a"];
  event->nTau                     = singleValuesInt["nTauGood"];
  event->nGenChargino = nGenChargino = singleValuesInt["nGenChargino"];
  
  for(uint i=0;i<event->nGenChargino;i++){
    event->genCharginos.push_back(Track(arrayValuesFloat["GenChargino_eta"][i],
                                        arrayValuesFloat["GenChargino_phi"][i],
                                        arrayValuesFloat["GenChargino_charge"][i],
                                        -1,
                                        arrayValuesFloat["GenChargino_pt"][i],
                                        arrayValuesFloat["GenChargino_mass"][i]));
  }
  
  
  event->metNoMuTrigger           = singleValuesInt["HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"]
                                 || singleValuesInt["HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"];
  event->flag_goodVertices        = singleValuesInt["Flag_goodVertices"];
  event->flag_badPFmuon           = singleValuesInt["Flag_BadPFMuonFilter"];
  event->flag_HBHEnoise           = singleValuesInt["Flag_HBHENoiseFilter"];
  event->flag_HBHEnoiseIso        = singleValuesInt["Flag_HBHENoiseIsoFilter"];
  event->flag_EcalDeadCell        = singleValuesInt["Flag_EcalDeadCellTriggerPrimitiveFilter"];
  event->flag_eeBadSc             = singleValuesInt["Flag_eeBadScFilter"];
  event->flag_badChargedCandidate = singleValuesInt["Flag_BadChargedCandidateFilter"];
  event->flag_ecalBadCalib        = singleValuesInt["Flag_ecalBadCalibFilter"];
  event->flag_globalTightHalo2016 = singleValuesInt["Flag_globalTightHalo2016Filter"];
  
  event->vertex      = make_unique<Point>(singleValuesFloat["vertex_x"],
                                          singleValuesFloat["vertex_y"],
                                          singleValuesFloat["vertex_z"]);
  
  event->xsec        = singleValuesFloat["xsec"];
  event->wgtsum      = singleValuesFloat["wgtsum"];
  event->genWeight   = singleValuesFloat["genWeight"];
  event->metSumEt    = singleValuesFloat["met_sumEt"];
  event->metPt       = singleValuesFloat["met_pt"];
  event->metGenPt    = singleValuesFloat["met_genPt"];
  event->metMass     = singleValuesFloat["met_mass"];
  event->metPhi      = singleValuesFloat["met_phi"];
  event->metEta      = singleValuesFloat["met_eta"];
  event->metNoMuPt          = singleValuesFloat["metNoMu_pt"];
  event->metNoMuPtJecUp     = singleValuesFloat["met_jecUp_pt"];
  event->metNoMuPtJecDown   = singleValuesFloat["met_jecDown_pt"];
  event->metNoMuMass = singleValuesFloat["metNoMu_mass"];
  event->metNoMuPhi  = singleValuesFloat["metNoMu_phi"];
  event->metNoMuEta  = singleValuesFloat["metNoMu_eta"];
  
  // Load additional info from friend tree
  if(!friendTree) return event;
  
  // uberhack to select needed entry directly...
  friendTree->Draw("Entry$>>hist(Entries$,0,Entries$)",
             ("lumiBlock=="+to_string(singleValuesUint["lumi"])+
              "&&runNumber=="+to_string(singleValuesUint["run"])+
              "&&eventNumber=="+to_string(singleValuesUlonglong["evt"])).c_str(),
             "goff");
  
  TH1I *hist = (TH1I*)gDirectory->Get("hist");
  Long64_t iEntry = hist->GetBinLowEdge(hist->FindFirstBinAbove(0));
  if(hist) delete hist;
  
  if(iEntry==0){
    Log(0)<<"No event with run: "<<singleValuesUint["run"]<<"\tlumi: "<<singleValuesUint["lumi"]<<"\tevNumber: "<<to_string(singleValuesUlonglong["evt"])<<" was found\n";
    return event;
  }
  
  friendTree->GetEntry(iEntry);
  event->hasFriendData = true;
  // end of uberhack
  
  map<int, string> subDetMap = {
    {0,  "PixelBarrel"},
    {1,  "PixelEndcap"},
    {2,  "TIB"},
    {3,  "TOB"},
    {4,  "TID"},
    {5,  "TEC"},
    {6,  "CSC"},
    {7,  "DT"},
    {8,  "RPCBarrel"},
    {9,  "RPCEndcap"},
    {10, "GEM"},
    {11, "ME0"},
    {12, "P2OTB"},
    {13, "P2OTEC"},
    {14, "P1PXB"},
    {15, "P1PXEC"},
    {16, "P2PXB"},
    {17, "P2PXEC"},
    {18, "TimingBarrel"},
    {19, "TimingEndcap"},
    {20, "invalidDet"}
  };
  
  if(friendTree->GetBranchStatus("generalTrack_px")){
    for(uint i=0;i<arrayValuesFriendFloat["generalTrack_px"]->size();i++){
      Point origin(0,0,0);
      Point momentum(arrayValuesFriendFloat["generalTrack_px"]->at(i),
                     arrayValuesFriendFloat["generalTrack_py"]->at(i),
                     arrayValuesFriendFloat["generalTrack_pz"]->at(i));
      
      int charge = arrayValuesFriendFloat["generalTrack_charge"]->at(i);
      
      
      Helix generalTrack(origin, momentum, charge);
      
      generalTrack.nLoops = arrayValuesFriendInt["generalTrack_nLoops"]->at(i);
      generalTrack.isLooper = arrayValuesFriendInt["generalTrack_isLooper"]->at(i);
      generalTrack.nRecPionHits = arrayValuesFriendInt["generalTrack_nPionHits"]->at(i);
      generalTrack.d0 = arrayValuesFriendFloat["generalTrack_d0"]->at(i);
      generalTrack.chi2 = arrayValuesFriendFloat["generalTrack_chi2"]->at(i);
      generalTrack.eta = arrayValuesFriendFloat["generalTrack_eta"]->at(i);
      generalTrack.phi = arrayValuesFriendFloat["generalTrack_phi"]->at(i);
      generalTrack.nRecHits = arrayValuesFriendFloat["generalTrack_nHits"]->at(i);
      generalTrack.nMissingHits = arrayValuesFriendFloat["generalTrack_nMissingHits"]->at(i);
      
      event->AddGeneralTrack(generalTrack);
    }
  }
  for(uint i=0;i<arrayValuesFriendFloat["pion_vx"]->size();i++){
    // change units from cm to mm and from GeV to MeV
    event->genPionHelices.emplace(event->genPionHelices.end(),
                                  Helix(Point(10*arrayValuesFriendFloat["pion_vx"]->at(i),
                                              10*arrayValuesFriendFloat["pion_vy"]->at(i),
                                              10*arrayValuesFriendFloat["pion_vz"]->at(i)),
                                        Point(1000*arrayValuesFriendFloat["pion_px"]->at(i),
                                              1000*arrayValuesFriendFloat["pion_py"]->at(i),
                                              1000*arrayValuesFriendFloat["pion_pz"]->at(i)),
                                        1*arrayValuesFriendInt["pion_charge"]->at(i)));
  }
  
  if(config.params["load_hits"]){
    for(uint i=0;i<arrayValuesFriendFloat["pion_simHits_x"]->size();i++){
      // convert cm to mm
      event->pionSimHits.push_back(make_shared<Point>(10*arrayValuesFriendFloat["pion_simHits_x"]->at(i),
                                                      10*arrayValuesFriendFloat["pion_simHits_y"]->at(i),
                                                      10*arrayValuesFriendFloat["pion_simHits_z"]->at(i),
                                                      0,
                                                      subDetMap[arrayValuesFriendInt["pion_simHits_subDet"]->at(i)],
                                                      0,0,0,0,-1,
                                                      arrayValuesFriendFloat["pion_simHits_t"]->at(i)));
    }
    
    for(uint i=0;i<arrayValuesFriendFloat["chargino_simHits_x"]->size();i++){
      // convert cm to mm
      event->charginoSimHits.push_back(make_shared<Point>(10*arrayValuesFriendFloat["chargino_simHits_x"]->at(i),
                                                          10*arrayValuesFriendFloat["chargino_simHits_y"]->at(i),
                                                          10*arrayValuesFriendFloat["chargino_simHits_z"]->at(i),
                                                          0,
                                                          subDetMap[arrayValuesFriendInt["chargino_simHits_subDet"]->at(i)]));
    }
    
    // Parameters for all hits in the pixel barrel
    for(uint i=0;i<arrayValuesFriendFloat["pixelCluster_x"]->size();i++){
      // convert cm to mm
      event->trackerClusters.push_back(make_shared<Point>(10*arrayValuesFriendFloat["pixelCluster_x"]->at(i),
                                                          10*arrayValuesFriendFloat["pixelCluster_y"]->at(i),
                                                          10*arrayValuesFriendFloat["pixelCluster_z"]->at(i),
                                                          arrayValuesFriendFloat["pixelCluster_charge"]->at(i),
                                                          subDetMap[arrayValuesFriendInt["pixelCluster_subDet"]->at(i)]));
    }
    
    for(uint i=0;i<arrayValuesFriendFloat["stripCluster_x"]->size();i++){
      string detName = subDetMap[arrayValuesFriendInt["stripCluster_subDet"]->at(i)];
      // convert cm to mm
      if(detName == "P1PXEC" || detName == "TID" || detName == "TEC"){
        event->trackerClusters.push_back(make_shared<Point>(10*arrayValuesFriendFloat["stripCluster_x"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_y"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_z"]->at(i),
                                                            arrayValuesFriendFloat["stripCluster_charge"]->at(i),
                                                            subDetMap[arrayValuesFriendInt["stripCluster_subDet"]->at(i)],
                                                            10*arrayValuesFriendFloat["stripCluster_ex"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_ez"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_ey"]->at(i)));
      }
      else{
        
        event->trackerClusters.push_back(make_shared<Point>(10*arrayValuesFriendFloat["stripCluster_x"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_y"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_z"]->at(i),
                                                            arrayValuesFriendFloat["stripCluster_charge"]->at(i),
                                                            subDetMap[arrayValuesFriendInt["stripCluster_subDet"]->at(i)],
                                                            10*arrayValuesFriendFloat["stripCluster_ex"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_ey"]->at(i),
                                                            10*arrayValuesFriendFloat["stripCluster_ez"]->at(i)));
      }
    }
    
    for(uint i=0;i<arrayValuesFriendFloat["pionCluster_x"]->size();i++){
      string detName = subDetMap[arrayValuesFriendInt["pionCluster_subDet"]->at(i)];
      // convert cm to mm
      if(detName == "P1PXEC" || detName == "TID" || detName == "TEC"){
        event->pionClusters.push_back(make_shared<Point>(10*arrayValuesFriendFloat["pionCluster_x"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_y"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_z"]->at(i),
                                                         arrayValuesFriendFloat["pionCluster_charge"]->at(i),
                                                         subDetMap[arrayValuesFriendInt["pionCluster_subDet"]->at(i)],
                                                         10*arrayValuesFriendFloat["pionCluster_ex"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_ez"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_ey"]->at(i)));
      }
      else{
        event->pionClusters.push_back(make_shared<Point>(10*arrayValuesFriendFloat["pionCluster_x"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_y"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_z"]->at(i),
                                                         arrayValuesFriendFloat["pionCluster_charge"]->at(i),
                                                         subDetMap[arrayValuesFriendInt["pionCluster_subDet"]->at(i)],
                                                         10*arrayValuesFriendFloat["pionCluster_ex"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_ey"]->at(i),
                                                         10*arrayValuesFriendFloat["pionCluster_ez"]->at(i)));
      }
    }
  }
  
  auto charginoSimHitsByLayer = pointsProcessor.SortByLayer(event->charginoSimHits);
  int maxCharginoLayer = -1;
  for(auto &hits : charginoSimHitsByLayer){
    for(auto &hit : hits){
      if(hit->GetLayer() > maxCharginoLayer) maxCharginoLayer = hit->GetLayer();
    }
  }
  int nCharginoLayers = event->charginoSimHits.size()==0 ? 0 : maxCharginoLayer+1;
  int nCharginos = (int)arrayValuesFriendFloat["chargino_eta"]->size();
  
  for(uint i=0;i<nCharginos;i++){
    event->genCharginoTrack.push_back(Track(arrayValuesFriendFloat["chargino_eta"]->at(i),
                                            arrayValuesFriendFloat["chargino_phi"]->at(i),
                                            arrayValuesFriendInt["chargino_charge"]->at(i),
                                            nCharginos==1 ? nCharginoLayers : -1, // works only if there's one chargino in the event!!
                                            /*arrayValuesFriendInt["chargino_nTrackerLayers"]->at(i),*/
                                            arrayValuesFriendFloat["chargino_pt"]->at(i),
                                            arrayValuesFriendFloat["chargino_mass"]->at(i)));
    
  }
  
  return event;
}

void EventProcessor::SaveEventToTree(shared_ptr<Event> event)
{
  singleValuesUint["lumi"]      = event->lumiSection;
  singleValuesUint["run"]       = event->runNumber;
  singleValuesUint["wasTagged"] = event->wasTagged;
  singleValuesUlonglong["evt"]  = event->eventNumber;
  
  singleValuesInt["nVert"]                                         = event->nVertices;
  singleValuesInt["nJet30"]                                        = event->nJet30;
  singleValuesInt["nJet30a"]                                       = event->nJet30a;
  singleValuesInt["nTauGood"]                                      = event->nTau;
  singleValuesInt["nGenChargino"]                                  = event->nGenChargino;
  singleValuesInt["HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"] = event->metNoMuTrigger;
  singleValuesInt["Flag_goodVertices"]                             = event->flag_goodVertices;
  singleValuesInt["Flag_BadPFMuonFilter"]                          = event->flag_badPFmuon;
  singleValuesInt["Flag_HBHENoiseFilter"]                          = event->flag_HBHEnoise;
  singleValuesInt["Flag_HBHENoiseIsoFilter"]                       = event->flag_HBHEnoiseIso;
  singleValuesInt["Flag_EcalDeadCellTriggerPrimitiveFilter"]       = event->flag_EcalDeadCell;
  singleValuesInt["Flag_eeBadScFilter"]                            = event->flag_eeBadSc;
  singleValuesInt["Flag_BadChargedCandidateFilter"]                = event->flag_badChargedCandidate;
  singleValuesInt["Flag_ecalBadCalibFilter"]                       = event->flag_ecalBadCalib;
  singleValuesInt["Flag_globalTightHalo2016Filter"]                = event->flag_globalTightHalo2016;
  
  singleValuesFloat["vertex_x"]     = event->vertex->GetX();
  singleValuesFloat["vertex_y"]     = event->vertex->GetY();
  singleValuesFloat["vertex_z"]     = event->vertex->GetZ();
  singleValuesFloat["xsec"]         = event->xsec;
  singleValuesFloat["wgtsum"]       = event->wgtsum;
  singleValuesFloat["genWeight"]    = event->genWeight;
  singleValuesFloat["met_sumEt"]    = event->metSumEt;
  singleValuesFloat["met_pt"]       = event->metPt;
  singleValuesFloat["met_GenPt"]    = event->metGenPt;
  singleValuesFloat["met_mass"]     = event->metMass;
  singleValuesFloat["met_phi"]      = event->metPhi;
  singleValuesFloat["met_eta"]      = event->metEta;
  singleValuesFloat["metNoMu_pt"]     = event->metNoMuPt;
  singleValuesFloat["met_jecUp_pt"]   = event->metNoMuPtJecUp;
  singleValuesFloat["met_jecDown_pt"] = event->metNoMuPtJecDown;
  singleValuesFloat["metNoMu_mass"] = event->metNoMuMass;
  singleValuesFloat["metNoMu_phi"]  = event->metNoMuPhi;
  singleValuesFloat["metNoMu_eta"]  = event->metNoMuEta;
  
  for(uint i=0;i<event->nGenChargino;i++){
    arrayValuesFloat["GenChargino_eta"][i]    = event->genCharginos[i].GetEta();
    arrayValuesFloat["GenChargino_phi"][i]    = event->genCharginos[i].GetPhi();
    arrayValuesFloat["GenChargino_charge"][i] = event->genCharginos[i].GetCharge();
    arrayValuesFloat["GenChargino_pt"][i]     = event->genCharginos[i].GetPt();
    arrayValuesFloat["GenChargino_mass"][i]   = event->genCharginos[i].GetMass();
   }
  
}

void EventProcessor::SetupBranchesForReading(TTree *tree, TTree *friendTree, TTree *prefireTree)
{
  for(string name : singleNamesFloat){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      singleValuesFloat[name] = 0;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &singleValuesFloat[name]);
  }
  
  for(string name : arrayNamesFloat){
//    arrayValuesFloat[name] = nullptr;
    
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesFloat[name]);
  }
  
  for(string name : singleNamesInt){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      singleValuesInt[name] = 0;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &singleValuesInt[name]);
  }
  
  for(string name : singleNamesUint){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      singleValuesUint[name] = 0;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &singleValuesUint[name]);
  }
  
  for(string name : singleNamesUlongLong){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      singleValuesUlonglong[name] = 0;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &singleValuesUlonglong[name]);
  }
  
  if(prefireTree){
    if(!prefireTree->GetBranchStatus("prefiringweight")){
      cout<<"WARNING -- no branch named prefiringweight"<<endl;
      prefireWeight = 1.0;
    }
    else{
      prefireTree->SetBranchAddress("prefiringweight", &prefireWeight);
    }
  }
  
  if(!friendTree) return;
  
  for(string name : arrayNamesFriendFloat){
    arrayValuesFriendFloat[name] = nullptr;
    
    if(!friendTree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    friendTree->SetBranchAddress(name.c_str(), &arrayValuesFriendFloat[name]);
  }
  
  for(string name : arrayNamesFriendInt){
    arrayValuesFriendInt[name] = nullptr;
    
    if(!friendTree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    friendTree->SetBranchAddress(name.c_str(), &arrayValuesFriendInt[name]);
  }
}

void EventProcessor::SetupBranchesForWriting(TTree *tree)
{
  for(string name : singleNamesFloat){
    tree->Branch(name.c_str(), &singleValuesFloat[name], Form("%s/F", name.c_str()));
  }
  
  for(string name : singleNamesInt){
    tree->Branch(name.c_str(), &singleValuesInt[name], Form("%s/I", name.c_str()));
  }
  
  for(string name : singleNamesUint){
    tree->Branch(name.c_str(), &singleValuesUint[name], Form("%s/i", name.c_str()));
  }
  
  for(string name : singleNamesUlongLong){
    tree->Branch(name.c_str(), &singleValuesUlonglong[name], Form("%s/l", name.c_str()));
  }
  
  for(string name : arrayNamesFloat){
     tree->Branch(name.c_str(), &arrayValuesFloat[name], Form("%s[nGenChargino]/F", name.c_str()));
   }
}
