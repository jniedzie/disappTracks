//
//  EventProcessor.cpp
//
//  Created by Jeremi Niedziela on 30/01/2019.
//

#include "EventProcessor.hpp"

EventProcessor::EventProcessor() :
trackProcessor(make_unique<TrackProcessor>()),
jetProcessor(make_unique<JetProcessor>()),
leptonProcessor(make_unique<LeptonProcessor>())
{
  singleNamesUint = {
    "run",
    "lumi"
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
    "met_mass",
    "met_phi",
    "met_eta",
    "metNoMu_pt",
    "metNoMu_mass",
    "metNoMu_phi",
    "metNoMu_eta"
  };
  
  singleNamesInt = {
    "nVert",
    "nJet30",
    "nJet30a",
    "nTauGood",
    "nGenChargino",
    "HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
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
  
}

EventProcessor::~EventProcessor()
{
  
}

void EventProcessor::ApplyTrackCut(shared_ptr<Event> event, const unique_ptr<TrackCut> &cut)
{
  for(auto track = event->tracks.begin(); track != event->tracks.end();){
    if(!trackProcessor->IsPassingCut(*track,cut)) track = event->tracks.erase(track);
    else                                          track++;
  }
}

void EventProcessor::ApplyJetCut(shared_ptr<Event> event, const unique_ptr<JetCut> &cut)
{
  for(auto jet = event->jets.begin(); jet != event->jets.end();){
    
    if(!jetProcessor->IsPassingCut(*jet,cut)){
      jet = event->jets.erase(jet);
    }
    else{
      // check separation with all tracks in the event
      bool overlapsWithTrack = false;
      double minTrackDeltaR = cut->GetTrackDeltaR().GetMin();
      
      if(minTrackDeltaR > 0){
        for(auto track : event->tracks){
          double deltaR_2 =   pow(track->GetPhi() - (*jet)->GetPhi(),2)
                            + pow(track->GetEta() - (*jet)->GetEta(),2);
          
          if(deltaR_2 < (minTrackDeltaR*minTrackDeltaR)){
            overlapsWithTrack = true;
            break;
          }
        }
      }
      
      if(overlapsWithTrack) jet = event->jets.erase(jet);
      else                  jet++;
    }
  }
  
}

void EventProcessor::ApplyLeptonCut(shared_ptr<Event> event, const unique_ptr<LeptonCut> &cut)
{
  for(auto lepton = event->leptons.begin(); lepton != event->leptons.end();){
    if(!leptonProcessor->IsPassingCut(*lepton, cut))  lepton = event->leptons.erase(lepton);
    else                                              lepton++;
  }
}


bool EventProcessor::IsPassingCut(const shared_ptr<Event> event, const EventCut &cut)
{
  // check the trigger
  if(cut.requiresMetNoMuTrigger && !event->metNoMuTrigger){
    cutReasons[0]++;
    return false;
  }
  
  // check MET filters
  if(cut.requiresPassingAllFilters){
    if(   !event->flag_goodVertices
       || !event->flag_badPFmuon
       || !event->flag_HBHEnoise
       || !event->flag_HBHEnoiseIso
       || !event->flag_EcalDeadCell
       || !event->flag_eeBadSc
     //|| flag_badChargedCandidate/
       || !event->flag_ecalBadCalib
       || !event->flag_globalTightHalo2016){
      cutReasons[1]++;
      return false;
    }
  }
  
  // check number of objects
  if(cut.nLeptons.IsOutside(event->GetNleptons())){
    cutReasons[2]++;
    return false;
  }
  if(cut.nTaus.IsOutside(event->nTau)){
    cutReasons[3]++;
    return false;
  }
  
  vector<shared_ptr<Lepton>> muons;
  for(auto l : event->leptons){
    if(abs(l->GetPid()) == 13) muons.push_back(l);
  }
  
  // check number of muons
  if(cut.nMuons.IsOutside((int)muons.size())){
    cutReasons[4]++;
    return false;
  }
  
  // make sure they have an opposite sign
  if(cut.requiresTwoOpositeMuons){
    if(muons.size() != 2){
      cutReasons[5]++;
      return false;
    }
    if(muons[0]->GetPid() != -muons[1]->GetPid()){
      cutReasons[6]++;
      return false;
    }
  }
  
  // apply tight muon cuts (tightID flag, pt > 20 GeV, isolation < 0.15
  if(cut.requiresTightMuon){
    auto tightMuonCut = LeptonCut();
    tightMuonCut.SetRelativeIsolation(range<double>(-inf,0.15));
    tightMuonCut.SetRequireTightID(true);
    tightMuonCut.SetPt(range<double>(20.0,inf));
    
    bool atLeastOneTightMuon = false;
    for(auto muon : muons){
      if(leptonProcessor->IsPassingCut(muon, make_unique<LeptonCut>(tightMuonCut))) atLeastOneTightMuon = true;
    }
    if(!atLeastOneTightMuon){
      cutReasons[7]++;
      return false;
    }
  }
  
  // check that invariant mass of muons is close to Z mass
  if(cut.requiresMuonsFromZ){
    TLorentzVector muon1vector, muon2vector, muonVectorSum;
    
    if(muons.size() != 2){
      cout<<"ERROR -- requested muons to come from Z decay, but there is "<<muons.size()<<" muons in the event!!"<<endl;
      cout<<"This event will be discarded!! Maybe you should require exactly two muons in the event?"<<endl;
      cutReasons[8]++;
      return false;
    }
    
    muon1vector.SetPtEtaPhiM(muons[0]->GetPt(),muons[0]->GetEta(),muons[0]->GetPhi(), 0.1057);
    muon2vector.SetPtEtaPhiM(muons[1]->GetPt(),muons[1]->GetEta(),muons[1]->GetPhi(), 0.1057);
    
    muonVectorSum += muon1vector;
    muonVectorSum += muon2vector;
    
    if(muonVectorSum.M() < 60. || muonVectorSum.M() > 120.){
      cutReasons[9]++;
      return false;
    }
  }
  
  if(cut.metPt.IsOutside(event->metPt)){
    cutReasons[10]++;
    return false;
  }
  if(cut.metNoMuPt.IsOutside(event->metNoMuPt)){
    cutReasons[11]++;
    return false;
  }
  
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
          if(cut.jetMuonDeltaPhi.IsOutside(fabs(jetVector.DeltaR(muonVector)))){
            return false;
          }
        }
      }
    }
  }
  
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
        if(cut.trackMuonDeltaPhi.IsOutside(fabs(trackVector.DeltaR(muonVector)))){
          return false;
        }
      }
    }
  }
  
  // check number of tracks and jets after removing those that are too close to muons
  if(cut.nJets.IsOutside(event->GetNcentralJets())){
   cutReasons[12]++;
    return false;
  }
  if(cut.nTracks.IsOutside(event->GetNtracks())){
    cutReasons[13]++;
    return false;
  }
  
  // find the jet with the highest pt meeting leading jet criteria
  shared_ptr<Jet> leadingJet = nullptr;
  double highestPt = -inf;
  
  for(auto jet : event->jets){
      if(jet->GetPt() > highestPt){
        highestPt = jet->GetPt();
        leadingJet = jet;
      }
  }
  
  if(cut.leadingJetPt.IsOutside(leadingJet->GetPt())        ||
     cut.leadingJetEta.IsOutside(leadingJet->GetEta())      ||
     cut.leadingJetChHEF.IsOutside(leadingJet->GetChHEF())  ||
     cut.leadingJetNeHEF.IsOutside(leadingJet->GetNeHEF())){
    leadingJet = nullptr;
  }
  
  // check if there is a leading jet in the event
  if(!leadingJet){
    cutReasons[14]++;
    return false;
  }

  if(cut.jetMetDeltaPhi.GetMin() > 0.0){
    TLorentzVector metVector, jetVector;
//    metVector.SetPtEtaPhiM(event->metPt, event->metEta, event->metPhi, event->metMass);
    metVector.SetPtEtaPhiM(event->metNoMuPt, event->metNoMuEta, event->metNoMuPhi, event->metNoMuMass);
    
    for(auto j : event->jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(cut.jetMetDeltaPhi.IsOutside(fabs(metVector.DeltaPhi(jetVector)) )){
        cutReasons[15]++;
        return false;
      }
    }
  }
  
  survivingEvents.push_back(event);
  return true;
}

shared_ptr<Event> EventProcessor::GetEventFromTree(xtracks::EDataType dataType, int setIter)
{
  for(auto &[name, val] : singleValuesInt ){
    if(val < -999999){
      cout<<"ERROR -- branch "<<name<<" was not read correctly!"<<endl;
    }
  }
  
  auto event = make_shared<Event>();
  
  double lumi = config.totalLuminosity * 1000.; // transform from fb^-1 to pb^-1
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
    weight *= 0.001 * (signalCrossSectionOneTrack[(ESignal)setIter] +
                       signalCrossSectionTwoTracks[(ESignal)setIter]);
    
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
  event->weight   = weight;
  event->dataType = dataType;
  event->setIter  = setIter;
  
  event->lumiSection = singleValuesUint["lumi"];
  event->runNumber   = singleValuesUint["run"];
  event->eventNumber = singleValuesUlonglong["evt"];
  
  event->nVertices                = singleValuesInt["nVert"];
  event->nJet30                   = singleValuesInt["nJet30"];
  event->nJet30a                  = singleValuesInt["nJet30a"];
  event->nTau                     = singleValuesInt["nTauGood"];
  event->nGenChargino             = singleValuesInt["nGenChargino"];
  event->metNoMuTrigger           = singleValuesInt["HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"];
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
  event->metMass     = singleValuesFloat["met_mass"];
  event->metPhi      = singleValuesFloat["met_phi"];
  event->metEta      = singleValuesFloat["met_eta"];
  event->metNoMuPt   = singleValuesFloat["metNoMu_pt"];
  event->metNoMuMass = singleValuesFloat["metNoMu_mass"];
  event->metNoMuPhi  = singleValuesFloat["metNoMu_phi"];
  event->metNoMuEta  = singleValuesFloat["metNoMu_eta"];
  
  return event;
}

void EventProcessor::SaveEventToTree(shared_ptr<Event> event)
{
  singleValuesUint["lumi"]      = event->lumiSection;
  singleValuesUint["run"]       = event->runNumber;
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
  singleValuesFloat["met_mass"]     = event->metMass;
  singleValuesFloat["met_phi"]      = event->metPhi;
  singleValuesFloat["met_eta"]      = event->metEta;
  singleValuesFloat["metNoMu_pt"]   = event->metNoMuPt;
  singleValuesFloat["metNoMu_mass"] = event->metNoMuMass;
  singleValuesFloat["metNoMu_phi"]  = event->metNoMuPhi;
  singleValuesFloat["metNoMu_eta"]  = event->metNoMuEta;
}

void EventProcessor::SetupBranchesForReading(TTree *tree)
{
  for(string name : singleNamesFloat){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      singleValuesFloat[name] = 0;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &singleValuesFloat[name]);
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
}
