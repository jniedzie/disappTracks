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
  auto track = event->tracks.begin();
  
  while(track != event->tracks.end()){
    if(!trackProcessor->IsPassingCut(*track,cut)) track = event->tracks.erase(track);
    else                                          track++;
  }
}

void EventProcessor::ApplyJetCut(shared_ptr<Event> event, const unique_ptr<JetCut> &cut)
{
  auto jet = event->jets.begin();
  
  while(jet != event->jets.end()){
    if(!jetProcessor->IsPassingCut(*jet,cut))
    jet = event->jets.erase(jet);
    else{
      // check separation with all tracks in the event
      bool overlapsWithTrack = false;
      double minTrackDeltaR = cut->GetTrackDeltaR().GetMin();
      
      if(minTrackDeltaR > 0){
        for(auto track : event->tracks){
          double deltaR_2 = pow(track->GetPhi() - (*jet)->GetPhi(),2)+pow(track->GetEta() - (*jet)->GetEta(),2);
          
          if(deltaR_2 < (minTrackDeltaR*minTrackDeltaR)){
            overlapsWithTrack = true;
            break;
          }
        }
      }
      
      if(overlapsWithTrack){
        jet = event->jets.erase(jet);
      }
      else{
        jet++;
      }
    }
  }
}

void EventProcessor::ApplyLeptonCut(shared_ptr<Event> event, const unique_ptr<LeptonCut> &cut)
{
  auto lepton = event->leptons.begin();
  
  while(lepton != event->leptons.end()){
    if(!leptonProcessor->IsPassingCut(*lepton, cut))
    lepton = event->leptons.erase(lepton);
    else
    lepton++;
  }
}


bool EventProcessor::IsPassingCut(const shared_ptr<Event> event, const unique_ptr<EventCut> &cut)
{
  // check the trigger
  if(cut->RequiresMetNoMuTrigger() && !event->metNoMuTrigger){
    return false;
  }
  // check filters
  if(cut->GetRequiresPassingAllFilters()){
    if(   !event->flag_goodVertices
       || !event->flag_goodVertices
       || !event->flag_badPFmuon
       || !event->flag_HBHEnoise
       || !event->flag_HBHEnoiseIso
       || !event->flag_EcalDeadCell
       || !event->flag_eeBadSc
    /*|| flag_badChargedCandidate*/
       || !event->flag_ecalBadCalib
       || !event->flag_globalTightHalo2016){
      return false;
    }
  }
  
  // check number of objects
  if(cut->GetNleptons().IsOutside(event->nLepton))  return false;
  if(cut->GetNtaus().IsOutside(event->nTau)) return false;
  
  vector<shared_ptr<Lepton>> muons;
  for(auto l : event->leptons){
    if(abs(l->GetPid()) == 13) muons.push_back(l);
  }
  
  // check number of muons
  if(cut->GetNmuons().IsOutside((int)muons.size())) return false;
  
  // make sure they have an opposite sign
  if(cut->RequiresTwoOppositeMuons()){
    if(muons.size() != 2) return false;
    if(muons[0]->GetPid() != -muons[1]->GetPid()) return false;
  }
  
  // apply tight muon cuts (tightID flag, pt > 20 GeV, isolation < 0.15
  if(cut->RequiresTightMuon()){
    unique_ptr<LeptonCut> tightMuonCut = unique_ptr<LeptonCut>(new LeptonCut());
    tightMuonCut->SetRelativeIsolation(range<double>(-inf,0.15));
    tightMuonCut->SetRequireTightID(true);
    tightMuonCut->SetPt(range<double>(20.0,inf));
    
    bool atLeastOneTightMuon = false;
    for(auto muon : muons){
      if(leptonProcessor->IsPassingCut(muon, tightMuonCut)) atLeastOneTightMuon = true;
    }
    if(!atLeastOneTightMuon) return false;
  }
  
  // check that invariant mass of muons is close to Z mass
  if(cut->RequiresMuonsFromZ()){
    TLorentzVector muon1vector, muon2vector, muonVectorSum;
    
    if(muons.size() != 2){
      cout<<"ERROR -- requested muons to come from Z decay, but there is "<<muons.size()<<" muons in the event!!"<<endl;
      cout<<"This event will be discarded!! Maybe you should require exactly two muons in the event?"<<endl;
      return false;
    }
    
    muon1vector.SetPtEtaPhiM(muons[0]->GetPt(),muons[0]->GetEta(),muons[0]->GetPhi(), 0.1057);
    muon2vector.SetPtEtaPhiM(muons[1]->GetPt(),muons[1]->GetEta(),muons[1]->GetPhi(), 0.1057);
    
    muonVectorSum += muon1vector;
    muonVectorSum += muon2vector;
    
    if(muonVectorSum.M() < 60. || muonVectorSum.M() > 120.) return false;
  }
  
  if(cut->GetMetPt().IsOutside(event->metPt))  return false;
  if(cut->GetMetNoMuPt().IsOutside(event->metNoMuPt))  return false;
  
  // Remove jets that are too close to muons (they will be permanently removed from the event)
  if(cut->RequiresMuJetR0p4()){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), 0.1057);
      muonVectors.push_back(mv);
    }
    
    for(int iJet=0;iJet<event->GetNcentralJets();iJet++){
      shared_ptr<Jet> j = event->GetJet(iJet);
      if(!j->IsForward()){
        
        TLorentzVector jetVector;
        jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
        
        for(auto muonVector : muonVectors){
          if(jetVector.DeltaR(muonVector) < 0.4){
            event->jets.erase(event->jets.begin()+iJet);
            iJet--;
            break;
          }
        }
      }
    }
  }
  
  // Remove tracks that are too close to muons (they will be permanently removed from the event)
  if(cut->RequiresMuTrackR0p4()){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), 0.1057);
      muonVectors.push_back(mv);
    }
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      shared_ptr<Track> t = event->tracks[iTrack];
      
      TLorentzVector trackVector;
      trackVector.SetPtEtaPhiM(t->GetPt(), t->GetEta(), t->GetPhi(), t->GetMass());
      
      for(auto muonVector : muonVectors){
        if(trackVector.DeltaR(muonVector) < 0.4){
          event->tracks.erase(event->tracks.begin()+iTrack);
          iTrack--;
          break;
        }
      }
    }
  }
  
  // check number of tracks and jets after removing those that are too close to muons
  if(cut->GetNjets().IsOutside(event->GetNcentralJets())) return false;
  if(cut->GetNtracks().IsOutside(event->GetNtracks())) return false;
  
  // find the jet with the highest pt
  shared_ptr<Jet> leadingJet = nullptr;
  double highestPt = -1.0;
  
  for(int iJet=0;iJet<event->GetNjets();iJet++){
    if(event->jets[iJet]->GetPt() > highestPt){
      highestPt = event->jets[iJet]->GetPt();
      leadingJet = event->jets[iJet];
    }
  }
  
  // check properties of the highest pt jet
  if(cut->RequiresHighJet() && !leadingJet) return false;
  if(cut->GetLeadingJetPt().IsOutside(leadingJet->GetPt())) return false;
  if(cut->GetLeadingJetEta().IsOutside(leadingJet->GetEta())) return false;
  if(cut->GetLeadingJetChHEF().IsOutside(leadingJet->GetChargedHadronEnergyFraction())) return false;
  if(cut->GetLeadingJetNeHEF().IsOutside(leadingJet->GetNeutralHadronEnergyFraction())) return false;
  
  if(cut->GetJetMetDeltaPhi().GetMin() > 0.0){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(event->metPt, event->metEta, event->metPhi, event->metMass);
    
    for(auto j : event->jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(cut->GetJetMetDeltaPhi().IsOutside(fabs(metVector.DeltaPhi(jetVector)) )) return false;
    }
  }
  
  if(cut->RequiresMetNoMuJetPhi0p5()){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(event->metNoMuPt, event->metNoMuEta, event->metNoMuPhi, event->metNoMuMass);
    
    for(auto j : event->jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(fabs(metVector.DeltaPhi(jetVector)) < 0.5) return false;
    }
  }
  
  return true;
}

shared_ptr<Event> EventProcessor::GetEventFromTree(xtracks::EDataType dataType, int setIter)
{
  auto event = make_shared<Event>();
  
  double lumi = config->totalLuminosity * 1000.; // transform from fb^-1 to pb^-1
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
      singleValuesFloat[name] = 0;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &singleValuesFloat[name]);
  }
  
  for(string name : singleNamesInt){
    tree->SetBranchAddress(name.c_str(), &singleValuesInt[name]);
  }
  
  for(string name : singleNamesUint){
    tree->SetBranchAddress(name.c_str(), &singleValuesUint[name]);
  }
  
  for(string name : singleNamesUlongLong){
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
