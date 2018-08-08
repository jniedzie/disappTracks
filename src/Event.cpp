//
//  Event.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Event.hpp"

#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
//#include <ROOT/TTreeProcessorMT.hxx>

Events::Events()
{
  
}

Events::Events(const Events &e)
{
  for(auto event : e.events){
    events.push_back(event);
  }
}

Events::Events(string fileName, int dataType)
{
  cout<<"Reading events from:"<<fileName<<endl;
  
  map<unsigned long long,Event*> eventsMap;
  
  TFile *inFile = TFile::Open(fileName.c_str());
  TTreeReader reader("tree", inFile);

  TTreeReaderValue<int>   nTracks(reader, "nIsoTrack");
  TTreeReaderValue<int>   _nVert(reader, "nVert");
  TTreeReaderValue<int>   nJets(reader, "nJet");
  TTreeReaderValue<int>   nJetsFwd(reader, "nJetFwd");
  TTreeReaderValue<int>   _nJet30(reader, "nJet30");
  TTreeReaderValue<int>   _nJet30a(reader, "nJet30a");
  TTreeReaderValue<int>   _nLepton(reader, "nLepGood");
  TTreeReaderValue<int>   _nTau(reader, "nTauGood");
  
  TTreeReaderValue<float> _xsec  (reader,(dataType==0 || dataType==1) ? "xsec" : "rho");
  TTreeReaderValue<float> _wgtsum(reader,(dataType==0 || dataType==1) ? "wgtsum" : "rho");
  TTreeReaderValue<float> _genwgt(reader,(dataType==0 || dataType==1) ? "genWeight" : "rho");

  TTreeReaderValue<float> _met_sumEt(reader, "met_sumEt");
  TTreeReaderValue<float> _met_pt(reader, "met_pt");
  TTreeReaderValue<float> _met_mass(reader, "met_mass");
  TTreeReaderValue<float> _met_phi(reader, "met_phi");
  TTreeReaderValue<float> _met_eta(reader, "met_eta");
  
  TTreeReaderValue<int>   _metNoMuTrigger(reader, "HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
  
  TTreeReaderValue<float> _metNoMu_pt(reader, "metNoMu_pt");
  TTreeReaderValue<float> _metNoMu_mass(reader, "metNoMu_mass");
  TTreeReaderValue<float> _metNoMu_phi(reader, "metNoMu_phi");
  TTreeReaderValue<float> _metNoMu_eta(reader, "metNoMu_eta");

  TTreeReaderArray<float> _eta(reader, "IsoTrack_eta");
  TTreeReaderArray<float> _phi(reader, "IsoTrack_phi");
  TTreeReaderArray<float> _caloEmEnergy(reader, "IsoTrack_caloEmEnergy");
  TTreeReaderArray<float> _caloHadEnergy(reader, "IsoTrack_caloHadEnergy");
  TTreeReaderArray<float> _dxyErr(reader, "IsoTrack_edxy");
  TTreeReaderArray<float> _dxy(reader, "IsoTrack_dxy");
  TTreeReaderArray<float> _dzErr(reader, "IsoTrack_edz");
  TTreeReaderArray<float> _dz(reader, "IsoTrack_dz");
  TTreeReaderArray<int>   _charge(reader, "IsoTrack_charge");
  TTreeReaderArray<float> _mass(reader, "IsoTrack_mass");
  TTreeReaderArray<float> _pt(reader, "IsoTrack_pt");
  TTreeReaderArray<int>   _pid(reader, "IsoTrack_pdgId");

  TTreeReaderArray<float> _lepton_pt(reader, "LepGood_pt");
  TTreeReaderArray<float> _lepton_phi(reader, "LepGood_phi");
  TTreeReaderArray<float> _lepton_eta(reader, "LepGood_eta");
  TTreeReaderArray<int>   _lepton_thight_pid(reader, "LepGood_tightId");
  TTreeReaderArray<float> _lepton_isolation(reader, "LepGood_relIso04");
  TTreeReaderArray<int>   _lepton_pid(reader, "LepGood_pdgId");
  
  TTreeReaderArray<float> _jet_pt(reader,  "Jet_pt");
  TTreeReaderArray<float> _jet_eta(reader, "Jet_eta");
  TTreeReaderArray<float> _jet_phi(reader, "Jet_phi");
  TTreeReaderArray<float> _jet_mass(reader, "Jet_mass");
  TTreeReaderArray<float> _jet_chHEF(reader, "Jet_chHEF");
  TTreeReaderArray<float> _jet_neHEF(reader, "Jet_neHEF");
  
  TTreeReaderArray<float> _jetFwd_pt(reader,  "JetFwd_pt");
  TTreeReaderArray<float> _jetFwd_eta(reader, "JetFwd_eta");
  TTreeReaderArray<float> _jetFwd_phi(reader, "JetFwd_phi");
  TTreeReaderArray<float> _jetFwd_mass(reader, "JetFwd_mass");
  TTreeReaderArray<float> _jetFwd_chHEF(reader, "JetFwd_chHEF");
  TTreeReaderArray<float> _jetFwd_neHEF(reader, "JetFwd_neHEF");

  
  TTreeReaderArray<float> *dedx[nLayers];
  TTreeReaderArray<int> *subDetId[nLayers];
  TTreeReaderArray<int> *sizeX[nLayers];
  TTreeReaderArray<int> *sizeY[nLayers];

  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx[iLayer] =      new TTreeReaderArray<float>(reader,Form("IsoTrack_dedxByLayer%i",iLayer));
    subDetId[iLayer] =  new TTreeReaderArray<int>(reader,Form("IsoTrack_subDetIdByLayer%i",iLayer));
    sizeX[iLayer] =     new TTreeReaderArray<int>(reader,Form("IsoTrack_sizeXbyLayer%i",iLayer));
    sizeY[iLayer] =     new TTreeReaderArray<int>(reader,Form("IsoTrack_sizeYbyLayer%i",iLayer));
  }
  int iter=0;
  while (reader.Next()){
//    if(iter>10000) break;
    iter++;
    
    Event *newEvent = new Event();

    for(int iTrack=0;iTrack<*nTracks;iTrack++){
      Track *track = new Track();
      track->SetEta(_eta[iTrack]);
      track->SetPhi(_phi[iTrack]);
      track->SetCaloEmEnergy(_caloEmEnergy[iTrack]);
      track->SetCaloHadEnergy(_caloHadEnergy[iTrack]);
      track->SetDxy(_dxy[iTrack],_dxyErr[iTrack]);
      track->SetDz(_dz[iTrack],_dzErr[iTrack]);
      track->SetCharge(_charge[iTrack]);
      track->SetMass(_mass[iTrack]);
      track->SetPt(_pt[iTrack]);
      track->SetPid(_pid[iTrack]);

      for(int iLayer=0;iLayer<nLayers;iLayer++){
        track->SetDeDxInLayer(iLayer, (*dedx[iLayer])[iTrack]);
        track->SetSubDetIdInLayer(iLayer, (*subDetId[iLayer])[iTrack]);
        track->SetSizeXinLayer(iLayer, (*sizeX[iLayer])[iTrack]);
        track->SetSizeYinLayer(iLayer, (*sizeY[iLayer])[iTrack]);
      }
      newEvent->AddTrack(track);
    }
    
    for(int iJet=0;iJet<*nJets;iJet++){
      Jet *jet = new Jet();
      jet->SetPt(_jet_pt[iJet]);
      jet->SetEta(_jet_eta[iJet]);
      jet->SetPhi(_jet_phi[iJet]);
      jet->SetMass(_jet_mass[iJet]);
      jet->SetChargedHadronEnergyFraction(_jet_chHEF[iJet]);
      jet->SetNeutralHadronEnergyFraction(_jet_neHEF[iJet]);
      jet->SetIsForward(false);
      newEvent->AddJet(jet);
    }
    
    for(int iJet=0;iJet<*nJetsFwd;iJet++){
      Jet *jet = new Jet();
      jet->SetPt(_jetFwd_pt[iJet]);
      jet->SetEta(_jetFwd_eta[iJet]);
      jet->SetPhi(_jetFwd_phi[iJet]);
      jet->SetMass(_jetFwd_mass[iJet]);
      jet->SetChargedHadronEnergyFraction(_jetFwd_chHEF[iJet]);
      jet->SetNeutralHadronEnergyFraction(_jetFwd_neHEF[iJet]);
      jet->SetIsForward(true);
      newEvent->AddJet(jet);
    }
    
    for(int iLepton=0;iLepton<*_nLepton;iLepton++){
      Lepton *lepton = new Lepton();
      lepton->SetPt(_lepton_pt[iLepton]);
      lepton->SetEta(_lepton_eta[iLepton]);
      lepton->SetPhi(_lepton_phi[iLepton]);
      lepton->SetTightID(_lepton_thight_pid[iLepton]);
      lepton->SetIsolation(_lepton_isolation[iLepton]);
      lepton->SetPid(_lepton_pid[iLepton]);
      newEvent->AddLepton(lepton);
    }
    
    double lumi = 41.37 * 1000.;
    double weight = lumi * (*_xsec) * (*_genwgt) / (*_wgtsum);

    if(dataType==1){
      weight = 182; // just invented some number to make S/B ~ 1
      weight *= 10000.0/reader.GetEntries(true); // correct for less entries in the tree than for background
    }
    else if(dataType==2){
      weight = 1.0;
    }
    
    newEvent->SetWeight(weight);
    
    newEvent->SetNvertices(*_nVert);
    newEvent->SetNjet30(*_nJet30);
    newEvent->SetNjet30a(*_nJet30a);
    newEvent->SetNlepton(*_nLepton);
    newEvent->SetNtau(*_nTau);
    
    newEvent->SetMetSumEt(*_met_sumEt);
    newEvent->SetMetPt(*_met_pt);
    newEvent->SetMetMass(*_met_mass);
    newEvent->SetMetEta(*_met_eta);
    newEvent->SetMetPhi(*_met_phi);
    
    newEvent->SetHasNoMuTrigger(*_metNoMuTrigger);
    newEvent->SetMetNoMuPt(*_metNoMu_pt);
    newEvent->SetMetNoMuMass(*_metNoMu_mass);
    newEvent->SetMetNoMuEta(*_metNoMu_eta);
    newEvent->SetMetNoMuPhi(*_metNoMu_phi);
    
    events.push_back(newEvent);
  }
}

Events::~Events()
{
  
}

Events* Events::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  Events *outputEvents = new Events(*this);
  
  if(trackCut)  outputEvents = outputEvents->ApplyTrackCut(trackCut);
  if(jetCut)    outputEvents = outputEvents->ApplyJetCut(jetCut);
  if(leptonCut) outputEvents = outputEvents->ApplyLeptonCut(leptonCut);
  if(eventCut)  outputEvents = outputEvents->ApplyEventCut(eventCut);
  
  return outputEvents;
}

Events* Events::ApplyEventCut(EventCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    if(events[iEvent]->IsPassingCut(cut)){
      outputEvents->AddEvent(events[iEvent]);
    }
  }
  return outputEvents;
}

Events* Events::ApplyTrackCut(TrackCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyTrackCut(cut));
  }
  return outputEvents;
}

Events* Events::ApplyJetCut(JetCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyJetCut(cut));
  }
  return outputEvents;
}

Events* Events::ApplyLeptonCut(LeptonCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyLeptonCut(cut));
  }
  return outputEvents;
}

double Events::WeightedSize(){
  if(events.size()==0) return 0;
  return events[0]->GetWeight()*events.size();
}

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

Event::Event()
{
  
}

Event::~Event()
{
  
}

void Event::Print(){
  for(auto t : tracks){ t->Print(); }
  for(auto j : jets){   j->Print(); }
}

Event* Event::CopyThisEventProperties()
{
  Event *outputEvent = new Event();
    
  outputEvent->SetWeight(weight);
  outputEvent->SetNvertices(nVertices);
  outputEvent->SetNjet30(nJet30);
  outputEvent->SetNjet30a(nJet30a);
  outputEvent->SetNlepton(nLepton);
  outputEvent->SetNtau(nTau);
  
  outputEvent->SetMetSumEt(metSumEt);
  outputEvent->SetMetPt(metPt);
  outputEvent->SetMetMass(metMass);
  outputEvent->SetMetPhi(metPhi);
  outputEvent->SetMetEta(metEta);
  
  outputEvent->SetMetNoMuPt(metNoMuPt);
  outputEvent->SetMetNoMuMass(metNoMuMass);
  outputEvent->SetMetNoMuPhi(metNoMuPhi);
  outputEvent->SetMetNoMuEta(metNoMuEta);
  outputEvent->SetHasNoMuTrigger(metNoMuTrigger);
  
  return outputEvent;
}

Event* Event::ApplyTrackCut(TrackCut *cut)
{
  Event *outputEvent = CopyThisEventProperties();
  
  for(auto j : jets){outputEvent->AddJet(j);}
  for(auto l : leptons){outputEvent->AddLepton(l);}
  
  vector<Track*> tracksPassingCut;
  
  for(auto track : tracks){
    if(track->IsPassingCut(cut)){
      outputEvent->AddTrack(track);
    }
  }
  return outputEvent;
}

Event* Event::ApplyJetCut(JetCut *cut)
{
  Event *outputEvent = CopyThisEventProperties();
  
  for(auto t : tracks){outputEvent->AddTrack(t);}
  for(auto l : leptons){outputEvent->AddLepton(l);}

  vector<Track*> jetPassingCuts;
  
  for(auto jet : jets){
    if(jet->IsPassingCut(cut)){
      outputEvent->AddJet(jet);
    }
  }
  return outputEvent;
}

Event* Event::ApplyLeptonCut(LeptonCut *cut)
{
  Event *outputEvent = CopyThisEventProperties();
  
  for(auto t : tracks){outputEvent->AddTrack(t);}
  for(auto j : jets){outputEvent->AddJet(j);}
  
  vector<Track*> leptonsPassingCuts;
  
  for(auto lepton : leptons){
    if(lepton->IsPassingCut(cut)){
      outputEvent->AddLepton(lepton);
    }
  }
  return outputEvent;
}

bool Event::IsPassingCut(EventCut *cut)
{
  // check the trigger
  if(cut->RequiresMetNoMuTrigger() && !metNoMuTrigger)  return false;
  
  // check number of objects
  if(nLepton < cut->GetMinNleptons() || nLepton > cut->GetMaxNleptons())  return false;
  if(nTau > cut->GetMaxNtau()) return false;
  
  vector<Lepton*> muons;
  for(auto l : leptons){
    if(abs(l->GetPid()) == 13) muons.push_back(l);
  }
  
  // check number of muons
  if(muons.size() < cut->GetMinNmuons() || muons.size() > cut->GetMaxNmuons()){
    return false;
  }
  
  // make sure they have an opposite sign
  if(cut->RequiresTwoOppositeMuons()){
    if(muons.size() != 2) return false;
    if(muons[0]->GetPid() != -muons[1]->GetPid()) return false;
  }
  
  // apply tight muon cuts (tightID flag, pt > 20 GeV, isolation < 0.15
  if(cut->RequiresTightMuon()){
    unsigned int leptonCutOptions = LeptonCut::kIsolated | LeptonCut::kTightID | LeptonCut::kPt20GeV;
    LeptonCut *tightMuonCut = new LeptonCut((LeptonCut::ECut)leptonCutOptions);

    bool atLeastOneTightMuon = false;
    if(muons[0]->IsPassingCut(tightMuonCut)) atLeastOneTightMuon = true;
    if(muons[1]->IsPassingCut(tightMuonCut)) atLeastOneTightMuon = true;
    if(!atLeastOneTightMuon)  return false;
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
  
  if(metPt < cut->GetMinMetPt())  return false;
  if(metNoMuPt < cut->GetMinMetNoMuPt())  return false;
  
  // Remove jets that are too close to muons (they will be permanently removed from the event)
  if(cut->RequiresMuJetR0p4()){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), 0.1057);
      muonVectors.push_back(mv);
    }
    
    for(int iJet=0;iJet<GetNcentralJets();iJet++){
      Jet *j = GetJet(iJet);
      if(!j->IsForward()){
        
        TLorentzVector jetVector;
        jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
        
        for(auto muonVector : muonVectors){
          if(jetVector.DeltaR(muonVector) < 0.4){
            jets.erase(jets.begin()+iJet);
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
    
    for(int iTrack=0;iTrack<GetNtracks();iTrack++){
      Track *t = tracks[iTrack];
      
      TLorentzVector trackVector;
      trackVector.SetPtEtaPhiM(t->GetPt(), t->GetEta(), t->GetPhi(), t->GetMass());
      
      for(auto muonVector : muonVectors){
        if(trackVector.DeltaR(muonVector) < 0.4){
          tracks.erase(tracks.begin()+iTrack);
          iTrack--;
          break;
        }
      }
    }
  }
  
  // check number of tracks and jets after removing those that are too close to muons
  if(GetNcentralJets() < cut->GetMinNjets()) return false;
  if(GetNtracks() < cut->GetMinNtracks()) return false;
  

  // find the jet with the highest pt
  Jet *highJet = nullptr;
  double highestPt = -1.0;

  for(int iJet=0;iJet<GetNjets();iJet++){
    if(jets[iJet]->GetPt() > highestPt){
      highestPt = jets[iJet]->GetPt();
      highJet = jets[iJet];
    }
  }

  // check properties of the highest pt jet
  if(cut->RequiresHighJet() && !highJet) return false;
  if(highJet->GetPt() < cut->GetHighJetMinPt()) return false;
  if(fabs(highJet->GetEta()) > cut->GetHighJetMaxEta()) return false;
  if(highJet->GetChargedHadronEnergyFraction() < cut->GetHighJetMinChHEF()) return false;
  if(highJet->GetNeutralHadronEnergyFraction() > cut->GetHighJetMaxNeHEF()) return false;


  
  if(cut->RequiresMetJetPhi0p5()){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(metPt, metEta, metPhi, metMass);

    for(auto j : jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(fabs(metVector.DeltaPhi(jetVector)) < 0.5) return false;
    }
  }

  if(cut->RequiresMetNoMuJetPhi0p5()){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(metNoMuPt, metNoMuEta, metNoMuPhi, metNoMuMass);

    for(auto j : jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(fabs(metVector.DeltaPhi(jetVector)) < 0.5) return false;
    }
  }
  

  
  return true;
}






