//
//  Event.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Event.hpp"

#include <TTreeReaderArray.h>
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

  TTreeReaderValue<int>   nJets(reader, "nJet");
  TTreeReaderValue<int>   nTracks(reader, "nIsoTrack");
  TTreeReaderValue<int>   _nVert(reader, "nVert");
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
  TTreeReaderArray<float> _lepton_thight_pid(reader, "LepGood_tightId");
  
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

  TTreeReaderArray<float> _jet_pt(reader,  "Jet_pt");
  TTreeReaderArray<float> _jet_eta(reader, "Jet_eta");
  TTreeReaderArray<float> _jet_phi(reader, "Jet_phi");

  while (reader.Next()){
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
      newEvent->AddJet(jet);
    }

    for(int iLepton=0;iLepton<*_nLepton;iLepton++){
      Lepton *lepton = new Lepton();
      lepton->SetPt(_lepton_pt[iLepton]);
      lepton->SetEta(_lepton_eta[iLepton]);
      lepton->SetPhi(_lepton_phi[iLepton]);
      lepton->SetPid(_lepton_thight_pid[iLepton]);
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

Events* Events::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut)
{
  Events *outputEvents = new Events(*this);
  
  if(trackCut)  outputEvents = outputEvents->ApplyTrackCut(trackCut);
  if(jetCut)    outputEvents = outputEvents->ApplyJetCut(jetCut);
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

Event* Event::ApplyTrackCut(TrackCut *cut)
{
  Event *outputEvent = new Event();
  for(auto j : jets){outputEvent->AddJet(j);}
  
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
  Event *outputEvent = new Event();
  for(auto t : tracks){outputEvent->AddTrack(t);}
  
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
  
  vector<Track*> jetPassingCuts;
  
  for(auto jet : jets){
    if(jet->IsPassingCut(cut)){
      outputEvent->AddJet(jet);
    }
  }
  return outputEvent;
}

bool Event::IsPassingCut(EventCut *cut)
{
  // check MET properties
  if(metPt < cut->GetMinMetPt()) return false;
  if(metNoMuPt < cut->GetMinMetNoMuPt()) return false;
  if(cut->RequiresMetNoMuTrigger() && !metNoMuTrigger) return false;
  
  // check number of objects
  if(GetNjets() < cut->GetMinNjets()) return false;
  if(GetNtracks() < cut->GetMinNtracks()) return false;
  if(nLepton > cut->GetMaxNlepton()) return false;
  if(nTau > cut->GetMaxNtau()) return false;
  
  return true;
}






