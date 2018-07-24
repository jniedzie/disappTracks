//
//  Event.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Event.hpp"

#include <TTreeReaderArray.h>

Events::Events()
{
  
}

Events::Events(string fileName)
{
  cout<<"Reading events from:"<<fileName<<endl;
  
  map<unsigned long long,Event*> eventsMap;
  
  TFile *inFile = TFile::Open(fileName.c_str());
  TTreeReader reader("tree", inFile);
  
  TTreeReaderValue<unsigned long long> eventNumber(reader, "evt");
  TTreeReaderValue<int> nJets(reader, "nJet");
  TTreeReaderValue<int> nTracks(reader, "nIsoTrack");
  
  TTreeReaderArray<float> *dedx[nLayers];
  TTreeReaderArray<int> *subDetId[nLayers];
  TTreeReaderArray<int> *sizeX[nLayers];
  TTreeReaderArray<int> *sizeY[nLayers];
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx[iLayer] =  new TTreeReaderArray<float>(reader,Form("IsoTrack_dedxByLayer%i",iLayer));
    subDetId[iLayer] =  new TTreeReaderArray<int>(reader,Form("IsoTrack_subDetIdByLayer%i",iLayer));
    sizeX[iLayer] =  new TTreeReaderArray<int>(reader,Form("IsoTrack_sizeXbyLayer%i",iLayer));
    sizeY[iLayer] =  new TTreeReaderArray<int>(reader,Form("IsoTrack_sizeYbyLayer%i",iLayer));
  }
  
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
  
  TTreeReaderArray<float> _jet_pt(reader,  "Jet_pt");
  TTreeReaderArray<float> _jet_eta(reader, "Jet_eta");
  TTreeReaderArray<float> _jet_phi(reader, "Jet_phi");
  
  TTreeReaderValue<int>   _nVert(reader, "nVert");
  TTreeReaderValue<int>   _nJet30(reader, "nJet30");
  TTreeReaderValue<int>   _nJet30a(reader, "nJet30a");
  TTreeReaderValue<float> _met_sumEt(reader, "met_sumEt");
  TTreeReaderValue<float> _met_pt(reader, "met_pt");
  TTreeReaderValue<float> _met_mass(reader, "met_mass");
  TTreeReaderValue<float> _met_phi(reader, "met_phi");
  TTreeReaderValue<float> _met_eta(reader, "met_eta");

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
    
    newEvent->SetNvertices(*_nVert);
    newEvent->SetNjetsFromTree(*nJets);
    newEvent->SetNisoTracksFromTree(*nTracks);
    newEvent->SetNjet30(*_nJet30);
    newEvent->SetNjet30a(*_nJet30a);
    newEvent->SetMetSumEt(*_met_sumEt);
    newEvent->SetMetPt(*_met_pt);
    newEvent->SetMetMass(*_met_mass);
    newEvent->SetMetEta(*_met_eta);
    newEvent->SetMetPhi(*_met_phi);
  
    events.push_back(newEvent);
    



//

  }
}

Events::~Events()
{
  
}

int Events::SizeNonEmpty()
{
  int n=0;
  for(auto event : events){
    if(event->GetNtracks() > 0 && event->GetNjets() > 0) n++;
  }
  return n;
}

int Events::GetNtracks()
{
  int nTracks = 0;
  for(int iEvent=0;iEvent < events.size();iEvent++){
    nTracks+= events[iEvent]->GetNtracks();
  }
  return nTracks;
}

int Events::GetNjets()
{
  int nJets = 0;
  for(int iEvent=0;iEvent < events.size();iEvent++){
    nJets += events[iEvent]->GetNjets();
  }
  return nJets;
}

Events* Events::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut)
{
  Events *outputEvents = new Events();
  if(eventCut){
    outputEvents = ApplyEventCut(eventCut);
    if(trackCut)  outputEvents = outputEvents->ApplyTrackCut(trackCut);
    if(jetCut)    outputEvents = outputEvents->ApplyJetCut(jetCut);
  }
  else{
    if(trackCut)  outputEvents = ApplyTrackCut(trackCut);
    if(jetCut)    outputEvents = outputEvents->ApplyJetCut(jetCut);
  }
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

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

void Event::Print(){
  for(auto t : tracks){ t->Print(); }
  for(auto j : jets){   j->Print(); }
}

Event* Event::ApplyTrackCut(TrackCut *cut)
{
  Event *outputEvent = new Event();
  for(auto j : jets){outputEvent->AddJet(j);}
  
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
  if(metPt < cut->GetMinMetPt()){
    return false;
  }
  
  return true;
}

