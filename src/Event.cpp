//
//  Event.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Event.hpp"

Events::Events()
{
  
}

Events::Events(string fileName)
{
  cout<<"Reading events from:"<<fileName<<endl;
  
  map<unsigned long long,Event*> eventsMap;
  
  TFile *inFile = TFile::Open(fileName.c_str());
  TTreeReader reader("tree", inFile);
  
  //  inFile->Get("tree")->Print();
  
  TTreeReaderValue<unsigned long long> eventNumber(reader, "evt");
  TTreeReaderValue<float> *dedx[nLayers];
  TTreeReaderValue<int> *subDetId[nLayers];
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx[iLayer] =  new TTreeReaderValue<float>(reader,Form("IsoTrack_dedxByLayer%i",iLayer));
    subDetId[iLayer] =  new TTreeReaderValue<int>(reader,Form("IsoTrack_subDetIdByLayer%i",iLayer));
  }
  
  TTreeReaderValue<float> _eta(reader, "IsoTrack_eta");
  TTreeReaderValue<float> _phi(reader, "IsoTrack_phi");
  TTreeReaderValue<float> _caloEmEnergy(reader, "IsoTrack_caloEmEnergy");
  TTreeReaderValue<float> _caloHadEnergy(reader, "IsoTrack_caloHadEnergy");
  TTreeReaderValue<float> _dxyErr(reader, "IsoTrack_edxy");
  TTreeReaderValue<float> _dxy(reader, "IsoTrack_dxy");
  TTreeReaderValue<float> _dzErr(reader, "IsoTrack_edz");
  TTreeReaderValue<float> _dz(reader, "IsoTrack_dz");
  TTreeReaderValue<int>   _charge(reader, "IsoTrack_charge");
  TTreeReaderValue<float> _mass(reader, "IsoTrack_mass");
  TTreeReaderValue<float> _pt(reader, "IsoTrack_pt");
  TTreeReaderValue<int>   _pid(reader, "IsoTrack_pdgId");
  
  while (reader.Next()){
    Track *track = new Track();
    
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      track->SetDeDxInLayer(iLayer, **dedx[iLayer]);
      track->SetSubDetIdInLayer(iLayer, **subDetId[iLayer]);
    }
    track->SetEta(*_eta);
    track->SetPhi(*_phi);
    track->SetCaloEmEnergy(*_caloEmEnergy);
    track->SetCaloHadEnergy(*_caloHadEnergy);
    track->SetDxy(*_dxy,*_dxyErr);
    track->SetDz(*_dz,*_dzErr);
    track->SetCharge(*_charge);
    track->SetMass(*_mass);
    track->SetPt(*_pt);
    track->SetPid(*_pid);
    
    int nDeDxPoints = 0;
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      if(track->GetDeDxInLayer(iLayer) > 0.000001) nDeDxPoints++;
    }
    if(nDeDxPoints <= shortTrackMaxNclusters) track->SetIsShort(true);
    
    if(eventsMap.find(*eventNumber) == eventsMap.end()){
      Event *newEvent = new Event();
      newEvent->AddTrack(track);
      eventsMap.insert(pair<unsigned long long,Event*>(*eventNumber,newEvent));
    }
    else{
      eventsMap[*eventNumber]->AddTrack(track);
    }
  }
  
  for(auto ev : eventsMap){
    events.push_back(ev.second);
  }
}

Events::~Events()
{
  
}

int Events::GetNtracks()
{
  int nTracks = 0;
  
  for(int iEvent=0;iEvent < events.size();iEvent++){
    nTracks+= events[iEvent]->GetNtracks();
  }
  return nTracks;
}

Events* Events::ApplyTrackCut(TrackCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyTrackCut(cut));
  }
  return outputEvents;
}

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

void Event::Print(){
  for(auto t : tracks){
    t->Print();
  }
}

Event* Event::ApplyTrackCut(TrackCut *cut)
{
  Event *outputEvent = new Event();
  
  vector<Track*> tracksPassingCut;
  
  for(auto track : tracks){
    if(track->IsPassingCut(cut)){
      outputEvent->AddTrack(track);
    }
  }
  return outputEvent;
}

vector<Track*> Event::GetTracksPassingCut(TrackCut *cut)
{
  vector<Track*> tracksPassingCut;
  
  for(auto track : tracks){
    if(track->IsPassingCut(cut)){
      tracksPassingCut.push_back(track);
    }
  }
  return tracksPassingCut;
}





