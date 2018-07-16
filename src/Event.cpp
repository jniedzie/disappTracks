//
//  Event.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Event.hpp"

Events::Events(string fileName)
{
  
  map<unsigned long long,Event*> eventsMap =  GetEventsFromFile(fileName);
  
  for(auto ev : eventsMap){
    events.push_back(ev.second);
  }
}

Events::~Events()
{
  
}

Event* Events::operator[](const int index)
{
  return events[index];
}

map<unsigned long long,Event*> Events::GetEventsFromFile(string fileName)
{
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
  
  map<unsigned long long,Event*> events;
  
  while (reader.Next()){
    Track *track = new Track();
    
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      track->SetDeDxInLayer(iLayer, **dedx[iLayer]);
      track->SetSubDetIdInLayer(iLayer, **subDetId[iLayer]);
    }
    track->SetEta(*_eta);
    track->SetPhi(*_phi);
    
    int nDeDxPoints = 0;
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      if(track->GetDeDxInLayer(iLayer) > 0.000001) nDeDxPoints++;
    }
    if(nDeDxPoints <= shortTrackMaxNclusters) track->SetIsShort(true);
    
    if(events.find(*eventNumber) == events.end()){
      Event *newEvent = new Event();
      newEvent->AddTrack(track);
      events.insert(pair<unsigned long long,Event*>(*eventNumber,newEvent));
    }
    else{
      events[*eventNumber]->AddTrack(track);
    }
  }
  return events;
}

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

void Event::Print(){
  for(auto t : tracks){
    t->Print();
  }
}

Event* Event::FilterShortTracksAboveThreshold(double threshold)
{
  Event *outputEvent = new Event();
  
  for(int iTrack=0;iTrack<GetNtracks();iTrack++){
    Track *track = GetTrack(iTrack);
    if(!track->GetIsShort()) continue;
    
    double totalDeDx = 0;
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      totalDeDx += track->GetDeDxInLayer(iLayer);
    }
    if(totalDeDx > threshold) outputEvent->AddTrack(track);
  }
  return outputEvent;
}

Event* Event::FilterShortTracks()
{
  Event *outputEvent = new Event();
  
  for(int iTrack=0;iTrack<GetNtracks();iTrack++){
    Track *track = GetTrack(iTrack);
    if(track->GetIsShort()) outputEvent->AddTrack(track);
  }
  return outputEvent;
  }






