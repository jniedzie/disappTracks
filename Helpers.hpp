//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <vector>
#include <iostream>
#include <map>

using namespace std;

const int nLayers = 14;
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };

const int shortTrackMaxNclusters = 3;

class Track{
public:
  Track(){
    isShort = false;
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      dedx.push_back(0.0);
      subDetId.push_back(-1);
    }
  };
  ~Track(){};
  
  void SetDeDxInLayer(int layer, float value){dedx[layer] = value;}
  void SetSubDetIdInLayer(int layer, int id){subDetId[layer] = id;}
  
  void SetEta(double _eta){eta=_eta;}
  void SetPhi(double _phi){phi=_phi;}
  
  void SetIsShort(bool _isShort){isShort = _isShort;}

  float GetDeDxInLayer(int layer){return dedx[layer];}
  int GetSubDetIdInLayer(int layer){return subDetId[layer];}
  
  double GetEta(){return eta;}
  double GetPhi(){return phi;}
  
  bool GetIsShort(){return isShort;}
  
  float GetTotalDedx(){return accumulate(dedx.begin(),dedx.end(),0.0);}
  
  int GetNclusters(){
    int nClusters=0;
    for(float d : dedx){
      if(d > 0.000001) nClusters++;
    }
    return nClusters;
  }
  
  void Print(){
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      cout<<"Layer:"<<iLayer<<"\tsub-det ID:"<<subDetId[iLayer]<<"\tdEdx:"<<dedx[iLayer]<<endl;
    }
  }
private:
  vector<float> dedx;   // dedx in consecutive layers
  vector<int> subDetId; // sub-detector IDs for each layer
  double eta;
  double phi;
  bool isShort; // track is short if it has max 3 dedx points
};

class Event{
public:
  Event(){};
  ~Event(){};
  
  void AddTrack(Track *track){tracks.push_back(track);}
  int GetNtracks(){return tracks.size(); }
  Track* GetTrack(int i){return tracks[i];}
  
  void Print(){
    for(auto t : tracks){
      t->Print();
    }
  }
private:
  vector<Track*> tracks; // vector of isolated tracks
  
};

Event* FilterShortTracksAboveThreshold(Event *inputEvent, double threshold)
{
  Event *outputEvent = new Event();
  
  for(int iTrack=0;iTrack<inputEvent->GetNtracks();iTrack++){
    Track *track = inputEvent->GetTrack(iTrack);
    if(!track->GetIsShort()) continue;
    
    double totalDeDx = 0;
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      totalDeDx += track->GetDeDxInLayer(iLayer);
    }
    if(totalDeDx > threshold) outputEvent->AddTrack(track);
  }
  return outputEvent;
}

Event* FilterShortTracks(Event *inputEvent)
{
  Event *outputEvent = new Event();
  
  for(int iTrack=0;iTrack<inputEvent->GetNtracks();iTrack++){
    Track *track = inputEvent->GetTrack(iTrack);
    if(track->GetIsShort()) outputEvent->AddTrack(track);
  }
  return outputEvent;
}

map<unsigned long long,Event*> GetEventsFromFile(const char *fileName)
{
  TFile *inFile = TFile::Open(fileName);
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

vector<Event*> GetEventsVectorFromFile(const char *fileName)
{
  vector<Event*> output;
  map<unsigned long long,Event*> events =  GetEventsFromFile(fileName);
  
  for(auto ev : events){
    output.push_back(ev.second);
  }
  return output;
}

TLegend* GetLegend(double legendW = 0.15, double legendH = 0.5, double legendX = 0.75, double legendY = 0.25,const char* header="")
{
  
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader(header);
  return leg;
}

#endif /* Helpers_h */
