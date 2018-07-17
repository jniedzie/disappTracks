#include "Helpers.hpp"
#include "Event.hpp"

#include <TSystem.h>
#include <TEveManager.h>
#include <TEveScene.h>
#include <TEvePointSet.h>
#include <TApplication.h>

const bool showUnderflowBins = false;
const bool showOverflowBins = true;

const int hitMarkerStyle = 20;
const int hitMarkerSize = 2.0;

const double dedxBinsMin = 1;
const double dedxBinsMax = 10;
const int dedxNbins = 5;

                                  //     Underflow                                       Overflow
const int dedxBinColors[dedxNbins+2] = { kGray,     kBlue, kCyan, kGreen, kYellow, kRed, kMagenta };

TEvePointSetArray* PreparePointsEventDisplay()
{
  TEvePointSetArray *points = new TEvePointSetArray("Hits");
  points->SetMarkerStyle(hitMarkerStyle);
  points->SetMarkerSize(hitMarkerSize);
  
  points->InitBins("Cluster Charge", dedxNbins, dedxBinsMin,dedxBinsMax);

  for(int i=0;i<(dedxNbins+2);i++){
    points->GetBin(i)->SetMainColor(dedxBinColors[i]);
  }
  
  points->GetBin(0)->SetRnrSelf(showUnderflowBins);
  points->GetBin(dedxNbins+1)->SetRnrSelf(showOverflowBins);
  
  return points;
}

void DrawPoints(TEvePointSetArray *points)
{
  points->SetRnrSelf(kTRUE);
  gEve->AddElement(points);
  gEve->Redraw3D();
}

void DrawEvent(Event *event, bool cleanView=true)
{
  gEve->GetEventScene()->DestroyElements();
  gSystem->ProcessEvents();
  
  TEvePointSetArray *points = PreparePointsEventDisplay();
  
  for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
    Track *track = event->GetTrack(iTrack);
    
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      double R = layerR[iLayer];
      double phi = track->GetPhi();
      double theta = 2*atan(exp(-track->GetEta()));
      
      double x = R*sin(theta)*cos(phi);
      double y = R*sin(theta)*sin(phi);
      double z = R*cos(theta);
      
      //          cout<<"Point R:"<<R<<"\teta:"<<eta<<"\tphi:"<<phi<<"\ttheta:"<<theta<<endl;
  //      cout<<"Point x:"<<x<<"\ty:"<<y<<"\tz:"<<z<<"\tvalue:"<<track->dedx[iLayer]<<endl;
      
      points->Fill(x,y,z,track->GetDeDxInLayer(iLayer));
    }
  }
  DrawPoints(points);
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  TEveManager::Create();
  
//  const char *inFileNameSignal = "../adish/Signal/tree.root";
//  const char *inFileNameBackground = "../adish/Background/tree.root";
  
    const char *inFileNameSignal = "../jniedzie/mcSignal/tree.root";
//    const char *inFileNameBackground = "../jniedzie/mcBackground/tree.root";
  
  Events *eventsSignal = new Events(inFileNameSignal);
//  Events *eventsBackground = new Events(inFileNameBackground);
  
  
  TrackCut *trackCut = new TrackCut(TrackCut::kShortAboveThreshold);
  
  Events *filteredSignalEvents = eventsSignal->ApplyTrackCut(trackCut);
  
  for(int iEvent=0;iEvent<filteredSignalEvents->size();iEvent++){
    Event* event = filteredSignalEvents->At(iEvent);
    
    if(event->GetNtracks() < 1) continue;

    cout<<"Event iter:"<<iEvent<<endl;
    DrawEvent(event);
    break;
  }
  
  theApp.Run();
  
  return 0;
}
