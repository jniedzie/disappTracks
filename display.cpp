#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"

#include <TSystem.h>
#include <TEveManager.h>
#include <TEveScene.h>
#include <TEvePointSet.h>
#include <TEveJetCone.h>
#include <TEveBox.h>
#include <TApplication.h>
#include <TGeoShape.h>
#include <TGeoTube.h>
#include <TEveGeoShape.h>

const double scale = 0.1;

const bool showUnderflowBins = false;
const bool showOverflowBins = true;

const int hitMarkerStyle = 20;
const int hitMarkerSize = 2.0;

const double dedxBinsMin = 1;
const double dedxBinsMax = 10;
const int dedxNbins = 5;

const double jetConeRadius = scale*0.4;

const double metRadius = scale * 2000;
const double metBoxSize = scale * 30;
const double metBoxAngularSize = 0.1;

const int geomTransparency = 90; // 30 - 100

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

void DrawEvent(shared_ptr<Event> event)
{
  gEve->GetEventScene()->DestroyElements();
  gSystem->ProcessEvents();
 
  for(int iJet=0;iJet<event->GetNjets();iJet++){
    shared_ptr<Jet> jet = event->GetJet(iJet);
    TEveJetCone *jetCone = new TEveJetCone();
    jetCone->SetCylinder(scale*2900, scale*5500);
    jetCone->AddCone(jet->GetEta(), jet->GetPhi(), jetConeRadius);
    jetCone->SetMainColorRGB((Float_t)1.0, 0.0, 0.0);
    jetCone->SetRnrSelf(kTRUE);
    gEve->AddElement(jetCone);
    gEve->Redraw3D();
  }
  
  
  TEvePointSetArray *points = PreparePointsEventDisplay();
  
  for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
    shared_ptr<Track> track = event->GetTrack(iTrack);
    
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      double R = layerR[iLayer];
      double phi = track->GetPhi();
      double theta = 2*atan(exp(-track->GetEta()));
      
      double x = R*sin(theta)*cos(phi);
      double y = R*sin(theta)*sin(phi);
      double z = R*cos(theta);
      
      //          cout<<"Point R:"<<R<<"\teta:"<<eta<<"\tphi:"<<phi<<"\ttheta:"<<theta<<endl;
  //      cout<<"Point x:"<<x<<"\ty:"<<y<<"\tz:"<<z<<"\tvalue:"<<track->dedx[iLayer]<<endl;
      
      points->Fill(scale*x,scale*y,scale*z,track->GetDeDxInLayer(iLayer));
    }
  }
  DrawPoints(points);
  
  
  // MET
  
  TEveBox *metBox = new TEveBox("MET");
  
  double metPhi = event->GetMetPhi();
  double metTheta = 2*atan(exp(-event->GetMetEta()));
  
//  double metX = metRadius*sin(metTheta)*cos(metPhi);
//  double metY = metRadius*sin(metTheta)*sin(metPhi);
//  double metZ = metRadius*cos(metTheta);
//
//  metBox->SetVertex(0, metX-metBoxSize, metY+metBoxSize, metZ+metBoxSize);
//  metBox->SetVertex(1, metX+metBoxSize, metY+metBoxSize, metZ+metBoxSize);
//  metBox->SetVertex(2, metX+metBoxSize, metY-metBoxSize, metZ+metBoxSize);
//  metBox->SetVertex(3, metX-metBoxSize, metY-metBoxSize, metZ+metBoxSize);
//
//  metBox->SetVertex(4, metX-metBoxSize, metY+metBoxSize, metZ-metBoxSize);
//  metBox->SetVertex(5, metX+metBoxSize, metY+metBoxSize, metZ-metBoxSize);
//  metBox->SetVertex(6, metX+metBoxSize, metY-metBoxSize, metZ-metBoxSize);
//  metBox->SetVertex(7, metX-metBoxSize, metY-metBoxSize, metZ-metBoxSize);
  
  
  metBox->SetVertex(0,
                    (metRadius-metBoxSize)*sin(metTheta+metBoxAngularSize)*cos(metPhi+metBoxAngularSize),
                    (metRadius-metBoxSize)*sin(metTheta+metBoxAngularSize)*sin(metPhi+metBoxAngularSize),
                    (metRadius-metBoxSize)*cos(metTheta+metBoxAngularSize));
  
  metBox->SetVertex(1,
                    (metRadius+metBoxSize)*sin(metTheta+metBoxAngularSize)*cos(metPhi+metBoxAngularSize),
                    (metRadius+metBoxSize)*sin(metTheta+metBoxAngularSize)*sin(metPhi+metBoxAngularSize),
                    (metRadius+metBoxSize)*cos(metTheta+metBoxAngularSize));
  
  metBox->SetVertex(2,
                    (metRadius+metBoxSize)*sin(metTheta-metBoxAngularSize)*cos(metPhi+metBoxAngularSize),
                    (metRadius+metBoxSize)*sin(metTheta-metBoxAngularSize)*sin(metPhi+metBoxAngularSize),
                    (metRadius+metBoxSize)*cos(metTheta-metBoxAngularSize));
  
  metBox->SetVertex(3,
                    (metRadius-metBoxSize)*sin(metTheta-metBoxAngularSize)*cos(metPhi+metBoxAngularSize),
                    (metRadius-metBoxSize)*sin(metTheta-metBoxAngularSize)*sin(metPhi+metBoxAngularSize),
                    (metRadius-metBoxSize)*cos(metTheta-metBoxAngularSize));

  
  metBox->SetVertex(4,
                    (metRadius-metBoxSize)*sin(metTheta+metBoxAngularSize)*cos(metPhi-metBoxAngularSize),
                    (metRadius-metBoxSize)*sin(metTheta+metBoxAngularSize)*sin(metPhi-metBoxAngularSize),
                    (metRadius-metBoxSize)*cos(metTheta+metBoxAngularSize));
  
  metBox->SetVertex(5,
                    (metRadius+metBoxSize)*sin(metTheta+metBoxAngularSize)*cos(metPhi-metBoxAngularSize),
                    (metRadius+metBoxSize)*sin(metTheta+metBoxAngularSize)*sin(metPhi-metBoxAngularSize),
                    (metRadius+metBoxSize)*cos(metTheta+metBoxAngularSize));
  
  metBox->SetVertex(6,
                    (metRadius+metBoxSize)*sin(metTheta-metBoxAngularSize)*cos(metPhi-metBoxAngularSize),
                    (metRadius+metBoxSize)*sin(metTheta-metBoxAngularSize)*sin(metPhi-metBoxAngularSize),
                    (metRadius+metBoxSize)*cos(metTheta-metBoxAngularSize));
  
  metBox->SetVertex(7,
                    (metRadius-metBoxSize)*sin(metTheta-metBoxAngularSize)*cos(metPhi-metBoxAngularSize),
                    (metRadius-metBoxSize)*sin(metTheta-metBoxAngularSize)*sin(metPhi-metBoxAngularSize),
                    (metRadius-metBoxSize)*cos(metTheta-metBoxAngularSize));
  
  
  metBox->SetMainColorRGB((Float_t)0.0, 1.0, 1.0);
  metBox->SetRnrSelf(true);
  
  gEve->AddElement(metBox);
  gEve->Redraw3D();
  
  // Geometry:
  
  TGeoTube *pixelTube = new TGeoTube(scale*0,scale*200, scale*1500);
  TEveGeoShape *pixel = new TEveGeoShape ("Pixel tracker","Pixel tracker");
  pixel->SetShape(pixelTube);
  pixel->SetMainTransparency(geomTransparency-30);
  pixel->SetMainColorRGB((Float_t)0.0, 1.0, 0.0);
  pixel->SetRnrSelf(true);
  
  TGeoTube *trackerTube = new TGeoTube(scale*230,scale*1100,scale*2800);
  TEveGeoShape *tracker = new TEveGeoShape ("Tracker","Tracker");
  tracker->SetShape(trackerTube);
  tracker->SetMainTransparency(geomTransparency-20);
  tracker->SetMainColorRGB((Float_t)1.0, 1.0, 0.0);
  tracker->SetRnrSelf(true);
  
  TGeoTube *emCalTube = new TGeoTube(scale*1100,scale*1800,scale*3700);
  TEveGeoShape *emCal = new TEveGeoShape ("EM calo","EM calo");
  emCal->SetShape(emCalTube);
  emCal->SetMainTransparency(geomTransparency-10);
  emCal->SetMainColorRGB((Float_t)0.0, 0.0, 1.0);
  emCal->SetRnrSelf(true);
  
  TGeoTube *hadCalTube = new TGeoTube(scale*1800,scale*2900,scale*5500);
  TEveGeoShape *hadCal = new TEveGeoShape ("Had calo","Had calo");
  hadCal->SetShape(hadCalTube);
  hadCal->SetMainTransparency(geomTransparency);
  hadCal->SetMainColorRGB((Float_t)1.0, 0.0, 0.5);
  hadCal->SetRnrSelf(true);
  
  
  gEve->AddElement(pixel);
  gEve->AddElement(tracker);
  gEve->AddElement(emCal);
  gEve->AddElement(hadCal);
  gEve->Redraw3D();
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  TEveManager::Create();
  
  auto events = shared_ptr<EventSet>(new EventSet());
  events->LoadEventsFromFiles("after_L2/");
  
  auto event = events->At(EventSet::kData, kElectron_Run2017B, 0);
  DrawEvent(event);
  event->Print();
  /*
  TrackCut *trackCut = new TrackCut();
  trackCut->SetNpixelLayers(range<int>(3,4));
  events->ApplyCuts(nullptr, trackCut, nullptr, nullptr);
  
  for(int iEvent=0;iEvent<events->size(EventSet::kSignal,kWino_M_300_cTau_30);iEvent++){
    shared_ptr<Event> event = events->At(EventSet::kSignal,kWino_M_300_cTau_30, iEvent);
    if(event->GetNtracks() < 1) continue;
    cout<<"Event iter:"<<iEvent<<endl;
    DrawEvent(event);
    break;
  }
  */
  
  theApp.Run();
  return 0;
}
