//
//  Display.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Display.hpp"

Display::Display()
{
  TEveManager::Create();
}

Display::~Display()
{
  
}

void Display::DrawSimplePoints(const shared_ptr<vector<Point>> points, const map<string,any> options)
{
  TEvePointSetArray *simplePoints = PreparePointsEventDisplay(options);
  
  for(auto p : *points){
    simplePoints->Fill(scale*p.GetX(),scale*p.GetY(),scale*p.GetZ(), p.GetValue());
  }
  
  simplePoints->SetRnrSelf(kTRUE);
  gEve->AddElement(simplePoints);
  gEve->Redraw3D();
}

void Display::DrawHelix(const unique_ptr<Helix> &helix, const map<string,any> options)
{
  TEvePointSetArray *helixPoints = PreparePointsEventDisplay(options);
  
  double tMin  = helix->GetTmin();
  double tMax  = helix->GetTmax();
  double tStep = helix->GetTstep();
  
  int zSign = sgn(helix->GetMomentum()->GetZ());
  
  auto fillPointForT = [&](double t){
    double x = helix->GetOrigin()->GetX();
    double y = helix->GetOrigin()->GetY();
    double z = helix->GetOrigin()->GetZ() + fabs(helix->GetSlope())*t;
    
    if(helix->GetCharge() > 0){
      x += helix->GetRadius()*cos(t);
      y += helix->GetRadius()*sin(t);
    }
    else{
      x += helix->GetRadius()*sin(t);
      y += helix->GetRadius()*cos(t);
    }
    helixPoints->Fill(scale*x,scale*y,scale*z, 0);
  };
  
  if(zSign < 0) for(double t = tMin; t > -tMax; t -= tStep){fillPointForT(t);}
  else          for(double t = tMin; t <  tMax; t += tStep){fillPointForT(t);}
  
  helixPoints->SetRnrSelf(kTRUE);
  gEve->AddElement(helixPoints);
  gEve->Redraw3D();
}

void Display::DrawEvent(const shared_ptr<Event> &event, const map<string,any> options)
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
  
  TEvePointSetArray *points = PreparePointsEventDisplay(options);
  
  for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
    shared_ptr<Track> track = event->GetTrack(iTrack);
    
    for(int iLayer=0;iLayer<nLayers;iLayer++){
      double R = layerR[iLayer];
      double phi = track->GetPhi();
      double theta = 2*atan(exp(-track->GetEta()));
      
      double x = R*sin(theta)*cos(phi);
      double y = R*sin(theta)*sin(phi);
      double z = R*cos(theta);
      
      points->Fill(scale*x,scale*y,scale*z,track->GetDeDxInLayer(iLayer));
    }
  }
  points->SetRnrSelf(kTRUE);
  gEve->AddElement(points);
  gEve->Redraw3D();
  
  // MET
  double metPhi = event->GetMetPhi();
  double metTheta = 2*atan(exp(-event->GetMetEta()));
  DrawMET(metPhi, metTheta);
  
  // Geometry:
  if(!showGeometry) return;
  
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

TEvePointSetArray* Display::PreparePointsEventDisplay(map<string,any> options)
{
  const char* title = any_cast<const char*>(options["title"]);
  TEvePointSetArray *points = new TEvePointSetArray(title);
  points->SetMarkerStyle(any_cast<int>(options["markerStyle"]));
  points->SetMarkerSize(any_cast<double>(options["markerSize"]));
  
  int nBins = options.count("nBins")   ? any_cast<int>(options["nBins"]) : 1;
  
  points->InitBins(title, nBins,
                   options.count("binsMin") ? any_cast<int>(options["binsMin"]) : -inf,
                   options.count("binsMax") ? any_cast<int>(options["binsMax"]) :  inf);
  
  for(int i=0;i<(nBins+2);i++){
    points->GetBin(i)->SetMainColor(options.count("color") ? any_cast<EColor>(options["color"]) : dedxBinColors[i]);
  }
  
  points->GetBin(0)->SetRnrSelf(showUnderflowBins);
  points->GetBin(nBins+1)->SetRnrSelf(showOverflowBins);
  
  return points;
}

void Display::DrawMET(double metPhi, double metTheta)
{
  TEveBox *metBox = new TEveBox("MET");
  
  vector<int> a = {-1, 1, 1,-1,-1, 1, 1,-1};
  vector<int> b = { 1, 1,-1,-1, 1, 1,-1,-1};
  vector<int> c = { 1, 1, 1, 1,-1,-1,-1,-1};
  
  // iterate over vertices of the box
  for(int i=0;i<8;i++){
    // calculate proper size of the box
    double R      = metRadius+a[i]*metBoxSize;
    double theta  = metTheta+b[i]*metBoxAngularSize;
    double phi    = metPhi+c[i]*metBoxAngularSize;
    // convert to XYZ coordinates
    metBox->SetVertex(i,R*sin(theta)*cos(phi),R*sin(theta)*sin(phi),R*cos(theta));
  }
  
  metBox->SetMainColorRGB((Float_t)0.0, 1.0, 1.0);
  metBox->SetRnrSelf(true);
  
  gEve->AddElement(metBox);
  gEve->Redraw3D();
}
