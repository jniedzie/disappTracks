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

void Display::DrawSimplePoints(const vector<shared_ptr<Point>> points, map<string,any> options)
{
  TEvePointSetArray *simplePoints = PreparePointsEventDisplay(options);
  auto stripClusters = new TEveElementList(any_cast<const char*>(options["title"]));
  
  for(auto &p : points){
    if(p->GetSubDetName() != "TOB" &&
       p->GetSubDetName() != "TIB"){
      simplePoints->Fill(scale*p->GetX(),scale*p->GetY(),scale*p->GetZ(), p->GetValue());
    }
    else{
      AddStripCluster(stripClusters, p, options);
    }
  }
  
  simplePoints->SetRnrSelf(kTRUE);
  
  gEve->AddElement(simplePoints);
  gEve->AddElement(stripClusters);
  gEve->Redraw3D();
}

void Display::AddStripCluster(TEveElementList *stripClusters,
                              const shared_ptr<Point> &point,
                              map<string,any> options)
{
  vector<int> a = {-1, 1, 1,-1,-1, 1, 1,-1};
  vector<int> b = { 1, 1,-1,-1, 1, 1,-1,-1};
  vector<int> c = { 1, 1, 1, 1,-1,-1,-1,-1};
  
  TEveBox *stripBox = new TEveBox("Strip");
  for(int i=0;i<8;i++){
    double x  = point->GetX() + a[i] * point->GetXerr();
    double y  = point->GetY() + b[i] * point->GetYerr();
    double z  = point->GetZ() + c[i] * point->GetZerr();
    stripBox->SetVertex(i,scale*x,scale*y,scale*z);
  }
  
  stripBox->SetMainColor(any_cast<EColor>(options["color"]));
  stripBox->SetRnrSelf(true);
  
  stripClusters->AddElement(stripBox);
}

void Display::DrawHelix(const unique_ptr<Helix> &helix, const map<string,any> options)
{
  TEvePointSetArray *helixPoints = PreparePointsEventDisplay(options);
  
  double tMin  = helix->GetTmin();
  double tMax  = helix->GetTmax();
  double tStep = helix->GetTstep();
  
  int zSign = sgn(helix->GetMomentum()->GetZ());
  
  auto fillPointForT = [&](double t){
    double x = helix->GetOrigin().GetX();
    double y = helix->GetOrigin().GetY();
    double z = helix->GetOrigin().GetZ() + fabs(helix->GetSlope())*t;
    
    if(helix->GetCharge()*zSign < 0){
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
//  gEve->GetEventScene()->DestroyElements();
//  gSystem->ProcessEvents();
	
  
  if(config->drawJets){
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
  }
    
  TEvePointSetArray *points = PreparePointsEventDisplay(options);
  
  for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
    shared_ptr<Track> track = event->GetTrack(iTrack);
    
    for(int iHit=0;iHit<nLayers;iHit++){
      
      int iLayer = track->GetLayerForHit(iHit) - 1;
      
      double phi    = track->GetPhi();
			double theta  = track->GetTheta();
      double decayR = layerR[iLayer] + 10*track->GetDxy();
      
      double x = decayR*cos(phi)              + 10*event->GetVertex()->GetX();
      double y = decayR*sin(phi)              + 10*event->GetVertex()->GetY();
      double z = decayR/sin(theta)*cos(theta) + 10*event->GetVertex()->GetZ() + 10*track->GetDz();
      
      points->Fill(scale*x,scale*y,scale*z,track->GetDeDxForHit(iHit));
    }
  }
  points->SetRnrSelf(kTRUE);
  gEve->AddElement(points);
  gEve->Redraw3D();
  
  // MET
  if(config->drawMET){
    double metPhi = event->GetMetPhi();
    double metTheta = 2*atan(exp(-event->GetMetEta()));
    DrawMET(metPhi, metTheta);
  }
    
  // Geometry:
  if(config->showGeometryPixel){
    for(int i=0;i<4;i++){
      TGeoTube *pixelTube = new TGeoTube(scale*layerR[i]-0.2,scale*layerR[i]+0.2, scale*pixelBarrelZsize);
      TEveGeoShape *pixel = new TEveGeoShape (Form("Pixel tracker %i",i),Form("Pixel tracker %i",i));
      pixel->SetShape(pixelTube);
      pixel->SetMainTransparency(geomTransparency-30);
      pixel->SetMainColorRGB((Float_t)0.0, 1.0, 0.0);
      pixel->SetRnrSelf(true);
      gEve->AddElement(pixel);
    }
  }
  
  if(config->showGeometryStrip){
    TGeoTube *trackerTube = new TGeoTube(scale*layerR[4],scale*layerR[nLayers-1],scale*2800);
    TEveGeoShape *tracker = new TEveGeoShape ("Tracker","Tracker");
    tracker->SetShape(trackerTube);
    tracker->SetMainTransparency(geomTransparency-20);
    tracker->SetMainColorRGB((Float_t)1.0, 1.0, 0.0);
    tracker->SetRnrSelf(true);
    gEve->AddElement(tracker);
  }
  
  if(config->showGeometryEcal){
    TGeoTube *emCalTube = new TGeoTube(scale*1100,scale*1800,scale*3700);
    TEveGeoShape *emCal = new TEveGeoShape ("EM calo","EM calo");
    emCal->SetShape(emCalTube);
    emCal->SetMainTransparency(geomTransparency-10);
    emCal->SetMainColorRGB((Float_t)0.0, 0.0, 1.0);
    emCal->SetRnrSelf(true);
    gEve->AddElement(emCal);
  }
  
  if(config->showGeometryHcal){
    TGeoTube *hadCalTube = new TGeoTube(scale*1800,scale*2900,scale*5500);
    TEveGeoShape *hadCal = new TEveGeoShape ("Had calo","Had calo");
    hadCal->SetShape(hadCalTube);
    hadCal->SetMainTransparency(geomTransparency);
    hadCal->SetMainColorRGB((Float_t)1.0, 0.0, 0.5);
    hadCal->SetRnrSelf(true);
    gEve->AddElement(hadCal);
  }

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
