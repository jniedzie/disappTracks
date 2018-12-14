#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "HelixFitter.hpp"
#include "Fitter.hpp"

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
#include <TF3.h>
#include <TH3F.h>

const double scale = 0.1;

const bool showUnderflowBins = false;
const bool showOverflowBins = true;

const bool showStipClusters = false;
const bool showGeometry = false;

const map<string,any> dedxOptions = {
  {"title", "dE/dx clusters"},
  {"binsMin" , 1},
  {"binsMax" , 10},
  {"nBins" , 5},
  {"markerStyle", 20},
  {"markerSize", 2.0}
};

const map<string,any> pixelHitsOptions = {
  {"title", "Pixel hits"},
  {"binsMin" , 0},
  {"binsMax" , 100000},
  {"nBins" , 5},
  {"markerStyle", 20},
  {"markerSize", 0.2}
};

const map<string,any> stripHitsOptions = {
  {"title", "Strip hits"},
  {"binsMin" , 0},
  {"binsMax" , 5000},
  {"nBins" , 5},
  {"markerStyle", 20},
  {"markerSize", 1.0}
};

const double jetConeRadius = scale*0.4;

const double metRadius = scale * 2000;
const double metBoxSize = scale * 30;
const double metBoxAngularSize = 0.1;

const int geomTransparency = 90; // 30 - 100

                            //     Underflow                                       Overflow
const vector<int> dedxBinColors = { kGray,     kBlue, kCyan, kGreen, kYellow, kRed, kMagenta };

// Parameters for all hits in the pixel barrel
const double chargeThreshold = 0; // 2000, 5000, 25000
const double minClusterSize = 0;
const double maxClusterSize = 100;

TEvePointSetArray* PreparePointsEventDisplay(map<string,any> options)
{
  TEvePointSetArray *points = new TEvePointSetArray(any_cast<const char*>(options["title"]));
  points->SetMarkerStyle(any_cast<int>(options["markerStyle"]));
  points->SetMarkerSize(any_cast<double>(options["markerSize"]));
  
  points->InitBins(any_cast<const char*>(options["title"]),
                   any_cast<int>(options["nBins"]),
                   any_cast<int>(options["binsMin"]),
                   any_cast<int>(options["binsMax"]));

  for(int i=0;i<(any_cast<int>(options["nBins"])+2);i++){
    points->GetBin(i)->SetMainColor(dedxBinColors[i]);
  }
  
  points->GetBin(0)->SetRnrSelf(showUnderflowBins);
  points->GetBin(any_cast<int>(options["nBins"])+1)->SetRnrSelf(showOverflowBins);
  
  return points;
}

void DrawMET(double metPhi, double metTheta)
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
  
  TEvePointSetArray *points = PreparePointsEventDisplay(dedxOptions);
  
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

vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber)
{
  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists.root");
//  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists_unfiltered.root");
  if(!inFile){
    cout<<"ERROR -- no file with all hits was found"<<endl;
    return vector<Point>();
  }
  TTree *tree = (TTree*)inFile->Get("hitsExtractor/hits");
  
  if(!tree){
    cout<<"ERROR -- no tree with all hits was found"<<endl;
    return vector<Point>();
  }
  
  vector<double> *hitX = nullptr;
  vector<double> *hitY = nullptr;
  vector<double> *hitZ = nullptr;
  vector<double> *hitCharge = nullptr;
  vector<double> *hitSizeX = nullptr;
  vector<double> *hitSizeY = nullptr;
  vector<double> *stripX = nullptr;
  vector<double> *stripY = nullptr;
  vector<double> *stripZ = nullptr;
  vector<double> *stripCharge = nullptr;
  
  uint run;
  uint lumi;
  unsigned long long event;
  
  tree->SetBranchAddress("hitX",&hitX);
  tree->SetBranchAddress("hitY",&hitY);
  tree->SetBranchAddress("hitZ",&hitZ);
  tree->SetBranchAddress("hitCharge",&hitCharge);
  tree->SetBranchAddress("hitSizeX",&hitSizeX);
  tree->SetBranchAddress("hitSizeY",&hitSizeY);
  tree->SetBranchAddress("stripX",&stripX);
  tree->SetBranchAddress("stripY",&stripY);
  tree->SetBranchAddress("stripZ",&stripZ);
  tree->SetBranchAddress("stripCharge",&stripCharge);
  
  tree->SetBranchAddress("runNumber",&run);
  tree->SetBranchAddress("lumiBlock",&lumi);
  tree->SetBranchAddress("eventNumber",&event);
  
  bool eventFound = false;
  
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    
    if(run == runNumber && lumi == lumiSection && event == eventNumber){
      eventFound = true;
      break;
    }
  }
  
  vector<Point> simplePoints;
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find all hits for requested event!\n\n"<<endl;
    return simplePoints;
  }
  
  TEvePointSetArray *pixelPoints = PreparePointsEventDisplay(pixelHitsOptions);
  
  for(int i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    
    // convert cm to mm
    pixelPoints->Fill(scale*10*hitX->at(i),
                      scale*10*hitY->at(i),
                      scale*10*hitZ->at(i),
                      hitCharge->at(i));
    
    simplePoints.push_back(Point(10*hitX->at(i),
                                              10*hitY->at(i),
                                              10*hitZ->at(i)));
  }
  
  pixelPoints->SetRnrSelf(kTRUE);
  gEve->AddElement(pixelPoints);
  gEve->Redraw3D();
  
  TEvePointSetArray *stripPoints = PreparePointsEventDisplay(stripHitsOptions);
  
  for(int i=0;i<stripX->size();i++){
    if(stripCharge->at(i) < chargeThreshold) continue;
    
    // convert cm to mm
    stripPoints->Fill(scale*10*stripX->at(i),
                      scale*10*stripY->at(i),
                      scale*10*stripZ->at(i),
                      stripCharge->at(i));
  }
  stripPoints->SetRnrSelf(showStipClusters);
  stripPoints->SetRnrChildren(showStipClusters);
  gEve->AddElement(stripPoints);
  gEve->Redraw3D();
  
  return simplePoints;
}

void DrawSimplePoints(vector<Point> points)
{
  TEvePointSetArray *simplePoints = PreparePointsEventDisplay(pixelHitsOptions);
  
  for(auto p : points){
    simplePoints->Fill(scale*p.x,scale*p.y,scale*p.z, 100000);
  }
  
  simplePoints->SetRnrSelf(kTRUE);
  gEve->AddElement(simplePoints);
  gEve->Redraw3D();
}

/// Draws a helix in given t parameter range
/// \param helix Object of type Helix
/// \param tMin t parameter minimum
/// \param tMax t parameter maximum
/// \param tStep t parameter step
void DrawHelix(Helix helix, double tMin=0, double tMax=5*2*TMath::Pi(), double tStep=0.01)
{
  TEvePointSetArray *helixPoints = PreparePointsEventDisplay(pixelHitsOptions);
  
  for(double t=tMin;t<tMax;t+=tStep){
    double x = helix.R*cos(t) + helix.x0;
    double y = helix.R*sin(t) + helix.y0;
    double z = helix.c*t      + helix.z0;
    
    helixPoints->Fill(scale*x,scale*y,scale*z, 10000);
  }
  
  helixPoints->SetRnrSelf(kTRUE);
  gEve->AddElement(helixPoints);
  gEve->Redraw3D();
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  TEveManager::Create();
  
  auto events = shared_ptr<EventSet>(new EventSet());
  events->LoadEventsFromFiles("after_L2/3layers/");
//  events->LoadEventsFromFiles("after_L0/");
  
//  auto event = events->At(EventSet::kSignal, kWino_M_300_cTau_3, 0);
  
  uint searchRun = 297100;
  uint searchLumi = 136;
  unsigned long long searchEvent = 245000232;

  auto event = events->GetEvent(EventSet::kData, searchRun, searchLumi, searchEvent);
  if(!event){
    cout<<"event not found"<<endl;
    exit(0);
  }

  DrawEvent(event);
  event->Print();

  
  // Helix fitting part
  vector<Point> allSimplePoints; // all hits in the event
  allSimplePoints = LoadAllHits(searchRun, searchLumi, searchEvent);
  

  double decayR = 140; // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
  double theta = 2*atan(exp(-event->GetTrack(0)->GetEta()));
  double phi = event->GetTrack(0)->GetPhi();
  
  double pionR = 40; // radius of the pion spiral
  double pionC = 4;  // helix slope in Z direction
  
  // location of the secondary vertex
  double decayX = decayR*sin(theta)*cos(phi);
  double decayY = decayR*sin(theta)*sin(phi);
  double decayZ = decayR*cos(theta);
  
  vector<Point> decayPoint = {Point(decayX,decayY,decayZ)};
  DrawSimplePoints(decayPoint);
  
  // true pion helix (has to be moved so that it begins in the secondary vertex, rather than has the center there)
  Point pionHelixCenter(decayX,decayY,decayZ);
  pionHelixCenter.PerpendicularShift(pionR);
  
  double helixThickness = 3.0; // gather only points within X mm
  
  Helix pionHelix(pionR,pionC,
                  pionHelixCenter.x,
                  pionHelixCenter.y,
                  pionHelixCenter.z,
                  helixThickness);
  
  DrawHelix(pionHelix,TMath::Pi(),30,0.01);
  
  vector<Point> pionPoints = pionHelix.GetPointsHittingSilicon(TMath::Pi(),5*2*TMath::Pi());
  for(auto &p : pionPoints){
    p.isPionHit = true;
    pionHelix.nPionPoints++;
  }
  DrawSimplePoints(pionPoints);
  pionHelix.nPoints = (int)pionPoints.size();
  
  // inject hits from pion into all points in the tracker
  allSimplePoints.insert(allSimplePoints.end(),pionPoints.begin(), pionPoints.end());
  
  // remove hits that for sure don't belong to the pion's helix
  cout<<"size before:"<<allSimplePoints.size()<<endl;
  for(int i=0;i<allSimplePoints.size();i++){
    Point p = allSimplePoints[i];
    
    if(   (p.z * decayZ < 0) // remove wrong Z points
//       || (sqrt(pow(p.x - decayX,2) + pow(p.y - decayY,2)) > 2*pionR)
       ){
      allSimplePoints.erase(allSimplePoints.begin()+i);
      i--;
    }
  }
  cout<<"size after:"<<allSimplePoints.size()<<endl;
  
  // Set fitter limits (all distances in [mm])
  HelixFitter::helixThickness = helixThickness; // gather only points within X mm
  
  // Pion helix parameters:
  double minR   = 25; // can't be smaller than the minimum to reach 2 different layers. Physically, from pion momentum it should be around 130 mm
  double maxR   = 70;
  double startR = pionR;
  double minC   = 3; // Fitter could try to make it very small, to gahter some additional points by chance while still fitting perfectly the real pion hits. This should be given more attention later...
  double maxC   = 70;
  HelixFitter::startC = pionC;
  
  // Position of the decay vertex along chargino's track (later should take into account that it's not a straight line)
  HelixFitter::minL   = layerR[2]; // minimum on the surface of 3rd layer
  HelixFitter::maxL   = layerR[3]; // maximum on the surface of 4th layer
  HelixFitter::startL = decayR;//(layerR[2]+layerR[3])/2.; // starting position between 3rd and 4th layer
  
  // Theta and phi of the chargino
  HelixFitter::theta = theta;
  HelixFitter::phi = phi;
  
  HelixFitter::InitHelixFitter();
//  HelixFitter::InitHelixFitterWithSeeds();
  HelixFitter::points = allSimplePoints;
  
  TH2D *pointsXY = new TH2D("pointsXY","pointsXY",50,-250,250,50,-250,250);
  pointsXY->GetXaxis()->SetTitle("X");
  pointsXY->GetYaxis()->SetTitle("Y");
  for(auto point : allSimplePoints){pointsXY->Fill(point.x,point.y);}
  
  vector<pair<double, double>> points2D;
  
  for(int binX=0;binX<pointsXY->GetNbinsX();binX++){
    for(int binY=0;binY<pointsXY->GetNbinsY();binY++){
      if(pointsXY->GetBinContent(binX,binY) < 2){
        pointsXY->SetBinContent(binX,binY, 0);
      }
      else{
        points2D.push_back(make_pair(pointsXY->GetXaxis()->GetBinCenter(binX),
                                     pointsXY->GetYaxis()->GetBinCenter(binY)));
      }
    }
  }
  
  auto chi2Function = [&](const double *par) {
    double f = 0;
    
    double x0 = par[0];
    double y0 = par[1];
    double R  = par[2];
    
    double xa = par[3];
    double ya = par[4];
    double xb = par[5];
    double yb = par[6];
    double xc = par[7];
    double yc = par[8];
    
    f  = pow(sqrt(pow(xa-x0,2)+pow(ya-y0,2))-R,2);
    f += pow(sqrt(pow(xb-x0,2)+pow(yb-y0,2))-R,2);
    f += pow(sqrt(pow(xc-x0,2)+pow(yc-y0,2))-R,2);
    
    return f;
  };
  
  double minDecayX = HelixFitter::minL*sin(theta)*cos(phi);
  double maxDecayX = HelixFitter::maxL*sin(theta)*cos(phi);
  
  double minDecayY = HelixFitter::minL*sin(theta)*sin(phi);
  double maxDecayY = HelixFitter::maxL*sin(theta)*sin(phi);
  
  double minX = minDecayX - HelixFitter::maxR;
  double maxX = maxDecayX + HelixFitter::maxR;
  double minY = minDecayY - HelixFitter::maxR;
  double maxY = maxDecayY + HelixFitter::maxR;
  
  Fitter *fitter = new Fitter(9);
  fitter->SetFitFunction(chi2Function);
  
  fitter->SetParameter(0, "x0", (maxX+minX)/2., minX, maxX);
  fitter->SetParameter(1, "y0", (maxY+minY)/2., minY, maxY);
  fitter->SetParameter(2, "R" , (maxR+minR)/2., minR, maxR);
  
  vector<Circle> circles;
  
  cout<<"Limits in X:"<<minX<<"\t"<<maxX<<endl;
  cout<<"Limits in Y:"<<minY<<"\t"<<maxY<<endl;
  
  int nPoints = 0;//(int)points2D.size();
  TH1D *chi2hist = new TH1D("chi2hist","chi2hist",1000,0,10);
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){
        fitter->FixParameter(3, "xa", points2D[i].first);
        fitter->FixParameter(4, "ya", points2D[i].second);
        fitter->FixParameter(5, "xb", points2D[j].first);
        fitter->FixParameter(6, "yb", points2D[j].second);
        fitter->FixParameter(7, "xc", points2D[k].first);
        fitter->FixParameter(8, "yc", points2D[k].second);

        if(fitter->RunFitting()) {
          auto result = fitter->GetResult();
          chi2hist->Fill(result.MinFcnValue());
          
          double x0 = result.GetParams()[0];
          double y0 = result.GetParams()[1];
          double R = result.GetParams()[2];
          
          Circle circle(x0,y0,R);
          
          int nCircleBins = circle.GetNbinsOverlappingWithHist(pointsXY);
          
          if(nCircleBins > 3 && nCircleBins < 7){
            circles.push_back(circle);
//            result.Print(std::cout);
          }
        }
      }
    }
  }
  
//  HelixFitter::fitFunctionWithSeeds(a, nullptr, chi2, par, 0);

  TH2D *hist = new TH2D("hist","hist",
                        (HelixFitter::maxR-HelixFitter::minR)+1,HelixFitter::minR,HelixFitter::maxR,
                        (HelixFitter::maxC-HelixFitter::minC)+1,HelixFitter::minC,HelixFitter::maxC);
  
  TH2D *hist2 = new TH2D("hist2","hist2",
                        (HelixFitter::maxR-HelixFitter::minR)+1,HelixFitter::minR,HelixFitter::maxR,
                        (HelixFitter::maxL-HelixFitter::minL)+1,HelixFitter::minL,HelixFitter::maxL);
  
  /*
  //  chi2            R     c     L
  map<double,tuple<double,double,double>> results;
  
  for(double R=HelixFitter::minR;R<HelixFitter::maxR;R+=1.0){
    for(double c=HelixFitter::minC;c<HelixFitter::maxC;c+=1.0){
      for(double L=HelixFitter::minL;L<HelixFitter::maxL;L+=1.0){
//      for(double L=140;L<141;L+=1.0){
        int a;
        double chi2;
        double par[3] = {R,c,L};
        
        HelixFitter::fitFunction(a, nullptr, chi2, par, 0);
        if(chi2 < 1000){
          hist->Fill(R,c,chi2);
          hist2->Fill(R,L,chi2);
          
          results[chi2] = make_tuple(R,c,L);
        }
        if(chi2 < bestChi2){
          bestChi2 = chi2;
          bestR = R;
          bestC = c;
          bestL = L;
        }
      }
    }
  }
  */
 
  TCanvas *c1 = new TCanvas("c1","c1",12800,800);
  
  c1->Divide(3,3);
  c1->cd(1);
  hist->DrawCopy("colz");
  c1->cd(2);
  chi2hist->Draw();
//  hist->ProjectionX()->DrawCopy();
  c1->cd(3);
  
  pointsXY->Draw("colz");
  
  TH1D *radii = new TH1D("radii","radii",
                         (HelixFitter::maxR-HelixFitter::minR)+1,HelixFitter::minR,HelixFitter::maxR);
  
  TH1D *Xs = new TH1D("Xs","Xs",
                      (maxX-minX)+1,minX,maxX);
  
  TH1D *Ys = new TH1D("Ys","Ys",
                      (maxY-minY)+1,minY,maxY);
  
  for(auto circle : circles){
    TArc *circleArc = circle.GetArc();
    circleArc->SetLineColor(kRed);
    circleArc->SetLineWidth(2);
    circleArc->ResetAttFill();
    circleArc->Draw("sameL");
    
    radii->Fill(circle.R);
    Xs->Fill(circle.x);
    Ys->Fill(circle.y);
  }
  c1->cd(5);
  radii->Draw();
  c1->cd(6);
  Xs->Draw();
  c1->cd(7);
  Ys->Draw();
  
  cout<<"N matching circles:"<<circles.size()<<endl;
  
  Circle trueCircle(pionHelixCenter.x,pionHelixCenter.y,pionR);
  int nTruePoints = trueCircle.GetNbinsOverlappingWithHist(pointsXY);
  cout<<"N points on true pion helix:"<<nTruePoints<<endl;
  
  TArc *trueCircleArc = trueCircle.GetArc();
  trueCircleArc->SetLineColor(kGreen);
  trueCircleArc->SetLineWidth(2);
  trueCircleArc->ResetAttFill();
  trueCircleArc->Draw("sameL");
  
  c1->cd(4);
  hist2->DrawCopy("colz");
  c1->Update();
  
//  HelixFitter::RunFitter();
//  Helix fittedHelix = HelixFitter::GetFittedHelix();
//  HelixFitter::CountPionPointsOnHelix(fittedHelix);
//  Helix fittedHelix = HelixFitter::bestSeedHelix;
  
  cout<<"\nPion helix:"<<endl;
  pionHelix.Print();
  
//  cout<<"\n\nFitted helix:"<<endl;
//  fittedHelix.Print();

//  DrawHelix(fittedHelix,TMath::Pi(),5*2*TMath::Pi());
  
  
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
  
  
  auto chi2Function3D = [&](const Double_t *par) {
    Helix helix(par[3],par[4], par[0], par[1], par[2],helixThickness);
    
    int nPoints=0;
    double chi2 = 0;
    
    for(Point p : allSimplePoints){
      Point q = helix.GetClosestPoint(p);
      double d = p.distance(q);
      
      if(d < HelixFitter::helixThickness){
        nPoints++;
        chi2 += pow(p.x-q.x,2)/fabs(q.x);
        chi2 += pow(p.y-q.y,2)/fabs(q.y);
        chi2 += pow(p.z-q.z,2)/fabs(q.z);
        
//        chi2 += 1+d*d;
      }
    }
    if(nPoints > 2) chi2 /= nPoints;
    else            chi2 = 99999;
//    cout<<"z0:"<<par[2]<<"\tc:"<<par[4]<<"\tchi2:"<<chi2<<endl;
    return chi2;
  };
  
  Fitter *fitter3D = new Fitter(5);
  fitter->SetFitFunction(chi2Function3D);
  
  circles.clear();
  circles.push_back(trueCircle);
  
  allSimplePoints.clear();
  allSimplePoints = pionPoints;
  
  Helix bestHelixFromCircles;
  double bestChi2fromCircles = 99999;
  
  ROOT::Fit::BinData data(100,2);
  
  for(auto p : allSimplePoints){
    data.Add(p.x,p.y,p.z);
  }
  
  bool fixCircles = true;
  
  for(auto circle : circles){
    fitter3D->SetParameter(0, "x0", circle.x, 0.95*circle.x, 1.05*circle.x, fixCircles);
    fitter3D->SetParameter(1, "y0", circle.y, 0.95*circle.y, 1.05*circle.y, fixCircles);
    fitter3D->SetParameter(3, "R" , circle.R, 0.95*circle.R, 1.05*circle.R, fixCircles);
    
    fitter3D->SetParameter(2, "z0", decayZ, decayZ-maxR, decayZ+maxR);
    fitter3D->SetParameter(4, "c" , pionC , minC, maxC);
    
    bool fitOk =true;// fitter3D.FitFCN();
//    fitter3D.CalculateHessErrors();
    
    auto result = fitter3D->GetResult();
    result.Print(std::cout);
    cout<<"Fit ok? - "<<fitOk<<endl;
    
    if(result.MinFcnValue() < bestChi2fromCircles && fitOk){
      bestChi2fromCircles = result.MinFcnValue();
      
      bestHelixFromCircles = Helix(result.GetParams()[3],
                                   result.GetParams()[4],
                                   result.GetParams()[0],
                                   result.GetParams()[1],
                                   result.GetParams()[2],
                                   helixThickness);
      
      cout<<result.Chi2()<<"\t"<<result.MinFcnValue()<<endl;
      
    }
  }
  
  DrawHelix(bestHelixFromCircles);
  
  cout<<"\n\ndone\n\n"<<endl;
  
  theApp.Run();
  return 0;
}
