#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "HelixFitter.hpp"

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

vector<HelixFitter::Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber)
{
  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists.root");
//  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists_unfiltered.root");
  if(!inFile){
    cout<<"ERROR -- no file with all hits was found"<<endl;
    return vector<HelixFitter::Point>();
  }
  TTree *tree = (TTree*)inFile->Get("hitsExtractor/hits");
  
  if(!tree){
    cout<<"ERROR -- no tree with all hits was found"<<endl;
    return vector<HelixFitter::Point>();
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
  
  vector<HelixFitter::Point> simplePoints;
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find all hits for requested event!\n\n"<<endl;
    return simplePoints;
  }
  
  TEvePointSetArray *pixelPoints = PreparePointsEventDisplay(pixelHitsOptions);
  
  for(int i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    
    pixelPoints->Fill(hitX->at(i),
                 hitY->at(i),
                 hitZ->at(i),
                 hitCharge->at(i));
    
    simplePoints.push_back(HelixFitter::Point(hitX->at(i),
                                              hitY->at(i),
                                              hitZ->at(i)));
  }
  
  pixelPoints->SetRnrSelf(kTRUE);
  gEve->AddElement(pixelPoints);
  gEve->Redraw3D();
  
  TEvePointSetArray *stripPoints = PreparePointsEventDisplay(stripHitsOptions);
  
  for(int i=0;i<stripX->size();i++){
    if(stripCharge->at(i) < chargeThreshold) continue;
    
    stripPoints->Fill(stripX->at(i),
                      stripY->at(i),
                      stripZ->at(i),
                      stripCharge->at(i));
  }
  stripPoints->SetRnrSelf(showStipClusters);
  stripPoints->SetRnrChildren(showStipClusters);
  gEve->AddElement(stripPoints);
  gEve->Redraw3D();
  
  return simplePoints;
}

void DrawSimplePoints(vector<HelixFitter::Point> points)
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
/// \param helix Object of type HelixFitter::Helix
/// \param tMin t parameter minimum
/// \param tMax t parameter maximum
/// \param tStep t parameter step
void DrawHelix(HelixFitter::Helix helix, double tMin=0, double tMax=5*2*TMath::Pi(), double tStep=0.01)
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

/// Returns vector of points along helix trajectory that hit the tracker
/// \param helix Object of type HelixFitter::Helix
/// \param tMin t parameter minimum
/// \param tMax t parameter maximum
/// \param tStep t parameter step
vector<HelixFitter::Point> GetHelixPointsHittingSilicon(HelixFitter::Helix helix,
                                                         double tMin=0, double tMax=5*2*TMath::Pi(), double tStep=0.01)
{
  vector<HelixFitter::Point> points;
  double threshold = scale*1.0; // how close to the tracker layer hits must be
  
  double x,y,z;
  for(double t=tMin;t<tMax;t+=tStep){
    x = helix.R*cos(t) + helix.x0;
    y = helix.R*sin(t) + helix.y0;
    z = helix.c*t      + helix.z0;
   
    for(int iLayer=0;iLayer<5;iLayer++){
      if(fabs(sqrt(x*x+y*y)-layerR[iLayer]) < threshold){
        points.push_back(HelixFitter::Point(x,y,z));
      }
    }
  }
  return points;
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
  vector<HelixFitter::Point> allSimplePoints; // all hits in the event
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
  
  vector<HelixFitter::Point> decayPoint = {HelixFitter::Point(decayX,decayY,decayZ)};
  DrawSimplePoints(decayPoint);
  
  // true pion helix (has to be moved so that it begins in the secondary vertex, rather than has the center there)
  HelixFitter::Helix pionHelix(pionR,pionC,
                               decayX+pionR,
                               decayY,
                               decayZ-pionC*TMath::Pi());
  DrawHelix(pionHelix,TMath::Pi(),30,0.01);
  
  vector<HelixFitter::Point> pionPoints = GetHelixPointsHittingSilicon(pionHelix,TMath::Pi(),5*2*TMath::Pi());
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
    HelixFitter::Point p = allSimplePoints[i];
    
    if(   (p.z * decayZ < 0) // remove wrong Z points
//       || (sqrt(pow(p.x - decayX,2) + pow(p.y - decayY,2)) > 2*pionR)
       ){
      allSimplePoints.erase(allSimplePoints.begin()+i);
      i--;
    }
  }
  cout<<"size after:"<<allSimplePoints.size()<<endl;
  
  // Set fitter limits (all distances in [mm])
  HelixFitter::helixThickness = 5.0; // gather only points within X mm
  
  // Pion helix parameters:
  HelixFitter::minR   = 25; // can't be smaller than the minimum to reach 2 different layers. Physically, from pion momentum it should be around 130 mm
  HelixFitter::maxR   = 120;
  HelixFitter::startR = pionR;
  HelixFitter::minC   = 3; // Fitter could try to make it very small, to gahter some additional points by chance while still fitting perfectly the real pion hits. This should be given more attention later...
  HelixFitter::maxC   = 70;
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
  
  
  
//  HelixFitter::fitFunctionWithSeeds(a, nullptr, chi2, par, 0);

  TH2D *hist = new TH2D("hist","hist",
                        100,HelixFitter::minR,HelixFitter::maxR,
                        100,HelixFitter::minC,HelixFitter::maxC);
  
  TH2D *hist2 = new TH2D("hist2","hist2",
                        100,HelixFitter::minR,HelixFitter::maxR,
                        100,HelixFitter::minL,HelixFitter::maxL);
  
  double bestChi2 = 9999;
  double bestR=-1, bestC=-1, bestL=-1;
  
  for(double R=HelixFitter::minR;R<HelixFitter::maxR;R+=1.0){
    for(double c=HelixFitter::minC;c<HelixFitter::maxC;c+=1.0){
//      for(double L=HelixFitter::minL;L<HelixFitter::maxL;L+=1.0){
      for(double L=140;L<141;L+=1.0){
        int a;
        double chi2;
        double par[3] = {R,c,L};
        
        HelixFitter::fitFunction(a, nullptr, chi2, par, 0);
        if(chi2 < 1000){
          hist->Fill(R,c,chi2);
          hist2->Fill(R,L,chi2);
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
  
  TCanvas *c1 = new TCanvas("c1","c1",12800,800);
  c1->Divide(2,2);
  c1->cd(1);
  hist->DrawCopy("colz");
  c1->cd(2);
  hist->ProjectionX()->DrawCopy();
  c1->cd(3);
  hist->ProjectionY()->DrawCopy();
  c1->cd(4);
  hist2->DrawCopy("colz");
  c1->Update();
  
  cout<<"\nBest chi2:"<<bestChi2<<"\tfor R:"<<bestR<<"\tc:"<<bestC<<"\tL:"<<bestL<<endl;
  
  double x0 = bestL*sin(theta)*cos(phi);
  double y0 = bestL*sin(theta)*sin(phi);
  double z0 = bestL*cos(theta);
  
  HelixFitter::Helix helix(bestR,bestC, x0+bestR, y0, z0-bestC*TMath::Pi());
  HelixFitter::CountPionPointsOnHelix(helix);
  
  cout<<"Best helix:"<<endl;
  helix.Print();
  
//  HelixFitter::RunFitter();
//  HelixFitter::Helix fittedHelix = HelixFitter::GetFittedHelix();
//  HelixFitter::CountPionPointsOnHelix(fittedHelix);
//  HelixFitter::Helix fittedHelix = HelixFitter::bestSeedHelix;
  
  cout<<"\nPion helix:"<<endl;
  pionHelix.Print();
  
//  cout<<"\n\nFitted helix:"<<endl;
//  fittedHelix.Print();
//
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
  
  theApp.Run();
  return 0;
}
