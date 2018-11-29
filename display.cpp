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

void LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber)
{
  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists.root");
  TTree *tree = (TTree*)inFile->Get("hitsExtractor/hits");
  
  vector<double> *hitX = nullptr;
  vector<double> *hitY = nullptr;
  vector<double> *hitZ = nullptr;
  vector<double> *hitCharge = nullptr;
  vector<double> *hitSizeX = nullptr;
  vector<double> *hitSizeY = nullptr;
  vector<double> *detX = nullptr;
  vector<double> *detY = nullptr;
  vector<double> *detZ = nullptr;
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
  tree->SetBranchAddress("detX",&detX);
  tree->SetBranchAddress("detY",&detY);
  tree->SetBranchAddress("detZ",&detZ);
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
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find all hits for requested event!\n\n"<<endl;
    return;
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
  stripPoints->SetRnrSelf(kTRUE);
  gEve->AddElement(stripPoints);
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
  
  LoadAllHits(searchRun, searchLumi, searchEvent);
  
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
