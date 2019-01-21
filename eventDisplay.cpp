#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "FitterConfig.hpp"
#include "HelixProcessor.hpp"

Display *display;
shared_ptr<FitterConfig> config;

bool showStipClusters = false;

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

const map<string,any> decayPointOptions = {
  {"title", "Decay Point"},
  {"markerStyle", 22},
  {"markerSize", 2.0},
  {"color", kGreen}
};

const map<string,any> pionPointsOptions = {
  {"title", "Pion points"},
  {"markerStyle", 21},
  {"markerSize", 1.6},
  {"color", kMagenta}
};

const map<string,any> helixOptions = {
  {"title", "Helix"},
  {"markerStyle", 20},
  {"markerSize", 0.2},
  {"color", kGreen}
};

const map<string,any> filteredPointsOptions = {
  {"title", "Filtered Points"},
  {"markerStyle", 20},
  {"markerSize", 1.0},
  {"color", kYellow}
};

// Parameters for all hits in the pixel barrel
const double chargeThreshold = 0; // 2000, 5000, 25000
const double minClusterSize = 0;
const double maxClusterSize = 100;

shared_ptr<vector<Point>> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber);

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  config = make_shared<FitterConfig>("configs/helixFitter.md");
  display = new Display();
  unique_ptr<HelixProcessor> helixProcessor = make_unique<HelixProcessor>(config);
  
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

  display->DrawEvent(event, dedxOptions);
  event->Print();

  // ------------------------------------------------------------------------------------------------------------
  // Helix fitting part
  // ------------------------------------------------------------------------------------------------------------
  
  cout<<"Preparing hits, track and pion's helix"<<endl;
  shared_ptr<vector<Point>> allSimplePoints; // all hits in the event
  allSimplePoints = LoadAllHits(searchRun, searchLumi, searchEvent);
  
  shared_ptr<Track> track = make_shared<Track>();
  track->FillRandomly(config->GetNTrackHits(), config->GetMaxTrackEta());
  
  // Draw decay point to make sure that it's correctly located
  shared_ptr<vector<Point>> decayPoint = make_shared<vector<Point>>();
  decayPoint->push_back(Point(track->GetDecayPoint()->GetX(),
                              track->GetDecayPoint()->GetY(),
                              track->GetDecayPoint()->GetZ()));
  display->DrawSimplePoints(decayPoint, decayPointOptions);
  
  // Draw true pion helix
//  double pionCharge = -1;
//  unique_ptr<Point> pionVector = make_unique<Point>(120,120,-800); // Total momentum ~200 MeV
  unique_ptr<Helix> pionHelix = helixProcessor->GetRandomPionHelix(track);
//   make_unique<Helix>(track->GetDecayPoint(), pionVector, pionCharge, config);
  display->DrawHelix(pionHelix,helixOptions);
  
  // Calculate and draw points along the helix that hit the silicon
  shared_ptr<vector<Point>> pionPoints = helixProcessor->GetPointsHittingSilicon(pionHelix);
  for(auto &p : *pionPoints){p.SetIsPionHit(true);}
  pionHelix->SetPoints(pionPoints);
  display->DrawSimplePoints(pionPoints, pionPointsOptions);
  
  // inject hits from pion into all points in the tracker
  allSimplePoints->insert(allSimplePoints->end(),pionPoints->begin(), pionPoints->end());
  
  // remove hits that for sure don't belong to the pion's helix
  cout<<"Fitting best helix"<<endl;
  unique_ptr<Fitter> fitter = make_unique<Fitter>(config);
  
  unique_ptr<Helix> bestHelix = fitter->GetBestFittingHelix(allSimplePoints, track);
//  display->DrawSimplePoints(allSimplePoints, filteredPointsOptions);
  
  if(bestHelix){
    map<string,any> bestHelixOptions = {
      {"title", "Best helix"},
      {"markerStyle", 20},
      {"markerSize", 0.2},
      {"color", kRed}
    };
    
    display->DrawHelix(bestHelix, bestHelixOptions);
    
    map<string,any> fitPointsOptions = {
      {"title", "Fit helix points"},
      {"markerStyle", 20},
      {"markerSize", 1.0},
      {"color", kRed}
    };
    display->DrawSimplePoints(bestHelix->GetPoints(), fitPointsOptions);
    
    cout<<"\n\nBest helix is:"<<endl;
    bestHelix->Print();
    
    cout<<"\nPion helix:"<<endl;
    pionHelix->Print();
  }
    
  
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
  
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}








shared_ptr<vector<Point>> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber)
{
  TFile *inFile = TFile::Open("pickhists.root");
//  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists.root");
  //  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists_unfiltered.root");
  if(!inFile){
    cout<<"ERROR -- no file with all hits was found"<<endl;
    return make_shared<vector<Point>>();
  }
  TTree *tree = (TTree*)inFile->Get("hitsExtractor/hits");
  
  if(!tree){
    cout<<"ERROR -- no tree with all hits was found"<<endl;
    return make_shared<vector<Point>>();
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
  
  shared_ptr<vector<Point>> pixelPoints = make_shared<vector<Point>>();
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find all hits for requested event!\n\n"<<endl;
    return pixelPoints;
  }
  
  for(uint i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    // convert cm to mm
    pixelPoints->push_back(Point(10*hitX->at(i),10*hitY->at(i),10*hitZ->at(i),hitCharge->at(i)));
  }
  
  display->DrawSimplePoints(pixelPoints, pixelHitsOptions);
  
  if(showStipClusters){
    shared_ptr<vector<Point>> stripPoints = make_shared<vector<Point>>();
    for(uint i=0;i<stripX->size();i++){
      if(stripCharge->at(i) < chargeThreshold) continue;
      stripPoints->push_back(Point(10*stripX->at(i),10*stripY->at(i),10*stripZ->at(i),stripCharge->at(i)));
    }
    display->DrawSimplePoints(stripPoints, stripHitsOptions);
  }
  
  return pixelPoints;
}
