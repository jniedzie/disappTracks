#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "TrackProcessor.hpp"

uint searchRun = 297100;
uint searchLumi = 136;
unsigned long long searchEvent = 245000232;

string configPath = "configs/eventDisplay.md";

Display *display;

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

map<string,any> filteredPointsOptions = {
  {"title", "Filtered Points"},
  {"markerStyle", 20},
  {"markerSize", 1.0},
  {"color", kYellow}
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  config = make_unique<ConfigManager>(configPath);
  display = new Display();
  auto helixProcessor = make_unique<HelixProcessor>();
  auto trackProcessor = make_unique<TrackProcessor>();
  
  auto events = make_shared<EventSet>();
//  events->LoadEventsFromFiles("/");
  events->LoadEventsFromFiles("after_L1/");
  
//  auto event = events->GetEvent(EventSet::kBackground, searchRun, searchLumi, searchEvent);
  auto event = events->At(xtracks::kSignal, kWino_M_300_cTau_10, 0);
  
  if(!event){
    cout<<"eventDisplay -- event not found"<<endl;
    exit(0);
  }

  auto jetCut = make_unique<JetCut>();
  
  jetCut->SetPt(range<double>(30.0, inf));
  jetCut->SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut->SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
  
  auto eventProcessor = make_unique<EventProcessor>();
  eventProcessor->ApplyJetCut(event, jetCut);
  
  display->DrawEvent(event, dedxOptions);
  event->Print();

  // ------------------------------------------------------------------------------------------------------------
  // Helix fitting part
  // ------------------------------------------------------------------------------------------------------------
  
  cout<<"Preparing hits, track and pion's helix"<<endl;
  
  shared_ptr<vector<Point>> allSimplePoints = event->GetTrackerHits();
  display->DrawSimplePoints(allSimplePoints, filteredPointsOptions);
  
  shared_ptr<vector<unique_ptr<Helix>>> truePionHelices = event->GetTruePionHelices();
  
  for(auto &helix : *truePionHelices){
    display->DrawHelix(helix,helixOptions);
    helix->SetPoints(allSimplePoints);
    auto helixPoints = helix->GetPoints();
    filteredPointsOptions["title"] = "true helix points";
    filteredPointsOptions["color"] = kRed;
    display->DrawSimplePoints(helixPoints, filteredPointsOptions);
    
    helix->Print();
  }
  
//  shared_ptr<Track> track = trackProcessor->GetRandomTrack(config->nTrackHits, config->maxEta);
  shared_ptr<Track> track = event->GetTrack(0);
  
  // Draw decay point to make sure that it's correctly located
//  shared_ptr<vector<Point>> decayPoint = make_shared<vector<Point>>();
//  decayPoint->push_back(Point(track->GetDecayPoint()->GetX(),
//                              track->GetDecayPoint()->GetY(),
//                              track->GetDecayPoint()->GetZ()));
//  display->DrawSimplePoints(decayPoint, decayPointOptions);
  
  // Draw true pion helix
  
  
//  unique_ptr<Helix> pionHelix = helixProcessor->GetRandomPionHelix(track);
//  display->DrawHelix(pionHelix,helixOptions);
  
  // Calculate and draw points along the helix that hit the silicon
//  shared_ptr<vector<Point>> pionPoints = helixProcessor->GetPointsHittingSilicon(pionHelix);
//  for(auto &p : *pionPoints){p.SetIsPionHit(true);}
//  pionHelix->SetPoints(pionPoints);
//  display->DrawSimplePoints(pionPoints, pionPointsOptions);
  
  // inject hits from pion into all points in the tracker
//  allSimplePoints->insert(allSimplePoints->end(),pionPoints->begin(), pionPoints->end());
  
  // remove hits that for sure don't belong to the pion's helix
  cout<<"Fitting best helix"<<endl;
  unique_ptr<Fitter> fitter = make_unique<Fitter>();
  
  unique_ptr<Helix> bestHelix = nullptr; //fitter->GetBestFittingHelix(allSimplePoints, track);
  
  
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
    
//    cout<<"\nPion helix:"<<endl;
//    pionHelix->Print();
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
