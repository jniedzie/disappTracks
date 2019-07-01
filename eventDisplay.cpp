#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "TrackProcessor.hpp"

uint searchRun = 1;
uint searchLumi = 1;
unsigned long long searchEvent = 2662;

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L1/all/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 2;

bool injectPion = false;
bool fitHelix = true;

bool pionHitsOnly = true;

// "after_L1/all/":
// ok: 0,

// no chance:
// no hits: 1

// "after_L2/all/":
// ok: 1, 2, 4, 8, 14, 15, 17, 19, 40, 41, 42, 44
//
// ok, but could be improved:
// 5: poor chi2
// 16: touches next layer at turn back, so missing hit is added and prevents turning back
//
// no way:
// not enough hits: 0, 9,
// missing seed hits: 3, 11, 12, 13
// endcaps only: 18, 20
//
// other:
// 6: chargino goes to endcaps
// 7: those hits got crazy
// 10: missing hit poorly placed, as it's the first one after the seed. Then cannot extend from there.
// 43: high middle Δφ, then scattered a bit


// "after_L2/4layers/":

// to improve:
// turns back to previous layer: 11, 22
// maybe different s(t) and r(t): 26
// seeds fitting issues: 0
// maybe with endcaps: 1, 21, 24, 29, 36
// missing seed hits: 6, 8 (+ endcaps needed), 12, 14
// special cases:
// - 34: turns back to the same layer, but probably missing hits in the next one
// - 2: this one is a mess and a bless, soo many hits that we have a huge number of possible helices

// "no way" events: 4, 5, 9, 10, 16, 17, 20, 23, 32, 40


Display *display;
shared_ptr<EventSet> events;

map<string,any> filteredPointsOptions = {
  {"title", "Filtered Points"},
	{"binsMin" , 0},
	{"binsMax" , 100000},
  {"markerStyle", 20},
  {"markerSize", 1.0},
  {"color", kYellow}
};

void DrawHitsOrClusters(const shared_ptr<Event> event, int pointsType)
{
  
  vector<shared_ptr<Point>> hitsOrClusters;
  
  map<string,any> drawingOptions = {
    {"markerStyle", (pointsType==2) ? 20 : 22},
    {"markerSize", (pointsType==2) ? 1.0 : 2.0},
  };
  string typeName;
  
  if(pointsType == 0){
    if(!config.drawPionSimHits) return;
    
    hitsOrClusters = event->GetPionSimHits();
    drawingOptions["color"] = kCyan;
    typeName = "Pions hits ";
  }
  else if(pointsType == 1){
    if(!config.drawCharginoSimHits) return;
    
    hitsOrClusters = event->GetCharginoSimHits();
    drawingOptions["color"] = kMagenta;
    typeName = "Charginos hits ";
  }
  else if(pointsType == 2){
    if(!config.drawTrackerClusters) return;
    
    hitsOrClusters = event->GetTrackerClusters();
    drawingOptions["color"] = kYellow;
    typeName = "Tracker clusters ";
  }
  else if(pointsType == 3){
    if(!config.drawPionClusters) return;
    
    hitsOrClusters = event->GetPionClusters();
    pointsProcessor.SetPointsLayers(hitsOrClusters);
    drawingOptions["color"] = kBlue;
    typeName = "Pion clusters ";
  }
  
  map<string, vector<shared_ptr<Point>>> hitsOrClustersBySubDet;
  
  map<int, string> subDetMap = {
    {0,  "PixelBarrel"},
    {1,  "PixelEndcap"},
    {2,  "TIB"},
    {3,  "TOB"},
    {4,  "TID"},
    {5,  "TEC"},
    {6,  "CSC"},
    {7,  "DT"},
    {8,  "RPCBarrel"},
    {9,  "RPCEndcap"},
    {10, "GEM"},
    {11, "ME0"},
    {12, "P2OTB"},
    {13, "P2OTEC"},
    {14, "P1PXB"},
    {15, "P1PXEC"},
    {16, "P2PXB"},
    {17, "P2PXEC"},
    {18, "TimingBarrel"},
    {19, "TimingEndcap"},
    {20, "invalidDet"}
  };
  
  for(auto &[iter, name] : subDetMap){
    hitsOrClustersBySubDet[name] = vector<shared_ptr<Point>>();
  }
  
  for(auto &hit : hitsOrClusters){
    hitsOrClustersBySubDet[hit->GetSubDetName()].push_back(hit);
  }
  
  for(auto &[name, hitsVector] : hitsOrClustersBySubDet){
    if(hitsVector.size() == 0) continue;
    
    drawingOptions["title"] = (typeName+name).c_str();
    display->DrawSimplePoints(hitsVector, drawingOptions);
  }
}

shared_ptr<Event> GetEvent()
{
  EventSet events;
	//  events->LoadEventsFromFiles("/");
//  events.LoadEventsFromFiles(cutLevel);
  events.LoadEventFromFiles(dataType, setIter, iEvent, cutLevel);
	
//  auto event = events.GetEvent(dataType, searchRun, searchLumi, searchEvent);
	auto event = events.At(dataType, setIter, 0);
	
	if(!event){
		cout<<"eventDisplay -- event not found"<<endl;
		exit(0);
	}
	
  JetCut jetCut;
	
	jetCut.SetPt(range<double>(30.0, inf));
	jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
	jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
	
  eventProcessor.ApplyJetCut(event, jetCut);
	
	return event;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  config = ConfigManager(configPath);
  display = new Display();
	
	auto event = GetEvent();
	
	const map<string,any> dedxOptions = {
		{"title", "dE/dx clusters"},
		{"binsMin" , 1},
		{"binsMax" , 10},
		{"nBins" , 5},
		{"markerStyle", 20},
		{"markerSize", 2.0}
	};
	
	display->DrawEvent(event, dedxOptions);
  event->Print();

  // -----------------------------------------------------------------------------------------------------
  // Helix fitting part
  // -----------------------------------------------------------------------------------------------------
  
  cout<<"Preparing hits, track and pion's helix"<<endl;
  
  DrawHitsOrClusters(event, 0); // pion sim hits
  DrawHitsOrClusters(event, 1); // chargino sim hits
  DrawHitsOrClusters(event, 2); // tracker clusters
  DrawHitsOrClusters(event, 3); // pion rec clusters
  
  
	const map<string,any> trueHelixOptions = {
		{"title", "True helix"},
		{"markerStyle", 20},
		{"markerSize", 0.2},
		{"color", kGreen}
	};
	
  auto allSimplePoints = event->GetTrackerClusters();
  
  vector<Helix> truePionHelices = event->GetGenPionHelices();
  
  for(auto &helix : truePionHelices){
    display->DrawHelix(helix,trueHelixOptions);
    //      helix.SetPoints(allSimplePoints);
    //      auto helixPoints = helix.GetPoints();
    //      filteredPointsOptions["title"] = "true helix points";
    //      filteredPointsOptions["color"] = kRed;
    //      display->DrawSimplePoints(helixPoints, filteredPointsOptions);
    
    cout<<"\n\nTrue pion helix:"<<endl;
    helix.Print();
  }
	
  for(auto p = allSimplePoints.begin(); p != allSimplePoints.end();){
    shared_ptr<Point> point = *p;
    if(   point->GetSubDetName() == "TID"
       || point->GetSubDetName() == "TEC"
       || point->GetSubDetName() == "P1PXEC")
      p = allSimplePoints.erase(p);
    else p++;
  }
  
	if(fitHelix){
		cout<<"Fitting best helix"<<endl;
		auto fitter = make_unique<Fitter>();
    
    auto pionClusters = event->GetPionClusters();
    auto pointsByLayer = pointsProcessor.SortByLayer(allSimplePoints);
  
    for(auto &track : event->GetTracks()){
      for(shared_ptr<Point> point : pointsByLayer[track->GetNtrackerLayers()]){
        auto trackPoint = pointsProcessor.GetPointOnTrack(layerR[track->GetNtrackerLayers()], *track, *event->GetVertex());
        
        if(pointsProcessor.distance(make_shared<Point>(trackPoint), point) < 100){
          if(find_if(pionClusters.begin(), pionClusters.end(), [&](const shared_ptr<Point> &p) {return *p == *point;}) == pionClusters.end()){
            pionClusters.push_back(point);
          }
        }
      }
    }
      
    map<string,any> pionClustersOptions = {
      {"title", "Pion hits"},
      {"binsMin" , 0},
      {"binsMax" , 100000},
      {"markerStyle", 20},
      {"markerSize", 2.0},
      {"color", kCyan}
    };

    vector<int> rndIndices = {};
    
    vector<shared_ptr<Point>> pointsNoEndcaps;
    
    for(auto &point : allSimplePoints){
      if(point->GetSubDetName() == "TID" || point->GetSubDetName() == "TEC" || point->GetSubDetName() == "P1PXEC") continue;
      
      if(point->GetSubDetName() != "TIB" && point->GetSubDetName() != "TOB" && point->GetSubDetName() != "P1PXB"){
        cout<<"Weird detector:"<<point->GetSubDetName()<<endl;
      }
      
      pointsNoEndcaps.push_back(point);
    }
    pointsProcessor.SetPointsLayers(pionClusters);
    display->DrawSimplePoints(pionHitsOnly ? pionClusters : pointsNoEndcaps, pionClustersOptions);
    
    auto start = now();
    vector<Helix> fittedHelices;
    
    for(auto &track : event->GetTracks()){
      bool trackLayersOk = false;
      bool chargeOk = false;
      
      for(auto chargino : event->GetGenCharginoTracks()){
        if(fabs(track->GetEta() - chargino.GetEta()) < 0.1){
          if(track->GetNtrackerLayers() == chargino.GetNtrackerLayers()) trackLayersOk = true;
          if(track->GetCharge() == chargino.GetCharge()) chargeOk = true;
        }
      }
      cout<<"\n\nTrack layers ok: "<<(trackLayersOk ? "YES" : "NO")<<endl;
      cout<<"\n\nTrack charge ok: "<<(chargeOk ? "YES" : "NO")<<endl;
      cout<<"Fitting for track: "; track->Print();
      
      vector<Helix> helices = fitter->FitHelices(pionHitsOnly ? pionClusters : pointsNoEndcaps, *track, *event->GetVertex());
      fittedHelices.insert(fittedHelices.end(), helices.begin(), helices.end());
    }
    auto end = now();
    
    cout<<"Fitting time: "<<duration(start, end)<<endl;

    map<string,any> bestHelixOptions = {
      {"title", "Best helix"},
      {"markerStyle", 20},
      {"markerSize", 0.2},
      {"color", kRed}
    };
    
    map<string,any> helixPointsOptions = {
      {"title", "Helix vertex"},
      {"binsMin" , 0},
      {"binsMax" , 100000},
      {"markerStyle", 20},
      {"markerSize", 2.0},
      {"color", kYellow}
    };
    
    for(auto helix : fittedHelices){
      cout<<endl; helix.Print();
      
      bestHelixOptions["title"] = ("Helix "+to_string(helix.GetSeedID())).c_str();
      
      display->DrawShrinkingHelix(helix, bestHelixOptions);
      display->DrawSimplePoints(helix.GetPoints(), helixPointsOptions);
    }
  }
  
     
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
