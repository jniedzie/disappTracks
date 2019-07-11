#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "TrackProcessor.hpp"

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L1/all/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 36;

// endcap track: 7, 18, 25, 40
// endcap hits: 11, 24, 25, 36, 40
// 17: nice endcap hits, but missing seed hits...

bool fitHelix = true;

bool pionHitsOnly = false;
bool removePionClusters = false;
bool removeEncapClusters = false;

unique_ptr<Display> display;
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
    if(!config.params["draw_pion_simhits"]) return;
    
    hitsOrClusters = event->GetPionSimHits();
    drawingOptions["color"] = kCyan;
    typeName = "Pions hits ";
  }
  else if(pointsType == 1){
    if(!config.params["draw_chargino_simhits"]) return;
    
    hitsOrClusters = event->GetCharginoSimHits();
    drawingOptions["color"] = kMagenta;
    typeName = "Charginos hits ";
  }
  else if(pointsType == 2){
    if(!config.params["draw_tracker_clusters"]) return;
    
    hitsOrClusters = event->GetTrackerClusters();
    drawingOptions["color"] = kYellow;
    typeName = "Tracker clusters ";
  }
  else if(pointsType == 3){
    if(!config.params["draw_pion_clusters"]) return;
    
    hitsOrClusters = event->GetPionClusters();
    drawingOptions["color"] = kBlue;
    typeName = "Pion clusters ";
  }
  pointsProcessor.SetPointsLayers(hitsOrClusters);
  pointsProcessor.SetPointsDisks(hitsOrClusters);

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
    display->DrawPoints(hitsVector, drawingOptions);
  }
}

shared_ptr<Event> GetEvent()
{
  EventSet events;
  events.LoadEventFromFiles(dataType, setIter, iEvent, cutLevel);
	
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
  display = make_unique<Display>();
	
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
	
  Helices truePionHelices = event->GetGenPionHelices();
  
  for(auto &helix : truePionHelices){
    display->DrawHelix(helix,trueHelixOptions);
    cout<<"\n\nTrue pion helix:"<<endl; helix.Print();
  }
  
	if(fitHelix){
		cout<<"Fitting best helix"<<endl;
		auto fitter = make_unique<Fitter>();
    
    Points pionClusters = event->GetPionClusters();
    Points eventClusters = event->GetClusters(removePionClusters, removeEncapClusters);
    vector<Points> clustersByLayer = pointsProcessor.SortByLayer(eventClusters);
  
    for(auto &track : event->GetTracks()){
      for(shared_ptr<Point> point : clustersByLayer[track->GetNtrackerLayers()]){
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

    pointsProcessor.SetPointsLayers(pionHitsOnly ? pionClusters : eventClusters);
    pointsProcessor.SetPointsDisks(pionHitsOnly ? pionClusters : eventClusters);
    display->DrawPoints(pionHitsOnly ? pionClusters : eventClusters, pionClustersOptions);
    
    auto start = now();
    Helices fittedHelices;
    
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
      
      Helices helices = fitter->FitHelices(pionHitsOnly ? pionClusters : eventClusters, *track, *event->GetVertex());
      
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
      bestHelixOptions["title"] = ("Helix "+to_string(helix.GetUniqueID())).c_str();
      display->DrawShrinkingHelix(helix, bestHelixOptions);
      display->DrawPoints(helix.GetPoints(), helixPointsOptions);
    }
  }
  
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
