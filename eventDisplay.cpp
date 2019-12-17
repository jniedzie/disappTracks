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
int setIter = kTaggerSignalNoPU;
int iEvent = 6;

// endcap track: 7, 18, 25, 40
// endcap hits: 11, 24, 25, 36, 40
// 17: nice endcap hits, but missing seed hits...

unique_ptr<Display> display;
shared_ptr<EventSet> events;

map<string,any> filteredPointsOptions = {
  {"title", "Filtered Points"},
	{"binsMin" , 0},
	{"binsMax" , 100000},
  {"markerStyle", 20},
  {"markerSize", 1.0},
  {"color", kMagenta+3}
};

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

/**
 Returns event specified by global variables `dataType`, `setIter`, `iEvent` and `cutLevel`.
 In case event could not be found, execution stops.
 */
shared_ptr<Event> GetEvent()
{
  EventSet events; events.LoadEventsFromFiles(dataType, setIter, cutLevel, iEvent);
  auto event = events.At(dataType, setIter, 2017, 0);
  
  if(!event){
    cout<<"eventDisplay -- event not found"<<endl;
    exit(0);
  }
  return event;
}

/**
 Draws pion/chargino sim hits, all tracker clusters and clusters associated with the pion,
 according to options specified in the config
 */
void DrawHitsAndClusters(const shared_ptr<Event> event)
{
  map<string,any> drawingOptions = {
    {"markerStyle", 22},
    {"markerSize" , 2.0},
  };
  string typeName;
  
  auto drawFunction = [&](Points points){
    pointsProcessor.SetPointsLayers(points);
    pointsProcessor.SetPointsDisks(points);
    
    map<string, Points> hitsOrClustersBySubDet;
    for(auto &[iter, name] : subDetMap) hitsOrClustersBySubDet[name] = Points();
    
    for(auto &point : points){
      hitsOrClustersBySubDet[point->GetSubDetName()].push_back(point);
    }
    for(auto &[name, hitsVector] : hitsOrClustersBySubDet){
      if(hitsVector.size() == 0) continue;
      
      drawingOptions["title"] = (typeName+name).c_str();
      display->DrawPoints(hitsVector, drawingOptions);
    }
  };
  
  if(config.params["draw_pion_simhits"]){
    drawingOptions["color"] = kCyan;
    typeName = "Pions hits ";
    drawFunction(event->GetPionSimHits());
  }
  if(config.params["draw_chargino_simhits"]){
    
    drawingOptions["color"] = kMagenta;
    typeName = "Charginos hits ";
    drawFunction(event->GetCharginoSimHits());
  }
  if(config.params["draw_tracker_clusters"]){
    drawingOptions["markerStyle"] = 20;
    drawingOptions["markerSize"] = 1.0;
    drawingOptions["color"] = (EColor)(kMagenta+2);
    typeName = "Tracker clusters ";
    drawFunction(event->GetTrackerClusters());
  }
  if(config.params["draw_pion_clusters"]){
    drawingOptions["color"] = kBlue;
    typeName = "Pion clusters ";
    drawFunction(event->GetPionClusters());
  }
}

/**
 Defines drawing options, draws and prints provided event.
 */
void DrawEvent(const shared_ptr<Event> &event)
{
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
}

/**
 Defines drawing options, prints and draws any gen-level pion helix found in the provided event.
 */
void DrawTruePionHelices(const shared_ptr<Event> &event)
{
  const map<string,any> trueHelixOptions = {
    {"title", "True helix"},
    {"markerStyle", 20},
    {"markerSize", 0.4},
    {"color", kGreen}
  };
  
  Helices truePionHelices = event->GetGenPionHelices();
  
  for(auto &helix : truePionHelices){
    display->DrawHelix(helix,trueHelixOptions);
    cout<<"\n\nTrue pion helix:"<<endl; helix.Print();
  }
}

/**
 Returns pion strip clusters + pixel clusters that were closer than 10 cm from the track.
 */
Points GetPionClusters(const shared_ptr<Event> &event)
{
  Points pionClusters = event->GetPionClusters();
  Points eventClusters = event->GetClusters();
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
  return pionClusters;
}

/**
 Returns points for fitting, which can be either pion clusters only, all tracker clusters (with or without
 end-cap clusters), or removing pion clusters form from the collection of tracker clusters, according to settings
 
 */
Points GetPointsToFit(const shared_ptr<Event> &event)
{
  Points pointsToFit;
  
  if(config.params["fit_pion_clusters_only"]) pointsToFit = GetPionClusters(event);
  else pointsToFit = event->GetClusters();
  
  pointsProcessor.SetPointsLayers(pointsToFit);
  pointsProcessor.SetPointsDisks(pointsToFit);
  
  return pointsToFit;
}

/**
 Defines drawing options and draws points used for fitting.
 */
void DrawFittingPoints(const Points &points)
{
  map<string,any> fittingPointsOptions = {
    {"title", "Fitting points"},
    {"color", (EColor)(kMagenta)},
    {"emphasis", false}
  };
  
  display->DrawPoints(points, fittingPointsOptions);
}

/**
 Defines drawing options and draws fitted helices and clusters that belong to them.
 */
void DrawFittedHelices(const Helices &helices)
{
  map<string,any> helixCurveOptions = {
    {"markerStyle", 20},
    {"markerSize", 0.4},
    {"color", kRed}
  };
  
  map<string,any> helixPointsOptions = {
    {"title", "Fitted helix points"},
    {"color", (EColor)(kGreen)},
    {"emphasis", true}
  };
  
  for(auto helix : helices){
    cout<<endl; helix.Print();
    helixCurveOptions["title"] = ("Fitted helix "+to_string(helix.GetUniqueID())).c_str();
    display->DrawShrinkingHelix(helix, helixCurveOptions);
    display->DrawPoints(helix.GetPoints(), helixPointsOptions);
  }
}

/**
 Prints information about a track, including whether or not number of layers and charge was correctly reconstructed.
 */
void PrintTrack(const shared_ptr<Track> &track, const shared_ptr<Event> &event)
{
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
}

/**
 The program execution starting point.
 */
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  display = make_unique<Display>();
	
	auto event = GetEvent();
  DrawEvent(event);
  DrawHitsAndClusters(event);
  DrawTruePionHelices(event);
  
	if(config.params["fit_helices"]){
		cout<<"Fitting best helix"<<endl;
		auto fitter = make_unique<Fitter>();
    
    Points pointsToFit = GetPointsToFit(event);
    DrawFittingPoints(pointsToFit);
    
    auto start = now();
    Helices fittedHelices;
    
    for(auto &track : event->GetTracks()){
      PrintTrack(track, event);
      Helices helices = fitter->FitHelices(pointsToFit, *track, *event->GetVertex());
      fittedHelices.insert(fittedHelices.end(), helices.begin(), helices.end());
    }
    auto end = now();
    cout<<"Fitting time: "<<duration(start, end)<<endl;

    DrawFittedHelices(fittedHelices);
  }
  
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
