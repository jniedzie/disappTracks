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
string cutLevel = "after_L2/4layers/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 13;

bool injectPion = false;
bool fitHelix = true;

// "after_L2/4layers/":

// 0 - bad hits
// 1 - bad hits
// 2 - missing second hit
// 4 - not enough hits
// 5 - scattered
// 6 - missing second hit
// 8 - high p_z (820)
// 9 - scattered
// 10 - not enough hits
// 12 - scattered
// 14 - missing (very scattered) first pion hit
// 15 - almost ok
// 16 - no hits
// 17 - bad hits
// 20 - scattered
// 21 - high p_z (1700 MeV)
// 23 - no hits
// 24 - high p_z
// 25 - high p_z
// 29 - high p_z
// 32 - no hits
// 35 - missing first hit
// 36 - high p_z/small R
// 40 - no hits
// 42 - high p_z


// "after_L1/4layers/":
// 0 - charge mismatch
// 1 - no gen pion
// 2 - no valid seeds (to be checked, pion turns back from chargino)
// 3 - track misreco
// 4 - ... Z+, Q+, only 4 points, looks like it hits fitter limits
// 5 - only 3 pion hits (2 of them in the same layer)
// 6 - no valid seeds, due to high pz (650 MeV)
// 7 - charge mismatch + only 2 layers of strip hits (then some endcap hits far away)
// 8 - ... no valid seeds, but actually this one doesn't look bad - to investirage
// 9 - no gen pion
// 10 - ... Z+, Q- (wrong helix bending)
// 11 - track misreco
// 12 - no gen pion from the reconstructed track
// 13 - ... Z+, Q+
// 14 - missing (very displaced) first strip hit
// 15 - charge mismatch
// 16 - no pion hits
// 17 - too high pz (790 MeV)
// 18 - track misreco
// 19 - too high pz (910 MeV)
// 20 - too high pz (630 MeV)
// 21 - high track eta, high vertex z (510 mm)
// 22 - missing first hit (?)
// 23 - no gen pion
// 24 - ... Z-, Q+
// 25 - track misreco
// 26 - track misreco
// 27 - OK, Z-, Q-
// 28 - only a few pion hits, scattered a lot
// 29 - track misreco
// 30 - ... Z+, Q+
// __31 - OK, Z+, Q+
// 32 - shitty hits
// 33 - no gen pion
// 34 - high pz (570 MeV)
// 35 - ... Z+, Q+, perfect case, idk why doesn't work
// 36 - ... Z-, Q-, very nice, should work
// 37 - ... no seeds. Just a few hits, but should work...
// 38 - too high pz (590 MeV)
// 39 - shitty hits
// 40 - ... Q+, Z-, no seeds, but looks good...

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
	shared_ptr<Track> track = event->GetTrack(0);
	
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
  
//    for(auto &point : pointsByLayer[track->GetNtrackerLayers()]){
//      auto trackPoint = pointsProcessor.GetPointOnTrack(layerR[track->GetNtrackerLayers()], *track, *event->GetVertex());
//
//      if(pointsProcessor.distance(make_shared<Point>(trackPoint), point) < 100){
//        pionClusters.push_back(point);
//      }
//    }
    
    map<string,any> pionClustersOptions = {
      {"title", "Layer 4"},
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
    
    for(int i=0; i<0; i++){
//      int r = rndIndices[i];
      int r = RandInt(0, (int)allSimplePoints.size()-1);
      auto point = allSimplePoints[r];

      if(find(pionClusters.begin(), pionClusters.end(), point) == pionClusters.end()) pionClusters.push_back(point);
      else                                                                            i--;
    }
    cout<<endl;
    
//    pionClusters = pointsProcessor.FilterNearbyPoints(pionClusters, 50);
    
    auto start = now();
//    vector<Helix> fittedHelices = fitter->FitHelices(allSimplePoints, *track, *event->GetVertex());

    display->DrawSimplePoints(pointsNoEndcaps, pionClustersOptions);
    vector<Helix> fittedHelices = fitter->FitHelices(pointsNoEndcaps, *track, *event->GetVertex());
    
//    pointsProcessor.SetPointsLayers(pionClusters);
//    display->DrawSimplePoints(pionClusters, pionClustersOptions);
//    vector<Helix> fittedHelices = fitter->FitHelices(pionClusters, *track, *event->GetVertex());
    
    auto end = now();
    
    cout<<"Fitting time: "<<duration(start, end)<<endl;
    
    
    map<string,any> bestHelixOptions = {
      {"title", "Best helix"},
      {"markerStyle", 20},
      {"markerSize", 0.2},
      {"color", kRed}
    };
    
    map<string,any> helixVertexOptions = {
      {"title", "Helix vertex"},
      {"binsMin" , 0},
      {"binsMax" , 100000},
      {"markerStyle", 20},
      {"markerSize", 2.0},
      {"color", kYellow}
    };
    
    for(auto helix : fittedHelices){

      auto helixPoints = helix.GetPoints();
      
      int nPionPoints = 0;
      for(auto &pionPoint : pionClusters){
        for(auto &helixPoint : helixPoints){
          if(*pionPoint == *helixPoint) nPionPoints++;
        }
      }
      
      int nFakePoints = (int)helixPoints.size()-1 - nPionPoints;
      
      cout<<"Helix pion points:"<<nPionPoints<<endl;
      cout<<"\tfake points:"<<nFakePoints<<endl;
      
      bestHelixOptions["title"] = ("Helix "+to_string(helix.GetUniqueID())).c_str();
      display->DrawShrinkingHelix(helix, bestHelixOptions);
      helix.Print();
      
      helixVertexOptions["markerStyle"] = 20;
      vector<shared_ptr<Point>> helixVertex = helix.GetPoints();
      display->DrawSimplePoints(helixVertex, helixVertexOptions);
    }
  }
  
     
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
