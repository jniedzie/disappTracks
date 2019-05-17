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
string cutLevel = "after_L1/4layers/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 0;

// 6  (q+, vz+, pz-) OK
// 10 (q+, vz+, pz-) OK - RECO
// 11 (q+, vz-, pz-) OK [strong breaking]
// 13 (q-, vz-, pz-) OK
// 14 (q-, vz+, pz+) OK
// 15 (q+, vz+, pz+) OK
// 18 (q-, vz+, pz+), (q+, vz+, pz+) OK [nice two helices]
// 19 (q-, vz-, pz-), (q+, vz-, pz-) OK [crazy stuff, probably 3rd soft particle there...]
// 20 (q+, vz+, pz+) OK [very high p_z]
// 21 (q-, vz+, pz-), (q+, vz+, pz+) OK [two big helices, one strongly breaking]
// 23 (q+, vz+, pz+) OK

bool injectPion = false;
bool fitHelix = true;

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
	events.LoadEventsFromFiles(cutLevel);
	
//  auto event = events.GetEvent(dataType, searchRun, searchLumi, searchEvent);
	auto event = events.At(dataType, setIter, iEvent);
	
	if(!event){
		cout<<"eventDisplay -- event not found"<<endl;
		exit(0);
	}
	
  JetCut jetCut;
	
	jetCut.SetPt(range<double>(30.0, inf));
	jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
	jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
	
  eventProcessor.ApplyJetCut(event, jetCut);
	
  event->LoadAdditionalInfo();
  
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
  
	if(injectPion){
		// Draw decay point to make sure that it's correctly located
    vector<shared_ptr<Point>> decayPoint;
    decayPoint.push_back(make_shared<Point>(track->GetDecayPoint().GetX(),
                                            track->GetDecayPoint().GetY(),
                                            track->GetDecayPoint().GetZ()));
		
		const map<string,any> decayPointOptions = {
			{"title", "Decay Point"},
			{"markerStyle", 22},
			{"markerSize", 2.0},
			{"color", kGreen}
		};
		
		display->DrawSimplePoints(decayPoint, decayPointOptions);
		
		// Draw true pion helix
    Helix pionHelix;
    helixProcessor.GetRandomPionHelix(track, pionHelix);
		display->DrawHelix(pionHelix,trueHelixOptions);
		cout<<"\n\nInjected pion helix:"<<endl;
		pionHelix.Print();
		
		// Calculate and draw points along the helix that hit the silicon
    auto pionPoints = helixProcessor.GetPointsHittingSilicon(pionHelix);
		for(auto &p : pionPoints){p->SetIsPionHit(true);}
		pionHelix.SetPoints(pionPoints);
		
		const map<string,any> pionPointsOptions = {
			{"title", "Pion points"},
			{"markerStyle", 21},
			{"markerSize", 1.6},
			{"color", kMagenta}
		};
		
		display->DrawSimplePoints(pionPoints, pionPointsOptions);
		
		// inject hits from pion into all points in the tracker
		allSimplePoints.insert(allSimplePoints.end(), pionPoints.begin(), pionPoints.end());
	}
	else{
    vector<Helix> truePionHelices;
    event->GetGenPionHelices(truePionHelices);
		
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
	}
	
  for(auto p = allSimplePoints.begin(); p != allSimplePoints.end();){
    shared_ptr<Point> point = *p;
    if(point->GetSubDetName() == "TID" || point->GetSubDetName() == "TEC") p = allSimplePoints.erase(p);
    else p++;
  }
  
	if(fitHelix){
		cout<<"Fitting best helix"<<endl;
		auto fitter = make_unique<Fitter>();
		
    auto pionClusters = event->GetPionClusters();

    for(auto &p : allSimplePoints){
      int layer = -1;
      double minDist = inf;
      double pointR = sqrt(pow(p->GetX(), 2) + pow(p->GetY(), 2));
      
      for(int iLayer=0; iLayer<nLayers; iLayer++){
        double pointLayerDist = fabs(layerR[iLayer] - pointR);
        
        if(pointLayerDist < minDist){
          minDist = pointLayerDist;
          layer = iLayer;
        }
      }
      
      if(layer == track->GetNtrackerLayers()){
        Point trackPoint(layerR[layer] * cos(track->GetPhi())    + 10*event->GetVertex()->GetX(),
                         layerR[layer] * sin(track->GetPhi())    + 10*event->GetVertex()->GetY(),
                         layerR[layer] / tan(track->GetTheta())  + 10*event->GetVertex()->GetZ());

        if(pointsProcessor.distance(make_shared<Point>(trackPoint), p) < 100){
          pionClusters.push_back(p);
        }
      }
    }
    
//    map<string,any> pionClustersOptions = {
//      {"title", "Pion clusters"},
//      {"binsMin" , 0},
//      {"binsMax" , 100000},
//      {"markerStyle", 20},
//      {"markerSize", 2.0},
//      {"color", kCyan}
//    };
//
//    display->DrawSimplePoints(pionClusters, pionClustersOptions);
    
    vector<int> rndIndices = { };
    
    // Turn this on to inject some noise
    cout<<endl;
    for(int i=0;i<0;i++){
//      int r = rndIndices[i];
      int r = RandInt(0, (int)allSimplePoints.size()-1);
      cout<<r<<",";
      pionClusters.insert(pionClusters.end(),allSimplePoints[r]);
    }
    cout<<endl;
    
//    display->DrawSimplePoints(pionClusters, filteredPointsOptions);
    
    vector<Helix> fittedHelices = fitter->FitHelices(pionClusters, *track, *event->GetVertex());
    		
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
    
    for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){

      bestHelixOptions["title"] = ("Helix "+to_string(fittedHelices[iHelix].uniqueID)).c_str();
      display->DrawShrinkingHelix(fittedHelices[iHelix], bestHelixOptions);
      fittedHelices[iHelix].Print();
      
      helixVertexOptions["markerStyle"] = 20;
      vector<shared_ptr<Point>> helixVertex = fittedHelices[iHelix].GetPoints();
      display->DrawSimplePoints(helixVertex, helixVertexOptions);
    }
  }
  
     
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
