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
string cutLevel = "after_L0/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 10;

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
    if(!config->drawPionSimHits) return;
    
    hitsOrClusters = event->GetPionSimHits();
    drawingOptions["color"] = kCyan;
    typeName = "Pions hits ";
  }
  else if(pointsType == 1){
    if(!config->drawCharginoSimHits) return;
    
    hitsOrClusters = event->GetCharginoSimHits();
    drawingOptions["color"] = kMagenta;
    typeName = "Charginos hits ";
  }
  else if(pointsType == 2){
    if(!config->drawTrackerClusters) return;
    
    hitsOrClusters = event->GetTrackerClusters();
    drawingOptions["color"] = kYellow;
    typeName = "Tracker clusters ";
  }
  else if(pointsType == 3){
    if(!config->drawPionClusters) return;
    
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
	events = make_shared<EventSet>();
	//  events->LoadEventsFromFiles("/");
	events->LoadEventsFromFiles(cutLevel);
	
//  auto event = events->GetEvent(dataType, searchRun, searchLumi, searchEvent);
	auto event = events->At(dataType, setIter, iEvent);
	
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
	
  event->LoadAdditionalInfo();
  
	return event;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  config = make_unique<ConfigManager>(configPath);
  display = new Display();
  auto helixProcessor = make_unique<HelixProcessor>();
  auto trackProcessor = make_unique<TrackProcessor>();
	
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

  // ------------------------------------------------------------------------------------------------------------
  // Helix fitting part
  // ------------------------------------------------------------------------------------------------------------
  
  cout<<"Preparing hits, track and pion's helix"<<endl;
  
  DrawHitsOrClusters(event, 0); // pions
  DrawHitsOrClusters(event, 1); // charginos
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
    decayPoint.push_back(make_shared<Point>(track->GetDecayPoint()->GetX(),
                                            track->GetDecayPoint()->GetY(),
                                            track->GetDecayPoint()->GetZ()));
		
		const map<string,any> decayPointOptions = {
			{"title", "Decay Point"},
			{"markerStyle", 22},
			{"markerSize", 2.0},
			{"color", kGreen}
		};
		
		display->DrawSimplePoints(decayPoint, decayPointOptions);
		
		// Draw true pion helix
		auto pionHelix = helixProcessor->GetRandomPionHelix(track);
		display->DrawHelix(pionHelix,trueHelixOptions);
		cout<<"\n\nInjected pion helix:"<<endl;
		pionHelix->Print();
		
		// Calculate and draw points along the helix that hit the silicon
		auto pionPoints = helixProcessor->GetPointsHittingSilicon(pionHelix);
		for(auto &p : pionPoints){p->SetIsPionHit(true);}
		pionHelix->SetPoints(pionPoints);
		
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
		const vector<unique_ptr<Helix>> &truePionHelices = event->GetGenPionHelices();
		
		for(auto &helix : truePionHelices){
			display->DrawHelix(helix,trueHelixOptions);
			helix->SetPoints(allSimplePoints);
			auto helixPoints = helix->GetPoints();
			filteredPointsOptions["title"] = "true helix points";
			filteredPointsOptions["color"] = kRed;
			display->DrawSimplePoints(helixPoints, filteredPointsOptions);
			
			cout<<"\n\nTrue pion helix:"<<endl;
			helix->Print();
		}
	}
	
	if(fitHelix){
		cout<<"Fitting best helix"<<endl;
		auto fitter = make_unique<Fitter>();
		
    auto bestHelix = fitter->FitHelix(event->GetPionSimHits(), track, event->GetVertex());
//    auto bestHelix = fitter->GetBestFittingHelix(allSimplePoints, track, event->GetVertex());
		
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
		}
    else{
      cout<<"\n\nCould not fit any helix..."<<endl;
    }
  }
	
  /*
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  TH1D *chargeHistClusters = new TH1D("chargeHistClusters","chargeHistClusters",100,0,1000);
  TH1D *chargeHistPion = new TH1D("chargeHistPion","chargeHistPion",100,0,1000);
  chargeHistClusters->Sumw2();
  chargeHistPion->Sumw2();
  chargeHistClusters->SetLineColor(kRed);
  chargeHistPion->SetLineColor(kBlue);
  
  for(int i=0;i<30;i++){
    auto event = events->At(dataType, setIter, i);
    event->LoadAdditionalInfo();
    
    auto clusters = event->GetTrackerClusters();
    auto pionHits = event->GetPionClusters();
    
    for(auto cluster : *clusters){
      chargeHistClusters->Fill(cluster.GetValue());
    }
    for(auto hit : *pionHits){
      chargeHistPion->Fill(hit.GetValue());
    }
  }
  
  chargeHistClusters->DrawNormalized();
  chargeHistPion->DrawNormalized("same");
  c1->Update();
  
  cout<<"tracker mean:"<<chargeHistClusters->GetMean()<<" +/- "<<chargeHistClusters->GetMeanError()<<endl;
  cout<<"pion mean:"<<chargeHistPion->GetMean()<<" +/- "<<chargeHistPion->GetMeanError()<<endl;
  
  */
  
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
