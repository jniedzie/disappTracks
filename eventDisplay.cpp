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
string cutLevel = "";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 26; // 8, 12, 20, 21, 26

bool injectPion = false;
bool fitHelix = false;

Display *display;

map<string,any> filteredPointsOptions = {
  {"title", "Filtered Points"},
	{"binsMin" , 0},
	{"binsMax" , 100000},
  {"markerStyle", 20},
  {"markerSize", 1.0},
  {"color", kYellow}
};

shared_ptr<Event> GetEvent()
{
	auto events = make_shared<EventSet>();
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
  auto allSimplePoints = event->GetTrackerClusters();
  display->DrawSimplePoints(allSimplePoints, filteredPointsOptions);
	
  
  auto charginoSimHits = event->GetCharginoSimHits();
  
  const map<string,any> charginoSimHitsOptions = {
    {"title", "Chargino sim hits"},
    {"markerStyle", 22},
    {"markerSize", 2.0},
    {"color", kMagenta}
  };
  display->DrawSimplePoints(charginoSimHits, charginoSimHitsOptions);
  
  auto pionSimHits = event->GetPionSimHits();
  
  const map<string,any> pionSimHitsOptions = {
    {"title", "Pion sim hits"},
    {"markerStyle", 22},
    {"markerSize", 2.0},
    {"color", kCyan}
  };
  display->DrawSimplePoints(pionSimHits, pionSimHitsOptions);
  
   
	const map<string,any> trueHelixOptions = {
		{"title", "True helix"},
		{"markerStyle", 20},
		{"markerSize", 0.2},
		{"color", kGreen}
	};
	
	if(injectPion){
		// Draw decay point to make sure that it's correctly located
		auto decayPoint = make_shared<vector<Point>>();
		decayPoint->push_back(Point(track->GetDecayPoint()->GetX(),
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
		for(auto &p : *pionPoints){p.SetIsPionHit(true);}
		pionHelix->SetPoints(pionPoints);
		
		const map<string,any> pionPointsOptions = {
			{"title", "Pion points"},
			{"markerStyle", 21},
			{"markerSize", 1.6},
			{"color", kMagenta}
		};
		
		display->DrawSimplePoints(pionPoints, pionPointsOptions);
		
		// inject hits from pion into all points in the tracker
		allSimplePoints->insert(allSimplePoints->end(), pionPoints->begin(), pionPoints->end());
	}
	else{
		auto truePionHelices = event->GetGenPionHelices();
		
		for(auto &helix : *truePionHelices){
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
		auto bestHelix = fitter->GetBestFittingHelix(allSimplePoints, track, event->GetVertex());
		
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
	
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
