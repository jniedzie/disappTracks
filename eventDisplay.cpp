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

vector<vector<double>> charginoPoints = {};
vector<vector<double>> pionPoints = {};

// event 8 (2646)
//vector<vector<double>> charginoPoints = {
//  // low hits
//  {2.459, -1.452, 6.62},
//  {6.839, -1.409, 6.624},
//  {10.95, -1.582, 6.628},
//  {2.459, -1.452, 6.62},
//  // high hits
//  // none
//};
//
//vector<vector<double>> pionPoints = {
//  // low hits
//  {2.718, 1.614, 6.579},
//  {6.774, 1.699, 6.529},
//  {10.64, 1.468, -0.08096},
//  {14.94, 6.188, -0.1679},
//  // high hits
//  // none
//};

// event 12 (2602)
//vector<vector<double>> charginoPoints = {
//  // low hits
//  {-0.008219, -3.375, 13.26},
//  {0.1146, -7.249, 13.29},
//  {-0.02884, -11.34, 13.31},
//  {-0.02884, -11.34, 13.31},
//  {0.1146, -7.249, 13.29},
//  // high hits
//  // none
//};
//
//vector<vector<double>> pionPoints = {
//  {-0.01882, -16.44, 19.94},
//};

// event 21 (1358)
//vector<vector<double>> charginoPoints = {
//  {2.715, 1.623, -0.06231},
//  {6.449, 1.522, -0.1197},
//  {5.373, 4.458, -0.1305},
//  {9.797, 4.397, -6.776},
//  {14.65, 6.018, -6.867},
//  {9.797, 4.397, -6.776}
//};
//vector<vector<double>> pionPoints = {
//};

// event 26 (2089)
//vector<vector<double>> charginoPoints = {
//  {2.457, -1.46, 6.618},
//  {2.764, -1.539, 6.617},
//  {5.207, -4.1, 6.61},
//  {10.06, -4.603, 6.565},
//  {13.14, -8.838, 6.538},
//  {5.207, -4.1, 6.611},
//};
//
//vector<vector<double>> pionPoints = {
//};


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
  auto allSimplePoints = event->GetTrackerHits();
  display->DrawSimplePoints(allSimplePoints, filteredPointsOptions);
	
  
  auto charginoSimHits = make_shared<vector<Point>>();
  for(int i=0;i<charginoPoints.size();i++){
    charginoSimHits->push_back(Point(10*charginoPoints[i][0],
                                     10*charginoPoints[i][1],
                                     10*charginoPoints[i][2]));
  }
  const map<string,any> charginoSimHitsOptions = {
    {"title", "Chargino sim hits"},
    {"markerStyle", 22},
    {"markerSize", 2.0},
    {"color", kMagenta}
  };
  display->DrawSimplePoints(charginoSimHits, charginoSimHitsOptions);
  
  auto pionSimHits = make_shared<vector<Point>>();
  for(int i=0;i<pionPoints.size();i++){
    pionSimHits->push_back(Point(10*pionPoints[i][0],
                                 10*pionPoints[i][1],
                                 10*pionPoints[i][2]));
  }
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
		auto truePionHelices = event->GetTruePionHelices();
		
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
