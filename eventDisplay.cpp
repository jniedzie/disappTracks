#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "FitterConfig.hpp"

Display *display;
shared_ptr<FitterConfig> config;

bool showStipClusters = false;

// assumptions about the pion
double decayR = 140; // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double pionCharge = 1;

unique_ptr<Point> pionVector = make_unique<Point>(115,115,115); // Total momentum ~200 MeV

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

vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber);

double trackTheta, trackPhi;

unique_ptr<Helix> GetBestFittingHelix(vector<Point> allSimplePoints)
{
  // Prepare 2D projections in XY
  vector<Point> points2D;
  vector<vector<Point>> pointsByLine = Point::SplitPointsIntoLines(allSimplePoints, config->GetLinesTolerance());
  
  for(vector<Point> line : pointsByLine){
    if(line.size() >= config->GetMinPointsAlongZ()){
      points2D.push_back(Point(line));
    }
  }
  
  // Create fitter to fit circles to 2D distribution
  unique_ptr<Fitter> fitter = unique_ptr<Fitter>(new Fitter(3));
  double helixThickness = config->GetHelixThickness();
  double minL = config->GetMinL();
  double maxL = config->GetMaxL();
  double minPx = config->GetMinPx();
  double maxPx = config->GetMaxPx();
  double minPy = config->GetMinPy();
  double maxPy = config->GetMaxPy();
  fitter->SetParameter(0, "L", (maxL+minL)/2., minL-helixThickness, maxL+helixThickness);
  fitter->SetParameter(1, "px", (maxPx-minPx)/2., minPx, maxPx);
  fitter->SetParameter(2, "py", (maxPy-minPy)/2., minPy, maxPy);
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  cout<<"Looking for candidate circles"<<endl;
  int nPoints = (int)points2D.size();
  double circleThickness = config->GetCircleThickness();
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){
        
        auto chi2Function = [&](const double *par) {
          double f = 0;
          
          double L = par[0];
          double px = par[1];
          double py = par[2];
          
          double x0 = L*sin(trackTheta)*cos(trackPhi);
          double y0 = L*sin(trackTheta)*sin(trackPhi);
          double z0 = L*cos(trackTheta);
          
          unique_ptr<Point> decayPoint  = make_unique<Point>(x0,y0,z0);
          unique_ptr<Point> momentum    = make_unique<Point>(px,py,0);
          Circle circle(decayPoint, momentum, pionCharge, circleThickness);
          
          f  = pow(circle.GetDistanceToPoint(points2D[i]),2);
          f += pow(circle.GetDistanceToPoint(points2D[j]),2);
          f += pow(circle.GetDistanceToPoint(points2D[k]),2);
          
          return f;
        };
        fitter->SetFitFunction(chi2Function);
        
        if(fitter->RunFitting()) {
          auto result = fitter->GetResult();
          
          double L = result.GetParams()[0];
          double px = result.GetParams()[1];
          double py = result.GetParams()[2];
          
          double x0 = L*sin(trackTheta)*cos(trackPhi);
          double y0 = L*sin(trackTheta)*sin(trackPhi);
          double z0 = L*cos(trackTheta);
          
          unique_ptr<Point> decayPoint = make_unique<Point>(x0,y0,z0);
          unique_ptr<Point> momentum = make_unique<Point>(px,py,0);
          unique_ptr<Circle> circle = make_unique<Circle>(decayPoint, momentum, pionCharge, circleThickness);
          
          int nPoints=0;
          for(Point p : points2D){
            if(circle->GetDistanceToPoint(p) < circleThickness) nPoints++;
          }
          if(nPoints > 3){
            circle->SetPoints(allSimplePoints);
            circles.push_back(move(circle));
          }
        }
      }
    }
  }
  cout<<"Found "<<circles.size()<<" circles"<<endl;
  
  unique_ptr<Helix> bestHelix = nullptr;
  int maxNregularPoints = 0;
  double maxFractionRegularPoints = 0;
  cout<<"Scanning candidate circles to find the best helix"<<endl;
  for(auto &circle : circles){
    vector<Point> points = circle->GetPoints();
    
    for(double pz = config->GetMaxPz(); pz >= config->GetMinPz() ; pz-=config->GetStepPz() ){
      
      double c = Point(circle->GetMomentum()->GetX(), circle->GetMomentum()->GetY(), pz).GetVectorSlopeC();
      unique_ptr<Helix> helix = make_unique<Helix>(c, circle, config);
      helix->SetPoints(points);
      
      int nRegularPoints = helix->GetNregularPoints();
      double fractionRegularPoints = nRegularPoints/(double)helix->GetNpoints();
      
      // Here is a condition to accept new solution as the best one
      // Accept as a new best solution if:
      // - it gives more reqular points than before or,
      // - it gives the same number of regular points, but they counstitute higher fraction of all points than before
      if(nRegularPoints < maxNregularPoints) continue;
      else if(nRegularPoints == maxNregularPoints){
        if(fractionRegularPoints - maxFractionRegularPoints < 0.001) continue;
      }
      
      // If we reach till this point, save this solution as the best one so far
      helix->SetPz(pz);
      bestHelix = move(helix);
      maxNregularPoints = nRegularPoints;
      maxFractionRegularPoints = fractionRegularPoints;
    }
  }
  
  return bestHelix;
}


int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  config = make_shared<FitterConfig>("configs/helixFitter.md");
  display = new Display();
  
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
  vector<Point> allSimplePoints; // all hits in the event
  allSimplePoints = LoadAllHits(searchRun, searchLumi, searchEvent);
  
  // what we can calculate from the assumptions
  double theta  = 2*atan(exp(-1.08));//2*atan(exp(-event->GetTrack(0)->GetEta()));
  double phi    = -2.16;//event->GetTrack(0)->GetPhi();
  double decayX = decayR*sin(theta)*cos(phi);
  double decayY = decayR*sin(theta)*sin(phi);
  double decayZ = decayR*cos(theta);
  
  trackPhi = phi;
  trackTheta = theta;
  
  // Draw decay point to make sure that it's correctly located
  vector<Point> decayPoint = {Point(decayX,decayY,decayZ)};
  display->DrawSimplePoints(decayPoint, decayPointOptions);
  
  // Draw true pion helix
  unique_ptr<Point> pionHelixCenter = unique_ptr<Point>(new Point(decayX,decayY,decayZ));
  unique_ptr<Helix> pionHelix = make_unique<Helix>(pionHelixCenter, pionVector, pionCharge, config);
  display->DrawHelix(pionHelix,helixOptions);
  
  // Calculate and draw points along the helix that hit the silicon
  vector<Point> pionPoints = pionHelix->GetPointsHittingSilicon();
  for(auto &p : pionPoints){p.SetIsPionHit(true);}
  pionHelix->SetPoints(pionPoints);
  display->DrawSimplePoints(pionPoints, pionPointsOptions);
  
  // inject hits from pion into all points in the tracker
  allSimplePoints.insert(allSimplePoints.end(),pionPoints.begin(), pionPoints.end());
  
  // remove hits that for sure don't belong to the pion's helix
  cout<<"Fitting best helix"<<endl;
  unique_ptr<Helix> bestHelix = GetBestFittingHelix(allSimplePoints);
  
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
    display->DrawSimplePoints(*bestHelix->GetPoints(), fitPointsOptions);
    
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
  
  cout<<"\n\ndone\n\n"<<endl;
//  c1->Update();
  theApp.Run();
  return 0;
}








vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber)
{
  TFile *inFile = TFile::Open("pickhists.root");
//  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists.root");
  //  TFile *inFile = TFile::Open("/afs/cern.ch/work/j/jniedzie/private/pickhists_unfiltered.root");
  if(!inFile){
    cout<<"ERROR -- no file with all hits was found"<<endl;
    return vector<Point>();
  }
  TTree *tree = (TTree*)inFile->Get("hitsExtractor/hits");
  
  if(!tree){
    cout<<"ERROR -- no tree with all hits was found"<<endl;
    return vector<Point>();
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
  
  vector<Point> pixelPoints;
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find all hits for requested event!\n\n"<<endl;
    return pixelPoints;
  }
  
  for(int i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    // convert cm to mm
    pixelPoints.push_back(Point(10*hitX->at(i),10*hitY->at(i),10*hitZ->at(i),hitCharge->at(i)));
  }
  
  display->DrawSimplePoints(pixelPoints, pixelHitsOptions);
  
  if(showStipClusters){
    vector<Point> stripPoints;
    for(int i=0;i<stripX->size();i++){
      if(stripCharge->at(i) < chargeThreshold) continue;
      stripPoints.push_back(Point(10*stripX->at(i),10*stripY->at(i),10*stripZ->at(i),stripCharge->at(i)));
    }
    display->DrawSimplePoints(stripPoints, stripHitsOptions);
  }
  
  return pixelPoints;
}
