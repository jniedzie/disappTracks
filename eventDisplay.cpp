#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"

Display *display;

bool showStipClusters = false;

 // determines how far points can be from helix to be assigned to it (in mm)
double helixThickness = 3.0;

// constants
double B = 3.7; // T

// assumptions about the pion
double decayR = 140; // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double pionMomentum = 200; // MeV
double pionCharge = 1;
double pionR = pionMomentum/B*0.3*10; // radius of the pion spiral in mm
double pionC = 4;  // helix slope in Z direction (should be properly calculated from momentum vector
double pionNturns = 5;

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

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  
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
  
  vector<Point> allSimplePoints; // all hits in the event
  allSimplePoints = LoadAllHits(searchRun, searchLumi, searchEvent);
  
  cout<<"Pion R:"<<pionR<<endl;
  
  // what we can calculate from the assumptions
  double theta  = 2*atan(exp(-event->GetTrack(0)->GetEta()));
  double phi    = event->GetTrack(0)->GetPhi();
  double decayX = decayR*sin(theta)*cos(phi);
  double decayY = decayR*sin(theta)*sin(phi);
  double decayZ = decayR*cos(theta);
  double tShift = acos(1/sqrt(pow(decayX/decayY,2)+1));
  
  // Draw decay point to make sure that it's correctly located
  vector<Point> decayPoint = {Point(decayX,decayY,decayZ)};
  display->DrawSimplePoints(decayPoint, decayPointOptions);
  
  // Draw true pion helix
  Point pionHelixCenter(decayX,decayY,decayZ);
  pionHelixCenter.PerpendicularShift(pionR,pionC,pionCharge); // shift helix to start in the decay point
  Helix pionHelix(pionR,pionC,pionHelixCenter,helixThickness);
  display->DrawHelix(pionHelix,helixOptions,-tShift,pionNturns*2*TMath::Pi());
  
  // Calculate and draw points along the helix that hit the silicon
  vector<Point> pionPoints = pionHelix.GetPointsHittingSilicon(-tShift,pionNturns*2*TMath::Pi());
  for(auto &p : pionPoints){p.isPionHit = true;}
  pionHelix.nPionPoints = pionHelix.nPoints = (int)pionPoints.size();
  display->DrawSimplePoints(pionPoints, pionPointsOptions);
  
  // inject hits from pion into all points in the tracker
  allSimplePoints.insert(allSimplePoints.end(),pionPoints.begin(), pionPoints.end());
  
  // remove hits that for sure don't belong to the pion's helix
  
  cout<<"size before:"<<allSimplePoints.size()<<endl;
  for(int i=0;i<allSimplePoints.size();i++){
    Point p = allSimplePoints[i];
    
    if((decayZ > 0 && p.z < (decayZ-helixThickness)) || // remove wrong Z points
       (decayZ <=0 && p.z > (decayZ+helixThickness))
//       || (sqrt(pow(p.x - decayX,2) + pow(p.y - decayY,2)) > 2*pionR)
       ){
      allSimplePoints.erase(allSimplePoints.begin()+i);
      i--;
    }
  }
  cout<<"size after:"<<allSimplePoints.size()<<endl;
  
//  display->DrawSimplePoints(allSimplePoints, filteredPointsOptions);
  
  
  // Set some limits on parameters (all distances in [mm])
  
  // Pion helix parameters:
  double minR   = 25; // can't be smaller than the minimum to reach 2 different layers. Physically, from pion momentum it should be around 130 mm
  double maxR   = 200;
  double minC   = 3; // Fitter could try to make it very small, to gahter some additional points by chance while still fitting perfectly the real pion hits. This should be given more attention later...
  double maxC   = 70;
  
  // Position of the decay vertex along chargino's track (later should take into account that it's not a straight line)
  double minL   = layerR[2]; // minimum on the surface of 3rd layer
  double maxL   = layerR[3]; // maximum on the surface of 4th layer
  
  int minNpointsAlongZ = 2; // minimum number of hits for each line parallel to Z axis
  
  // Prepare 2D projections in XY
  TH2D *pointsXY = new TH2D("pointsXY","pointsXY",500/helixThickness,-250,250,500/helixThickness,-250,250);
  pointsXY->GetXaxis()->SetTitle("X");
  pointsXY->GetYaxis()->SetTitle("Y");
  for(auto point : allSimplePoints){pointsXY->Fill(point.x,point.y);}
  
  vector<pair<double, double>> points2D;
  
  for(int binX=0;binX<pointsXY->GetNbinsX();binX++){
    for(int binY=0;binY<pointsXY->GetNbinsY();binY++){
      if(pointsXY->GetBinContent(binX,binY) < minNpointsAlongZ){
        pointsXY->SetBinContent(binX,binY, 0);
      }
      else{
        points2D.push_back(make_pair(pointsXY->GetXaxis()->GetBinCenter(binX),
                                     pointsXY->GetYaxis()->GetBinCenter(binY)));
      }
    }
  }
  
  // Create fitter to fit circles to 2D distribution
  Fitter *fitter = new Fitter(2);
  fitter->SetParameter(0, "L", (maxL+minL)/2., minL-helixThickness, maxL+helixThickness);
  fitter->SetParameter(1, "R", (maxR+minR)/2., minR, maxR);
  
  // Store fitted circles for each triplet of points
  vector<Circle> circles;
  
  int nPoints = (int)points2D.size();
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){

        auto chi2Function = [&](const double *par) {
          double f = 0;
          
          double L = par[0];
          double R  = par[1];
          
          double x0 = L*sin(theta)*cos(phi);
          double y0 = L*sin(theta)*sin(phi);
          
          Point c(x0,y0,0.0);
          c.PerpendicularShift(R, 0.0);
          
          double xa = points2D[i].first;
          double ya = points2D[i].second;
          double xb = points2D[j].first;
          double yb = points2D[j].second;
          double xc = points2D[k].first;
          double yc = points2D[k].second;
          
          f  = pow(sqrt(pow(xa-c.x,2)+pow(ya-c.y,2))-R,2);
          f += pow(sqrt(pow(xb-c.x,2)+pow(yb-c.y,2))-R,2);
          f += pow(sqrt(pow(xc-c.x,2)+pow(yc-c.y,2))-R,2);
          
          return f;
        };
        fitter->SetFitFunction(chi2Function);
   
        if(fitter->RunFitting()) {
          auto result = fitter->GetResult();
          
          double x0 = result.GetParams()[0]*sin(theta)*cos(phi);
          double y0 = result.GetParams()[0]*sin(theta)*sin(phi);
          double z0 = result.GetParams()[0]*cos(theta);
          
          double R  = result.GetParams()[1];
          
          Point cc(x0,y0,z0);
          double tShift = cc.PerpendicularShift(R, 0.0);
          Circle circle(cc.x,cc.y,R);
          circle.z = cc.z;
          circle.tShift = tShift;
          
          int nCircleBins = circle.GetNbinsOverlappingWithHist(pointsXY);
          
          if(nCircleBins > 3){
            circles.push_back(circle);
//            result.Print(std::cout);
          }
        }
      }
    }
  }
  
  // Plot the results
  TCanvas *c1 = new TCanvas("c1","c1",1280,1800);
  c1->Divide(2,2);
  c1->cd(1);
  pointsXY->DrawCopy("colz");
  
  for(auto circle : circles){
    TArc *circleArc = circle.GetArc();
    circleArc->SetLineColor(kRed);
    circleArc->SetLineWidth(2);
    circleArc->ResetAttFill();
    circleArc->Draw("sameL");
  }
  Circle pionCircle = Circle(pionHelixCenter.x, pionHelixCenter.y, pionR);
  TArc *pionCircleArc = pionCircle.GetArc();
  
  pionCircleArc->SetLineColor(kGreen);
  pionCircleArc->SetLineWidth(2);
  pionCircleArc->ResetAttFill();
  pionCircleArc->Draw("sameL");
  
  cout<<"N circles found:"<<circles.size()<<endl;
  
  //
  // verified up to this point
  //
  
  cout<<"Assigning points to circles"<<endl;
  for(auto &circle : circles){
    for(auto point : allSimplePoints){
      Point q = circle.GetClosestPoint(point);
      q.z = point.z;
      double d = q.distance(point);
      if(d < helixThickness){
        circle.points.push_back(point);
      }
    }
  }
  
  cout<<"n circles before:"<<circles.size()<<endl;
  for(int i=0;i<circles.size();i++){
    if(circles[i].points.size() < 8){
      circles.erase(circles.begin()+i);
      i--;
    }
  }
  cout<<"n circles after:"<<circles.size()<<endl;
 
  // Create fitter to fit helices to 3D points assigned to circles
  /*
  double minZ = minL*cos(theta);
  double maxZ = maxL*cos(theta);
  
  cout<<"Fitter 3D min z:"<<minZ-helixThickness<<"\tmax z:"<<maxZ+helixThickness<<endl;
  cout<<"\tmin c:"<<minC<<"\tmax c:"<<maxC<<endl;
  
  Fitter *fitter3D = new Fitter(1);
//  fitter3D->SetParameter(0, "z0", (maxZ+minZ)/2., minZ-helixThickness, maxZ+helixThickness);
  fitter3D->SetParameter(0, "c" , (maxC+minC)/2., minC, maxC);
  
  for(auto &circle : circles){
    vector<Point> points = circle.points;
    
    double x0 = circle.x;
    double y0 = circle.y;
    double z0 = circle.z;
    double R = circle.R;
    double tShift = acos(1/sqrt(pow(x0/y0,2)+1));
    
    auto chi2Function = [&](const double *par){
      double f = 0;
      
      double c  = par[0];
      
      
      Helix helix(R, c, x0, y0, z0+tShift*c, helixThickness);
      
      for(auto point : points){
        Point q = helix.GetClosestPoint(point);
        double d = q.distance(point);
        
        f+= d*d;
      }
      f /= points.size();
      return f;
    };
    fitter3D->SetFitFunction(chi2Function);
    bool fitterStatus = fitter3D->RunFitting();
    
    if(1){
      auto result = fitter3D->GetResult();
      result.Print(std::cout);
      
      double c = result.GetParams()[0];
      
      Helix helix(R,c,x0,y0,z0+tShift*c,helixThickness);
      helix.CountMatchingPoints(points);
      
      circle.helix = helix;
      circle.chi2 = result.MinFcnValue();
    }
  }
  */
  
  Helix bestHelix;
  Circle bestCircle;
  int maxNpoints = 0;
  double bestChi2 = 9999999;
  vector<Point> bestPointsOnHelix;
  
  for(auto &circle : circles){
    vector<Point> points = circle.points;
    
    double x0 = circle.x;
    double y0 = circle.y;
    double R = circle.R;
    
    int xSign=1, ySign=1;
    if(x0> 0 && y0> 0){xSign= 1; ySign=-1;}
    if(x0<=0 && y0> 0){xSign= 1; ySign= 1;}
    if(x0<=0 && y0<=0){xSign=-1; ySign= 1;}
    if(x0> 0 && y0<=0){xSign=-1; ySign=-1;}
    double x = x0 - xSign * R/sqrt(pow(x0/y0,2)+1);
    double y = y0 - ySign * R/sqrt(pow(y0/x0,2)+1);
    double tShift = acos(1/sqrt(pow(x/y,2)+1));
  
    for(double c=maxC;c>=minC;c-=0.1){
      double z0 = circle.z + c*tShift;
      Helix helix(R, c, x0, y0, z0, helixThickness);
      helix.CountMatchingPoints(points);
      
      double chi2=0;
      vector<Point> pointsOnHelix;
      for(auto point : points){
        Point q = helix.GetClosestPoint(point);
        pointsOnHelix.push_back(q);
        double d = q.distance(point);
        chi2 += d*d;
      }
      chi2 /= points.size();
      
      if(helix.nPoints > maxNpoints){
        bestCircle = circle;
        maxNpoints = helix.nPoints;
        bestChi2 = chi2;
        bestHelix = helix;
        bestPointsOnHelix = pointsOnHelix;
      }
      else if(helix.nPoints == maxNpoints){
        if(chi2 < bestChi2){
          bestCircle = circle;
          bestChi2 = chi2;
          bestHelix = helix;
          bestHelix.tShift = circle.tShift;
          bestPointsOnHelix = pointsOnHelix;
        }
      }
    }
  }
  
  // Try to increase c and check if we don't loose points
  tShift = acos(1/sqrt(pow(bestHelix.x0/bestHelix.y0,2)+1));
  
  map<string,any> bestHelixOptions = {
    {"title", "Best helix"},
    {"markerStyle", 20},
    {"markerSize", 0.2},
    {"color", kRed}
  };

  display->DrawHelix(bestHelix, bestHelixOptions,-bestHelix.tShift,pionNturns*2*TMath::Pi());
  
  map<string,any> fitPointsOptions = {
    {"title", "Fit helix points"},
    {"markerStyle", 20},
    {"markerSize", 1.0},
    {"color", kRed}
  };
  display->DrawSimplePoints(bestHelix.points, fitPointsOptions);
  
  map<string,any> closestPointsOptions = {
    {"title", "Closest points"},
    {"markerStyle", 20},
    {"markerSize", 1.0},
    {"color", kOrange}
  };
  display->DrawSimplePoints(bestPointsOnHelix, closestPointsOptions);
  
  cout<<"\n\nBest helix is:"<<endl;
  bestHelix.Print();
  
  cout<<"\nPion helix:"<<endl;
  pionHelix.Print();
  
  
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
  c1->Update();
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
