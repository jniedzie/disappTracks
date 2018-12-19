//
//  helixFitter.cpp
//
//  Created by Jeremi Niedziela on 17/12/2018.
//

#include "Helpers.hpp"
#include "Helix.hpp"
#include "Circle.hpp"
#include "Point.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"

// determines how far points can be from helix to be assigned to it (in mm)
double helixThickness = 3.0;

bool injectPionHits = true;
int nTests = 100;
const char* outFileName = "tests.root";

// constants
double B = 3.7; // T

// assumptions about the pion
double decayR = 140; // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double pionMomentum = 200; // MeV
double pionCharge = 1;
double pionR = pionMomentum/B*0.3*10; // radius of the pion spiral in mm
double pionC = 4;  // helix slope in Z direction (should be properly calculated from momentum vector
double pionNturns = 5;

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


// what we can calculate from the assumptions
double eta = 1.8;
double theta  = 2*atan(exp(-eta));
double phi    = 0.45*TMath::Pi();

unique_ptr<Helix> GetBestFittingHelix(const unique_ptr<Helix> &pionHelix);
vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber);

vector<Point> originalPixelPoints;

inline double RandDouble(double min, double max)
{
  return min + static_cast<double>(rand()) /( static_cast<double>(RAND_MAX/(max-min)));
}

void FillRandomPoints(int nPoints)
{
  originalPixelPoints.clear();
  double phi, R;
  int layerIndex;
  
  for(int i=0;i<nPoints;i++){
    phi = RandDouble(0, 2*TMath::Pi());
    layerIndex = RandDouble(0, 4);
    R = layerR[layerIndex];
    Point p(R*cos(phi), R*sin(phi), RandDouble(-500, 500));
    originalPixelPoints.push_back(p);
  }
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  // load hits from an event (could be replaced by random points)
//  originalPixelPoints = LoadAllHits(297100, 136, 245000232);
//  Display *display = new Display();
//
//  map<string,any> options = {
//    {"title", "points"},
//    {"color", kRed},
//    {"markerStyle", 20},
//    {"markerSize", 1.0}
//  };
//
//  display->DrawSimplePoints(originalPixelPoints, options);
  
  TH2D *radiusResponse = new TH2D("radiusResponse","radiusResponse",100,minR,maxR,100,minR,maxR);
  radiusResponse->GetXaxis()->SetTitle("R_{true}");
  radiusResponse->GetYaxis()->SetTitle("R_{fit}");
  
  TH2D *cResponse = new TH2D("cResponse","cResponse",100,minC,maxC,100,minC,maxC);
  cResponse->GetXaxis()->SetTitle("c_{true}");
  cResponse->GetYaxis()->SetTitle("c_{fit}");
  
  TH2D *xResponse = new TH2D("xResponse","xResponse",500,-250,250,500,-250,250);
  xResponse->GetXaxis()->SetTitle("x_{true}");
  xResponse->GetYaxis()->SetTitle("x_{fit}");
  
  TH2D *yResponse = new TH2D("yResponse","yResponse",500,-250,250,500,-250,250);
  yResponse->GetXaxis()->SetTitle("y_{true}");
  yResponse->GetYaxis()->SetTitle("y_{fit}");
  
  TH2D *zResponse = new TH2D("zResponse","zResponse",500,-500,500,500,-500,500);
  zResponse->GetXaxis()->SetTitle("z_{true}");
  zResponse->GetYaxis()->SetTitle("z_{fit}");
  
  TH1D *nPointsOnHelix = new TH1D("nPointsOnHelix","nPointsOnHelix",100,0,100);
  nPointsOnHelix->GetXaxis()->SetTitle("N point on helix");
  
  TH1D *chi2ofHelix = new TH1D("chi2ofHelix","chi2ofHelix",500,0,500);
  chi2ofHelix->GetXaxis()->SetTitle("Chi2");
  
  TH1D *nPionPoints = new TH1D("nPionPoints","nPionPoints",1000,0,10);
  nPionPoints->GetXaxis()->SetTitle("N_{fit}/N_{true}");
  
  for(int i=0;i<nTests;i++){
    FillRandomPoints(500);
    
    pionR = RandDouble(minR, maxR);
    pionC = RandDouble(minC, maxC);
    decayR = RandDouble(minL, maxL);
    cout<<"R iter:"<<i<<"\tR:"<<pionR<<"\tc:"<<pionC<<"\tdecayR:"<<decayR<<endl;
    
    // Create true pion helix
    double decayX = decayR*sin(theta)*cos(phi);
    double decayY = decayR*sin(theta)*sin(phi);
    double decayZ = decayR*cos(theta);
    
    Point pionHelixCenter(decayX,decayY,decayZ);
    unique_ptr<Helix> pionHelix = unique_ptr<Helix>(new Helix(pionR,pionC,pionHelixCenter,pionNturns, helixThickness));
    pionHelix->Shift(pionCharge); // shift helix to start in the decay point
    
    unique_ptr<Helix> bestHelix = GetBestFittingHelix(pionHelix);
    
    radiusResponse->Fill(pionHelix->R, bestHelix->R);
    cResponse->Fill(pionHelix->c, bestHelix->c);
    xResponse->Fill(pionHelix->x0, bestHelix->x0);
    yResponse->Fill(pionHelix->y0, bestHelix->y0);
    zResponse->Fill(pionHelix->z0, bestHelix->z0);
    nPointsOnHelix->Fill(bestHelix->nPoints);
    chi2ofHelix->Fill(bestHelix->chi2 < 500 ? bestHelix->chi2 : 499);
    nPionPoints->Fill(bestHelix->nPionPoints/(double)pionHelix->nPionPoints);
    
    cout<<"Pion helix:"<<endl;
    pionHelix->Print();
    cout<<"Fitted helix:"<<endl;
    bestHelix->Print();
  }
  
  // Plot the results
  TCanvas *c1 = new TCanvas("c1","c1",1280,1000);
  c1->Divide(3,3);
  c1->cd(1);
  radiusResponse->Draw("colz");
  c1->cd(2);
  cResponse->Draw("colz");
  c1->cd(3);
  xResponse->Draw("colz");
  c1->cd(4);
  yResponse->Draw("colz");
  c1->cd(5);
  zResponse->Draw("colz");
  c1->cd(6);
  nPointsOnHelix->Draw();
  c1->cd(7);
  chi2ofHelix->Draw();
  c1->cd(8);
  nPionPoints->Draw();
  
  TFile *outFile = new TFile(outFileName,"recreate");
  outFile->cd();
  radiusResponse->Write();
  cResponse->Write();
  xResponse->Write();
  yResponse->Write();
  zResponse->Write();
  nPointsOnHelix->Write();
  chi2ofHelix->Write();
  nPionPoints->Write();
  outFile->Close();
  
  c1->Update();
  theApp.Run();
  return 0;
}

unique_ptr<Helix> GetBestFittingHelix(const unique_ptr<Helix> &pionHelix)
{
  // Calculate points along the helix that hit the silicon and inject them into all points in the tracker
  vector<Point> pionPoints = pionHelix->GetPointsHittingSilicon();
  for(auto &p : pionPoints){p.isPionHit = true;}
  pionHelix->nPionPoints = pionHelix->nPoints = (int)pionPoints.size();
  vector<Point> allSimplePoints = originalPixelPoints;
  if(injectPionHits) allSimplePoints.insert(allSimplePoints.end(),pionPoints.begin(), pionPoints.end());
  
  // remove hits that for sure don't belong to the pion's helix
  for(int i=0;i<allSimplePoints.size();i++){
    Point p = allSimplePoints[i];
    
    // --------------------------------
    // TODO: Here I use truth info to remove some htis -> should be independent of pionHelix and decay point params!!!
    // --------------------------------
    
    if((pionHelix->z0 > 0 && p.z < (pionHelix->z0-helixThickness)) || // remove wrong Z points
       (pionHelix->z0 <=0 && p.z > (pionHelix->z0+helixThickness))){
      allSimplePoints.erase(allSimplePoints.begin()+i);
      i--;
    }
  }
  static int plotIter=0;
  // Prepare 2D projections in XY
  TH2D *pointsXY = new TH2D(Form("pointsXY%i",plotIter),Form("pointsXY%i",plotIter),500/helixThickness,-250,250,500/helixThickness,-250,250);
  plotIter++;
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
  unique_ptr<Fitter> fitter = unique_ptr<Fitter>(new Fitter(2));
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
          
          Circle c(x0,y0,R);
          c.Shift();
          
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
          
          double L = result.GetParams()[0];
          double R  = result.GetParams()[1];
          
          double x0 = L*sin(theta)*cos(phi);
          double y0 = L*sin(theta)*sin(phi);
          double z0 = L*cos(theta);
          
          Circle circle(x0,y0,R);
          circle.Shift();
          circle.z = z0;
          
          if(circle.GetNbinsOverlappingWithHist(pointsXY) > 3) circles.push_back(circle);
        }
      }
    }
  }
  
  for(auto &circle : circles){
    for(auto point : allSimplePoints){
      Point q = circle.GetClosestPoint(point);
      q.z = point.z;
      if(q.distance(point) < helixThickness){
        circle.points.push_back(point);
      }
    }
  }
  
  for(int i=0;i<circles.size();i++){
    if(circles[i].points.size() < 8){
      circles.erase(circles.begin()+i);
      i--;
    }
  }
  
  unique_ptr<Helix> bestHelix = unique_ptr<Helix>(new Helix());
  Circle bestCircle(0,0,0);
  int maxNpoints = 0;
  double bestChi2 = 9999999;
  
  for(auto &circle : circles){
    vector<Point> points = circle.points;
    
    for(double c=maxC;c>=minC;c-=0.1){
      unique_ptr<Helix> helix = unique_ptr<Helix>(new Helix(c, circle, pionNturns, helixThickness));
      helix->CountMatchingPoints(points);
      
      double chi2=0;
      vector<Point> pointsOnHelix;
      for(auto point : points){
        Point q = helix->GetClosestPoint(point);
        pointsOnHelix.push_back(q);
        double d = q.distance(point);
        chi2 += d*d;
      }
      chi2 /= points.size();
      
      if(helix->nPoints > maxNpoints){
        bestCircle = circle;
        maxNpoints = helix->nPoints;
        bestChi2 = chi2;
        bestHelix = move(helix);
      }
      else if(helix->nPoints == maxNpoints){
        if(chi2 < bestChi2){
          bestCircle = circle;
          bestChi2 = chi2;
          bestHelix = move(helix);
        }
      }
    }
  }
  return bestHelix;
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
  
  // Parameters for all hits in the pixel barrel
  const double chargeThreshold = 0; // 2000, 5000, 25000
  const double minClusterSize = 0;
  const double maxClusterSize = 100;
  
  for(int i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    // convert cm to mm
    pixelPoints.push_back(Point(10*hitX->at(i),10*hitY->at(i),10*hitZ->at(i),hitCharge->at(i)));
  }
  
  return pixelPoints;
}
