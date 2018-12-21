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


// Tunable parameters of the algorithm
const double helixThickness = 1.0;  // determines how far points can be from helix to be assigned to it (in mm)
const double circleThickness = 1.0; // determines how far points can be from circle to be assigned to it (in mm)
const double stepPz = 0.5;
const double zRegularityTolerance = 1.0;
const int minNpointsAlongZ = 2; // minimum number of hits for each line parallel to Z axis

// Settings
bool injectPionHits = true;
int nTests = 5;
const char* outFileName = "tests.root";

// Parameters to tune
double pionNturns = 5;    // How many turns the pion does (should be removed and it should just go until the end of the pixel barrel

double minPx = 50;
double minPy = 50;
double minPz = 50;

double maxPx = 250;
double maxPy = 250;
double maxPz = 250;

double minL   = layerR[2]; // minimum on the surface of 3rd layer
double maxL   = layerR[3]; // maximum on the surface of 4th layer

// Conditinos to accept fitted helix as a proper solution
const double toleranceR = 10; // mm
const double toleranceC = 0.1;
const double toleranceX = 10; // mm
const double toleranceY = 10; // mm
const double toleranceZ = 10; // mm

// Will be calculated automatically or don't matter
double trackEta, trackTheta, trackPhi; // parameters of the chargino track
double decayR;  // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double pionR;   // radius of the pion spiral in mm
double pionC;   // helix slope in Z direction (should be properly calculated from momentum vector
double pionCharge = 1;

unique_ptr<Point> pionVector = nullptr;

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

void SetRandomTrack()
{
  trackEta = 1.08;//RandDouble(-1.2, 1.2);
  trackTheta = 2*atan(exp(-trackEta));
  trackPhi = -2.16;//RandDouble(0, 2*TMath::Pi());
}

/// Checks if input and output helices are identical.
/// \return Returns zero if identical, otherwise returns failure reason code
vector<int> AreHelicesIdentical(const unique_ptr<Helix> &h1, const unique_ptr<Helix> &h2)
{
  vector<int> reasons;
  
  if(fabs(h1->GetOrigin()->x - h2->GetOrigin()->x) > toleranceX) reasons.push_back(1);
  if(fabs(h1->GetOrigin()->y - h2->GetOrigin()->y) > toleranceY) reasons.push_back(2);
  if(fabs(h1->GetOrigin()->z - h2->GetOrigin()->z) > toleranceZ) reasons.push_back(3);
  if(fabs(h1->GetRadius()    - h2->GetRadius()   ) > toleranceR) reasons.push_back(4);
  if(fabs(h1->GetSlope()     - h2->GetSlope()    ) > toleranceC) reasons.push_back(5);
  
  return reasons;
}

double minR = GetRadiusInMagField(minPx, minPy, solenoidField);
double maxR = GetRadiusInMagField(maxPx, maxPy, solenoidField);

double minC = Point(maxPx, maxPy, minPz).GetVectorSlopeC();
double maxC = Point(maxPx, maxPy, maxPz).GetVectorSlopeC();

vector<tuple<const char*,int,double,double>> monitors1Dparams = {
  {"nPointsOnHelix",100 ,0,100},
  {"chi2ofHelix",   50  ,0,50 },
  {"nPionPoints",   100, 0,1.2},
  {"nFakeHits",     1000,0,10 },
  {"failReason",    10,  0,10 },
};

vector<tuple<const char*,int,double,double,int,double,double>> monitors2Dparams = {
  {"radiusResponse",100,minR,maxR,100,minR,maxR},
  {"cResponse",     100,minC,maxC,100,minC,maxC},
  {"xResponse",     500,-250,250, 500,-250,250 },
  {"yResponse",     500,-250,250, 500,-250,250 },
  {"zResponse",     500,-250,250, 500,-250,250 },
  {"pxResponse",    200,-maxPx,maxPx, 200,-maxPx,maxPx },
  {"pyResponse",    200,-maxPy,maxPy, 200,-maxPy,maxPy },
  {"pzResponse",    200,-maxPz,maxPz, 200,-maxPz,maxPz },
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  // load hits from an event (could be replaced by random points)
//  originalPixelPoints = LoadAllHits(297100, 136, 245000232);
  
  map<string, TH1D*> monitors1D;
  map<string, TH2D*> monitors2D;
  
  for(auto params : monitors1Dparams){
    monitors1D[get<0>(params)] = new TH1D(get<0>(params),get<0>(params),get<1>(params),get<2>(params),get<3>(params));
  }
  
  for(auto params : monitors2Dparams){
    monitors2D[get<0>(params)] = new TH2D(get<0>(params),get<0>(params),
                                          get<1>(params),get<2>(params),get<3>(params),
                                          get<4>(params),get<5>(params),get<6>(params));
  }
  
  int nSuccess = 0;
  
  for(int i=0;i<nTests;i++){
    SetRandomTrack();         // Randomly generate chargino's track
    FillRandomPoints(500);    // Randomly fill in pixel barrel with noise hits
    
    pionVector = make_unique<Point>(RandDouble(minPx, maxPx),
                                    RandDouble(minPy, maxPy),
                                    RandDouble(minPz, maxPz));
    
    cout<<"True pion momentum vector:"; pionVector->Print();
    
    pionR = GetRadiusInMagField(pionVector->x, pionVector->y, solenoidField);
    pionC = pionVector->GetVectorSlopeC();
    decayR = RandDouble(minL, maxL);
    cout<<"R iter:"<<i<<"\tR:"<<pionR<<"\tc:"<<pionC<<"\tdecayR:"<<decayR<<endl;
    
    // Create true pion helix
    double decayX = decayR*sin(trackTheta)*cos(trackPhi);
    double decayY = decayR*sin(trackTheta)*sin(trackPhi);
    double decayZ = decayR*cos(trackTheta);
    
    unique_ptr<Point> pionHelixCenter = make_unique<Point>(decayX,decayY,decayZ);
    unique_ptr<Helix> pionHelix = make_unique<Helix>(pionHelixCenter, pionVector, pionCharge,
                                                     pionNturns, helixThickness, zRegularityTolerance);
    
    unique_ptr<Helix> bestHelix = GetBestFittingHelix(pionHelix);
    
    monitors2D["radiusResponse"]->Fill(pionHelix->GetRadius(), bestHelix->GetRadius());
    monitors2D["cResponse"]->Fill(pionHelix->GetSlope(), bestHelix->GetSlope());
    monitors2D["xResponse"]->Fill(pionHelix->GetOrigin()->x, bestHelix->GetOrigin()->x);
    monitors2D["yResponse"]->Fill(pionHelix->GetOrigin()->y, bestHelix->GetOrigin()->y);
    monitors2D["zResponse"]->Fill(pionHelix->GetOrigin()->z, bestHelix->GetOrigin()->z);
    monitors2D["pxResponse"]->Fill(pionVector->x, bestHelix->GetMomentum()->x);
    monitors2D["pyResponse"]->Fill(pionVector->y, bestHelix->GetMomentum()->y);
    monitors2D["pzResponse"]->Fill(pionVector->z, bestHelix->GetMomentum()->z);
    monitors1D["nPointsOnHelix"]->Fill(bestHelix->GetNpoints());
    monitors1D["chi2ofHelix"]->Fill(bestHelix->GetChi2() < 50 ? bestHelix->GetChi2() : 49);
    monitors1D["nPionPoints"]->Fill(bestHelix->GetNpionPoints()/(double)pionHelix->GetNpionPoints());
    monitors1D["nFakeHits"]->Fill((bestHelix->GetNpoints()-bestHelix->GetNpionPoints())/(double)bestHelix->GetNpoints());
    
    vector<int> failureCodes = AreHelicesIdentical(pionHelix, bestHelix);
    if(failureCodes.size()==0) nSuccess++;
    else{
      for(int f : failureCodes) monitors1D["failReason"]->Fill(f);
    }
    
    cout<<"Pion helix:"<<endl;
    pionHelix->Print();
    cout<<"Fitted helix:"<<endl;
    bestHelix->Print();
  }
  
  // Plot the results
  TCanvas *c1 = new TCanvas("c1","c1",1280,1000);
  c1->Divide(4,4);
  TFile *outFile = new TFile(outFileName,"recreate");
  outFile->cd();
  
  int i=1;
  for(auto &[title, hist] : monitors2D){
    c1->cd(i++);
    hist->Draw("colz");
    hist->Write();
  }

  for(auto &[title, hist] : monitors1D){
    c1->cd(i++);
    hist->Draw();
    hist->Write();
  }
  outFile->Close();
  
  cout<<"Percentage of successful fits:"<<nSuccess/(double)nTests<<endl;
  
  c1->Update();
  theApp.Run();
  return 0;
}

unique_ptr<Helix> GetBestFittingHelix(const unique_ptr<Helix> &pionHelix)
{
  // Calculate points along the helix that hit the silicon and inject them into all points in the tracker
  vector<Point> pionPoints = pionHelix->GetPointsHittingSilicon();
  for(auto &p : pionPoints){p.isPionHit = true;}
  pionHelix->SetPoints(pionPoints);
  vector<Point> allSimplePoints = originalPixelPoints;
  if(injectPionHits) allSimplePoints.insert(allSimplePoints.end(),pionPoints.begin(), pionPoints.end());
  
  // remove hits that for sure don't belong to the pion's helix
  for(int i=0;i<allSimplePoints.size();i++){
    Point p = allSimplePoints[i];
    
    // --------------------------------
    // TODO: Here I use truth info to remove some htis -> should be independent of pionHelix and decay point params!!!
    // --------------------------------
    
    if((pionHelix->GetOrigin()->z > 0 && p.z < (pionHelix->GetOrigin()->z - helixThickness)) || // remove wrong Z points
       (pionHelix->GetOrigin()->z <=0 && p.z > (pionHelix->GetOrigin()->z + helixThickness))){
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
  
  vector<Point> points2D;
  
  for(int binX=0;binX<pointsXY->GetNbinsX();binX++){
    for(int binY=0;binY<pointsXY->GetNbinsY();binY++){
      if(pointsXY->GetBinContent(binX,binY) < minNpointsAlongZ){
        pointsXY->SetBinContent(binX,binY, 0);
      }
      else{
        points2D.push_back(Point(pointsXY->GetXaxis()->GetBinCenter(binX),
                                 pointsXY->GetYaxis()->GetBinCenter(binY),
                                 0.0));
      }
    }
  }
  
  // Create fitter to fit circles to 2D distribution
  unique_ptr<Fitter> fitter = unique_ptr<Fitter>(new Fitter(3));
  fitter->SetParameter(0, "L", (maxL+minL)/2., minL-helixThickness, maxL+helixThickness);
  fitter->SetParameter(1, "px", (maxPx-minPx)/2., minPx, maxPx);
  fitter->SetParameter(2, "py", (maxPy-minPy)/2., minPy, maxPy);
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  
  int nPoints = (int)points2D.size();
  
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
          
          if(circle->GetNbinsOverlappingWithHist(pointsXY) > 3) circles.push_back(move(circle));
        }
      }
    }
  }
  
  cout<<"N circles found:"<<circles.size()<<endl;
  for(int i=0;i<circles.size();i++){
    circles[i]->SetPoints(allSimplePoints);
    if(circles[i]->GetNpoints() < 8){
      circles.erase(circles.begin()+i);
      i--;
    }
  }
  cout<<"N circles after cleaning:"<<circles.size()<<endl;
  
  unique_ptr<Helix> bestHelix = nullptr;
  int maxNregularPoints = 0;
  double maxFractionRegularPoints = 0;
  
  for(auto &circle : circles){
    vector<Point> points = circle->GetPoints();
    
    for(double pz = maxPz; pz >= minPz ; pz-=stepPz ){
      
      double c = Point(circle->GetMomentum()->x, circle->GetMomentum()->y, pz).GetVectorSlopeC();
      unique_ptr<Helix> helix = make_unique<Helix>(c, circle, pionNturns, helixThickness, zRegularityTolerance);
      helix->SetPoints(points);
      
      int nRegularPoints = helix->GetNregularPoints();
      double fractionRegularPoints = nRegularPoints/(double)helix->GetNpoints();
      
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
