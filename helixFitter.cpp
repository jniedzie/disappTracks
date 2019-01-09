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
#include "FitterConfig.hpp"

double pionNturns = 5;    // How many turns the pion does (should be removed and it should just go until the end of the pixel barrel

unique_ptr<FitterConfig> config;

// Will be calculated automatically
double trackEta, trackTheta, trackPhi; // parameters of the chargino track
double decayR;  // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double pionCharge = 1;

unique_ptr<Point> pionVector = nullptr;

unique_ptr<Helix> GetBestFittingHelix(const unique_ptr<Helix> &pionHelix);
vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber);

vector<Point> originalPixelPoints;

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
  trackEta = RandDouble(-config->GetMaxTrackEta(), config->GetMaxTrackEta());
  trackTheta = 2*atan(exp(-trackEta));
  trackPhi = RandDouble(0, 2*TMath::Pi());
}

/// Checks if input and output helices are identical.
/// \return Returns zero if identical, otherwise returns failure reason code
vector<int> AreHelicesIdentical(const unique_ptr<Helix> &h1, const unique_ptr<Helix> &h2)
{
  vector<int> reasons;
  
  if(fabs(h1->GetOrigin()->GetX() - h2->GetOrigin()->GetX()) > config->GetToleranceX()) reasons.push_back(1);
  if(fabs(h1->GetOrigin()->GetY() - h2->GetOrigin()->GetY()) > config->GetToleranceY()) reasons.push_back(2);
  if(fabs(h1->GetOrigin()->GetZ() - h2->GetOrigin()->GetZ()) > config->GetToleranceZ()) reasons.push_back(3);
  if(fabs(h1->GetMomentum()->GetX() - h2->GetMomentum()->GetX()) > config->GetTolerancePx()) reasons.push_back(4);
  if(fabs(h1->GetMomentum()->GetY() - h2->GetMomentum()->GetY()) > config->GetTolerancePy()) reasons.push_back(5);
  if(fabs(h1->GetMomentum()->GetZ() - h2->GetMomentum()->GetZ()) > config->GetTolerancePz()) reasons.push_back(6);
  
  return reasons;
}

Point AveragePoints(vector<Point> points)
{
  double x=0, y=0, z=0;
  for(Point p : points){
    x += p.GetX();
    y += p.GetY();
    z += p.GetZ();
  }
  x /= points.size();
  y /= points.size();
  z /= points.size();
  
  return Point(x,y,z);
}

vector<vector<Point>> SplitPointsIntoLines(vector<Point> points)
{
  vector<vector<Point>> pointsByLines;
  bool addedToExisting;
  
  for(Point p : points){
    addedToExisting = false;
    
    // loop over existing lines and check if this point belongs to one of them
    for(vector<Point> &line : pointsByLines){
      // if distance to this line is small enough, just add the point to this line and go to next point
      if(AveragePoints(line).distanceXY(p) < config->GetLinesTolerance()){
        line.push_back(p);
        addedToExisting = true;
        break;
      }
    }
    if(addedToExisting) continue;
    
    // If the point was not added to any line, create a new line for it
    vector<Point> line;
    line.push_back(p);
    pointsByLines.push_back(line);
  }
  
  return pointsByLines;
}

double minPx, minPy, minPz, maxPx, maxPy, maxPz, minL, maxL;

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = make_unique<FitterConfig>("configs/helixFitter.md");
  
  // load hits from an event (could be replaced by random points)
//  originalPixelPoints = LoadAllHits(297100, 136, 245000232);
  
  map<string, TH1D*> monitors1D;
  map<string, TH2D*> monitors2D;
  
  const vector<tuple<const char*,int,double,double>> monitors1Dparams = {
    {"nPointsOnHelix",100 ,0,100},
    {"chi2ofHelix",   50  ,0,50 },
    {"nPionPoints",   100, 0,1.2},
    {"nFakeHits",     1000,0,10 },
    {"failReason",    10,  0,10 },
  };
  
  minPx = config->GetMinPx();
  minPy = config->GetMinPy();
  minPz = config->GetMinPz();
  maxPx = config->GetMaxPx();
  maxPy = config->GetMaxPy();
  maxPz = config->GetMaxPz();
  
  const vector<tuple<const char*,int,double,double,int,double,double>> monitors2Dparams = {
    {"xResponse",     500,-250,250, 500,-250,250 },
    {"yResponse",     500,-250,250, 500,-250,250 },
    {"zResponse",     500,-250,250, 500,-250,250 },
    {"pxResponse",    200,-maxPx,maxPx, 200,-maxPx,maxPx },
    {"pyResponse",    200,-maxPy,maxPy, 200,-maxPy,maxPy },
    {"pzResponse",    200,-maxPz,maxPz, 200,-maxPz,maxPz },
  };
  
  for(auto params : monitors1Dparams){
    monitors1D[get<0>(params)] = new TH1D(get<0>(params),get<0>(params),get<1>(params),get<2>(params),get<3>(params));
  }
  
  for(auto params : monitors2Dparams){
    monitors2D[get<0>(params)] = new TH2D(get<0>(params),get<0>(params),
                                          get<1>(params),get<2>(params),get<3>(params),
                                          get<4>(params),get<5>(params),get<6>(params));
  }
  
  int nSuccess = 0;
  int nTests = config->GetNtests();
  for(int i=0;i<nTests;i++){
    SetRandomTrack();         // Randomly generate chargino's track
    FillRandomPoints(500);    // Randomly fill in pixel barrel with noise hits
    
    pionVector = make_unique<Point>(/*RandSign()*/RandDouble(minPx, maxPx),
                                    /*RandSign()*/RandDouble(minPy, maxPy),
                                    /*RandSign()*/RandDouble(minPz, maxPz));
    
    cout<<"True pion momentum vector:"; pionVector->Print();
    
    minL = config->GetMinL();
    maxL = config->GetMaxL();
    decayR = RandDouble(minL, maxL);
    cout<<"R iter:"<<i<<"\tdecayR:"<<decayR<<endl;
    
    // Create true pion helix
    double decayX = decayR*sin(trackTheta)*cos(trackPhi);
    double decayY = decayR*sin(trackTheta)*sin(trackPhi);
    double decayZ = decayR*cos(trackTheta);
    
    unique_ptr<Point> pionHelixCenter = make_unique<Point>(decayX,decayY,decayZ);
    unique_ptr<Helix> pionHelix = make_unique<Helix>(pionHelixCenter, pionVector, pionCharge,
                                                     pionNturns, config->GetHelixThickness(), config->GetZregularityTolerance());
    
    unique_ptr<Helix> bestHelix = GetBestFittingHelix(pionHelix);
    
    if(!bestHelix){
      monitors1D["failReason"]->Fill(7);
      continue;
    }
    
    monitors2D["xResponse"]->Fill(pionHelix->GetOrigin()->GetX(), bestHelix->GetOrigin()->GetX());
    monitors2D["yResponse"]->Fill(pionHelix->GetOrigin()->GetY(), bestHelix->GetOrigin()->GetY());
    monitors2D["zResponse"]->Fill(pionHelix->GetOrigin()->GetZ(), bestHelix->GetOrigin()->GetZ());
    monitors2D["pxResponse"]->Fill(pionVector->GetX(), bestHelix->GetMomentum()->GetX());
    monitors2D["pyResponse"]->Fill(pionVector->GetY(), bestHelix->GetMomentum()->GetY());
    monitors2D["pzResponse"]->Fill(pionVector->GetZ(), bestHelix->GetMomentum()->GetZ());
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
  TFile *outFile = new TFile(config->GetOutputPath(),"recreate");
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
  for(auto &p : pionPoints){p.SetIsPionHit(true);}
  pionHelix->SetPoints(pionPoints);
  vector<Point> allSimplePoints = originalPixelPoints;
  if(config->GetInjectPionHits()) allSimplePoints.insert(allSimplePoints.end(),pionPoints.begin(), pionPoints.end());
  
  // Prepare 2D projections in XY
  vector<Point> points2D;
  vector<vector<Point>> pointsByLine = SplitPointsIntoLines(allSimplePoints);
  
  for(vector<Point> line : pointsByLine){
    if(line.size() >= config->GetMinPointsAlongZ()){
      points2D.push_back(AveragePoints(line));
    }
  }
  
  // Create fitter to fit circles to 2D distribution
  unique_ptr<Fitter> fitter = unique_ptr<Fitter>(new Fitter(3));
  double helixThickness = config->GetHelixThickness();
  fitter->SetParameter(0, "L", (maxL+minL)/2., minL-helixThickness, maxL+helixThickness);
  fitter->SetParameter(1, "px", (maxPx-minPx)/2., minPx, maxPx);
  fitter->SetParameter(2, "py", (maxPy-minPy)/2., minPy, maxPy);
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  
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
  
  unique_ptr<Helix> bestHelix = nullptr;
  int maxNregularPoints = 0;
  double maxFractionRegularPoints = 0;
  
  for(auto &circle : circles){
    vector<Point> points = circle->GetPoints();
    
    for(double pz = maxPz; pz >= minPz ; pz-=config->GetStepPz() ){
      
      double c = Point(circle->GetMomentum()->GetX(), circle->GetMomentum()->GetY(), pz).GetVectorSlopeC();
      unique_ptr<Helix> helix = make_unique<Helix>(c, circle, pionNturns, helixThickness, config->GetZregularityTolerance());
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
