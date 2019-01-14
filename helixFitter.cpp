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

string configPath = "configs/helixFitter.md";
shared_ptr<FitterConfig> config;

// Will be calculated automatically
double trackEta, trackTheta, trackPhi; // parameters of the chargino track
double decayR;  // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double pionCharge = 1;

double minPx, minPy, minPz, maxPx, maxPy, maxPz, minL, maxL;

// Monitoring histograms
map<string, TH1D*> monitors1D;
map<string, TH2D*> monitors2D;

unique_ptr<Helix> GetBestFittingHelix(vector<Point> allSimplePoints);
vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber);

void SetRandomTrack()
{
  trackEta = RandDouble(-config->GetMaxTrackEta(), config->GetMaxTrackEta());
  trackTheta = 2*atan(exp(-trackEta));
  trackPhi = RandDouble(0, 2*TMath::Pi());
  minL = config->GetMinL();
  maxL = config->GetMaxL();
  decayR = RandDouble(minL, maxL);
}

void SetupMonitors()
{
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
}

bool FillMonitors(const unique_ptr<Helix> &fittedHelix, const unique_ptr<Helix> &trueHelix)
{
  bool success = false;
  if(!fittedHelix){
    monitors1D["failReason"]->Fill(7);
    return success;
  }
  
  monitors2D["xResponse"]->Fill(trueHelix->GetOrigin()->GetX(), fittedHelix->GetOrigin()->GetX());
  monitors2D["yResponse"]->Fill(trueHelix->GetOrigin()->GetY(), fittedHelix->GetOrigin()->GetY());
  monitors2D["zResponse"]->Fill(trueHelix->GetOrigin()->GetZ(), fittedHelix->GetOrigin()->GetZ());
  monitors2D["pxResponse"]->Fill(trueHelix->GetMomentum()->GetX(), fittedHelix->GetMomentum()->GetX());
  monitors2D["pyResponse"]->Fill(trueHelix->GetMomentum()->GetY(), fittedHelix->GetMomentum()->GetY());
  monitors2D["pzResponse"]->Fill(trueHelix->GetMomentum()->GetZ(), fittedHelix->GetMomentum()->GetZ());
  monitors1D["nPointsOnHelix"]->Fill(fittedHelix->GetNpoints());
  monitors1D["chi2ofHelix"]->Fill(fittedHelix->GetChi2() < 50 ? fittedHelix->GetChi2() : 49);
  monitors1D["nPionPoints"]->Fill(fittedHelix->GetNpionPoints()/(double)trueHelix->GetNpionPoints());
  monitors1D["nFakeHits"]->Fill((fittedHelix->GetNpoints()-fittedHelix->GetNpionPoints())/(double)fittedHelix->GetNpoints());
  
  vector<int> failureCodes = Helix::AreHelicesIdentical(fittedHelix, trueHelix);
  
  if(failureCodes.size()==0) success = true;
  else  for(int f : failureCodes) monitors1D["failReason"]->Fill(f);
  
  return success;
}

unique_ptr<Helix> GetPionHelix(const unique_ptr<Point> &pionVector)
{
  // Create true pion helix
  double decayX = decayR*sin(trackTheta)*cos(trackPhi);
  double decayY = decayR*sin(trackTheta)*sin(trackPhi);
  double decayZ = decayR*cos(trackTheta);
  
  unique_ptr<Point> pionHelixCenter = make_unique<Point>(decayX,decayY,decayZ);
  unique_ptr<Helix> pionHelix = make_unique<Helix>(pionHelixCenter, pionVector, pionCharge, config);
  vector<Point> pionPoints = pionHelix->GetPointsHittingSilicon();
  for(auto &p : pionPoints){p.SetIsPionHit(true);}
  pionHelix->SetPoints(pionPoints);
  
  return pionHelix;
}

void InjectPionPointsToCollectionOfPoints(const unique_ptr<Helix> &pionHelix, vector<Point> &pixelPoints)
{
  vector<Point> *pionPoints = pionHelix->GetPoints();
   pixelPoints.insert(pixelPoints.end(),pionPoints->begin(), pionPoints->end());
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = make_shared<FitterConfig>(configPath);
  SetupMonitors();
  
  // load hits from an event (could be replaced by random points)
//  originalPixelPoints = LoadAllHits(297100, 136, 245000232);
  
  int nSuccess = 0;
  int nTests = config->GetNtests();
  
  for(int i=0;i<nTests;i++){
    cout<<"\n========================================================"<<endl;
    cout<<"Test iter:"<<i<<endl;
    
    SetRandomTrack(); // Randomly generate chargino's track
    vector<Point> pixelPoints = Point::GetRandomPoints(config->GetNnoiseHits());
    
    unique_ptr<Point> pionVector = make_unique<Point>(/*RandSign()*/RandDouble(minPx, maxPx),
                                                      /*RandSign()*/RandDouble(minPy, maxPy),
                                                      /*RandSign()*/RandDouble(minPz, maxPz));
    
    unique_ptr<Helix> pionHelix = GetPionHelix(pionVector);
    if(config->GetInjectPionHits()){
      InjectPionPointsToCollectionOfPoints(pionHelix, pixelPoints);
    }
    
    unique_ptr<Helix> bestHelix = GetBestFittingHelix(pixelPoints);
    
    bool success = FillMonitors(bestHelix, pionHelix);
    
    if(success){
      nSuccess++;
      cout<<"Pion helix:"; pionHelix->Print();
      cout<<"Fitted helix:"; bestHelix->Print();
    }
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

vector<unique_ptr<Circle>> FitCirclesToPoints(vector<Point> allSimplePoints)
{
  // Prepare 2D projections in XY
  vector<Point> points2D;
  vector<vector<Point>> pointsByLine = Point::SplitPointsIntoLines(allSimplePoints, config->GetLinesToleranceForCircles());
  
  for(vector<Point> line : pointsByLine){
    if(line.size() >= config->GetMinPointsAlongZ()){
      points2D.push_back(Point(line));
    }
  }
  cout<<"N 2D points:"<<points2D.size()<<endl;
  
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
  if(circles.size() == 0) return circles;
  
  cout<<"N circles:"<<circles.size()<<endl;
  
  sort(circles.begin(), circles.end(),
       [](const auto &c1, const auto &c2) -> bool {return c1->GetRadius() < c2->GetRadius();});
  
  for(int i=0; i<circles.size()-1; i++){
    if(fabs(circles[i]->GetRadius() - circles[i+1]->GetRadius()) < config->GetCircleThickness()){
      circles.erase(circles.begin()+i);
      i--;
    }
  }
  cout<<"N circles after:"<<circles.size()<<endl;
  
  return circles;
}

unique_ptr<Helix> GetBestFittingHelix(vector<Point> allSimplePoints)
{
  vector<unique_ptr<Circle>> circles = FitCirclesToPoints(allSimplePoints);
  if(circles.size() == 0){
    cout<<"No circles were found"<<endl;
    return nullptr;
  }
  
  unique_ptr<Helix> bestHelix = nullptr;
  int maxNregularPoints = 0;
  double maxFractionRegularPoints = 0;
  
  for(auto &circle : circles){
    vector<Point> points = circle->GetPoints();
    
    for(double pz = maxPz; pz >= minPz ; pz-=config->GetStepPz()){
      double c = Point(circle->GetMomentum()->GetX(), circle->GetMomentum()->GetY(), pz).GetVectorSlopeC();
      unique_ptr<Helix> helix = make_unique<Helix>(c, circle, config);
      helix->SetPoints(points);
      helix->CalculateNregularPoints();
      
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
    /*
    for(double pz = -maxPz; pz <= -minPz ; pz+=config->GetStepPz()){
      double c = Point(circle->GetMomentum()->GetX(), circle->GetMomentum()->GetY(), pz).GetVectorSlopeC();
      unique_ptr<Helix> helix = make_unique<Helix>(c, circle, config);
      helix->SetPoints(points);
      helix->CalculateNregularPoints();
      
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
     */
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
