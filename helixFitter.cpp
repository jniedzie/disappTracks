//
//  helixFitter.cpp
//
//  Created by Jeremi Niedziela on 17/12/2018.
//

#include "Helpers.hpp"
#include "Helix.hpp"
#include "Circle.hpp"
#include "Point.hpp"
#include "PointsProcessor.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "FitterConfig.hpp"
#include "MonitorsManager.hpp"

int verbosityLevel = 2;

string configPath = "configs/helixFitter.md";
shared_ptr<FitterConfig> config;
unique_ptr<PointsProcessor> pointsProcessor;
unique_ptr<MonitorsManager> monitorsManager;

// Will be calculated automatically
double trackEta, trackTheta, trackPhi; // parameters of the chargino track
double decayR;  // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer

vector<Point> LoadAllHits(uint runNumber, uint lumiSection, unsigned long long eventNumber);

void SetRandomTrack()
{
  trackEta = RandDouble(-config->GetMaxTrackEta(), config->GetMaxTrackEta());
  trackTheta = 2*atan(exp(-trackEta));
  trackPhi = RandDouble(0, 2*TMath::Pi());
  decayR = RandDouble(config->GetMinL(), config->GetMaxL());
}

unique_ptr<Helix> GetRandomPionHelix()
{
  unique_ptr<Point> pionVector = make_unique<Point>(RandSign()*RandDouble(config->GetMinPx(), config->GetMaxPx()),
                                                    RandSign()*RandDouble(config->GetMinPy(), config->GetMaxPy()),
                                                    RandSign()*RandDouble(config->GetMinPz(), config->GetMaxPz()));
  int pionCharge = RandSign();
  
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



void PerformTests(int &nSuccess, int &nFullSuccess)
{
  int nTests = config->GetNtests();
  unique_ptr<Fitter> fitter = make_unique<Fitter>(config);
    
  for(int i=0;i<nTests;i++){
    if(verbosityLevel >= 2){
      cout<<"\n========================================================"<<endl;
      cout<<"Test iter:"<<i<<endl;
    }
      
    SetRandomTrack(); // Randomly generate chargino's track
    vector<Point> pixelPoints = pointsProcessor->GetRandomPoints(config->GetNnoiseHits());
    //  vector<Point> pixelPoints = LoadAllHits(297100, 136, 245000232);
    
    unique_ptr<Helix> pionHelix = GetRandomPionHelix();
    if(config->GetInjectPionHits()) InjectPionPointsToCollectionOfPoints(pionHelix, pixelPoints);
    
    unique_ptr<Helix> bestHelix = fitter->GetBestFittingHelix(pixelPoints, trackTheta, trackPhi);
    monitorsManager->FillMonitors(bestHelix, pionHelix);
    
    auto successCode = monitorsManager->GetFittingStatus(bestHelix, pionHelix);
    
    if(successCode == MonitorsManager::kFail)         continue;
    if(successCode == MonitorsManager::kSuccess)      nSuccess++;
    if(successCode == MonitorsManager::kFullSuccess)  nFullSuccess++;
    
    if(verbosityLevel >= 2){
      cout<<"Pion helix:"; pionHelix->Print();
      cout<<"Fitted helix:"; bestHelix->Print();
      
      if(bestHelix->GetCharge() != pionHelix->GetCharge()){
        cout<<"\n\nwrong charge\n\n"<<endl;
        cout<<"best charge:"<<bestHelix->GetCharge()<<endl;
        cout<<"true charge:"<<pionHelix->GetCharge()<<endl;
        cout<<"best t shift:"<<bestHelix->GetTmin()<<endl;
        cout<<"true t shift:"<<pionHelix->GetTmin()<<endl;
      }
    }
  }
}

void ScanParameter()
{
  string paramName = "n_noise_hits";
  double paramMin = 0;
  double paramMax = 1500;
  double paramStep = 100;
  
  TH1D *eff_vs_param = new TH1D(("eff_vs_"+paramName).c_str(),
                                ("eff_vs_"+paramName).c_str(),
                                (paramMax-paramMin)/paramStep+1,paramMin,paramMax);
  
  TH1D *full_eff_vs_param = new TH1D(("full_eff_vs_"+paramName).c_str(),
                                     ("full_eff_vs_"+paramName).c_str(),
                                     (paramMax-paramMin)/paramStep+1,paramMin,paramMax);
  
  for(double param=paramMin;param<paramMax;param+=paramStep){
    cout<<"param:"<<param<<endl;
    config->nNoiseHits = param;
    
    int nTests = config->GetNtests();
    int nFullSuccess = 0;
    int nSuccess = 0;
    PerformTests(nSuccess, nFullSuccess);
    
    eff_vs_param->SetBinContent(eff_vs_param->GetXaxis()->FindFixBin(param), nSuccess/(double)nTests);
    full_eff_vs_param->SetBinContent(full_eff_vs_param->GetXaxis()->FindFixBin(param), nFullSuccess/(double)nTests);
  }
  
  TCanvas *c2 = new TCanvas("c2","c2",1280,800);
  c2->Divide(2,2);
  c2->cd(1);
  eff_vs_param->Draw();
  c2->cd(2);
  full_eff_vs_param->Draw();
  c2->Update();
  eff_vs_param->SaveAs(("eff_vs_"+paramName+".root").c_str());
  full_eff_vs_param->SaveAs(("full_eff_vs_"+paramName+".root").c_str());
  
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = make_shared<FitterConfig>(configPath);
  pointsProcessor = make_unique<PointsProcessor>();
  monitorsManager = make_unique<MonitorsManager>(config);
  
  auto startTime = now();
  
  int nTests = config->GetNtests();
  int nSuccess, nFullSuccess;
  nFullSuccess = nSuccess = 0;
  PerformTests(nSuccess, nFullSuccess);
  cout<<"Percentage of successful fits:"<<nSuccess/(double)nTests<<endl;
  cout<<"Percentage of fully successful fits:"<<nFullSuccess/(double)nTests<<endl;
  
  monitorsManager->PlotAndSaveMonitors();
  
//  ScanParameter();
  config->Print();
  
  double timeElapsed = duration(startTime, now());
  cout<<"Average time per event:"<<timeElapsed/nTests<<" seconds"<<endl;
  
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
  
  // Parameters for all hits in the pixel barrel
  const double chargeThreshold = 0; // 2000, 5000, 25000
  const double minClusterSize = 0;
  const double maxClusterSize = 100;
  
  for(uint i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    // convert cm to mm
    pixelPoints.push_back(Point(10*hitX->at(i),10*hitY->at(i),10*hitZ->at(i),hitCharge->at(i)));
  }
  
  return pixelPoints;
}
