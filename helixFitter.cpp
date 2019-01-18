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

int verbosityLevel = 2;

string configPath = "configs/helixFitter.md";
shared_ptr<FitterConfig> config;
unique_ptr<PointsProcessor> pointsProcessor;

// Will be calculated automatically
double trackEta, trackTheta, trackPhi; // parameters of the chargino track
double decayR;  // secondary vertex R (from 0,0,0) just somewhere between 3rd and 4th layer
double minPx, minPy, minPz, maxPx, maxPy, maxPz, minL, maxL;

// Monitoring histograms
map<string, TH1D*> monitors1D;
map<string, TH2D*> monitors2D;

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

/// Returns 0 in case no helix was fitter, 1 it there was a helix fitted but it was far from the true one
/// and 2 if a perfectly correct helix was found
int FillMonitors(const unique_ptr<Helix> &fittedHelix, const unique_ptr<Helix> &trueHelix)
{
  int success = 0;
  if(!fittedHelix){
    monitors1D["failReason"]->Fill(8);
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
  
  if(failureCodes.size()==0) success = 2;
  else{
    success = 1;
    for(int f : failureCodes) monitors1D["failReason"]->Fill(f);
  }
  return success;
}

unique_ptr<Helix> GetRandomPionHelix()
{
  unique_ptr<Point> pionVector = make_unique<Point>(RandSign()*RandDouble(minPx, maxPx),
                                                    RandSign()*RandDouble(minPy, maxPy),
                                                    RandSign()*RandDouble(minPz, maxPz));
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

void PlotAndSaveMonitors()
{
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
  c1->Update();
}

void PerformTests(int &nSuccess, int &nFullSuccess)
{
  int nTests = config->GetNtests();
  unique_ptr<Fitter> fitter = make_unique<Fitter>(config);
  
  TH1D *success_vs_pz = new TH1D("success_vs_pz","success_vs_pz",config->GetMaxPz()-config->GetMinPz()+1,config->GetMinPz(),config->GetMaxPz());
  TH1D *full_success_vs_pz = new TH1D("full_success_vs_pz","full_success_vs_pz",config->GetMaxPz()-config->GetMinPz()+1,config->GetMinPz(),config->GetMaxPz());
  
  success_vs_pz->Sumw2();
  TH1D *den_vs_pz = new TH1D("den_vs_pz","den_vs_pz",config->GetMaxPz()-config->GetMinPz()+1,config->GetMinPz(),config->GetMaxPz());
  den_vs_pz->Sumw2();
  
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
    int successCode = FillMonitors(bestHelix, pionHelix);
    
    den_vs_pz->Fill(fabs(pionHelix->GetMomentum()->GetZ()));
    
    if(successCode == 1){
      nSuccess++;
      success_vs_pz->Fill(fabs(pionHelix->GetMomentum()->GetZ()));
    }
    else if( successCode == 2){
      nSuccess++;
      nFullSuccess++;
      success_vs_pz->Fill(fabs(pionHelix->GetMomentum()->GetZ()));
      full_success_vs_pz->Fill(fabs(pionHelix->GetMomentum()->GetZ()));
    }
    
    if(bestHelix){
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
  
  success_vs_pz->Divide(den_vs_pz);
  full_success_vs_pz->Divide(den_vs_pz);
  
  TCanvas *c2 = new TCanvas("c2","c2",1280,800);
  c2->Divide(2,2);
  c2->cd(1);
  success_vs_pz->Draw();
  c2->cd(2);
  full_success_vs_pz->Draw();
  c2->Update();
  success_vs_pz->SaveAs("success_vs_pz.root");
  full_success_vs_pz->SaveAs("full_success_vs_pz.root");
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
  SetupMonitors();
  
  auto startTime = now();
  
  int nTests = config->GetNtests();
  int nSuccess, nFullSuccess;
  nFullSuccess = nSuccess = 0;
  PerformTests(nSuccess, nFullSuccess);
  cout<<"Percentage of successful fits:"<<nSuccess/(double)nTests<<endl;
  cout<<"Percentage of fully successful fits:"<<nFullSuccess/(double)nTests<<endl;
  PlotAndSaveMonitors();
  
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
