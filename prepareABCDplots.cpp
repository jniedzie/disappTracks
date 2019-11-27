//  prepareABCDplots.cpp
//
//  Created by Jeremi Niedziela on 20/12/2019.

#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"
#include "Logger.hpp"

typedef tuple<vector<double>, vector<double>> binning;

// Desired number of MET and dE/dx bins and limits of those
int nDedxBins = 3, nMetBins  = 3;
const double minMet  = 200 , maxMet  = 500 , stepMet  = 10;
//const double minDedx = 2.0 , maxDedx = 5.1 , stepDedx = 0.1; // for min dE/dx
const double minDedx = 2.5 , maxDedx = 11.0 , stepDedx = 0.1; // for dE/dx likelihood

const double limitMet = 1000.0;
const double limitDedx = 20.0;

string configPath = "configs/analysis.md";
string outputPath = "results/abcd_optimization.txt";

string afsPath = "/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/combine/CMSSW_10_2_13/src/";

bool simulateTagger = false;
//double taggerEfficiency = 0.595152; // no PU
//double taggerFakeRate   = 0.119221;
double taggerEfficiency = 0.852; // with PU
double taggerFakeRate   = 0.23;


const int ratioRebin = 2;
string sampleTag = "_Wmunu";
string backgroundHistNams = "background";
string dataHistNams = "data";


//------------------------------------------------
// 2x2, 3 layers
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
  { kWino_M_300_cTau_3    ,  {{430}, {2.4}}},
  { kWino_M_300_cTau_10   ,  {{370}, {2.4}}},
  { kWino_M_300_cTau_30   ,  {{320}, {2.8}}},
  { kWino_M_500_cTau_10   ,  {{450}, {2.9}}},
  { kWino_M_500_cTau_20   ,  {{450}, {3.0}}},
  { kWino_M_650_cTau_10   ,  {{440}, {2.9}}},
  { kWino_M_650_cTau_20   ,  {{450}, {3.0}}},
  { kWino_M_800_cTau_10   ,  {{450}, {3.0}}}, // BEST
  { kWino_M_800_cTau_20   ,  {{450}, {3.0}}},
  { kWino_M_1000_cTau_10  ,  {{450}, {3.6}}},
  { kWino_M_1000_cTau_20  ,  {{450}, {3.6}}},
};
*/
//------------------------------------------------
// 2x2, 4 layers
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
  { kWino_M_300_cTau_3    ,  {{400}, {3.8}}},
  { kWino_M_300_cTau_10   ,  {{320}, {2.1}}},
  { kWino_M_300_cTau_30   ,  {{320}, {2.1}}},
  { kWino_M_500_cTau_10   ,  {{410}, {2.6}}},
  { kWino_M_500_cTau_20   ,  {{400}, {2.8}}},
  { kWino_M_650_cTau_10   ,  {{490}, {2.4}}},
  { kWino_M_650_cTau_20   ,  {{400}, {2.6}}},
  { kWino_M_800_cTau_10   ,  {{490}, {2.4}}},
  { kWino_M_800_cTau_20   ,  {{440}, {2.5}}}, // BEST
  { kWino_M_1000_cTau_10  ,  {{490}, {2.4}}},
  { kWino_M_1000_cTau_20  ,  {{490}, {2.4}}},
};
*/


//------------------------------------------------
// 2x3, 3 layers
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
 { kWino_M_300_cTau_3   , {{300, 430}, {2.4}}},
 { kWino_M_300_cTau_10  , {{300, 450}, {2.5}}},
 { kWino_M_300_cTau_30  , {{300, 410}, {2.8}}},
 { kWino_M_500_cTau_10  , {{320, 450}, {2.9}}},
 { kWino_M_500_cTau_20  , {{300, 450}, {3.0}}},
 { kWino_M_650_cTau_10  , {{300, 440}, {2.9}}},
 { kWino_M_650_cTau_20  , {{300, 450}, {3.0}}},
 { kWino_M_800_cTau_10  , {{300, 450}, {3.0}}},
 { kWino_M_800_cTau_20  , {{300, 450}, {3.0}}}, // BEST
 { kWino_M_1000_cTau_10 , {{300, 450}, {3.6}}},
 { kWino_M_1000_cTau_20 , {{300, 450}, {3.5}}},
};
*/
//------------------------------------------------
// 2x3, 4 layers
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
 { kWino_M_300_cTau_3    , {{400, 470}, {3.8}}},
 { kWino_M_300_cTau_10   , {{300, 400}, {2.1}}},
 { kWino_M_300_cTau_30   , {{320, 400}, {2.1}}},
 { kWino_M_500_cTau_10   , {{310, 440}, {2.5}}},
 { kWino_M_500_cTau_20   , {{320, 440}, {2.6}}},
 { kWino_M_650_cTau_10   , {{310, 490}, {2.5}}},
 { kWino_M_650_cTau_20   , {{320, 440}, {2.5}}},
 { kWino_M_800_cTau_10   , {{400, 490}, {2.4}}}, // BEST
 { kWino_M_800_cTau_20   , {{320, 490}, {2.5}}},
 { kWino_M_1000_cTau_10  , {{320, 490}, {2.4}}},
 { kWino_M_1000_cTau_20  , {{320, 490}, {2.4}}},
};
*/



//------------------------------------------------
// 3x2, 3 layers
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
 { kWino_M_300_cTau_3   , {{430}, {2.4, 5.1}}},
 { kWino_M_300_cTau_10  , {{370}, {2.1, 3.6}}},
 { kWino_M_300_cTau_30  , {{390}, {2.8, 3.7}}},
 { kWino_M_500_cTau_10  , {{450}, {3.0, 5.0}}}, // BEST
 { kWino_M_500_cTau_20  , {{450}, {3.0, 4.0}}},
 { kWino_M_650_cTau_10  , {{440}, {4.8, 5.0}}},
 { kWino_M_650_cTau_20  , {{450}, {3.0, 5.1}}},
 { kWino_M_800_cTau_10  , {{450}, {3.6, 4.0}}},
 { kWino_M_800_cTau_20  , {{460}, {4.8, 5.1}}},
 { kWino_M_1000_cTau_10 , {{440}, {4.8, 5.0}}},
 { kWino_M_1000_cTau_20 , {{460}, {4.8, 5.1}}},
};
*/
//------------------------------------------------
// 3x2, 4 layers
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
 { kWino_M_300_cTau_3    , {{400}, {2.1, 3.8}}},
 { kWino_M_300_cTau_10   , {{320}, {2.1, 5.1}}},
 { kWino_M_300_cTau_30   , {{320}, {2.1, 3.9}}},
 { kWino_M_500_cTau_10   , {{390}, {2.6, 2.7}}},
 { kWino_M_500_cTau_20   , {{320}, {2.6, 4.4}}},
 { kWino_M_650_cTau_10   , {{440}, {2.6, 2.8}}},
 { kWino_M_650_cTau_20   , {{490}, {2.6, 2.8}}},
 { kWino_M_800_cTau_10   , {{450}, {2.6, 3.2}}},
 { kWino_M_800_cTau_20   , {{480}, {2.5, 3.0}}},
 { kWino_M_1000_cTau_10  , {{440}, {2.5, 3.6}}}, // BEST
 { kWino_M_1000_cTau_20  , {{440}, {2.5, 3.6}}},
};
*/
 
 

//------------------------------------------------
// 3x3, 3 layers, likelihood
//------------------------------------------------

map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
//  { kWino_M_300_cTau_3    , {{330, 440}, {4.6}}},
//  { kWino_M_300_cTau_10   , {{300, 370}, {4.0}}},
//  { kWino_M_300_cTau_30   , {{320, 410}, {4.0}}},
//  { kWino_M_500_cTau_10   , {{320, 450}, {4.0}}},
//  { kWino_M_500_cTau_20   , {{330, 450}, {4.0}}},
//  { kWino_M_650_cTau_10   , {{300, 440}, {4.0}}},
//  { kWino_M_650_cTau_20   , {{320, 450}, {4.0}}},
//  { kWino_M_800_cTau_10   , {{300, 440}, {4.0}}},
//  { kWino_M_800_cTau_20   , {{350, 450}, {4.0}}},
//  { kWino_M_1000_cTau_10  , {{350, 450}, {4.0}}},
  { kWino_M_1000_cTau_20  , {{390, 490}, {4.6}}}, // BEST
};

//------------------------------------------------
// 3x3, 4 layers, likelihood
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
//  { kWino_M_300_cTau_3    , {{300, 470}, {3.9, 5.3}}},
//  { kWino_M_300_cTau_10   , {{320, 400}, {4.0}}},
//  { kWino_M_300_cTau_30   , {{320, 400}, {4.0}}},
//  { kWino_M_500_cTau_10   , {{300, 440}, {4.0}}},
//  { kWino_M_500_cTau_20   , {{310, 440}, {4.0}}},
//  { kWino_M_650_cTau_10   , {{310, 490}, {4.0}}},
//  { kWino_M_650_cTau_20   , {{320, 490}, {4.0}}},
//  { kWino_M_800_cTau_10   , {{320, 410}, {3.7, 6.3}}},
  { kWino_M_800_cTau_20   , {{320, 470}, {4.1, 7.1}}},  // BEST
//  { kWino_M_1000_cTau_10  , {{310, 350}, {4.1, 6.9}}},
//  { kWino_M_1000_cTau_20  , {{320, 490}, {5.0}}},
};
*/
//------------------------------------------------
// 3x3, 5-6 layers, likelihood
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
//  { kWino_M_300_cTau_3    , {{300, 390}, {3.7, 4.5}}},
//  { kWino_M_300_cTau_10   , {{300, 340}, {4.0}}},
//  { kWino_M_300_cTau_30   , {{300, 340}, {4.0}}},
//  { kWino_M_500_cTau_10   , {{340, 470}, {4.0}}},
//  { kWino_M_500_cTau_20   , {{300, 340}, {4.0}}},
//  { kWino_M_650_cTau_10   , {{340, 450}, {3.7, 5.7}}},
//  { kWino_M_650_cTau_20   , {{340, 470}, {4.0}}},
  { kWino_M_800_cTau_10   , {{350, 450}, {3.7, 5.3}}}, // BEST
//  { kWino_M_800_cTau_20   , {{340, 470}, {4.0}}},
//  { kWino_M_1000_cTau_10  , {{300, 360}, {4.3, 5.3}}},
//  { kWino_M_1000_cTau_20  , {{300, 350}, {4.3, 5.3}}},
};
*/

/**
 Returns number of counts in ABCD... regions determined by criticalMet and criticalDedx values.
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param bestValues Positions of bins border on MET and dE/dx axes
 */
vector<vector<double>> GetABCD(const TH2D *metVsDedxHist, const binning bestValues, double scale)
{
  auto &[criticalMet, criticalDedx] = bestValues;
  vector<double> binsMet = { 0 };
  vector<double> binsDedx = { 0 };
  
  for(float met : criticalMet) binsMet.push_back(met);
  for(float dedx : criticalDedx) binsDedx.push_back(dedx);
  
  binsMet.push_back(inf);
  binsDedx.push_back(inf);
  
  vector<vector<double>> abcd(binsMet.size()-1, vector<double>(binsDedx.size()-1));
    
  for(int iMet=0; iMet<binsMet.size()-1; iMet++){
    for(int iDedx=0; iDedx<binsDedx.size()-1; iDedx++){
      
      int binX1 = metVsDedxHist->GetXaxis()->FindFixBin(binsDedx[iDedx]+stepDedx/2.);
      int binX2 = metVsDedxHist->GetXaxis()->FindFixBin(binsDedx[iDedx+1]-stepDedx/2.);
      
      int binY1 = metVsDedxHist->GetYaxis()->FindFixBin(binsMet[iMet]+stepMet/2.);
      int binY2 = metVsDedxHist->GetYaxis()->FindFixBin(binsMet[iMet+1]-stepMet/2.);
      
      abcd[iMet][iDedx] = scale*metVsDedxHist->Integral(binX1, binX2, binY1, binY2);
    }
  }
  
  return abcd;
}

/**
 Creates ABCD... histogram. In case of backgrounds, all samples are merged into the same histogram.
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param bestValues Positions of bins border on MET and dE/dx axes
 \param dataType Specifies whether background or signal events should be analyzed
 \param setIter For signal, specified which samples to use
 */
TH2D* GetABCDplot(const TH2D* metVsDedxHist, const binning bestValues,
                  xtracks::EDataType dataType, int setIter=0)
{
  auto &[criticalMet, criticalDedx] = bestValues;
  if(criticalMet.size()==0 || criticalDedx.size()==0){
    cout<<"No critical MET or dE/dx value were privided!"<<endl;
    return nullptr;
  }
  
  string title = "ABCD_"+to_string_with_precision(criticalMet[0], 0)+"_"+to_string_with_precision(criticalDedx[0], 1);
  double scale = 1.0;
  
  if(dataType == xtracks::kSignal){
    title += ("_"+signalName.at((ESignal)setIter));
    if(simulateTagger) scale = taggerEfficiency;
  }
  else if(dataType == kData){
    title += "_Data";
    if(simulateTagger) scale = taggerFakeRate;
  }
  else{
    title += "_Backgrounds";
    if(simulateTagger) scale = taggerFakeRate;
  }
  
  vector<float> binsMet = { 0 };
  vector<float> binsDedx = { 0 };
  
  for(float met : criticalMet) binsMet.push_back(met);
  for(float dedx : criticalDedx) binsDedx.push_back(dedx);
  
  binsMet.push_back(limitMet);
  binsDedx.push_back(limitDedx);
  
  float binsMetArray[100];
  float binsDedxArray[100];
  
  for(int i=0; i<binsMet.size(); i++) binsMetArray[i] = binsMet[i];
  for(int i=0; i<binsDedx.size(); i++) binsDedxArray[i] = binsDedx[i];
  
  TH2D *hist = new TH2D(title.c_str(), title.c_str(),
                        (int)binsDedx.size()-1, binsDedxArray,
                        (int)binsMet.size()-1, binsMetArray);
 
  auto abcd = GetABCD(metVsDedxHist, bestValues, scale);
  
  for(int iMet=0; iMet<binsMet.size()-1; iMet++){
    for(int iDedx=0; iDedx<binsDedx.size()-1; iDedx++){
      float midMet = (binsMet[iMet+1] + binsMet[iMet]) / 2.;
      float midDedx = (binsDedx[iDedx+1] + binsDedx[iDedx]) / 2.;
      hist->Fill(midDedx, midMet, abcd[iMet][iDedx]);
    }
  }
  return hist;
}

double GetSignificance(const TH2D *metVsDedxHistBackground, const TH2D *metVsDedxHistSignal,
                       const binning bestValues)
{
  auto abcdBackground = GetABCD(metVsDedxHistBackground, bestValues, simulateTagger ? taggerFakeRate : 1.0);
  auto abcdSignal     = GetABCD(metVsDedxHistSignal, bestValues, simulateTagger ? taggerEfficiency : 1.0);
  
  double significance = 0;
  
  for(int x=0; x<abcdBackground.size(); x++){
    for(int y=0; y<abcdBackground[x].size(); y++){
      if(abcdBackground[x][y] == 0) return 0; // protection agains bins with no background events
      if(abcdSignal[x][y] + abcdBackground[x][y] == 0) continue;
      
      // Simplified version
      double signif = abcdSignal[x][y] / sqrt(abcdSignal[x][y] + abcdBackground[x][y]);
      
      // Full version (from Andrea)
      double signifFull = sqrt(2*( (abcdSignal[x][y]+abcdBackground[x][y]) * log(1+abcdSignal[x][y]/abcdBackground[x][y])-abcdSignal[x][y]));
      
      significance += pow(signifFull, 2);
    }
  }
  return sqrt(significance);
}

/// Fills correlation histograms from background events
TH2D* GetMetVsDedxHist(const EventSet &events, xtracks::EDataType dataType, int setIter=0)
{
  TH2D *hist = new TH2D("metVsDedx","metVsDedx",1000, 0.0, 100.0, 1000, 0, 10000);
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    if(dataType == xtracks::kBackground){
      for(EBackground iBck : backgrounds){
        if(!config.runBackground[iBck]) continue;
        
        for(int iEvent=0; iEvent<events.size(dataType, iBck, year); iEvent++){
          auto event = events.At(dataType, iBck, year, iEvent);
          
          for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
            auto track = event->GetTrack(iTrack);
            //          hist->Fill(track->GetMinDedx(), event->GetMetNoMuPt(), event->GetWeight());
            hist->Fill(track->GetDedxLikelihood(), event->GetMetNoMuPt(), event->GetWeight());
          }
        }
      }
      
    }
    else if(dataType == kData){
      for(EData iData : datas){
        if(!config.runData[iData]) continue;
        
        for(int iEvent=0; iEvent<events.size(dataType, iData, year); iEvent++){
          auto event = events.At(dataType, iData, year, iEvent);
          
          for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
            auto track = event->GetTrack(iTrack);
            //          hist->Fill(track->GetMinDedx(), event->GetMetNoMuPt(), event->GetWeight());
            hist->Fill(track->GetDedxLikelihood(), event->GetMetNoMuPt(), event->GetWeight());
          }
        }
      }
    }
    else if(dataType == xtracks::kSignal){
      for(int iEvent=0;iEvent<events.size(dataType, setIter, year);iEvent++){
        auto event = events.At(dataType, setIter, year, iEvent);
        
        for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
          auto track = event->GetTrack(iTrack);
          //        hist->Fill(track->GetMinDedx(), event->GetMetNoMuPt(), event->GetWeight());
          hist->Fill(track->GetDedxLikelihood(), event->GetMetNoMuPt(), event->GetWeight());
        }
      }
      
    }
  }
  return hist;
}

TGraphErrors* GetRatioGraph(const TH2D *metVsDedxHistBackground, const binning bestValues)
{
  auto &[criticalMet, criticalDedx] = bestValues;
  
  TH2D *metVsDedxRebinned = new TH2D(*metVsDedxHistBackground);
  metVsDedxRebinned->RebinX(ratioRebin);
  
  TGraphErrors *abcdBackgroundRatio = new TGraphErrors();
  int iPoint=0;
  
  double min = metVsDedxRebinned->GetYaxis()->GetBinLowEdge(1);
  double max = metVsDedxRebinned->GetYaxis()->GetBinUpEdge(metVsDedxRebinned->GetNbinsY());
  
  for(int iMet=0; iMet<criticalMet.size(); iMet++){
    double minMet, middleMet, maxMet;
    
    if(iMet==0) minMet = min;
    else        minMet = criticalMet[iMet-1];
    
    if(iMet==criticalMet.size())  middleMet = max;
    else                          middleMet = criticalMet[iMet];
    
    if((iMet+1)==criticalMet.size())  maxMet = max;
    else                              maxMet = criticalMet[iMet+1];
    
    range<double> metRangeNum(minMet, middleMet);
    range<double> metRangeDen(middleMet, maxMet);
    
    TH1D *histNum = new TH1D("histNum", "histNum",
                             metVsDedxRebinned->GetNbinsX(),
                             metVsDedxRebinned->GetXaxis()->GetXmin(),
                             metVsDedxRebinned->GetXaxis()->GetXmax());
    histNum->Sumw2();
    TH1D *histDen = new TH1D(*histNum);
    
    for(int binX=1; binX<=metVsDedxRebinned->GetNbinsX(); binX++){
      double dedx = metVsDedxRebinned->GetXaxis()->GetBinCenter(binX);
      
      for(int binY=1; binY<=metVsDedxRebinned->GetNbinsY(); binY++){
        double met = metVsDedxRebinned->GetYaxis()->GetBinCenter(binY);
        
        double nEvents = metVsDedxRebinned->GetBinContent(binX, binY);
        
        if(metRangeNum.IsInside(met)) histNum->Fill(dedx, nEvents);
        if(metRangeDen.IsInside(met)) histDen->Fill(dedx, nEvents);
      }
    }
    
    histNum->Divide(histDen);
    
    for(int binX=1; binX<=histNum->GetNbinsX(); binX++){
      double dedx = histNum->GetXaxis()->GetBinCenter(binX);

      if(histNum->GetBinContent(binX) != 0){
        abcdBackgroundRatio->SetPoint(iPoint, dedx, histNum->GetBinContent(binX));
        abcdBackgroundRatio->SetPointError(iPoint, 0, histNum->GetBinError(binX));
        iPoint++;
      }
    }
  }
  return abcdBackgroundRatio;
}


/**
 Draws and saves ABCD plots for given values of critical MET and critical dE/dx
 \param metVsDedxHistBackground MET vs. dE/dx histogram with backgrounds
 \param metVsDedxHistsSignal map of MET vs. dE/dx histograms for each signal
 \param bins Tuple of vectors with bin possitions for MET and dE/dx, respectively
 \param outputPath Path to which all histograms will be saved
 */
void drawAndSaveABCDplots(const TH2D *metVsDedxHistBackground,
                          const map<int, TH2D*> &metVsDedxHistsSignal,
                          const binning bins, string outputPath,
                          const TH2D *metVsDedxHistData = nullptr)
{
  TCanvas *abcdCanvas = new TCanvas("ABCD", "ABCD", 1000, 1500);
  abcdCanvas->Divide(4,3);
  TFile *outFile = new TFile(outputPath.c_str(),"recreate");
  
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".1f");
  
  cout<<"Plotting background ABCD"<<endl;
  metVsDedxHistBackground->Write();
  TH2D *abcdPlotBackgrounds = GetABCDplot(metVsDedxHistBackground, bins, kBackground);
  abcdPlotBackgrounds->SetMarkerSize(3.0);
  
  TCanvas *ratioCanvas = new TCanvas("Ratio", "Ratio", 1200, 600);
  ratioCanvas->Divide(2, 1);
  
  if(metVsDedxHistData){
    metVsDedxHistData->Write();
    TH2D *abcdPlotData = GetABCDplot(metVsDedxHistData, bins, kData);
    outFile->cd();
    abcdPlotData->SetName(dataHistNams.c_str());
    abcdPlotData->SetTitle(dataHistNams.c_str());
    abcdPlotData->Write();
  
    ratioCanvas->cd(1);
    
    TGraphErrors *abcdDataRatio = GetRatioGraph(metVsDedxHistData, bins);
    abcdDataRatio->SetMarkerStyle(20);
    abcdDataRatio->SetMarkerSize(1.0);
    abcdDataRatio->SetMarkerColor(kViolet);
    abcdDataRatio->GetXaxis()->SetTitle("dE/dx likelihood");
    abcdDataRatio->GetYaxis()->SetTitle("Ratio");
    
    TF1 *linearFunction = new TF1("linearFunction", "[0]", 3, 10);
    linearFunction->SetParameter(0, 1.0);
    linearFunction->SetLineColor(kGreen);
    
    TF1 *linearFunctionWithTilt = new TF1("linearFunctionWithTilt", "[1]+x*[0]", 3, 10);
    linearFunctionWithTilt->SetParameter(0, 1.0);
    linearFunctionWithTilt->SetParameter(1, 0.0);
    linearFunctionWithTilt->SetLineColor(kRed);
    abcdDataRatio->Draw("APE");
    
    abcdDataRatio->Fit(linearFunction, "", "", 3, 10);
    linearFunction->Draw("sameL");
    
    abcdDataRatio->Fit(linearFunctionWithTilt, "", "", 3, 10);
    linearFunctionWithTilt->Draw("sameL");
    
  }
  
  ratioCanvas->cd(2);
  
  TGraphErrors *abcdBackgroundRatio = GetRatioGraph(metVsDedxHistBackground, bins);
  abcdBackgroundRatio->SetMarkerStyle(20);
  abcdBackgroundRatio->SetMarkerSize(1.0);
  abcdBackgroundRatio->SetMarkerColor(kViolet);
  abcdBackgroundRatio->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
  abcdBackgroundRatio->GetYaxis()->SetTitle("Ratio");
  
  TF1 *linearFunction = new TF1("linearFunction", "[0]", 3, 10);
  linearFunction->SetParameter(0, 1.0);
  linearFunction->SetLineColor(kGreen);
  
  TF1 *linearFunctionWithTilt = new TF1("linearFunctionWithTilt", "[1]+x*[0]", 3, 10);
  linearFunctionWithTilt->SetParameter(0, 1.0);
  linearFunctionWithTilt->SetParameter(1, 0.0);
  linearFunctionWithTilt->SetLineColor(kRed);
  abcdBackgroundRatio->Draw("APE");
  
  abcdBackgroundRatio->Fit(linearFunction, "", "", 3, 10);
  linearFunction->Draw("sameL");
  
  abcdBackgroundRatio->Fit(linearFunctionWithTilt, "", "", 3, 10);
  linearFunctionWithTilt->Draw("sameL");
  
  ratioCanvas->Update();
  ratioCanvas->SaveAs("plots/abcdRatio.pdf");
  
  int iPad=1;
  abcdCanvas->cd(iPad++);
  abcdPlotBackgrounds->Draw("colzText");
  
  outFile->cd();
  abcdPlotBackgrounds->SetName(backgroundHistNams.c_str());
  abcdPlotBackgrounds->SetTitle(backgroundHistNams.c_str());
  abcdPlotBackgrounds->Write();
  abcdBackgroundRatio->Write();
  
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    cout<<"Plotting "<<signalTitle.at(iSig)<<" ABCD"<<endl;
    TH2D *abcdPlot = GetABCDplot(metVsDedxHistsSignal.at((ESignal)iSig), bins, xtracks::kSignal, iSig);
    abcdCanvas->cd(iPad++);
    abcdPlot->SetMarkerSize(3.0);
    abcdPlot->Draw("colzText");
    outFile->cd();
    abcdPlot->SetName(signalName.at(iSig).c_str());
    abcdPlot->SetTitle(signalName.at(iSig).c_str());
    abcdPlot->Write();
  }
  
  abcdCanvas->Update();
  abcdCanvas->Write();
  abcdCanvas->SaveAs("plots/abcd.pdf");
  outFile->Close();
}

/**
 Adds all combinations of numbers between the last element of each of the vectors in `groups` and `max`.
 \param groups Vector of vectors of lines, which will be extended by all possible combinations.
 \param max Maximum value of the numbers to be put in `gorups` vector
 \param step Step between the numbers to be added to `groups`
 \param nLines Desired total number of values in each vector
 */
void AddValuesCombinations(vector<vector<double>> &groups, double max, double step, int nLines)
{
  if(groups[0].size()==nLines) return;
  
  vector<vector<double>> newGroups;
  
  for(vector<double> &group : groups){
    vector<vector<double>> replacementForGroup;
    
    for(double value=group[group.size()-1]+step; value<max; value+=step){
      vector<double> newGroup;
      for(double element : group) newGroup.push_back(element);
      newGroup.push_back(value);
      replacementForGroup.push_back(newGroup);
    }
    for(auto element : replacementForGroup) newGroups.push_back(element);
  }
  groups = newGroups;
  
  AddValuesCombinations(groups, max, step, nLines);
}

/**
 Finds the best positions of dE/dx bin borders, optimizing significance
 \param metVsDedxHistBackground MET vs. dE/dx histogram for background
 \param metVsDedxHistSignal MET vs. dE/dx histogram for signal
 */
tuple<vector<double>,vector<double>> findBestBinning(const TH2D *metVsDedxHistBackground,
                                                     const TH2D *metVsDedxHistSignal,
                                                     const vector<vector<double>> &groupsDedx,
                                                     const vector<vector<double>> &groupsMet)
{
 
  
  cout<<"Combinations in vector: "<<groupsDedx.size()*groupsDedx[0].size()*groupsMet.size()*groupsMet[0].size()<<endl;
 
  // Find the best combination
  vector<double> bestMet, bestDedx;
  double bestSignificance = -inf;
  
  Log(2)<<"\t|\t";
  for(auto groupDedx : groupsDedx){ for(double dedx : groupDedx) Log(2)<<dedx<<"\t";}
  cout<<endl;
  
  for(auto groupMet : groupsMet){
    for(double met : groupMet) cout<<met<<"\t";
    Log(2)<<"|\t";
  
    for(auto groupDedx : groupsDedx){
      double significance = GetSignificance(metVsDedxHistBackground, metVsDedxHistSignal, {groupMet, groupDedx});
      Log(2)<<significance<<" ";
      
      if(significance > bestSignificance){
        bestSignificance = significance;
        bestMet = groupMet;
        bestDedx = groupDedx;
      }
    }
    cout<<endl;
  }
  return make_tuple(bestMet, bestDedx);
}

/// Creates Combine datacard using Andrea's tool
void createDatacard(string outFileName, string outputPath)
{
  string command = "/usr/local/bin/python3 ../DatacardCreatorABCD/mkDatacards.py  --inputHistoFile "+outputPath;
  command += " --dataHistoName  "+backgroundHistNams;
  command += " --sigHistoNameTemplate  Wino";
  command += " --nuisancesFile   test/nuisances.py";
  
  cout<<"Executing command:\n"<<command<<endl;
  
  string commandOutput = exec(command.c_str());
  cout<<commandOutput<<endl;
  
  command = "cp test/datacard_mytest.txt ";
  command += "results/datacards/datacard_"+outFileName+".txt";
  exec(command.c_str());
}

/// Transfers datacard created by this app to a hardcoded lxplus location using scp
void copyDatacardToLxplus(string outFileName)
{
  string command = "scp results/datacards/datacard_"+outFileName+".txt";
  command += " jniedzie@lxplus.cern.ch:/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/combine/datacards/andrea/";
  exec(command.c_str());
}

/// Runs Combine on lxplus and transfers results back to local machine
void runCombine(string outFileName)
{
  string command = "ssh jniedzie@lxplus7.cern.ch './runCombine.sh datacard_"+outFileName+"'";
  string output = exec(command.c_str());
  cout<<output<<endl;
  
  command = "cp "+afsPath+"limits_datacard_"+outFileName+".txt";
  command += " macros/limitsData/combineOutput/";
  exec(command.c_str());
}

void convertRtoLimits(string outFileName)
{
  string command = "/usr/local/Cellar/root/6.18.04/bin/root -q -b -l ";
  command += "\"macros/getLimitsFromR.C(\\\"macros/limitsData/combineOutput/limits_datacard_"+outFileName+".txt\\\", ";
  command += "\\\"macros/limitsData/cms_short_disappearing_"+outFileName+".txt\\\")\"";
  exec(command.c_str());
}

map<int, TH2D*> loadSignalHists(const EventSet &events)
{
  map<int, TH2D*> metVsDedxHistsSignal;
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    metVsDedxHistsSignal[iSig] = GetMetVsDedxHist(events, xtracks::kSignal, iSig);
  }
  return metVsDedxHistsSignal;
}

void runBinningScan(const TH2D *metVsDedxHistBackground, const map<int, TH2D*> &metVsDedxHistsSignal)
{
  // Find all combinations of binning in dE/dx and MET
  vector<vector<double>> groupsDedx;
  for(double startingDedx=minDedx; startingDedx<maxDedx; startingDedx+=stepDedx) groupsDedx.push_back({startingDedx});
  AddValuesCombinations(groupsDedx, maxDedx, stepDedx, nDedxBins-1);

  vector<vector<double>> groupsMet;
  for(double startingMet=minMet; startingMet<maxMet; startingMet+=stepMet) groupsMet.push_back({startingMet});
  AddValuesCombinations(groupsMet, maxMet, stepMet, nMetBins-1);
  
  ofstream outFile(outputPath);
  
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    auto result = findBestBinning(metVsDedxHistBackground, metVsDedxHistsSignal.at(iSig), groupsDedx, groupsMet);
    double significance = GetSignificance(metVsDedxHistBackground, metVsDedxHistsSignal.at(iSig), result);

    Log(0)<<"Category: "<<config.category<<"\n";
    Log(0)<<"Sample: "<<signalTitle.at(iSig)<<"\n";
    Log(0)<<"MET bins: "; for(double met : get<0>(result)) Log(0)<<met<<"\t"; Log(0)<<"\n";
    Log(0)<<"dE/dx bins: "; for(double dedx : get<1>(result)) Log(0)<<dedx<<"\t"; Log(0)<<"\n";
    Log(0)<<"significance: "<<significance<<"\n";
    
    outFile<<"Category: "<<config.category<<"\n";
    outFile<<"Sample: "<<signalTitle.at(iSig)<<"\n";
    outFile<<"MET bins: "; for(double met : get<0>(result)) outFile<<met<<"\t"; outFile<<"\n";
    outFile<<"dE/dx bins: "; for(double dedx : get<1>(result)) outFile<<dedx<<"\t"; outFile<<"\n";
    outFile<<"significance: "<<significance<<"\n";
  }
  outFile.close();
}

void getLimitsForSignal(string outFileName, string outputPath)
{
  cout<<"\n\n--------------------------------------------------------"<<endl;
  cout<<"Creating datacard"<<endl;
  createDatacard(outFileName, outputPath);
  cout<<"\n\n--------------------------------------------------------"<<endl;
  cout<<"Transferring card to lxplus"<<endl;
  copyDatacardToLxplus(outFileName);
  cout<<"\n\n--------------------------------------------------------"<<endl;
  cout<<"Running Combine and copying results back to local machine"<<endl;
  runCombine(outFileName);
  cout<<"\n\n--------------------------------------------------------"<<endl;
  cout<<"Converting signal strength R to limits in mass-ct"<<endl;
  convertRtoLimits(outFileName);
  cout<<"\n\n--------------------------------------------------------"<<endl;
  cout<<"Done"<<endl;
}

void produceLimits(const TH2D *metVsDedxHistBackground, const map<int, TH2D*> &metVsDedxHistsSignal)
{
  vector<string> producedLimits;
  
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    if(bestValues.find((ESignal)iSig) == bestValues.end()) continue;
    
    string outFileName = to_string_with_precision(nDedxBins, 0)+"x"+to_string_with_precision(nMetBins, 0)+"_"+config.category+"_"+sampleTag+"_"+signalShortName.at(iSig);
    string outputPath = "results/abcd_plots_"+outFileName+".root";
    
    cout<<"\n\n--------------------------------------------------------"<<endl;
    cout<<"Drawing plots"<<endl;
    drawAndSaveABCDplots(metVsDedxHistBackground, metVsDedxHistsSignal,
                         bestValues.at((ESignal)iSig), outputPath);
    getLimitsForSignal(outFileName, outputPath);
    
    producedLimits.push_back("cms_short_disappearing_"+outFileName+".txt");
  }
  
  
  cout<<"\nFiles ready to add to limits plot:"<<endl;
  for(string path : producedLimits) cout<<path<<endl;
}

/// Starting point of the application
int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 6){
    cout<<"No or 2 argument expected: optimizer_output_path sample_index category nDedxBins nMetBins"<<endl;
    exit(0);
  }
  
  srand((uint)time(0));
  
  cout.imbue(locale("de_DE"));
  
  config = ConfigManager(configPath);
  
  if(argc == 6){
    outputPath = argv[1];
    for(ESignal iSig : signals) config.runSignal[iSig] = false;
    config.runSignal[atoi(argv[2])] = true;
    config.category = argv[3];
    nDedxBins = atoi(argv[4]);
    nMetBins = atoi(argv[5]);
  }
  cout<<"Output will be stored in "<<outputPath<<endl;
  
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // Load sll events with initial cuts only
  EventSet events;
  string prefix = "";
  if(config.secondaryCategory == "Wmunu") prefix += "Wmunu/";
  if(config.secondaryCategory == "Zmumu") prefix += "Zmumu/";
  
  prefix += "after_L"+to_string_with_precision(config.params["cuts_level"], 0)+"/"+config.category+"/";
  events.LoadEventsFromFiles(prefix);
  
  // Create histograms with number of events for each MET-dE/dx bin
  TH2D *metVsDedxHistBackground = GetMetVsDedxHist(events, xtracks::kBackground);
  map<int, TH2D*> metVsDedxHistsSignal = loadSignalHists(events);
  
//  runBinningScan(metVsDedxHistBackground, metVsDedxHistsSignal);
  
//  produceLimits(metVsDedxHistBackground, metVsDedxHistsSignal);

//  drawAndSaveABCDplots(metVsDedxHistBackground, metVsDedxHistsSignal, {{400},{3.0, 3.1}} , "results/abcd_plots_debug.root");
  
  
  TH2D *metVsDedxHistData = GetMetVsDedxHist(events, kData);
  
  drawAndSaveABCDplots(metVsDedxHistBackground, metVsDedxHistsSignal, {{250},{3.0}} ,
                       "results/abcd_plots_Wmunu.root", metVsDedxHistData);
  
  
  theApp->Run();
  return 0;
}

