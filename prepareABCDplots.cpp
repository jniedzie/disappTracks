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
const int nDedxBins = 3, nMetBins  = 3;
const double minMet  = 300 , maxMet  = 500 , stepMet  = 10;
const double minDedx = 2.0 , maxDedx = 5.1 , stepDedx = 0.1;

const int ratioRebin = 1;
string sampleTag = "notag";
string backgroundHistNams = "background";


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
// 3x3, 3 layers (fixed)
//------------------------------------------------
/*
map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
  { kWino_M_300_cTau_3    ,  {{330, 430},{2.4, 5.1}}},
  { kWino_M_300_cTau_10   ,  {{300, 450},{2.1, 3.6}}},
  { kWino_M_300_cTau_30   ,  {{300, 410},{2.8, 3.7}}},
  { kWino_M_500_cTau_10   ,  {{320, 450},{3.0, 5.0}}}, // BEST
  { kWino_M_500_cTau_20   ,  {{450, 480},{4.3, 5.0}}},
  { kWino_M_650_cTau_10   ,  {{300, 440},{4.8, 5.0}}},
  { kWino_M_650_cTau_20   ,  {{300, 450},{3.0, 4.7}}},
  { kWino_M_800_cTau_10   ,  {{430, 480},{3.0, 4.7}}},
  { kWino_M_800_cTau_20   ,  {{370, 460},{4.8, 5.1}}},
  { kWino_M_1000_cTau_10  ,  {{330, 440},{4.8, 5.0}}},
  { kWino_M_1000_cTau_20  ,  {{360, 460},{4.8, 5.1}}},
};
*/
//------------------------------------------------
// 3x3, 4 layers (fixed)
//------------------------------------------------

map<ESignal, binning> bestValues = { // best MET and dE/dx bins for each signal
  { kWino_M_300_cTau_3    , {{310, 400}, {2.1, 3.8}}},
  { kWino_M_300_cTau_10   , {{300, 400}, {2.1, 3.1}}},
  { kWino_M_300_cTau_30   , {{320, 400}, {2.1, 3.9}}},
  { kWino_M_500_cTau_10   , {{310, 390}, {2.6, 2.7}}},
  { kWino_M_500_cTau_20   , {{320, 490}, {2.5, 2.7}}},
  { kWino_M_650_cTau_10   , {{310, 440}, {2.6, 2.8}}},
  { kWino_M_650_cTau_20   , {{320, 490}, {2.6, 2.8}}},
  { kWino_M_800_cTau_10   , {{440, 460}, {2.5, 3.2}}},
  { kWino_M_800_cTau_20   , {{400, 480}, {2.5, 3.0}}},
  { kWino_M_1000_cTau_10  , {{320, 440}, {2.5, 3.6}}}, // BEST
  { kWino_M_1000_cTau_20  , {{320, 440}, {2.5, 3.6}}},
};




/**
 Returns number of counts in ABCD... regions determined by criticalMet and criticalDedx values.
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param bestValues Positions of bins border on MET and dE/dx axes
 */
vector<vector<double>> GetABCD(const TH2D *metVsDedxHist, const binning bestValues)
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
      
      abcd[iMet][iDedx] = metVsDedxHist->Integral(binX1, binX2, binY1, binY2);
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
TH2D* GetABCDplot(const TH2D* metVsDedxHist, const binning bestValues, xtracks::EDataType dataType, int setIter=0)
{
  auto &[criticalMet, criticalDedx] = bestValues;
  if(criticalMet.size()==0 || criticalDedx.size()==0){
    cout<<"No critical MET or dE/dx value were privided!"<<endl;
    return nullptr;
  }
  
  string title = "ABCD_"+to_string_with_precision(criticalMet[0], 0)+"_"+to_string_with_precision(criticalDedx[0], 1);
  if(dataType == xtracks::kSignal)  title += ("_"+signalName[setIter]);
  else                              title += "_Backgrounds";
  
  vector<float> binsMet = { 0 };
  vector<float> binsDedx = { 0 };
  
  for(float met : criticalMet) binsMet.push_back(met);
  for(float dedx : criticalDedx) binsDedx.push_back(dedx);
  
  binsMet.push_back(1000);
  binsDedx.push_back(6.0);
  
  float binsMetArray[100];
  float binsDedxArray[100];
  
  for(int i=0; i<binsMet.size(); i++) binsMetArray[i] = binsMet[i];
  for(int i=0; i<binsDedx.size(); i++) binsDedxArray[i] = binsDedx[i];
  
  TH2D *hist = new TH2D(title.c_str(), title.c_str(),
                        (int)binsDedx.size()-1, binsDedxArray,
                        (int)binsMet.size()-1, binsMetArray);
 
  auto abcd = GetABCD(metVsDedxHist, bestValues);
  
  for(int iMet=0; iMet<binsMet.size()-1; iMet++){
    for(int iDedx=0; iDedx<binsDedx.size()-1; iDedx++){
      float midMet = (binsMet[iMet+1] + binsMet[iMet]) / 2.;
      float midDedx = (binsDedx[iDedx+1] + binsDedx[iDedx]) / 2.;
      hist->Fill(midDedx, midMet, abcd[iMet][iDedx]);
    }
  }
  return hist;
}

/**
 Returns variance calculated from predicted/expected ratio in each 2x2 sub-histogram of the full ABCD...
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param bestValues Positions of bins border on MET and dE/dx axes
 */
tuple<double, double> GetVariance(const TH2D *metVsDedxHist, const binning bestValues)
{
  auto abcd = GetABCD(metVsDedxHist, bestValues);
  double variance = 0;
  double varianceError = 0;
  int nEntries = 0;
  for(int x=0; x<abcd.size()-1; x++){
    for(int y=0; y<abcd[x].size()-1; y++){
      if(abcd[x][y+1]==0  || abcd[x][y]==0    ||
         abcd[x+1][y]==0  || abcd[x+1][y+1]==0) return make_tuple(inf, inf);
      
      double ratio = abcd[x][y+1]/abcd[x][y] * abcd[x+1][y]/abcd[x+1][y+1];
      double ratioError = sqrt(1/(abcd[x][y+1]*abcd[x+1][y])+1/(abcd[x][y]*abcd[x+1][y+1])) * ratio;
      variance += pow(ratio-1, 2);
      varianceError += pow((2*ratio-2)*ratioError, 2);
      nEntries++;
    }
  }
  variance /= nEntries;
  
  varianceError = sqrt(varianceError);
  varianceError /= nEntries;
  
  return make_tuple(variance, varianceError);
}

double GetSignificance(const TH2D *metVsDedxHistBackground, const TH2D *metVsDedxHistSignal,
                       const binning bestValues)
{
  auto abcdBackground = GetABCD(metVsDedxHistBackground, bestValues);
  auto abcdSignal     = GetABCD(metVsDedxHistSignal, bestValues);
  
  double significance = 0;
  
  for(int x=0; x<abcdBackground.size(); x++){
    for(int y=0; y<abcdBackground[x].size(); y++){
      if(abcdSignal[x][y] + abcdBackground[x][y] == 0) continue;
      if(abcdBackground[x][y] == 0) return 0; // protection agains bins with no background events
      
      double signif = abcdSignal[x][y] / sqrt(abcdSignal[x][y] + abcdBackground[x][y]);
      significance += pow(signif, 2);
    }
  }
  return sqrt(significance);
}

/// Fills correlation histograms from background events
TH2D* GetMetVsDedxHist(const EventSet &events, xtracks::EDataType dataType, int setIter=0)
{
  TH2D *hist = new TH2D("metVsDedx","metVsDedx",1000, 0.0, 100.0, 1000, 0, 10000);
  
  if(dataType == xtracks::kBackground){
    for(int iBck=0; iBck<kNbackgrounds; iBck++){
      if(!config.runBackground[iBck]) continue;
      
      for(int iEvent=0;iEvent<events.size(dataType, iBck);iEvent++){
        auto event = events.At(dataType, iBck, iEvent);
        
        for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
          auto track = event->GetTrack(iTrack);
          hist->Fill(track->GetMinDedx(), event->GetMetNoMuPt(), event->GetWeight());
        }
      }
    }
  }
  else if(dataType == xtracks::kSignal){
    for(int iEvent=0;iEvent<events.size(dataType, setIter);iEvent++){
      auto event = events.At(dataType, setIter, iEvent);
      
      for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
        auto track = event->GetTrack(iTrack);
        hist->Fill(track->GetMinDedx(), event->GetMetNoMuPt(), event->GetWeight());
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

      abcdBackgroundRatio->SetPoint(iPoint, dedx, histNum->GetBinContent(binX));
      abcdBackgroundRatio->SetPointError(iPoint, 0, histNum->GetBinError(binX));
      iPoint++;
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
void drawAndSaveABCDplots(const TH2D *metVsDedxHistBackground, const map<int, TH2D*> &metVsDedxHistsSignal,
                          const binning bins, string outputPath)
{
  TCanvas *abcdCanvas = new TCanvas("ABCD", "ABCD", 1000, 1500);
  abcdCanvas->Divide(4,3);
  TFile *outFile = new TFile(outputPath.c_str(),"recreate");
  
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".1f");
  
  cout<<"Plotting background ABCD"<<endl;
  metVsDedxHistBackground->Write();
  TH2D *abcdPlotBackgrounds = GetABCDplot(metVsDedxHistBackground, bins, xtracks::kBackground);
  abcdPlotBackgrounds->SetMarkerSize(3.0);
  
  TCanvas *ratioCanvas = new TCanvas("Ratio", "Ratio", 800, 600);
  ratioCanvas->cd();
  
  TGraphErrors *abcdBackgroundRatio = GetRatioGraph(metVsDedxHistBackground, bins);
  abcdBackgroundRatio->SetMarkerStyle(20);
  abcdBackgroundRatio->SetMarkerSize(1.0);
  abcdBackgroundRatio->SetMarkerColor(kViolet);
  abcdBackgroundRatio->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
  abcdBackgroundRatio->GetYaxis()->SetTitle("Ratio");
  
  TF1 *linearFunction = new TF1("linearFunction", "[0]", 0, 7.0);
  linearFunction->SetParameter(0, 1.0);
  linearFunction->SetLineColor(kGreen);
  
  TF1 *linearFunctionWithTilt = new TF1("linearFunctionWithTilt", "[1]+x*[0]", 0, 7.0);
  linearFunctionWithTilt->SetParameter(0, 1.0);
  linearFunctionWithTilt->SetParameter(1, 0.0);
  linearFunctionWithTilt->SetLineColor(kRed);
  abcdBackgroundRatio->Draw("APE");
  
  abcdBackgroundRatio->Fit(linearFunction, "", "", 0, 7.0);
  linearFunction->Draw("sameL");
  
  abcdBackgroundRatio->Fit(linearFunctionWithTilt, "", "", 0, 7.0);
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
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    cout<<"Plotting "<<signalTitle[iSig]<<" ABCD"<<endl;
    TH2D *abcdPlot = GetABCDplot(metVsDedxHistsSignal.at((ESignal)iSig), bins, xtracks::kSignal, iSig);
    abcdCanvas->cd(iPad++);
    abcdPlot->SetMarkerSize(3.0);
    abcdPlot->Draw("colzText");
    outFile->cd();
    abcdPlot->SetName(signalName[iSig].c_str());
    abcdPlot->SetTitle(signalName[iSig].c_str());
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
  Log(1)<<"\n";
  
  for(auto groupMet : groupsMet){
    for(double met : groupMet) Log(1)<<met<<"\t";
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
    Log(1)<<"\n";
  }
  return make_tuple(bestMet, bestDedx);
}

/// Creates Combine datacard using Andrea's tool
void createDatacard(string outFileName, string outputPath)
{
  string command = "python ../DatacardCreatorABCD/mkDatacards.py  --inputHistoFile "+outputPath;
  command += " --dataHistoName  "+backgroundHistNams;
  command += " --sigHistoNameTemplate  Wino";
  
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
  
  command = "cp /afs/cern.ch/work/j/jniedzie/private/CMSSW_9_4_6_patch1/src/limits_datacard_"+outFileName+".txt";
  command += " macros/limitsData/combineOutput/";
  exec(command.c_str());
}

void convertRtoLimits(string outFileName)
{
  string command = "/Applications/root_v6.16.00/bin/root -q -b -l ";
  command += "\"macros/getLimitsFromR.C(\\\"macros/limitsData/combineOutput/limits_datacard_"+outFileName+".txt\\\", ";
  command += "\\\"macros/limitsData/cms_short_disappearing_"+outFileName+".txt\\\")\"";
  exec(command.c_str());
}

map<int, TH2D*> loadSignalHists(const EventSet &events)
{
  map<int, TH2D*> metVsDedxHistsSignal;
  for(int iSig=0; iSig<kNsignals; iSig++){
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
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    auto result = findBestBinning(metVsDedxHistBackground, metVsDedxHistsSignal.at(iSig), groupsDedx, groupsMet);
    double significance = GetSignificance(metVsDedxHistBackground, metVsDedxHistsSignal.at(iSig), result);

    Log(0)<<"Sample: "<<signalTitle[iSig]<<"\n";
    Log(0)<<"MET bins: "; for(double met : get<0>(result)) Log(0)<<met<<"\t"; Log(0)<<"\n";
    Log(0)<<"dE/dx bins: "; for(double dedx : get<1>(result)) Log(0)<<dedx<<"\t"; Log(0)<<"\n";
    Log(0)<<"significance: "<<significance<<"\n";
  }
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
  
  for(int iSig=0; iSig<kNsignals; iSig++){
     if(!config.runSignal[iSig]) continue;
     
     string outFileName = to_string_with_precision(nDedxBins, 0)+"x"+to_string_with_precision(nMetBins, 0)+"_"+config.category+"_"+sampleTag+"_"+signalShortName[iSig];
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
  srand((uint)time(0));
  
  cout.imbue(locale("de_DE"));
  TApplication *theApp = new TApplication("App", &argc, argv);
  config = ConfigManager("configs/analysis.md");
  
 
  
  // Load sll events with initial cuts only
  EventSet events;
  string prefix = "after_L"+to_string_with_precision(config.params["cuts_level"], 0)+"/"+config.category+"/";
  events.LoadEventsFromFiles(prefix);
  
  // Create histograms with number of events for each MET-dE/dx bin
  TH2D *metVsDedxHistBackground = GetMetVsDedxHist(events, xtracks::kBackground);
  map<int, TH2D*> metVsDedxHistsSignal = loadSignalHists(events);
  
  
//  runBinningScan(metVsDedxHistBackground, metVsDedxHistsSignal);
  
  produceLimits(metVsDedxHistBackground, metVsDedxHistsSignal);

//  drawAndSaveABCDplots(metVsDedxHistBackground, metVsDedxHistsSignal, {{400},{3.0, 3.1}} , "results/abcd_plots_debug.root");
  
   
  theApp->Run();
  return 0;
}




