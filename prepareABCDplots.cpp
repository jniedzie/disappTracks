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

string outputPath = "results/abcd_plots_3x3_3layers.root";

// Desired number of MET and dE/dx bins and limits of those
const int nMetBins  = 3, nDedxBins = 3;
const double minMet  = 300 , maxMet  = 500 , stepMet  = 50;
const double minDedx = 2.0 , maxDedx = 5.1 , stepDedx = 0.2;

/**
 Returns number of counts in ABCD... regions determined by criticalMet and criticalDedx values.
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param criticalMet Positions of bins border on MET axis
 \param criticalDedx Positions of bins border on dE/dx axis
 */
vector<vector<double>> GetABCD(const TH2D *metVsDedxHist,
                               const vector<double> &criticalMet,
                               const vector<double> &criticalDedx)
{
  vector<float> binsMet = { 0 };
  vector<float> binsDedx = { 0 };
  
  for(float met : criticalMet) binsMet.push_back(met);
  for(float dedx : criticalDedx) binsDedx.push_back(dedx);
  
  binsMet.push_back(inf);
  binsDedx.push_back(inf);
  
  vector<vector<double>> abcd(binsMet.size()-1, vector<double>(binsDedx.size()-1));
  
  for(int iMet=0; iMet<binsMet.size()-1; iMet++){
    for(int iDedx=0; iDedx<binsDedx.size()-1; iDedx++){
      
      int binX1 = metVsDedxHist->GetXaxis()->FindFixBin(binsDedx[iDedx]);
      int binX2 = metVsDedxHist->GetXaxis()->FindFixBin(binsDedx[iDedx+1])-1;
      
      int binY1 = metVsDedxHist->GetYaxis()->FindFixBin(binsMet[iMet]);
      int binY2 = metVsDedxHist->GetYaxis()->FindFixBin(binsMet[iMet+1])-1;
      
      abcd[iMet][iDedx] = metVsDedxHist->Integral(binX1, binX2, binY1, binY2);
    }
  }
  
  return abcd;
}

/**
 Creates ABCD... histogram. In case of backgrounds, all samples are merged into the same histogram.
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param criticalMet Positions of bins border on MET axis
 \param criticalDedx Positions of bins border on dE/dx axis
 \param dataType Specifies whether background or signal events should be analyzed
 \param setIter For signal, specified which samples to use
 */
TH2D* GetABCDplot(TH2D* metVsDedxHist, vector<double> criticalMet, vector<double> criticalDedx,
                  xtracks::EDataType dataType, int setIter=0)
{
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
 
  auto abcd = GetABCD(metVsDedxHist, criticalMet, criticalDedx);
  
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
 \param criticalMet Positions of bins border on MET axis
 \param criticalDedx Positions of bins border on dE/dx axis
 */
tuple<double, double> GetVariance(const TH2D *metVsDedxHist, const vector<double> criticalMet, const vector<double> criticalDedx)
{
  auto abcd = GetABCD(metVsDedxHist, criticalMet, criticalDedx);
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
                       const vector<double> &criticalMet, const vector<double> &criticalDedx)
{
  auto abcdBackground = GetABCD(metVsDedxHistBackground, criticalMet, criticalDedx);
  auto abcdSignal     = GetABCD(metVsDedxHistSignal, criticalMet, criticalDedx);
  
  double significance = 0;
  
  for(int x=0; x<abcdBackground.size(); x++){
    for(int y=0; y<abcdBackground[x].size(); y++){
      if(abcdSignal[x][y] + abcdBackground[x][y] == 0) continue;
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

/**
 Searches for criticalMet and criticalDedx for which relative difference between predicted and expected number
 of events in the signal bin is the smallers. Uses background events to find those values.
 \param metVsDedxHist Histogram containing number of events for each MET-dE/dx bin
 \param minMet Starting value of MET pt
 \param maxMet Maximum value of MET pt
 \param stepMet Step of MET pt
 \param minDedx Starting value of dE/dx
 \param maxDedx Maximum value of dE/dx
 \param stepDedx Step of dE/dx
 \return returns tuple containing: best critical MET, best critical dE/dx
 */
tuple<vector<double>,vector<double>> GetBestMetAndDedx(TH2D *metVsDedxHist,
                                                       int nMetBins, double minMet, double maxMet, double stepMet,
                                                       int nDedxBins, double minDedx, double maxDedx, double stepDedx,
                                                       vector<double> initialMet, vector<double> initialDedx)
{
  double bestVariance = inf;
  vector<double> bestMet, bestDedx;
  
  if(initialMet.size()==0){
    double initialMetBinWidth = (maxMet-minMet)/nMetBins;
    for(int metBin=0; metBin<nMetBins-1; metBin++) bestMet.push_back(minMet+metBin*initialMetBinWidth);
  }
  else bestMet = initialMet;
  
  if(initialDedx.size()==0){
    double initialDedxBinWidth = (maxDedx-minDedx)/nDedxBins;
    for(int dedxBin=0; dedxBin<nDedxBins-1; dedxBin++)  bestDedx.push_back(minDedx+dedxBin*initialDedxBinWidth);
  }
  else bestDedx = initialDedx;
  
  cout<<"\nInitial bins:"<<endl;
  cout<<"MET bins: "; for(double met : bestMet) cout<<met<<"\t"; cout<<endl;
  cout<<"dE/dx bins: "; for(double dedx : bestDedx) cout<<dedx<<"\t"; cout<<"\n"<<endl;
  
  vector<double> foundMet, foundDedx;
  
  for(int metBin=0; metBin<nMetBins-1; metBin++){
    
    double _minMet = metBin == 0 ? minMet : bestMet[metBin-1];
    double _maxMet = metBin == nMetBins-2 ? maxMet : bestMet[metBin+1];
    
    for(int dedxBin=0; dedxBin<nDedxBins-1; dedxBin++){
      Log(1)<<"Testing "<<metBin<<" MET and "<<dedxBin<<" dE/dx bin border\n";
      
      double _minDedx = dedxBin == 0 ? minDedx : bestDedx[dedxBin-1];
      double _maxDedx = dedxBin == nDedxBins-2 ? maxDedx : bestDedx[dedxBin+1];
      
      for(bestMet[metBin] = _minMet; bestMet[metBin] <= _maxMet; bestMet[metBin] += stepMet){
        for(bestDedx[dedxBin] = _minDedx; bestDedx[dedxBin] <= _maxDedx; bestDedx[dedxBin] += stepDedx){
          Log(2)<<"\tcurrent MET: "<<bestMet[metBin]<<"\tdE/dx: "<<bestDedx[dedxBin];
          auto [variance, varianceError] = GetVariance(metVsDedxHist, bestMet, bestDedx);
          Log(2)<<"\tvariance:"<<variance<<" +/- "<<varianceError<<"\n";
          
          if(variance < bestVariance){
            bestVariance = variance;
            foundMet = bestMet;
            foundDedx = bestDedx;
          }
        }
      }
    }
  }
  return make_tuple(foundMet, foundDedx);
}

tuple<vector<double>,vector<double>> GetBestMetAndDedxRandomly(TH2D *metVsDedxHist,
                                                               int nMetBins, double minMet, double maxMet, double stepMet,
                                                               int nDedxBins, double minDedx, double maxDedx, double stepDedx,
                                                               int nTrials)
{
  double bestVariance = inf;
  vector<double> bestMet, bestDedx;
  vector<double> foundMet, foundDedx;
  
  for(int metBin=0; metBin<nMetBins-1; metBin++) bestMet.push_back(0);
  for(int dedxBin=0; dedxBin<nDedxBins-1; dedxBin++)  bestDedx.push_back(0);
  
  double currentMet, currentDedx;
  
  for(int iTrial=0; iTrial<nTrials; iTrial++){
    if(iTrial%100==0) Log(1)<<"Trial "<<iTrial<<"\n";
  
    for(int metBin=0; metBin<nMetBins-1; metBin++){
      do{ currentMet = (int)(RandDouble(minMet, maxMet)/stepMet)*stepMet;}
      while(find(bestMet.begin(), bestMet.end(), currentMet) != bestMet.end());
      bestMet[metBin] = currentMet;
    }
    for(int dedxBin=0; dedxBin<nDedxBins-1; dedxBin++){
      do{ currentDedx = (int)(RandDouble(minDedx, maxDedx)/stepDedx)*stepDedx;}
      while(find(bestDedx.begin(), bestDedx.end(), currentDedx) != bestDedx.end());
      bestDedx[dedxBin] = currentDedx;
    }
    sort(bestMet.begin(), bestMet.end());
    sort(bestDedx.begin(), bestDedx.end());
    
    auto [variance, varianceError] = GetVariance(metVsDedxHist, bestMet, bestDedx);
    
    if(variance < bestVariance){
      bestVariance = variance;
      foundMet = bestMet;
      foundDedx = bestDedx;
    }
  }
  
  return make_tuple(foundMet, foundDedx);
}

TGraphErrors* GetRatioGraph(const TH2D *metVsDedxHistBackground,
                            const vector<double> &criticalMet, const vector<double> &criticalDedx)
{
  TH2D *metVsDedxRebinned = new TH2D(*metVsDedxHistBackground);
  metVsDedxRebinned->RebinX(2);
  
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


/// Draws and saves ABCD plots for given values of critical MET and critical dE/dx
void DrawAndSaveABCDplots(TH2D *metVsDedxHistBackground,
                          map<int, TH2D*> metVsDedxHistsSignal,
                          vector<double> criticalMet, vector<double> criticalDedx)
{
  TCanvas *abcdCanvas = new TCanvas("ABCD", "ABCD", 1000, 1500);
  abcdCanvas->Divide(4,3);
  TFile *outFile = new TFile(outputPath.c_str(),"recreate");
  
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".1f");
  
  cout<<"Plotting background ABCD"<<endl;
  TH2D *abcdPlotBackgrounds = GetABCDplot(metVsDedxHistBackground, criticalMet, criticalDedx, xtracks::kBackground);
  abcdPlotBackgrounds->SetMarkerSize(3.0);
  
  TCanvas *ratioCanvas = new TCanvas("Ratio", "Ratio", 800, 600);
  ratioCanvas->cd();
  
  TGraphErrors *abcdBackgroundRatio = GetRatioGraph(metVsDedxHistBackground, criticalMet, criticalDedx);
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
  abcdPlotBackgrounds->SetName("background");
  abcdPlotBackgrounds->SetTitle("background");
  abcdPlotBackgrounds->Write();
  abcdBackgroundRatio->Write();
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    cout<<"Plotting "<<signalTitle[iSig]<<" ABCD"<<endl;
    TH2D *abcdPlot = GetABCDplot(metVsDedxHistsSignal[iSig], criticalMet, criticalDedx, xtracks::kSignal, iSig);
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
      double significance = GetSignificance(metVsDedxHistBackground, metVsDedxHistSignal, groupMet, groupDedx);
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
  map<int, TH2D*> metVsDedxHistsSignal;
  
  // Find all combinations of binning in dE/dx and MET
  vector<vector<double>> groupsDedx;
  for(double startingDedx=minDedx; startingDedx<maxDedx; startingDedx+=stepDedx) groupsDedx.push_back({startingDedx});
  AddValuesCombinations(groupsDedx, maxDedx, stepDedx, nDedxBins-1);

  vector<vector<double>> groupsMet;
  for(double startingMet=minMet; startingMet<maxMet; startingMet+=stepMet) groupsMet.push_back({startingMet});
  AddValuesCombinations(groupsMet, maxMet, stepMet, nMetBins-1);
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    metVsDedxHistsSignal[iSig] = GetMetVsDedxHist(events, xtracks::kSignal, iSig);
    
//    auto result = findBestBinning(metVsDedxHistBackground, metVsDedxHistsSignal[iSig],
//                                  groupsDedx, groupsMet);
//
//    vector<double> bestMet = get<0>(result);
//    vector<double> bestDedx = get<1>(result);
//
//    double significance = GetSignificance(metVsDedxHistBackground,
//                                          metVsDedxHistsSignal[iSig],
//                                          bestMet, bestDedx);
//
//    Log(0)<<"Sample: "<<signalTitle[iSig]<<"\n";
//    Log(0)<<"MET bins: "; for(double met : bestMet) Log(0)<<met<<"\t"; Log(0)<<"\n";
//    Log(0)<<"dE/dx bins: "; for(double dedx : bestDedx) Log(0)<<dedx<<"\t"; Log(0)<<"\n";
//    Log(0)<<"significance: "<<significance<<"\n";
  }
  
  vector<double> bestMet={300, 450};
  vector<double> bestDedx={3.0, 4.5};

//  vector<double> bestMet={500};
//  vector<double> bestDedx={2.3};

  DrawAndSaveABCDplots(metVsDedxHistBackground, metVsDedxHistsSignal, bestMet, bestDedx);

  theApp->Run();
  return 0;
}




