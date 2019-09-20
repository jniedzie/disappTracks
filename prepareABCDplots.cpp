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
  
  vector<vector<double>> abcd;
  for(int iMet=0; iMet<binsMet.size()-1; iMet++){
    vector<double> tmp;
    for(int iDedx=0; iDedx<binsDedx.size()-1; iDedx++){
      tmp.push_back(0.0);
    }
    abcd.push_back(tmp);
  }
  
  for(int iMet=0; iMet<binsMet.size()-1; iMet++){
    for(int iDedx=0; iDedx<binsDedx.size()-1; iDedx++){
      
      int binX1 = metVsDedxHist->GetXaxis()->FindFixBin(binsDedx[iDedx]);
      int binX2 = metVsDedxHist->GetXaxis()->FindFixBin(binsDedx[iDedx+1]);
      
      int binY1 = metVsDedxHist->GetYaxis()->FindFixBin(binsMet[iMet]);
      int binY2 = metVsDedxHist->GetYaxis()->FindFixBin(binsMet[iMet+1]);
      
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
double GetVariance(const TH2D *metVsDedxHist, const vector<double> criticalMet, const vector<double> criticalDedx)
{
  auto abcd = GetABCD(metVsDedxHist, criticalMet, criticalDedx);
  double variance = 0;
  int nEntries = 0;
  for(int x=0; x<abcd.size()-1; x++){
    for(int y=0; y<abcd[x].size()-1; y++){
      if(abcd[x][y+1]==0  || abcd[x][y]==0    ||
         abcd[x+1][y]==0  || abcd[x+1][y+1]==0) return inf;
      
      double ratio = abcd[x][y+1]/abcd[x][y] * abcd[x+1][y]/abcd[x+1][y+1];
      variance += pow(ratio-1, 2);
      nEntries++;
    }
  }
  return variance/nEntries;
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
          hist->Fill(track->GetAverageDedx(), event->GetMetNoMuPt());
        }
      }
    }
  }
  else if(dataType == xtracks::kSignal){
    for(int iEvent=0;iEvent<events.size(dataType, setIter);iEvent++){
      auto event = events.At(dataType, setIter, iEvent);
      
      for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
        auto track = event->GetTrack(iTrack);
        hist->Fill(track->GetAverageDedx(), event->GetMetNoMuPt());
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
          double variance = GetVariance(metVsDedxHist, bestMet, bestDedx);
          Log(2)<<"\tvariance:"<<variance<<"\n";
          
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
    
    double variance = GetVariance(metVsDedxHist, bestMet, bestDedx);
    
    if(variance < bestVariance){
      bestVariance = variance;
      foundMet = bestMet;
      foundDedx = bestDedx;
    }
  }
  
  return make_tuple(foundMet, foundDedx);
}


/// Draws and saves ABCD plots for given values of critical MET and critical dE/dx
void DrawAndSaveABCDplots(TH2D *metVsDedxHistBackground,
                          map<int, TH2D*> metVsDedxHistsSignal,
                          vector<double> criticalMet, vector<double> criticalDedx)
{
  TCanvas *abcdCanvas = new TCanvas("ABCD", "ABCD", 1000, 1500);
  abcdCanvas->Divide(4,3);
  
  TFile *outFile = new TFile("results/abcd_plots.root","recreate");
  
  gStyle->SetOptStat(0);
  
  cout<<"Plotting background ABCD"<<endl;
  TH2D *abcdPlotBackgrounds = GetABCDplot(metVsDedxHistBackground, criticalMet, criticalDedx, xtracks::kBackground);
  int iPad=1;
  abcdCanvas->cd(iPad++);
  abcdPlotBackgrounds->SetMarkerSize(3.0);
  abcdPlotBackgrounds->Draw("colzText");
  outFile->cd();
  abcdPlotBackgrounds->Write();
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    cout<<"Plotting "<<signalTitle[iSig]<<" ABCD"<<endl;
    TH2D *abcdPlot = GetABCDplot(metVsDedxHistsSignal[iSig], criticalMet, criticalDedx, xtracks::kSignal, iSig);
    abcdCanvas->cd(iPad++);
    abcdPlot->SetMarkerSize(3.0);
    abcdPlot->Draw("colzText");
    outFile->cd();
    abcdPlot->Write();
  }
  
  abcdCanvas->Update();
  abcdCanvas->Write();
  outFile->Close();
}

/// Starting point of the application
int main(int argc, char* argv[])
{
  srand((uint)time(0));
  
  cout.imbue(locale("de_DE"));
  TApplication *theApp = new TApplication("App", &argc, argv);
  config = ConfigManager("configs/analysis.md");
  
  // All events with initial cuts only
  EventSet events;
  string prefix = "after_L"+to_string_with_precision(config.params["cuts_level"], 0)+"/"+config.category+"/";
  events.LoadEventsFromFiles(prefix);
  
  // Find the best values of critical MET and critical dE/dx
  int nMetBins  = 3, nDedxBins = 3;
  
  double minMet  = 200 , maxMet  = 300 , stepMet  = 10;
  double minDedx = 3.0 , maxDedx = 4.5 , stepDedx = 0.1;
  
  vector<double> bestMet, bestDedx;
  
  TH2D *metVsDedxHistBackground = GetMetVsDedxHist(events, xtracks::kBackground);
  
  map<int, TH2D*> metVsDedxHistsSignal;
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    metVsDedxHistsSignal[iSig] = GetMetVsDedxHist(events, xtracks::kSignal, iSig);
  }
  
  double bestVariance = inf;
  
//  for(int i=0; i<10; i++){
//    Log(0)<<"\n\nStarting "<<i+1<<" iteration\n";
//    auto result = GetBestMetAndDedx(metVsDedxHistBackground,
//                                    nMetBins, minMet, maxMet, stepMet,
//                                    nDedxBins, minDedx, maxDedx, stepDedx,
//                                    bestMet, bestDedx);
//
//    vector<double> foundMet = get<0>(result);
//    vector<double> foundDedx = get<1>(result);
//
//    double variance = GetVariance(metVsDedxHistBackground, foundMet, foundDedx);
//    Log(0)<<"Variance: "<<variance<<"\n";
//
//    if(variance < bestVariance){
//      bestVariance = variance;
//      bestMet = foundMet;
//      bestDedx = foundDedx;
//    }
//    else{
//      break;
//    }
//  }
  
//  auto result = GetBestMetAndDedxRandomly(metVsDedxHistBackground,
//                                          nMetBins, minMet, maxMet, stepMet,
//                                          nDedxBins, minDedx, maxDedx, stepDedx,
//                                          1000);
//
//  bestMet = get<0>(result);
//  bestDedx = get<1>(result);
  
//  Log(0)<<"MET bins: "; for(double met : bestMet) Log(0)<<met<<"\t"; Log(0)<<"\n";
//  Log(0)<<"dE/dx bins: "; for(double dedx : bestDedx) Log(0)<<dedx<<"\t"; Log(0)<<"\n";
//  Log(0)<<"variance: "<<GetVariance(metVsDedxHistBackground, bestMet, bestDedx)<<"\n";
  bestMet={200, 250, 300};
  bestDedx={4.2, 4.6};
  
  // Print restuls
 
  DrawAndSaveABCDplots(metVsDedxHistBackground, metVsDedxHistsSignal, bestMet, bestDedx);
  
  theApp->Run();
  return 0;
}




