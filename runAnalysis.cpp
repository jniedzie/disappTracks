#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

string configPath = "configs/analysis.md";

void saveEvents(const EventSet &events, const string &suffix)
{
  if(!config.params["save_events"]) return;
    
  string prefix = "after_L"+to_string((int)config.params["cuts_level"]);
  prefix = prefix + "/" + suffix + "/";
  if(config.params["cuts_level"]==10) prefix = "adish_cuts";
  events.SaveEventsToFiles(prefix);
  
}

void plotEvents(const EventSet &events)
{
  if(config.params["draw_standard_plots"])  events.DrawStandardPlots();
  if(config.params["draw_per_layer_plots"])  events.DrawPerLayerPlots();
}

void printDetails(const EventSet &events)
{
  if(config.params["print_background_details"]){
    for(EBackground iBck : backgrounds){
      if(!config.runBackground[iBck]) continue;
      for(int year : years){
        cout<<"Background events in "<<backgroundTitle.at(iBck)<<":"<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kBackground, iBck, year); iEvent++){
          events.At(xtracks::kBackground, iBck, year, iEvent)->Print();
        }
      }
    }
  }
  if(config.params["print_data_details"]){
    for(EData iData : datas){
      if(!config.runData[iData]) continue;
      for(int year : years){
        cout<<"Data events in "<<dataTitle[iData]<<":"<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kData, iData, year); iEvent++){
          events.At(xtracks::kData, iData, year, iEvent)->Print();
        }
      }
    }
  }
  if(config.params["print_signal_details"]){
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      for(int year : years){
        cout<<"Signal events in "<<signalTitle.at(iSig)<<":"<<endl;
        for(int iEvent=0;iEvent<events.size(xtracks::kSignal, iSig, year); iEvent++){
          events.At(xtracks::kSignal, iSig, year, iEvent)->Print();
        }
      }
    }
  }
}

void saveSurvivingEvents(const EventSet &events)
{
  bool updateData = false;
  
  for(EData iData : datas){
    if(config.runData[iData]){ updateData = true; break; }
  }
  
  string suffix = "L"+to_string_with_precision(config.params["cuts_level"],0);
  if(config.params["cuts_level"]==1) suffix+= "_"+config.category;
  
  if(updateData){
    ofstream dataSurvivingFile;
    dataSurvivingFile.open("results/survivingDataEventsAfter"+suffix+".txt");
    
    ofstream dataSurvivingFileByLS;
    dataSurvivingFileByLS.open("results/survivingDataEventsAfter"+suffix+"_byLS.txt");
    
    set<int> lumiSections;
    
    for(EData iData : datas){
      if(!config.runData[iData]) continue;
      for(int year : years){
        cout<<"Data events surviving cuts in "<<dataTitle[iData]<<":"<<events.size(xtracks::kData, iData, year)<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kData, iData, year); iEvent++){
          auto event = events.At(xtracks::kData, iData, year, iEvent);
          int runNumber = event->GetRunNumber();
          int lumiSection = event->GetLumiSection();
          long long int eventNumber = event->GetEventNumber();
          
          dataSurvivingFile<<runNumber<<":"<<lumiSection<<":"<<eventNumber<<"\n";
          
          if(lumiSections.find(lumiSection) == lumiSections.end()){
            dataSurvivingFileByLS<<runNumber<<" "<<lumiSection<<"\n";
            lumiSections.insert(lumiSection);
          }
        }
      }
    }
    dataSurvivingFile.close();
  }
  
  bool updateSignals = false;
  
  for(ESignal iSig : signals){
    if(config.runSignal[iSig]){ updateSignals = true; break; }
  }
  
  if(updateSignals){
    ofstream signalSurvivingFile;
    signalSurvivingFile.open ("results/survivingSignalEventsAfter"+suffix+".txt");
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      for(int year : years){
        cout<<"Signal events surviving cuts in "<<signalTitle.at(iSig)<<":"<<events.size(xtracks::kSignal, iSig, year)<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kSignal, iSig, year); iEvent++){
          auto event = events.At(xtracks::kSignal, iSig, year, iEvent);
          int runNumber = event->GetRunNumber();
          int lumiSection = event->GetLumiSection();
          long long int eventNumber = event->GetEventNumber();
          signalSurvivingFile<<runNumber<<":"<<lumiSection<<":"<<eventNumber<<"\n";
        }
      }
    }
    signalSurvivingFile.close();
  }
}

void processCuts(EventSet &events,
                 const EventCut &eventCut, const TrackCut &trackCut,
                 const JetCut &jetCut, const LeptonCut &leptonCut,
                 string suffix = "")
{
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(config.params["print_yields"]){
    cout<<"\n\nYields after level "<<config.params["cuts_level"]<<" cuts"<<endl;
    events.PrintYields();
  }
  saveEvents(events, suffix);
  plotEvents(events);
  printDetails(events);
  saveSurvivingEvents(events);
}

/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix;
  
  if(config.params["cuts_level"]==0 || config.params["cuts_level"]==10) prefix = "";
  if(config.params["cuts_level"]==1) prefix = "after_L0/";
  if(config.params["cuts_level"]==2) prefix = "after_L1/"+config.category+"/";
  if(config.params["cuts_level"]==20) prefix = "afterHelixTagging/";
  
  return prefix;
}

void runMETbinning(const EventSet &events,
                   EventCut &eventCut, const TrackCut &trackCut,
                   const JetCut &jetCut, const LeptonCut &leptonCut)
{
  cout<<"Combined significances from MET bins:"<<endl;
  
  int nMetBins = 10;
  
  vector<vector<double>> significances; // [metBin][iSig]
  
  for(int iMetBin=0; iMetBin<nMetBins; iMetBin++){
    
    double binMin = 200+iMetBin*100;
    double binMax = 200+(iMetBin+1)*100;
    eventCut.SetMetPt(range<double>(binMin, binMax));
    
    EventSet eventsForMetBin(events);
    processCuts(eventsForMetBin, eventCut, trackCut, jetCut, leptonCut);
    significances.push_back(eventsForMetBin.GetSignificance());
  }
  
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    
    double combinedSignificance = 0;
    
    for(int iMetBin=0; iMetBin<nMetBins; iMetBin++){
      if(significances[iMetBin][iSig] < 0 || !isnormal(significances[iMetBin][iSig])) continue;
      
      combinedSignificance += pow(significances[iMetBin][iSig],2);
    }
    combinedSignificance = sqrt(combinedSignificance);
    
    cout<<signalTitle.at(iSig)<<"\t"<<combinedSignificance<<endl;
  }
}

void scanMETbinning(const EventSet &events,
                    EventCut &eventCut, const TrackCut &trackCut,
                    const JetCut &jetCut, const LeptonCut &leptonCut)
{
  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("significance","significance",800,500);
  c1->cd();
  
  TH1D hists[signals.size()];
  
  
  double min = 201, max = 1001, step = 20;
  int nBins = (max-min)/step + 1;
  
  for(ESignal iSig : signals){
    hists[iSig] = TH1D(signalTitle.at(iSig).c_str(),signalTitle.at(iSig).c_str(),nBins,min,max);
  }
  
  for(int i=0;i<nBins;i++){
    double split = min+i*step;
    cout<<"\n\nSplit:"<<split<<"\n";
    EventSet events1(events);
    EventSet events2(events);
    
    eventCut.SetMetPt(range<double>(200,split));
    processCuts(events1, eventCut, trackCut, jetCut, leptonCut);
    
    eventCut.SetMetPt(range<double>(split,inf));
    processCuts(events2, eventCut, trackCut, jetCut, leptonCut);
    
    vector<double> significances1 = events1.GetSignificance();
    vector<double> significances2 = events2.GetSignificance();
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      
      double combinedSignificance = sqrt(pow(significances1[iSig],2) + pow(significances2[iSig],2));
      cout<<signalTitle.at(iSig)<<"\t"<<((isnormal(combinedSignificance) && combinedSignificance < 1000) ? combinedSignificance : 0.0)<<endl;
      hists[iSig].SetBinContent(i+1, (isnormal(combinedSignificance) && combinedSignificance < 1000) ? combinedSignificance : 0.0);
    }
  }
  
  bool first = true;
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    hists[iSig].SetMarkerSize(2.0);
    
    hists[iSig].SetLineColor(SignalColor((ESignal)iSig));
    hists[iSig].SetMarkerStyle(signalMarkers[iSig]);
    hists[iSig].SetMarkerColor(SignalColor((ESignal)iSig));
    
    
    double mean = 0;
    double minY = inf;
    
    for(int i=0;i<hists[iSig].GetNbinsX();i++){
      mean += hists[iSig].GetBinContent(i+1);
      if(hists[iSig].GetBinContent(i+1) < minY) minY = hists[iSig].GetBinContent(i+1);
    }
    mean /= nBins;
    
    for(int i=0;i<hists[iSig].GetNbinsX();i++){
      hists[iSig].SetBinContent(i+1,hists[iSig].GetBinContent(i+1)-minY);
    }
    hists[iSig].Scale(1/hists[iSig].Integral());
    hists[iSig].Sumw2(false);
    if(first){
      
      hists[iSig].Draw("PL");
      first = false;
    }
    else{
      hists[iSig].Draw("samePL");
    }
  }
  c1->Update();
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  
  string initPrefix = getPathPrefix();

  EventSet events;
  events.LoadEventsFromFiles(initPrefix);
  cout<<"\n\nInitial yields"<<endl;
  events.PrintYields();
  
  CutsManager cutsManager;
  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(config.params["cuts_level"] == 0){
    processCuts(events, eventCut, trackCut, jetCut, leptonCut);
  }
	
  if(config.params["cuts_level"] == 1){
    if(config.params["scan_MET_binning"]){
      scanMETbinning(events, eventCut, trackCut, jetCut, leptonCut);
    }
    else if(config.params["do_MET_binning"]){
      // for 2 tracks we don't have enough stats to do MET binning, just get a regular S/B
      if(config.category != "2-tracks") runMETbinning(events, eventCut, trackCut, jetCut, leptonCut);
      else                              processCuts(events, eventCut, trackCut, jetCut, leptonCut);
    }
    else{
      processCuts(events, eventCut, trackCut, jetCut, leptonCut, config.category);
    }
  }
  
  if(config.params["cuts_level"] == 2){
    processCuts(events, eventCut, trackCut, jetCut, leptonCut, config.category);
  }
  
  //---------------------------------------------------------------------------
  // Draw plots after helix tagging
  //---------------------------------------------------------------------------
  if(config.params["cuts_level"] == 20) events.DrawStandardPlots();
  
  cout<<"Done"<<endl;
  
  if(config.params["draw_standard_plots"]  ||
     config.params["draw_per_layer_plots"] ||
     config.params["scan_MET_binning"])
    theApp.Run();
  
  return 0;
}



