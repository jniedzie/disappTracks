#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

#include "TGraph.h"

#include <TApplication.h>

string configPath = "configs/analysis.md";

void ProcessCuts(shared_ptr<EventSet> events,
                 const EventCut   &eventCut,
                 const TrackCut   &trackCut,
                 const JetCut     &jetCut,
                 const LeptonCut  &leptonCut,
                 string suffix = "")
{
  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(config.printYields){
    cout<<"\n\nYields after level "<<config.performCutsLevel<<" cuts"<<endl;
    events->PrintYields();
  }
  if(config.saveEvents){
    string prefix = "after_L"+to_string(config.performCutsLevel);
    prefix = prefix + "/" + suffix + "/";
    if(config.performCutsLevel==10) prefix = "adish_cuts";
    events->SaveEventsToFiles(prefix);
  }
  if(config.drawStandardPlots){
    events->DrawStandardPlots();
  }
  if(config.drawPerLayerPlots){
    events->DrawPerLayerPlots();
  }
  if(config.printBackgroundDetails){
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!config.runBackground[iBck]) continue;
      cout<<"Background events in "<<backgroundTitle[iBck]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(xtracks::kBackground,iBck);iEvent++){
        events->At(xtracks::kBackground,iBck,iEvent)->Print();
      }
    }
  }
  if(config.printDataDetails){
    for(int iData=0;iData<kNdata;iData++){
      if(!config.runData[iData]) continue;
      cout<<"Data events in "<<dataTitle[iData]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(xtracks::kData,iData);iEvent++){
        events->At(xtracks::kData,iData,iEvent)->Print();
      }
    }
  }
  if(config.printSignalDetails){
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!config.runSignal[iSig]) continue;
      cout<<"Signal events in "<<signalTitle[iSig]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(xtracks::kSignal,iSig);iEvent++){
        events->At(xtracks::kSignal,iSig,iEvent)->Print();
      }
    }
  }
  
  ofstream dataSurvivingFile;
  dataSurvivingFile.open ("dataAfterL1.txt");
  
  for(int iData=0;iData<kNdata;iData++){
    if(!config.runData[iData]) continue;
    cout<<"Data events surviving cuts in "<<dataTitle[iData]<<":"<<events->size(xtracks::kData,iData)<<endl;
    for(int iEvent=0;iEvent<events->size(xtracks::kData,iData);iEvent++){
      int runNumber = events->At(xtracks::kData,iData,iEvent)->GetRunNumber();
      int lumiSection = events->At(xtracks::kData,iData,iEvent)->GetLumiSection();
      long long int eventNumber = events->At(xtracks::kData,iData,iEvent)->GetEventNumber();
      dataSurvivingFile<<runNumber<<":"<<lumiSection<<":"<<eventNumber<<"\n";
    }
  }
  dataSurvivingFile.close();
  
  ofstream signalSurvivingFile;
  signalSurvivingFile.open ("signalAfterL1.txt");
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config.runSignal[iSig]) continue;
    cout<<"Signal events surviving cuts in "<<signalTitle[iSig]<<":"<<events->size(xtracks::kSignal,iSig)<<endl;
    for(int iEvent=0;iEvent<events->size(xtracks::kSignal,iSig);iEvent++){
      int runNumber = events->At(xtracks::kSignal,iSig,iEvent)->GetRunNumber();
      int lumiSection = events->At(xtracks::kSignal,iSig,iEvent)->GetLumiSection();
      long long int eventNumber = events->At(xtracks::kSignal,iSig,iEvent)->GetEventNumber();
      signalSurvivingFile<<runNumber<<":"<<lumiSection<<":"<<eventNumber<<"\n";
    }
  }
  signalSurvivingFile.close();
  
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  config = ConfigManager(configPath);
  shared_ptr<EventSet> events = shared_ptr<EventSet>(new EventSet);
  
  string initPrefix = "after_L"+to_string(config.performCutsLevel-1)+"/";
  if(config.performCutsLevel==0 || config.performCutsLevel==10) initPrefix = "";
  if(config.performCutsLevel==20) initPrefix = "afterHelixTagging/";
  
  events->LoadEventsFromFiles(initPrefix);
  cout<<"\n\nInitial yields"<<endl;
  events->PrintYields();
  
  CutsManager cutsManager;
  
  //---------------------------------------------------------------------------
  // Level 0
  //---------------------------------------------------------------------------

  
  EventCut eventCut;
  TrackCut trackCut;
  JetCut jetCut;
  LeptonCut leptonCut;
  
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  
  if(config.performCutsLevel == 0){
    ProcessCuts(events, eventCut, trackCut, jetCut, leptonCut);
  }
	
  if(config.performCutsLevel == 1){
    if(config.scanMETbinning){
      gStyle->SetOptStat(0);
      TCanvas *c1 = new TCanvas("significance","significance",800,500);
      c1->cd();
      
      TH1D *hists[kNsignals];
      
      
      double min = 201, max = 1001, step = 20;
      int nBins = (max-min)/step + 1;
      
      for(int iSig=0;iSig<kNsignals;iSig++){
        hists[iSig] = new TH1D(signalTitle[iSig].c_str(),signalTitle[iSig].c_str(),nBins,min,max);
      }
      
      for(int i=0;i<nBins;i++){
        double split = min+i*step;
        cout<<"\n\nSplit:"<<split<<"\n";
        auto events1 = shared_ptr<EventSet>(new EventSet(*events));
        auto events2 = shared_ptr<EventSet>(new EventSet(*events));
        
        eventCut.SetMetPt(range<double>(200,split));
        ProcessCuts(events1, eventCut, trackCut, jetCut, leptonCut);
        
        eventCut.SetMetPt(range<double>(split,inf));
        ProcessCuts(events2, eventCut, trackCut, jetCut, leptonCut);
        
        vector<double> significances1 = events1->GetSignificance();
        vector<double> significances2 = events2->GetSignificance();
        
        for(int iSig=0;iSig<kNsignals;iSig++){
          if(!config.runSignal[iSig]) continue;
          
          double combinedSignificance = sqrt(pow(significances1[iSig],2) + pow(significances2[iSig],2));
          cout<<signalTitle[iSig]<<"\t"<<((isnormal(combinedSignificance) && combinedSignificance < 1000) ? combinedSignificance : 0.0)<<endl;
          hists[iSig]->SetBinContent(i+1, (isnormal(combinedSignificance) && combinedSignificance < 1000) ? combinedSignificance : 0.0);
        }
      }
      
      bool first = true;
      for(int iSig=0;iSig<kNsignals;iSig++){
        if(!config.runSignal[iSig]) continue;
        hists[iSig]->SetMarkerSize(2.0);
        
        hists[iSig]->SetLineColor(SignalColor((ESignal)iSig));
        hists[iSig]->SetMarkerStyle(signalMarkers[iSig]);
        hists[iSig]->SetMarkerColor(SignalColor((ESignal)iSig));
        
        
        double mean = 0;
        double minY = inf;
        
        for(int i=0;i<hists[iSig]->GetNbinsX();i++){
          mean += hists[iSig]->GetBinContent(i+1);
          if(hists[iSig]->GetBinContent(i+1) < minY) minY = hists[iSig]->GetBinContent(i+1);
        }
        mean /= nBins;
        
        for(int i=0;i<hists[iSig]->GetNbinsX();i++){
          hists[iSig]->SetBinContent(i+1,hists[iSig]->GetBinContent(i+1)-minY);
        }
        hists[iSig]->Scale(1/hists[iSig]->Integral());
        hists[iSig]->Sumw2(false);
        if(first){
          
          hists[iSig]->Draw("PL");
          first = false;
        }
        else{
          hists[iSig]->Draw("samePL");
        }
      }
      c1->Update();
    }
    else if(config.doMETbinning){
      if(config.category != "2-tracks"){
        
        cout<<"Combined significances from MET bins:"<<endl;
        
        
        int nMetBins = 10;
        
        vector<vector<double>> significances; // [metBin][iSig]
        
        for(int iMetBin=0; iMetBin<nMetBins; iMetBin++){
          
          double binMin = 200+iMetBin*100;
          double binMax = 200+(iMetBin+1)*100;
          eventCut.SetMetPt(range<double>(binMin, binMax));
          
          auto eventsForMetBin = make_shared<EventSet>(*events);
          ProcessCuts(eventsForMetBin, eventCut, trackCut, jetCut, leptonCut);
          significances.push_back(eventsForMetBin->GetSignificance());
        }

        for(int iSig=0;iSig<kNsignals;iSig++){
          if(!config.runSignal[iSig]) continue;
          
          double combinedSignificance = 0;
          
          for(int iMetBin=0; iMetBin<nMetBins; iMetBin++){
            if(significances[iMetBin][iSig] < 0 || !isnormal(significances[iMetBin][iSig])) continue;
            
            combinedSignificance += pow(significances[iMetBin][iSig],2);
          }
          combinedSignificance = sqrt(combinedSignificance);
          
          cout<<signalTitle[iSig]<<"\t"<<combinedSignificance<<endl;
        }
      }
      else{
        // for 2 tracks we don't have enough stats to do MET binning, just get a regular S/B
        ProcessCuts(events, eventCut, trackCut, jetCut, leptonCut);
      }
    }
    else{
      string suffix = "";
      if(config.category == "2-tracks") suffix = "2tracks";
      if(config.category == "3-layers") suffix = "3layers";
      if(config.category == "4-layers") suffix = "4layers";
      ProcessCuts(events, eventCut, trackCut, jetCut, leptonCut, suffix);
    }
  }
  
  //---------------------------------------------------------------------------
  // Draw plots after helix tagging
  //---------------------------------------------------------------------------
  if(config.performCutsLevel == 20){
    events->DrawStandardPlots();
  }
  
  if(config.drawStandardPlots || config.drawPerLayerPlots || config.scanMETbinning)  theApp->Run();
  return 0;
}



