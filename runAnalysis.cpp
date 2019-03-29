#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"

#include "TGraph.h"

#include <TApplication.h>

string configPath = "configs/analysis.md";

void ProcessCuts(shared_ptr<EventSet> events,
                 const unique_ptr<EventCut> &eventCut,const  unique_ptr<TrackCut> &trackCut,
                 const unique_ptr<JetCut> &jetCut,const unique_ptr<LeptonCut> &leptonCut,
                 string suffix = "")
{
  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(config->printYields){
    cout<<"\n\nYields after level "<<config->performCutsLevel<<" cuts"<<endl;
    events->PrintYields();
  }
  if(config->saveEvents){
    string prefix = "after_L"+to_string(config->performCutsLevel);
    prefix = prefix + "/" + suffix + "/";
    if(config->performCutsLevel==10) prefix = "adish_cuts";
    events->SaveEventsToFiles(prefix);
  }
  if(config->drawStandardPlots){
    events->DrawStandardPlots();
  }
  if(config->drawPerLayerPlots){
    events->DrawPerLayerPlots();
  }
  if(config->printBackgroundDetails){
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!config->runBackground[iBck]) continue;
      cout<<"Background events in "<<backgroundTitle[iBck]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(xtracks::kBackground,iBck);iEvent++){
        events->At(xtracks::kBackground,iBck,iEvent)->Print();
      }
    }
  }
  if(config->printDataDetails){
    for(int iData=0;iData<kNdata;iData++){
      if(!config->runData[iData]) continue;
      cout<<"Data events in "<<dataTitle[iData]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(xtracks::kData,iData);iEvent++){
        events->At(xtracks::kData,iData,iEvent)->Print();
      }
    }
  }
  if(config->printSignalDetails){
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!config->runSignal[iSig]) continue;
      cout<<"Signal events in "<<signalTitle[iSig]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(xtracks::kSignal,iSig);iEvent++){
        events->At(xtracks::kSignal,iSig,iEvent)->Print();
      }
    }
  }
  
  ofstream dataSurvivingFile;
  dataSurvivingFile.open ("dataAfterL1.txt");
  
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
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
    if(!config->runSignal[iSig]) continue;
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
  config = make_unique<ConfigManager>(configPath);
  shared_ptr<EventSet> events = shared_ptr<EventSet>(new EventSet);
  
  string initPrefix = "after_L"+to_string(config->performCutsLevel-1)+"/";
  if(config->performCutsLevel==0 || config->performCutsLevel==10) initPrefix = "";
  if(config->performCutsLevel==20) initPrefix = "afterHelixTagging/";
  
  events->LoadEventsFromFiles(initPrefix);
  cout<<"\n\nInitial yields"<<endl;
  events->PrintYields();
  
  //---------------------------------------------------------------------------
  // Level 0
  //---------------------------------------------------------------------------

  
    auto eventCut_L0 = unique_ptr<EventCut>(new EventCut());
    auto trackCut_L0 = unique_ptr<TrackCut>(new TrackCut());
    auto jetCut_L0   = unique_ptr<JetCut>(new JetCut());
    auto leptonCut_L0= unique_ptr<LeptonCut>(new LeptonCut());
		
		// Remove bad jets
		jetCut_L0->SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
		jetCut_L0->SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
		jetCut_L0->SetPt(range<double>(30.0, inf));
		
		// Remove bad tracks
		trackCut_L0->SetNmissingInnerPixel(range<int>(0, 0));
		trackCut_L0->SetNmissingMiddleTracker(range<int>(0, 0));
		trackCut_L0->SetRelativeIsolation(range<double>(0.0, 0.5));
		trackCut_L0->SetNlayers(range<int>(2, inf));
		trackCut_L0->SetEta(range<double>(-2.1, 2.1));
		
		// Check MET properties
    eventCut_L0->SetMetNoMuPt(range<double>(200,inf));
    eventCut_L0->SetRequireMetNoMuTrigger(true);
    eventCut_L0->SetRequirePassingAllFilters(true);
		eventCut_L0->SetJetMetDeltaPhi(range<double>(0.5,inf));
		
		// Check leading jet properties
    eventCut_L0->SetLeadingJetPt(range<double>(100,inf));
    eventCut_L0->SetLeadingJetEta(range<double>(-2.4,2.4));
    eventCut_L0->SetLeadingJetNeHEF(range<double>(-inf,0.8));
    eventCut_L0->SetLeadingJetChHEF(range<double>(0.1,inf));
		
		// Check number of objects after cuts
		eventCut_L0->SetNtracks(range<int>(1,inf));
		eventCut_L0->SetNjets(range<int>(1,inf));
		eventCut_L0->SetNmuons(range<int>(0,0));
		eventCut_L0->SetNtaus(range<int>(0,0));
		eventCut_L0->SetNleptons(range<int>(0,0));
  
  if(config->performCutsLevel == 0){
    ProcessCuts(events,eventCut_L0, trackCut_L0, jetCut_L0, leptonCut_L0);
  }
	
  //---------------------------------------------------------------------------
  // Level 1
  //---------------------------------------------------------------------------
  if(config->performCutsLevel == 1){
    auto eventCut_L1 = make_unique<EventCut>(*eventCut_L0);
    auto trackCut_L1 = make_unique<TrackCut>(*trackCut_L0);
    auto jetCut_L1   = make_unique<JetCut>(*jetCut_L0);
    auto leptonCut_L1= make_unique<LeptonCut>(*leptonCut_L0);
    
//    trackCut_L1->SetRequireMcMatch(true);
    
    // pick category
    if(config->category == "2-tracks"){
      trackCut_L1->SetNmissingOuterTracker(range<int>(1, inf));
      trackCut_L1->SetCaloEmEnergy(range<double>(0.0,8.0));
      
      eventCut_L1->SetNtracks(range<int>(2,2));
    }
    else if(config->category == "3-layers"){
      trackCut_L1->SetNpixelLayers(range<int>(3, 3));
      trackCut_L1->SetRelativeIsolation(range<double>(0,0.11));
      trackCut_L1->SetNmissingOuterTracker(range<int>(3, inf));
      trackCut_L1->SetCaloEmEnergy(range<double>(0.0,3.0));
      trackCut_L1->SetDedxPerCluster(range<double>(2.0,inf));
      trackCut_L1->SetTrackMetDeltaPhi(range<double>(-2.3,2.3));
      
      eventCut_L1->SetNtracks(range<int>(1,1));
      eventCut_L1->SetJetMetDeltaPhi(range<double>(0.7,inf));
    }
    else if(config->category == "4-layers"){
      trackCut_L1->SetNpixelLayers(range<int>(4, 4));
      
//      trackCut_L1->SetNmissingOuterTracker(range<int>(7, inf));
      trackCut_L1->SetCaloEmEnergy(range<double>(0.0,0.4));
      trackCut_L1->SetDedxPerCluster(range<double>(2.0,inf));
//      trackCut_L1->SetPt(range<double>(100,inf));
      trackCut_L1->SetTrackMetDeltaPhi(range<double>(-2.3,2.3));
      
      
      eventCut_L1->SetNtracks(range<int>(1,1));
    }
    
    if(config->scanMETbinning){
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
        
        eventCut_L1->SetMetPt(range<double>(200,split));
        ProcessCuts(events1, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
        
        eventCut_L1->SetMetPt(range<double>(split,inf));
        ProcessCuts(events2, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
        
        vector<double> significances1 = events1->GetSignificance();
        vector<double> significances2 = events2->GetSignificance();
        
        for(int iSig=0;iSig<kNsignals;iSig++){
          if(!config->runSignal[iSig]) continue;
          
          double combinedSignificance = sqrt(pow(significances1[iSig],2) + pow(significances2[iSig],2));
          cout<<signalTitle[iSig]<<"\t"<<((isnormal(combinedSignificance) && combinedSignificance < 1000) ? combinedSignificance : 0.0)<<endl;
          hists[iSig]->SetBinContent(i+1, (isnormal(combinedSignificance) && combinedSignificance < 1000) ? combinedSignificance : 0.0);
        }
      }
      
      bool first = true;
      for(int iSig=0;iSig<kNsignals;iSig++){
        if(!config->runSignal[iSig]) continue;
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
    else if(config->doMETbinning){
      if(config->category != "2-tracks"){
        for(int iSig=0;iSig<kNsignals;iSig++){
          if(!config->runSignal[iSig]) continue;
          auto events1 = shared_ptr<EventSet>(new EventSet(*events));
          auto events2 = shared_ptr<EventSet>(new EventSet(*events));
          double split = -1;
          
          if(config->category == "3-layers"){
            if(iSig == kWino_M_300_cTau_3  || iSig == kWino_M_300_cTau_10 || iSig == kWino_M_300_cTau_30 ||
               iSig == kWino_M_500_cTau_10 || iSig == kWino_M_500_cTau_20){
              split = 380;
            }
            else if(iSig == kWino_M_650_cTau_10 || iSig == kWino_M_650_cTau_20 ||
                    iSig == kWino_M_800_cTau_20){
              split = 440;
            }
            else if(iSig == kWino_M_800_cTau_10  ||
                    iSig == kWino_M_1000_cTau_10 || iSig == kWino_M_1000_cTau_20){
              split = 580;
            }
          }
          else if(config->category == "4-layers"){
            if(iSig == kWino_M_300_cTau_3  || iSig == kWino_M_300_cTau_10 || iSig == kWino_M_300_cTau_30 ||
               iSig == kWino_M_500_cTau_10 || iSig == kWino_M_500_cTau_20 ||
               iSig == kWino_M_650_cTau_10 || iSig == kWino_M_650_cTau_20 ||
               iSig == kWino_M_800_cTau_20){
              split = 290;
            }
            else if(iSig == kWino_M_800_cTau_10  ||
                    iSig == kWino_M_1000_cTau_10 || iSig == kWino_M_1000_cTau_20){
              split = 440;
            }
          }
          
          eventCut_L1->SetMetPt(range<double>(200,split));
          ProcessCuts(events1, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
          
          eventCut_L1->SetMetPt(range<double>(split,inf));
          ProcessCuts(events2, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
          
          vector<double> significances1 = events1->GetSignificance();
          vector<double> significances2 = events2->GetSignificance();
          
          double combinedSignificance = sqrt(pow(significances1[iSig],2) + pow(significances2[iSig],2));
          cout<<signalTitle[iSig]<<"\t"<<combinedSignificance<<endl;
        }
        
        for(int iData=0;iData<kNdata;iData++){
          if(!config->runData[iData]) continue;
          
          vector<double> splits;
          
          if(config->category == "3-layers"){
            splits = { 380, 440, 580};
          }
          else if(config->category == "4-layers"){
            splits = { 290, 440 };
          }
          
          cout<<dataTitle[iData]<<"\n";
          
          for(double split : splits){
            auto events1 = shared_ptr<EventSet>(new EventSet(*events));
            auto events2 = shared_ptr<EventSet>(new EventSet(*events));
            
            eventCut_L1->SetMetPt(range<double>(200,split));
            ProcessCuts(events1, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
            
            eventCut_L1->SetMetPt(range<double>(split,inf));
            ProcessCuts(events2, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
            
            vector<double> significances1 = events1->GetSignificance(true);
            vector<double> significances2 = events2->GetSignificance(true);
            
            double combinedSignificance = sqrt(pow(significances1[iData],2) + pow(significances2[iData],2));
            cout<<"split:"<<split<<"\t"<<combinedSignificance<<endl;
          }
        }
      }
      else{
        // for 2 tracks we don't have enough stats to do MET binning, just get a regular S/B
        ProcessCuts(events, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
      }
    }
    else{
      string suffix = "";
      if(config->category == "2-tracks") suffix = "2tracks";
      if(config->category == "3-layers") suffix = "3layers";
      if(config->category == "4-layers") suffix = "4layers";
      ProcessCuts(events, eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1, suffix);
    }
  }
  
  //---------------------------------------------------------------------------
  // Adish cuts
  //---------------------------------------------------------------------------
  if(config->performCutsLevel == 10){
    auto eventCut_adish = unique_ptr<EventCut>(new EventCut());
    auto trackCut_adish = unique_ptr<TrackCut>(new TrackCut());
    auto jetCut_adish   = unique_ptr<JetCut>(new JetCut());
    auto leptonCut_adish= unique_ptr<LeptonCut>(new LeptonCut());
    
    // adish cuts
    eventCut_adish->SetRequireMetNoMuTrigger(true);
    eventCut_adish->SetMetNoMuPt(range<double>(200,inf));
    
    eventCut_adish->SetNjets(range<int>(1,inf));

    eventCut_adish->SetLeadingJetPt(range<double>(100,inf));
    eventCut_adish->SetLeadingJetEta(range<double>(-2.4,2.4));
    eventCut_adish->SetLeadingJetNeHEF(range<double>(-inf,0.8));
    eventCut_adish->SetLeadingJetChHEF(range<double>(0.1,inf));

    eventCut_adish->SetJetMetDeltaPhi(range<double>(0.5,inf));
    jetCut_adish->SetPt(range<double>(30, inf));
    
    eventCut_adish->SetNmuons(range<int>(0,0));
    eventCut_adish->SetNtaus(range<int>(0,0));
    eventCut_adish->SetNleptons(range<int>(0,0));
    
    ProcessCuts(events,eventCut_adish, trackCut_adish, jetCut_adish, leptonCut_adish);
  }
  
  //---------------------------------------------------------------------------
  // Draw plots after helix tagging
  //---------------------------------------------------------------------------
  if(config->performCutsLevel == 20){
    events->DrawStandardPlots();
  }
  
  if(config->drawStandardPlots || config->drawPerLayerPlots || config->scanMETbinning)  theApp->Run();
  return 0;
}



