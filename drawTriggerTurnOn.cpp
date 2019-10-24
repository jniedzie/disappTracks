#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

#include <TApplication.h>

string configPath = "configs/analysis.md";

int main(int argc, char* argv[])
{
  // All events with initial cuts only
  config = ConfigManager(configPath);
  EventSet events;
  events.LoadEventsFromFiles();
  
  float bins[25] = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 240, 260, 300, 350, 400, 500};
  
  TH1D *triggerHistBackgroundsNum[backgrounds.size()];
  TH1D *triggerHistBackgroundsDen[backgrounds.size()];
  
  TH1D *triggerHistSignalsNum[signals.size()];
  TH1D *triggerHistSignalsDen[signals.size()];
  
  string title = "trigger_turnon_backgrounds";
  TH1D *triggerHistBackgroundsAllNum = new TH1D(title.c_str(), title.c_str(), 24, bins);
  triggerHistBackgroundsAllNum->Sumw2();
  TH1D *triggerHistBackgroundsAllDen = new TH1D(*triggerHistBackgroundsAllNum);
  
  title = "trigger_turnon_signals";
  TH1D *triggerHistSignalsAllNum = new TH1D(title.c_str(), title.c_str(), 24, bins);
  triggerHistSignalsAllNum->Sumw2();
  TH1D *triggerHistSignalsAllDen = new TH1D(*triggerHistSignalsAllNum);
  
  TFile *outFile = new TFile("results/triggerTurnOn.root", "recreate");
  outFile->cd();
  
  int nBackground = 0;
  int nBackgroundPassing = 0;
  
  int nSignal = 0;
  int nSignalPassing = 0;
  
  for(auto iBck : backgrounds){
    if(!config.runBackground[iBck]) continue;
    
    string title = "trigger_turnon_"+backgroundName.at(iBck);
    triggerHistBackgroundsNum[iBck] = new TH1D(title.c_str(), title.c_str(), 24, bins);
    triggerHistBackgroundsNum[iBck]->Sumw2();
    triggerHistBackgroundsDen[iBck] = new TH1D(*triggerHistBackgroundsNum[iBck]);
  
  
    for(int iEvent=0; iEvent<events.size(kBackground, iBck, 2017); iEvent++){
      auto event = events.At(kBackground, iBck, 2017, iEvent);
      
      if(event->HasMetNoMuTrigger()){
        triggerHistBackgroundsNum[iBck]->Fill(event->GetMetNoMuPt());
        triggerHistBackgroundsAllNum->Fill(event->GetMetNoMuPt());
        
        nBackground++;
        if(event->GetMetNoMuPt() > 200) nBackgroundPassing++;
        
      }
      triggerHistBackgroundsDen[iBck]->Fill(event->GetMetNoMuPt());
      triggerHistBackgroundsAllDen->Fill(event->GetMetNoMuPt());
    }
    
    triggerHistBackgroundsNum[iBck]->Divide(triggerHistBackgroundsNum[iBck], triggerHistBackgroundsDen[iBck], 1, 1, "B");
    triggerHistBackgroundsNum[iBck]->Write();
  }
  
  triggerHistBackgroundsAllNum->Divide(triggerHistBackgroundsAllNum, triggerHistBackgroundsAllDen, 1, 1, "B");
  triggerHistBackgroundsAllNum->Write();
  
  for(auto iSig : signals){
    if(!config.runSignal[iSig]) continue;
    
    string title = "trigger_turnon_"+signalName.at(iSig);
    triggerHistSignalsNum[iSig] = new TH1D(title.c_str(), title.c_str(), 24, bins);
    triggerHistSignalsNum[iSig]->Sumw2();
    triggerHistSignalsDen[iSig] = new TH1D(*triggerHistSignalsNum[iSig]);
  
  
    for(int iEvent=0; iEvent<events.size(kSignal, iSig, 2017); iEvent++){
      auto event = events.At(kSignal, iSig, 2017, iEvent);
      
      if(event->HasMetNoMuTrigger()){
        triggerHistSignalsNum[iSig]->Fill(event->GetMetNoMuPt());
        triggerHistSignalsAllNum->Fill(event->GetMetNoMuPt());
        
        nSignal++;
        if(event->GetMetNoMuPt() > 200) nSignalPassing++;
      }
      triggerHistSignalsDen[iSig]->Fill(event->GetMetNoMuPt());
      triggerHistSignalsAllDen->Fill(event->GetMetNoMuPt());
    }
    
    triggerHistSignalsNum[iSig]->Divide(triggerHistSignalsNum[iSig], triggerHistSignalsDen[iSig], 1, 1, "B");
    triggerHistSignalsNum[iSig]->Write();
  }
  
  
  triggerHistSignalsAllNum->Divide(triggerHistSignalsAllNum, triggerHistSignalsAllDen, 1, 1, "B");
  triggerHistSignalsAllNum->Write();
  
  outFile->Close();
  
  cout<<"N background: "<<nBackground<<"\tpassing: "<<nBackgroundPassing<<endl;
  cout<<"Fraction background passing: "<<nBackgroundPassing/(double)nBackground<<endl;
  cout<<"N signal: "<<nSignal<<"\tpassing: "<<nSignalPassing<<endl;
  cout<<"Fraction signal passing: "<<nSignalPassing/(double)nSignal<<endl;
  
  return 0;
}
