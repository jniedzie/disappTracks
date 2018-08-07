#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

const int nHTbins = 4;

string basePath = "../ZnnStudy";

string ZmmFilePaths[nHTbins] = {
//  "Zmm/DYJetsM50_HT100to200",
//  "Zmm/DYJetsM50_HT200to400",
//  "Zmm/DYJetsM50_HT400to600",
  "Zmm/DYJetsM50_HT600to800",
  "Zmm/DYJetsM50_HT800to1200",
  "Zmm/DYJetsM50_HT1200to2500",
  "Zmm/DYJetsM50_HT2500toInf"
};

string ZvvFilePaths[nHTbins] = {
//  "Zvv/ZvvJets_HT100to200",
//  "Zvv/ZvvJets_HT200to400",
//  "Zvv/ZvvJets_HT400to600",
  "Zvv/ZvvJets_HT600to800",
  "Zvv/ZvvJets_HT800to1200",
  "Zvv/ZvvJets_HT1200to2500",
  "Zvv/ZvvJets_HT2500toInf"
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  Events *ZmmData[nHTbins];
  Events *ZvvData[nHTbins];
  
  TH1D *nEventsZvv[nHTbins];
  TH1D *nEventsZmm[nHTbins];
  
  TH1D *nEventsZvvHTsum = nullptr;
  TH1D *nEventsZmmHTsum = nullptr;
  TH1D *fFactorSum = nullptr;
  
  TCanvas *canvas = new TCanvas("F factor","F factor",2880,1800);
  canvas->Divide(3,3);
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  unsigned int eventCutOptions =
//  EventCut::kEmpty;
  EventCut::kOneTrack
  | EventCut::kOneJet
  | EventCut::kMetNoMu100GeV
  | EventCut::kMetNoMuTrigger
  | EventCut::kNoLepton
  | EventCut::kNoTau;
  
  unsigned int trackCutOptions =
//  TrackCut::kEmpty;
  TrackCut::kPt50GeV
  | TrackCut::kEta2p4;
//  | TrackCut::kLowDEdx;
  //  | TrackCut::kHighPt;
  //  | TrackCut::kMedium;
  
  unsigned int jetCutOptions =
    JetCut::kEmpty;
//  JetCut::kPt100GeV;
  
  EventCut  *eventCut = new EventCut((EventCut::ECut)eventCutOptions);
  TrackCut  *trackCut = new TrackCut((TrackCut::ECut)trackCutOptions);
  JetCut    *jetCut   = new JetCut((JetCut::ECut)jetCutOptions);
  
  for(int iHT=0;iHT<nHTbins;iHT++){
  
    ZmmData[iHT] = new Events(basePath+"/"+ZmmFilePaths[iHT]+"/tree.root",0);
    ZvvData[iHT] = new Events(basePath+"/"+ZvvFilePaths[iHT]+"/tree.root",0);
  
    ZmmData[iHT] = ZmmData[iHT]->ApplyCuts(eventCut, trackCut, jetCut);
    ZvvData[iHT] = ZvvData[iHT]->ApplyCuts(eventCut, trackCut, jetCut);
  
    nEventsZvv[iHT] = new TH1D(Form("N events Zvv (HT bin %i)",iHT),Form("N events Zvv (HT bin %i)",iHT), 10, 0, 1000);
    nEventsZmm[iHT] = new TH1D(Form("N events Zmm (HT bin %i)",iHT),Form("N events Zmm (HT bin %i)",iHT), 10, 0, 1000);
    
    nEventsZvv[iHT]->Sumw2();
    nEventsZmm[iHT]->Sumw2();
    
    for(int iEvent=0;iEvent<ZvvData[iHT]->size();iEvent++){
      nEventsZvv[iHT]->Fill(ZvvData[iHT]->At(iEvent)->GetMetNoMuPt(), ZvvData[iHT]->At(iEvent)->GetWeight());
    }
    
    for(int iEvent=0;iEvent<ZmmData[iHT]->size();iEvent++){
      nEventsZmm[iHT]->Fill(ZmmData[iHT]->At(iEvent)->GetMetNoMuPt(), ZmmData[iHT]->At(iEvent)->GetWeight());
    }
    
    if(iHT==0){
      nEventsZvvHTsum = new TH1D(*nEventsZvv[iHT]);
      nEventsZmmHTsum = new TH1D(*nEventsZmm[iHT]);
      
      nEventsZvvHTsum->SetTitle("N events Zvv (sum over all HT bins)");
      nEventsZmmHTsum->SetTitle("N events Zmm (sum over all HT bins)");
    }
    else{
      nEventsZvvHTsum->Add(nEventsZvv[iHT]);
      nEventsZmmHTsum->Add(nEventsZmm[iHT]);
    }
 
  
    canvas->cd(1);
    nEventsZvv[iHT]->DrawCopy(iHT==0 ? "" : "same");
    
    canvas->cd(2);
    nEventsZmm[iHT]->DrawCopy(iHT==0 ? "" : "same");
    
    canvas->cd(3);
    nEventsZvv[iHT]->Divide(nEventsZmm[iHT]);
    
    if(iHT==0){
      fFactorSum = new TH1D(*nEventsZvv[iHT]);
      fFactorSum->SetTitle("f factor (sum over all HT bins)");
    }
    else{
      fFactorSum->Add(nEventsZvv[iHT]);
    }
    
    nEventsZvv[iHT]->Draw(iHT==0 ? "" : "same");
  }
  
  canvas->cd(5);
  nEventsZmmHTsum->DrawCopy();
  
  canvas->cd(6);
  nEventsZvvHTsum->DrawCopy();
  
  canvas->cd(4);
  nEventsZvvHTsum->Divide(nEventsZmmHTsum);
  nEventsZvvHTsum->SetTitle("f factor (sum of num / sum of den)");
  nEventsZvvHTsum->Draw();
  
  canvas->cd(7);
  fFactorSum->Draw();
  
  cout<<"Finished"<<endl;
  
  canvas->Update();
  theApp.Run();
  return 0;
}



