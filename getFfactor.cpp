#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "LeptonCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

string basePath = "../SR_MC";

vector<string> ZmmFilePaths = {
//  "DYJetsM50_HT100to200",
//  "DYJetsM50_HT200to400",
//  "DYJetsM50_HT400to600",
//  "DYJetsM50_HT600to800",
//  "DYJetsM50_HT800to1200",
//  "DYJetsM50_HT1200to2500",
//  "DYJetsM50_HT2500toInf"
};

vector<string> WvlFilePaths = {
  "WJets_HT100to200",
  "WJets_HT200to400",
  "WJets_HT400to600",
  "WJets_HT600to800",
  "WJets_HT800to1200",
  "WJets_HT1200to2500",
  "WJets_HT2500toInf"
};

vector<string> ZvvFilePaths = {
  "ZvvJets_HT100to200",
  "ZvvJets_HT200to400",
  "ZvvJets_HT400to600",
  "ZvvJets_HT600to800",
  "ZvvJets_HT800to1200",
  "ZvvJets_HT1200to2500",
  "ZvvJets_HT2500toInf"
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  vector<Events*> ZmmData;
  vector<Events*> ZvvData;
  vector<Events*> WvlData;
  
  const int nBins = 11;
  double bins[] = {200., 225., 250., 275., 300., 350., 400., 450., 500., 600., 800., 1000.};
  
  TH1D *nEventsZvv = new TH1D("N events Zvv","N events Zvv", nBins, bins);
  TH1D *nEventsZmm = new TH1D("N events Zmm","N events Zmm", nBins, bins);
  TH1D *nEventsWvl = new TH1D("N events Wvl","N events Wvl", nBins, bins);
  
  nEventsZvv->Sumw2();
  nEventsZmm->Sumw2();
  nEventsWvl->Sumw2();
  
  TCanvas *canvas = new TCanvas("F factor","F factor",2880,1800);
  canvas->Divide(2,3);
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  unsigned int eventCutOptionsZmm =
    EventCut::kOneTrack
  | EventCut::kMetNoMuTrigger
  | EventCut::kNoTau
  | EventCut::kHighJetPt100GeV
  | EventCut::kHighJetChHEF0p1
  | EventCut::kHighJetNeHEF0p8
  | EventCut::kHighJetEta2p4
  | EventCut::kHighJet
  | EventCut::kMetNoMu200GeV
  | EventCut::kMetNoMuJetPhi0p5
  | EventCut::kMuonsFromZ
  | EventCut::kMuJetR0p4
  | EventCut::kMuTrackR0p4
  | EventCut::kTwoOppositeMuons
  | EventCut::kTightMuon
  ;
  
  unsigned int eventCutOptionsWvl =
    EventCut::kOneTrack
  | EventCut::kMetNoMuTrigger
  | EventCut::kNoTau
  | EventCut::kHighJetPt100GeV
  | EventCut::kHighJetChHEF0p1
  | EventCut::kHighJetNeHEF0p8
  | EventCut::kHighJetEta2p4
  | EventCut::kHighJet
  | EventCut::kMetNoMu200GeV
  | EventCut::kMetNoMuJetPhi0p5
//  | EventCut::kMuonsFromZ
  | EventCut::kMuJetR0p4
  | EventCut::kMuTrackR0p4
//  | EventCut::kTwoOppositeMuons
  | EventCut::kTightMuon
  ;
  
  unsigned int eventCutOptionsZvv =
    EventCut::kOneTrack
  | EventCut::kMetNoMuTrigger
  | EventCut::kNoTau
  | EventCut::kHighJetPt100GeV
  | EventCut::kHighJetChHEF0p1
  | EventCut::kHighJetNeHEF0p8
  | EventCut::kHighJetEta2p4
  | EventCut::kHighJet
  | EventCut::kMet200GeV
  | EventCut::kMetJetPhi0p5
  | EventCut::kNoLepton
  ;
  
  unsigned int trackCutOptions =
    TrackCut::kPt50GeV
  | TrackCut::kEta2p4
//  | TrackCut::kShort
//  | TrackCut::kLowCalo
//  | TrackCut::kLowDEdx
  ;

  unsigned int jetCutOptions =
    JetCut::kPt30GeV
  | JetCut::kFwdEta4p7
  ;

  unsigned int leptonCutOptions =
    LeptonCut::kEmpty;
  
  
  EventCut  *eventCutZmm  = new EventCut((EventCut::ECut)eventCutOptionsZmm);
  EventCut  *eventCutZvv  = new EventCut((EventCut::ECut)eventCutOptionsZvv);
  EventCut  *eventCutWvl  = new EventCut((EventCut::ECut)eventCutOptionsWvl);

  TrackCut  *trackCut     = new TrackCut((TrackCut::ECut)trackCutOptions);
  JetCut    *jetCut       = new JetCut((JetCut::ECut)jetCutOptions);
  LeptonCut *leptonCut    = new LeptonCut((LeptonCut::ECut)leptonCutOptions);
  
  for(int iHT=0;iHT<ZmmFilePaths.size();iHT++){
    ZmmData.push_back(new Events(basePath+"/"+ZmmFilePaths[iHT]+"/tree.root",Events::kBackground));
    ZmmData[iHT] = ZmmData[iHT]->ApplyCuts(eventCutZmm, trackCut, jetCut, leptonCut);
    
    for(int iEvent=0;iEvent<ZmmData[iHT]->size();iEvent++){
      nEventsZmm->Fill(ZmmData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       ZmmData[iHT]->At(iEvent)->GetWeight());
    }
  }
  
  for(int iHT=0;iHT<WvlFilePaths.size();iHT++){
    WvlData.push_back(new Events(basePath+"/"+WvlFilePaths[iHT]+"/tree.root",Events::kBackground));
    WvlData[iHT] = WvlData[iHT]->ApplyCuts(eventCutWvl, trackCut, jetCut, leptonCut);
    
    for(int iEvent=0;iEvent<WvlData[iHT]->size();iEvent++){
      nEventsWvl->Fill(WvlData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       WvlData[iHT]->At(iEvent)->GetWeight());
    }
  }
    
  for(int iHT=0;iHT<ZvvFilePaths.size();iHT++){
    ZvvData.push_back(new Events(basePath+"/"+ZvvFilePaths[iHT]+"/tree.root",Events::kBackground));
    ZvvData[iHT] = ZvvData[iHT]->ApplyCuts(eventCutZvv, trackCut, jetCut, leptonCut);

    for(int iEvent=0;iEvent<ZvvData[iHT]->size();iEvent++){
      nEventsZvv->Fill(ZvvData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       ZvvData[iHT]->At(iEvent)->GetWeight());
    }
  }
  
  TFile *outFile = new TFile("fFactor.root","recreate");
  outFile->cd();
  
  canvas->cd(1);
  nEventsZvv->DrawCopy();
  nEventsZvv->Write("Zvv");
  
  canvas->cd(2);
  nEventsZmm->DrawCopy();
  nEventsZmm->Write("Zmm");
  
  canvas->cd(3);
  nEventsWvl->DrawCopy();
  nEventsWvl->Write("Wvl");
  
  canvas->cd(4);
  TH1D *fFactorZmm = new TH1D(*nEventsZvv);
  fFactorZmm->Divide(nEventsZmm);
  fFactorZmm->Draw();
  fFactorZmm->Write("f factor (ZvvZmm)");
  
  canvas->cd(4);
  TH1D *fFactorWvl = new TH1D(*nEventsZvv);
  fFactorWvl->Divide(nEventsWvl);
  fFactorWvl->Draw();
  fFactorWvl->Write("f factor (ZvvWvl)");
  
  cout<<"Finished"<<endl;
  
  canvas->Update();
  outFile->Close();
  
  theApp.Run();
  return 0;
}



