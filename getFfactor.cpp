#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "LeptonCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

string basePath = "../ZnnStudy";

vector<string> ZmmFilePaths = {
  "Zmm/DYJetsM50_HT100to200",
  "Zmm/DYJetsM50_HT200to400",
  "Zmm/DYJetsM50_HT400to600",
  "Zmm/DYJetsM50_HT600to800",
  "Zmm/DYJetsM50_HT800to1200",
  "Zmm/DYJetsM50_HT1200to2500",
  "Zmm/DYJetsM50_HT2500toInf"
};

vector<string> ZvvFilePaths = {
  "Zvv/ZvvJets_HT100to200",
  "Zvv/ZvvJets_HT200to400",
  "Zvv/ZvvJets_HT400to600",
  "Zvv/ZvvJets_HT600to800",
  "Zvv/ZvvJets_HT800to1200",
  "Zvv/ZvvJets_HT1200to2500",
  "Zvv/ZvvJets_HT2500toInf"
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  vector<Events*> ZmmData;
  vector<Events*> ZvvData;
  
  const int nBins = 11;
  double bins[] = {200., 225., 250., 275., 300., 350., 400., 450., 500., 600., 800., 1000.};
  
  TH1D *nEventsZvv = new TH1D("N events Zvv","N events Zvv", nBins, bins);
  TH1D *nEventsZmm = new TH1D("N events Zmm","N events Zmm", nBins, bins);
  
  nEventsZvv->Sumw2();
  nEventsZmm->Sumw2();
  
  TCanvas *canvas = new TCanvas("F factor","F factor",2880,1800);
  canvas->Divide(2,2);
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  unsigned int eventCutOptionsZmm =
//  EventCut::kEmpty
    EventCut::kOneTrack
//    EventCut::kOneJet // ok
  | EventCut::kMetNoMu200GeV // ok
  | EventCut::kMetNoMuTrigger // ok
  | EventCut::kNoTau  // ok
  | EventCut::kMetNoMuJetPhi0p5 // ok
  | EventCut::kMuonsFromZ // ok
  | EventCut::kMuJetR0p4
  | EventCut::kMuTrackR0p4
  | EventCut::kTwoMuon
  | EventCut::kTightMuon
  | EventCut::kHighJetPt100GeV // ok
  | EventCut::kHighJetChHEF0p1 // ok
  | EventCut::kHighJetNeHEF0p8 // ok
  | EventCut::kHighJetEta2p4   // ok
  | EventCut::kHighJet
  ;
  
  unsigned int eventCutOptionsZvv =
//  EventCut::kEmpty
//    EventCut::kOneTrack
    EventCut::kOneJet // ok
  | EventCut::kMet200GeV // ok
  | EventCut::kMetNoMuTrigger // ok
  | EventCut::kNoTau // ok
  | EventCut::kMetJetPhi0p5 // ok
  | EventCut::kNoLepton // ok
  | EventCut::kHighJetPt100GeV // ok
  | EventCut::kHighJetChHEF0p1 // ok
  | EventCut::kHighJetNeHEF0p8 // ok
  | EventCut::kHighJetEta2p4 // ok
  | EventCut::kHighJet
  ;
  
  unsigned int trackCutOptions =
//    TrackCut::kEmpty
    TrackCut::kPt50GeV // ok
  | TrackCut::kEta2p4  // ok
  ;

  unsigned int jetCutOptions =
//    JetCut::kEmpty
    JetCut::kPt30GeV // ok
  | JetCut::kFwdEta4p7 // ok
  ;

  unsigned int leptonCutOptions =
    LeptonCut::kEmpty;
  
  
  EventCut  *eventCutZmm  = new EventCut((EventCut::ECut)eventCutOptionsZmm);
  EventCut  *eventCutZvv  = new EventCut((EventCut::ECut)eventCutOptionsZvv);

  TrackCut  *trackCut     = new TrackCut((TrackCut::ECut)trackCutOptions);
  JetCut    *jetCut       = new JetCut((JetCut::ECut)jetCutOptions);
  LeptonCut *leptonCut    = new LeptonCut((LeptonCut::ECut)leptonCutOptions);
  
  for(int iHT=0;iHT<ZmmFilePaths.size();iHT++){
    ZmmData.push_back(new Events(basePath+"/"+ZmmFilePaths[iHT]+"/tree.root",0));
    ZmmData[iHT] = ZmmData[iHT]->ApplyCuts(eventCutZmm, trackCut, jetCut, leptonCut);
    
    for(int iEvent=0;iEvent<ZmmData[iHT]->size();iEvent++){
      nEventsZmm->Fill(ZmmData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       ZmmData[iHT]->At(iEvent)->GetWeight());
    }
  }
    
  for(int iHT=0;iHT<ZvvFilePaths.size();iHT++){
    ZvvData.push_back(new Events(basePath+"/"+ZvvFilePaths[iHT]+"/tree.root",0));
    ZvvData[iHT] = ZvvData[iHT]->ApplyCuts(eventCutZvv, trackCut, jetCut, leptonCut);

    for(int iEvent=0;iEvent<ZvvData[iHT]->size();iEvent++){
      nEventsZvv->Fill(ZvvData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       ZvvData[iHT]->At(iEvent)->GetWeight());
    }
  }
  
  canvas->cd(1);
  nEventsZvv->DrawCopy();
  
  canvas->cd(2);
  nEventsZmm->DrawCopy();
  
  canvas->cd(3);
  nEventsZvv->Divide(nEventsZmm);
  nEventsZvv->Draw();
  
  cout<<"Finished"<<endl;
  
  canvas->Update();
  theApp.Run();
  return 0;
}



