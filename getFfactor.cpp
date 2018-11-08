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
  
  vector<shared_ptr<Events>> ZmmData;
  vector<shared_ptr<Events>> ZvvData;
  vector<shared_ptr<Events>> WvlData;
  
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
  
  EventCut  *eventCutZmm  = new EventCut();
  EventCut  *eventCutWvl  = new EventCut();
  EventCut  *eventCutZvv  = new EventCut();
  
  
  eventCutZmm->SetNtracks(range<int>(1,inf));
  eventCutZmm->SetRequireMetNoMuTrigger(true);
  eventCutZmm->SetNtaus(range<int>(0,0));
  eventCutZmm->SetLeadingJetPt(range<double>(100,inf));
  eventCutZmm->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutZmm->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutZmm->SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutZmm->SetRequireHighJet(true);
  eventCutZmm->SetMetNoMuPt(range<double>(200,inf));
  eventCutZmm->SetRequireMetNoMuJetPhi0p5(true);
  eventCutZmm->SetRequireMuonsFromZ(true);
  eventCutZmm->SetRequireMuJetR0p4(true);
  eventCutZmm->SetRequireMuTrackR0p4(true);
  eventCutZmm->SetRequireTwoOppositeMuons(true);
  eventCutZmm->SetRequireTightMuon(true);
  
  eventCutWvl->SetNtracks(range<int>(1,inf));
  eventCutWvl->SetRequireMetNoMuTrigger(true);
  eventCutWvl->SetNtaus(range<int>(0,0));
  eventCutWvl->SetLeadingJetPt(range<double>(100,inf));
  eventCutWvl->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutWvl->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutWvl->SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutWvl->SetRequireHighJet(true);
  eventCutWvl->SetMetNoMuPt(range<double>(200,inf));
  eventCutWvl->SetRequireMetNoMuJetPhi0p5(true);
  eventCutWvl->SetRequireMuJetR0p4(true);
  eventCutWvl->SetRequireMuTrackR0p4(true);
  eventCutWvl->SetRequireTightMuon(true);
  
 
  eventCutZvv->SetNtracks(range<int>(1,inf));
  eventCutZvv->SetRequireMetNoMuTrigger(true);
  eventCutZvv->SetNtaus(range<int>(0,0));
  eventCutZvv->SetLeadingJetPt(range<double>(100,inf));
  eventCutZvv->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutZvv->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutZvv->SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutZvv->SetRequireHighJet(true);
  eventCutZvv->SetMetPt(range<double>(200,inf));
  eventCutZvv->SetRequireMetJetPhi0p5(true);
  eventCutZvv->SetNleptons(range<int>(0,0));
  
  TrackCut *trackCut = new TrackCut();
  trackCut->SetPt(range<double>(50,inf));
  trackCut->SetEta(range<double>(-2.1, 2.1));

  JetCut *jetCut = new JetCut();
  jetCut->SetPt(range<double>(30,inf));
  jetCut->SetEtaForward(range<double>(-4.7, 4.7));
  
  LeptonCut *leptonCut = new LeptonCut();
  
  for(int iHT=0;iHT<(int)ZmmFilePaths.size();iHT++){
    ZmmData.push_back(make_shared<Events>(basePath+"/"+ZmmFilePaths[iHT]+"/tree.root",Events::kBackground));
    ZmmData[iHT] = ZmmData[iHT]->ApplyCuts(eventCutZmm, trackCut, jetCut, leptonCut);
    
    for(int iEvent=0;iEvent<ZmmData[iHT]->size();iEvent++){
      nEventsZmm->Fill(ZmmData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       ZmmData[iHT]->At(iEvent)->GetWeight());
    }
  }
  
  for(int iHT=0;iHT<(int)WvlFilePaths.size();iHT++){
    WvlData.push_back(make_shared<Events>(basePath+"/"+WvlFilePaths[iHT]+"/tree.root",Events::kBackground));
    WvlData[iHT] = WvlData[iHT]->ApplyCuts(eventCutWvl, trackCut, jetCut, leptonCut);
    
    for(int iEvent=0;iEvent<WvlData[iHT]->size();iEvent++){
      nEventsWvl->Fill(WvlData[iHT]->At(iEvent)->GetMetNoMuPt(),
                       WvlData[iHT]->At(iEvent)->GetWeight());
    }
  }
    
  for(int iHT=0;iHT<(int)ZvvFilePaths.size();iHT++){
    ZvvData.push_back(make_shared<Events>(basePath+"/"+ZvvFilePaths[iHT]+"/tree.root",Events::kBackground));
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



