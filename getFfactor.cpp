#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "LeptonCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  // Load W->μν, Z->νν and Z->μμ events
  shared_ptr<EventSet> eventsZmm = shared_ptr<EventSet>(new EventSet);
  eventsZmm->LoadEventsFromFiles(xtracks::kBackground,kZmumuJets);
  
  shared_ptr<EventSet> eventsWmv = shared_ptr<EventSet>(new EventSet);
  eventsWmv->LoadEventsFromFiles(xtracks::kBackground,kWmunuJets);
  
  shared_ptr<EventSet> eventsZvv = shared_ptr<EventSet>(new EventSet);
  eventsZvv->LoadEventsFromFiles(xtracks::kBackground,kZnunuJets);
  
  // Prepare histograms
  const int nBins = 11;
  double bins[] = {200., 225., 250., 275., 300., 350., 400., 450., 500., 600., 800., 1000.};
  
  TH1D *nEventsZvv = new TH1D("N events Zvv","N events Zvv", nBins, bins);
  TH1D *nEventsZmm = new TH1D("N events Zmm","N events Zmm", nBins, bins);
  TH1D *nEventsWmv = new TH1D("N events Wvl","N events Wvl", nBins, bins);
  
  nEventsZvv->Sumw2();
  nEventsZmm->Sumw2();
  nEventsWmv->Sumw2();
  
  TCanvas *canvas = new TCanvas("F factor","F factor",2880,1800);
  canvas->Divide(2,3);
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  
  auto eventCutZmm = EventCut();
  auto eventCutWmv = EventCut();
  auto eventCutZvv = EventCut();
  
  eventCutZmm.SetNtracks(range<int>(1,inf));
  eventCutZmm.SetRequireMetNoMuTrigger(true);
  eventCutZmm.SetNtaus(range<int>(0,0));
  eventCutZmm.SetLeadingJetPt(range<double>(100,inf));
  eventCutZmm.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutZmm.SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutZmm.SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutZmm.SetMetNoMuPt(range<double>(200,inf));
  eventCutZmm.SetJetMetDeltaPhi(range<double>(0.5,inf));
  eventCutZmm.SetRequireMuonsFromZ(true);
  eventCutZmm.SetJetMuonDeltaPhi(range<double>(0.4, inf));
  eventCutZmm.SetTrackMuonDeltaPhi(range<double>(0.4, inf));
  eventCutZmm.SetRequireTwoOppositeMuons(true);
  eventCutZmm.SetRequireTightMuon(true);
  
  eventCutWmv.SetNtracks(range<int>(1,inf));
  eventCutWmv.SetRequireMetNoMuTrigger(true);
  eventCutWmv.SetNtaus(range<int>(0,0));
  eventCutWmv.SetLeadingJetPt(range<double>(100,inf));
  eventCutWmv.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutWmv.SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutWmv.SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutWmv.SetMetNoMuPt(range<double>(200,inf));
  eventCutWmv.SetJetMetDeltaPhi(range<double>(0.5,inf));
  eventCutWmv.SetJetMuonDeltaPhi(range<double>(0.4, inf));
  eventCutWmv.SetTrackMuonDeltaPhi(range<double>(0.4, inf));
  eventCutWmv.SetRequireTightMuon(true);
  
  eventCutZvv.SetNtracks(range<int>(1,inf));
  eventCutZvv.SetRequireMetNoMuTrigger(true);
  eventCutZvv.SetNtaus(range<int>(0,0));
  eventCutZvv.SetLeadingJetPt(range<double>(100,inf));
  eventCutZvv.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutZvv.SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutZvv.SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutZvv.SetMetPt(range<double>(200,inf));
  eventCutZvv.SetJetMetDeltaPhi(range<double>(0.5,inf));
  eventCutZvv.SetNleptons(range<int>(0,0));
  
  TrackCut trackCut;
  trackCut.SetPt(range<double>(50,inf));
  trackCut.SetEta(range<double>(-2.1, 2.1));

  JetCut jetCut;
  jetCut.SetPt(range<double>(30,inf));
  jetCut.SetEtaForward(range<double>(-4.7, 4.7));
  
  LeptonCut leptonCut;
  
  eventsZmm->ApplyCuts(eventCutZmm, trackCut, jetCut, leptonCut);
  eventsWmv->ApplyCuts(eventCutWmv, trackCut, jetCut, leptonCut);
  eventsZvv->ApplyCuts(eventCutZvv, trackCut, jetCut, leptonCut);
  
  for(int iEvent=0;iEvent<eventsZmm->size(xtracks::kBackground,kZmumuJets);iEvent++){
    shared_ptr<Event> event = eventsZmm->At(xtracks::kBackground,kZmumuJets,iEvent);
    nEventsZmm->Fill(event->GetMetNoMuPt(), event->GetWeight());
  }
  
  for(int iEvent=0;iEvent<eventsWmv->size(xtracks::kBackground,kWmunuJets);iEvent++){
    shared_ptr<Event> event = eventsWmv->At(xtracks::kBackground,kWmunuJets,iEvent);
    nEventsWmv->Fill(event->GetMetNoMuPt(), event->GetWeight());
  }
  
  for(int iEvent=0;iEvent<eventsZvv->size(xtracks::kBackground,kZnunuJets);iEvent++){
    shared_ptr<Event> event = eventsZvv->At(xtracks::kBackground,kZnunuJets,iEvent);
    nEventsZvv->Fill(event->GetMetNoMuPt(), event->GetWeight());
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
  nEventsWmv->DrawCopy();
  nEventsWmv->Write("Wvl");
  
  canvas->cd(4);
  TH1D *fFactorZmm = new TH1D(*nEventsZvv);
  fFactorZmm->Divide(nEventsZmm);
  fFactorZmm->Draw();
  fFactorZmm->Write("f factor (ZvvZmm)");
  
  canvas->cd(4);
  TH1D *fFactorWvl = new TH1D(*nEventsZvv);
  fFactorWvl->Divide(nEventsWmv);
  fFactorWvl->Draw();
  fFactorWvl->Write("f factor (ZvvWvl)");
  
  cout<<"Finished"<<endl;
  
  canvas->Update();
  outFile->Close();
  
  theApp.Run();
  return 0;
}



