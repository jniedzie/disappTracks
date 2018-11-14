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
  eventsZmm->LoadEventsFromFiles(EventSet::kBackground,kZmumuJets);
  
  shared_ptr<EventSet> eventsWmv = shared_ptr<EventSet>(new EventSet);
  eventsWmv->LoadEventsFromFiles(EventSet::kBackground,kWmunuJets);
  
  shared_ptr<EventSet> eventsZvv = shared_ptr<EventSet>(new EventSet);
  eventsZvv->LoadEventsFromFiles(EventSet::kBackground,kZnunuJets);
  
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
  
  auto eventCutZmm = unique_ptr<EventCut>(new EventCut());
  auto eventCutWmv = unique_ptr<EventCut>(new EventCut());
  auto eventCutZvv = unique_ptr<EventCut>(new EventCut());
  
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
  
  eventCutWmv->SetNtracks(range<int>(1,inf));
  eventCutWmv->SetRequireMetNoMuTrigger(true);
  eventCutWmv->SetNtaus(range<int>(0,0));
  eventCutWmv->SetLeadingJetPt(range<double>(100,inf));
  eventCutWmv->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCutWmv->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCutWmv->SetLeadingJetChHEF(range<double>(0.1,inf));
  eventCutWmv->SetRequireHighJet(true);
  eventCutWmv->SetMetNoMuPt(range<double>(200,inf));
  eventCutWmv->SetRequireMetNoMuJetPhi0p5(true);
  eventCutWmv->SetRequireMuJetR0p4(true);
  eventCutWmv->SetRequireMuTrackR0p4(true);
  eventCutWmv->SetRequireTightMuon(true);
  
 
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
  
  auto trackCut = unique_ptr<TrackCut>(new TrackCut());
  trackCut->SetPt(range<double>(50,inf));
  trackCut->SetEta(range<double>(-2.1, 2.1));

  auto jetCut = unique_ptr<JetCut>(new JetCut());
  jetCut->SetPt(range<double>(30,inf));
  jetCut->SetEtaForward(range<double>(-4.7, 4.7));
  
  auto leptonCut = unique_ptr<LeptonCut>(new LeptonCut());
  
  eventsZmm->ApplyCuts(eventCutZmm, trackCut, jetCut, leptonCut);
  eventsWmv->ApplyCuts(eventCutWmv, trackCut, jetCut, leptonCut);
  eventsZvv->ApplyCuts(eventCutZvv, trackCut, jetCut, leptonCut);
  
  for(int iEvent=0;iEvent<eventsZmm->size(EventSet::kBackground,kZmumuJets);iEvent++){
    shared_ptr<Event> event = eventsZmm->At(EventSet::kBackground,kZmumuJets,iEvent);
    nEventsZmm->Fill(event->GetMetNoMuPt(), event->GetWeight());
  }
  
  for(int iEvent=0;iEvent<eventsWmv->size(EventSet::kBackground,kWmunuJets);iEvent++){
    shared_ptr<Event> event = eventsWmv->At(EventSet::kBackground,kWmunuJets,iEvent);
    nEventsWmv->Fill(event->GetMetNoMuPt(), event->GetWeight());
  }
  
  for(int iEvent=0;iEvent<eventsZvv->size(EventSet::kBackground,kZnunuJets);iEvent++){
    shared_ptr<Event> event = eventsZvv->At(EventSet::kBackground,kZnunuJets,iEvent);
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



