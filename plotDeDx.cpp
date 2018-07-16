#include "Helpers.hpp"
#include "Event.hpp"
#include "TrackCut.hpp"

#include <TApplication.h>

const bool drawPlots = true;
const bool analyzeData = false;

TrackCut* GetShortTrackCut();
TrackCut* GetShortAboveThresholdTrackCut();

vector<TH1D*> GetDedxPerLayerHists(Events *events)
{
  vector<TH1D*> hists;
  static int iter=0;
  
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    TH1D *dedxHist = new TH1D(Form("dedx_layer[%i]_%i",iLayer,iter),Form("dedx_layer[%i]_%i",iLayer,iter),50,0,10);
    hists.push_back(dedxHist);
  }
  
  for(int iEvent=0;iEvent<events->size();iEvent++){
    for(int iTrack=0;iTrack<events->At(iEvent)->GetNtracks();iTrack++){
      Track *track = events->At(iEvent)->GetTrack(iTrack);
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0001){
          hists[iLayer]->Fill(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  iter++;
  return hists;
}


int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  string inFileNameSignal = "../jniedzie/mcSignal/tree.root";
  string inFileNameBackground = "../jniedzie/mcBackground/tree.root";
//  const char *inFileNameSignal = "../adish/Signal/tree.root";
//  const char *inFileNameBackground = "../adish/Background/tree.root";
  string inFileNameData = "../adish/Data/tree.root";
  
  
  TrackCut *shortTrackCut = GetShortTrackCut();
  TrackCut *shortAboveTrasholdTrackCut = GetShortAboveThresholdTrackCut();
  
  Events *eventsSignal = new Events(inFileNameSignal);
  int nTracksSignal = eventsSignal->GetNtracks();
  int nShortTracksSignal = eventsSignal->ApplyTrackCut(shortTrackCut)->GetNtracks();
  int nShortTracksAboveThresholdSignal = eventsSignal->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  
  Events *eventsBackground = new Events(inFileNameBackground);
  int nTracksBackground = eventsBackground->GetNtracks();
  int nShortTracksBackground = eventsBackground->ApplyTrackCut(shortTrackCut)->GetNtracks();
  int nShortTracksAboveThresholdBackground = eventsBackground->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  
  int nTracksData, nShortTracksData, nShortTracksAboveThresholdData;
  if(analyzeData){
    Events *eventsData = new Events(inFileNameData);
    nTracksData = eventsData->GetNtracks();
    nShortTracksData = eventsData->ApplyTrackCut(shortTrackCut)->GetNtracks();
    nShortTracksAboveThresholdData = eventsData->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  }
  
  cout<<"Total number of tracks (signal):"<<nTracksSignal<<endl;
  cout<<"Total number of tracks (background):"<<nTracksBackground<<endl;
  if(analyzeData) cout<<"Total number of tracks (data):"<<nTracksData<<endl;
  
  cout<<"Total number of short tracks (signal):"<<nShortTracksSignal<<endl;
  cout<<"Total number of short tracks (background):"<<nShortTracksBackground<<endl;
  if(analyzeData) cout<<"Total number of short tracks (data):"<<nShortTracksData<<endl;
  
  cout<<"Total number of short tracks above threshold (signal):"<<nShortTracksAboveThresholdSignal<<endl;
  cout<<"Total number of short tracks above threshold (background):"<<nShortTracksAboveThresholdBackground<<endl;
  if(analyzeData) cout<<"Total number of short tracks above threshold (data):"<<nShortTracksAboveThresholdData<<endl;
  
  cout<<"% or short tracks (signal):"<<nShortTracksSignal/(double)nTracksSignal<<endl;
  cout<<"% or short tracks (background):"<<nShortTracksBackground/(double)nTracksBackground<<endl;
  if(analyzeData) cout<<"% or short tracks (data):"<<nShortTracksData/(double)nTracksData<<endl;
  
  cout<<"% or short tracks above threshold (signal):"<<nShortTracksAboveThresholdSignal/(double)nTracksSignal<<endl;
  cout<<"% or short tracks above threshold (background):"<<nShortTracksAboveThresholdBackground/(double)nTracksBackground<<endl;
  if(analyzeData) cout<<"% or short tracks above threshold (data):"<<nShortTracksAboveThresholdData/(double)nTracksData<<endl;
  
  cout<<"\t\tchi:"<<1-nShortTracksAboveThresholdSignal/(double)nTracksSignal+2*nShortTracksAboveThresholdBackground/(double)nTracksBackground;
  cout<<"\t\tsignal:"<<nShortTracksAboveThresholdSignal/(double)nTracksSignal<<"\t\tback:"<<nShortTracksAboveThresholdBackground/(double)nTracksBackground<<endl;
  
  if(!drawPlots) return 0;
  
  // Per layer plots
  vector<TH1D*> dedxPerLayerSignal = GetDedxPerLayerHists(eventsSignal);
  vector<TH1D*> dedxPerLayerBackground = GetDedxPerLayerHists(eventsBackground);
  
  TCanvas *c1 = new TCanvas("c1","c1",2880,1800);
  c1->Divide(4,4);
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    c1->cd(iLayer+1);
    
    dedxPerLayerSignal[iLayer]->SetLineColor(kGreen+1);
    dedxPerLayerSignal[iLayer]->Draw();
    
    dedxPerLayerBackground[iLayer]->SetLineColor(kRed+1);
    dedxPerLayerBackground[iLayer]->Draw("same");
  }
  
  // Global plots
  TH1D *nClustersPerTrackSignal = new TH1D("N clusters per track (signal)","N clusters per track (signal)",20,0,20);
  TH1D *nClustersPerTrackBackground = new TH1D("N clusters per track","N clusters per track",20,0,20);
  nClustersPerTrackSignal->Sumw2();
  nClustersPerTrackBackground->Sumw2();
  
  TH1D *totalDeDxSignal = new TH1D("total dedx per track (signal)","total dedx per track (signal)",50,0,200);
  TH1D *totalDeDxBackground = new TH1D("total dedx per track","total dedx per track",50,0,200);
  totalDeDxSignal->Sumw2();
  totalDeDxBackground->Sumw2();
  
  TH1D *dedxPerClusterSignal = new TH1D("dedx per cluster (signal)","dedx per cluster (signal)",50,0,20);
  TH1D *dedxPerClusterBackground = new TH1D("dedx per cluster","dedx per cluster",50,0,20);
  dedxPerClusterSignal->Sumw2();
  dedxPerClusterBackground->Sumw2();
  
  TH1D *totalDeDxByNclustersSignal = new TH1D("total dedx per track / n clusters (signal)","total dedx per track / n clusters (signal)",50,0,20);
  TH1D *totalDeDxByNclustersBackground = new TH1D("total dedx per track / n clusters","total dedx per track / n clusters",50,0,20);
  totalDeDxByNclustersSignal->Sumw2();
  totalDeDxByNclustersBackground->Sumw2();
  
  TH1D *ptSignal = new TH1D("pt dist (signal)","pt dist (signal)",50,0,1000);
  TH1D *ptBackground = new TH1D("pt dist","pt dist",50,0,1000);
  ptSignal->Sumw2();
  ptBackground->Sumw2();

  TH1D *etaSignal = new TH1D("eta dist (signal)","eta dist (signal)",50,-3,3);
  TH1D *etaBackground = new TH1D("eta dist","eta dist",50,-3,3);
  etaSignal->Sumw2();
  etaBackground->Sumw2();
  
  TH1D *phiSignal = new TH1D("phi dist (signal)","phi dist (signal)",50,-3.5,3.5);
  TH1D *phiBackground = new TH1D("phi dist","phi dist",50,-3.5,3.5);
  phiSignal->Sumw2();
  phiBackground->Sumw2();
  
  // Fill signal histograms
  for(int iEvent=0;iEvent<eventsSignal->size();iEvent++){
    Event *event = eventsSignal->At(iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrackSignal->Fill(track->GetNclusters());
      totalDeDxSignal->Fill(track->GetTotalDedx());
      totalDeDxByNclustersSignal->Fill(track->GetTotalDedx()/track->GetNclusters());
      ptSignal->Fill(track->GetPt());
      etaSignal->Fill(track->GetEta());
      phiSignal->Fill(track->GetPhi());
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerClusterSignal->Fill(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  // Fill background histograms
  for(int iEvent=0;iEvent<eventsBackground->size();iEvent++){
    Event *event = eventsBackground->At(iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrackBackground->Fill(track->GetNclusters());
      totalDeDxBackground->Fill(track->GetTotalDedx());
      totalDeDxByNclustersBackground->Fill(track->GetTotalDedx()/track->GetNclusters());
      ptBackground->Fill(track->GetPt());
      etaBackground->Fill(track->GetEta());
      phiBackground->Fill(track->GetPhi());
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerClusterBackground->Fill(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  // Plot histograms
  TCanvas *c2 = new TCanvas("c2","c2",2880,1800);
  c2->Divide(3,3);
  TLegend *leg = GetLegend(0.15,0.5,0.75,0.25,"Data type");
  
  c2->cd(1);
  nClustersPerTrackSignal->Draw();
  nClustersPerTrackBackground->SetLineColor(kRed);
  nClustersPerTrackBackground->Draw("same");
  leg->AddEntry(nClustersPerTrackSignal,"Signal","lp");
  leg->AddEntry(nClustersPerTrackBackground,"Background","lp");
  leg->Draw();
  
  c2->cd(2);
  totalDeDxSignal->Draw();
  totalDeDxSignal->GetXaxis()->SetTitle("#sum dE/dx (MeV)");
  totalDeDxBackground->SetLineColor(kRed);
  totalDeDxBackground->Draw("same");
  leg->Draw();
  
  c2->cd(3);
  dedxPerClusterSignal->Draw();
  dedxPerClusterSignal->GetXaxis()->SetTitle("dE/dx (MeV)");
  dedxPerClusterBackground->SetLineColor(kRed);
  dedxPerClusterBackground->Draw("same");
  leg->Draw();
  
  c2->cd(4);
  totalDeDxByNclustersSignal->Draw();
  totalDeDxByNclustersSignal->GetXaxis()->SetTitle("dE/dx (MeV)");
  totalDeDxByNclustersBackground->SetLineColor(kRed);
  totalDeDxByNclustersBackground->Draw("same");
  leg->Draw();
  
  c2->cd(5);
  ptSignal->Draw();
  ptSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  ptBackground->SetLineColor(kRed);
  ptBackground->Draw("same");
  leg->Draw();
  
  c2->cd(6);
  etaSignal->Draw();
  etaSignal->GetXaxis()->SetTitle("#eta");
  etaBackground->SetLineColor(kRed);
  etaBackground->Draw("same");
  leg->Draw();
  
  c2->cd(7);
  phiSignal->Draw();
  phiSignal->GetXaxis()->SetTitle("#phi");
  phiBackground->SetLineColor(kRed);
  phiBackground->Draw("same");
  leg->Draw();
  
  c1->Update();
  c2->Update();
  theApp.Run();
  
  return 0;
}




// Get different track cuts

TrackCut* GetShortTrackCut()
{
  TrackCut *cut = new TrackCut();
  cut->SetNdedxClusters(3, 3);
  return cut;
}

TrackCut* GetShortAboveThresholdTrackCut()
{
  TrackCut *cut = new TrackCut();
  cut->SetNdedxClusters(3, 3);
  cut->SetMinDedxPerCluster(5.0);
  return cut;
}



