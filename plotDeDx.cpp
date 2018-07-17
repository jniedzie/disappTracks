#include "Helpers.hpp"
#include "Event.hpp"
#include "TrackCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

const bool drawPlots = true;
const bool analyzeData = false;

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

void GetSizeXperLayerHists(Events *events, vector<TH1D*> &histsX, vector<TH1D*> &histsY)
{
//  vector<TH1D*> hists;
  static int iter=0;
  
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    TH1D *sizeXhist = new TH1D(Form("sizeX_layer[%i]_%i",iLayer,iter),Form("sizeX_layer[%i]_%i",iLayer,iter),10,0,10);
    histsX.push_back(sizeXhist);
    
    TH1D *sizeYhist = new TH1D(Form("sizeY_layer[%i]_%i",iLayer,iter),Form("sizeY_layer[%i]_%i",iLayer,iter),10,0,10);
    histsY.push_back(sizeYhist);
  }
  
  for(int iEvent=0;iEvent<events->size();iEvent++){
    for(int iTrack=0;iTrack<events->At(iEvent)->GetNtracks();iTrack++){
      Track *track = events->At(iEvent)->GetTrack(iTrack);
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetSizeXinLayer(iLayer) > 0.0001){
          histsX[iLayer]->Fill(track->GetSizeXinLayer(iLayer));
        }
        if(track->GetSizeYinLayer(iLayer) > 0.0001){
          histsY[iLayer]->Fill(track->GetSizeYinLayer(iLayer));
        }
      }
    }
  }
  iter++;
}


int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  string inFileNameSignal = "../jniedzie/mcSignal/tree.root";
  string inFileNameBackground = "../jniedzie/mcBackground/tree.root";
//  const char *inFileNameSignal = "../adish/Signal/tree.root";
//  const char *inFileNameBackground = "../adish/Background/tree.root";
  string inFileNameData = "../adish/Data/tree.root";
  
  
  TrackCut *shortTrackCut = new TrackCut(TrackCut::kShort);
  TrackCut *shortAboveTrasholdTrackCut = new TrackCut(TrackCut::kShortAboveThreshold);
  
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
  
  vector<TH1D*> sizeXperLayerSignal, sizeYperLayerSignal;
  vector<TH1D*> sizeXperLayerBackground, sizeYperLayerBackground;
  
  GetSizeXperLayerHists(eventsSignal, sizeXperLayerSignal, sizeYperLayerSignal);
  GetSizeXperLayerHists(eventsBackground, sizeXperLayerBackground, sizeYperLayerBackground);
  
  TCanvas *c1 = new TCanvas("c1","c1",2880,1800);
  c1->Divide(4,4);
  
  TCanvas *c4 = new TCanvas("c4","c4",2880,1800);
  c4->Divide(4,2);
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    c1->cd(iLayer+1);
    dedxPerLayerSignal[iLayer]->SetLineColor(kGreen+1);
    dedxPerLayerSignal[iLayer]->Draw();
    
    dedxPerLayerBackground[iLayer]->SetLineColor(kRed+1);
    dedxPerLayerBackground[iLayer]->Draw("same");
    
    if(iLayer > 3) continue;
    
    c4->cd(iLayer+1);
    sizeXperLayerSignal[iLayer]->SetLineColor(kGreen+1);
    sizeXperLayerSignal[iLayer]->Draw();
    sizeXperLayerBackground[iLayer]->SetLineColor(kRed+1);
    sizeXperLayerBackground[iLayer]->Draw("same");
    
    c4->cd(iLayer+1+4);
    sizeYperLayerSignal[iLayer]->SetLineColor(kGreen+1);
    sizeYperLayerSignal[iLayer]->Draw();
    sizeYperLayerBackground[iLayer]->SetLineColor(kRed+1);
    sizeYperLayerBackground[iLayer]->Draw("same");
    
  }
  c1->Update();
  c4->Update();
  
  // Global plots
  HistSet *nClustersPerTrack = new HistSet("N clusters per track",20,0,20);
  HistSet *totalDeDx = new HistSet("total dedx per track",50,0,200);
  HistSet *dedxPerCluster = new HistSet("dedx per cluster",50,0,20);
  HistSet *totalDeDxByNclusters = new HistSet("total dedx per track / n clusters",50,0,20);
  HistSet *pt = new HistSet("pt dist",50,0,1000);
  HistSet *eta = new HistSet("eta dist",50,-3,3);
  HistSet *phi = new HistSet("phi dist",50,-3.5,3.5);
  HistSet *caloEm = new HistSet("EM calo energy",100,0,500);
  HistSet *caloHad = new HistSet("Hadron calo energy",100,0,500);
  
  HistSet *dxy = new HistSet("Displacement in XY",100,-0.02,0.02);
  HistSet *dz = new HistSet("Displacement in Z",100,-0.02,0.02);
  HistSet *charge = new HistSet("Charge",100,-10,10);
  HistSet *mass = new HistSet("Mass",500,0,0.25);
  HistSet *pid = new HistSet("PDG PID",441,-220,220);
  
  // Fill signal histograms
  for(int iEvent=0;iEvent<eventsSignal->size();iEvent++){
    Event *event = eventsSignal->At(iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrack->FillSignal(track->GetNclusters());
      totalDeDx->FillSignal(track->GetTotalDedx());
      totalDeDxByNclusters->FillSignal(track->GetTotalDedx()/track->GetNclusters());
      pt->FillSignal(track->GetPt());
      eta->FillSignal(track->GetEta());
      phi->FillSignal(track->GetPhi());
      caloEm->FillSignal(track->GetCaloEmEnergy());
      caloHad->FillSignal(track->GetCaloHadEnergy());
      dxy->FillSignal(track->GetDxy());
      dz->FillSignal(track->GetDz());
      charge->FillSignal(track->GetCharge());
      mass->FillSignal(track->GetMass());
      pid->FillSignal(track->GetPid());
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerCluster->FillSignal(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  // Fill background histograms
  for(int iEvent=0;iEvent<eventsBackground->size();iEvent++){
    Event *event = eventsBackground->At(iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrack->FillBackground(track->GetNclusters());
      totalDeDx->FillBackground(track->GetTotalDedx());
      totalDeDxByNclusters->FillBackground(track->GetTotalDedx()/track->GetNclusters());
      pt->FillBackground(track->GetPt());
      eta->FillBackground(track->GetEta());
      phi->FillBackground(track->GetPhi());
      caloEm->FillBackground(track->GetCaloEmEnergy());
      caloHad->FillBackground(track->GetCaloHadEnergy());
      dxy->FillBackground(track->GetDxy());
      dz->FillBackground(track->GetDz());
      charge->FillBackground(track->GetCharge());
      mass->FillBackground(track->GetMass());
      pid->FillBackground(track->GetPid());
      
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerCluster->FillBackground(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  // Plot histograms
  TCanvas *c2 = new TCanvas("c2","c2",2880,1800);
  c2->Divide(3,3);
  
  nClustersPerTrack->Draw(c2,1);
  totalDeDx->Draw(c2,2);
  dedxPerCluster->Draw(c2,3);
  totalDeDxByNclusters->Draw(c2,4);
  pt->Draw(c2,5);
  eta->Draw(c2,6);
  phi->Draw(c2,7);
  caloEm->Draw(c2,8);
  caloHad->Draw(c2,9);
  
  TCanvas *c3 = new TCanvas("c3","c3",2880,1800);
  c3->Divide(2,3);
  dxy->Draw(c3,1);
  dz->Draw(c3,2);
  charge->Draw(c3,3);
  mass->Draw(c3,4);
  pid->Draw(c3,5);
  
  theApp.Run();
  
  return 0;
}



