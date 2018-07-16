#include "Helpers.hpp"


const bool drawPlots = true;
const bool analyzeData = false;

const int nColors = 7;
const int colors[nColors] = {kGreen, kBlue, kRed, kOrange, kBlack, kMagenta, kCyan};

//double minDeDxOfChargino = 120000; // random
//double minDeDxOfChargino = 89500;  // minimizes 1-frac_signal+frac_back
double minDeDxOfChargino = 163000;   // minimizes 1-frac_signal+2*frac_back

vector<TH1D*> GetHistsFromEvents(vector<Event*> events)
{
  static int iter = 0;
  
  vector<TH1D*> hists;
  int iEvent = 0;
  for(auto event : events){
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      TH1D *dedxHist = new TH1D(Form("%i_dedx_ev%i",iter,iEvent),Form("%i_dedx_ev%i",iter,iEvent),nLayers,0,nLayers-1);
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        dedxHist->SetBinContent(iLayer,track->GetDeDxInLayer(iLayer));
      }
      
      dedxHist->SetLineColor(colors[iEvent % nColors]);
      hists.push_back(dedxHist);
    }
    iEvent++;
  }
  iter++;
  return hists;
}

vector<TH1D*> GetDedxPerLayerHists(vector<Event*> events)
{
  vector<TH1D*> hists;
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    TH1D *dedxHist = new TH1D(Form("dedx_layer[%i]",iLayer),Form("dedx_layer[%i]",iLayer),50,0,10);
    hists.push_back(dedxHist);
  }
  
  for(auto event : events){
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0001){
//          cout<<"Layer:"<<iLayer<<"\t"<<track->dedx[iLayer]<<endl;
          hists[iLayer]->Fill(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  return hists;
}

int GetNtracks(vector<Event*> events, int maxEvent=9999999)
{
  int nTracks = 0;
  int iEvent = 0;
  for(auto event : events){
    if(iEvent > maxEvent) break;
    nTracks+= event->GetNtracks();
    iEvent++;
  }
  return nTracks;
}

int GetNshortTracks(vector<Event*> events, int maxEvent=9999999)
{
  int nShortTracks = 0;
  int iEvent = 0;
  for(auto event : events){
    if(iEvent > maxEvent) break;
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      if(track->GetIsShort()) nShortTracks++;
    }
    iEvent++;
  }
  return nShortTracks;
}

int GetNshortTracksAboveThreshold(vector<Event*> events, int maxEvent=9999999)
{
  int nShortTracksAboveThreshold = 0;
  int iEvent = 0;
  for(auto event : events){
    if(iEvent > maxEvent) break;
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      if(!track->GetIsShort()) continue;
      
      int nPointsAboveThreshold=0;
    
      double totalDeDx = 0;
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        totalDeDx += track->GetDeDxInLayer(iLayer);
//        if(track->dedx[iLayer] > minDeDxOfChargino) nPointsAboveThreshold++;
      }
//      if(nPointsAboveThreshold == 3) nShortTracksAboveThreshold++;
      if(totalDeDx > minDeDxOfChargino) nShortTracksAboveThreshold++;
    }
    iEvent++;
  }
  return nShortTracksAboveThreshold;
}

int plotDeDx()
//int main()
{

  
  TLegend *leg = GetLegend(0.15,0.5,0.75,0.25,"Data type");
  
  const char *inFileNameSignal = "../jniedzie/mcSignal/tree.root";
  const char *inFileNameBackground = "../jniedzie/mcBackground/tree.root";

//  const char *inFileNameSignal = "../adish/Signal/tree.root";
//  const char *inFileNameBackground = "../adish/Background/tree.root";
  const char *inFileNameData = "../adish/Data/tree.root";
  
  cout<<"Reading signal events"<<endl;
  vector<Event*> eventsSignal = GetEventsVectorFromFile(inFileNameSignal);
  
  cout<<"Reading background events"<<endl;
  vector<Event*> eventsBackground = GetEventsVectorFromFile(inFileNameBackground);
  
  eventsBackground[0]->Print();
  
//  for(minDeDxOfChargino = 161000;minDeDxOfChargino<166000;minDeDxOfChargino+=500){
  
  cout<<"minDeDxOfChargino:"<<minDeDxOfChargino<<endl;
  
  
  int nTracksSignal = GetNtracks(eventsSignal);
  int nShortTracksSignal = GetNshortTracks(eventsSignal);
  int nShortTracksAboveThresholdSignal = GetNshortTracksAboveThreshold(eventsSignal);
  
  int nTracksBackground = GetNtracks(eventsBackground);
  int nShortTracksBackground = GetNshortTracks(eventsBackground);
  int nShortTracksAboveThresholdBackground = GetNshortTracksAboveThreshold(eventsBackground);
  
  int nTracksData, nShortTracksData, nShortTracksAboveThresholdData;
  
  if(analyzeData){
    cout<<"Reading data events"<<endl;
    vector<Event*> eventsData = GetEventsVectorFromFile(inFileNameData);

    nTracksData = GetNtracks(eventsData);
    nShortTracksData = GetNshortTracks(eventsData);
    nShortTracksAboveThresholdData = GetNshortTracksAboveThreshold(eventsData);
  }
    
  cout<<"Total number of tracks in signal events:"<<nTracksSignal<<endl;
  cout<<"Total number of tracks in background events:"<<nTracksBackground<<endl;
  if(analyzeData) cout<<"Total number of tracks in data events:"<<nTracksData<<endl;
  
  cout<<"Total number of short (max 3 dedx points) tracks in signal events:"<<nShortTracksSignal<<endl;
  cout<<"Total number of short (max 3 dedx points) tracks in background events:"<<nShortTracksBackground<<endl;
  if(analyzeData) cout<<"Total number of short (max 3 dedx points) tracks in data events:"<<nShortTracksData<<endl;
  
  cout<<"Total number of short (max 3 dedx points) tracks above threshold in signal events:"<<nShortTracksAboveThresholdSignal<<endl;
  cout<<"Total number of short (max 3 dedx points) tracks above threshold in background events:"<<nShortTracksAboveThresholdBackground<<endl;
  if(analyzeData) cout<<"Total number of short (max 3 dedx points) tracks above threshold in data events:"<<nShortTracksAboveThresholdData<<endl;
  
  cout<<"% or short tracks in signal:"<<nShortTracksSignal/(double)nTracksSignal<<endl;
  cout<<"% or short tracks in background:"<<nShortTracksBackground/(double)nTracksBackground<<endl;
  if(analyzeData) cout<<"% or short tracks in data:"<<nShortTracksData/(double)nTracksData<<endl;
  
  cout<<"% or short tracks above threshold in signal:"<<nShortTracksAboveThresholdSignal/(double)nTracksSignal<<endl;
  cout<<"% or short tracks above threshold in background:"<<nShortTracksAboveThresholdBackground/(double)nTracksBackground<<endl;
  if(analyzeData) cout<<"% or short tracks above threshold in data:"<<nShortTracksAboveThresholdData/(double)nTracksData<<endl;
  
  cout<<"\t\tchi:"<<1-nShortTracksAboveThresholdSignal/(double)nTracksSignal+2*nShortTracksAboveThresholdBackground/(double)nTracksBackground;
  cout<<"\t\tsignal:"<<nShortTracksAboveThresholdSignal/(double)nTracksSignal<<"\t\tback:"<<nShortTracksAboveThresholdBackground/(double)nTracksBackground<<endl;
  
//  }
//  cout<<"S/B corected (short tracks) = "<<eventsBackground.GetNtracks()/(double)eventsSignal.size()*nShortTracksSignal/(double)nShortTracksBackground<<endl;
//  cout<<"S/B corected (short tracks above threshold) = "<<eventsBackground.size()/(double)eventsSignal.size()*nShortTracksAboveThresholdSignal/(double)nShortTracksAboveThresholdBackground<<endl;

  
  // Plotting part
  
  if(!drawPlots) return 0;

//  vector<TH1D*> histsSignal = GetHistsFromEvents(eventsSignal);
//  vector<TH1D*> histsBackground = GetHistsFromEvents(eventsBackground);
//  if(analyzeData) vector<TH1D*> histsData = GetHistsFromEvents(eventsBackground);
  
//  cout<<"Dedx:"<<eventsSignal[0]->GetTrack(0)->dedx[0]<<endl;
  
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
  
  TCanvas *c2 = new TCanvas("c2","c2",2880,1800);
  c2->Divide(2,2);
  
  TH1D *nClustersPerTrackSignal = new TH1D("N clusters per track","N clusters per track",20,0,20);
  TH1D *nClustersPerTrackBackground = new TH1D("N clusters per track","N clusters per track",20,0,20);
  nClustersPerTrackSignal->Sumw2();
  nClustersPerTrackBackground->Sumw2();
  
  TH1D *totalDeDxSignal = new TH1D("total dedx per track","total dedx per track",50,0,200);
  TH1D *totalDeDxBackground = new TH1D("total dedx per track","total dedx per track",50,0,200);
  totalDeDxSignal->Sumw2();
  totalDeDxBackground->Sumw2();
  
  TH1D *dedxPerClusterSignal = new TH1D("dedx per cluster","dedx per cluster",50,0,20);
  TH1D *dedxPerClusterBackground = new TH1D("dedx per cluster","dedx per cluster",50,0,20);
  dedxPerClusterSignal->Sumw2();
  dedxPerClusterBackground->Sumw2();
  
  TH1D *totalDeDxByNclustersSignal = new TH1D("total dedx per track / n clusters","total dedx per track / n clusters",50,0,20);
  TH1D *totalDeDxByNclustersBackground = new TH1D("total dedx per track / n clusters","total dedx per track / n clusters",50,0,20);
  totalDeDxByNclustersSignal->Sumw2();
  totalDeDxByNclustersBackground->Sumw2();
  
  for(Event *event : eventsSignal){
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrackSignal->Fill(track->GetNclusters());
      totalDeDxSignal->Fill(track->GetTotalDedx());
      totalDeDxByNclustersSignal->Fill(track->GetTotalDedx()/track->GetNclusters());
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerClusterSignal->Fill(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  for(Event *event : eventsBackground){
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrackBackground->Fill(track->GetNclusters());
      totalDeDxBackground->Fill(track->GetTotalDedx());
      totalDeDxByNclustersBackground->Fill(track->GetTotalDedx()/track->GetNclusters());
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerClusterBackground->Fill(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  c2->cd(1);
  nClustersPerTrackSignal->Draw();
  nClustersPerTrackBackground->SetLineColor(kRed);
  nClustersPerTrackBackground->Draw("same");
  leg->AddEntry(nClustersPerTrackSignal,"Signal","lp");
  leg->AddEntry(nClustersPerTrackBackground,"Background","lp");
  leg->Draw();
  
  c2->cd(2);
  totalDeDxSignal->Draw();
  totalDeDxBackground->SetLineColor(kRed);
  totalDeDxBackground->Draw("same");
  leg->Draw();
  
  c2->cd(3);
  dedxPerClusterSignal->Draw();
  dedxPerClusterBackground->SetLineColor(kRed);
  dedxPerClusterBackground->Draw("same");
  leg->Draw();
  
  c2->cd(4);
  totalDeDxByNclustersSignal->Draw();
  totalDeDxByNclustersBackground->SetLineColor(kRed);
  totalDeDxByNclustersBackground->Draw("same");
  leg->Draw();
  
  
//  TH1D *dedxMeanSignal = new TH1D("dedxMeanSignal","dedxMeanSignal",20,0,2.5);
//  dedxMeanSignal->Sumw2();
//  TH1D *dedxWidthSignal = new TH1D("dedxWidthSignal","dedxWidthSignal",20,0,2.5);
//  dedxWidthSignal->Sumw2();
//  TH1D *dedxIntegralSignal = new TH1D("dedxIntegralSignal","dedxIntegralSignal",20,0,300e3);
//  dedxIntegralSignal->Sumw2();
//
//  TH1D *dedxMeanBackgrund = new TH1D("dedxMeanBackgrund","dedxMeanBackgrund",20,0,2.5);
//  dedxMeanBackgrund->Sumw2();
//  TH1D *dedxWidthBackgrund = new TH1D("dedxWidthBackgrund","dedxWidthBackgrund",20,0,2.5);
//  dedxWidthBackgrund->Sumw2();
//  TH1D *dedxIntegralBackgrund = new TH1D("dedxIntegralBackgrund","dedxIntegralBackgrund",20,0,300e3);
//  dedxIntegralBackgrund->Sumw2();
//
//  bool firstHist = true;
//  c1->cd(1);
//  for(TH1D *hist : histsSignal){
//    dedxMeanSignal->Fill(hist->GetMean());
//    dedxWidthSignal->Fill(hist->GetStdDev());
//    dedxIntegralSignal->Fill(hist->Integral());
//
//    if(firstHist){
//      hist->Draw();
//      firstHist = false;
//    }
//    else{
//      hist->Draw("same");
//    }
//
//  }
//
//  for(TH1D *hist : histsBackground){
//    dedxMeanBackgrund->Fill(hist->GetMean());
//    dedxWidthBackgrund->Fill(hist->GetStdDev());
//    dedxIntegralBackgrund->Fill(hist->Integral());
//  }
//
//
//  c1->cd(2);
//  dedxMeanSignal->DrawNormalized();
//  dedxMeanBackgrund->SetLineColor(kRed);
//  dedxMeanBackgrund->DrawNormalized("same");
//  leg->AddEntry(dedxMeanSignal,"Signal","lp");
//  leg->AddEntry(dedxMeanBackgrund,"Background","lp");
//  leg->Draw();
//
//  c1->cd(3);
//  dedxWidthSignal->DrawNormalized();
//  dedxWidthBackgrund->SetLineColor(kRed);
//  dedxWidthBackgrund->DrawNormalized("same");
//  leg->Draw();
//
//  c1->cd(4);
//  dedxIntegralSignal->DrawNormalized();
//  dedxIntegralBackgrund->SetLineColor(kRed);
//  dedxIntegralBackgrund->DrawNormalized("same");
//  leg->Draw();
  
  return 0;
}
