//  getMetWeights.cpp
//  Created by Jeremi Niedziela on 16/01/2020.

#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

string configPath = "configs/analysis.md";
 
/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix = "";
   
  if(config.secondaryCategory == "Zmumu")   prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu")   prefix += "Wmunu/";
  if(config.secondaryCategory == "LowMET")   prefix += "LowMET/";
  
  if(config.params["cuts_level"]==0) prefix += "after_L0/";
  if(config.params["cuts_level"]==1) prefix += "after_L1/"+config.category+"/";
  
  return prefix;
}

void fillHistogram(EventSet &events, TH1D *hist, ESignal iSig, string mode)
{
  double sumW = 0;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(int iEvent=0; iEvent<events.size(kSignal, iSig, year); iEvent++){
      auto event = events.At(kSignal, iSig, year, iEvent);
      
      if(mode == "charginos"){
        
      }
      else if(mode == "gen_met") hist->Fill(event->GetMetGenPt(), event->GetWeight());
      else if(mode == "met")     hist->Fill(event->GetMetNoMuPt(), event->GetWeight());
      else{
        cout<<"ERROR -- unknown mode: "<<mode<<endl;
        return;
      }
      
      sumW += event->GetWeight();
    }
  }
  hist->Scale(1./sumW);
}

void reweightHistogram(EventSet &events, TH1D *hist, TH1D *ratio, ESignal iSig)
{
  double sumW = 0;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(int iEvent=0; iEvent<events.size(kSignal, iSig, year); iEvent++){
      auto event = events.At(kSignal, iSig, year, iEvent);
      
      double met = event->GetMetNoMuPt();
      double scale = ratio->GetBinContent(ratio->GetXaxis()->FindFixBin(met));
      
      hist->Fill(event->GetMetNoMuPt(), event->GetWeight()*scale);
      sumW += event->GetWeight()*scale;
    }
  }
  hist->Scale(1./sumW);
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  gStyle->SetOptStat(0);
  
  // All events with initial cuts only
  config = ConfigManager(configPath);
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  const int nBins = 50;
  
  TH1D *charginosPtGiovanni = new TH1D("charginos' pt" , "charginos' pt"  , nBins, 0, 1000);
  TH1D *charginosPtFilip    = new TH1D("charginos' pt ", "charginos' pt " , nBins, 0, 1000);
  TH1D *genMetGiovanni      = new TH1D("gen-MET"       , "gen-MET"        , nBins, 0, 1000);
  TH1D *genMetFilip         = new TH1D("gen-Met "      , "gen-Met "       , nBins, 0, 1000);
  TH1D *metGiovanni         = new TH1D("MET"           , "MET"            , nBins, 0, 1000);
  TH1D *metFilip            = new TH1D("Met "          , "Met "           , nBins, 0, 1000);
  
  charginosPtGiovanni->SetLineColor(kBlack);
  genMetGiovanni->SetLineColor(kBlack);
  metGiovanni->SetLineColor(kBlack);
  charginosPtFilip->SetLineColor(kRed);
  genMetFilip->SetLineColor(kRed);
  metFilip->SetLineColor(kRed);
  
  fillHistogram(events, charginosPtGiovanni , kWino_M_500_cTau_10 , "charginos" );
  fillHistogram(events, charginosPtFilip    , kChargino500_10     , "charginos" );
  fillHistogram(events, genMetGiovanni      , kWino_M_500_cTau_10 , "gen_met"   );
  fillHistogram(events, genMetFilip         , kChargino500_10     , "gen_met"   );
  fillHistogram(events, metGiovanni         , kWino_M_500_cTau_10 , "met"       );
  fillHistogram(events, metFilip            , kChargino500_10     , "met"       );

  TCanvas *c1 = new TCanvas("c1", "c1", 2880, 1800);
  c1->Divide(3,2);
  
  
  //
  //  Charginos' pt
  //
  
  TLegend *leg = new TLegend(0.85, 1.0, 0.85, 1.0);
  leg->AddEntry(charginosPtGiovanni , "charginos' pt MadGraph"  , "epl");
  leg->AddEntry(charginosPtFilip    , "charignos' pt Pythia"    , "epl");
  
  c1->cd(1);
  charginosPtGiovanni->Draw();
  charginosPtFilip->Draw("same");
  leg->Draw();
  
  c1->cd(4);
  TH1D *charginosPtRatio = new TH1D(*charginosPtGiovanni);
  charginosPtRatio->SetTitle("charginos' pt ratio");
  charginosPtRatio->SetLineColor(kBlue);
  charginosPtRatio->Sumw2();
  charginosPtRatio->Divide(charginosPtFilip);
  charginosPtRatio->Draw();
  
  //
  //  gen-MET
  //
  
  
  TLegend *leg3 = new TLegend(0.85, 1.0, 0.85, 1.0);
  leg3->AddEntry(genMetGiovanni      , "gen MET MadGraph"        , "epl");
  leg3->AddEntry(genMetFilip         , "gen MET Pythia"          , "epl");
  
  c1->cd(2);
  genMetGiovanni->Draw();
  genMetFilip->Draw("same");
  leg3->Draw();
  
  c1->cd(5);
  TH1D *genMetRatio = new TH1D(*genMetGiovanni);
  genMetRatio->SetTitle("gen-MET ratio");
  genMetRatio->SetLineColor(kBlue);
  genMetRatio->Sumw2();
  genMetRatio->Divide(genMetFilip);
  genMetRatio->Draw();
  
  
  //
  // MET
  //
  

  c1->cd(3);
  TH1D *metFilipReweighted = new TH1D(*metFilip);
  reweightHistogram(events, metFilipReweighted, genMetRatio, kChargino500_10);
  metGiovanni->Draw();
  metFilipReweighted->Draw("same");
  
  TLegend *leg2 = new TLegend(0.85, 1.0, 0.85, 1.0);
  leg2->AddEntry(metFilipReweighted , "reweighted MET Pythia", "epl");
  leg2->AddEntry(metGiovanni        , "MET MadGraph"        , "epl");
  
  leg2->Draw();
  
  c1->cd(6);
  TH1D *metReweightedRatio = new TH1D(*metGiovanni);
  metReweightedRatio->SetLineColor(kBlue);
  metReweightedRatio->Sumw2();
  metReweightedRatio->Divide(metFilipReweighted);
  metReweightedRatio->Draw();
  
  
  //
  // Save to file
  //
  
//  charginosPtRatio->SetTitle("metRatio");
//  charginosPtRatio->SetName("metRatio");
//  metRatio->SaveAs("../data/SIG-SR/metWeights.root");
  
  c1->Update();
  theApp.Run();
  return 0;
}

