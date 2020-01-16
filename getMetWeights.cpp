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

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  gStyle->SetOptStat(0);
  
  // All events with initial cuts only
  config = ConfigManager(configPath);
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  TH1D *metGiovanni = new TH1D("metGiovanni", "metGiovanni", 100, 0, 1000);
  TH1D *metFilip    = new TH1D("metFilip", "metFilip", 100, 0, 1000);
  metFilip->SetLineColor(kRed);
  
  double sumWgiovanni = 0;
  double sumWfilip = 0;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    cout<<"Year: "<<year<<endl;
    
    ESignal iSigGiovanni = kWino_M_500_cTau_10;
    ESignal iSigFilip = kChargino500_10;
    
    for(int iEvent=0; iEvent<events.size(kSignal, iSigGiovanni, year); iEvent++){
      auto event = events.At(kSignal, iSigGiovanni, year, iEvent);
      metGiovanni->Fill(event->GetMetNoMuPt(), event->GetWeight());
      sumWgiovanni += event->GetWeight();
    }
    
    for(int iEvent=0; iEvent<events.size(kSignal, iSigFilip, year); iEvent++){
      auto event = events.At(kSignal, iSigFilip, year, iEvent);
      metFilip->Fill(event->GetMetNoMuPt(), event->GetWeight());
      sumWfilip += event->GetWeight();
    }
  }
  cout<<"Done"<<endl;
  
  metGiovanni->Scale(1./sumWgiovanni);
  metFilip->Scale(1./sumWfilip);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 2880, 1800);
  c1->Divide(2,2);
  
  TLegend *leg = new TLegend(0.85, 1.0, 0.85, 1.0);
  leg->AddEntry(metGiovanni, "MET Giovanni", "epl");
  leg->AddEntry(metFilip, "MET Filip", "epl");
  
  c1->cd(1);
  
  metGiovanni->Draw();
  metFilip->Draw("same");
  leg->Draw();
  
  
  
  TH1D *metRatio = new TH1D(*metGiovanni);
  metRatio->Sumw2();
  metRatio->Divide(metFilip);
  c1->cd(3); metRatio->Draw();
  
  TH1D *metFilipReweighted = new TH1D("metFilipReweighted", "metFilipReweighted", 100, 0, 1000);
  metFilipReweighted->SetLineColor(kRed);
  
  double sumWfilipReweighted = 0;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    cout<<"Year: "<<year<<endl;
    
    ESignal iSigFilip = kChargino500_10;
    
    for(int iEvent=0; iEvent<events.size(kSignal, iSigFilip, year); iEvent++){
      auto event = events.At(kSignal, iSigFilip, year, iEvent);
      
      double met = event->GetMetNoMuPt();
      double scale = metRatio->GetBinContent(metRatio->GetXaxis()->FindFixBin(met));
      
      metFilipReweighted->Fill(event->GetMetNoMuPt(), event->GetWeight()*scale);
      sumWfilipReweighted += event->GetWeight()*scale;
    }
  }
  metFilipReweighted->Scale(1./sumWfilipReweighted);
  
  c1->cd(2);
  metGiovanni->Draw();
  metFilipReweighted->Draw("same");
  
  TH1D *metRatioReweighted = new TH1D(*metGiovanni);
  metRatioReweighted->Sumw2();
  metRatioReweighted->Divide(metFilipReweighted);
  c1->cd(4); metRatioReweighted->Draw();
  
  metRatio->SetTitle("metRatio");
  metRatio->SetName("metRatio");
//  metRatio->SaveAs("../data/SIG-SR/metWeights.root");
  
  c1->Update();
  theApp.Run();
  return 0;
}

