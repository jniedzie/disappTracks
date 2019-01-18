//
//  MonitorsManager.cpp
//
//  Created by Jeremi Niedziela on 18/01/2019.
//

#include "MonitorsManager.hpp"

MonitorsManager::MonitorsManager(const shared_ptr<FitterConfig> &_config) :
config(_config)
{
  const vector<tuple<const char*,int,double,double>> monitors1Dparams = {
    {"nPointsOnHelix",100 ,0,100},
    {"chi2ofHelix",   50  ,0,50 },
    {"nPionPoints",   100, 0,1.2},
    {"nFakeHits",     1000,0,10 },
    {"failReason",    10,  0,10 },
  };
  
  double maxPx = config->GetMaxPx();
  double maxPy = config->GetMaxPy();
  double maxPz = config->GetMaxPz();
  
  const vector<tuple<const char*,int,double,double,int,double,double>> monitors2Dparams = {
    {"xResponse",     500,-250,250, 500,-250,250 },
    {"yResponse",     500,-250,250, 500,-250,250 },
    {"zResponse",     500,-250,250, 500,-250,250 },
    {"pxResponse",    200,-maxPx,maxPx, 200,-maxPx,maxPx },
    {"pyResponse",    200,-maxPy,maxPy, 200,-maxPy,maxPy },
    {"pzResponse",    200,-maxPz,maxPz, 200,-maxPz,maxPz },
  };
  
  const vector<tuple<const char*,int,double,double>> fractionMonitorsParams = {
    {"successVsPz",     config->GetMaxPz()-config->GetMinPz()+1, config->GetMinPz(), config->GetMaxPz()},
    {"fullSuccessVsPz", config->GetMaxPz()-config->GetMinPz()+1, config->GetMinPz(), config->GetMaxPz()},
  };
  
  for(auto params : monitors1Dparams){
    monitors1D[get<0>(params)] = new TH1D(get<0>(params),get<0>(params),get<1>(params),get<2>(params),get<3>(params));
  }
  
  for(auto params : monitors2Dparams){
    monitors2D[get<0>(params)] = new TH2D(get<0>(params),get<0>(params),
                                          get<1>(params),get<2>(params),get<3>(params),
                                          get<4>(params),get<5>(params),get<6>(params));
  }
  
  for(auto params : fractionMonitorsParams){
    fractionMonitors[get<0>(params)].first  = new TH1D(get<0>(params),get<0>(params),
                                                       get<1>(params),get<2>(params),get<3>(params));
    
    fractionMonitors[get<0>(params)].second = new TH1D(get<0>(params),get<0>(params),
                                                       get<1>(params),get<2>(params),get<3>(params));
    
    fractionMonitors[get<0>(params)].first->Sumw2();
    fractionMonitors[get<0>(params)].second->Sumw2();
  }
  
}

MonitorsManager::~MonitorsManager()
{
  
}

void MonitorsManager::FillMonitors(const unique_ptr<Helix> &fittedHelix, const unique_ptr<Helix> &trueHelix)
{
  fractionMonitors["successVsPz"].second->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  fractionMonitors["fullSuccessVsPz"].second->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  
  if(!fittedHelix){
    monitors1D["failReason"]->Fill(8);
    return;
  }
  
  // Here (full) success part starts:
  fractionMonitors["successVsPz"].first->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  
  monitors2D["xResponse"]->Fill(trueHelix->GetOrigin()->GetX(), fittedHelix->GetOrigin()->GetX());
  monitors2D["yResponse"]->Fill(trueHelix->GetOrigin()->GetY(), fittedHelix->GetOrigin()->GetY());
  monitors2D["zResponse"]->Fill(trueHelix->GetOrigin()->GetZ(), fittedHelix->GetOrigin()->GetZ());
  monitors2D["pxResponse"]->Fill(trueHelix->GetMomentum()->GetX(), fittedHelix->GetMomentum()->GetX());
  monitors2D["pyResponse"]->Fill(trueHelix->GetMomentum()->GetY(), fittedHelix->GetMomentum()->GetY());
  monitors2D["pzResponse"]->Fill(trueHelix->GetMomentum()->GetZ(), fittedHelix->GetMomentum()->GetZ());
  monitors1D["nPointsOnHelix"]->Fill(fittedHelix->GetNpoints());
  monitors1D["chi2ofHelix"]->Fill(fittedHelix->GetChi2() < 50 ? fittedHelix->GetChi2() : 49);
  monitors1D["nPionPoints"]->Fill(fittedHelix->GetNpionPoints()/(double)trueHelix->GetNpionPoints());
  monitors1D["nFakeHits"]->Fill((fittedHelix->GetNpoints()-fittedHelix->GetNpionPoints())/(double)fittedHelix->GetNpoints());
  
  vector<int> failureCodes = Helix::AreHelicesIdentical(fittedHelix, trueHelix);
  
  if(failureCodes.size()!=0){ // if helix was fitted, but there was some discrepancy wrt the true helix
    for(int f : failureCodes) monitors1D["failReason"]->Fill(f);
  }
  else{ // full success case
    fractionMonitors["fullSuccessVsPz"].first->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  }
}

MonitorsManager::EFitStatus MonitorsManager::GetFittingStatus(const unique_ptr<Helix> &fittedHelix,
                                                              const unique_ptr<Helix> &trueHelix)
{
  if(!fittedHelix) return kFail;
  vector<int> failureCodes = Helix::AreHelicesIdentical(fittedHelix, trueHelix);
  if(failureCodes.size()==0) return kFullSuccess;
  return kSuccess;
}

void MonitorsManager::PlotAndSaveMonitors()
{
  // Regular 1D and 2D monitors
  TCanvas *c1 = new TCanvas("Monitors","Monitors",2880,1800);
  c1->Divide(4,4);
  TFile *outFile = new TFile(config->GetOutputPath(),"recreate");
  outFile->cd();
  
  int i=1;
  for(auto &[title, hist] : monitors2D){
    c1->cd(i++);
    hist->Draw("colz");
    hist->Write();
  }
  
  for(auto &[title, hist] : monitors1D){
    c1->cd(i++);
    hist->Draw();
    hist->Write();
  }
  
  // Fraction monitors:
  for(auto hists : fractionMonitors){
    hists.second.first->Divide(hists.second.second);
  }
  
  TCanvas *c2 = new TCanvas("Fraction monitors","Fraction monitors",2880,1800);
  c2->Divide(2,2);
  i=1;
  for(auto &[title, hist] : fractionMonitors){
    c2->cd(i++);
    hist.first->Draw();
    hist.first->Write();
  }
  
  outFile->Close();
  c1->Update();
  c2->Update();
}
