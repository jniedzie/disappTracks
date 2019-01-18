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

  double maxX = 250;
  double maxY = 250;
  double maxZ = 250;
  double maxEta = config->GetMaxTrackEta();
  
  double minPx = config->GetMinPx();
  double minPy = config->GetMinPy();
  double minPz = config->GetMinPz();
  double maxPx = config->GetMaxPx();
  double maxPy = config->GetMaxPy();
  double maxPz = config->GetMaxPz();
  
  
  const vector<tuple<const char*,int,double,double,int,double,double>> monitors2Dparams = {
    {"xResponse",     500, -maxX,  maxX,  500, -maxX,  maxX},
    {"yResponse",     500, -maxY,  maxY,  500, -maxY,  maxY},
    {"zResponse",     500, -maxZ,  maxZ,  500, -maxZ,  maxZ},
    {"pxResponse",    200, -maxPx, maxPx, 200, -maxPx, maxPx},
    {"pyResponse",    200, -maxPy, maxPy, 200, -maxPy, maxPy},
    {"pzResponse",    200, -maxPz, maxPz, 200, -maxPz, maxPz},
  };
  
  const vector<tuple<string,int,double,double>> fractionMonitorsParams = {
    {"successVsX",      maxX+1, 0, maxX},
    {"fullSuccessVsX",  maxX+1, 0, maxX},
    {"successVsY",      maxY+1, 0, maxY},
    {"fullSuccessVsY",  maxY+1, 0, maxY},
    {"successVsZ",      maxZ+1, 0, maxZ},
    {"fullSuccessVsZ",  maxZ+1, 0, maxZ},
    {"successVsPx",     maxPx-minPx+1, minPx, maxPx},
    {"fullSuccessVsPx", maxPx-minPx+1, minPx, maxPx},
    {"successVsPy",     maxPy-minPy+1, minPy, maxPy},
    {"fullSuccessVsPy", maxPy-minPy+1, minPy, maxPy},
    {"successVsPz",     maxPz-minPz+1, minPz, maxPz},
    {"fullSuccessVsPz", maxPz-minPz+1, minPz, maxPz},
    {"successVsEta",    100, 0, maxEta},
    {"fullSuccessVsEta",100, 0, maxEta},
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
    fractionMonitors[get<0>(params)].first  = new TH1D(get<0>(params).c_str(),get<0>(params).c_str(),
                                                       get<1>(params),get<2>(params),get<3>(params));
    
    fractionMonitors[get<0>(params)].second = new TH1D((get<0>(params)+"_den").c_str(),(get<0>(params)+"_den").c_str(),
                                                       get<1>(params),get<2>(params),get<3>(params));
    
    fractionMonitors[get<0>(params)].first->Sumw2();
    fractionMonitors[get<0>(params)].second->Sumw2();
  }
  
}

MonitorsManager::~MonitorsManager()
{
  
}

void MonitorsManager::FillMonitors(const unique_ptr<Helix> &fittedHelix,
                                   const unique_ptr<Helix> &trueHelix,
                                   const unique_ptr<Track> &track)
{
  fractionMonitors["successVsX"].second->Fill(fabs(trueHelix->GetOrigin()->GetX()));
  fractionMonitors["fullSuccessVsX"].second->Fill(fabs(trueHelix->GetOrigin()->GetX()));
  fractionMonitors["successVsY"].second->Fill(fabs(trueHelix->GetOrigin()->GetY()));
  fractionMonitors["fullSuccessVsY"].second->Fill(fabs(trueHelix->GetOrigin()->GetY()));
  fractionMonitors["successVsZ"].second->Fill(fabs(trueHelix->GetOrigin()->GetZ()));
  fractionMonitors["fullSuccessVsZ"].second->Fill(fabs(trueHelix->GetOrigin()->GetZ()));
  fractionMonitors["successVsPx"].second->Fill(fabs(trueHelix->GetMomentum()->GetX()));
  fractionMonitors["fullSuccessVsPx"].second->Fill(fabs(trueHelix->GetMomentum()->GetX()));
  fractionMonitors["successVsPy"].second->Fill(fabs(trueHelix->GetMomentum()->GetY()));
  fractionMonitors["fullSuccessVsPy"].second->Fill(fabs(trueHelix->GetMomentum()->GetY()));
  fractionMonitors["successVsPz"].second->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  fractionMonitors["fullSuccessVsPz"].second->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  fractionMonitors["successVsEta"].second->Fill(fabs(track->GetEta()));
  fractionMonitors["fullSuccessVsEta"].second->Fill(fabs(track->GetEta()));
  
  if(!fittedHelix){
    monitors1D["failReason"]->Fill(8);
    return;
  }
  
  // Here (full) success part starts:
  fractionMonitors["successVsX"].first->Fill(fabs(trueHelix->GetOrigin()->GetX()));
  fractionMonitors["successVsY"].first->Fill(fabs(trueHelix->GetOrigin()->GetY()));
  fractionMonitors["successVsZ"].first->Fill(fabs(trueHelix->GetOrigin()->GetZ()));
  fractionMonitors["successVsPx"].first->Fill(fabs(trueHelix->GetMomentum()->GetX()));
  fractionMonitors["successVsPy"].first->Fill(fabs(trueHelix->GetMomentum()->GetY()));
  fractionMonitors["successVsPz"].first->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
  fractionMonitors["successVsEta"].first->Fill(fabs(track->GetEta()));
  
  
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
    fractionMonitors["fullSuccessVsX"].first->Fill(fabs(trueHelix->GetOrigin()->GetX()));
    fractionMonitors["fullSuccessVsY"].first->Fill(fabs(trueHelix->GetOrigin()->GetY()));
    fractionMonitors["fullSuccessVsZ"].first->Fill(fabs(trueHelix->GetOrigin()->GetZ()));
    fractionMonitors["fullSuccessVsPx"].first->Fill(fabs(trueHelix->GetMomentum()->GetX()));
    fractionMonitors["fullSuccessVsPy"].first->Fill(fabs(trueHelix->GetMomentum()->GetY()));
    fractionMonitors["fullSuccessVsPz"].first->Fill(fabs(trueHelix->GetMomentum()->GetZ()));
    fractionMonitors["fullSuccessVsEta"].first->Fill(fabs(track->GetEta()));
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
  c2->Divide(4,4);
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
