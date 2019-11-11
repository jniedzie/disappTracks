//  drawHelixPlots.cpp
//  Created by Jeremi Niedziela on 11/11/2019.

#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "EventProcessor.hpp"
#include "EventSet.hpp"
#include "HelixProcessor.hpp"

string configPath = "configs/taggerPlotting.md";
string cutLevel = "afterHelixTagging";
string suffix = "";

xtracks::EDataType dataType = xtracks::kSignal;

double getPz(double R0, double s0, double q)
{
  double c = cos(TMath::Pi()/2-atan2(s0, q*R0));
  return c/sqrt(c*c-1)*3*R0*solenoidField/10;
}

void prepareHist(TH1D *hist, int color, int rebin)
{
  hist->Sumw2();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(0.4);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  
  hist->Rebin(rebin);
  hist->Scale(1./(rebin*hist->GetEntries()));
}

void FillMonitors(const EventSet &events)
{
  TH2D *vertexX = new TH2D("pionVertexX", "pionVertexX", 100, -1000, 1000, 100, -1000, 1000);
  TH2D *vertexY = new TH2D("pionVertexY", "pionVertexY", 100, -1000, 1000, 100, -1000, 1000);
  TH2D *vertexZ = new TH2D("pionVertexZ", "pionVertexZ", 100, -1000, 1000, 100, -1000, 1000);
  
  TH2D *vertexXbest = new TH2D("pionVertexXbest", "pionVertexXbest", 100, -1000, 1000, 100, -1000, 1000);
  TH2D *vertexYbest = new TH2D("pionVertexYbest", "pionVertexYbest", 100, -1000, 1000, 100, -1000, 1000);
  TH2D *vertexZbest = new TH2D("pionVertexZbest", "pionVertexZbest", 100, -1000, 1000, 100, -1000, 1000);
  
  TH2D *ptHist      = new TH2D("ptHist"     ,  "ptHist"   , 100, 0, 1000, 100, 0, 1000);
  TH2D *ptHistBest  = new TH2D("ptHistBest" , "ptHistBest", 100, 0, 1000, 100, 0, 1000);
  
  TH2D *pzHist      = new TH2D("pzHist"     ,  "pzHist"   , 100, -1000, 1000, 100, -1000, 1000);
  TH2D *pzHistBest  = new TH2D("pzHistBest" , "pzHistBest", 100, -1000, 1000, 100, -1000, 1000);
  
  TH1D *r0      = new TH1D("r0"     , "r0"    , 100, 0, 1000);
  TH1D *r0best  = new TH1D("r0best" , "r0best", 100, 0, 1000);
  
  TH1D *ar      = new TH1D("ar"     , "ar"    , 100, 0, 500);
  TH1D *arBest  = new TH1D("arBest" , "arBest", 100, 0, 500);
  
  TH1D *s0      = new TH1D("s0"     , "s0"    , 100, 0, 1000);
  TH1D *s0best  = new TH1D("s0best" , "s0best", 100, 0, 1000);
  
  TH1D *bs      = new TH1D("bs"     , "bs"    , 100, -1000, 0);
  TH1D *bsBest  = new TH1D("bsBest" , "bsBest", 100, -1000, 0);
  
  TH1D *vertexXRes      = new TH1D("vertexXRes"     , "vertexXRes"    , 40, -6, 6);
  TH1D *vertexXResBest  = new TH1D("vertexXResBest" , "vertexXResBest", 40, -6, 6);
  
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    
    for(int year : years){
      if(!config.params["load_"+to_string(year)]) continue;
      
      for(int iEvent=0; iEvent<events.size(dataType, iSig, year); iEvent++){
        auto event = events.At(dataType, iSig, year, iEvent);
        if(!event->WasTagged()) continue;
        
        for(int iHelix=0; iHelix<event->GetNhelices(); iHelix++){
          auto helix = event->GetHelix(iHelix);
          
          double radius = helix.GetRadius(helix.GetTmin());
          double slope = helix.GetSlope(helix.GetTmin());
          double pt = GetPtInMagField(radius, solenoidField);
          double pz = getPz(radius, slope, helix.GetCharge());
          double vx = helix.GetVertex()->GetX();
          double vy = helix.GetVertex()->GetY();
          double vz = helix.GetVertex()->GetZ();
          
          r0->Fill(radius);
          ar->Fill(helix.GetRadiusFactor());
          s0->Fill(slope);
          bs->Fill(helix.GetSlopeFactor());
        }
        
        for(int iPion=0; iPion<event->GetGenPionHelices().size(); iPion++){
          auto pionHelix = event->GetGenPionHelices()[iPion];
          
          double bestHelixVertexX=inf;
          double bestHelixVertexY=inf;
          double bestHelixVertexZ=inf;
          double bestPt=inf, bestPz=inf;
          double bestR0=inf, bestAr=inf, bestS0=inf, bestBs=inf;
          double bestVertexXres=inf;
          
          double pionPt = sqrt(pow(pionHelix.GetMomentum().GetX(), 2) +
                               pow(pionHelix.GetMomentum().GetY(), 2));
          double pionPz = pionHelix.GetMomentum().GetZ();
          
          int maxNhitsOnHelix = -inf;
          
          for(int iHelix=0; iHelix<event->GetNhelices(); iHelix++){
            auto helix = event->GetHelix(iHelix);
            
            double radius = helix.GetRadius(helix.GetTmin());
            double slope = helix.GetSlope(helix.GetTmin());
            double pt = GetPtInMagField(radius, solenoidField);
            double pz = getPz(radius, slope, helix.GetCharge());
            double vx = helix.GetVertex()->GetX();
            double vy = helix.GetVertex()->GetY();
            double vz = helix.GetVertex()->GetZ();
            
            vertexX->Fill(pionHelix.GetVertex()->GetX(), vx);
            vertexY->Fill(pionHelix.GetVertex()->GetY(), vy);
            vertexZ->Fill(pionHelix.GetVertex()->GetZ(), vz);
            ptHist->Fill(pionPt, pt);
            pzHist->Fill(pionPz, pz);
   
            vertexXRes->Fill((pionHelix.GetVertex()->GetX()-vx)/pionHelix.GetVertex()->GetX());
            
            if((int)helix.GetPoints().size() > maxNhitsOnHelix){
              maxNhitsOnHelix = (int)helix.GetPoints().size();
              
              bestHelixVertexX = vx;
              bestHelixVertexY = vy;
              bestHelixVertexZ = vz;
              bestPt = pt;
              bestPz = pz;
              bestR0 = helix.GetRadius(helix.GetTmin());
              bestAr = helix.GetRadiusFactor();
              bestS0 = helix.GetSlope(helix.GetTmin());
              bestBs = helix.GetSlopeFactor();
              bestVertexXres = (pionHelix.GetVertex()->GetX()-vx)/pionHelix.GetVertex()->GetX();
            }
          }
          
          vertexXbest->Fill(pionHelix.GetVertex()->GetX(), bestHelixVertexX);
          vertexYbest->Fill(pionHelix.GetVertex()->GetY(), bestHelixVertexY);
          vertexZbest->Fill(pionHelix.GetVertex()->GetZ(), bestHelixVertexZ);
          ptHistBest->Fill(pionPt, bestPt);
          pzHistBest->Fill(pionPz, bestPz);
          r0best->Fill(bestR0);
          arBest->Fill(bestAr);
          s0best->Fill(bestS0);
          bsBest->Fill(bestBs);
          vertexXResBest->Fill(bestVertexXres);
        }
      }
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",2880,1800);
  c1->Divide(4,4);
  
  c1->cd(1); vertexX->Draw("colz");
  c1->cd(2); vertexXbest->Draw("colz");
  
  c1->cd(3); vertexY->Draw("colz");
  c1->cd(4); vertexYbest->Draw("colz");
  
  c1->cd(5); vertexZ->Draw("colz");
  c1->cd(6); vertexZbest->Draw("colz");
  
  c1->cd(7);
  prepareHist(r0, kBlack, 1);
  r0->Draw();
  prepareHist(r0best, kMagenta, 5);
  r0best->Draw("same");
  
  c1->cd(8);
  prepareHist(ar, kBlack, 1);
  ar->DrawNormalized();
  prepareHist(arBest, kMagenta, 5);
  arBest->DrawNormalized("same");
  
  c1->cd(9);
  prepareHist(s0, kBlack, 1);
  s0->DrawNormalized();
  prepareHist(s0best, kMagenta, 5);
  s0best->DrawNormalized("same");
  
  c1->cd(10);
  prepareHist(bs, kBlack, 1);
  bs->DrawNormalized();
  prepareHist(bsBest, kMagenta, 5);
  bsBest->DrawNormalized("same");
  
  c1->cd(11); ptHist->Draw("colz");
  c1->cd(12); ptHistBest->Draw("colz");
  
  c1->cd(13); pzHist->Draw("colz");
  c1->cd(14); pzHistBest->Draw("colz");

  c1->cd(15);
  prepareHist(vertexXRes, kBlack, 1);
  vertexXRes->Fit("gaus");
  vertexXRes->Draw();
  
  prepareHist(vertexXResBest, kMagenta, 1);
  vertexXResBest->Draw("same");
  
  c1->Update();
}



/// The program execution starting point.
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events; events.LoadEventsFromFiles(cutLevel+suffix+"/");
  FillMonitors(events);

  theApp.Run();
  return 0;
}

