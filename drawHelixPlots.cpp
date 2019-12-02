//  drawHelixPlots.cpp
//  Created by Jeremi Niedziela on 11/11/2019.

#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "EventProcessor.hpp"
#include "EventSet.hpp"
#include "HelixProcessor.hpp"

string configPath = "configs/taggerPlotting.md";
string cutLevel = "afterHelixTagging";
string suffix = "_noHighPtHits";

xtracks::EDataType dataType = xtracks::kSignal;

TH2D *vertexX, *vertexY, *vertexZ, *vertexXbest, *vertexYbest, *vertexZbest;
TH2D *ptHist, *ptHistBest, *pzHist, *pzHistBest;

TH1D *r0, *r0best, *ar, *arBest, *s0, *s0best, *bs, *bsBest ;

TH1D *vertexXRes, *vertexXResBest, *vertexYRes, *vertexYResBest, *vertexZRes, *vertexZResBest, *momentumTRes, *momentumTResBest, *momentumZRes, *momentumZResBest, *fractionTrueHitsHist, *fractionTrueHitsHistBest;

TH1D *pionPtHist = new TH1D("pionPt", "pionPt", 100, 0, 1000);

/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix = "";
   
  if(config.secondaryCategory == "Zmumu") prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu") prefix += "Wmunu/";
  
  if(config.params["cuts_level"]==0) prefix += "after_L0/";
  if(config.params["cuts_level"]==1) prefix += "after_L1/"+config.category+"/";
  
  prefix += "afterHelixTagging"+suffix+"/";
  
  return prefix;
}

double getPz(double R0, double s0, double q)
{
  double c = cos(TMath::Pi()/2-atan2(s0, q*R0));
  return c/sqrt(c*c-1)*3*R0*solenoidField/10;
}

void fillSingleHists(const Helix &helix, bool best)
{
  double radius = helix.GetRadius(helix.GetTmin());
  double slope = helix.GetSlope(helix.GetTmin());
  double trueHitsFraction = helix.GetNrecPionHits()/(double)(helix.GetNrecHits()-1);
  
  if(best){
    fractionTrueHitsHistBest->Fill(trueHitsFraction);
    r0best->Fill(radius);
    arBest->Fill(helix.GetRadiusFactor());
    s0best->Fill(slope);
  }
  
  fractionTrueHitsHist->Fill(trueHitsFraction);
  r0->Fill(radius);
  ar->Fill(helix.GetRadiusFactor());
  s0->Fill(slope);
  bs->Fill(helix.GetSlopeFactor());
  
}

void fillDoubleHists(const Helix &recHelix, const Helix &genHelix, bool best)
{
  if(recHelix.GetRadius(recHelix.GetTmin()) > 800 ||
     recHelix.GetRadius(recHelix.GetTmin()) < 50) return;
  
  double genPt = sqrt(pow(genHelix.GetMomentum().GetX(), 2) + pow(genHelix.GetMomentum().GetY(), 2));
  double genPz = genHelix.GetMomentum().GetZ();
  double genVx = genHelix.GetVertex()->GetX();
  double genVy = genHelix.GetVertex()->GetY();
  double genVz = genHelix.GetVertex()->GetZ();
  
  
  double recSlope = recHelix.GetSlope(recHelix.GetTmin());
  
  double recPt = GetPtInMagField(recHelix.GetRadius(recHelix.GetTmin()), solenoidField);
//  double recPt = GetPtInMagField(recHelix.GetRadius(-recHelix.GetTmin()), solenoidField);
  
  double recPz = getPz(recHelix.GetRadius(recHelix.GetTmin()), recSlope, recHelix.GetCharge());
  double recVx = recHelix.GetVertex()->GetX();
  double recVy = recHelix.GetVertex()->GetY();
  double recVz = recHelix.GetVertex()->GetZ();
  
  
  if(best){
    vertexXbest->Fill(genVx, recVx);
    vertexYbest->Fill(genVy, recVy);
    vertexZbest->Fill(genVz, recVz);
    ptHistBest->Fill(genPt, recPt);
    pzHistBest->Fill(genPz, recPz);
    
    vertexXResBest->Fill((recVx-genVx));
    vertexYResBest->Fill((recVy-genVy));
    vertexZResBest->Fill((recVz-genVz));
    momentumTResBest->Fill((recPt-genPt)/genPt);
    momentumZResBest->Fill((recPz-genPz)/genPz);
  }
  else{
    vertexX->Fill(genVx, recVx);
    vertexY->Fill(genVy, recVy);
    vertexZ->Fill(genVz, recVz);
    ptHist->Fill(genPt, recPt);
    pzHist->Fill(genPz, recPz);
    
    vertexXRes->Fill((recVx-genVx));
    vertexYRes->Fill((recVy-genVy));
    vertexZRes->Fill((recVz-genVz));
    momentumTRes->Fill((recPt-genPt)/genPt);
    momentumZRes->Fill((recPz-genPz)/genPz);
  }
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
  vertexX = new TH2D("pionVertexX", "pionVertexX", 100, -1000, 1000, 100, -1000, 1000);
  vertexY = new TH2D("pionVertexY", "pionVertexY", 100, -1000, 1000, 100, -1000, 1000);
  vertexZ = new TH2D("pionVertexZ", "pionVertexZ", 100, -1000, 1000, 100, -1000, 1000);
  
  vertexXbest = new TH2D("pionVertexXbest", "pionVertexXbest", 100, -1000, 1000, 100, -1000, 1000);
  vertexYbest = new TH2D("pionVertexYbest", "pionVertexYbest", 100, -1000, 1000, 100, -1000, 1000);
  vertexZbest = new TH2D("pionVertexZbest", "pionVertexZbest", 100, -1000, 1000, 100, -1000, 1000);
  
  ptHist      = new TH2D("ptHist"     ,  "ptHist"   , 100, 0, 1000, 100, 0, 1000);
  ptHistBest  = new TH2D("ptHistBest" , "ptHistBest", 100, 0, 1000, 100, 0, 1000);
  
  pzHist      = new TH2D("pzHist"     ,  "pzHist"   , 100, -1000, 1000, 100, -1000, 1000);
  pzHistBest  = new TH2D("pzHistBest" , "pzHistBest", 100, -1000, 1000, 100, -1000, 1000);
  
  r0      = new TH1D("r0"     , "r0"    , 100, 0, 1000);
  r0best  = new TH1D("r0best" , "r0best", 100, 0, 1000);
  
  ar      = new TH1D("ar"     , "ar"    , 100, 0, 500);
  arBest  = new TH1D("arBest" , "arBest", 100, 0, 500);
  
  s0      = new TH1D("s0"     , "s0"    , 100, 0, 1000);
  s0best  = new TH1D("s0best" , "s0best", 100, 0, 1000);
  
  bs      = new TH1D("bs"     , "bs"    , 100, -1000, 0);
  bsBest  = new TH1D("bsBest" , "bsBest", 100, -1000, 0);
  
  vertexXRes      = new TH1D("vertexXRes"     , "vertexXRes"    , 50, -50, 50);
  vertexXResBest  = new TH1D("vertexXResBest" , "vertexXResBest", 50, -50, 50);
  
  vertexYRes      = new TH1D("vertexYRes"     , "vertexYRes"    , 50, -50, 50);
  vertexYResBest  = new TH1D("vertexYResBest" , "vertexYResBest", 50, -50, 50);
  
  vertexZRes      = new TH1D("vertexZRes"     , "vertexZRes"    , 50, -50, 50);
  vertexZResBest  = new TH1D("vertexZResBest" , "vertexZResBest", 50, -50, 50);
  
  momentumTRes      = new TH1D("momentumTRes"     , "momentumTRes"    , 60, -10, 10);
  momentumTResBest  = new TH1D("momentumTResBest" , "momentumTResBest", 60, -10, 10);
  
  momentumZRes      = new TH1D("momentumZRes"     , "momentumZRes"    , 40, -60, 60);
  momentumZResBest  = new TH1D("momentumZResBest" , "momentumZResBest", 40, -60, 60);
  
  pionPtHist = new TH1D("pionPt", "pionPt", 100, 0, 1000);
  
  fractionTrueHitsHist      = new TH1D("fractionTrueHitsHist", "fractionTrueHitsHist", 100, 0, 1);
  fractionTrueHitsHistBest  = new TH1D("fractionTrueHitsHistBest", "fractionTrueHitsHistBest", 100, 0, 1);
  
  map<double, int> nEventsWithTrueHelix;
  map<double, int> nEventsWithTrueHelix2;
  int nEvents = 0;
  vector<double> trueHelixThresholds = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    
    for(int year : years){
      if(!config.params["load_"+to_string(year)]) continue;
      
      for(int iEvent=0; iEvent<events.size(dataType, iSig, year); iEvent++){
        auto event = events.At(dataType, iSig, year, iEvent);
        if(!event->WasTagged()) continue;
      
        map<double, bool> trueHelixFound;
        map<double, bool> trueHelixFound2;
        
        
        Helix bestHelix(Point(), Point(), 0);
        
        int maxNhitsOnHelix = -inf;
        
        int nPionHits = event->GetPionClusters().size();
//        if(nPionHits==0) continue;
        
        nEvents++;
        
        for(int iHelix=0; iHelix<event->GetNhelices(); iHelix++){
          auto helix = event->GetHelix(iHelix);
          
          double trueHitsFraction = helix.GetNrecPionHits()/(double)(helix.GetNrecHits()-1);
          double pionHitsFraction = helix.GetNrecPionHits()/(double)nPionHits;
          
          if(helix.GetNrecHits()==0) trueHitsFraction=0;
          if(nPionHits == 0) pionHitsFraction=0;
          
          for(double threshold : trueHelixThresholds){
            
            if(!trueHelixFound[threshold]){
              if(trueHitsFraction >= threshold){
                nEventsWithTrueHelix[threshold]++;
                trueHelixFound[threshold] = true;
              }
            }
            if(!trueHelixFound2[threshold]){
              if(pionHitsFraction >= threshold){
                nEventsWithTrueHelix2[threshold]++;
                trueHelixFound2[threshold] = true;
              }
            }
          }
          
          fillSingleHists(helix, false);
          
          if(helix.GetNrecHits() > maxNhitsOnHelix){
//          if((int)helix.GetPoints().size() > maxNhitsOnHelix){
            maxNhitsOnHelix = helix.GetNrecHits();
            bestHelix = helix;
          }
        }
        
        fillSingleHists(bestHelix, true);
        
        for(int iPion=0; iPion<event->GetGenPionHelices().size(); iPion++){
          auto pionHelix = event->GetGenPionHelices()[iPion];
          fillDoubleHists(bestHelix, pionHelix, true);
          
          for(int iHelix=0; iHelix<event->GetNhelices(); iHelix++){
            auto helix = event->GetHelix(iHelix);
            fillDoubleHists(helix, pionHelix, false);
          }
        }
      }
      
      cout<<"N events: "<<nEvents<<endl;
      
      for(double threshold : trueHelixThresholds){
//        cout<<"N events with true helix above "<<threshold<<": "<<nEventsWithTrueHelix[threshold]<<endl;
        cout<<"Fraction of events with true helix above "<<threshold<<":\t"<<nEventsWithTrueHelix[threshold]/(double)nEvents<<"\t\t";
    cout<<"+/-\t\t"<<nEventsWithTrueHelix[threshold]/(double)nEvents*sqrt(1./nEventsWithTrueHelix[threshold]+1/nEvents)<<endl;
        
        cout<<"Fraction of events with true helix above "<<threshold<<" (second method):\t"<<nEventsWithTrueHelix2[threshold]/(double)nEvents<<"\t\t";
        cout<<"+\-\t\t"<<nEventsWithTrueHelix2[threshold]/(double)nEvents*sqrt(1./nEventsWithTrueHelix2[threshold]+1/nEvents)<<endl;
        
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
  
  c1->cd(15); fractionTrueHitsHist->Draw();
  c1->cd(16); fractionTrueHitsHistBest->Draw();

  
  TCanvas *c2 = new TCanvas("c2","c2",2880,1800);
  c2->Divide(2,3);
  
  c2->cd(1);
  prepareHist(vertexXRes, kBlack, 1);
  vertexXRes->Fit("gaus");
  vertexXRes->Draw();
  prepareHist(vertexXResBest, kMagenta, 1);
  vertexXResBest->Draw("same");
  
  c2->cd(2);
  prepareHist(vertexYRes, kBlack, 1);
  vertexYRes->Fit("gaus");
  vertexYRes->Draw();
  prepareHist(vertexYResBest, kMagenta, 1);
  vertexYResBest->Draw("same");
  
  c2->cd(3);
  prepareHist(vertexZRes, kBlack, 1);
  vertexZRes->Fit("gaus");
  vertexZRes->Draw();
  prepareHist(vertexZResBest, kMagenta, 1);
  vertexZResBest->Draw("same");
  
  c2->cd(4);
  prepareHist(momentumTRes, kBlack, 1);
  momentumTRes->Fit("gaus");
  momentumTRes->Draw();
  prepareHist(momentumTResBest, kMagenta, 1);
  momentumTResBest->Draw("same");
  
  c2->cd(5);
  prepareHist(momentumZRes, kBlack, 1);
  momentumZRes->Fit("gaus");
  momentumZRes->Draw();
  prepareHist(momentumZResBest, kMagenta, 1);
  momentumZResBest->Draw("same");
  
  c2->cd(6);
  pionPtHist->Draw();

  
  
  c1->Update();
  c2->Update();
}



/// The program execution starting point.
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  FillMonitors(events);

  theApp.Run();
  return 0;
}

