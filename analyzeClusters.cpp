//  analyzeClusters.cpp
//
//  Created by Jeremi Niedziela on 10/05/2019.

#include "Helpers.hpp"
#include "EventSet.hpp"

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L1/3layers/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 3;

struct ComparePointByZ{
  bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
    return (p1->GetZ() < p2->GetZ());
  }
};

struct ComparePointByX{
  bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
    return (p1->GetX() < p2->GetX());
  }
};

struct ComparePointByY{
  bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
    return (p1->GetY() < p2->GetY());
  }
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
 
  config = ConfigManager(configPath);
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  
  
  TH1D *lumiBlockHist = new TH1D("Lumi", "Lumi", 300, 0, 300);
  
  TH1D *pionPxHist    = new TH1D("Pion px", "Pion px", 200, 0, 1000);
  pionPxHist->GetXaxis()->SetTitle("|p_{x}| (MeV)");
  TH1D *pionPyHist    = new TH1D("Pion py", "Pion py", 200, 0, 1000);
  pionPyHist->GetXaxis()->SetTitle("|p_{y}| (MeV)");
  TH1D *pionPzHist    = new TH1D("Pion pz", "Pion pz", 400, 0, 2000);
  pionPzHist->GetXaxis()->SetTitle("|p_{z}| (MeV)");
  TH1D *pionPtHist    = new TH1D("Pion pt", "Pion pt", 100,  0, 1000);
  pionPtHist->GetXaxis()->SetTitle("p_{T} (MeV)");
  
  TH1D *nPionClusters = new TH1D("N pion clusters", "N pion clusters", 100,  0, 100);
  
  TH1D *initialRadius = new TH1D("Initial pion helix radius", "Initial pion helix radius", 100,  0, 1000);
  initialRadius->GetXaxis()->SetTitle("R_{max} (mm)");
  TH1D *finalRadius   = new TH1D("Final pion helix radius", "Final pion helix radius", 100,  0, 1000);
  finalRadius->GetXaxis()->SetTitle("R_{min} (mm)");
  
  TH1D *zRangeHist    = new TH1D("Z range", "Z range", 100,  0, 10000);
  zRangeHist->GetXaxis()->SetTitle("#Delta z (mm)");
  TH1D *pionVertexZ   = new TH1D("Pion vertex Z", "Pion vertex Z", 100,  -1000, 1000);
  pionVertexZ->GetXaxis()->SetTitle("v_{z} (mm)");
  TH1D *pionVertexR   = new TH1D("Pion vertex radius", "Pion vertex radius", 100,  0, 1000);
  pionVertexR->GetXaxis()->SetTitle("v_{xy} (mm)");
  
  TH1D *nPionClustersInDet[20];
  TH1D *pionClusterChargeHist[20];
  TH1D *trackerClusterChargeHist[20];
  
  for(int iDet=0; iDet<20; iDet++){
    pionClusterChargeHist[iDet]   = new TH1D(Form("Pion cluster charge %s", subDetMap[iDet].c_str()),
                                             Form("Pion cluster charge %s", subDetMap[iDet].c_str()),  100,  0, 600);
    
    trackerClusterChargeHist[iDet] = new TH1D(Form("Tracker cluster charge %s", subDetMap[iDet].c_str()),
                                              Form("Tracker cluster charge %s", subDetMap[iDet].c_str()),
                                              100,  0, 600);
    
    nPionClustersInDet[iDet] = new TH1D(Form("N pion clusters %s", subDetMap[iDet].c_str()),
                                        Form("N pion clusters %s", subDetMap[iDet].c_str()),
                                              100,  0, 100);
  }
  
  int nHitsInDet[20] = {0};
  int nLastHitsInDet[20] = {0};
  int nPionHits = 0;
  
  cout<<"Filling histograms"<<endl;
  
  for(int iEvent=0; iEvent<events.size(dataType, setIter); iEvent++){
    auto event = events.At(dataType, setIter, iEvent);
    
    if(!event){
      cout<<"eventDisplay -- event not found"<<endl;
      exit(0);
    }
    
    // Fill basic info about events
    lumiBlockHist->Fill(event->GetLumiSection());
    vector<shared_ptr<Point>> pionClusters = event->GetPionClusters();
    nPionClusters->Fill(pionClusters.size());
    nPionHits += pionClusters.size();
    
    // Fill gen-level info about pions
    for(auto pion : event->GetGenPionHelices()){
      pionPxHist->Fill(fabs(pion.GetMomentum()->GetX()));
      pionPyHist->Fill(fabs(pion.GetMomentum()->GetY()));
      pionPzHist->Fill(fabs(pion.GetMomentum()->GetZ()));
      
      pionPtHist->Fill(sqrt(pow(pion.GetMomentum()->GetX(), 2)
                            +pow(pion.GetMomentum()->GetY(), 2)));
      
      pionVertexZ->Fill(pion.GetVertex()->GetZ());
      pionVertexR->Fill(sqrt(pow(pion.GetVertex()->GetX(), 2)
                             + pow(pion.GetVertex()->GetY(), 2)));
      
    }
    
    // Fill info about pion rec clusters
    if(pionClusters.size() >= 2){
      auto eventVertex = event->GetVertex();
      
      sort(pionClusters.begin(), pionClusters.end(), ComparePointByZ());
      
      shared_ptr<Point> farthestPoint;
      if(fabs(pionClusters.front()->GetZ()) > fabs(pionClusters.back()->GetZ()))  farthestPoint = pionClusters.front();
      else                                                                        farthestPoint = pionClusters.back();
      
      zRangeHist->Fill(fabs(eventVertex->GetZ() - farthestPoint->GetZ()));
      
      sort(pionClusters.begin(), pionClusters.end(), ComparePointByX());
      double minX = pionClusters.front()->GetX();
      double maxX = pionClusters.back()->GetX();
      
      sort(pionClusters.begin(), pionClusters.end(), ComparePointByY());
      double minY = pionClusters.front()->GetY();
      double maxY = pionClusters.back()->GetY();
      
      Point center((minX+maxX)/2., (minY+maxY)/2., 0);
      double minDistanceToCenter = inf;
      double maxDistanceToCenter = -inf;
      
      for(auto &cluster : pionClusters){
        double dist = pointsProcessor.distanceXY(center, *cluster);
        if(dist < minDistanceToCenter) minDistanceToCenter = dist;
        if(dist > maxDistanceToCenter) maxDistanceToCenter = dist;
      }
      
      initialRadius->Fill(maxDistanceToCenter);
      finalRadius->Fill(minDistanceToCenter);
    }
    
    // Fill per-detector info
    double maxZ = -999;
    int maxDet = -1;
    int nClustersInDet[20] = {0};
    
    for(int iDet=0; iDet<20; iDet++){
      for(auto &cluster : pionClusters){
        if(cluster->GetSubDetName() == subDetMap[iDet]){
          pionClusterChargeHist[iDet]->Fill(cluster->GetValue());
          nHitsInDet[iDet]++;
          nClustersInDet[iDet]++;
          
          if(fabs(cluster->GetZ()) > maxZ){
            maxZ = fabs(cluster->GetZ());
            maxDet = iDet;
          }
        }
      }
      
      auto trackerClusters = event->GetTrackerClusters();
      
      for(auto &cluster : trackerClusters){
        if(cluster->GetSubDetName() == subDetMap[iDet]){
          trackerClusterChargeHist[iDet]->Fill(cluster->GetValue());
        }
      }
    }
    
    for(int iDet=0; iDet<20; iDet++){
      nPionClustersInDet[iDet]->Fill(nClustersInDet[iDet]);
    }
    
    nLastHitsInDet[maxDet]++;
  }
  
  TCanvas *c1 = new TCanvas("c1", "c1", 2880, 1800);
  c1->Divide(4,3);
  
  c1->cd(1);
  lumiBlockHist->Draw();
  c1->cd(2);
  pionPzHist->Draw();
  c1->cd(3);
  pionPxHist->Draw();
  c1->cd(4);
  pionPyHist->Draw();
  c1->cd(5);
  pionPtHist->Draw();
  c1->cd(6);
  nPionClusters->Draw();
  c1->cd(7);
  initialRadius->Draw();
  c1->cd(8);
  finalRadius->Draw();
  c1->cd(9);
  zRangeHist->Draw();
  c1->cd(10);
  pionVertexZ->Draw();
  c1->cd(11);
  pionVertexR->Draw();
  
  
  TCanvas *chargeCanvas = new TCanvas("Charge", "Charge", 2880, 1800);
  chargeCanvas->Divide(2,3);
  
  TCanvas *nClustersCanvas = new TCanvas("N clusters", "N clusters", 2880, 1800);
  nClustersCanvas->Divide(2,3);
  
  int iPad=1;
  
  int nLastHits=0;
  
  for(int iDet=0; iDet<20; iDet++){
    nLastHits += nLastHitsInDet[iDet];
  }
  
  for(int iDet=0; iDet<20; iDet++){
    cout<<endl;
    cout<<"Number of pion hits in "<<subDetMap[iDet]<<": "<<nHitsInDet[iDet]<<"\t("<<nHitsInDet[iDet]/(double)nPionHits<<" %)"<<endl;
    cout<<"Number of last hits in "<<subDetMap[iDet]<<": "<<nLastHitsInDet[iDet]<<"\t("<<nLastHitsInDet[iDet]/(double)nLastHits<<" %)"<<endl;
    
    if(subDetMap[iDet] != "TIB" && subDetMap[iDet] != "TOB" &&
       subDetMap[iDet] != "TID" && subDetMap[iDet] != "TEC" &&
       subDetMap[iDet] != "P1PXB" && subDetMap[iDet] != "P1PXEC") continue;
    
    chargeCanvas->cd(iPad);
    pionClusterChargeHist[iDet]->SetLineColor(kRed);
    trackerClusterChargeHist[iDet]->SetLineColor(kBlue);
    
    THStack *stack = new THStack(subDetMap[iDet].c_str(),subDetMap[iDet].c_str());
    
    if(pionClusterChargeHist[iDet]->GetEntries() != 0) pionClusterChargeHist[iDet]->Scale(1/pionClusterChargeHist[iDet]->Integral());
    if(trackerClusterChargeHist[iDet]->GetEntries() != 0) trackerClusterChargeHist[iDet]->Scale(1/trackerClusterChargeHist[iDet]->Integral());
    
    stack->Add(pionClusterChargeHist[iDet]);
    stack->Add(trackerClusterChargeHist[iDet]);
    stack->Draw("nostack");
    
    nClustersCanvas->cd(iPad++);
    nPionClustersInDet[iDet]->Draw();
  }
  
  c1->Update();
  chargeCanvas->Update();
  nClustersCanvas->Update();
  
  theApp.Run();
}
