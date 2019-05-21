//  analyzeClusters.cpp
//
//  Created by Jeremi Niedziela on 10/05/2019.

#include "Helpers.hpp"
#include "EventSet.hpp"

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L1/3layers/";//after_L1/";
string outfileName = "genPionHists3layers";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 3;


vector<tuple<string, int, double, double, string>> histOptions1D = {
// title                        nBins  min     max    titleX
  {"lumi"                       , 300, 0      , 300   , "lumi block"              },
  {"pion_px"                    , 50 , 0      , 1000  , "|p_{x}| (MeV)"           },
  {"pion_py"                    , 50 , 0      , 1000  , "|p_{y}| (MeV)"           },
  {"pion_pz"                    , 50 , 0      , 2000  , "|p_{z}| (MeV)"           },
  {"pion_pt"                    , 50 , 0      , 1000  , " p_{t}  (MeV)"           },
  {"pion_n_clusters"            , 50 , 0      , 100   , "# clusters"              },
  {"pion_n_clusters_TIB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TOB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TID"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TEC"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_PXB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_PXE"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_initial_radius"        , 50 , 0      , 1000  , "R_{max} (mm)"            },
  {"pion_final_radius"          , 50 , 0      , 1000  , "R_{min} (mm)"            },
  {"pion_range_z"               , 100, 0      , 3500  , "#Delta z (mm)"           },
  {"pion_vertex_xy"             , 100, 0      , 500   , "v_{xy} (mm)"             },
  {"pion_vertex_z"              , 50 , -1000  , 1000  , "v_{z} (mm)"              },
  {"delta_phi_pion_chargino"    , 50 , 0      , 3.2   , "#Delta #phi (#pi, #chi)" },
  {"next_point_delta_phi"       , 100, 0      , 3.2   , "#Delta #phi"             },
  {"next_point_delta_z"         , 100, 0      , 500   , "#Delta z"                },
  {"middle_seed_hit_delta_phi"  , 32 , 0      , 3.2   , "#Delta #phi"             },
  {"middle_seed_hit_delta_z"    , 50 , 0      , 500   , "#Delta z"                },
  {"last_seed_hit_delta_phi"    , 32 , 0      , 3.2   , "#Delta #phi"             },
  {"last_seed_hit_delta_z"      , 50 , 0      , 500   , "#Delta z"                },
  {"pion_cluster_charge_TIB"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_TOB"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_TID"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_TEC"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_PXB"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_PXE"    , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TIB" , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TOB" , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TID" , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TEC" , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_PXB" , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_PXE" , 100, 0      , 600   , ""                        },
};

vector<tuple<string, int, double, double, string, int, double, double, string>> histOptions2D = {
  // title                        nBinsX  minX maxX    titleX nBinsY  minY maxY    titleY
  {"next_point_delta_phi_pion_pt" , 50, 0, 3.14 , "#Delta #phi" , 50, 0, 1000, "pion p_{t} (MeV)" },
  {"next_point_delta_z_pion_pz"   , 50, 0, 1000 , "#Delta z"    , 50, 0, 1000, "pion p_{z} (MeV)" },
};
map<string, TH1D *> hists1D;
map<string, TH2D *> hists2D;

void SetupHists()
{
  for(auto &[title, nBins, min, max, titleX] : histOptions1D){
    hists1D[title] = new TH1D(title.c_str(), title.c_str(), nBins, min, max);
    hists1D[title]->GetXaxis()->SetTitle(titleX.c_str());
  }
  
  for(auto &[title, nBinsX, minX, maxX, titleX, nBinsY, minY, maxY, titleY] : histOptions2D){
    hists2D[title] = new TH2D(title.c_str(), title.c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    hists2D[title]->GetXaxis()->SetTitle(titleX.c_str());
    hists2D[title]->GetYaxis()->SetTitle(titleY.c_str());
  }
}

range<double> GetPointsRingSize(vector<shared_ptr<Point>> points)
{
  sort(points.begin(), points.end(), PointsProcessor::ComparePointByX());
  double minX = points.front()->GetX();
  double maxX = points.back()->GetX();

  sort(points.begin(), points.end(), PointsProcessor::ComparePointByY());
  double minY = points.front()->GetY();
  double maxY = points.back()->GetY();

  Point center((minX+maxX)/2., (minY+maxY)/2., 0);
  double minDistanceToCenter = inf;
  double maxDistanceToCenter = -inf;

  for(auto &point : points){
    double dist = pointsProcessor.distanceXY(center, *point);
    if(dist < minDistanceToCenter) minDistanceToCenter = dist;
    if(dist > maxDistanceToCenter) maxDistanceToCenter = dist;
  }
  
  return range<double>(minDistanceToCenter, maxDistanceToCenter);
}
  
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  
  EventSet events; events.LoadEventsFromFiles(cutLevel);
  SetupHists();
  
  map<string, int> nHitsInDet;
  map<string, int> nLastHitsInDet;
  int nPionHits = 0;
  
  cout<<"Filling histograms"<<endl;
  
  for(int iEvent=0; iEvent<events.size(dataType, setIter); iEvent++){
    auto event = events.At(dataType, setIter, iEvent);
    
    if(!event){
      cout<<"Event not found"<<endl;
      exit(0);
    }
    
    // Pick only events with one track and one gen pion and resonable number of pion clusters
    if(event->GetGenPionHelices().size() != 1) continue;
    if(event->GetNtracks() != 1) continue;
    if(event->GetPionClusters().size() < 3) continue;
    if(event->GetPionSimHits().size() < 3) continue;
    
    // Get what's needed from the event
    auto pionHelix       = event->GetGenPionHelices()[0];
    auto track           = event->GetTrack(0);
    auto eventVertex     = event->GetVertex();
    auto pionSimHits     = event->GetPionSimHits();
    auto pionClusters    = event->GetPionClusters();
    auto trackerClusters = event->GetTrackerClusters();
    
    // Fill basic info about events
    hists1D["lumi"]->Fill(event->GetLumiSection());
    
    // Fill gen-pion histograms
    double pionPt = sqrt(pow(pionHelix.GetMomentum()->GetX(), 2) + pow(pionHelix.GetMomentum()->GetY(), 2));
    double pionPz = pionHelix.GetMomentum()->GetZ();
    double phiPion = atan2(pionHelix.GetMomentum()->GetY(), pionHelix.GetMomentum()->GetX());
    
    hists1D["pion_px"]->Fill(fabs(pionHelix.GetMomentum()->GetX()));
    hists1D["pion_py"]->Fill(fabs(pionHelix.GetMomentum()->GetY()));
    hists1D["pion_pz"]->Fill(fabs(pionPz));
    hists1D["pion_pt"]->Fill(pionPt);
    hists1D["pion_vertex_z"]->Fill(pionHelix.GetVertex()->GetZ());
    hists1D["pion_vertex_xy"]->Fill(sqrt(pow(pionHelix.GetVertex()->GetX(), 2) + pow(pionHelix.GetVertex()->GetY(), 2)));
    
    for(auto &track : event->GetTracks()){
      hists1D["delta_phi_pion_chargino"]->Fill(fabs(phiPion - track->GetPhi()));
    }
    
    // Fill pion sim hit histograms
    sort(pionSimHits.begin(), pionSimHits.end(), PointsProcessor::ComparePointByZ());
    
    for(int iHit=1; iHit<pionSimHits.size()-1; iHit++){
      double deltaPhi = pointsProcessor.GetPointingAngleXY(*pionSimHits[iHit-1],
                                                           *pionSimHits[iHit],
                                                           *pionSimHits[iHit+1]);
      
      double deltaZ = fabs(pionSimHits[iHit]->GetZ() - pionSimHits[iHit+1]->GetZ());
      
      hists1D["next_point_delta_phi"]->Fill(deltaPhi);
      hists2D["next_point_delta_phi_pion_pt"]->Fill(deltaPhi, pionPt);
      hists1D["next_point_delta_z"]->Fill(deltaZ);
      hists2D["next_point_delta_z_pion_pz"]->Fill(deltaZ, pionPz);
    }
    
    shared_ptr<Point> farthestPoint;
    if(fabs(pionSimHits.front()->GetZ()) > fabs(pionSimHits.back()->GetZ()))  farthestPoint = pionSimHits.front();
    else                                                                        farthestPoint = pionSimHits.back();
    
    hists1D["pion_range_z"]->Fill(fabs(eventVertex->GetZ() - farthestPoint->GetZ()));
    
    range<double> pionRingSize = GetPointsRingSize(pionSimHits);
    
    hists1D["pion_initial_radius"]->Fill(pionRingSize.GetMax());
    hists1D["pion_final_radius"]->Fill(pionRingSize.GetMin());
    
    double maxZ = -999;
    string maxDet;
    map<string, int> nClustersInDet;
    
    // Fill pion rec clusters histograms
    for(auto &cluster : pionClusters)
    {
      string detName = cluster->GetSubDetName();
      string title = "pion_cluster_charge_" + detName;
      if(hists1D.find(title) != hists1D.end()){
        hists1D[title]->Fill(cluster->GetValue());
      }
      
      nClustersInDet[detName]++;
      nHitsInDet[detName]++;
      
      if(fabs(cluster->GetZ()) > maxZ){
        maxZ = fabs(cluster->GetZ());
        maxDet = detName;
      }
    }
    
    nPionHits += pionSimHits.size();
    hists1D["pion_n_clusters"]->Fill(pionClusters.size());
    hists1D["pion_n_clusters_TIB"]->Fill(nClustersInDet["TIB"]);
    hists1D["pion_n_clusters_TOB"]->Fill(nClustersInDet["TOB"]);
    hists1D["pion_n_clusters_TID"]->Fill(nClustersInDet["TID"]);
    hists1D["pion_n_clusters_TEC"]->Fill(nClustersInDet["TEC"]);
    hists1D["pion_n_clusters_PXB"]->Fill(nClustersInDet["PXB"]);
    hists1D["pion_n_clusters_PXE"]->Fill(nClustersInDet["PXE"]);
    
    
    // Fill seeds tracking histograms
    Point decayVertexMin = pointsProcessor.GetPointOnTrack(layerR[track->GetNtrackerLayers()-1], *track, *eventVertex);
    Point decayVertexMax = pointsProcessor.GetPointOnTrack(layerR[track->GetNtrackerLayers()],   *track, *eventVertex);
    Point middleSeedHit  = *pionSimHits[0];
    Point lastSeedHit    = *pionSimHits[1];
    
    double middleHitDeltaPhi = min( pointsProcessor.GetPointingAngleXY(Point(0,0,0), decayVertexMin, middleSeedHit),
                                    pointsProcessor.GetPointingAngleXY(Point(0,0,0), decayVertexMax, middleSeedHit));
    
    double middleHitDeltaZ = min( fabs(middleSeedHit.GetZ() - decayVertexMin.GetZ()),
                                  fabs(middleSeedHit.GetZ() - decayVertexMax.GetZ()));
    
    double lastHitDeltaPhi = min( pointsProcessor.GetPointingAngleXY(decayVertexMin, middleSeedHit, lastSeedHit),
                                  pointsProcessor.GetPointingAngleXY(decayVertexMax, middleSeedHit, lastSeedHit));
    
    
    double lastPointDeltaZ = fabs(middleSeedHit.GetZ() - lastSeedHit.GetZ());
    
    
    hists1D["middle_seed_hit_delta_phi"]->Fill(middleHitDeltaPhi);
    hists1D["middle_seed_hit_delta_z"]->Fill(middleHitDeltaZ);
    hists1D["last_seed_hit_delta_phi"]->Fill(lastHitDeltaPhi);
    hists1D["last_seed_hit_delta_z"]->Fill(lastPointDeltaZ);
  
    // Fill noise clusters histograms
    for(auto &cluster : trackerClusters){
      string title = "tracker_cluster_charge_" + cluster->GetSubDetName();
      
      if(hists1D.find(title) != hists1D.end()){
        hists1D[title]->Fill(cluster->GetValue());
      }
    }
    
    nLastHitsInDet[maxDet]++;
  }
  
  TCanvas *genPionCanvas = new TCanvas("gen_pion", "gen_pion", 2880, 1800);
  genPionCanvas->Divide(4,3);
  
  genPionCanvas->cd(1);   hists1D["lumi"]->Draw();
  genPionCanvas->cd(2);   hists1D["pion_pz"]->Draw();
  genPionCanvas->cd(3);   hists1D["pion_px"]->Draw();
  genPionCanvas->cd(4);   hists1D["pion_py"]->Draw();
  genPionCanvas->cd(5);   hists1D["pion_pt"]->Draw();
  genPionCanvas->cd(6);   hists1D["pion_n_clusters"]->Draw();
  genPionCanvas->cd(7);   hists1D["pion_initial_radius"]->Draw();
  genPionCanvas->cd(8);   hists1D["pion_final_radius"]->Draw();
  genPionCanvas->cd(9);   hists1D["pion_range_z"]->Draw();
  genPionCanvas->cd(10);  hists1D["pion_vertex_z"]->Draw();
  genPionCanvas->cd(11);  hists1D["pion_vertex_xy"]->Draw();
  genPionCanvas->cd(12);  hists1D["delta_phi_pion_chargino"]->Draw();
  
  TCanvas *trackingCanvas = new TCanvas("tracking", "tracking", 2880, 1800);
  trackingCanvas->Divide(4,2);
  
  trackingCanvas->cd(1);  hists1D["middle_seed_hit_delta_phi"]->Draw();
  trackingCanvas->cd(2);  hists1D["middle_seed_hit_delta_z"]->Draw();
  trackingCanvas->cd(3);  hists1D["last_seed_hit_delta_phi"]->Draw();
  trackingCanvas->cd(4);  hists1D["last_seed_hit_delta_z"]->Draw();
  trackingCanvas->cd(5);  hists1D["next_point_delta_phi"]->Draw();
  trackingCanvas->cd(6);  hists2D["next_point_delta_phi_pion_pt"]->Draw("colz");
  trackingCanvas->cd(7);  hists1D["next_point_delta_z"]->Draw();
  trackingCanvas->cd(8);  hists2D["next_point_delta_z_pion_pz"]->Draw("colz");
  
  TCanvas *chargeCanvas   = new TCanvas("charge", "charge", 2880, 1800);
  TCanvas *clustersCanvas = new TCanvas("clusters", "clusters", 2880, 1800);
  
  chargeCanvas->Divide(2,3);
  clustersCanvas->Divide(2,3);
  
  vector<string> detNames = {"TIB", "TOB", "TID", "TEC", "PXB", "PXE" };
  int nLastHits=0;
  for(auto &[key, val] : nLastHitsInDet) nLastHits += val;
  
  for(int i=0; i<detNames.size(); i++){
    chargeCanvas->cd(i+1);
    hists1D["pion_cluster_charge_"+detNames[i]]->SetLineColor(kRed);
    hists1D["tracker_cluster_charge_"+detNames[i]]->SetLineColor(kBlue);
    hists1D["pion_cluster_charge_"+detNames[i]]->DrawNormalized();
    hists1D["tracker_cluster_charge_"+detNames[i]]->DrawNormalized("same");
    
    clustersCanvas->cd(i+1);
    hists1D["pion_n_clusters_"+detNames[i]]->Draw();
    
    cout<<"Number of pion hits in "<<detNames[i]<<": "<<nHitsInDet[detNames[i]]<<"\t("<<nHitsInDet[detNames[i]]/(double)nPionHits<<" %)"<<endl;
    cout<<"Number of last hits in "<<detNames[i]<<": "<<nLastHitsInDet[detNames[i]]<<"\t("<<nLastHitsInDet[detNames[i]]/(double)nLastHits<<" %)"<<endl;
  }
  
  genPionCanvas->Update();
  trackingCanvas->Update();
  chargeCanvas->Update();
  clustersCanvas->Update();
  
  TFile *outFile = new TFile(outfileName.c_str(), "recreate");
  outFile->cd();
  
  for(auto &[name, hist] : hists1D) hist->Write();
  for(auto &[name, hist] : hists2D) hist->Write();
  
  outFile->Close();
  
  theApp.Run();
}
