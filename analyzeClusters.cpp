//  analyzeClusters.cpp
//
//  Created by Jeremi Niedziela on 10/05/2019.

#include "Helpers.hpp"
#include "EventSet.hpp"

string configPath  = "configs/eventDisplay.md";
string cutLevel    = "after_L1/all/";// "after_L1/4layers/"; //"after_L1/4layers/";//after_L1/";
string outfileName = "results/tmp.root";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

vector<tuple<string, int, double, double, string>> histOptions1D = {
// title                        nBins  min     max    titleX
  {"lumi"                       , 300, 0      , 300   , "lumi block"              },
  {"pion_px"                    , 50 , 0      , 1000  , "|p_{x}| (MeV)"           },
  {"pion_py"                    , 50 , 0      , 1000  , "|p_{y}| (MeV)"           },
  {"pion_pz"                    , 50 , 0      , 1500  , "|p_{z}| (MeV)"           },
  {"pion_pt"                    , 50 , 0      , 1000  , " p_{t}  (MeV)"           },
  {"pion_n_clusters"            , 50 , 0      , 100   , "# clusters"              },
  {"pion_n_clusters_TIB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TOB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TID"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TEC"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_P1PXB"      , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_P1PXE"      , 50 , 0      , 50    , "# clusters"              },
  {"pion_initial_radius"        , 50 , 0      , 1000  , "R_{max} (mm)"            },
  {"pion_final_radius"          , 50 , 0      , 1000  , "R_{min} (mm)"            },
  {"pion_range_z"               , 100, 0      , 3500  , "#Delta z (mm)"           },
  {"pion_vertex_xy"             , 100, 0      , 500   , "v_{xy} (mm)"             },
  {"pion_vertex_z"              , 50 , -1000  , 1000  , "v_{z} (mm)"              },
  {"delta_phi_pion_chargino"    , 50 , 0      , 3.2   , "#Delta #phi (#pi, #chi)" },
  {"delta_theta_pion_chargino"  , 100, -6.3   , 6.3   , "#Delta #theta(#pi, #chi)"},
  {"pion_turning_layer"         , 20 , -1     , 19    , "Layer number"            },
  
  {"next_point_delta_phi"       , 100, -3.2   , 3.2   , "#Delta #phi"             },
  {"next_point_delta_phi_plu"   , 100, -3.2   , 3.2   , "#Delta #phi"             },
  {"next_point_delta_phi_min"   , 100, -3.2   , 3.2   , "#Delta #phi"             },
  {"next_point_delta_z"         , 100, 0      , 500   , "#Delta z"                },
  
  {"middle_seed_hit_delta_phi"      , 100 , -3.2      , 3.2   , "#Delta #phi"      },
  {"middle_seed_hit_delta_phi_plu"  , 100 , -3.2      , 3.2   , "#Delta #phi"      },
  {"middle_seed_hit_delta_phi_min"  , 100 , -3.2      , 3.2   , "#Delta #phi"      },
  
  {"middle_seed_hit_delta_z"        , 50  , 0         , 500   , "#Delta z"         },
  
  {"last_seed_hit_delta_phi"        , 100 , -3.2      , 3.2   , "#Delta #phi"      },
  {"last_seed_hit_delta_phi_plu"    , 100 , -3.2      , 3.2   , "#Delta #phi"      },
  {"last_seed_hit_delta_phi_min"    , 100 , -3.2      , 3.2   , "#Delta #phi"      },
  
  {"last_seed_hit_delta_z"      , 50 , 0      , 500   , "#Delta z"                },
  
  {"pion_cluster_charge_TIB"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_TOB"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_TID"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_TEC"    , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_P1PXB"  , 100, 0      , 600   , ""                        },
  {"pion_cluster_charge_P1PXE"  , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TIB"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TOB"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TID"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TEC"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_P1PXB" , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_P1PXE" , 100, 0      , 600   , ""                        },
  
  {"num_chargino_q_efficiency_vs_gen_pt"  , 100, 0      , 1000  , "p_T (GeV)"           },
  {"den_chargino_q_efficiency_vs_gen_pt"  , 100, 0      , 1000  , "p_T (GeV)"           },
  
  {"num_chargino_q_efficiency_vs_track_pt"  , 100, 0      , 1000  , "p_T (GeV)"           },
  {"den_chargino_q_efficiency_vs_track_pt"  , 100, 0      , 1000  , "p_T (GeV)"           },
  
  // special hists
  {"pion_pz_noTIB"              , 50 , 0      , 2000  , "|p_{z}| (MeV)"           },
  {"pion_pt_noTIB"              , 50 , 0      , 1000  , " p_{t}  (MeV)"           },
  {"pion_vertex_xy_noTIB"       , 100, 0      , 500   , "v_{xy} (mm)"             },
  {"pion_vertex_z_noTIB"        , 50 , -1000  , 1000  , "v_{z} (mm)"              },
  {"noise_n_clusters"           , 100, 0      , 10000 , "# clusters"              },
  {"tracker_clusters_r"         , 1100,0      , 1500  , "R (mm)"                  },
  
  {"last_to_avg_dedx_ratio"     , 500 ,0      , 10    , "dE/dx last / dE/dx avg"  },
};

vector<tuple<string, int, double, double, string, int, double, double, string>> histOptions2D = {
  // title                        nBinsX  minX maxX    titleX nBinsY  minY maxY    titleY
  {"next_point_delta_phi_pion_pt" , 50, 0, 3.14 , "#Delta #phi" , 50, 0, 1000 , "pion p_{t} (MeV)"  },
  {"next_point_delta_z_pion_pz"   , 50, 0, 1000 , "#Delta z"    , 50, 0, 1000 , "pion p_{z} (MeV)"  },
  {"charge_gen_vs_rec"            , 3 ,-1, 2    , "q_{rec}"     , 3 ,-1, 2    , "q_{gen}"           },
  {"layers_gen_vs_rec_0_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_0_5"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_0_7"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_0_8"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_0_9"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_1_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_1_2"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_2_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_3_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_4_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_5_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_6_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"layers_gen_vs_rec_7_0"        ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"pt_gen_vs_rec"                ,100 ,0 , 1000   , "pt_{rec}", 100 , 0, 1000, "pt_{gen}" },
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

shared_ptr<Point> FindClusterForHit(const Point &hit, const vector<shared_ptr<Point>> &clusters)
{
  for(auto &cluster : clusters){
    if(fabs(hit.GetX() - cluster->GetX()) < cluster->GetXerr() &&
       fabs(hit.GetY() - cluster->GetY()) < cluster->GetYerr() &&
       fabs(hit.GetZ() - cluster->GetZ()) < cluster->GetZerr()){
      return cluster;
    }
  }
  return nullptr;
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
    if(event->GetGenCharginoTracks().size() != 1) continue;
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
    auto chargino        = event->GetGenCharginoTracks()[0];
    auto charginoSimHits = event->GetCharginoSimHits();
    
    
    
    vector<double> trackDedx;
    
    double avgTrackDedxWithoutLast=0;
    
    for(int iLayer=1; iLayer<(track->GetLastBarrelLayer()+1); iLayer++){
      avgTrackDedxWithoutLast += track->GetDedxInBarrelLayer(iLayer);
    }
    avgTrackDedxWithoutLast /= (track->GetNtrackerLayers()-1);
    
    double lastToAvgTrackDedxRatio = track->GetDedxInBarrelLayer((track->GetLastBarrelLayer()+1)) / avgTrackDedxWithoutLast;
    
    hists1D["last_to_avg_dedx_ratio"]->Fill(lastToAvgTrackDedxRatio);
    
    auto charginoSimHitsByLayer = pointsProcessor.SortByLayer(charginoSimHits);
    int maxCharginoLayer = -1;
    for(auto &hits : charginoSimHitsByLayer){
      for(auto &hit : hits){
        if(hit->GetLayer() > maxCharginoLayer) maxCharginoLayer = hit->GetLayer();
      }
    }
    int nCharginoLayers = maxCharginoLayer+1;
    
    // Fill basic info about events
    hists1D["lumi"]->Fill(event->GetLumiSection());
    hists1D["noise_n_clusters"]->Fill(trackerClusters.size()-pionClusters.size());
    hists2D["charge_gen_vs_rec"]->Fill(track->GetCharge(), chargino.GetCharge());
    
    hists2D["layers_gen_vs_rec_0_0"]->Fill(track->GetNtrackerLayers(), nCharginoLayers);
    hists2D["layers_gen_vs_rec_0_5"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 0.5 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_0_7"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 0.7 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_0_8"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 0.8 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_0_9"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 0.9 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_1_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 1.0 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_1_2"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 1.2 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_2_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 2.0 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_3_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 2.0 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_4_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 2.0 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_5_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 2.0 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_6_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 2.0 ? 1 : 0), nCharginoLayers);
    hists2D["layers_gen_vs_rec_7_0"]->Fill(track->GetNtrackerLayers()-(lastToAvgTrackDedxRatio < 2.0 ? 1 : 0), nCharginoLayers);

    hists2D["pt_gen_vs_rec"]->Fill(track->GetPt(), chargino.GetPt());
    
    if(track->GetCharge()==chargino.GetCharge()){
      hists1D["num_chargino_q_efficiency_vs_gen_pt"]->Fill(chargino.GetPt());
      hists1D["num_chargino_q_efficiency_vs_track_pt"]->Fill(track->GetPt());
    }
    
    hists1D["den_chargino_q_efficiency_vs_gen_pt"]->Fill(chargino.GetPt());
    hists1D["den_chargino_q_efficiency_vs_track_pt"]->Fill(track->GetPt());
    
    // Fill gen-pion histograms
    double pionPt = sqrt(pow(pionHelix.GetMomentum()->GetX(), 2) + pow(pionHelix.GetMomentum()->GetY(), 2));
    double pionPz = pionHelix.GetMomentum()->GetZ();
    double phiPion = atan2(pionHelix.GetMomentum()->GetY(), pionHelix.GetMomentum()->GetX());
    double pionVertexXY = sqrt(pow(pionHelix.GetVertex()->GetX(), 2) + pow(pionHelix.GetVertex()->GetY(), 2));
    double pionVertexZ = pionHelix.GetVertex()->GetZ();
    
    hists1D["pion_px"]->Fill(fabs(pionHelix.GetMomentum()->GetX()));
    hists1D["pion_py"]->Fill(fabs(pionHelix.GetMomentum()->GetY()));
    hists1D["pion_pz"]->Fill(fabs(pionPz));
    hists1D["pion_pt"]->Fill(pionPt);
    hists1D["pion_vertex_z"]->Fill(pionVertexZ);
    hists1D["pion_vertex_xy"]->Fill(pionVertexXY);
    
    double pionTheta = TMath::Pi()/2-pionHelix.GetMomentum()->GetVectorSlopeC();
    
    hists1D["delta_phi_pion_chargino"]->Fill(fabs(phiPion - track->GetPhi()));
    hists1D["delta_theta_pion_chargino"]->Fill(track->GetTheta() - pionTheta);
    
    // Fill pion sim hit histograms
    pointsProcessor.SetPointsLayers(pionSimHits);
//    pointsProcessor.FilterNearbyPoints(pionSimHits, 100);
    
    sort(pionSimHits.begin(), pionSimHits.end(), PointsProcessor::ComparePointByZ());
    
    int turningLayer = -1;
    
    for(int iHit=0; iHit<pionSimHits.size()-1; iHit++){
//      auto pionCluster = FindClusterForHit(*pionSimHits[iHit],  pionClusters);
    
      auto thisHit = pionSimHits[iHit];
      auto nextHit = pionSimHits[iHit+1];
      
      if(thisHit->GetLayer() == thisHit->GetLayer()){
        turningLayer = pionSimHits[iHit]->GetLayer();
        break;
      }
    }
    hists1D["pion_turning_layer"]->Fill(turningLayer);
    
    Point trackPointMid = pointsProcessor.GetPointOnTrack((layerR[track->GetNtrackerLayers()-1]+layerR[track->GetNtrackerLayers()])/2., *track, *eventVertex);
    
    auto middleSeedCluster = FindClusterForHit(*pionSimHits[0],  pionClusters);// pionSimHits[0];//
    auto lastSeedCluster   = FindClusterForHit(*pionSimHits[1],  pionClusters);// pionSimHits[1];//
    
    if(!middleSeedCluster || !lastSeedCluster) continue;
    
    double middleHitDeltaPhi = pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMid, *middleSeedCluster);
    double middleHitDeltaZ = fabs(middleSeedCluster->GetZ() - trackPointMid.GetZ());
    
    
    double lastHitDeltaPhi = pointsProcessor.GetPointingAngleXY(trackPointMid, *middleSeedCluster, *lastSeedCluster);
    double lastPointDeltaZ = fabs(middleSeedCluster->GetZ() - lastSeedCluster->GetZ());
    
    hists1D["middle_seed_hit_delta_phi"]->Fill(middleHitDeltaPhi);
    hists1D["last_seed_hit_delta_phi"]->Fill(lastHitDeltaPhi);
    
    if(pionHelix.GetCharge() > 0){
      hists1D["middle_seed_hit_delta_phi_plu"]->Fill(middleHitDeltaPhi);
      hists1D["last_seed_hit_delta_phi_plu"]->Fill(lastHitDeltaPhi);
    }
    else{
      hists1D["middle_seed_hit_delta_phi_min"]->Fill(middleHitDeltaPhi);
      hists1D["last_seed_hit_delta_phi_min"]->Fill(lastHitDeltaPhi);
    }
    
    hists1D["middle_seed_hit_delta_z"]->Fill(middleHitDeltaZ);
    hists1D["last_seed_hit_delta_z"]->Fill(lastPointDeltaZ);
    
    for(int iHit=2; iHit<pionSimHits.size()-1; iHit++){
      auto previousCluster = FindClusterForHit(*pionSimHits[iHit-1],  pionClusters);
      auto thisCluster     = FindClusterForHit(*pionSimHits[iHit],    pionClusters);
      auto nextCluster     = FindClusterForHit(*pionSimHits[iHit+1],  pionClusters);
      
      if(!previousCluster || !thisCluster || !nextCluster){
//        cout<<"Could not find cluster for sim hit"<<endl;
        continue;
      }
      
      double deltaPhi = pointsProcessor.GetPointingAngleXY(*previousCluster, *thisCluster, *nextCluster);
      double deltaZ = fabs(thisCluster->GetZ() - nextCluster->GetZ());

      hists1D["next_point_delta_phi"]->Fill(deltaPhi);
      if(pionHelix.GetCharge() > 0) hists1D["next_point_delta_phi_plu"]->Fill(deltaPhi);
      else                          hists1D["next_point_delta_phi_min"]->Fill(deltaPhi);
      
      hists2D["next_point_delta_phi_pion_pt"]->Fill(deltaPhi, pionPt);
      hists1D["next_point_delta_z"]->Fill(deltaZ);
      hists2D["next_point_delta_z_pion_pz"]->Fill(deltaZ, pionPz);
    }
    
//    for(int iCluster=1; iCluster<pionClusters.size()-1; iCluster++){
//      double deltaPhi = pointsProcessor.GetPointingAngleXY(*pionClusters[iCluster-1],
//                                                           *pionClusters[iCluster],
//                                                           *pionClusters[iCluster+1]);
//
//      double deltaZ = fabs(pionClusters[iCluster]->GetZ() - pionClusters[iCluster+1]->GetZ());
//
//      hists1D["next_point_delta_phi"]->Fill(deltaPhi);
//      hists2D["next_point_delta_phi_pion_pt"]->Fill(deltaPhi, pionPt);
//      hists1D["next_point_delta_z"]->Fill(deltaZ);
//      hists2D["next_point_delta_z_pion_pz"]->Fill(deltaZ, pionPz);
//    }
    
    
    shared_ptr<Point> farthestPoint;
    if(fabs(pionSimHits.front()->GetZ()) > fabs(pionSimHits.back()->GetZ()))  farthestPoint = pionSimHits.front();
    else                                                                      farthestPoint = pionSimHits.back();
    
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
    hists1D["pion_n_clusters_P1PXB"]->Fill(nClustersInDet["P1PXB"]);
    hists1D["pion_n_clusters_P1PXE"]->Fill(nClustersInDet["P1PXE"]);
    
    // Fill special histograms
    if(nClustersInDet["TIB"] == 0){
      hists1D["pion_pt_noTIB"]->Fill(pionPt);
      hists1D["pion_pz_noTIB"]->Fill(pionPz);
      hists1D["pion_vertex_xy_noTIB"]->Fill(pionVertexXY);
      hists1D["pion_vertex_z_noTIB"]->Fill(pionVertexZ);
    }
    
    // Fill seeds tracking histograms
  
  
    // Fill noise clusters histograms
    for(auto &cluster : trackerClusters){
      string title = "tracker_cluster_charge_" + cluster->GetSubDetName();
      
      if(hists1D.find(title) != hists1D.end()){
        hists1D[title]->Fill(cluster->GetValue());
      }
      
      double r = sqrt(pow(cluster->GetX(), 2) + pow(cluster->GetY(), 2));
      
      if(cluster->GetSubDetName() == "TIB" ||
         cluster->GetSubDetName() == "TOB" ||
         cluster->GetSubDetName() == "P1PXB"){
        hists1D["tracker_clusters_r"]->Fill(r);
      }
      
    }
    
    nLastHitsInDet[maxDet]++;
  }
  
  TCanvas *genPionCanvas = new TCanvas("gen_pion", "gen_pion", 2880, 1800);
  genPionCanvas->Divide(4,4);
  
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
  genPionCanvas->cd(13);  hists1D["delta_theta_pion_chargino"]->Draw();
  genPionCanvas->cd(14);  hists1D["pion_turning_layer"]->Draw();
  
  TCanvas *trackingCanvas = new TCanvas("tracking", "tracking", 2880, 1800);
  trackingCanvas->Divide(4,4);
  
  trackingCanvas->cd(1);  hists1D["middle_seed_hit_delta_phi"]->Draw();
  trackingCanvas->cd(2);  hists1D["middle_seed_hit_delta_z"]->Draw();
  trackingCanvas->cd(3);  hists1D["last_seed_hit_delta_phi"]->Draw();
  trackingCanvas->cd(4);  hists1D["last_seed_hit_delta_z"]->Draw();
  trackingCanvas->cd(5);  hists1D["next_point_delta_phi"]->Draw();
  trackingCanvas->cd(6);  hists2D["next_point_delta_phi_pion_pt"]->Draw("colz");
  trackingCanvas->cd(7);  hists1D["next_point_delta_z"]->Draw();
  trackingCanvas->cd(8);  hists2D["next_point_delta_z_pion_pz"]->Draw("colz");
  trackingCanvas->cd(9);  hists1D["next_point_delta_phi_plu"]->Draw();
  trackingCanvas->cd(10);  hists1D["next_point_delta_phi_min"]->Draw();
  trackingCanvas->cd(11);  hists1D["middle_seed_hit_delta_phi_plu"]->Draw();
  trackingCanvas->cd(12);  hists1D["middle_seed_hit_delta_phi_min"]->Draw();
  trackingCanvas->cd(13);  hists1D["last_seed_hit_delta_phi_plu"]->Draw();
  trackingCanvas->cd(14);  hists1D["last_seed_hit_delta_phi_min"]->Draw();
  
  TCanvas *charginoCanvas = new TCanvas("chargino", "chargino", 2880, 1800);
  charginoCanvas->Divide(2,3);
  hists2D["charge_gen_vs_rec"]->SetMarkerSize(3.0);
  charginoCanvas->cd(1);  hists2D["charge_gen_vs_rec"]->DrawNormalized("text colz");
  hists2D["layers_gen_vs_rec_0_0"]->SetMarkerSize(3.0);
  charginoCanvas->cd(2);  hists2D["layers_gen_vs_rec_0_0"]->DrawNormalized("text colz");
  
  hists1D["num_chargino_q_efficiency_vs_gen_pt"]->Divide(hists1D["den_chargino_q_efficiency_vs_gen_pt"]);
  charginoCanvas->cd(3);  hists1D["num_chargino_q_efficiency_vs_gen_pt"]->Draw();
  
  hists1D["num_chargino_q_efficiency_vs_track_pt"]->Divide(hists1D["den_chargino_q_efficiency_vs_track_pt"]);
  charginoCanvas->cd(4);  hists1D["num_chargino_q_efficiency_vs_track_pt"]->Draw();
  
  charginoCanvas->cd(5); hists2D["pt_gen_vs_rec"]->DrawNormalized("colz");
  
  TCanvas *chargeCanvas   = new TCanvas("charge", "charge", 2880, 1800);
  TCanvas *clustersCanvas = new TCanvas("clusters", "clusters", 2880, 1800);
  
  chargeCanvas->Divide(2,3);
  clustersCanvas->Divide(2,3);
  
  vector<string> detNames = {"TIB", "TOB", "TID", "TEC", "P1PXB", "P1PXE" };
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
  
  TCanvas *specialCanvas = new TCanvas("specialCanvas", "specialCanvas", 2880, 1800);
  specialCanvas->Divide(3,3);

  specialCanvas->cd(1);  hists1D["pion_pt_noTIB"]->Draw();
  specialCanvas->cd(2);  hists1D["pion_pz_noTIB"]->Draw();
  specialCanvas->cd(3);  hists1D["pion_vertex_xy_noTIB"]->Draw();
  specialCanvas->cd(4);  hists1D["pion_vertex_z_noTIB"]->Draw();
  specialCanvas->cd(5);  hists1D["noise_n_clusters"]->Draw();
  specialCanvas->cd(6);  hists1D["tracker_clusters_r"]->Draw();
  specialCanvas->cd(7);  hists1D["last_to_avg_dedx_ratio"]->Draw();
  
  genPionCanvas->Update();
  trackingCanvas->Update();
  chargeCanvas->Update();
  clustersCanvas->Update();
  specialCanvas->Update();
  charginoCanvas->Update();
  
  TFile *outFile = new TFile(outfileName.c_str(), "recreate");
  outFile->cd();
  
  for(auto &[name, hist] : hists1D) hist->Write();
  for(auto &[name, hist] : hists2D) hist->Write();
  
  outFile->Close();
  
  theApp.Run();
}
