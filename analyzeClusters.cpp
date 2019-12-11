//  analyzeClusters.cpp
//
//  Created by Jeremi Niedziela on 10/05/2019.

#include "Helpers.hpp"
#include "EventSet.hpp"

string configPath  = "configs/clusters.md";
string sampleTag = "500_10_test";

vector<tuple<string, int, double, double, string>> histOptions1D = {
// title                        nBins  min     max    titleX
  {"lumi"                       , 300, 0      , 300   , "lumi block"              },
  {"pion_px"                    , 50 , 0      , 1000  , "|p_{x}| (MeV)"           },
  {"pion_py"                    , 50 , 0      , 1000  , "|p_{y}| (MeV)"           },
  {"pion_pz"                    , 50 , 0      , 1500  , "|p_{z}^{gen}| (MeV)"     },
  {"pion_pt"                    , 50 , 0      , 1000  , "p_{t}^{gen}  (MeV)"      },
  {"pion_n_clusters"            , 50 , 0      , 100   , "# clusters"              },
  {"pion_n_clusters_TIB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TOB"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TID"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_TEC"        , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_P1PXB"      , 50 , 0      , 50    , "# clusters"              },
  {"pion_n_clusters_P1PXEC"     , 50 , 0      , 50    , "# clusters"              },
  {"pion_initial_radius"        , 50 , 0      , 1000  , "R_{max} (mm)"            },
  {"pion_final_radius"          , 50 , 0      , 1000  , "R_{min} (mm)"            },
  {"pion_simhits_z"             , 100, 0      , 3500  , "Pion simHit z (mm)"      },
  {"pion_simhits_range_z"       , 100, 0      , 3500  , "Pion simHit range z (mm)"},
  {"pion_vertex_xy"             , 100, 0      , 500   , "v_{xy} (mm)"             },
  {"pion_vertex_z"              , 50 , -1000  , 1000  , "v_{z} (mm)"              },
  {"delta_phi_pion_chargino"    , 50 , 0      , 3.2   , "#Delta #phi^{gen} (#pi, #chi)"   },
  {"delta_theta_pion_chargino"  , 100, -6.3   , 6.3   , "#Delta #theta^{gen} (#pi, #chi)" },
  {"pion_turning_layer"         , 20 , -1     , 19    , "Layer number"            },
  {"chargino_abs_eta"           , 20 ,  0     , 4.0   , "#chi^{#pm} |#eta^{gen}|"     },
  {"chargino_n_layers_rec"      , 20 ,  0     , 20    , "#chi^{#pm} N^{rec}_{layers}" },
  {"chargino_n_layers_gen"      , 20 ,  0     , 20    , "#chi^{#pm} N^{gen}_{layers}" },
  
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
  {"pion_cluster_charge_P1PXEC"  , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TIB"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TOB"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TID"     , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_TEC"     , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_P1PXB"   , 100, 0      , 600   , ""                        },
  {"tracker_cluster_charge_P1PXEC"  , 100, 0      , 600   , ""                        },
  {"pion_cluster_z"                 , 100, 0      , 3500   , "Pion cluster z"       },
  {"pion_cluster_range_z"           , 100, 0      , 3500   , "Pion cluster range z" },
  
  {"num_chargino_q_efficiency_vs_gen_pt"    , 100, 0  , 1000  , "p_T (GeV)" },
  {"den_chargino_q_efficiency_vs_gen_pt"    , 100, 0  , 1000  , "p_T (GeV)" },
  
  {"num_chargino_q_efficiency_vs_track_pt"  , 100, 0  , 1000  , "p_T (GeV)" },
  {"den_chargino_q_efficiency_vs_track_pt"  , 100, 0  , 1000  , "p_T (GeV)" },
  
  // special hists
  {"pion_pz_noTIB"              , 50 , 0      , 2000  , "|p_{z}| (MeV)"           },
  {"pion_pt_noTIB"              , 50 , 0      , 1000  , " p_{t}  (MeV)"           },
  {"pion_vertex_xy_noTIB"       , 100, 0      , 500   , "v_{xy} (mm)"             },
  {"pion_vertex_z_noTIB"        , 50 , -1000  , 1000  , "v_{z} (mm)"              },
  {"noise_n_clusters"           , 50 , 0      , 50000 , "# clusters"              },
  {"tracker_clusters_r"         , 1100,0      , 1500  , "R (mm)"                  },
  {"tracker_clusters_z_PXE"     , 100 ,-600   , 600   , "z (mm)"                  },
  {"tracker_clusters_z_TID"     , 400 ,-1200  , 1200  , "z (mm)"                  },
  {"tracker_clusters_z_TEC"     , 300 ,-3000  , 3000  , "z (mm)"                  },
  
  {"last_to_avg_dedx_ratio"     , 500 ,0      , 10    , "dE/dx last / dE/dx avg"  },
};

vector<tuple<string, int, double, double, string, int, double, double, string>> histOptions2D = {
  // title                        nBinsX  minX maxX    titleX nBinsY  minY maxY    titleY
  {"next_point_delta_phi_pion_pt" , 50, 0, 3.14 , "#Delta #phi" , 50, 0, 1000 , "pion p_{t} (MeV)"  },
  {"next_point_delta_z_pion_pz"   , 50, 0, 1000 , "#Delta z"    , 50, 0, 1000 , "pion p_{z} (MeV)"  },
  {"charge_gen_vs_rec"            , 3 ,-1, 2    , "q_{rec}"     , 3 ,-1, 2    , "q_{gen}"           },
  {"layers_gen_vs_rec"            ,10 ,0 , 10   , "N^{layers}_{rec}", 10 , 0, 10, "q^{layer}_{gen}" },
  {"pt_gen_vs_rec"                ,100 ,0 , 1000   , "pt_{rec}", 100 , 0, 1000, "pt_{gen}" },
  {"pt_gen_vs_rec_good_charge"    ,100 ,0 , 1000   , "pt_{rec}", 100 , 0, 1000, "pt_{gen}" },
  {"pion_pt_vs_n_layers"          ,100 ,0 , 1000   , "pion pt (MeV)", 20 , 0, 20, "n chargino layers" },
  {"eta_vs_n_layers"              ,50  ,0 , 2.5    , "track |#eta|" , 20 , 0, 20, "n chargino layers" },
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

range<double> GetPointsRingSize(Points points)
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

shared_ptr<Point> FindClusterForHit(const Point &hit, const Points &clusters)
{
  for(auto &cluster : clusters){
    if(fabs(hit.GetX() - cluster->GetX()) < 1 &&
       fabs(hit.GetY() - cluster->GetY()) < 1 &&
       fabs(hit.GetZ() - cluster->GetZ()) < cluster->GetZerr()){
      return cluster;
    }
  }
  return nullptr;
}

void fillGenPionHists(const shared_ptr<Event> &event)
{
  auto pionHelix = event->GetGenPionHelices()[0];
  
  // Fill gen-pion histograms
  double pionPt = sqrt(pow(pionHelix.GetMomentum().GetX(), 2) + pow(pionHelix.GetMomentum().GetY(), 2));
  double pionPz = pionHelix.GetMomentum().GetZ();
  double phiPion = atan2(pionHelix.GetMomentum().GetY(), pionHelix.GetMomentum().GetX());
  double pionVertexXY = sqrt(pow(pionHelix.GetVertex()->GetX(), 2) + pow(pionHelix.GetVertex()->GetY(), 2));
  double pionVertexZ = pionHelix.GetVertex()->GetZ();
  double pionTheta = TMath::Pi()/2-pionHelix.GetMomentum().GetVectorSlopeC();
  
  hists1D["pion_px"]->Fill(fabs(pionHelix.GetMomentum().GetX()));
  hists1D["pion_py"]->Fill(fabs(pionHelix.GetMomentum().GetY()));
  hists1D["pion_pz"]->Fill(fabs(pionPz));
  hists1D["pion_pt"]->Fill(pionPt);
  hists1D["pion_vertex_z"]->Fill(pionVertexZ);
  hists1D["pion_vertex_xy"]->Fill(pionVertexXY);
  
  if(event->GetNtracks() == 1){
    auto track = event->GetTrack(0);
    hists2D["pion_pt_vs_n_layers"]->Fill(pionHelix.GetMomentum().GetTransverse(), track->GetNtrackerLayers());
  }
  if(event->GetGenCharginoTracks().size() == 1){
    auto chargino = event->GetGenCharginoTracks()[0];
    hists1D["delta_phi_pion_chargino"]->Fill(fabs(chargino.GetPhi() - phiPion));
    hists1D["delta_theta_pion_chargino"]->Fill(chargino.GetTheta() - pionTheta);
  }
  
}

void fillGenCharginoHists(const shared_ptr<Event> &event)
{
  auto chargino = event->GetGenCharginoTracks()[0];
  auto charginoSimHits = event->GetCharginoSimHits();
  auto charginoSimHitsByLayer = pointsProcessor.SortByLayer(charginoSimHits);
  
  hists1D["den_chargino_q_efficiency_vs_gen_pt"]->Fill(chargino.GetPt());
  hists1D["chargino_abs_eta"]->Fill(fabs(chargino.GetEta()));
  hists1D["chargino_n_layers_gen"]->Fill(chargino.GetNtrackerLayers());
  
  int maxCharginoLayer = -1;
  for(auto &hits : charginoSimHitsByLayer){
    for(auto &hit : hits){
      if(hit->GetLayer() > maxCharginoLayer) maxCharginoLayer = hit->GetLayer();
    }
  }
  int nCharginoLayers = maxCharginoLayer+1;
  
  if(event->GetNtracks() == 1){
    auto track = event->GetTrack(0);
    
    hists2D["charge_gen_vs_rec"]->Fill(track->GetCharge(), chargino.GetCharge());
    hists2D["layers_gen_vs_rec"]->Fill(track->GetNtrackerLayers(), nCharginoLayers);
    hists2D["pt_gen_vs_rec"]->Fill(track->GetPt(), chargino.GetPt());
    
    if(track->GetCharge()==chargino.GetCharge()){
      hists1D["num_chargino_q_efficiency_vs_gen_pt"]->Fill(chargino.GetPt());
      hists1D["num_chargino_q_efficiency_vs_track_pt"]->Fill(track->GetPt());
      hists2D["pt_gen_vs_rec_good_charge"]->Fill(track->GetPt(), chargino.GetPt());
    }
  }
}

void fillTrackHists(const shared_ptr<Event> &event)
{
  auto track = event->GetTrack(0);
  
  hists2D["eta_vs_n_layers"]->Fill(fabs(track->GetEta()), track->GetNtrackerLayers());
  hists1D["den_chargino_q_efficiency_vs_track_pt"]->Fill(track->GetPt());
  hists1D["chargino_n_layers_rec"]->Fill(track->GetNtrackerLayers());
  
  double avgTrackDedxWithoutLast=0;
  
  for(int iLayer=1; iLayer<(track->GetLastBarrelLayer()+1); iLayer++){
    avgTrackDedxWithoutLast += track->GetDedxInBarrelLayer(iLayer);
  }
  avgTrackDedxWithoutLast /= (track->GetNtrackerLayers()-1);
  
  
  double lastToAvgTrackDedxRatio = track->GetDedxInBarrelLayer((track->GetLastBarrelLayer()+1)) / avgTrackDedxWithoutLast;
  hists1D["last_to_avg_dedx_ratio"]->Fill(lastToAvgTrackDedxRatio);
  
}

void fillClusterHists(const shared_ptr<Event> &event)
{
  if(!event->HasFriendData()) return;
  
  auto trackerClusters = event->GetTrackerClusters();
  auto trackerClustersPerLayer = pointsProcessor.SortByLayer(trackerClusters);
  
  auto pionClusters = event->GetPionClusters();
  auto pionSimHits  = event->GetPionSimHits();
  pointsProcessor.SetPointsLayers(pionSimHits);
  sort(pionSimHits.begin(), pionSimHits.end(), PointsProcessor::ComparePointByTime());
  
  unsigned long nNoise = trackerClusters.size()-pionClusters.size();
  hists1D["noise_n_clusters"]->Fill(nNoise);
  
  if(pionClusters.size() >= 3 && pionSimHits.size() >= 3){
    if(event->GetNtracks() == 1){
      auto track        = event->GetTrack(0);
      auto eventVertex  = event->GetVertex();
      
      int turningLayer = -1;
      
      for(int iHit=0; iHit<pionSimHits.size()-1; iHit++){
        auto thisHit = pionSimHits[iHit];
        auto nextHit = pionSimHits[iHit+1];
        
        if(thisHit->GetLayer() == thisHit->GetLayer()){
          turningLayer = pionSimHits[iHit]->GetLayer();
          break;
        }
      }
      hists1D["pion_turning_layer"]->Fill(turningLayer);
      
      Point trackPointMid = pointsProcessor.GetPointOnTrack((layerR[track->GetNtrackerLayers()-1]+
                                                             layerR[track->GetNtrackerLayers()])/2.,
                                                            *track, *eventVertex);
      
      auto middleSeedCluster = FindClusterForHit(*pionSimHits[0],  pionClusters);// pionSimHits[0];//
      auto lastSeedCluster   = FindClusterForHit(*pionSimHits[1],  pionClusters);// pionSimHits[1];//
      
      if(middleSeedCluster && lastSeedCluster){
        
        double middleHitDeltaPhi = pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMid, *middleSeedCluster);
        double middleHitDeltaZ = fabs(middleSeedCluster->GetZ() - trackPointMid.GetZ());
        
        
        double lastHitDeltaPhi = pointsProcessor.GetPointingAngleXY(trackPointMid, *middleSeedCluster, *lastSeedCluster);
        double lastPointDeltaZ = fabs(middleSeedCluster->GetZ() - lastSeedCluster->GetZ());
        
        hists1D["middle_seed_hit_delta_phi"]->Fill(middleHitDeltaPhi);
        hists1D["last_seed_hit_delta_phi"]->Fill(lastHitDeltaPhi);
        
        if(event->GetGenPionHelices().size()==1){
          auto pionHelix = event->GetGenPionHelices()[0];
          
          if(pionHelix.GetCharge() > 0){
            hists1D["middle_seed_hit_delta_phi_plu"]->Fill(middleHitDeltaPhi);
            hists1D["last_seed_hit_delta_phi_plu"]->Fill(lastHitDeltaPhi);
          }
          else{
            hists1D["middle_seed_hit_delta_phi_min"]->Fill(middleHitDeltaPhi);
            hists1D["last_seed_hit_delta_phi_min"]->Fill(lastHitDeltaPhi);
          }
        }
        
        hists1D["middle_seed_hit_delta_z"]->Fill(middleHitDeltaZ);
        hists1D["last_seed_hit_delta_z"]->Fill(lastPointDeltaZ);
      }
    }
  }
  
  
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
    else if(cluster->GetSubDetName() == "TID")    hists1D["tracker_clusters_z_TID"]->Fill(cluster->GetZ());
    else if(cluster->GetSubDetName() == "TEC")    hists1D["tracker_clusters_z_TEC"]->Fill(cluster->GetZ());
    else if(cluster->GetSubDetName() == "P1PXEC") hists1D["tracker_clusters_z_PXE"]->Fill(cluster->GetZ());
  }
}

void fillPionPointHists(const shared_ptr<Event> &event)
{
  auto pionClusters = event->GetPionClusters();
  auto pionSimHits  = event->GetPionSimHits();
  
  pointsProcessor.SetPointsLayers(pionSimHits);
  pointsProcessor.SetPointsLayers(pionClusters);
  
  vector<Points> pionSimHitsRegrouped = pointsProcessor.RegroupNerbyPoints(pionSimHits);
  
  Points pionSimHitsNoDouble;
  
  for(Points hits : pionSimHitsRegrouped) pionSimHitsNoDouble.push_back(hits.front());
  
//  sort(pionSimHits.begin(), pionSimHits.end(), PointsProcessor::ComparePointByZ());
//  sort(pionClusters.begin(), pionClusters.end(), PointsProcessor::ComparePointByZ());
  
  sort(pionSimHitsNoDouble.begin(), pionSimHitsNoDouble.end(), PointsProcessor::ComparePointByTime());
  sort(pionClusters.begin(), pionClusters.end(), PointsProcessor::ComparePointByTime());
  
  for(int iHit=2; iHit<pionSimHitsNoDouble.size()-1; iHit++){
    auto previousCluster = FindClusterForHit(*pionSimHitsNoDouble[iHit-1],  pionClusters);
    auto thisCluster     = FindClusterForHit(*pionSimHitsNoDouble[iHit],    pionClusters);
    auto nextCluster     = FindClusterForHit(*pionSimHitsNoDouble[iHit+1],  pionClusters);
    
    if(!previousCluster || !thisCluster || !nextCluster) continue;
    
    double deltaPhi = pointsProcessor.GetPointingAngleXY(*previousCluster, *thisCluster, *nextCluster);
    double deltaZ = fabs(thisCluster->GetZ() - nextCluster->GetZ());
    
    if(event->GetGenPionHelices().size() == 1){
      auto pionHelix = event->GetGenPionHelices()[0];
      
      // Fill gen-pion histograms
      double pionPt = sqrt(pow(pionHelix.GetMomentum().GetX(), 2) + pow(pionHelix.GetMomentum().GetY(), 2));
      double pionPz = pionHelix.GetMomentum().GetZ();

      hists1D["next_point_delta_phi"]->Fill(deltaPhi);
      if(pionHelix.GetCharge() > 0) hists1D["next_point_delta_phi_plu"]->Fill(deltaPhi);
      else                          hists1D["next_point_delta_phi_min"]->Fill(deltaPhi);
      
      hists2D["next_point_delta_phi_pion_pt"]->Fill(deltaPhi, pionPt);
      hists1D["next_point_delta_z"]->Fill(deltaZ);
      hists2D["next_point_delta_z_pion_pz"]->Fill(deltaZ, pionPz);
    }
  }
  
  for(auto &hit : pionSimHitsNoDouble){
    hists1D["pion_simhits_z"]->Fill(fabs(hit->GetZ()));
  }
  
  shared_ptr<Point> farthestPoint;
  if(fabs(pionSimHitsNoDouble.front()->GetZ()) > fabs(pionSimHitsNoDouble.back()->GetZ()))  farthestPoint = pionSimHitsNoDouble.front();
  else                                                                      farthestPoint = pionSimHitsNoDouble.back();
  
  hists1D["pion_simhits_range_z"]->Fill(fabs(farthestPoint->GetZ()));
  
  range<double> pionRingSize = GetPointsRingSize(pionSimHitsNoDouble);
  
  hists1D["pion_initial_radius"]->Fill(pionRingSize.GetMax());
  hists1D["pion_final_radius"]->Fill(pionRingSize.GetMin());
  
  double maxZ = -inf;
  string maxDet;
  map<string, int> nClustersInDet;
  
  // Fill pion rec clusters histograms
  for(auto &cluster : pionClusters){
    string detName = cluster->GetSubDetName();
    string title = "pion_cluster_charge_" + detName;
    if(hists1D.find(title) != hists1D.end()) hists1D[title]->Fill(cluster->GetValue());

    nClustersInDet[detName]++;
    hists1D["pion_cluster_z"]->Fill(fabs(cluster->GetZ()));
  }
  
  
  hists1D["pion_n_clusters"]->Fill(pionClusters.size());
  hists1D["pion_n_clusters_TIB"]->Fill(nClustersInDet["TIB"]);
  hists1D["pion_n_clusters_TOB"]->Fill(nClustersInDet["TOB"]);
  hists1D["pion_n_clusters_TID"]->Fill(nClustersInDet["TID"]);
  hists1D["pion_n_clusters_TEC"]->Fill(nClustersInDet["TEC"]);
  hists1D["pion_n_clusters_P1PXB"]->Fill(nClustersInDet["P1PXB"]);
  hists1D["pion_n_clusters_P1PXEC"]->Fill(nClustersInDet["P1PXEC"]);
  
  
  
  shared_ptr<Point> farthestCluster;
  if(fabs(pionClusters.front()->GetZ()) > fabs(pionClusters.back()->GetZ())) farthestCluster = pionClusters.front();
  else                                                                       farthestCluster = pionClusters.back();
  
  hists1D["pion_cluster_range_z"]->Fill(fabs(farthestCluster->GetZ()));
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  
  int part = 5;
  int nEvents = 500;
  
  if(argc==3){
    nEvents = atoi(argv[1]);
    part = atoi(argv[2]);
  }
  cout<<"Running part "<<part<<" with "<<nEvents<<" events"<<endl;
  
  int eventOffset = part*nEvents;
  string outfileName = "results/clusters_"+sampleTag+"_p"+to_string(part)+".root";
  
  string cutLevel = "";
  if(config.params["cuts_level"]==0) cutLevel += "after_L0/";
  if(config.params["cuts_level"]==1) cutLevel += "after_L1/"+config.category+"/";
  
  EventSet events; //events.LoadEventsFromFiles(cutLevel);
  SetupHists();
  
  cout<<"Filling histograms"<<endl;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
  
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      
      for(auto iEvent=eventOffset; iEvent<nEvents+eventOffset; iEvent++){
        if(iEvent%10==0) cout<<"\niEvent: "<<iEvent<<"\n"<<endl;
        events.LoadEventsFromFiles(kSignal, iSig, cutLevel, iEvent);
        auto event = events.At(kSignal, iSig, year, iEvent-eventOffset);
        
//      for(int iEvent=0; iEvent<events.size(kSignal, iSig, year); iEvent++){
//        auto event = events.At(kSignal, iSig, year, iEvent);
        
        if(!event){
          cout<<"Event not found"<<endl;
          exit(0);
        }
        
        // Pick only events with one track and one gen pion and resonable number of pion clusters
        if(event->GetGenPionHelices().size() == 1)    fillGenPionHists(event);
        if(event->GetNtracks() == 1)                  fillTrackHists(event);
        if(event->GetGenCharginoTracks().size() == 1) fillGenCharginoHists(event);
        if(event->GetPionClusters().size() >= 3 &&
           event->GetPionSimHits().size()  >= 3)      fillPionPointHists(event);
        fillClusterHists(event);
        
        hists1D["lumi"]->Fill(event->GetLumiSection());
      }
    }
  }
  
//  TCanvas *genPionCanvas = new TCanvas("gen_pion", "gen_pion", 2880, 1800);
//  genPionCanvas->Divide(4,4);
//
//  genPionCanvas->cd(1);   hists1D["lumi"]->Draw();
//  genPionCanvas->cd(2);   hists1D["pion_pz"]->Draw();
//  genPionCanvas->cd(3);   hists1D["pion_px"]->Draw();
//  genPionCanvas->cd(4);   hists1D["pion_py"]->Draw();
//  genPionCanvas->cd(5);   hists1D["pion_pt"]->Draw();
//  genPionCanvas->cd(6);   hists1D["pion_n_clusters"]->Draw();
//  genPionCanvas->cd(7);   hists1D["pion_initial_radius"]->Draw();
//  genPionCanvas->cd(8);   hists1D["pion_final_radius"]->Draw();
//  genPionCanvas->cd(9);   hists1D["pion_simhits_range_z"]->Draw();
//  genPionCanvas->cd(10);  hists1D["pion_simhits_z"]->Draw();
//  genPionCanvas->cd(11);  hists1D["pion_vertex_z"]->Draw();
//  genPionCanvas->cd(12);  hists1D["pion_vertex_xy"]->Draw();
//  genPionCanvas->cd(13);  hists1D["delta_phi_pion_chargino"]->Draw();
//  genPionCanvas->cd(14);  hists1D["delta_theta_pion_chargino"]->Draw();
//  genPionCanvas->cd(15);  hists1D["pion_turning_layer"]->Draw();
//
//  TCanvas *trackingCanvas = new TCanvas("tracking", "tracking", 2880, 1800);
//  trackingCanvas->Divide(4,4);
//
//  trackingCanvas->cd(1);  hists1D["middle_seed_hit_delta_phi"]->Draw();
//  trackingCanvas->cd(2);  hists1D["middle_seed_hit_delta_z"]->Draw();
//  trackingCanvas->cd(3);  hists1D["last_seed_hit_delta_phi"]->Draw();
//  trackingCanvas->cd(4);  hists1D["last_seed_hit_delta_z"]->Draw();
//  trackingCanvas->cd(5);  hists1D["next_point_delta_phi"]->Draw();
//  trackingCanvas->cd(6);  hists2D["next_point_delta_phi_pion_pt"]->Draw("colz");
//  trackingCanvas->cd(7);  hists1D["next_point_delta_z"]->Draw();
//  trackingCanvas->cd(8);  hists2D["next_point_delta_z_pion_pz"]->Draw("colz");
//  trackingCanvas->cd(9);  hists1D["next_point_delta_phi_plu"]->Draw();
//  trackingCanvas->cd(10);  hists1D["next_point_delta_phi_min"]->Draw();
//  trackingCanvas->cd(11);  hists1D["middle_seed_hit_delta_phi_plu"]->Draw();
//  trackingCanvas->cd(12);  hists1D["middle_seed_hit_delta_phi_min"]->Draw();
//  trackingCanvas->cd(13);  hists1D["last_seed_hit_delta_phi_plu"]->Draw();
//  trackingCanvas->cd(14);  hists1D["last_seed_hit_delta_phi_min"]->Draw();
  
//  TCanvas *charginoCanvas = new TCanvas("chargino", "chargino", 2880, 1800);
//  charginoCanvas->Divide(3,3);
  hists2D["charge_gen_vs_rec"]->SetMarkerSize(3.0);
//  charginoCanvas->cd(1);  hists2D["charge_gen_vs_rec"]->DrawNormalized("text colz");
  hists2D["layers_gen_vs_rec"]->SetMarkerSize(3.0);
//  charginoCanvas->cd(2);  hists2D["layers_gen_vs_rec"]->DrawNormalized("text colz");
  
  hists1D["num_chargino_q_efficiency_vs_gen_pt"]->Divide(hists1D["den_chargino_q_efficiency_vs_gen_pt"]);
//  charginoCanvas->cd(3);  hists1D["num_chargino_q_efficiency_vs_gen_pt"]->Draw();
  
  hists1D["num_chargino_q_efficiency_vs_track_pt"]->Divide(hists1D["den_chargino_q_efficiency_vs_track_pt"]);
//  charginoCanvas->cd(4);  hists1D["num_chargino_q_efficiency_vs_track_pt"]->Draw();
//
//  charginoCanvas->cd(5); hists2D["pt_gen_vs_rec"]->DrawNormalized("colz");
//  charginoCanvas->cd(6); hists2D["pt_gen_vs_rec_good_charge"]->DrawNormalized("colz");
//
//  charginoCanvas->cd(7); hists2D["pion_pt_vs_n_layers"]->DrawNormalized("colz");
//  charginoCanvas->cd(8); hists2D["eta_vs_n_layers"]->DrawNormalized("colz");
//
//  TCanvas *chargeCanvas   = new TCanvas("charge", "charge", 2880, 1800);
//  TCanvas *clustersCanvas = new TCanvas("clusters", "clusters", 2880, 1800);
  
//  chargeCanvas->Divide(2,3);
//  clustersCanvas->Divide(2,3);
//
//  vector<string> detNames = {"TIB", "TOB", "TID", "TEC", "P1PXB", "P1PXEC" };
//
//  for(int i=0; i<detNames.size(); i++){
//    chargeCanvas->cd(i+1);
//    hists1D["pion_cluster_charge_"+detNames[i]]->SetLineColor(kRed);
//    hists1D["tracker_cluster_charge_"+detNames[i]]->SetLineColor(kBlue);
//    hists1D["pion_cluster_charge_"+detNames[i]]->DrawNormalized();
//    hists1D["tracker_cluster_charge_"+detNames[i]]->DrawNormalized("same");
//
//    clustersCanvas->cd(i+1);
//    hists1D["pion_n_clusters_"+detNames[i]]->Draw();
//  }
//
//  TCanvas *specialCanvas = new TCanvas("specialCanvas", "specialCanvas", 2880, 1800);
//  specialCanvas->Divide(3,4);

//  specialCanvas->cd(1);  hists1D["pion_pt_noTIB"]->Draw();
//  specialCanvas->cd(2);  hists1D["pion_pz_noTIB"]->Draw();
//  specialCanvas->cd(3);  hists1D["pion_vertex_xy_noTIB"]->Draw();
//  specialCanvas->cd(4);  hists1D["pion_vertex_z_noTIB"]->Draw();
//  specialCanvas->cd(5);  hists1D["noise_n_clusters"]->Draw();
//  specialCanvas->cd(6);  hists1D["tracker_clusters_r"]->Draw();
//  specialCanvas->cd(7);  hists1D["tracker_clusters_z_PXE"]->Draw();
//  specialCanvas->cd(8);  hists1D["tracker_clusters_z_TID"]->Draw();
//  specialCanvas->cd(9);  hists1D["tracker_clusters_z_TEC"]->Draw();
//  specialCanvas->cd(10);  hists1D["pion_cluster_z"]->Draw();
//  specialCanvas->cd(11);  hists1D["pion_cluster_range_z"]->Draw();
////  specialCanvas->cd(8);  hists1D["last_to_avg_dedx_ratio"]->Draw();
//
//  genPionCanvas->Update();
//  trackingCanvas->Update();
//  chargeCanvas->Update();
//  clustersCanvas->Update();
//  specialCanvas->Update();
//  charginoCanvas->Update();
  
  TFile *outFile = new TFile(outfileName.c_str(), "recreate");
  outFile->cd();
  
  for(auto &[name, hist] : hists1D) hist->Write();
  for(auto &[name, hist] : hists2D) hist->Write();
  
  outFile->Close();
  
//  theApp.Run();
}
