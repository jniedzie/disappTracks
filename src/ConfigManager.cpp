//
//  FitterConfig.cpp
//
//  Created by Jeremi Niedziela on 09/01/2019.
//

#include "ConfigManager.hpp"

ConfigManager config("init");

ConfigManager::ConfigManager(string _path)
{
  if(_path=="init") return;
  
  configFile = make_unique<TEnv>();
  
  if(configFile->ReadFile(_path.c_str(), kEnvUser) < 0){
    cout<<"ERROR - could not load config file:"<<_path<<". Will use default values"<<endl;
    return;
  }
  else{
    cout<<"INFO -- Successfully read config file: "<<_path<<endl;
  }
  
  verbosity   = configFile->GetValue("verbosity_level", 0);
  
  doubleHitsMaxDistance     = configFile->GetValue("double_hit_max_distance", 0.0);
  doAsymmetricConstraints   = configFile->GetValue("do_asymmetric_constraints", 0);
  
  seedMaxChi2               = configFile->GetValue("seed_max_chi2",100.0);
  seedMiddleHitMaxDeltaZ    = configFile->GetValue("seed_middle_hit_max_delta_z",1000.0);
  seedLastHitMaxDeltaZ      = configFile->GetValue("seed_last_hit_max_delta_z",1000.0);
  trackMaxChi2              = configFile->GetValue("track_max_chi2",100.0);
  
  if(doAsymmetricConstraints){
    seedMiddleHitDeltaPhi     = range<double>(configFile->GetValue("seed_middle_hit_min_delta_phi",-3.14),
                                              configFile->GetValue("seed_middle_hit_max_delta_phi", 3.14));
    
    seedLastHitDeltaPhi       = range<double>(configFile->GetValue("seed_last_hit_min_delta_phi",-3.14),
                                              configFile->GetValue("seed_last_hit_max_delta_phi", 3.14));
    
    nextPointDeltaPhi         = range<double>(configFile->GetValue("next_point_min_delta_phi",-3.14),
                                              configFile->GetValue("next_point_max_delta_phi", 3.14));
  }
  else{
    seedMiddleHitDeltaPhi     = range<double>(0,
                                              fabs(configFile->GetValue("seed_middle_hit_max_delta_phi", 3.14)));
    
    seedLastHitDeltaPhi       = range<double>(0,
                                              fabs(configFile->GetValue("seed_last_hit_max_delta_phi", 3.14)));
    
    nextPointDeltaPhi         = range<double>(0,
                                              fabs(configFile->GetValue("next_point_max_delta_phi",3.14)));
  }
  
  nextPointMaxDeltaZ        = configFile->GetValue("next_point_max_delta_z",1000.0);
  nextPointMaxDeltaXY       = configFile->GetValue("next_point_max_delta_xy",1000.0);
  nextPointMaxDeltaT        = configFile->GetValue("next_point_max_delta_t",1000.0);
  
  mergingMaxDifferentPoints = configFile->GetValue("merging_max_different_point", 0);
  candidateMinNpoints       = configFile->GetValue("candidate_min_n_points", 0);
  trackMinNpoints           = configFile->GetValue("track_min_n_points", 0);
  trackMinNlayers           = configFile->GetValue("track_min_n_layers", 0);
  
  maxNmissingHits           = configFile->GetValue("max_n_missing_hits", 0);
  maxNmissingHitsInRow      = configFile->GetValue("max_n_missing_hits_in_raw", 0);
  mergeAtTurnBack           = configFile->GetValue("merge_at_turn_back", 0);
  mergeFinalHelices         = configFile->GetValue("merge_final_helices", 0);
  allowTurningBack          = configFile->GetValue("allow_turning_back", 0);
  requireGoodStartingValues = configFile->GetValue("require_good_starting_values", 0);
  expRadiusFunction         = configFile->GetValue("exp_radius_function", 0);
  expSlopeFunction          = configFile->GetValue("exp_slope_function", 0);
  
  allowOneLessLayer          = configFile->GetValue("allow_one_less_layer", 0);
  allowOneMoreLayer          = configFile->GetValue("allow_one_more_layer", 0);
  checkOppositeChargeBelowNlayers = configFile->GetValue("check_opposite_charge_below_Nlayers", 0);
  minLayersForDeltaXY         = configFile->GetValue("min_layers_for_delta_xy", inf);
  
  startR0 = configFile->GetValue("start_R0", 320.0);
  minR0 = configFile->GetValue("min_R0", 50.0);
  maxR0 = configFile->GetValue("max_R0", 1000.0);
  minRslope = configFile->GetValue("min_Rslope", 0.0);
  maxRslope = configFile->GetValue("max_Rslope", 10000.0);
  minS0 = configFile->GetValue("min_S0", -10000.0);
  maxS0 = configFile->GetValue("max_S0", 10000.0);
  minSslope = configFile->GetValue("min_Sslope", -10000.0);
  maxSslope = configFile->GetValue("max_Sslope", 0.0);
  
  minX0 = configFile->GetValue("min_X0", -2000.0);
  maxX0 = configFile->GetValue("max_X0",  2000.0);
  minY0 = configFile->GetValue("min_Y0", -2000.0);
  maxY0 = configFile->GetValue("max_Y0",  2000.0);
  minZ0 = configFile->GetValue("min_Z0", -2000.0);
  maxZ0 = configFile->GetValue("max_Z0",  2000.0);
  
  maxEta                      = configFile->GetValue("max_eta",10.0);
  nTrackHits                  = configFile->GetValue("n_track_hits",3);
  minPx                       = configFile->GetValue("min_px",50.0);
  minPy                       = configFile->GetValue("min_py",50.0);
  minPz                       = configFile->GetValue("min_pz",50.0);
  maxPx                       = configFile->GetValue("max_px",250.0);
  maxPy                       = configFile->GetValue("max_py",250.0);
  maxPz                       = configFile->GetValue("max_pz",250.0);
  injectPionHits              = configFile->GetValue("inject_pion_hits",0);
  nTests                      = configFile->GetValue("n_tests",1);
  outputPath                  = configFile->GetValue("output_path","unnamed.root");
  toleranceX                  = configFile->GetValue("tolerance_x",10.0);
  toleranceY                  = configFile->GetValue("tolerance_y",10.0);
  toleranceZ                  = configFile->GetValue("tolerance_z",10.0);
  tolerancePx                 = configFile->GetValue("tolerance_px",30.0);
  tolerancePy                 = configFile->GetValue("tolerance_py",30.0);
  tolerancePz                 = configFile->GetValue("tolerance_pz",30.0);
  nNoiseHits                  = configFile->GetValue("n_noise_hits",500);
  nTrackerLayers              = configFile->GetValue("n_tracker_layers",4);
  
 
  performCutsLevel       = configFile->GetValue("cuts_level",2);
  
  category               = configFile->GetValue("analysis_category","");
  scanMETbinning         = configFile->GetValue("scan_MET_binning",0);
  doMETbinning           = configFile->GetValue("do_MET_binning",0);
  saveEvents             = configFile->GetValue("save_events",0);
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){runBackground.push_back(false);}
  for(int iSig=0;iSig<kNsignals;iSig++){runSignal.push_back(false);}
  for(int iData=0;iData<kNdata;iData++){runData.push_back(false);}

  runBackground[kQCD]       = configFile->GetValue("do_QCD",0);
  runBackground[kZmumuJets] = configFile->GetValue("do_Zmm",0);
  runBackground[kTT]        = configFile->GetValue("do_tops",0);
  runBackground[kVV]        = configFile->GetValue("do_dibosons",0);
  runBackground[kWmunuJets] = configFile->GetValue("do_Wmv",0);
  runBackground[kZnunuJets] = configFile->GetValue("do_Zvv",0);
  
  runSignal[kWino_M_300_cTau_3]   = configFile->GetValue("do_300_3",0);
  runSignal[kWino_M_300_cTau_10]  = configFile->GetValue("do_300_10",0);
  runSignal[kWino_M_300_cTau_30]  = configFile->GetValue("do_300_30",0);
  runSignal[kWino_M_500_cTau_10]  = configFile->GetValue("do_500_10",0);
  runSignal[kWino_M_500_cTau_20]  = configFile->GetValue("do_500_20",0);
  runSignal[kWino_M_650_cTau_10]  = configFile->GetValue("do_650_10",0);
  runSignal[kWino_M_650_cTau_20]  = configFile->GetValue("do_650_20",0);
  runSignal[kWino_M_800_cTau_10]  = configFile->GetValue("do_800_10",0);
  runSignal[kWino_M_800_cTau_20]  = configFile->GetValue("do_800_20",0);
  runSignal[kWino_M_1000_cTau_10] = configFile->GetValue("do_1000_10",0);
  runSignal[kWino_M_1000_cTau_20] = configFile->GetValue("do_1000_20",0);
  
  runData[kElectron_Run2017B] = configFile->GetValue("do_2017",0);
  
  printYields            = configFile->GetValue("print_yields",1);
  printBackgroundDetails = configFile->GetValue("print_background_details",0);
  printDataDetails       = configFile->GetValue("print_data_details",0);
  printSignalDetails     = configFile->GetValue("print_signal_details",0);
  
  drawStandardPlots      = configFile->GetValue("draw_standard_plots",1);
  drawPerLayerPlots      = configFile->GetValue("draw_per_layer_plots",0);
  showLegends            = configFile->GetValue("show_legends",1);
  
  
  
  maxNeventsBackground   = configFile->GetValue("max_N_events_background",-1);
  maxNeventsSignal       = configFile->GetValue("max_N_events_signal",-1);
  maxNeventsData         = configFile->GetValue("max_N_events_data",-1);
  
  loadFriendTree         = configFile->GetValue("load_friend_tree",0);
  
  totalLuminosity        = configFile->GetValue("total_luminosity",146.91);
  
  showGeometryPixel      = configFile->GetValue("show_geometry_pixel",0);
  showGeometryStrip      = configFile->GetValue("show_geometry_strip",0);
  showGeometryEcal       = configFile->GetValue("show_geometry_ecal",0);
  showGeometryHcal       = configFile->GetValue("show_geometry_hcal",0);
  
  drawTrackerClusters = configFile->GetValue("draw_tracker_clusters",1);
  drawMET             = configFile->GetValue("draw_met",1);
  drawJets            = configFile->GetValue("draw_jets",1);
  drawPionSimHits     = configFile->GetValue("draw_pion_simhits",1);
  drawPionClusters    = configFile->GetValue("draw_pion_clusters",1);
  drawCharginoSimHits = configFile->GetValue("draw_chargino_simhits",1);
}


void ConfigManager::Print()
{
  cout<<"\n\n========================================================"<<endl;
  cout<<"Basic configuration:"<<endl;
  cout<<"\tN tests:"<<nTests<<endl;
  cout<<"\tN noise hits:"<<nNoiseHits<<endl;
  cout<<"\tOutput file:"<<outputPath<<endl;
  cout<<"\tTrack eta:"<<maxEta<<"\t n hits:"<<nTrackHits<<endl;
  cout<<"\tN tracker layers:"<<nTrackerLayers<<endl;
  
  cout<<"Pion properties:"<<endl;
  cout<<"\tInjecting pion hits: "<<(injectPionHits ? "yes" : "no")<<endl;
  cout<<"\tpx:"<<minPx<<" - "<<maxPx<<endl;
  cout<<"\tpy:"<<minPy<<" - "<<maxPy<<endl;
  cout<<"\tpz:"<<minPz<<" - "<<maxPz<<endl;
  
  cout<<"Benchmark options:"<<endl;
  cout<<"\tTolerance X:"<<toleranceX<<endl;
  cout<<"\tTolerance Y:"<<toleranceY<<endl;
  cout<<"\tTolerance Z:"<<toleranceZ<<endl;
  cout<<"\tTolerance Px:"<<tolerancePx<<endl;
  cout<<"\tTolerance Py:"<<tolerancePy<<endl;
  cout<<"\tTolerance Pz:"<<tolerancePz<<endl;
  
  cout<<"========================================================\n\n"<<endl;
}
