//  FitterConfig.cpp
//
//  Created by Jeremi Niedziela on 09/01/2019.

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
  
  
  // using a map
  
  params["verbosity_level"]             = configFile->GetValue("verbosity_level", 0);
  params["double_hit_max_distance"]     = configFile->GetValue("double_hit_max_distance", 0.0);
  params["do_asymmetric_constraints"]   = configFile->GetValue("do_asymmetric_constraints", 0);
  
  params["seed_max_chi2"]                 = configFile->GetValue("seed_max_chi2",100.0);
  params["seed_middle_hit_min_delta_phi"] = configFile->GetValue("seed_middle_hit_min_delta_phi",-3.14);
  params["seed_middle_hit_max_delta_phi"] = configFile->GetValue("seed_middle_hit_max_delta_phi", 3.14);
  params["seed_middle_hit_max_delta_z"]   = configFile->GetValue("seed_middle_hit_max_delta_z",1000.0);
  params["seed_last_hit_min_delta_phi"]   = configFile->GetValue("seed_last_hit_min_delta_phi",-3.14);
  params["seed_last_hit_max_delta_phi"]   = configFile->GetValue("seed_last_hit_max_delta_phi", 3.14);
  params["seed_last_hit_max_delta_z"]     = configFile->GetValue("seed_last_hit_max_delta_z",1000.0);

  params["track_max_chi2"]           = configFile->GetValue("track_max_chi2",100.0);
  params["next_point_min_delta_phi"] = configFile->GetValue("next_point_min_delta_phi",-3.14);
  params["next_point_max_delta_phi"] = configFile->GetValue("next_point_max_delta_phi", 3.14);
  params["next_point_max_delta_z"]   = configFile->GetValue("next_point_max_delta_z",1000.0);
  params["next_point_max_delta_xy"]  = configFile->GetValue("next_point_max_delta_xy",1000.0);
  params["next_point_max_delta_t"]   = configFile->GetValue("next_point_max_delta_t",1000.0);
  
  params["merging_max_different_point"] = configFile->GetValue("merging_max_different_point", 0);
  params["candidate_min_n_points"] = configFile->GetValue("candidate_min_n_points", 0);
  params["track_min_n_points"] = configFile->GetValue("track_min_n_points", 0);
  params["track_min_n_layers"] = configFile->GetValue("track_min_n_layers", 0);
  
  params["max_n_missing_hits"] = configFile->GetValue("max_n_missing_hits", 0);
  params["max_n_missing_hits_in_raw"] = configFile->GetValue("max_n_missing_hits_in_raw", 0);
  params["merge_at_turn_back"] = configFile->GetValue("merge_at_turn_back", 0);
  params["merge_final_helices"] = configFile->GetValue("merge_final_helices", 0);
  params["allow_turning_back"] = configFile->GetValue("allow_turning_back", 0);
  params["require_good_starting_values"] = configFile->GetValue("require_good_starting_values", 0);
  params["exp_radius_function"] = configFile->GetValue("exp_radius_function", 0);
  params["exp_slope_function"] = configFile->GetValue("exp_slope_function", 0);
  
  params["allow_one_less_layer"] = configFile->GetValue("allow_one_less_layer", 0);
  params["allow_one_more_layer"] = configFile->GetValue("allow_one_more_layer", 0);
  params["check_opposite_charge_below_Nlayers"] = configFile->GetValue("check_opposite_charge_below_Nlayers", 0);
  params["min_layers_for_delta_xy"] = configFile->GetValue("min_layers_for_delta_xy", inf);
  
  params["start_R0"] = configFile->GetValue("start_R0", 320.0);
  params["min_R0"] = configFile->GetValue("min_R0", 50.0);
  params["max_R0"] = configFile->GetValue("max_R0", 1000.0);
  params["min_Rslope"] = configFile->GetValue("min_Rslope", 0.0);
  params["max_Rslope"] = configFile->GetValue("max_Rslope", 10000.0);
  params["min_S0"] = configFile->GetValue("min_S0", -10000.0);
  params["max_S0"] = configFile->GetValue("max_S0", 10000.0);
  params["min_Sslope"] = configFile->GetValue("min_Sslope", -10000.0);
  params["max_Sslope"] = configFile->GetValue("max_Sslope", 0.0);
  
  params["min_X0"] = configFile->GetValue("min_X0", -2000.0);
  params["max_X0"] = configFile->GetValue("max_X0",  2000.0);
  params["min_Y0"] = configFile->GetValue("min_Y0", -2000.0);
  params["max_Y0"] = configFile->GetValue("max_Y0",  2000.0);
  params["min_Z0"] = configFile->GetValue("min_Z0", -2000.0);
  params["max_Z0"] = configFile->GetValue("max_Z0",  2000.0);
  
  params["max_eta"] = configFile->GetValue("max_eta",10.0);
  params["n_track_hits"] = configFile->GetValue("n_track_hits",3);
  params["min_px"] = configFile->GetValue("min_px",50.0);
  params["min_py"] = configFile->GetValue("min_py",50.0);
  params["min_pz"] = configFile->GetValue("min_pz",50.0);
  params["max_px"] = configFile->GetValue("max_px",250.0);
  params["max_py"] = configFile->GetValue("max_py",250.0);
  params["max_pz"] = configFile->GetValue("max_pz",250.0);
  params["inject_pion_hits"] = configFile->GetValue("inject_pion_hits",0);
  params["n_tests"] = configFile->GetValue("n_tests",1);
  
  params["tolerance_x"] = configFile->GetValue("tolerance_x",10.0);
  params["tolerance_y"] = configFile->GetValue("tolerance_y",10.0);
  params["tolerance_z"] = configFile->GetValue("tolerance_z",10.0);
  params["tolerance_px"] = configFile->GetValue("tolerance_px",30.0);
  params["tolerance_py"] = configFile->GetValue("tolerance_py",30.0);
  params["tolerance_pz"] = configFile->GetValue("tolerance_pz",30.0);
  params["n_noise_hits"] = configFile->GetValue("n_noise_hits",500);
  params["n_tracker_layers"] = configFile->GetValue("n_tracker_layers",4);
  
  
  params["cuts_level"] = configFile->GetValue("cuts_level",2);
  
  params["scan_MET_binning"] = configFile->GetValue("scan_MET_binning",0);
  params["do_MET_binning"] = configFile->GetValue("do_MET_binning",0);
  params["save_events"] = configFile->GetValue("save_events",0);
  
  params["print_yields"]             = configFile->GetValue("print_yields",1);
  params["print_background_details"] = configFile->GetValue("print_background_details",0);
  params["print_data_details"]       = configFile->GetValue("print_data_details",0);
  params["print_signal_details"]     = configFile->GetValue("print_signal_details",0);
  
  params["draw_standard_plots"]      = configFile->GetValue("draw_standard_plots",1);
  params["draw_per_layer_plots"]     = configFile->GetValue("draw_per_layer_plots",0);
  params["show_legends"]             = configFile->GetValue("show_legends",1);
  
  params["max_N_events_background"] = configFile->GetValue("max_N_events_background",-1);
  params["max_N_events_signal"]     = configFile->GetValue("max_N_events_signal",-1);
  params["max_N_events_data"]       = configFile->GetValue("max_N_events_data",-1);
  
  params["load_friend_tree"]        = configFile->GetValue("load_friend_tree",0);
  
  params["total_luminosity"]        = configFile->GetValue("total_luminosity",146.91);
  
  params["show_geometry_pixel"]     = configFile->GetValue("show_geometry_pixel",0);
  params["show_geometry_strip"]     = configFile->GetValue("show_geometry_strip",0);
  params["show_geometry_ecal"]      = configFile->GetValue("show_geometry_ecal",0);
  params["show_geometry_hcal"]      = configFile->GetValue("show_geometry_hcal",0);
  
  params["draw_tracker_clusters"]   = configFile->GetValue("draw_tracker_clusters",1);
  params["draw_met"]                = configFile->GetValue("draw_met",1);
  params["draw_jets"]               = configFile->GetValue("draw_jets",1);
  params["draw_pion_simhits"]       = configFile->GetValue("draw_pion_simhits",1);
  params["draw_pion_clusters"]      = configFile->GetValue("draw_pion_clusters",1);
  params["draw_chargino_simhits"]   = configFile->GetValue("draw_chargino_simhits",1);
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
