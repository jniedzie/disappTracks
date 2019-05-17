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
  
  helixThickness              = configFile->GetValue("helix_thickness",1.0);
  circleThickness             = configFile->GetValue("circle_thickness",1.0);
  linesToleranceForCircles    = configFile->GetValue("lines_tolerance_for_circles",1.0);
  linesToleranceForRegularity = configFile->GetValue("lines_tolerance_for_regularity",10.0);
  stepPz                      = configFile->GetValue("step_pz",1.0);
  zRegularityTolerance        = configFile->GetValue("z_regularity_tolerance",1.0);
  minNpointsAlongZ            = configFile->GetValue("min_n_points_along_z",2);
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
  
  cout<<"Algorithm parameters:"<<endl;
  cout<<"\tHelix thickness:"<<helixThickness<<endl;
  cout<<"\tCircle thickness:"<<circleThickness<<endl;
  cout<<"\tLines tolerance for circles:"<<linesToleranceForCircles<<endl;
  cout<<"\tLines tolerance for circles:"<<linesToleranceForRegularity<<endl;
  cout<<"\tStep z:"<<stepPz<<endl;
  cout<<"\tZ regularity tolerance:"<<zRegularityTolerance<<endl;
  cout<<"\tMin points along Z:"<<minNpointsAlongZ<<endl;
  
  cout<<"Benchmark options:"<<endl;
  cout<<"\tTolerance X:"<<toleranceX<<endl;
  cout<<"\tTolerance Y:"<<toleranceY<<endl;
  cout<<"\tTolerance Z:"<<toleranceZ<<endl;
  cout<<"\tTolerance Px:"<<tolerancePx<<endl;
  cout<<"\tTolerance Py:"<<tolerancePy<<endl;
  cout<<"\tTolerance Pz:"<<tolerancePz<<endl;
  
  cout<<"========================================================\n\n"<<endl;
}
