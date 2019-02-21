//
//  FitterConfig.cpp
//
//  Created by Jeremi Niedziela on 09/01/2019.
//

#include "ConfigManager.hpp"

unique_ptr<ConfigManager> config = nullptr;

ConfigManager::ConfigManager(string _path)
{
  config = make_unique<TEnv>();
  
  if(config->ReadFile(_path.c_str(), kEnvUser) < 0){
    cout<<"ERROR - could not load config file:"<<_path<<endl;
    cout<<"Will use default values"<<endl;
    return;
  }
  
  helixThickness              = config->GetValue("helix_thickness",1.0);
  circleThickness             = config->GetValue("circle_thickness",1.0);
  linesToleranceForCircles    = config->GetValue("lines_tolerance_for_circles",1.0);
  linesToleranceForRegularity = config->GetValue("lines_tolerance_for_regularity",10.0);
  stepPz                      = config->GetValue("step_pz",1.0);
  zRegularityTolerance        = config->GetValue("z_regularity_tolerance",1.0);
  minNpointsAlongZ            = config->GetValue("min_n_points_along_z",2);
  maxEta                      = config->GetValue("max_eta",10.0);
  nTrackHits                  = config->GetValue("n_track_hits",3);
  minPx                       = config->GetValue("min_px",50.0);
  minPy                       = config->GetValue("min_py",50.0);
  minPz                       = config->GetValue("min_pz",50.0);
  maxPx                       = config->GetValue("max_px",250.0);
  maxPy                       = config->GetValue("max_py",250.0);
  maxPz                       = config->GetValue("max_pz",250.0);
  injectPionHits              = config->GetValue("inject_pion_hits",0);
  nTests                      = config->GetValue("n_tests",1);
  outputPath                  = config->GetValue("output_path","unnamed.root");
  toleranceX                  = config->GetValue("tolerance_x",10.0);
  toleranceY                  = config->GetValue("tolerance_y",10.0);
  toleranceZ                  = config->GetValue("tolerance_z",10.0);
  tolerancePx                 = config->GetValue("tolerance_px",30.0);
  tolerancePy                 = config->GetValue("tolerance_py",30.0);
  tolerancePz                 = config->GetValue("tolerance_pz",30.0);
  nNoiseHits                  = config->GetValue("n_noise_hits",500);
  nTrackerLayers              = config->GetValue("n_tracker_layers",4);
  
 
  performCutsLevel       = config->GetValue("cuts_level",2);
  
  category               = config->GetValue("analysis_category","");
  scanMETbinning         = config->GetValue("scan_MET_binning",0);
  doMETbinning           = config->GetValue("do_MET_binning",0);
  saveEvents             = config->GetValue("save_events",0);
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){runBackground.push_back(false);}
  for(int iSig=0;iSig<kNsignals;iSig++){runSignal.push_back(false);}
  for(int iData=0;iData<kNdata;iData++){runData.push_back(false);}

  runBackground[kQCD]       = config->GetValue("do_QCD",0);
  runBackground[kZmumuJets] = config->GetValue("do_Zmm",0);
  runBackground[kTT]        = config->GetValue("do_tops",0);
  runBackground[kVV]        = config->GetValue("do_dibosons",0);
  runBackground[kWmunuJets] = config->GetValue("do_Wmv",0);
  runBackground[kZnunuJets] = config->GetValue("do_Zvv",0);
  
  runSignal[kWino_M_300_cTau_3]   = config->GetValue("do_300_3",0);
  runSignal[kWino_M_300_cTau_10]  = config->GetValue("do_300_10",0);
  runSignal[kWino_M_300_cTau_30]  = config->GetValue("do_300_30",0);
  runSignal[kWino_M_500_cTau_10]  = config->GetValue("do_500_10",0);
  runSignal[kWino_M_500_cTau_20]  = config->GetValue("do_500_20",0);
  runSignal[kWino_M_650_cTau_10]  = config->GetValue("do_650_10",0);
  runSignal[kWino_M_650_cTau_20]  = config->GetValue("do_650_20",0);
  runSignal[kWino_M_800_cTau_10]  = config->GetValue("do_800_10",0);
  runSignal[kWino_M_800_cTau_20]  = config->GetValue("do_800_20",0);
  runSignal[kWino_M_1000_cTau_10] = config->GetValue("do_1000_10",0);
  runSignal[kWino_M_1000_cTau_20] = config->GetValue("do_1000_20",0);
  
  runData[kElectron_Run2017B] = config->GetValue("do_2017",0);
  
  printYields            = config->GetValue("print_yields",1);
  printBackgroundDetails = config->GetValue("print_background_details",0);
  printDataDetails       = config->GetValue("print_data_details",0);
  printSignalDetails     = config->GetValue("print_signal_details",0);
  
  drawStandardPlots      = config->GetValue("draw_standard_plots",1);
  drawPerLayerPlots      = config->GetValue("draw_per_layer_plots",0);
  showLegends            = config->GetValue("show_legends",1);
  
  
  
  maxNeventsBackground   = config->GetValue("max_N_events_background",-1);
  maxNeventsSignal       = config->GetValue("max_N_events_signal",-1);
  maxNeventsData         = config->GetValue("max_N_events_data",-1);
  
  totalLuminosity        = config->GetValue("total_luminosity",146.91);
  
  showGeometryPixel      = config->GetValue("show_geometry_pixel",0);
  showGeometryStrip      = config->GetValue("show_geometry_strip",0);
  showGeometryEcal       = config->GetValue("show_geometry_ecal",0);
  showGeometryHcal       = config->GetValue("show_geometry_hcal",0);
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
