//
//  FitterConfig.cpp
//
//  Created by Jeremi Niedziela on 09/01/2019.
//

#include "FitterConfig.hpp"

FitterConfig::FitterConfig(string _path)
{
  config = make_unique<TEnv>();
  
  if(config->ReadFile(_path.c_str(), kEnvUser) < 0){
    cout<<"ERROR - could not load config file:"<<_path<<endl;
    exit(0);
  }
  
  helixThickness = config->GetValue("helix_thickness",1.0);
  circleThickness = config->GetValue("circle_thickness",1.0);
  linesToleranceForCircles = config->GetValue("lines_tolerance_for_circles",1.0);
  linesToleranceForRegularity = config->GetValue("lines_tolerance_for_regularity",10.0);
  stepPz = config->GetValue("step_pz",1.0);
  zRegularityTolerance = config->GetValue("z_regularity_tolerance",1.0);
  minNpointsAlongZ = config->GetValue("min_n_points_along_z",2);
  maxEta = config->GetValue("max_eta",10.0);
  nTrackHits = config->GetValue("n_track_hits",3);
  minPx = config->GetValue("min_px",50.0);
  minPy = config->GetValue("min_py",50.0);
  minPz = config->GetValue("min_pz",50.0);
  maxPx = config->GetValue("max_px",250.0);
  maxPy = config->GetValue("max_py",250.0);
  maxPz = config->GetValue("max_pz",250.0);
  injectPionHits = config->GetValue("inject_pion_hits",0);
  nTests = config->GetValue("n_tests",1);
  outputPath = config->GetValue("output_path","unnamed.root");
  toleranceX = config->GetValue("tolerance_x",10.0);
  toleranceY = config->GetValue("tolerance_y",10.0);
  toleranceZ = config->GetValue("tolerance_z",10.0);
  tolerancePx = config->GetValue("tolerance_px",30.0);
  tolerancePy = config->GetValue("tolerance_py",30.0);
  tolerancePz = config->GetValue("tolerance_pz",30.0);
  nNoiseHits = config->GetValue("n_noise_hits",500);
}


void FitterConfig::Print()
{
  cout<<"\n\n========================================================"<<endl;
  cout<<"Basic configuration:"<<endl;
  cout<<"\tN tests:"<<nTests<<endl;
  cout<<"\tN noise hits:"<<nNoiseHits<<endl;
  cout<<"\tOutput file:"<<outputPath<<endl;
  cout<<"\tTrack eta:"<<maxEta<<"\t n hits:"<<config->GetValue("max_l",-1)<<endl;
  
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
