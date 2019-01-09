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
}

double FitterConfig::GetHelixThickness()
{
  return config->GetValue("helix_thickness",1.0);
}

double FitterConfig::GetCircleThickness()
{
  return config->GetValue("circle_thickness",1.0);
}

double FitterConfig::GetLinesTolerance()
{
  return config->GetValue("lines_tolerance",1.0);
}

double FitterConfig::GetStepPz()
{
  return config->GetValue("step_pz",1.0);
}

double FitterConfig::GetZregularityTolerance()
{
  return config->GetValue("z_regularity_tolerance",1.0);
}

int FitterConfig::GetMinPointsAlongZ()
{
  return config->GetValue("min_n_points_along_z",2);
}

double FitterConfig::GetMaxTrackEta()
{
  return config->GetValue("max_eta",10.0);
}

double FitterConfig::GetMinL()
{
  return layerR[config->GetValue("min_l",2)];
}

double FitterConfig::GetMaxL()
{
  return layerR[config->GetValue("max_l",3)];
}

double FitterConfig::GetMinPx()
{
  return config->GetValue("min_px",50.0);
}

double FitterConfig::GetMinPy()
{
  return config->GetValue("min_py",50.0);
}

double FitterConfig::GetMinPz()
{
  return config->GetValue("min_pz",50.0);
}

double FitterConfig::GetMaxPx()
{
  return config->GetValue("max_px",250.0);
}

double FitterConfig::GetMaxPy()
{
  return config->GetValue("max_py",250.0);
}

double FitterConfig::GetMaxPz()
{
  return config->GetValue("max_pz",250.0);
}

bool FitterConfig::GetInjectPionHits()
{
  return config->GetValue("inject_pion_hits",0);
}

int FitterConfig::GetNtests()
{
  return config->GetValue("n_tests",1);
}

const char* FitterConfig::GetOutputPath()
{
  return config->GetValue("output_path","unnamed.root");
}

double FitterConfig::GetToleranceX()
{
  return config->GetValue("tolerance_x",10.0);
}

double FitterConfig::GetToleranceY()
{
  return config->GetValue("tolerance_y",10.0);
}

double FitterConfig::GetToleranceZ()
{
  return config->GetValue("tolerance_z",10.0);
}

double FitterConfig::GetTolerancePx()
{
  return config->GetValue("tolerance_px",30.0);
}

double FitterConfig::GetTolerancePy()
{
  return config->GetValue("tolerance_py",30.0);
}

double FitterConfig::GetTolerancePz()
{
  return config->GetValue("tolerance_pz",30.0);
}

int FitterConfig::GetNnoiseHits()
{
  return config->GetValue("n_noise_hits",500);
}
