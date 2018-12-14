//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter(int _nPar) : nPar(_nPar)
{
  fitter = new ROOT::Fit::Fitter();
}

Fitter::~Fitter()
{
  
}

void Fitter::SetParameter(int i, string name, double start, double min, double max, bool fix)
{
  fitter->Config().ParSettings(i).SetName(name);
  fitter->Config().ParSettings(i).SetValue(start);
  fitter->Config().ParSettings(i).SetLimits((min < max) ? min : max,(min < max) ? max : min);
  if(fix) fitter->Config().ParSettings(i).Fix();
}

void Fitter::FixParameter(int i, string name, double val)
{
  fitter->Config().ParSettings(i).SetName(name);
  fitter->Config().ParSettings(i).SetValue(val);
  fitter->Config().ParSettings(i).Fix();
}

bool Fitter::RunFitting()
{
  return fitter->FitFCN();
}

const ROOT::Fit::FitResult& Fitter::GetResult()
{
  return fitter->Result();
}
