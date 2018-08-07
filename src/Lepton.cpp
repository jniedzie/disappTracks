//
//  Lepton.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#include "Lepton.hpp"

Lepton::Lepton() :
eta(99999),
phi(99999),
pt(99999),
pid(99999),
isolation(0)
{

};

bool Lepton::IsPassingCut(LeptonCut *cut)
{
  // check pt
  if(pt < cut->GetMinPt() || pt > cut->GetMaxPt()) return false;

  // check isolation
  if(isolation >= cut->GetMaxIsolation()) return false;
  
  // check tight id requirement
  if(cut->RequiresTightID() && !tightID) return false;
  
  return true;
}



