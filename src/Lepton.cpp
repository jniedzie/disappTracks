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
relativeIsolation(0)
{

};

bool Lepton::IsPassingCut(LeptonCut *cut)
{
  // check pt
  if(cut->GetPt().IsOutside(pt)) return false;

  // check isolation
  if(cut->GetRelativeIsolation().IsOutside(relativeIsolation)) return false;
  
  // check tight id requirement
  if(cut->RequiresTightID() && !tightID) return false;
  
  return true;
}



