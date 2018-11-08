//
//  Lepton.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#include "Lepton.hpp"

Lepton::Lepton() :
pt(99999),
eta(99999),
phi(99999),
relativeIsolation(0),
tightID(false),
pid(99999)
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



