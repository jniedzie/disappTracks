//
//  Jet.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 17/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Jet.hpp"

#include <iostream>

using namespace std;

Jet::Jet() :
pt(0),
eta(0),
phi(0)
{

};

Jet::~Jet()
{
  
}

void Jet::Print()
{
  cout<<"Jet pt:"<<pt<<"\teta:"<<eta<<"\tphi:"<<phi<<endl;
}

bool Jet::IsPassingCut(JetCut *cut)
{
  // check jet's pT
  if(pt < cut->GetMinPt() || pt > cut->GetMaxPt()){
    return false;
  }
  
  return true;
}
