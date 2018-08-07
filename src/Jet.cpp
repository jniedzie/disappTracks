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
pt(99999),
eta(99999),
phi(99999),
mass(999999)
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
  
  // check eta
  if(eta > cut->GetMaxEta()) return false;
  
  // check hadron energy fractions
  if(chargedHadronEnergyFraction < cut->GetMinChargedHadronEnergyFraction()) return false;
  
  if(neutralHadronEnergyFraction > cut->GetMaxNeutralHadronEnergyFraction()) return false;
  
  return true;
}
