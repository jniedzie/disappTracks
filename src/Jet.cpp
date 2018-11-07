//
//  Jet.cpp
//
//  Created by Jeremi Niedziela on 17/07/2018.
//

#include "Jet.hpp"

#include <TMath.h>

#include <iostream>

using namespace std;

Jet::Jet() :
pt(99999),
eta(99999),
phi(99999),
mass(999999),
isForward(false)
{

};

Jet::~Jet()
{
  
}

void Jet::Print()
{
  cout<<"Jet pt:"<<pt<<"\teta:"<<eta<<"\tphi:"<<phi<<endl;
  cout<<"f_CH:"<<chargedHadronEnergyFraction<<"\tf_NE:"<<neutralHadronEnergyFraction<<endl;
  cout<<"forward:"<<isForward<<endl;
}

bool Jet::IsPassingCut(JetCut *cut)
{
  // check jet's pT
  if(pt < cut->GetMinPt() || pt > cut->GetMaxPt())return false;
  
  // check eta
  if(fabs(eta) > cut->GetMaxEta()) return false;
  if(isForward && fabs(eta) > cut->GetMaxEtaFwd())return false;
  
  // check hadron energy fractions
  if(chargedHadronEnergyFraction < cut->GetMinChargedHadronEnergyFraction()) return false;
  if(neutralHadronEnergyFraction > cut->GetMaxNeutralHadronEnergyFraction()) return false;
  
  return true;
}
