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
  // check jet's kinematical properties
  if(cut->GetPt().IsOutside(pt))  return false;
  if(!isForward && cut->GetEta().IsOutside(eta))  return false;
  if( isForward && cut->GetEtaForward().IsOutside(eta)) return false;
  
  // check hadron energy fractions
  if(cut->GetChargedHadronEnergyFraction().IsOutside(chargedHadronEnergyFraction)) return false;
  if(cut->GetNeutralHadronEnergyFraction().IsOutside(neutralHadronEnergyFraction)) return false;
  
  return true;
}
