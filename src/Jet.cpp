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
pt(inf),
eta(inf),
phi(inf),
mass(inf),
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
