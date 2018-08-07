//
//  JetCut.cpp
//
//  Created by Jeremi Niedziela on 20/07/2018.
//

#include "JetCut.hpp"

#include <iostream>

using namespace std;

JetCut::JetCut(ECut cutType) :
minPt(0.0),
maxPt(9999999),
maxEta(9999999),
minChargedHadronEnergyFraction(0),
maxNeutralHadronEnergyFraction(9999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kPt100GeV) minPt = 100.0;
  if(cutType&kPt200GeV) minPt = 200.0;
  if(cutType&kEta2p4)   maxEta = 2.4;
  if(cutType&kChHEF0p1) minChargedHadronEnergyFraction = 0.1;
  if(cutType&kNeHEF0p8) maxNeutralHadronEnergyFraction = 0.8;
}

JetCut::~JetCut()
{
  
}
