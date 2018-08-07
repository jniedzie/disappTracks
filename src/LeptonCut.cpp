//
//  LeptonCut.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#include "LeptonCut.hpp"

LeptonCut::LeptonCut(ECut cutType) :
minPt(0.0),
maxPt(9999999),
maxIsolation(9999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kTightID)  requireTightID = true;
  if(cutType&kPt20GeV)  minPt = 20.0;
  if(cutType&kIsolated) maxIsolation = 0.15;
}

LeptonCut::~LeptonCut()
{
  
}
