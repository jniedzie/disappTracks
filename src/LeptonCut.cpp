//
//  LeptonCut.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#include "LeptonCut.hpp"

LeptonCut::LeptonCut(ECut cutType) :
pt(range<double>()),
relativeIsolation(range<double>()),
requireTightID(false)
{
  if(cutType&kEmpty)    return;
  if(cutType&kTightID)  requireTightID = true;
  if(cutType&kPt20GeV)  pt = range<double>(20.0, 999999);
  if(cutType&kIsolated) relativeIsolation = range<double>(-999999, 0.15);
}

LeptonCut::~LeptonCut()
{
  
}
