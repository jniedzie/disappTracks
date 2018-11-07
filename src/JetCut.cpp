//
//  JetCut.cpp
//
//  Created by Jeremi Niedziela on 20/07/2018.
//

#include "JetCut.hpp"

JetCut::JetCut(ECut cutType) :
pt(range<double>()),
eta(range<double>()),
etaForward(range<double>()),
chargedHadronEnergyFraction(range<double>()),
neutralHadronEnergyFraction(range<double>()),
trackDeltaR(range<double>())
{
  if(cutType&kEmpty)      return;
  if(cutType&kPt30GeV)    pt = range<double>(30.0, 999999);
  if(cutType&kPt100GeV)   pt = range<double>(100.0, 999999);
  if(cutType&kPt200GeV)   pt = range<double>(200.0, 999999);
  if(cutType&kEta2p4)     eta = range<double>(-2.4, 2.4);
  if(cutType&kFwdEta4p7)  etaForward = range<double>(-4.7, 4.7);
  if(cutType&kChHEF0p1)   chargedHadronEnergyFraction = range<double>(0.1, 999999);
  if(cutType&kNeHEF0p8)   neutralHadronEnergyFraction = range<double>(-999999, 0.8);
}

JetCut::~JetCut()
{
  
}
