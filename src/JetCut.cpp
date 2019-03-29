//
//  JetCut.cpp
//
//  Created by Jeremi Niedziela on 20/07/2018.
//

#include "JetCut.hpp"

JetCut::JetCut() :
pt(range<double>()),
eta(range<double>()),
etaForward(range<double>()),
chargedHadronEnergyFraction(range<double>()),
neutralHadronEnergyFraction(range<double>()),
trackDeltaR(range<double>())
{

}

JetCut::JetCut(const JetCut &c)
{
  pt                          = c.pt;
  eta                         = c.eta;
  etaForward                  = c.etaForward;
  chargedHadronEnergyFraction = c.chargedHadronEnergyFraction;
  neutralHadronEnergyFraction = c.neutralHadronEnergyFraction;
  trackDeltaR                 = c.trackDeltaR;
}

JetCut::~JetCut()
{
  
}
