//  JetCut.cpp
//
//  Created by Jeremi Niedziela on 20/07/2018.

#include "JetCut.hpp"

JetCut::JetCut() :

pt(range<double>()),
eta(range<double>()),
etaForward(range<double>()),
ChHEF(range<double>()),
NeHEF(range<double>())
{

}

JetCut::JetCut(const JetCut &c)
{
  pt          = c.pt;
  eta         = c.eta;
  etaForward  = c.etaForward;
  ChHEF       = c.ChHEF;
  NeHEF       = c.NeHEF;
}

JetCut::~JetCut()
{
  
}
