//  LeptonCut.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.

#include "LeptonCut.hpp"

LeptonCut::LeptonCut() :

pt(range<double>()),
relativeIsolation(range<double>()),
requiresTightID(false)
{

}

LeptonCut::LeptonCut(const LeptonCut &c)
{
  pt                = c.pt;
  relativeIsolation = c.relativeIsolation;
  requiresTightID   = c.requiresTightID;
}

LeptonCut::~LeptonCut()
{
  
}
