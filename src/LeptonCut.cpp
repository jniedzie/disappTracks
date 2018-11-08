//
//  LeptonCut.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#include "LeptonCut.hpp"

LeptonCut::LeptonCut() :
pt(range<double>()),
relativeIsolation(range<double>()),
requireTightID(false)
{

}

LeptonCut::~LeptonCut()
{
  
}
