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
maxPt(9999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kPt100GeV) minPt = 100.0;
  if(cutType&kPt200GeV) minPt = 200.0;
}

JetCut::~JetCut()
{
  
}
