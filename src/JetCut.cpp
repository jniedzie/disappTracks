//
//  JetCut.cpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 20/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "JetCut.hpp"

#include <iostream>

using namespace std;

JetCut::JetCut(ECut cutType) :
minPt(0.0),
maxPt(9999999)
{
  switch (cutType) {
    case kEmpty:
      break;
    case kHighPt:
      minPt = 100.0;
      break;
    default:
      cout<<"ERROR -- no jet cut specified... in case you want a blank cut to customize, use ECut::kEmpty."<<endl;
      exit(0);
      break;
  }
}

JetCut::~JetCut()
{
  
}
