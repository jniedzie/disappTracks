//
//  EventCut.cpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "EventCut.hpp"

#include <iostream>

using namespace std;

EventCut::EventCut(ECut cutType) :
minMetPt(0)
{
  switch (cutType) {
    case kEmpty:
      break;
    case kMet100GeV:
      minMetPt = 100.0;
      break;
    default:
      cout<<"ERROR -- no event cut specified... in case you want a blank cut to customize, use ECut::kEmpty."<<endl;
      exit(0);
      break;
  }
}

EventCut::~EventCut()
{
  
}
