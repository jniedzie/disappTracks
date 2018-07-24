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
minMetPt(0),
minNjets(0),
minNtracks(0)
{
  switch (cutType) {
    case kEmpty:
      break;
    case kOneTrack:
      minNtracks = 1;
      break;
    case kOneJet:
      minNjets = 1;
      break;
    case kOneTrackOneJet:
      minNtracks = 1;
      minNjets = 1;
      break;
    case kMet100GeV:
      minMetPt = 100.0;
      break;
    case kMet100GeVOneJet:
      minMetPt = 100.0;
      minNjets = 1;
      break;
    case kMet100GeVOneTrack:
      minMetPt = 100.0;
      minNtracks = 1;
      break;
    case kMet100GeVOneTrackOneJet:
      minMetPt = 100.0;
      minNtracks = 1;
      minNjets = 1;
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
