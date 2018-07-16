//
//  Track.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "Track.hpp"

Track::Track()
{
  isShort = false;
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    dedx.push_back(0.0);
    subDetId.push_back(-1);
  }
};


int Track::GetNclusters()
{
  int nClusters=0;
  for(float d : dedx){
    if(d > 0.000001) nClusters++;
  }
  return nClusters;
}

void Track::Print()
{
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    cout<<"Layer:"<<iLayer<<"\tsub-det ID:"<<subDetId[iLayer]<<"\tdEdx:"<<dedx[iLayer]<<endl;
  }
}
