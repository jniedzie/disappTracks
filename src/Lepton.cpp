//
//  Lepton.cpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#include "Lepton.hpp"

Lepton::Lepton() :
eta(99999),
phi(99999),
pt(99999),
pid(99999)
{

};

bool Lepton::IsPassingCut(/*LeptonCut *cut*/)
{
  // check number of dedx clusters
//  if(GetNclusters() < cut->GetMinDedxClusters() || GetNclusters() > cut->GetMaxDedxClusters()){
//    return false;
//  }
//
//  // check values of dedx along the Lepton
//  if(GetTotalDedx() < cut->GetMinTotalDedx() || GetTotalDedx() > cut->GetMaxTotalDedx()){
//    return false;
//  }
//
//  for(int iCluster=0;iCluster<GetNclusters();iCluster++){
//    if(dedx[iCluster] < cut->GetMinDedxPerCluster()) return false;
//  }
//
//  // check pt
//  if(pt < cut->GetMinPt()) return false;
//
//  // check calo energy
//  if(caloEmEnergy > cut->GetMaxEmCalo()) return false;
//  if(caloHadEnergy > cut->GetMaxHadCalo()) return false;
//
//  // check eta
//  if(fabs(eta) > cut->GetMaxEta()) return false;
//
  return true;
}



