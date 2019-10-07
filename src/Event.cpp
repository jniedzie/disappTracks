//  Event.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.

#include "Event.hpp"

#include <TLorentzVector.h>

Event::Event() :
vertex(make_unique<Point>(0,0,0)),
hasFriendData(false),
wasTagged(false)
{

}

Event::Event(const Event &e)
{
  for(auto t : e.tracks){ tracks.push_back(t);}
  for(auto j : e.jets){   jets.push_back(j);}
  for(auto l : e.leptons){leptons.push_back(l);}
  for(auto h : e.helices){helices.push_back(h);}
  
  SetWeight(e.weight);
  SetNvertices(e.nVertices);
  vertex = make_unique<Point>(*e.vertex);
  SetNjet30(e.nJet30);
  SetNjet30a(e.nJet30a);
  SetNtau(e.nTau);
  
  SetMetSumEt(e.metSumEt);
  SetMetPt(e.metPt);
  SetMetMass(e.metMass);
  SetMetPhi(e.metPhi);
  SetMetEta(e.metEta);
  
  SetMetNoMuPt(e.metNoMuPt);
  SetMetNoMuMass(e.metNoMuMass);
  SetMetNoMuPhi(e.metNoMuPhi);
  SetMetNoMuEta(e.metNoMuEta);
  SetHasNoMuTrigger(e.metNoMuTrigger);
  
  SetGoodVerticesFlag(e.flag_goodVertices);
  SetBadPFmuonFlag(e.flag_badPFmuon);
  SetHBHEnoiseFlag(e.flag_HBHEnoise);
  SetHBHEnoiseIsoFlag(e.flag_HBHEnoiseIso);
  SetEcalDeadCellFlag(e.flag_EcalDeadCell);
  SetEeBadScFlag(e.flag_eeBadSc);
  SetBadChargedCandidateFlag(e.flag_badChargedCandidate);
  SetEcalBadCalibFlag(e.flag_ecalBadCalib);
  SetGlobalTightHalo2016Flag(e.flag_globalTightHalo2016);
  
  SetNgenChargino(e.nGenChargino);
  SetXsec(e.xsec);
  SetWgtSum(e.wgtsum);
  SetGenWeight(e.genWeight);
}


Event::~Event()
{
  
}

void Event::Print(){
  cout<<"\n\n================================================"<<endl;
  cout<<"Event:"<<endl;
  cout<<"\t lumi section:"<<lumiSection<<"\trun number:"<<runNumber<<"\tevent number:"<<eventNumber<<endl;
  cout<<"\t n vertices:"<<nVertices<<endl;
  cout<<"\t primary vertex:("<<vertex->GetX()<<", "<<vertex->GetY()<<", "<<vertex->GetZ()<<")"<<endl;
  cout<<"\t MET pT:"<<metPt<<"\teta:"<<metEta<<"\tphi:"<<metPhi<<endl;
  cout<<"weigth:"<<genWeight<<endl;
  cout<<"xsec:"<<xsec<<endl;
  
  cout<<"\nTracks:"<<endl;
  for(auto t : tracks){ t->Print(); }
  
  cout<<"\nJets:"<<endl;
  for(auto j : jets){   j->Print(); }
  
  cout<<"================================================\n\n"<<endl;
}

Points Event::GetClusters()
{
  Points resultClusters;
  
  for(auto &point : trackerClusters){
    if(!config.params["include_endcaps"] &&
       (point->GetSubDetName() == "TID" ||
        point->GetSubDetName() == "TEC" ||
        point->GetSubDetName() == "P1PXEC")) continue;
    
    if(point->GetSubDetName() != "TIB" &&
       point->GetSubDetName() != "TOB" &&
       point->GetSubDetName() != "P1PXB" &&
       point->GetSubDetName() != "TID" &&
       point->GetSubDetName() != "TEC" &&
       point->GetSubDetName() != "P1PXEC"){
      cout<<"Unknown detector:"<<point->GetSubDetName()<<endl;
    }
    
    if(config.params["fit_noise_clusters_only"]){
      bool isPionHit = false;
      for(auto &pionCluster : pionClusters){
        if(*pionCluster == *point){
          isPionHit = true;
          break;
        }
      }
      if(isPionHit) continue;
    }
    
    resultClusters.push_back(point);
  }
  
  return resultClusters;
}
