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
  
  dataType = e.dataType;
  setIter = e.setIter;
  lumiSection = e.lumiSection;
  runNumber = e.runNumber;
  
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
  
  hasFriendData     = e.hasFriendData;
  wasTagged         = e.wasTagged;
  trackerClusters   = e.trackerClusters;
  pionClusters      = e.pionClusters;
  pionSimHits       = e.pionSimHits;
  charginoSimHits   = e.charginoSimHits;
  genPionHelices    = e.genPionHelices;
  genCharginoTrack  = e.genCharginoTrack;
  friendTree        = e.friendTree;
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
  
  map<string, double> noiseThreshold = {
    { "TIB" , 70.0}, { "TOB" , 120.0}, { "TID" , 65.0}, { "TEC" , 80.0},
  };
  
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
    
    for(auto &pionCluster : pionClusters){
      if(*pionCluster == *point){
        point->SetIsPionHit(true);
        break;
      }
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
    
    if(config.params["cut_noise_hits"]){
      
      if(noiseThreshold.find(point->GetSubDetName()) != noiseThreshold.end()){
        if(point->GetValue() < noiseThreshold[point->GetSubDetName()]) continue;
      }
    }
    resultClusters.push_back(point);
  }
  
  return resultClusters;
}
