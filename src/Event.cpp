//
//  Event.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "Event.hpp"

#include <TLorentzVector.h>

Event::Event() :
vertex(make_unique<Point>(0,0,0)),
trackProcessor(make_unique<TrackProcessor>())
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
  SetNlepton(e.nLepton);
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


void Event::LoadAdditionalInfo()
{
  string basePath;
  
  if(dataType == xtracks::kSignal)      basePath = inFileNameSignal[setIter];
  if(dataType == xtracks::kBackground)  basePath = inFileNameBackground[setIter][0];
  if(dataType == xtracks::kData)        basePath = inFileNameData[setIter][0];
  
  TFile *inFile = TFile::Open(Form("%s/tree_friend.root",basePath.c_str()));
  
  if(!inFile){
    cout<<"ERROR -- no file with all hits was found in path: "<<basePath<<"/tree_friend.root"<<endl;
    return;
  }
  TTree *tree = (TTree*)inFile->Get("CharginoAnalyzer/tree");
  
  if(!tree){
    cout<<"ERROR -- no tree with additional information was found"<<endl;
    return;
  }
  
  vector<double> *pionVx            = nullptr;
  vector<double> *pionVy            = nullptr;
  vector<double> *pionVz            = nullptr;
  vector<double> *pionPx            = nullptr;
  vector<double> *pionPy            = nullptr;
  vector<double> *pionPz            = nullptr;
  vector<int>    *pionCharge        = nullptr;
  vector<double> *pionSimHitsX      = nullptr;
  vector<double> *pionSimHitsY      = nullptr;
  vector<double> *pionSimHitsZ      = nullptr;
  vector<int>    *pionSimHitsSubDet = nullptr;
  
  vector<double> *charginoSimHitsX      = nullptr;
  vector<double> *charginoSimHitsY      = nullptr;
  vector<double> *charginoSimHitsZ      = nullptr;
  vector<int>    *charginoSimHitsSubDet = nullptr;
  
  vector<double> *pixelClusterX      = nullptr;
  vector<double> *pixelClusterY      = nullptr;
  vector<double> *pixelClusterZ      = nullptr;
  vector<double> *pixelClusterCharge = nullptr;
  vector<int>    *pixelClusterSubDet = nullptr;
  
  vector<double> *stripClusterX      = nullptr;
  vector<double> *stripClusterY      = nullptr;
  vector<double> *stripClusterZ      = nullptr;
  vector<double> *stripClusterXerr   = nullptr;
  vector<double> *stripClusterYerr   = nullptr;
  vector<double> *stripClusterZerr   = nullptr;
  vector<double> *stripClusterCharge = nullptr;
  vector<int>    *stripClusterSubDet = nullptr;
  
  vector<double> *pionClusterX      = nullptr;
  vector<double> *pionClusterY      = nullptr;
  vector<double> *pionClusterZ      = nullptr;
  vector<double> *pionClusterXerr   = nullptr;
  vector<double> *pionClusterYerr   = nullptr;
  vector<double> *pionClusterZerr   = nullptr;
  vector<double> *pionClusterCharge = nullptr;
  vector<int>    *pionClusterSubDet = nullptr;
  
  uint run;
  uint lumi;
  unsigned long long event;
  
  tree->SetBranchAddress("runNumber",&run);
  tree->SetBranchAddress("lumiBlock",&lumi);
  tree->SetBranchAddress("eventNumber",&event);
  
  tree->SetBranchAddress("pion_vx",&pionVx);
  tree->SetBranchAddress("pion_vy",&pionVy);
  tree->SetBranchAddress("pion_vz",&pionVz);
  tree->SetBranchAddress("pion_px",&pionPx);
  tree->SetBranchAddress("pion_py",&pionPy);
  tree->SetBranchAddress("pion_pz",&pionPz);
  tree->SetBranchAddress("pion_charge",&pionCharge);
  
  tree->SetBranchAddress("pion_simHits_x",&pionSimHitsX);
  tree->SetBranchAddress("pion_simHits_y",&pionSimHitsY);
  tree->SetBranchAddress("pion_simHits_z",&pionSimHitsZ);
  tree->SetBranchAddress("pion_simHits_subDet",&pionSimHitsSubDet);
  
  tree->SetBranchAddress("chargino_simHits_x",&charginoSimHitsX);
  tree->SetBranchAddress("chargino_simHits_y",&charginoSimHitsY);
  tree->SetBranchAddress("chargino_simHits_z",&charginoSimHitsZ);
  tree->SetBranchAddress("chargino_simHits_subDet",&charginoSimHitsSubDet);
  
  tree->SetBranchAddress("pixelCluster_x",&pixelClusterX);
  tree->SetBranchAddress("pixelCluster_y",&pixelClusterY);
  tree->SetBranchAddress("pixelCluster_z",&pixelClusterZ);
  tree->SetBranchAddress("pixelCluster_charge",&pixelClusterCharge);
  tree->SetBranchAddress("pixelCluster_subDet",&pixelClusterSubDet);
  
  tree->SetBranchAddress("stripCluster_x",&stripClusterX);
  tree->SetBranchAddress("stripCluster_y",&stripClusterY);
  tree->SetBranchAddress("stripCluster_z",&stripClusterZ);
  tree->SetBranchAddress("stripCluster_ex",&stripClusterXerr);
  tree->SetBranchAddress("stripCluster_ey",&stripClusterYerr);
  tree->SetBranchAddress("stripCluster_ez",&stripClusterZerr);
  tree->SetBranchAddress("stripCluster_charge",&stripClusterCharge);
  tree->SetBranchAddress("stripCluster_subDet",&stripClusterSubDet);
  
  tree->SetBranchAddress("pionCluster_x",&pionClusterX);
  tree->SetBranchAddress("pionCluster_y",&pionClusterY);
  tree->SetBranchAddress("pionCluster_z",&pionClusterZ);
  tree->SetBranchAddress("pionCluster_ex",&pionClusterXerr);
  tree->SetBranchAddress("pionCluster_ey",&pionClusterYerr);
  tree->SetBranchAddress("pionCluster_ez",&pionClusterZerr);
  tree->SetBranchAddress("pionCluster_charge",&pionClusterCharge);
  tree->SetBranchAddress("pionCluster_subDet",&pionClusterSubDet);
  
  bool eventFound = false;
  
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(run == runNumber && lumi == lumiSection && event == eventNumber){
      eventFound = true;
      break;
    }
  }
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find all hits for requested event!\n\n"<<endl;
    return;
  }

  for(uint i=0;i<pionVx->size();i++){
    // change units from cm to mm and from GeV to MeV
    auto helix = make_unique<Helix>(Point(10*pionVx->at(i),
                                          10*pionVy->at(i),
                                          10*pionVz->at(i)),
                                    make_unique<Point>(1000*pionPx->at(i),
                                                       1000*pionPy->at(i),
                                                       1000*pionPz->at(i)),
                                    1*pionCharge->at(i));
    
    genPionHelices.push_back(move(helix));
  }
  
  for(uint i=0;i<pionSimHitsX->size();i++){
    // convert cm to mm
    pionSimHits.push_back(make_shared<Point>(10*pionSimHitsX->at(i),
                                             10*pionSimHitsY->at(i),
                                             10*pionSimHitsZ->at(i),
                                             0,
                                             subDetMap[pionSimHitsSubDet->at(i)]));
  }
  
  for(uint i=0;i<charginoSimHitsX->size();i++){
    // convert cm to mm
    charginoSimHits.push_back(make_shared<Point>(10*charginoSimHitsX->at(i),
                                                 10*charginoSimHitsY->at(i),
                                                 10*charginoSimHitsZ->at(i),
                                                 0,
                                                 subDetMap[charginoSimHitsSubDet->at(i)]));
  }
  
  // Parameters for all hits in the pixel barrel
  for(uint i=0;i<pixelClusterX->size();i++){
    // convert cm to mm
    trackerClusters.push_back(make_shared<Point>(10*pixelClusterX->at(i),
                                                 10*pixelClusterY->at(i),
                                                 10*pixelClusterZ->at(i),
                                                 pixelClusterCharge->at(i),
                                                 subDetMap[pixelClusterSubDet->at(i)]));
  }
  
 
  
  for(uint i=0;i<stripClusterX->size();i++){
    // convert cm to mm
    trackerClusters.push_back(make_shared<Point>(10*stripClusterX->at(i),
                                                 10*stripClusterY->at(i),
                                                 10*stripClusterZ->at(i),
                                                 stripClusterCharge->at(i),
                                                 subDetMap[stripClusterSubDet->at(i)],
                                                 10*stripClusterXerr->at(i),
                                                 10*stripClusterYerr->at(i),
                                                 10*stripClusterZerr->at(i)));
    
    
  }
  
  for(uint i=0;i<pionClusterX->size();i++){
    // convert cm to mm
    pionClusters.push_back(make_shared<Point>(10*pionClusterX->at(i),
                                              10*pionClusterY->at(i),
                                              10*pionClusterZ->at(i),
                                              pionClusterCharge->at(i),
                                              subDetMap[pionClusterSubDet->at(i)],
                                              10*pionClusterXerr->at(i),
                                              10*pionClusterYerr->at(i),
                                              10*pionClusterZerr->at(i)));
  }
}
