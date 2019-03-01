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
  vertex = make_unique<Point>(e.vertex);
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
  cout<<"\t MET pT:"<<metPt<<"\teta:"<<metEta<<"\tphi:"<<metPhi<<endl;
  cout<<"weigth:"<<genWeight<<endl;
  cout<<"xsec:"<<xsec<<endl;
  
  cout<<"\nTracks:"<<endl;
  for(auto t : tracks){ t->Print(); }
  
  cout<<"\nJets:"<<endl;
  for(auto j : jets){   j->Print(); }
  
  cout<<"================================================\n\n"<<endl;
}


shared_ptr<vector<Point>> Event::GetTrackerHits()
{
  shared_ptr<vector<Point>> trackerPoints = make_shared<vector<Point>>();
  
  string basePath;
  
  if(dataType == xtracks::kSignal)      basePath = inFileNameSignal[setIter];
  if(dataType == xtracks::kBackground)  basePath = inFileNameBackground[setIter][0];
  if(dataType == xtracks::kData)        basePath = inFileNameData[setIter][0];
  
  TFile *inFile = TFile::Open(Form("%s/charginoAllHits.root",basePath.c_str()));
  
  if(!inFile){
    cout<<"ERROR -- no file with all hits was found"<<endl;
    return trackerPoints;
  }
  TTree *tree = (TTree*)inFile->Get("hitsExtractor/hits");
  
  if(!tree){
    cout<<"ERROR -- no tree with all hits was found"<<endl;
    return trackerPoints;
  }
  
  vector<double> *hitX = nullptr;
  vector<double> *hitY = nullptr;
  vector<double> *hitZ = nullptr;
  vector<double> *hitCharge = nullptr;
  vector<double> *hitSizeX = nullptr;
  vector<double> *hitSizeY = nullptr;
  vector<double> *stripX = nullptr;
  vector<double> *stripY = nullptr;
  vector<double> *stripZ = nullptr;
  vector<double> *stripCharge = nullptr;
  
  uint run;
  uint lumi;
  unsigned long long event;
  
  tree->SetBranchAddress("runNumber",&run);
  tree->SetBranchAddress("lumiBlock",&lumi);
  tree->SetBranchAddress("eventNumber",&event);
  
  tree->SetBranchAddress("hitX",&hitX);
  tree->SetBranchAddress("hitY",&hitY);
  tree->SetBranchAddress("hitZ",&hitZ);
  tree->SetBranchAddress("hitCharge",&hitCharge);
  tree->SetBranchAddress("hitSizeX",&hitSizeX);
  tree->SetBranchAddress("hitSizeY",&hitSizeY);
  tree->SetBranchAddress("stripX",&stripX);
  tree->SetBranchAddress("stripY",&stripY);
  tree->SetBranchAddress("stripZ",&stripZ);
  tree->SetBranchAddress("stripCharge",&stripCharge);
  
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
    return trackerPoints;
  }

  // Parameters for all hits in the pixel barrel
  const double chargeThreshold = 0; // 2000, 5000, 25000
  const double minClusterSize = 0;
  const double maxClusterSize = inf;
  
  for(uint i=0;i<hitX->size();i++){
    if(hitCharge->at(i) < chargeThreshold) continue;
    double clusterSize = sqrt(pow(hitSizeX->at(i),2)+pow(hitSizeY->at(i),2));
    if(clusterSize < minClusterSize || clusterSize > maxClusterSize) continue;
    // convert cm to mm
    trackerPoints->push_back(Point(10*hitX->at(i),10*hitY->at(i),10*hitZ->at(i),hitCharge->at(i)));
  }
  return trackerPoints;
}


shared_ptr<vector<unique_ptr<Helix>>> Event::GetTruePionHelices()
{
  auto trueHelices = make_shared<vector<unique_ptr<Helix>>>();
  
  string basePath;
  
  if(dataType == xtracks::kSignal)      basePath = inFileNameSignal[setIter];
  if(dataType == xtracks::kBackground)  basePath = inFileNameBackground[setIter][0];
  if(dataType == xtracks::kData)        basePath = inFileNameData[setIter][0];
  
  TFile *inFile = TFile::Open(Form("%s/charginoAnalysisPions.root",basePath.c_str()));
  
  if(!inFile){
    cout<<"ERROR -- no file with true pion helices was found"<<endl;
    return trueHelices;
  }
  TTree *tree = (TTree*)inFile->Get("CharginoAnalyzer/pions");
  
  if(!tree){
    cout<<"ERROR -- no tree with gen-level pion helices was found"<<endl;
    return trueHelices;
  }
  
  vector<double> *pionVx = nullptr;
  vector<double> *pionVy = nullptr;
  vector<double> *pionVz = nullptr;
  vector<double> *pionPx = nullptr;
  vector<double> *pionPy = nullptr;
  vector<double> *pionPz = nullptr;
  vector<int> *pionCharge = nullptr;
  
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
  
  bool eventFound = false;
  
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(run == runNumber && lumi == lumiSection && event == eventNumber){
      eventFound = true;
      break;
    }
  }
  
  if(!eventFound){
    cout<<"\n\nERROR - could not find gen pions for requested event!\n\n"<<endl;
    return trueHelices;
  }
  
  
  for(uint i=0;i<pionVx->size();i++){
    // change units from cm to mm and from GeV to MeV
    auto helix = make_unique<Helix>(make_unique<Point>(10*pionVx->at(i), 10*pionVy->at(i), 10*pionVz->at(i)),
                                    make_unique<Point>(1000*pionPx->at(i), 1000*pionPy->at(i), 1000*pionPz->at(i)),
                                    pionCharge->at(i));
    
    trueHelices->push_back(move(helix));
  }
  return trueHelices;
}
