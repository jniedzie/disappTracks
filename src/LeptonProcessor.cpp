//
//  LeptonProcessor.cpp
//
//  Created by Jeremi Niedziela on 31/01/2019.
//

#include "LeptonProcessor.hpp"

LeptonProcessor::LeptonProcessor()
{
  arrayNamesFloat = {
    "LepGood_pt",
    "LepGood_eta",
    "LepGood_phi",
    "LepGood_relIso04"
  };
  
  arrayNamesInt = {
    "LepGood_tightId",
    "LepGood_pdgId",
  };
}

LeptonProcessor::~LeptonProcessor()
{
  
}

bool LeptonProcessor::IsPassingCut(const shared_ptr<Lepton> lepton,
                                   const unique_ptr<LeptonCut> &cut)
{
  // check pt
  if(cut->GetPt().IsOutside(lepton->pt)) return false;
  
  // check isolation
  if(cut->GetRelativeIsolation().IsOutside(lepton->relativeIsolation)) return false;
  
  // check tight id requirement
  if(cut->RequiresTightID() && !lepton->tightID) return false;
  
  return true;
}

vector<shared_ptr<Lepton>> LeptonProcessor::GetLeptonsFromTree()
{
  vector<shared_ptr<Lepton>> leptons = vector<shared_ptr<Lepton>>();
  
  for(int iLepton=0;iLepton<nLeptons;iLepton++){
    auto lepton = make_shared<Lepton>();
    
    // float array variables
    lepton->pt = arrayValuesFloat["LepGood_pt"][iLepton];
    lepton->eta = arrayValuesFloat["LepGood_eta"][iLepton];
    lepton->phi = arrayValuesFloat["LepGood_phi"][iLepton];
    lepton->relativeIsolation = arrayValuesFloat["LepGood_relIso04"][iLepton];
    
    // int array variables
    lepton->tightID = arrayValuesInt["LepGood_tightId"][iLepton];
    lepton->pid = arrayValuesInt["LepGood_pdgId"][iLepton];
    
    leptons.push_back(lepton);
  }

  return leptons;
}

void LeptonProcessor::SaveLeptonsToTree(vector<shared_ptr<Lepton>> leptons)
{
  nLeptons = (int)leptons.size();
  
  for(int iLepton=0;iLepton<nLeptons;iLepton++){
    arrayValuesFloat["LepGood_pt"][iLepton]       = leptons[iLepton]->GetPt();
    arrayValuesFloat["LepGood_eta"][iLepton]      = leptons[iLepton]->GetEta();
    arrayValuesFloat["LepGood_phi"][iLepton]      = leptons[iLepton]->GetPhi();
    arrayValuesFloat["LepGood_relIso04"][iLepton] = leptons[iLepton]->GetRelativeIsolation();
    arrayValuesInt["LepGood_tightId"][iLepton]    = leptons[iLepton]->GetTightID();
    arrayValuesInt["LepGood_pdgId"][iLepton]      = leptons[iLepton]->GetPid();
  }
}

void LeptonProcessor::SetupBranchesForReading(TTree *tree)
{
  // single int variables
  tree->SetBranchAddress("nLepGood",    &nLeptons);
  
  for(string name : arrayNamesFloat){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesFloat[name]);
  }
  
  for(string name : arrayNamesInt){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesInt[name]);
  }
}

void LeptonProcessor::SetupBranchesForWriting(TTree *tree)
{
  tree->Branch("nLepGood",    &nLeptons,     "nLepGood/I");
  
  for(string name : arrayNamesFloat){
    tree->Branch(name.c_str(), &arrayValuesFloat[name], Form("%s[nLepGood]/F", name.c_str()));
  }
  
  for(string name : arrayNamesInt){
    tree->Branch(name.c_str(), &arrayValuesInt[name], Form("%s[nLepGood]/I", name.c_str()));
  }
}
