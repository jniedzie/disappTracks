//
//  LeptonProcessor.cpp
//
//  Created by Jeremi Niedziela on 31/01/2019.
//

#include "LeptonProcessor.hpp"

LeptonProcessor::LeptonProcessor()
{
  
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


void LeptonProcessor::SetupBranches(TTree *tree)
{
  // single int variables
  tree->SetBranchAddress("nLepGood",    &nLeptons);
  
  // float array variables
  tree->SetBranchAddress("LepGood_pt",        &arrayValuesFloat["LepGood_pt"]);
  tree->SetBranchAddress("LepGood_eta",       &arrayValuesFloat["LepGood_eta"]);
  tree->SetBranchAddress("LepGood_phi",       &arrayValuesFloat["LepGood_phi"]);
  tree->SetBranchAddress("LepGood_relIso04",  &arrayValuesFloat["LepGood_relIso04"]);
  
  // int array variables
  tree->SetBranchAddress("LepGood_tightId",   &arrayValuesInt["LepGood_tightId"]);
  tree->SetBranchAddress("LepGood_pdgId",     &arrayValuesInt["LepGood_pdgId"]);
}
