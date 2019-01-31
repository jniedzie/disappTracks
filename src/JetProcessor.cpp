//
//  JetProcessor.cpp
//
//  Created by Jeremi Niedziela on 31/01/2019.
//

#include "JetProcessor.hpp"

JetProcessor::JetProcessor()
{
  
}

JetProcessor::~JetProcessor()
{
  
}

vector<shared_ptr<Jet>> JetProcessor::GetJetsFromTree()
{
  vector<shared_ptr<Jet>> jets = vector<shared_ptr<Jet>>();
  
  for(int iJet=0;iJet<nJets;iJet++){
    auto jet = make_shared<Jet>();

    jet->pt = arrayValuesFloat["Jet_pt"][iJet];
    jet->eta = arrayValuesFloat["Jet_eta"][iJet];
    jet->phi = arrayValuesFloat["Jet_phi"][iJet];
    jet->mass = arrayValuesFloat["Jet_mass"][iJet];
    jet->chargedHadronEnergyFraction = arrayValuesFloat["Jet_chHEF"][iJet];
    jet->neutralHadronEnergyFraction = arrayValuesFloat["Jet_neHEF"][iJet];
    jet->isForward = false;
    
    jets.push_back(jet);
  }
  
  for(int iJet=0;iJet<nJetsFwd;iJet++){
    auto jet = make_shared<Jet>();
    
    jet->pt = arrayValuesFloat["JetFwd_pt"][iJet];
    jet->eta = arrayValuesFloat["JetFwd_eta"][iJet];
    jet->phi = arrayValuesFloat["JetFwd_phi"][iJet];
    jet->mass = arrayValuesFloat["JetFwd_mass"][iJet];
    jet->chargedHadronEnergyFraction = arrayValuesFloat["JetFwd_chHEF"][iJet];
    jet->neutralHadronEnergyFraction = arrayValuesFloat["JetFwd_neHEF"][iJet];
    jet->isForward = true;
    
    jets.push_back(jet);
  }
  
  return jets;
}


void JetProcessor::SetupBranches(TTree *tree)
{
  // single int variables
  tree->SetBranchAddress("nJet",    &nJets);
  tree->SetBranchAddress("nJetFwd", &nJetsFwd);
  
  // float array variables
  tree->SetBranchAddress("Jet_pt",        &arrayValuesFloat["Jet_pt"]);
  tree->SetBranchAddress("Jet_eta",       &arrayValuesFloat["Jet_eta"]);
  tree->SetBranchAddress("Jet_phi",       &arrayValuesFloat["Jet_phi"]);
  tree->SetBranchAddress("Jet_mass",      &arrayValuesFloat["Jet_mass"]);
  tree->SetBranchAddress("Jet_chHEF",     &arrayValuesFloat["Jet_chHEF"]);
  tree->SetBranchAddress("Jet_neHEF",     &arrayValuesFloat["Jet_neHEF"]);
  
  tree->SetBranchAddress("JetFwd_pt",     &arrayValuesFloat["JetFwd_pt"]);
  tree->SetBranchAddress("JetFwd_eta",    &arrayValuesFloat["JetFwd_eta"]);
  tree->SetBranchAddress("JetFwd_phi",    &arrayValuesFloat["JetFwd_phi"]);
  tree->SetBranchAddress("JetFwd_mass",   &arrayValuesFloat["JetFwd_mass"]);
  tree->SetBranchAddress("JetFwd_chHEF",  &arrayValuesFloat["JetFwd_chHEF"]);
  tree->SetBranchAddress("JetFwd_neHEF",  &arrayValuesFloat["JetFwd_neHEF"]);
}
