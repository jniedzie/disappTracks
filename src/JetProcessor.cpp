//
//  JetProcessor.cpp
//
//  Created by Jeremi Niedziela on 31/01/2019.
//

#include "JetProcessor.hpp"

JetProcessor::JetProcessor()
{
  arrayNamesFloat = {
    "Jet_pt",
    "Jet_eta",
    "Jet_phi",
    "Jet_mass",
    "Jet_chHEF",
    "Jet_neHEF",
    "JetFwd_pt",
    "JetFwd_eta",
    "JetFwd_phi",
    "JetFwd_mass",
    "JetFwd_chHEF",
    "JetFwd_neHEF"
  };
}

JetProcessor::~JetProcessor()
{
  
}

bool JetProcessor::IsPassingCut(const shared_ptr<Jet> jet,
                                const unique_ptr<JetCut> &cut)
{
  // check jet's kinematical properties
  if(cut->GetPt().IsOutside(jet->pt))  return false;
  if(!jet->isForward && cut->GetEta().IsOutside(jet->eta))  return false;
  if( jet->isForward && cut->GetEtaForward().IsOutside(jet->eta)) return false;
  
  // check hadron energy fractions
  if(cut->GetChargedHadronEnergyFraction().IsOutside(jet->chargedHadronEnergyFraction)) return false;
  if(cut->GetNeutralHadronEnergyFraction().IsOutside(jet->neutralHadronEnergyFraction)) return false;
  
  return true;
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

void JetProcessor::SaveJetsToTree(vector<shared_ptr<Jet>> jets)
{
  nJets = (int)jets.size();
  nJetsFwd = 0;
  
  for(int iJet=0;iJet<nJets;iJet++){
    arrayValuesFloat["Jet_pt"][iJet]    = jets[iJet]->GetPt();
    arrayValuesFloat["Jet_eta"][iJet]   = jets[iJet]->GetEta();
    arrayValuesFloat["Jet_phi"][iJet]   = jets[iJet]->GetPhi();
    arrayValuesFloat["Jet_mass"][iJet]  = jets[iJet]->GetMass();
    arrayValuesFloat["Jet_chHEF"][iJet] = jets[iJet]->GetChargedHadronEnergyFraction();
    arrayValuesFloat["Jet_neHEF"][iJet] = jets[iJet]->GetNeutralHadronEnergyFraction();
  }
}

void JetProcessor::SetupBranchesForReading(TTree *tree)
{
  // single int variables
  tree->SetBranchAddress("nJet",    &nJets);
  tree->SetBranchAddress("nJetFwd", &nJetsFwd);
  
  for(string name : arrayNamesFloat){
    tree->SetBranchAddress(name.c_str(), &arrayValuesFloat[name]);
  }
}

void JetProcessor::SetupBranchesForWriting(TTree *tree)
{
  tree->Branch("nJet",    &nJets,     "nJet/I");
  tree->Branch("nJetFwd", &nJetsFwd,  "nJetFwd/I");
  
  for(string name : arrayNamesFloat){
    tree->Branch(name.c_str(), &arrayValuesFloat[name], Form("%s[nIsoTrack]/F", name.c_str()));
  }
}
