//  JetProcessor.cpp
//
//  Created by Jeremi Niedziela on 31/01/2019.

#include "JetProcessor.hpp"

JetProcessor jetProcessor = JetProcessor();

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

bool JetProcessor::IsPassingCut(const shared_ptr<Jet> jet, const JetCut &cut)
{
  // check jet's kinematical properties
  if(cut.pt.IsOutside(jet->pt))  return false;
  if(!jet->isForward && cut.eta.IsOutside(jet->eta))  return false;
  if( jet->isForward && cut.etaForward.IsOutside(jet->eta)) return false;
  
  // check hadron energy fractions
  if(cut.ChHEF.IsOutside(jet->chHEF)) return false;
  if(cut.NeHEF.IsOutside(jet->neHEF)) return false;
  
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
    jet->chHEF = arrayValuesFloat["Jet_chHEF"][iJet];
    jet->neHEF = arrayValuesFloat["Jet_neHEF"][iJet];
    jet->isForward = false;
    
    jets.push_back(jet);
  }
  
  for(int iJet=0;iJet<nJetsFwd;iJet++){
    auto jet = make_shared<Jet>();
    
    jet->pt = arrayValuesFloat["JetFwd_pt"][iJet];
    jet->eta = arrayValuesFloat["JetFwd_eta"][iJet];
    jet->phi = arrayValuesFloat["JetFwd_phi"][iJet];
    jet->mass = arrayValuesFloat["JetFwd_mass"][iJet];
    jet->chHEF = arrayValuesFloat["JetFwd_chHEF"][iJet];
    jet->neHEF = arrayValuesFloat["JetFwd_neHEF"][iJet];
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
    arrayValuesFloat["Jet_chHEF"][iJet] = jets[iJet]->GetChHEF();
    arrayValuesFloat["Jet_neHEF"][iJet] = jets[iJet]->GetNeHEF();
  }
}

void JetProcessor::SetupBranchesForReading(TTree *tree)
{
  // single int variables
  tree->SetBranchAddress("nJet",    &nJets);
  tree->SetBranchAddress("nJetFwd", &nJetsFwd);
  
  for(string name : arrayNamesFloat){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesFloat[name]);
  }
}

void JetProcessor::SetupBranchesForWriting(TTree *tree)
{
  tree->Branch("nJet",    &nJets,     "nJet/I");
  tree->Branch("nJetFwd", &nJetsFwd,  "nJetFwd/I");
  
  for(string name : arrayNamesFloat){
    tree->Branch(name.c_str(), &arrayValuesFloat[name], Form("%s[nJet]/F", name.c_str()));
  }
}
