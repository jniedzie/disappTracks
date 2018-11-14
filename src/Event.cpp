//
//  Event.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "Event.hpp"

#include <TLorentzVector.h>

Event::Event()
{
  
}

Event::Event(const Event &e)
{
  for(auto t : e.tracks){
    tracks.push_back(t);
  }
  
  for(auto j : e.jets){
    jets.push_back(j);
  }
  
  for(auto l : e.leptons){
    leptons.push_back(l);
  }
  
  SetWeight(e.weight);
  SetNvertices(e.nVertices);
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

void Event::ApplyTrackCut(const unique_ptr<TrackCut> &cut)
{
  auto track = tracks.begin();

  while(track != tracks.end()){
    if(!(*track)->IsPassingCut(cut))
      track = tracks.erase(track);
    else
      track++;
  }
}

void Event::ApplyJetCut(const unique_ptr<JetCut> &cut)
{
  auto jet = jets.begin();
  
  while(jet != jets.end()){
    if(!(*jet)->IsPassingCut(cut))
      jet = jets.erase(jet);
    else{
      // check separation with all tracks in the event
      bool overlapsWithTrack = false;
      double minTrackDeltaR = cut->GetTrackDeltaR().GetMin();
      
      if(minTrackDeltaR > 0){
        for(auto track : tracks){
          double deltaR_2 = pow(track->GetPhi() - (*jet)->GetPhi(),2)+pow(track->GetEta() - (*jet)->GetEta(),2);
          
          if(deltaR_2 < (minTrackDeltaR*minTrackDeltaR)){
            overlapsWithTrack = true;
            break;
          }
        }
      }
      
      if(overlapsWithTrack){
        jet = jets.erase(jet);
      }
      else{
        jet++;
      }
    }
  }
}

void Event::ApplyLeptonCut(const unique_ptr<LeptonCut> &cut)
{
  auto lepton = leptons.begin();
  
  while(lepton != leptons.end()){
    if(!(*lepton)->IsPassingCut(cut))
      lepton = leptons.erase(lepton);
    else
      lepton++;
  }
}


bool Event::IsPassingCut(const unique_ptr<EventCut> &cut)
{
  // check the trigger
  if(cut->RequiresMetNoMuTrigger() && !metNoMuTrigger){
    return false;
  }
  // check filters
  if(cut->GetRequiresPassingAllFilters()){
    if(   !flag_goodVertices  || !flag_goodVertices         || !flag_badPFmuon
       || !flag_HBHEnoise     || !flag_HBHEnoiseIso         || !flag_EcalDeadCell
       || !flag_eeBadSc       /*|| flag_badChargedCandidate*/   || !flag_ecalBadCalib
       || !flag_globalTightHalo2016){
      return false;
    }
  }
  
  // check number of objects
  if(cut->GetNleptons().IsOutside(nLepton))  return false;
  if(cut->GetNtaus().IsOutside(nTau)) return false;
  
  vector<Lepton*> muons;
  for(auto l : leptons){
    if(abs(l->GetPid()) == 13) muons.push_back(l);
  }
  
  // check number of muons
  if(cut->GetNmuons().IsOutside((int)muons.size())) return false;
  
  // make sure they have an opposite sign
  if(cut->RequiresTwoOppositeMuons()){
    if(muons.size() != 2) return false;
    if(muons[0]->GetPid() != -muons[1]->GetPid()) return false;
  }
  
  // apply tight muon cuts (tightID flag, pt > 20 GeV, isolation < 0.15
  if(cut->RequiresTightMuon()){
    unique_ptr<LeptonCut> tightMuonCut = unique_ptr<LeptonCut>(new LeptonCut());
    tightMuonCut->SetRelativeIsolation(range<double>(-inf,0.15));
    tightMuonCut->SetRequireTightID(true);
    tightMuonCut->SetPt(range<double>(20.0,inf));
    
    bool atLeastOneTightMuon = false;
    for(auto muon : muons){
      if(muon->IsPassingCut(tightMuonCut)) atLeastOneTightMuon = true;
    }
    if(!atLeastOneTightMuon) return false;
  }
  
  // check that invariant mass of muons is close to Z mass
  if(cut->RequiresMuonsFromZ()){
    TLorentzVector muon1vector, muon2vector, muonVectorSum;

    if(muons.size() != 2){
      cout<<"ERROR -- requested muons to come from Z decay, but there is "<<muons.size()<<" muons in the event!!"<<endl;
      cout<<"This event will be discarded!! Maybe you should require exactly two muons in the event?"<<endl;
      return false;
    }
    
    muon1vector.SetPtEtaPhiM(muons[0]->GetPt(),muons[0]->GetEta(),muons[0]->GetPhi(), 0.1057);
    muon2vector.SetPtEtaPhiM(muons[1]->GetPt(),muons[1]->GetEta(),muons[1]->GetPhi(), 0.1057);
    
    muonVectorSum += muon1vector;
    muonVectorSum += muon2vector;

    if(muonVectorSum.M() < 60. || muonVectorSum.M() > 120.) return false;
  }
  
  if(cut->GetMetPt().IsOutside(metPt))  return false;
  if(cut->GetMetNoMuPt().IsOutside(metNoMuPt))  return false;
  
  // Remove jets that are too close to muons (they will be permanently removed from the event)
  if(cut->RequiresMuJetR0p4()){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), 0.1057);
      muonVectors.push_back(mv);
    }
    
    for(int iJet=0;iJet<GetNcentralJets();iJet++){
      Jet *j = GetJet(iJet);
      if(!j->IsForward()){
        
        TLorentzVector jetVector;
        jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
        
        for(auto muonVector : muonVectors){
          if(jetVector.DeltaR(muonVector) < 0.4){
            jets.erase(jets.begin()+iJet);
            iJet--;
            break;
          }
        }
      }
    }
  }
  
  // Remove tracks that are too close to muons (they will be permanently removed from the event)
  if(cut->RequiresMuTrackR0p4()){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), 0.1057);
      muonVectors.push_back(mv);
    }
    
    for(int iTrack=0;iTrack<GetNtracks();iTrack++){
      Track *t = tracks[iTrack];
      
      TLorentzVector trackVector;
      trackVector.SetPtEtaPhiM(t->GetPt(), t->GetEta(), t->GetPhi(), t->GetMass());
      
      for(auto muonVector : muonVectors){
        if(trackVector.DeltaR(muonVector) < 0.4){
          tracks.erase(tracks.begin()+iTrack);
          iTrack--;
          break;
        }
      }
    }
  }
  
  // check number of tracks and jets after removing those that are too close to muons
  if(cut->GetNjets().IsOutside(GetNcentralJets())) return false;
  if(cut->GetNtracks().IsOutside(GetNtracks())) return false;

  // find the jet with the highest pt
  Jet *leadingJet = nullptr;
  double highestPt = -1.0;

  for(int iJet=0;iJet<GetNjets();iJet++){
    if(jets[iJet]->GetPt() > highestPt){
      highestPt = jets[iJet]->GetPt();
      leadingJet = jets[iJet];
    }
  }

  // check properties of the highest pt jet
  if(cut->RequiresHighJet() && !leadingJet) return false;
  if(cut->GetLeadingJetPt().IsOutside(leadingJet->GetPt())) return false;
  if(cut->GetLeadingJetEta().IsOutside(leadingJet->GetEta())) return false;
  if(cut->GetLeadingJetChHEF().IsOutside(leadingJet->GetChargedHadronEnergyFraction())) return false;
  if(cut->GetLeadingJetNeHEF().IsOutside(leadingJet->GetNeutralHadronEnergyFraction())) return false;

  if(cut->GetJetMetDeltaPhi().GetMin() > 0.0){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(metPt, metEta, metPhi, metMass);

    for(auto j : jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(cut->GetJetMetDeltaPhi().IsOutside(fabs(metVector.DeltaPhi(jetVector)) )) return false;
    }
  }

  if(cut->RequiresMetNoMuJetPhi0p5()){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(metNoMuPt, metNoMuEta, metNoMuPhi, metNoMuMass);

    for(auto j : jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(fabs(metVector.DeltaPhi(jetVector)) < 0.5) return false;
    }
  }
  
  return true;
}






