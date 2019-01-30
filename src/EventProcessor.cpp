//
//  EventProcessor.cpp
//
//  Created by Jeremi Niedziela on 30/01/2019.
//

#include "EventProcessor.hpp"

EventProcessor::EventProcessor() :
trackProcessor(make_unique<TrackProcessor>())
{
  
}

EventProcessor::~EventProcessor()
{
  
}


void EventProcessor::ApplyTrackCut(shared_ptr<Event> event, const unique_ptr<TrackCut> &cut)
{
  auto track = event->tracks.begin();
  
  while(track != event->tracks.end()){
    if(!trackProcessor->IsPassingCut(*track,cut)) track = event->tracks.erase(track);
    else                                          track++;
  }
}

void EventProcessor::ApplyJetCut(shared_ptr<Event> event, const unique_ptr<JetCut> &cut)
{
  auto jet = event->jets.begin();
  
  while(jet != event->jets.end()){
    if(!(*jet)->IsPassingCut(cut))
    jet = event->jets.erase(jet);
    else{
      // check separation with all tracks in the event
      bool overlapsWithTrack = false;
      double minTrackDeltaR = cut->GetTrackDeltaR().GetMin();
      
      if(minTrackDeltaR > 0){
        for(auto track : event->tracks){
          double deltaR_2 = pow(track->GetPhi() - (*jet)->GetPhi(),2)+pow(track->GetEta() - (*jet)->GetEta(),2);
          
          if(deltaR_2 < (minTrackDeltaR*minTrackDeltaR)){
            overlapsWithTrack = true;
            break;
          }
        }
      }
      
      if(overlapsWithTrack){
        jet = event->jets.erase(jet);
      }
      else{
        jet++;
      }
    }
  }
}

void EventProcessor::ApplyLeptonCut(shared_ptr<Event> event, const unique_ptr<LeptonCut> &cut)
{
  auto lepton = event->leptons.begin();
  
  while(lepton != event->leptons.end()){
    if(!(*lepton)->IsPassingCut(cut))
    lepton = event->leptons.erase(lepton);
    else
    lepton++;
  }
}


bool EventProcessor::IsPassingCut(const shared_ptr<Event> event, const unique_ptr<EventCut> &cut)
{
  // check the trigger
  if(cut->RequiresMetNoMuTrigger() && !event->metNoMuTrigger){
    return false;
  }
  // check filters
  if(cut->GetRequiresPassingAllFilters()){
    if(   !event->flag_goodVertices
       || !event->flag_goodVertices
       || !event->flag_badPFmuon
       || !event->flag_HBHEnoise
       || !event->flag_HBHEnoiseIso
       || !event->flag_EcalDeadCell
       || !event->flag_eeBadSc
    /*|| flag_badChargedCandidate*/
       || !event->flag_ecalBadCalib
       || !event->flag_globalTightHalo2016){
      return false;
    }
  }
  
  // check number of objects
  if(cut->GetNleptons().IsOutside(event->nLepton))  return false;
  if(cut->GetNtaus().IsOutside(event->nTau)) return false;
  
  vector<shared_ptr<Lepton>> muons;
  for(auto l : event->leptons){
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
  
  if(cut->GetMetPt().IsOutside(event->metPt))  return false;
  if(cut->GetMetNoMuPt().IsOutside(event->metNoMuPt))  return false;
  
  // Remove jets that are too close to muons (they will be permanently removed from the event)
  if(cut->RequiresMuJetR0p4()){
    vector<TLorentzVector> muonVectors;
    
    for(auto m : muons){
      TLorentzVector mv;
      mv.SetPtEtaPhiM(m->GetPt(),m->GetEta(),m->GetPhi(), 0.1057);
      muonVectors.push_back(mv);
    }
    
    for(int iJet=0;iJet<event->GetNcentralJets();iJet++){
      shared_ptr<Jet> j = event->GetJet(iJet);
      if(!j->IsForward()){
        
        TLorentzVector jetVector;
        jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
        
        for(auto muonVector : muonVectors){
          if(jetVector.DeltaR(muonVector) < 0.4){
            event->jets.erase(event->jets.begin()+iJet);
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
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      shared_ptr<Track> t = event->tracks[iTrack];
      
      TLorentzVector trackVector;
      trackVector.SetPtEtaPhiM(t->GetPt(), t->GetEta(), t->GetPhi(), t->GetMass());
      
      for(auto muonVector : muonVectors){
        if(trackVector.DeltaR(muonVector) < 0.4){
          event->tracks.erase(event->tracks.begin()+iTrack);
          iTrack--;
          break;
        }
      }
    }
  }
  
  // check number of tracks and jets after removing those that are too close to muons
  if(cut->GetNjets().IsOutside(event->GetNcentralJets())) return false;
  if(cut->GetNtracks().IsOutside(event->GetNtracks())) return false;
  
  // find the jet with the highest pt
  shared_ptr<Jet> leadingJet = nullptr;
  double highestPt = -1.0;
  
  for(int iJet=0;iJet<event->GetNjets();iJet++){
    if(event->jets[iJet]->GetPt() > highestPt){
      highestPt = event->jets[iJet]->GetPt();
      leadingJet = event->jets[iJet];
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
    metVector.SetPtEtaPhiM(event->metPt, event->metEta, event->metPhi, event->metMass);
    
    for(auto j : event->jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(cut->GetJetMetDeltaPhi().IsOutside(fabs(metVector.DeltaPhi(jetVector)) )) return false;
    }
  }
  
  if(cut->RequiresMetNoMuJetPhi0p5()){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(event->metNoMuPt, event->metNoMuEta, event->metNoMuPhi, event->metNoMuMass);
    
    for(auto j : event->jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(fabs(metVector.DeltaPhi(jetVector)) < 0.5) return false;
    }
  }
  
  return true;
}
