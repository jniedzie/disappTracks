//
//  Event.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "Event.hpp"

#include <TTreeReaderArray.h>
#include <TLorentzVector.h>
//#include <ROOT/TTreeProcessorMT.hxx>

Events::Events()
{
  
}

Events::Events(string fileName, EDataType dataType, int maxNevents, ESignal iSig)
{
  AddEventsFromFile(fileName,dataType,maxNevents, iSig);
}

Events::Events(const Events &e)
{
  for(auto event : e.events){
    events.push_back(event);
  }
}

Events::~Events()
{
  
}

double Events::weightedSize(){
  double sum=0;
  for(Event *ev : events){
    sum += ev->GetWeight();
  }
  return sum;
}


void Events::AddEventsFromFile(std::string fileName, EDataType dataType, int maxNevents, ESignal iSig)
{
  cout<<"Reading events from:"<<fileName<<endl;
  TFile *inFile = TFile::Open(fileName.c_str());
  TTreeReader reader("tree", inFile);

  TTreeReaderValue<int>   _nTracks(reader, "nIsoTrack");
  TTreeReaderValue<int>   _nVert(reader, "nVert");
  TTreeReaderValue<int>   _nJets(reader, "nJet");
  TTreeReaderValue<int>   _nJetsFwd(reader, "nJetFwd");
  TTreeReaderValue<int>   _nJet30(reader, "nJet30");
  TTreeReaderValue<int>   _nJet30a(reader, "nJet30a");
  TTreeReaderValue<int>   _nLepton(reader, "nLepGood");
  TTreeReaderValue<int>   _nTau(reader, "nTauGood");
  TTreeReaderValue<int>   _nGenChargino(reader, "nGenChargino");
  
  TTreeReaderValue<float> _xSec  (reader,(dataType==kBackground || dataType==kSignal) ? "xsec" : "rho");
  TTreeReaderValue<float> _sumWgt(reader,(dataType==kBackground || dataType==kSignal) ? "wgtsum" : "rho");
  TTreeReaderValue<float> _genWgt(reader,(dataType==kBackground || dataType==kSignal) ? "genWeight" : "rho");

  TTreeReaderValue<float> _metSumEt(reader, "met_sumEt");
  TTreeReaderValue<float> _metPt(reader, "met_pt");
  TTreeReaderValue<float> _metMass(reader, "met_mass");
  TTreeReaderValue<float> _metPhi(reader, "met_phi");
  TTreeReaderValue<float> _metEta(reader, "met_eta");
  
  TTreeReaderValue<int>   _metNoMuTrigger(reader, "HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight");
  
  TTreeReaderValue<float> _metNoMuPt(reader, "metNoMu_pt");
  TTreeReaderValue<float> _metNoMuMass(reader, "metNoMu_mass");
  TTreeReaderValue<float> _metNoMuPhi(reader, "metNoMu_phi");
  TTreeReaderValue<float> _metNoMuEta(reader, "metNoMu_eta");

  TTreeReaderArray<float> _trackEta(reader, "IsoTrack_eta");
  TTreeReaderArray<float> _trackPhi(reader, "IsoTrack_phi");
  TTreeReaderArray<float> _trackCaloEmEnergy(reader, "IsoTrack_caloEmEnergy");
  TTreeReaderArray<float> _trackCaloHadEnergy(reader, "IsoTrack_caloHadEnergy");
  TTreeReaderArray<float> _trackDxyErr(reader, "IsoTrack_edxy");
  TTreeReaderArray<float> _trackDxy(reader, "IsoTrack_dxy");
  TTreeReaderArray<float> _trackDzErr(reader, "IsoTrack_edz");
  TTreeReaderArray<float> _trackDz(reader, "IsoTrack_dz");
  TTreeReaderArray<int>   _trackCharge(reader, "IsoTrack_charge");
  TTreeReaderArray<float> _trackMass(reader, "IsoTrack_mass");
  TTreeReaderArray<float> _trackPt(reader, "IsoTrack_pt");
  TTreeReaderArray<int>   _trackPid(reader, "IsoTrack_pdgId");
  
  TTreeReaderArray<int>   _trackTrackerLayers(reader, "IsoTrack_trackerLayers");
  TTreeReaderArray<int>   _trackPixelLayers(reader, "IsoTrack_pixelLayers");
  TTreeReaderArray<int>   _trackTrackerHits(reader, "IsoTrack_trackerHits");
  TTreeReaderArray<int>   _trackPixelHits(reader, "IsoTrack_pixelHits");
  TTreeReaderArray<int>   _trackMissingInnerPixelHits(reader, "IsoTrack_missingInnerPixelHits");
  TTreeReaderArray<int>   _trackMissingOuterPixelHits(reader, "IsoTrack_missingOuterPixelHits");
  TTreeReaderArray<int>   _trackMissingInnerStripHits(reader, "IsoTrack_missingInnerStripHits");
  TTreeReaderArray<int>   _trackMissingOuterStripHits(reader, "IsoTrack_missingOuterStripHits");
  TTreeReaderArray<int>   _trackMissingInnerTrackerHits(reader, "IsoTrack_missingInnerTrackerHits");
  TTreeReaderArray<int>   _trackMissingOuterTrackerHits(reader, "IsoTrack_missingOuterTrackerHits");
  TTreeReaderArray<int>   _trackMissingMiddleTrackerHits(reader, "IsoTrack_missingMiddleTrackerHits");
  
  TTreeReaderArray<float> _leptonPt(reader, "LepGood_pt");
  TTreeReaderArray<float> _leptonPhi(reader, "LepGood_phi");
  TTreeReaderArray<float> _leptonEta(reader, "LepGood_eta");
  TTreeReaderArray<int>   _leptonThightId(reader, "LepGood_tightId");
  TTreeReaderArray<float> _leptonIsolation(reader, "LepGood_relIso04");
  TTreeReaderArray<int>   _leptonPid(reader, "LepGood_pdgId");
  
  TTreeReaderArray<float> _jetPt(reader,  "Jet_pt");
  TTreeReaderArray<float> _jetEta(reader, "Jet_eta");
  TTreeReaderArray<float> _jetPhi(reader, "Jet_phi");
  TTreeReaderArray<float> _jetMass(reader, "Jet_mass");
  TTreeReaderArray<float> _jetChHEF(reader, "Jet_chHEF");
  TTreeReaderArray<float> _jetNeHEF(reader, "Jet_neHEF");
  
  TTreeReaderArray<float> _jetFwdPt(reader,  "JetFwd_pt");
  TTreeReaderArray<float> _jetFwdEta(reader, "JetFwd_eta");
  TTreeReaderArray<float> _jetFwdPhi(reader, "JetFwd_phi");
  TTreeReaderArray<float> _jetFwdMass(reader, "JetFwd_mass");
  TTreeReaderArray<float> _jetFwdChHEF(reader, "JetFwd_chHEF");
  TTreeReaderArray<float> _jetFwdNeHEF(reader, "JetFwd_neHEF");

  TTreeReaderArray<float> *_dedx[nLayers];
  TTreeReaderArray<int> *_subDetId[nLayers];
  TTreeReaderArray<int> *_sizeX[nLayers];
  TTreeReaderArray<int> *_sizeY[nLayers];

  for(int iLayer=0;iLayer<nLayers;iLayer++){
    _dedx[iLayer] =      new TTreeReaderArray<float>(reader,Form("IsoTrack_dedxByLayer%i",iLayer));
    _subDetId[iLayer] =  new TTreeReaderArray<int>(reader,Form("IsoTrack_subDetIdByLayer%i",iLayer));
    _sizeX[iLayer] =     new TTreeReaderArray<int>(reader,Form("IsoTrack_sizeXbyLayer%i",iLayer));
    _sizeY[iLayer] =     new TTreeReaderArray<int>(reader,Form("IsoTrack_sizeYbyLayer%i",iLayer));
  }
  int iter=-1;
  while (reader.Next()){
    iter++;
    if(maxNevents>0 && iter>maxNevents) break;
    
    Event *newEvent = new Event();

    for(int iTrack=0;iTrack<*_nTracks;iTrack++){
      Track *track = new Track();
      track->SetEta(_trackEta[iTrack]);
      track->SetPhi(_trackPhi[iTrack]);
      track->SetCaloEmEnergy(_trackCaloEmEnergy[iTrack]);
      track->SetCaloHadEnergy(_trackCaloHadEnergy[iTrack]);
      track->SetDxy(_trackDxy[iTrack],_trackDxyErr[iTrack]);
      track->SetDz(_trackDz[iTrack],_trackDzErr[iTrack]);
      track->SetCharge(_trackCharge[iTrack]);
      track->SetMass(_trackMass[iTrack]);
      track->SetPt(_trackPt[iTrack]);
      track->SetPid(_trackPid[iTrack]);
      
      track->SetNtrackerLayers(_trackTrackerLayers[iTrack]);
      track->SetNpixelLayers(_trackPixelLayers[iTrack]);
      track->SetNtrackerHits(_trackTrackerHits[iTrack]);
      track->SetNpixelHits(_trackPixelHits[iTrack]);
      track->SetNmissingInnerPixelHits(_trackMissingInnerPixelHits[iTrack]);
      track->SetNmissingOuterPixelHits(_trackMissingOuterPixelHits[iTrack]);
      track->SetNmissingInnerStripHits(_trackMissingInnerStripHits[iTrack]);
      track->SetNmissingOuterStripHits(_trackMissingOuterStripHits[iTrack]);
      track->SetNmissingInnerTrackerHits(_trackMissingInnerTrackerHits[iTrack]);
      track->SetNmissingOuterTrackerHits(_trackMissingOuterTrackerHits[iTrack]);
      track->SetNmissingMiddleTrackerHits(_trackMissingMiddleTrackerHits[iTrack]);
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        track->SetDeDxInLayer(iLayer, (*_dedx[iLayer])[iTrack]);
        track->SetSubDetIdInLayer(iLayer, (*_subDetId[iLayer])[iTrack]);
        track->SetSizeXinLayer(iLayer, (*_sizeX[iLayer])[iTrack]);
        track->SetSizeYinLayer(iLayer, (*_sizeY[iLayer])[iTrack]);
      }
      newEvent->AddTrack(track);
    }
    
    for(int iJet=0;iJet<*_nJets;iJet++){
      Jet *jet = new Jet();
      jet->SetPt(_jetPt[iJet]);
      jet->SetEta(_jetEta[iJet]);
      jet->SetPhi(_jetPhi[iJet]);
      jet->SetMass(_jetMass[iJet]);
      jet->SetChargedHadronEnergyFraction(_jetChHEF[iJet]);
      jet->SetNeutralHadronEnergyFraction(_jetNeHEF[iJet]);
      jet->SetIsForward(false);
      newEvent->AddJet(jet);
    }
    
    for(int iJet=0;iJet<*_nJetsFwd;iJet++){
      Jet *jet = new Jet();
      jet->SetPt(_jetFwdPt[iJet]);
      jet->SetEta(_jetFwdEta[iJet]);
      jet->SetPhi(_jetFwdPhi[iJet]);
      jet->SetMass(_jetFwdMass[iJet]);
      jet->SetChargedHadronEnergyFraction(_jetFwdChHEF[iJet]);
      jet->SetNeutralHadronEnergyFraction(_jetFwdNeHEF[iJet]);
      jet->SetIsForward(true);
      newEvent->AddJet(jet);
    }
    
    for(int iLepton=0;iLepton<*_nLepton;iLepton++){
      Lepton *lepton = new Lepton();
      lepton->SetPt(_leptonPt[iLepton]);
      lepton->SetEta(_leptonEta[iLepton]);
      lepton->SetPhi(_leptonPhi[iLepton]);
      lepton->SetTightID(_leptonThightId[iLepton]);
      lepton->SetIsolation(_leptonIsolation[iLepton]);
      lepton->SetPid(_leptonPid[iLepton]);
      newEvent->AddLepton(lepton);
    }
    
    double lumi = 41.37 * 1000.;
    double weight = lumi * (*_genWgt) / (*_sumWgt);

    if(dataType==kBackground){
      weight *= (*_xSec);
    }
    if(dataType==kSignal){
      if(*_nGenChargino == 1){
        weight *= 0.001 * signalCrossSectionOneTrack[iSig]; // cross section for given signal (stored in fb, here transformed to pb to match background units
      }
      else if(*_nGenChargino == 2){
        weight *= 0.001 * signalCrossSectionTwoTracks[iSig];
      }
      else{
        cout<<"WARNING -- number of generator-level charginos different than 1 or 2"<<endl;
      }
      weight *= 100.0; // scale up to make it visible
    }
    else if(dataType==kData){
      weight = 1;
    }
    
    newEvent->SetWeight(weight);
    
    newEvent->SetNvertices(*_nVert);
    newEvent->SetNjet30(*_nJet30);
    newEvent->SetNjet30a(*_nJet30a);
    newEvent->SetNlepton(*_nLepton);
    newEvent->SetNtau(*_nTau);
    
    newEvent->SetMetSumEt(*_metSumEt);
    newEvent->SetMetPt(*_metPt);
    newEvent->SetMetMass(*_metMass);
    newEvent->SetMetEta(*_metEta);
    newEvent->SetMetPhi(*_metPhi);
    
    newEvent->SetHasNoMuTrigger(*_metNoMuTrigger);
    newEvent->SetMetNoMuPt(*_metNoMuPt);
    newEvent->SetMetNoMuMass(*_metNoMuMass);
    newEvent->SetMetNoMuEta(*_metNoMuEta);
    newEvent->SetMetNoMuPhi(*_metNoMuPhi);
    
    events.push_back(newEvent);
  }
}

Events* Events::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  Events *outputEvents = new Events(*this);
  
  if(trackCut)  outputEvents = outputEvents->ApplyTrackCut(trackCut);
  if(jetCut)    outputEvents = outputEvents->ApplyJetCut(jetCut);
  if(leptonCut) outputEvents = outputEvents->ApplyLeptonCut(leptonCut);
  if(eventCut)  outputEvents = outputEvents->ApplyEventCut(eventCut);
  
  return outputEvents;
}

Events* Events::ApplyEventCut(EventCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    if(events[iEvent]->IsPassingCut(cut)){
      outputEvents->AddEvent(events[iEvent]);
    }
  }
  return outputEvents;
}

Events* Events::ApplyTrackCut(TrackCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyTrackCut(cut));
  }
  return outputEvents;
}

Events* Events::ApplyJetCut(JetCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyJetCut(cut));
  }
  return outputEvents;
}

Events* Events::ApplyLeptonCut(LeptonCut *cut)
{
  Events *outputEvents = new Events();
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyLeptonCut(cut));
  }
  return outputEvents;
}

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

Event::Event()
{
  
}

Event::~Event()
{
  
}

void Event::Print(){
  for(auto t : tracks){ t->Print(); }
  for(auto j : jets){   j->Print(); }
}

Event* Event::CopyThisEventProperties()
{
  Event *outputEvent = new Event();
    
  outputEvent->SetWeight(weight);
  outputEvent->SetNvertices(nVertices);
  outputEvent->SetNjet30(nJet30);
  outputEvent->SetNjet30a(nJet30a);
  outputEvent->SetNlepton(nLepton);
  outputEvent->SetNtau(nTau);
  
  outputEvent->SetMetSumEt(metSumEt);
  outputEvent->SetMetPt(metPt);
  outputEvent->SetMetMass(metMass);
  outputEvent->SetMetPhi(metPhi);
  outputEvent->SetMetEta(metEta);
  
  outputEvent->SetMetNoMuPt(metNoMuPt);
  outputEvent->SetMetNoMuMass(metNoMuMass);
  outputEvent->SetMetNoMuPhi(metNoMuPhi);
  outputEvent->SetMetNoMuEta(metNoMuEta);
  outputEvent->SetHasNoMuTrigger(metNoMuTrigger);
  
  return outputEvent;
}

Event* Event::ApplyTrackCut(TrackCut *cut)
{
  Event *outputEvent = CopyThisEventProperties();
  
  for(auto j : jets){outputEvent->AddJet(j);}
  for(auto l : leptons){outputEvent->AddLepton(l);}
  
  vector<Track*> tracksPassingCut;
  
  for(auto track : tracks){
    if(track->IsPassingCut(cut)){
      outputEvent->AddTrack(track);
    }
  }
  return outputEvent;
}

Event* Event::ApplyJetCut(JetCut *cut)
{
  Event *outputEvent = CopyThisEventProperties();
  
  for(auto t : tracks){outputEvent->AddTrack(t);}
  for(auto l : leptons){outputEvent->AddLepton(l);}

  vector<Track*> jetPassingCuts;
  
  for(auto jet : jets){
    if(jet->IsPassingCut(cut)){
      outputEvent->AddJet(jet);
    }
  }
  return outputEvent;
}

Event* Event::ApplyLeptonCut(LeptonCut *cut)
{
  Event *outputEvent = CopyThisEventProperties();
  
  for(auto t : tracks){outputEvent->AddTrack(t);}
  for(auto j : jets){outputEvent->AddJet(j);}
  
  vector<Track*> leptonsPassingCuts;
  
  for(auto lepton : leptons){
    if(lepton->IsPassingCut(cut)){
      outputEvent->AddLepton(lepton);
    }
  }
  return outputEvent;
}

bool Event::IsPassingCut(EventCut *cut)
{
  // check the trigger
  if(cut->RequiresMetNoMuTrigger() && !metNoMuTrigger)  return false;
  
  // check number of objects
  if(nLepton < cut->GetMinNleptons() || nLepton > cut->GetMaxNleptons())  return false;
  if(nTau > cut->GetMaxNtau()) return false;
  
  vector<Lepton*> muons;
  for(auto l : leptons){
    if(abs(l->GetPid()) == 13) muons.push_back(l);
  }
  
  // check number of muons
  if(muons.size() < cut->GetMinNmuons() || muons.size() > cut->GetMaxNmuons()){
    return false;
  }
  
  // make sure they have an opposite sign
  if(cut->RequiresTwoOppositeMuons()){
    if(muons.size() != 2) return false;
    if(muons[0]->GetPid() != -muons[1]->GetPid()) return false;
  }
  
  // apply tight muon cuts (tightID flag, pt > 20 GeV, isolation < 0.15
  if(cut->RequiresTightMuon()){
    unsigned int leptonCutOptions = LeptonCut::kIsolated | LeptonCut::kTightID | LeptonCut::kPt20GeV;
    LeptonCut *tightMuonCut = new LeptonCut((LeptonCut::ECut)leptonCutOptions);

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
  
  if(metPt < cut->GetMinMetPt())  return false;
  if(metNoMuPt < cut->GetMinMetNoMuPt())  return false;
  
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
  if(GetNcentralJets() < cut->GetMinNjets()) return false;
  if(GetNtracks() < cut->GetMinNtracks() || GetNtracks() > cut->GetMaxNtracks()) return false;
  

  // find the jet with the highest pt
  Jet *highJet = nullptr;
  double highestPt = -1.0;

  for(int iJet=0;iJet<GetNjets();iJet++){
    if(jets[iJet]->GetPt() > highestPt){
      highestPt = jets[iJet]->GetPt();
      highJet = jets[iJet];
    }
  }

  // check properties of the highest pt jet
  if(cut->RequiresHighJet() && !highJet) return false;
  if(highJet->GetPt() < cut->GetHighJetMinPt()) return false;
  if(fabs(highJet->GetEta()) > cut->GetHighJetMaxEta()) return false;
  if(highJet->GetChargedHadronEnergyFraction() < cut->GetHighJetMinChHEF()) return false;
  if(highJet->GetNeutralHadronEnergyFraction() > cut->GetHighJetMaxNeHEF()) return false;


  
  if(cut->GetMinJetMetPhi() > 0.0){
    TLorentzVector metVector, jetVector;
    metVector.SetPtEtaPhiM(metPt, metEta, metPhi, metMass);

    for(auto j : jets){
      jetVector.SetPtEtaPhiM(j->GetPt(), j->GetEta(), j->GetPhi(), j->GetMass());
      if(fabs(metVector.DeltaPhi(jetVector)) < cut->GetMinJetMetPhi()) return false;
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






