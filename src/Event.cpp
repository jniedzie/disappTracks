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
  for(auto &event : e.events){
    events.push_back(make_shared<Event>(*event));
  }
}

Events::~Events()
{
  
}

double Events::weightedSize(){
  double sum=0;
  for(auto &ev : events){
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
  TTreeReaderValue<int>   _flag_goodVertices(reader, "Flag_goodVertices");
  TTreeReaderValue<int>   _flag_badPFmuon(reader, "Flag_BadPFMuonFilter");
  TTreeReaderValue<int>   _flag_HBHEnoise(reader, "Flag_HBHENoiseFilter");
  TTreeReaderValue<int>   _flag_HBHEnoiseIso(reader, "Flag_HBHENoiseIsoFilter");
  TTreeReaderValue<int>   _flag_EcalDeadCell(reader, "Flag_EcalDeadCellTriggerPrimitiveFilter");
  TTreeReaderValue<int>   _flag_eeBadSc(reader, "Flag_eeBadScFilter");
  TTreeReaderValue<int>   _flag_badChargedCandidate(reader, "Flag_BadChargedCandidateFilter");
  TTreeReaderValue<int>   _flag_ecalBadCalib(reader, "Flag_ecalBadCalibFilter");
  TTreeReaderValue<int>   _flag_globalTightHalo2016(reader, "Flag_globalTightHalo2016Filter");
  
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
  TTreeReaderArray<float> _trackRelIso03(reader, "IsoTrack_relIso03");
  
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
    
    unique_ptr<Event> newEvent = unique_ptr<Event>(new Event());

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
      track->SetRelativeIsolation(_trackRelIso03[iTrack]);
      
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
      lepton->SetRelativeIsolation(_leptonIsolation[iLepton]);
      lepton->SetPid(_leptonPid[iLepton]);
      newEvent->AddLepton(lepton);
    }
    
    double lumi = 41.37 * 1000.;
    double weight = lumi * (*_genWgt) / (*_sumWgt);

//    static map<string,set<double>> wgts;
//
//    if(*_genWgt != 1.0 && wgts[fileName].find(*_genWgt) == wgts[fileName].end() ){
//      wgts[fileName].insert(*_genWgt);
//      cout<<*_genWgt<<"\t"<<fileName<<endl;
//    }
    
    if(dataType==kBackground){
      weight *= (*_xSec);
    }
    if(dataType==kSignal){
      // it's not clear how to calculate weights for the signal...
      
      // cross section for given signal (stored in fb, here transformed to pb to match background units
      weight *= 0.001 * (signalCrossSectionOneTrack[iSig] + signalCrossSectionTwoTracks[iSig]);
      
//      if(*_nGenChargino == 1){
//        weight *= 0.001 * signalCrossSectionOneTrack[iSig]; // cross section for given signal (stored in fb, here transformed to pb to match background units
//      }
//      else if(*_nGenChargino == 2){
//        weight *= 0.001 * signalCrossSectionTwoTracks[iSig];
//      }
//      else{
//        cout<<"WARNING -- number of generator-level charginos different than 1 or 2"<<endl;
//      }
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
    
    newEvent->SetGoodVerticesFlag(*_flag_goodVertices);
    newEvent->SetBadPFmuonFlag(*_flag_badPFmuon);
    newEvent->SetHBHEnoiseFlag(*_flag_HBHEnoise);
    newEvent->SetHBHEnoiseIsoFlag(*_flag_HBHEnoiseIso);
    newEvent->SetEcalDeadCellFlag(*_flag_EcalDeadCell);
    newEvent->SetEeBadScFlag(*_flag_eeBadSc);
    newEvent->SetBadChargedCandidateFlag(*_flag_badChargedCandidate);
    newEvent->SetEcalBadCalibFlag(*_flag_ecalBadCalib);
    newEvent->SetGlobalTightHalo2016Flag(*_flag_globalTightHalo2016);
    
    newEvent->SetNgenChargino(*_nGenChargino);
    newEvent->SetXsec(*_xSec);
    newEvent->SetWgtSum(*_sumWgt);
    newEvent->SetGenWeight(*_genWgt);
    
    events.push_back(move(newEvent));
  }
}

void Events::SaveToTree(string fileName)
{
  TFile outFile(fileName.c_str(),"RECREATE");
  outFile.cd();
  TTree *tree = new TTree("tree","tree");
  
  const int nTracks = 100;
  const int nJets = 100;
  const int nLeptons = 100;
  const int nJetsFwd = 0;
  
  int nIsoTracks, nVert, nJet, nJetFwd, nJet30, nJet30a, nLepGood, nTauGood, nGenChargino;
  float xsec, wgtsum, genWeight, met_sumEt, met_pt, met_mass, met_phi, met_eta;
  int metNoMuTrigger, flag_goodVertices, flag_badPFmuon, flag_HBHEnoise, flag_HBHEnoiseIso, flag_EcalDeadCell, flag_eeBadSc, flag_badChargedCandidate, flag_ecalBadCalib, flag_globalTightHalo2016;
  float metNoMu_pt, metNoMu_mass, metNoMu_phi, metNoMu_eta;
  
  float IsoTrack_eta[nTracks], IsoTrack_phi[nTracks], IsoTrack_caloEmEnergy[nTracks], IsoTrack_caloHadEnergy[nTracks], IsoTrack_edxy[nTracks], IsoTrack_dxy[nTracks], IsoTrack_edz[nTracks], IsoTrack_dz[nTracks], IsoTrack_mass[nTracks], IsoTrack_pt[nTracks], IsoTrack_relIso03[nTracks];
  
  int IsoTrack_charge[nTracks], IsoTrack_pdgId[nTracks];
  int IsoTrack_trackerLayers[nTracks], IsoTrack_pixelLayers[nTracks], IsoTrack_trackerHits[nTracks], IsoTrack_pixelHits[nTracks], IsoTrack_missingInnerPixelHits[nTracks], IsoTrack_missingOuterPixelHits[nTracks], IsoTrack_missingInnerStripHits[nTracks], IsoTrack_missingOuterStripHits[nTracks], IsoTrack_missingInnerTrackerHits[nTracks], IsoTrack_missingOuterTrackerHits[nTracks], IsoTrack_missingMiddleTrackerHits[nTracks];
  
  float LepGood_pt[nLeptons], LepGood_phi[nLeptons], LepGood_eta[nLeptons], LepGood_tightId[nLeptons], LepGood_relIso04[nLeptons], LepGood_pdgId[nLeptons];
  float Jet_pt[nJets], Jet_eta[nJets], Jet_phi[nJets], Jet_mass[nJets], Jet_chHEF[nJets], Jet_neHEF[nJets];
  
  
  float JetFwd_pt[nJetsFwd], JetFwd_eta[nJetsFwd], JetFwd_phi[nJetsFwd], JetFwd_mass[nJetsFwd], JetFwd_chHEF[nJetsFwd], JetFwd_neHEF[nJetsFwd];
  
  tree->Branch("nIsoTrack", &nIsoTracks, "nIsoTrack/I");
  tree->Branch("nVert", &nVert, "nVert/I");
  tree->Branch("nJet", &nJet, "nJet/I");
  tree->Branch("nJetFwd", &nJetFwd, "nJetFwd/I");
  tree->Branch("nJet30", &nJet30, "nJet30/I");
  tree->Branch("nJet30a", &nJet30a, "nJet30a/I");
  tree->Branch("nLepGood", &nLepGood, "nLepGood/I");
  tree->Branch("nTauGood", &nTauGood, "nTauGood/I");
  tree->Branch("nGenChargino", &nGenChargino, "nGenChargino/I");
  
  tree->Branch("xsec", &xsec, "xsec/F");
  tree->Branch("wgtsum", &wgtsum, "wgtsum/F");
  tree->Branch("genWeight", &genWeight, "genWeight/F");
  
  tree->Branch("met_sumEt", &met_sumEt, "met_sumEt/F");
  tree->Branch("met_pt", &met_pt, "met_pt/F");
  tree->Branch("met_mass", &met_mass, "met_mass/F");
  tree->Branch("met_phi", &met_phi, "met_phi/F");
  tree->Branch("met_eta", &met_eta, "met_eta/F");
  
  tree->Branch("HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &metNoMuTrigger, "HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/I");
  tree->Branch("Flag_goodVertices", &flag_goodVertices, "Flag_goodVertices/I");
  tree->Branch("Flag_BadPFMuonFilter", &flag_badPFmuon, "Flag_BadPFMuonFilter/I");
  tree->Branch("Flag_HBHENoiseFilter", &flag_HBHEnoise, "Flag_HBHENoiseFilter/I");
  tree->Branch("Flag_HBHENoiseIsoFilter", &flag_HBHEnoiseIso, "Flag_HBHENoiseIsoFilter/I");
  tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &flag_EcalDeadCell, "Flag_EcalDeadCellTriggerPrimitiveFilter/I");
  tree->Branch("Flag_eeBadScFilter", &flag_eeBadSc, "Flag_eeBadScFilter/I");
  tree->Branch("Flag_BadChargedCandidateFilter", &flag_badChargedCandidate, "Flag_BadChargedCandidateFilter/I");
  tree->Branch("Flag_ecalBadCalibFilter", &flag_ecalBadCalib, "Flag_ecalBadCalibFilter/I");
  tree->Branch("Flag_globalTightHalo2016Filter", &flag_globalTightHalo2016, "Flag_globalTightHalo2016Filter/I");
  
  tree->Branch("metNoMu_pt", &metNoMu_pt, "metNoMu_pt/F");
  tree->Branch("metNoMu_mass", &metNoMu_mass, "metNoMu_mass/F");
  tree->Branch("metNoMu_phi", &metNoMu_phi, "metNoMu_phi/F");
  tree->Branch("metNoMu_eta", &metNoMu_eta, "metNoMu_eta/F");
  

  tree->Branch("IsoTrack_eta", &IsoTrack_eta, "IsoTrack_eta[nIsoTrack]/F");
  tree->Branch("IsoTrack_phi", &IsoTrack_phi, "IsoTrack_phi[nIsoTrack]/F");
  tree->Branch("IsoTrack_caloEmEnergy", &IsoTrack_caloEmEnergy, "IsoTrack_caloEmEnergy[nIsoTrack]/F");
  tree->Branch("IsoTrack_caloHadEnergy", &IsoTrack_caloHadEnergy, "IsoTrack_caloHadEnergy[nIsoTrack]/F");
  tree->Branch("IsoTrack_edxy", &IsoTrack_edxy, "IsoTrack_edxy[nIsoTrack]/F");
  tree->Branch("IsoTrack_dxy", &IsoTrack_dxy, "IsoTrack_dxy[nIsoTrack]/F");
  tree->Branch("IsoTrack_edz", &IsoTrack_edz, "IsoTrack_edz[nIsoTrack]/F");
  tree->Branch("IsoTrack_dz", &IsoTrack_dz, "IsoTrack_dz[nIsoTrack]/F");
  tree->Branch("IsoTrack_charge", &IsoTrack_charge, "IsoTrack_charge[nIsoTrack]/I");
  tree->Branch("IsoTrack_mass", &IsoTrack_mass, "IsoTrack_mass[nIsoTrack]/F");
  tree->Branch("IsoTrack_pt", &IsoTrack_pt, "IsoTrack_pt[nIsoTrack]/F");
  tree->Branch("IsoTrack_pdgId", &IsoTrack_pdgId, "IsoTrack_pdgId[nIsoTrack]/I");
  tree->Branch("IsoTrack_relIso03", &IsoTrack_relIso03, "IsoTrack_relIso03[nIsoTrack]/F");
  
  tree->Branch("IsoTrack_trackerLayers", &IsoTrack_trackerLayers, "IsoTrack_trackerLayers[nIsoTrack]/I");
  tree->Branch("IsoTrack_pixelLayers", &IsoTrack_pixelLayers, "IsoTrack_pixelLayers[nIsoTrack]/I");
  tree->Branch("IsoTrack_trackerHits", &IsoTrack_trackerHits, "IsoTrack_trackerHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_pixelHits", &IsoTrack_pixelHits, "IsoTrack_pixelHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingInnerPixelHits", &IsoTrack_missingInnerPixelHits, "IsoTrack_missingInnerPixelHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingOuterPixelHits", &IsoTrack_missingOuterPixelHits, "IsoTrack_missingOuterPixelHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingInnerStripHits", &IsoTrack_missingInnerStripHits, "IsoTrack_missingInnerStripHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingOuterStripHits", &IsoTrack_missingOuterStripHits, "IsoTrack_missingOuterStripHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingInnerTrackerHits", &IsoTrack_missingInnerTrackerHits, "IsoTrack_missingInnerTrackerHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingOuterTrackerHits", &IsoTrack_missingOuterTrackerHits, "IsoTrack_missingOuterTrackerHits[nIsoTrack]/I");
  tree->Branch("IsoTrack_missingMiddleTrackerHits", &IsoTrack_missingMiddleTrackerHits, "IsoTrack_missingMiddleTrackerHits[nIsoTrack]/I");
  

  
  tree->Branch("LepGood_pt", &LepGood_pt, "LepGood_pt[nLepGood]/F");
  tree->Branch("LepGood_phi", &LepGood_phi, "LepGood_phi[nLepGood]/F");
  tree->Branch("LepGood_eta", &LepGood_eta, "LepGood_eta[nLepGood]/F");
  tree->Branch("LepGood_tightId", &LepGood_tightId, "LepGood_tightId[nLepGood]/I");
  tree->Branch("LepGood_relIso04", &LepGood_relIso04, "LepGood_relIso04[nLepGood]/F");
  tree->Branch("LepGood_pdgId", &LepGood_pdgId, "LepGood_pdgId[nLepGood]/I");
  
  tree->Branch("Jet_pt", &Jet_pt, "Jet_pt[nJet]/F");
  tree->Branch("Jet_eta", &Jet_eta, "Jet_eta[nJet]/F");
  tree->Branch("Jet_phi", &Jet_phi, "Jet_phi[nJet]/F");
  tree->Branch("Jet_mass", &Jet_mass, "Jet_mass[nJet]/F");
  tree->Branch("Jet_chHEF", &Jet_chHEF, "Jet_chHEF[nJet]/F");
  tree->Branch("Jet_neHEF", &Jet_neHEF, "Jet_neHEF[nJet]/F");
  
  tree->Branch("JetFwd_pt", &JetFwd_pt, "JetFwd_pt[nJetFwd]/F");
  tree->Branch("JetFwd_eta", &JetFwd_eta, "JetFwd_eta[nJetFwd]/F");
  tree->Branch("JetFwd_phi", &JetFwd_phi, "JetFwd_phi[nJetFwd]/F");
  tree->Branch("JetFwd_mass", &JetFwd_mass, "JetFwd_mass[nJetFwd]/F");
  tree->Branch("JetFwd_chHEF", &JetFwd_chHEF, "JetFwd_chHEF[nJetFwd]/F");
  tree->Branch("JetFwd_neHEF", &JetFwd_neHEF, "JetFwd_neHEF[nJetFwd]/F");
  
  float dedx[nLayers];
  int subDetId[nLayers], sizeX[nLayers], sizeY[nLayers];
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    
    tree->Branch(Form("IsoTrack_dedxByLayer%i",iLayer), &dedx,
                 Form("IsoTrack_dedxByLayer%i[nIsoTrack]/F",iLayer));
    
    tree->Branch(Form("IsoTrack_subDetIdByLayer%i",iLayer), &subDetId,
                 Form("IsoTrack_subDetIdByLayer%i[nIsoTrack]/I",iLayer));
    
    tree->Branch(Form("IsoTrack_sizeXbyLayer%i",iLayer), &sizeX,
                 Form("IsoTrack_sizeXbyLayer%i[nIsoTrack]/I",iLayer));
    
    tree->Branch(Form("IsoTrack_sizeYbyLayer%i",iLayer), &sizeY,
                 Form("IsoTrack_sizeYbyLayer%i[nIsoTrack]/I",iLayer));
  }
  
  for(auto &event : events){
    nVert = event->GetNvertices();
    nIsoTracks = (int)event->GetNtracks();
    nJet = (int)event->GetNjets();
    nJetFwd = 0;
    nJet30 = event->GetNjet30();
    nJet30a = event->GetNjet30a();
    nLepGood = event->GetNlepton();
    nTauGood = event->GetNtau();
    nGenChargino = event->GetNgenChargino();
    xsec = event->GetXsec();
    wgtsum = event->GetWgtSum();
    genWeight = event->GetGenWeight();
    met_sumEt = event->GetMetSumEt();
    met_pt = event->GetMetPt();
    met_mass = event->GetMetMass();
    met_phi = event->GetMetPhi();
    met_eta = event->GetMetEta();
    
    metNoMu_pt = event->GetMetNoMuPt();
    metNoMu_mass = event->GetMetNoMuMass();
    metNoMu_phi = event->GetMetNoMuPhi();
    metNoMu_eta = event->GetMetNoMuEta();
    
    metNoMuTrigger    = event->HetMetNoMuTrigger();
    flag_goodVertices = event->GetGoodVerticesFlag();
    flag_badPFmuon    = event->GetBadPFmuonFlag();
    flag_HBHEnoise    = event->GetHBHEnoiseFlag();
    flag_HBHEnoiseIso = event->GetHBHEnoiseIsoFlag();
    flag_EcalDeadCell = event->GetEcalDeadCellFlag();
    flag_eeBadSc      = event->GetEeBadScFlag();
    flag_ecalBadCalib = event->GetEcalBadCalibFlag();
    flag_badChargedCandidate = event->GetBadChargedCandidateFlag();
    flag_globalTightHalo2016 = event->GetGlobalTightHalo2016Flag();
    
    for(int iTrack=0;iTrack<nIsoTracks;iTrack++){
      IsoTrack_eta[iTrack] = event->GetTrack(iTrack)->GetEta();
      IsoTrack_phi[iTrack] = event->GetTrack(iTrack)->GetPhi();
      IsoTrack_caloEmEnergy[iTrack] = event->GetTrack(iTrack)->GetCaloEmEnergy();
      IsoTrack_caloHadEnergy[iTrack] = event->GetTrack(iTrack)->GetCaloHadEnergy();
      IsoTrack_edxy[iTrack] = event->GetTrack(iTrack)->GetDxyErr();
      IsoTrack_dxy[iTrack] = event->GetTrack(iTrack)->GetDxy();
      IsoTrack_edz[iTrack] = event->GetTrack(iTrack)->GetDzErr();
      IsoTrack_dz[iTrack] = event->GetTrack(iTrack)->GetDz();
      IsoTrack_mass[iTrack] = event->GetTrack(iTrack)->GetMass();
      IsoTrack_pt[iTrack] = event->GetTrack(iTrack)->GetPt();
      IsoTrack_relIso03[iTrack] = event->GetTrack(iTrack)->GetRelativeIsolation();

      IsoTrack_charge[iTrack] = event->GetTrack(iTrack)->GetCharge();
      IsoTrack_pdgId[iTrack] = event->GetTrack(iTrack)->GetPid();
      
      IsoTrack_trackerLayers[iTrack] = event->GetTrack(iTrack)->GetNtrackerLayers();
      IsoTrack_pixelLayers[iTrack] = event->GetTrack(iTrack)->GetNpixelLayers();
      IsoTrack_trackerHits[iTrack] = event->GetTrack(iTrack)->GetNtrackerHits();
      IsoTrack_pixelHits[iTrack] = event->GetTrack(iTrack)->GetNpixelHits();
      IsoTrack_missingInnerPixelHits[iTrack] = event->GetTrack(iTrack)->GetNmissingInnerPixelHits();
      IsoTrack_missingOuterPixelHits[iTrack] = event->GetTrack(iTrack)->GetNmissingOuterPixelHits();
      IsoTrack_missingInnerStripHits[iTrack] = event->GetTrack(iTrack)->GetNmissingInnerStripHits();
      IsoTrack_missingOuterStripHits[iTrack] = event->GetTrack(iTrack)->GetNmissingOuterStripHits();
      IsoTrack_missingInnerTrackerHits[iTrack] = event->GetTrack(iTrack)->GetNmissingInnerTrackerHits();
      IsoTrack_missingOuterTrackerHits[iTrack] = event->GetTrack(iTrack)->GetNmissingOuterTrackerHits();
      IsoTrack_missingMiddleTrackerHits[iTrack] = event->GetTrack(iTrack)->GetNmissingMiddleTrackerHits();
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        dedx[iLayer] = event->GetTrack(iTrack)->GetDeDxInLayer(iLayer);
        subDetId[iLayer] = event->GetTrack(iTrack)->GetSubDetIdInLayer(iLayer);
        sizeX[iLayer] = event->GetTrack(iTrack)->GetSizeXinLayer(iLayer);
        sizeY[iLayer] = event->GetTrack(iTrack)->GetSizeYinLayer(iLayer);
      }
    }
    
    for(int iLep=0;iLep<nLepGood;iLep++){
      LepGood_pt[iLep] = event->GetLepton(iLep)->GetPt();
      
      LepGood_phi[iLep] = event->GetLepton(iLep)->GetPt();
      LepGood_eta[iLep] = event->GetLepton(iLep)->GetEta();
      LepGood_tightId[iLep] = event->GetLepton(iLep)->GetTightID();
      LepGood_relIso04[iLep] = event->GetLepton(iLep)->GetRelativeIsolation();
      LepGood_pdgId[iLep] = event->GetLepton(iLep)->GetPid();
    }
    
    for(int iJet=0;iJet<nJet;iJet++){
      Jet_pt[iJet] = event->GetJet(iJet)->GetPt();
      Jet_eta[iJet] = event->GetJet(iJet)->GetEta();
      Jet_phi[iJet] = event->GetJet(iJet)->GetPhi();
      Jet_mass[iJet] = event->GetJet(iJet)->GetMass();
      Jet_chHEF[iJet] = event->GetJet(iJet)->GetChargedHadronEnergyFraction();
      Jet_neHEF[iJet] = event->GetJet(iJet)->GetNeutralHadronEnergyFraction();
    }
    
    tree->Fill();
  }
  tree->Write();
//  outFile.Close();
}

void Events::LoadEventsFromFiles(vector<shared_ptr<Events>> &eventsSignal,
                                 vector<shared_ptr<Events>> &eventsBackground,
                                 vector<shared_ptr<Events>> &eventsData,
                                 string prefix)
{
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]){
      eventsData.push_back(nullptr);
    }
    else{
      eventsData.push_back(make_shared<Events>((inFileNameData[iData]+prefix+"tree.root"), Events::kData, maxNeventsData));
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]){
      eventsSignal.push_back(nullptr);
    }
    else{
      eventsSignal.push_back(make_shared<Events>((inFileNameSignal[iSig]+prefix+"tree.root"), Events::kSignal, maxNeventsSignal,(ESignal)iSig));
    }
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]){
      eventsBackground.push_back(nullptr);
    }
    else{
      eventsBackground.push_back(make_shared<Events>());
      
      if(prefix==""){
        for(string path : inFileNameBackground[iBck]){
          eventsBackground[iBck]->AddEventsFromFile((path+prefix+"tree.root"),
                                                    Events::kBackground, maxNeventsBackground);
        }
      }
      else{
        string path = inFileNameBackground[iBck][0];
        eventsBackground[iBck]->AddEventsFromFile((path+prefix+"tree.root"),
                                                  Events::kBackground, maxNeventsBackground);
      }
    }
  }
}

void Events::SaveEventsToFiles(vector<shared_ptr<Events>> &eventsSignal,
                               vector<shared_ptr<Events>> &eventsBackground,
                               vector<shared_ptr<Events>> &eventsData,
                               string prefix)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    system(("mkdir -p "+inFileNameSignal[iSig]+prefix).c_str());
    eventsSignal[iSig]->SaveToTree((inFileNameSignal[iSig]+prefix+"tree.root").c_str());
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    
    // merged events will be stored in the first directory for given background
    string path = inFileNameBackground[iBck][0];
    system(("mkdir -p "+path+prefix).c_str());
    eventsBackground[iBck]->SaveToTree((path+prefix+"tree.root").c_str());
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    system(("mkdir -p "+inFileNameData[iData]+prefix).c_str());
    eventsData[iData]->SaveToTree((inFileNameData[iData]+prefix+"tree.root").c_str());
  }
}

void Events::ApplyCuts(vector<shared_ptr<Events>> &eventsSignal,
               vector<shared_ptr<Events>> &eventsBackground,
               vector<shared_ptr<Events>> &eventsData,
               EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    eventsSignal[iSig] = eventsSignal[iSig]->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    eventsBackground[iBck] = eventsBackground[iBck]->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    eventsData[iData] = eventsData[iData]->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  }
}

void Events::PrintYields(vector<shared_ptr<Events>> &eventsSignal,
                         vector<shared_ptr<Events>> &eventsBackground,
                         vector<shared_ptr<Events>> &eventsData)
{
  double nBackgroundTotal=0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck] || !eventsBackground[iBck]) continue;
    nBackgroundTotal += eventsBackground[iBck]->weightedSize();
    
    if(printYields){
      cout<<backgroundTitle[iBck]<<"\t";
      cout<<eventsBackground[iBck]->weightedSize();
      cout<<"\t("<<eventsBackground[iBck]->size()<<")"<<endl;
    }
  }
  
  if(printYields){
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      cout<<signalTitle[iSig]<<"\tN events:\t";
      cout<<eventsSignal[iSig]->weightedSize();
      cout<<"\t("<<eventsSignal[iSig]->size()<<")"<<endl;
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    cout<<signalTitle[iSig]<<"\tS/sqrt(B):\t";
    cout<<eventsSignal[iSig]->weightedSize()/sqrt(nBackgroundTotal+eventsSignal[iSig]->weightedSize())<<endl;
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    cout<<dataTitle[iData]<<"\tS/sqrt(B):\t";
    cout<<eventsData[iData]->weightedSize()/sqrt(nBackgroundTotal+eventsData[iData]->weightedSize())<<endl;
  }
}

shared_ptr<Events> Events::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
//  shared_ptr<Events> outputEvents = shared_ptr<Events>(new Events(*this));
  shared_ptr<Events> outputEvents = make_shared<Events>(*this);
  
  if(trackCut)  outputEvents = outputEvents->ApplyTrackCut(trackCut);
  if(jetCut)    outputEvents = outputEvents->ApplyJetCut(jetCut);
  if(leptonCut) outputEvents = outputEvents->ApplyLeptonCut(leptonCut);
  if(eventCut)  outputEvents = outputEvents->ApplyEventCut(eventCut);
  
  return outputEvents;
}

shared_ptr<Events> Events::ApplyEventCut(EventCut *cut)
{
  shared_ptr<Events> outputEvents = shared_ptr<Events>(new Events());
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    if(events[iEvent]->IsPassingCut(cut)){
      outputEvents->AddEvent(events[iEvent]);
    }
  }
  return outputEvents;
}

shared_ptr<Events> Events::ApplyTrackCut(TrackCut *cut)
{
  shared_ptr<Events> outputEvents = shared_ptr<Events>(new Events());
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyTrackCut(cut));
  }
  return outputEvents;
}

shared_ptr<Events> Events::ApplyJetCut(JetCut *cut)
{
  shared_ptr<Events> outputEvents = shared_ptr<Events>(new Events());
  
  for(int iEvent=0;iEvent<events.size();iEvent++){
    outputEvents->AddEvent(events[iEvent]->ApplyJetCut(cut));
  }
  return outputEvents;
}

shared_ptr<Events> Events::ApplyLeptonCut(LeptonCut *cut)
{
  shared_ptr<Events> outputEvents = shared_ptr<Events>(new Events());
  
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

unique_ptr<Event> Event::CopyThisEventProperties()
{
  unique_ptr<Event> outputEvent = unique_ptr<Event>(new Event());
  
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
  
  outputEvent->SetGoodVerticesFlag(flag_goodVertices);
  outputEvent->SetBadPFmuonFlag(flag_badPFmuon);
  outputEvent->SetHBHEnoiseFlag(flag_HBHEnoise);
  outputEvent->SetHBHEnoiseIsoFlag(flag_HBHEnoiseIso);
  outputEvent->SetEcalDeadCellFlag(flag_EcalDeadCell);
  outputEvent->SetEeBadScFlag(flag_eeBadSc);
  outputEvent->SetBadChargedCandidateFlag(flag_badChargedCandidate);
  outputEvent->SetEcalBadCalibFlag(flag_ecalBadCalib);
  outputEvent->SetGlobalTightHalo2016Flag(flag_globalTightHalo2016);
  
  outputEvent->SetNgenChargino(nGenChargino);
  outputEvent->SetXsec(xsec);
  outputEvent->SetWgtSum(wgtsum);
  outputEvent->SetGenWeight(genWeight);
  
  return outputEvent;
}

unique_ptr<Event> Event::ApplyTrackCut(TrackCut *cut)
{
  unique_ptr<Event> outputEvent = CopyThisEventProperties();
  
  for(auto j : jets){outputEvent->AddJet(j);}
  for(auto l : leptons){outputEvent->AddLepton(l);}
  
  for(auto track : tracks){
    if(track->IsPassingCut(cut)){
      outputEvent->AddTrack(track);
    }
  }
  return outputEvent;
}

unique_ptr<Event> Event::ApplyJetCut(JetCut *cut)
{
  unique_ptr<Event> outputEvent = CopyThisEventProperties();
  
  for(auto t : tracks){outputEvent->AddTrack(t);}
  for(auto l : leptons){outputEvent->AddLepton(l);}
  
  for(auto jet : jets){
    if(jet->IsPassingCut(cut)){
      
      // check separation with all tracks in the event
      bool overlapsWithTrack = false;
      double minTrackDeltaR = cut->GetTrackDeltaR().GetMin();
      
      if(minTrackDeltaR > 0){
        for(auto track : tracks){
          double deltaR_2 = pow(track->GetPhi() - jet->GetPhi(),2)+pow(track->GetEta() - jet->GetEta(),2);
          
          if(deltaR_2 < (minTrackDeltaR*minTrackDeltaR)){
            overlapsWithTrack = true;
            break;
          }
        }
      }
      
      if(!overlapsWithTrack)  outputEvent->AddJet(jet);
    }
  }
  return outputEvent;
}

unique_ptr<Event> Event::ApplyLeptonCut(LeptonCut *cut)
{
  unique_ptr<Event> outputEvent = CopyThisEventProperties();
  
  for(auto t : tracks){outputEvent->AddTrack(t);}
  for(auto j : jets){outputEvent->AddJet(j);}
  
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






