//
//  EventSet.cpp
//
//  Created by Jeremi Niedziela on 08/11/2018.
//

#include "EventSet.hpp"
#include "HistSet.hpp"

#include <TTreeReaderArray.h>
#include <TLorentzVector.h>

EventSet::EventSet()
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    eventsSignal.push_back(vector<shared_ptr<Event>>());
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    eventsBackground.push_back(vector<shared_ptr<Event>>());
  }
  for(int iData=0;iData<kNdata;iData++){
    eventsData.push_back(vector<shared_ptr<Event>>());
  }
}

EventSet::EventSet(string fileName, EDataType dataType, int maxNevents, ESignal iSig)
{
  AddEventsFromFile(fileName,dataType,maxNevents,iSig);
}

EventSet::EventSet(const EventSet &e)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    eventsSignal.push_back(vector<shared_ptr<Event>>());
    for(auto &event : e.eventsSignal[iSig]){
      eventsSignal[iSig].push_back(make_shared<Event>(*event));
    }
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    eventsBackground.push_back(vector<shared_ptr<Event>>());
    for(auto &event : e.eventsBackground[iBck]){
      eventsBackground[iBck].push_back(make_shared<Event>(*event));
    }
  }
  for(int iData=0;iData<kNdata;iData++){
    eventsData.push_back(vector<shared_ptr<Event>>());
    for(auto &event : e.eventsData[iData]){
      eventsData[iData].push_back(make_shared<Event>(*event));
    }
  }
}

EventSet::~EventSet()
{
  
}

void EventSet::SaveToTree(string fileName, EDataType dataType, int setIter)
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
  
  
  function<void(shared_ptr<Event>, TTree*)> func = [&](shared_ptr<Event> event, TTree *tree) -> void {
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
  };
  
  if(dataType == kSignal){
    for(auto &event : eventsSignal[(ESignal)setIter]){func(event, tree);}
  }
  else if(dataType == kBackground){
    for(auto &event : eventsBackground[(EBackground)setIter]){func(event, tree);}
  }
  else if(dataType == kData){
    for(auto &event : eventsData[(EData)setIter]){func(event, tree);}
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  
  tree->Write();
  //  outFile.Close();
}

void EventSet::LoadEventsFromFiles(string prefix)
{
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    
    if(prefix==""){
      for(string path : inFileNameBackground[iBck]){
        AddEventsFromFile((path+prefix+"tree.root"),kBackground, maxNeventsBackground, iBck);
      }
    }
    else{
      string path = inFileNameBackground[iBck][0];
      AddEventsFromFile((path+prefix+"tree.root"),kBackground, maxNeventsBackground, iBck);
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    AddEventsFromFile((inFileNameSignal[iSig]+prefix+"tree.root"),kSignal,maxNeventsSignal,iSig);
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    AddEventsFromFile((inFileNameData[iData]+prefix+"tree.root"),kData,maxNeventsData,iData);
  }
}

void EventSet::LoadEventsFromFiles(EDataType dataType, int setIter, string prefix)
{
  if(dataType == kBackground){
    if(prefix==""){
      for(string path : inFileNameBackground[setIter]){
        AddEventsFromFile((path+prefix+"tree.root"),kBackground, maxNeventsBackground, setIter);
      }
    }
    else{
      string path = inFileNameBackground[setIter][0];
      AddEventsFromFile((path+prefix+"tree.root"),kBackground, maxNeventsBackground, setIter);
    }
  }
  else if(dataType == kSignal){
    AddEventsFromFile((inFileNameSignal[setIter]+prefix+"tree.root"),kSignal,maxNeventsSignal,setIter);
  }
  else if(dataType == kData){
    AddEventsFromFile((inFileNameData[setIter]+prefix+"tree.root"),kData,maxNeventsData,setIter);
  }
}

void EventSet::SaveEventsToFiles(string prefix)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    system(("mkdir -p "+inFileNameSignal[iSig]+prefix).c_str());
    
    SaveToTree((inFileNameSignal[iSig]+prefix+"tree.root").c_str(), kSignal, (int)iSig);
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    
    // merged events will be stored in the first directory for given background
    string path = inFileNameBackground[iBck][0];
    system(("mkdir -p "+path+prefix).c_str());
    SaveToTree((path+prefix+"tree.root").c_str(), kBackground, (int)iBck);
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    system(("mkdir -p "+inFileNameData[iData]+prefix).c_str());
    SaveToTree((inFileNameData[iData]+prefix+"tree.root").c_str(), kData, (int)iData);
  }
}

void EventSet::PrintYields()
{
  double nBackgroundTotal=0;
  int nBackgroundTotalRaw=0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    nBackgroundTotal += weightedSize(kBackground,iBck);
    nBackgroundTotalRaw += size(kBackground,iBck);
    
    if(printYields){
      cout<<backgroundTitle[iBck]<<"\t";
      cout<<weightedSize(kBackground, (int)iBck);
      cout<<"\t("<<size(kBackground,(int)iBck)<<")"<<endl;
    }
  }
  
  if(printYields){
    cout<<"Background total:\t"<<nBackgroundTotal<<"\t("<<nBackgroundTotalRaw<<")"<<endl;
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      cout<<signalTitle[iSig]<<"\tN events:\t";
      cout<<weightedSize(kSignal,iSig);
      cout<<"\t("<<size(kSignal,iSig)<<")"<<endl;
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    cout<<signalTitle[iSig]<<"\tS/sqrt(B):\t";
    cout<<weightedSize(kSignal,iSig)/sqrt(nBackgroundTotal+weightedSize(kSignal,iSig))<<endl;
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    cout<<dataTitle[iData]<<"\tS/sqrt(B):\t";
    cout<<weightedSize(kData,iData)/sqrt(nBackgroundTotal+weightedSize(kData,iData))<<endl;
  }
}

void EventSet::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    eventsSignal[iSig] = ApplyCuts(eventCut, trackCut, jetCut, leptonCut, kSignal, iSig);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    eventsBackground[iBck] = ApplyCuts(eventCut, trackCut, jetCut, leptonCut, kBackground, iBck);
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    eventsData[iData] = ApplyCuts(eventCut, trackCut, jetCut, leptonCut, kData, iData);
  }
}

void EventSet::ApplyCutsInPlace(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    for(int iEvent=0;iEvent<eventsSignal[iSig].size();){
      
      eventsSignal[iSig][iEvent]->ApplyTrackCutInPlace(trackCut);
      eventsSignal[iSig][iEvent]->ApplyJetCutInPlace(jetCut);
      eventsSignal[iSig][iEvent]->ApplyLeptonCutInPlace(leptonCut);
      
      if(!eventsSignal[iSig][iEvent]->IsPassingCut(eventCut)){
        eventsSignal[iSig].erase(eventsSignal[iSig].begin()+iEvent);
      }
      else{
        iEvent++;
      }
      
    }
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    for(int iEvent=0;iEvent<eventsBackground[iBck].size();){
      
      eventsBackground[iBck][iEvent]->ApplyTrackCutInPlace(trackCut);
      eventsBackground[iBck][iEvent]->ApplyJetCutInPlace(jetCut);
      eventsBackground[iBck][iEvent]->ApplyLeptonCutInPlace(leptonCut);
      
      if(!eventsBackground[iBck][iEvent]->IsPassingCut(eventCut)){
        eventsBackground[iBck].erase(eventsBackground[iBck].begin()+iEvent);
      }
      else{
        iEvent++;
      }
      
    }
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    for(int iEvent=0;iEvent<eventsData[iData].size();){
      
      eventsData[iData][iEvent]->ApplyTrackCutInPlace(trackCut);
      eventsData[iData][iEvent]->ApplyJetCutInPlace(jetCut);
      eventsData[iData][iEvent]->ApplyLeptonCutInPlace(leptonCut);
      
      if(!eventsData[iData][iEvent]->IsPassingCut(eventCut)){
        eventsData[iData].erase(eventsData[iData].begin()+iEvent);
      }
      else{
        iEvent++;
      }
      
    }
  }
}

void EventSet::DrawStandardPlots(string prefix)
{
  // Create standard per event, per track and per jet plots
  map<string, HistSet*> hists;
  
  hists["nVertices"]  = new HistSet(kNvertices);
  hists["nIsoTrack"]  = new HistSet(kNisoTracks);
  hists["nJet"]       = new HistSet(kNjets);
  hists["nJet30"]     = new HistSet(kNjets30);
  hists["nJet30a"]    = new HistSet(kNjets30a);
  hists["nMetSumEt"]  = new HistSet(kMetSumEt);
  hists["nMetPt"]     = new HistSet(kMetPt);
  hists["nMetMass"]   = new HistSet(kMetMass);
  hists["nMetEta"]    = new HistSet(kMetEta);
  hists["nMetPhi"]    = new HistSet(kMetPhi);
  hists["nMetJetDphi"]= new HistSet(kMetJetDphi);
  
  hists["nClustersPerTrack"]  = new HistSet(kTrackNclusters);
  hists["totalDeDx"]          = new HistSet(kTrackTotalDedx);
  hists["totalDeDxByNclusters"] = new HistSet(kTrackDedxPerCluster);
  hists["missingOuterTracker"]  = new HistSet(kTrackMissingOuterTrackerHits);
  
  hists["pt"]           = new HistSet(kTrackPt);
  hists["eta"]          = new HistSet(kTrackEta);
  hists["phi"]          = new HistSet(kTrackPhi);
  hists["caloEm"]       = new HistSet(kTrackCaloEm);
  hists["caloHad"]      = new HistSet(kTrackCaloHad);
  hists["pixelHits"]    = new HistSet(kTrackPixelHits);
  hists["trackerHits"]  = new HistSet(kTrackTrackerHits);
  hists["isolation"]    = new HistSet(kTrackRelativeIsolation);
  hists["absIsolation"] = new HistSet(kTrackAbsoluteIsolation);
  hists["trackMetDphi"] = new HistSet(kTrackMetDphi);
  hists["dedx"]         = new HistSet(kTrackDedxPerHit);
  
  hists["dxy"]    = new HistSet(kTrackDxy);
  hists["dz"]     = new HistSet(kTrackDz);
  hists["charge"] = new HistSet(kTrackCharge);
  hists["mass"]   = new HistSet(kTrackMass);
  hists["pid"]    = new HistSet(kTrackPid);
  
  hists["jet_pt"]     = new HistSet(kJetPt);
  hists["jet_eta"]    = new HistSet(kJetEta);
  hists["jet_phi"]    = new HistSet(kJetPhi);
  hists["jetTrackDr"] = new HistSet(kJetTrackDr);
  hists["jetCHF"] = new HistSet(kJetCHF);
  hists["jetNHF"] = new HistSet(kJetNHF);
  
  for(pair<string, HistSet*> hist : hists){
    hist.second->FillFromEvents(make_shared<EventSet>(*this));
  }
  
  // Plot histograms
  TCanvas *canvasEvents = new TCanvas((prefix+"Events").c_str(),(prefix+"Events").c_str(),2880,1800);
  canvasEvents->Divide(3,2);
  
  hists["nVertices"]->Draw(canvasEvents,1);
  hists["nIsoTrack"]->Draw(canvasEvents,2);
  hists["nJet"]->Draw(canvasEvents,3);
  hists["jet_pt"]->Draw(canvasEvents, 4);
  hists["nMetPt"]->Draw(canvasEvents,5);
  hists["nMetJetDphi"]->Draw(canvasEvents,6);
  
  TCanvas *canvasTrack = new TCanvas((prefix+"Tracks").c_str(),(prefix+"Tracks").c_str(),2880,1800);
  canvasTrack->Divide(4,3);
  
  hists["pt"]->Draw(canvasTrack,1);
  hists["caloEm"]->Draw(canvasTrack,2);
  hists["caloHad"]->Draw(canvasTrack,3);
  hists["missingOuterTracker"]->Draw(canvasTrack,4);
  hists["pixelHits"]->Draw(canvasTrack,5);
  hists["trackerHits"]->Draw(canvasTrack,6);
  hists["isolation"]->Draw(canvasTrack,7);
  hists["trackMetDphi"]->Draw(canvasTrack,8);
  hists["eta"]->Draw(canvasTrack,9);
  hists["dedx"]->Draw(canvasTrack,10);
  hists["absIsolation"]->Draw(canvasTrack,11);
  hists["dz"]->Draw(canvasTrack,12);
  
  TCanvas *canvasJets = new TCanvas("Jets","Jets",2880,1800);
  canvasJets->Divide(3,2);
  
  hists["jetTrackDr"]->Draw(canvasJets, 1);
  hists["jet_eta"]->Draw(canvasJets, 2);
  hists["jet_phi"]->Draw(canvasJets, 3);
  hists["jetCHF"]->Draw(canvasJets, 4);
  hists["jetNHF"]->Draw(canvasJets, 5);
}

void EventSet::DrawPerLayerPlots()
{
  HistSet *dedxPerLayer = new HistSet(kDedx);
  dedxPerLayer->FillFromEvents(make_shared<EventSet>(*this));
  dedxPerLayer->DrawPerLayer();
  
  //  HistSet *sizeXperLayer = new HistSet(kSizeX);
  //  sizeXperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  //  sizeXperLayer->DrawPerLayer();
  //
  //  HistSet *sizeYperLayer = new HistSet(kSizeY);
  //  sizeYperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  //  sizeYperLayer->DrawPerLayer();
}


int EventSet::size(EDataType dataType, int setIter)
{
  if(dataType == kSignal){
    return (int)eventsSignal[(ESignal)setIter].size();
  }
  else if(dataType == kBackground){
    return (int)eventsBackground[(EBackground)setIter].size();
  }
  else if(dataType == kData){
    return (int)eventsData[(EData)setIter].size();
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  return 1;
}

double EventSet::weightedSize(EDataType dataType, int setIter)
{
  double sum=0;
  
  if(dataType == kSignal){
    for(auto &ev : eventsSignal[(ESignal)setIter]){sum += ev->GetWeight();}
  }
  else if(dataType == kBackground){
    for(auto &ev : eventsBackground[(EBackground)setIter]){sum += ev->GetWeight();}
  }
  else if(dataType == kData){
    for(auto &ev : eventsData[(EData)setIter]){sum += ev->GetWeight();}
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  
  return sum;
}

shared_ptr<Event> EventSet::At(EDataType dataType, int setIter, int index)
{
  if(dataType == kSignal){
    return eventsSignal[(ESignal)setIter][index];
  }
  else if(dataType == kBackground){
    return eventsBackground[(EBackground)setIter][index];
  }
  else if(dataType == kData){
    return eventsData[(EData)setIter][index];
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
}

void EventSet::AddEvent(shared_ptr<Event> event, EDataType dataType, int setIter)
{
  
  if(dataType == kSignal){
    eventsSignal[(ESignal)setIter].push_back(event);
  }
  else if(dataType == kBackground){
    eventsBackground[(EBackground)setIter].push_back(event);
  }
  else if(dataType == kData){
    eventsData[(EData)setIter].push_back(event);
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
}

void EventSet::AddEventsFromFile(std::string fileName, EDataType dataType, int maxNevents, int setIter)
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
    
    shared_ptr<Event> newEvent = shared_ptr<Event>(new Event());
    
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
      weight *= 0.001 * (signalCrossSectionOneTrack[(ESignal)setIter] +
                         signalCrossSectionTwoTracks[(ESignal)setIter]);
      
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
    
    if(dataType == kSignal){
      eventsSignal[setIter].push_back(newEvent);
    }
    else if(dataType == kBackground){
      eventsBackground[setIter].push_back(newEvent);
    }
    else if(dataType == kData){
      eventsData[setIter].push_back(newEvent);
    }
    else{
      throw out_of_range("Unknown data type provided");
    }
    
  }
}


vector<shared_ptr<Event>> EventSet::ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut, EDataType dataType, int setIter)
{
  vector<shared_ptr<Event>> outputEvents;
  
  if(dataType == kSignal){
    outputEvents = eventsSignal[(ESignal)setIter];
  }
  else if(dataType == kBackground){
    outputEvents = eventsBackground[(EBackground)setIter];
  }
  else if(dataType == kData){
    outputEvents = eventsData[(EData)setIter];
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  
  if(trackCut)  outputEvents = ApplyTrackCut(outputEvents, trackCut);
  if(jetCut)    outputEvents = ApplyJetCut(outputEvents, jetCut);
  if(leptonCut) outputEvents = ApplyLeptonCut(outputEvents, leptonCut);
  if(eventCut)  outputEvents = ApplyEventCut(outputEvents, eventCut);
  
  return outputEvents;
}

vector<shared_ptr<Event>> EventSet::ApplyEventCut(vector<shared_ptr<Event>> events, EventCut *cut)
{
  vector<shared_ptr<Event>> outputEvents;
  
  for(int iEvent=0;iEvent<(int)events.size();iEvent++){
    if(events[iEvent]->IsPassingCut(cut)){
      outputEvents.push_back(events[iEvent]);
    }
  }
  
  return outputEvents;
}

vector<shared_ptr<Event>> EventSet::ApplyTrackCut(vector<shared_ptr<Event>> events, TrackCut *cut)
{
  vector<shared_ptr<Event>> outputEvents;
  
  for(int iEvent=0;iEvent<(int)events.size();iEvent++){
    outputEvents.push_back(events[iEvent]->ApplyTrackCut(cut));
  }
  return outputEvents;
}

vector<shared_ptr<Event>> EventSet::ApplyJetCut(vector<shared_ptr<Event>> events, JetCut *cut)
{
  vector<shared_ptr<Event>> outputEvents;
  
  for(int iEvent=0;iEvent<(int)events.size();iEvent++){
    outputEvents.push_back(events[iEvent]->ApplyJetCut(cut));
  }
  return outputEvents;
}

vector<shared_ptr<Event>> EventSet::ApplyLeptonCut(vector<shared_ptr<Event>> events, LeptonCut *cut)
{
  vector<shared_ptr<Event>> outputEvents;
  
  for(int iEvent=0;iEvent<(int)events.size();iEvent++){
    outputEvents.push_back(events[iEvent]->ApplyLeptonCut(cut));
  }
  return outputEvents;
}
