//
//  EventSet.cpp
//
//  Created by Jeremi Niedziela on 08/11/2018.
//

#include "EventSet.hpp"
#include "HistSet.hpp"

EventSet::EventSet() :
trackProcessor(make_unique<TrackProcessor>()),
jetProcessor(make_unique<JetProcessor>()),
helixProcessor(make_unique<HelixProcessor>()),
eventProcessor(make_unique<EventProcessor>()),
leptonProcessor(make_unique<LeptonProcessor>())
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

EventSet::EventSet(string fileName, xtracks::EDataType dataType, int maxNevents, ESignal iSig) :
trackProcessor(make_unique<TrackProcessor>()),
eventProcessor(make_unique<EventProcessor>())
{
  AddEventsFromFile(fileName,dataType,maxNevents,iSig);
}

EventSet::EventSet(const EventSet &e)
{
  eventProcessor = make_unique<EventProcessor>();
  trackProcessor = make_unique<TrackProcessor>();
  
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

void EventSet::SaveToTree(string fileName, xtracks::EDataType dataType, int setIter)
{
  TFile outFile(fileName.c_str(),"RECREATE");
  outFile.cd();
  TTree *tree = new TTree("tree","tree");
	
  const int nHelices = 100;
  
  unsigned long long evt;
  uint lumi, run;
  int nVert, nJet30, nJet30a, nTauGood, nGenChargino;
  float vertex_x, vertex_y, vertex_z;
  float xsec, wgtsum, genWeight, met_sumEt, met_pt, met_mass, met_phi, met_eta;
  int metNoMuTrigger, flag_goodVertices, flag_badPFmuon, flag_HBHEnoise, flag_HBHEnoiseIso, flag_EcalDeadCell, flag_eeBadSc, flag_badChargedCandidate, flag_ecalBadCalib, flag_globalTightHalo2016;
  float metNoMu_pt, metNoMu_mass, metNoMu_phi, metNoMu_eta;
	
  int nFittedHelices;
  float helix_x[nHelices],  helix_y[nHelices],  helix_z[nHelices],
        helix_px[nHelices], helix_py[nHelices], helix_pz[nHelices];
  
  int helix_charge[nHelices];
  
	trackProcessor->SetupBranchesForWriting(tree);
  jetProcessor->SetupBranchesForWriting(tree);
  leptonProcessor->SetupBranchesForWriting(tree);
	
  tree->Branch("lumi", &lumi, "lumi/i");
  tree->Branch("run", &run, "run/i");
  tree->Branch("evt", &evt, "evt/l");

  tree->Branch("nFittedHelices", &nFittedHelices, "nFittedHelices/I");
  tree->Branch("nVert", &nVert, "nVert/I");
  tree->Branch("vertex_x", &vertex_x, "vertex_x/F");
  tree->Branch("vertex_y", &vertex_y, "vertex_y/F");
  tree->Branch("vertex_z", &vertex_z, "vertex_z/F");
  tree->Branch("nJet30", &nJet30, "nJet30/I");
  tree->Branch("nJet30a", &nJet30a, "nJet30a/I");
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
	
  tree->Branch("helix_x", &helix_x, "helix_x[nFittedHelices]/F");
  tree->Branch("helix_y", &helix_y, "helix_y[nFittedHelices]/F");
  tree->Branch("helix_z", &helix_z, "helix_z[nFittedHelices]/F");
  tree->Branch("helix_px", &helix_px, "helix_px[nFittedHelices]/F");
  tree->Branch("helix_py", &helix_py, "helix_py[nFittedHelices]/F");
  tree->Branch("helix_pz", &helix_pz, "helix_pz[nFittedHelices]/F");
  
  tree->Branch("helix_charge", &helix_charge, "helix_charge[nFittedHelices]/I");
  
  function<void(shared_ptr<Event>, TTree*)> func = [&](shared_ptr<Event> event, TTree *tree) -> void {
    lumi = event->GetLumiSection();
    run = event->GetRunNumber();
    evt = event->GetEventNumber();
    
    nVert = event->GetNvertices();
    vertex_x = event->GetVertex()->GetX();
    vertex_y = event->GetVertex()->GetY();
    vertex_z = event->GetVertex()->GetZ();
    nFittedHelices = (int)event->GetNhelices();
    nJet30 = event->GetNjet30();
    nJet30a = event->GetNjet30a();
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
		
		trackProcessor->SaveTracksToTree(event->GetTracks());
		
    for(int iHelix=0;iHelix<nFittedHelices;iHelix++){
      auto helix = event->GetHelix(iHelix);
      if(helix){
        helix_x[iHelix]      = helix->GetOrigin()->GetX();
        helix_y[iHelix]      = helix->GetOrigin()->GetY();
        helix_z[iHelix]      = helix->GetOrigin()->GetZ();
        helix_px[iHelix]     = helix->GetMomentum()->GetX();
        helix_py[iHelix]     = helix->GetMomentum()->GetY();
        helix_pz[iHelix]     = helix->GetMomentum()->GetZ();
        helix_charge[iHelix] = helix->GetCharge();
      }
    }
    
    leptonProcessor->SaveLeptonsToTree(event->GetLeptons());
    jetProcessor->SaveJetsToTree(event->GetJets());
    
    tree->Fill();
  };
  
  if(dataType == xtracks::kSignal){
    for(auto &event : eventsSignal[(ESignal)setIter]){func(event, tree);}
  }
  else if(dataType == xtracks::kBackground){
    for(auto &event : eventsBackground[(EBackground)setIter]){func(event, tree);}
  }
  else if(dataType == xtracks::kData){
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
    if(!config->runBackground[iBck]) continue;
    
    if(prefix==""){
      for(string path : inFileNameBackground[iBck]){
        AddEventsFromFile((path+prefix+"tree.root"),xtracks::kBackground, config->maxNeventsBackground, iBck);
      }
    }
    else{
      string path = inFileNameBackground[iBck][0];
      AddEventsFromFile((path+prefix+"tree.root"),xtracks::kBackground, config->maxNeventsBackground, iBck);
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config->runSignal[iSig]) continue;
    AddEventsFromFile((inFileNameSignal[iSig]+prefix+"tree.root"),xtracks::kSignal,config->maxNeventsSignal,iSig);
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
    
    if(prefix==""){
      for(string path : inFileNameData[iData]){
        AddEventsFromFile((path+prefix+"tree.root"),xtracks::kData, config->maxNeventsData, iData);
      }
    }
    else{
      string path = inFileNameData[iData][0];
      AddEventsFromFile((path+prefix+"tree.root"),xtracks::kData, config->maxNeventsData, iData);
    }
  }
}

void EventSet::LoadEventsFromFiles(xtracks::EDataType dataType, int setIter, string prefix)
{
  if(dataType == xtracks::kBackground){
    if(prefix==""){
      for(string path : inFileNameBackground[setIter]){
        AddEventsFromFile((path+prefix+"tree.root"),xtracks::kBackground, config->maxNeventsBackground, setIter);
      }
    }
    else{
      string path = inFileNameBackground[setIter][0];
      AddEventsFromFile((path+prefix+"tree.root"),xtracks::kBackground, config->maxNeventsBackground, setIter);
    }
  }
  else if(dataType == xtracks::kSignal){
    AddEventsFromFile((inFileNameSignal[setIter]+prefix+"tree.root"),xtracks::kSignal,config->maxNeventsSignal,setIter);
  }
  else if(dataType == xtracks::kData){
    if(prefix==""){
      for(string path : inFileNameData[setIter]){
        AddEventsFromFile((path+prefix+"tree.root"),xtracks::kData, config->maxNeventsData, setIter);
      }
    }
    else{
      string path = inFileNameData[setIter][0];
      AddEventsFromFile((path+prefix+"tree.root"),xtracks::kData, config->maxNeventsData, setIter);
    }
  }
}

void EventSet::SaveEventsToFiles(string prefix)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config->runSignal[iSig]) continue;
    system(("mkdir -p "+inFileNameSignal[iSig]+prefix).c_str());
    
    SaveToTree((inFileNameSignal[iSig]+prefix+"tree.root").c_str(), xtracks::kSignal, iSig);
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    
    // merged events will be stored in the first directory for given background
    string path = inFileNameBackground[iBck][0];
    system(("mkdir -p "+path+prefix).c_str());
    SaveToTree((path+prefix+"tree.root").c_str(), xtracks::kBackground, iBck);
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
    string path = inFileNameData[iData][0];
    system(("mkdir -p "+path+prefix).c_str());
    SaveToTree((path+prefix+"tree.root").c_str(), xtracks::kData, iData);
  }
}

void EventSet::PrintYields()
{
  double nBackgroundTotal=0;
  int nBackgroundTotalRaw=0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    nBackgroundTotal += weightedSize(xtracks::kBackground,iBck);
    nBackgroundTotalRaw += size(xtracks::kBackground,iBck);
    
    if(config->printYields){
      cout<<backgroundTitle[iBck]<<"\t";
      cout<<weightedSize(xtracks::kBackground, (int)iBck);
      cout<<"\t("<<size(xtracks::kBackground,(int)iBck)<<")"<<endl;
    }
  }
  
  if(config->printYields){
    cout<<"Background total:\t"<<nBackgroundTotal<<"\t("<<nBackgroundTotalRaw<<")"<<endl;
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!config->runSignal[iSig]) continue;
      cout<<signalTitle[iSig]<<"\tN events:\t";
      cout<<weightedSize(xtracks::kSignal,iSig);
      cout<<"\t("<<size(xtracks::kSignal,iSig)<<")"<<endl;
    }
    
    for(int iData=0;iData<kNdata;iData++){
      if(!config->runData[iData]) continue;
      cout<<dataTitle[iData]<<"\tsize:\t";
      cout<<weightedSize(xtracks::kData,iData)<<"\n";
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config->runSignal[iSig]) continue;
    cout<<signalTitle[iSig]<<"\tS/sqrt(S+B):\t";
    cout<<weightedSize(xtracks::kSignal,iSig)/sqrt(nBackgroundTotal+weightedSize(xtracks::kSignal,iSig))<<endl;
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
    cout<<dataTitle[iData]<<"\t(M-B)/sqrt(M):\t";
    cout<<(weightedSize(xtracks::kData,iData)-nBackgroundTotal)/sqrt(weightedSize(xtracks::kData,iData))<<endl;
  }
}

vector<double> EventSet::GetSignificance(bool inData)
{
  double nBackgroundTotal=0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    nBackgroundTotal += weightedSize(xtracks::kBackground,iBck);
  }
  
  vector<double> results;
  
  if(inData){
    for(int iData=0;iData<kNdata;iData++){
      if(!config->runData[iData]){
        results.push_back(inf);
        continue;
      }
      results.push_back((weightedSize(xtracks::kData,iData)-nBackgroundTotal)/sqrt(weightedSize(xtracks::kData,iData)));
    }
  }
  else{
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!config->runSignal[iSig]){
        results.push_back(inf);
        continue;
      }
      results.push_back(weightedSize(xtracks::kSignal,iSig)/sqrt(nBackgroundTotal+weightedSize(xtracks::kSignal,iSig)));
    }
  }
  return results;
}

void EventSet::ApplyCuts(const unique_ptr<EventCut> &eventCut,const unique_ptr<TrackCut> &trackCut,
                         const unique_ptr<JetCut> &jetCut,const unique_ptr<LeptonCut> &leptonCut)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config->runSignal[iSig]) continue;
    for(int iEvent=0;iEvent<(int)eventsSignal[iSig].size();){
      
      eventProcessor->ApplyTrackCut(eventsSignal[iSig][iEvent], trackCut);
      eventProcessor->ApplyJetCut(eventsSignal[iSig][iEvent], jetCut);
      eventProcessor->ApplyLeptonCut(eventsSignal[iSig][iEvent], leptonCut);
      
      if(!eventProcessor->IsPassingCut(eventsSignal[iSig][iEvent], eventCut)){
        EraseFast(eventsSignal[iSig], iEvent);
      }
      else{
        iEvent++;
      }
      
    }
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    for(int iEvent=0;iEvent<(int)eventsBackground[iBck].size();){
      
      eventProcessor->ApplyTrackCut(eventsBackground[iBck][iEvent], trackCut);
      eventProcessor->ApplyJetCut(eventsBackground[iBck][iEvent], jetCut);
      eventProcessor->ApplyLeptonCut(eventsBackground[iBck][iEvent], leptonCut);
      
      if(!eventProcessor->IsPassingCut(eventsBackground[iBck][iEvent], eventCut)){
        EraseFast(eventsBackground[iBck], iEvent);
      }
      else{
        iEvent++;
      }
      
    }
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
    for(int iEvent=0;iEvent<(int)eventsData[iData].size();){
      
      eventProcessor->ApplyTrackCut(eventsData[iData][iEvent], trackCut);
      eventProcessor->ApplyJetCut(eventsData[iData][iEvent], jetCut);
      eventProcessor->ApplyLeptonCut(eventsData[iData][iEvent], leptonCut);
      
      if(!eventProcessor->IsPassingCut(eventsData[iData][iEvent], eventCut)){
        EraseFast(eventsData[iData], iEvent);
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
  hists["nHelices"]   = new HistSet(kNhelices);
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
  hists["trackerLayers"]= new HistSet(kTrackTrackerLayers);
  hists["isolation"]    = new HistSet(kTrackRelativeIsolation);
  hists["absIsolation"] = new HistSet(kTrackAbsoluteIsolation);
  hists["trackMetDphi"] = new HistSet(kTrackMetDphi);
  hists["dedx"]         = new HistSet(kTrackDedxPerHit);
  
  hists["helixX"]       = new HistSet(kHelixX);
  hists["helixY"]       = new HistSet(kHelixY);
  hists["helixZ"]       = new HistSet(kHelixZ);
  hists["helixPx"]      = new HistSet(kHelixPx);
  hists["helixPy"]      = new HistSet(kHelixPy);
  hists["helixPz"]      = new HistSet(kHelixPz);
  hists["helixCharge"]  = new HistSet(kHelixCharge);
  
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
  
  TCanvas *helixCanvas = new TCanvas((prefix+"Helix").c_str(),(prefix+"Helix").c_str(),2880,1800);
  helixCanvas->Divide(3,3);
  hists["nHelices"]->Draw(helixCanvas,1);
  hists["helixX"]->Draw(helixCanvas,2);
  hists["helixY"]->Draw(helixCanvas,3);
  hists["helixZ"]->Draw(helixCanvas,4);
  hists["helixPx"]->Draw(helixCanvas,5);
  hists["helixPy"]->Draw(helixCanvas,6);
  hists["helixPz"]->Draw(helixCanvas,7);
  hists["helixCharge"]->Draw(helixCanvas,8);
  
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


int EventSet::size(xtracks::EDataType dataType, int setIter)
{
  if(dataType == xtracks::kSignal){
    return (int)eventsSignal[(ESignal)setIter].size();
  }
  else if(dataType == xtracks::kBackground){
    return (int)eventsBackground[(EBackground)setIter].size();
  }
  else if(dataType == xtracks::kData){
    return (int)eventsData[(EData)setIter].size();
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  return 1;
}

double EventSet::weightedSize(xtracks::EDataType dataType, int setIter)
{
  double sum=0;
  
  if(dataType == xtracks::kSignal){
    for(auto &ev : eventsSignal[(ESignal)setIter]){sum += ev->GetWeight();}
  }
  else if(dataType == xtracks::kBackground){
    for(auto &ev : eventsBackground[(EBackground)setIter]){sum += ev->GetWeight();}
  }
  else if(dataType == xtracks::kData){
    for(auto &ev : eventsData[(EData)setIter]){sum += ev->GetWeight();}
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  
  return sum;
}

shared_ptr<Event> EventSet::At(xtracks::EDataType dataType, int setIter, int index)
{
  if(dataType == xtracks::kSignal){
    return eventsSignal[(ESignal)setIter][index];
  }
  else if(dataType == xtracks::kBackground){
    return eventsBackground[(EBackground)setIter][index];
  }
  else if(dataType == xtracks::kData){
    return eventsData[(EData)setIter][index];
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
}

shared_ptr<Event> EventSet::GetEvent(xtracks::EDataType dataType, uint run, uint lumi, unsigned long long eventNumber)
{
  if(dataType == xtracks::kSignal){
    for(auto events : eventsSignal){
      for(auto event : events){
        if(event->GetRunNumber() == run &&
           event->GetLumiSection() == lumi &&
           event->GetEventNumber() == eventNumber){
          return event;
        }
      }
    }
  }
  else if(dataType == xtracks::kBackground){
    for(auto events : eventsBackground){
      for(auto event : events){
        if(event->GetRunNumber() == run &&
           event->GetLumiSection() == lumi &&
           event->GetEventNumber() == eventNumber){
          return event;
        }
      }
    }
  }
  else if(dataType == xtracks::kData){
    for(auto events : eventsData){
      for(auto event : events){
        if(event->GetRunNumber() == run &&
           event->GetLumiSection() == lumi &&
           event->GetEventNumber() == eventNumber){
          return event;
        }
      }
    }
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  return nullptr;
}

void EventSet::AddEvent(shared_ptr<Event> event, xtracks::EDataType dataType, int setIter)
{
  
  if(dataType == xtracks::kSignal){
    eventsSignal[(ESignal)setIter].push_back(event);
  }
  else if(dataType == xtracks::kBackground){
    eventsBackground[(EBackground)setIter].push_back(event);
  }
  else if(dataType == xtracks::kData){
    eventsData[(EData)setIter].push_back(event);
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
}

void EventSet::AddEventsFromFile(std::string fileName, xtracks::EDataType dataType, int maxNevents, int setIter)
{
  cout<<"Reading events from:"<<fileName<<endl;
  TFile *inFile = TFile::Open(fileName.c_str());
  TTree *tree = (TTree*)inFile->Get("tree");
  TTreeReader reader("tree", inFile);
  
  trackProcessor->SetupBranchesForReading(tree);
  jetProcessor->SetupBranchesForReading(tree);
  leptonProcessor->SetupBranchesForReading(tree);
  helixProcessor->SetupBranches(tree);
  
  TTreeReaderValue<uint>   _run(reader, "run");
  TTreeReaderValue<uint>   _lumi(reader, "lumi");
  TTreeReaderValue<unsigned long long>   _evt(reader, "evt");
  
  TTreeReaderValue<int>   _nTracks(reader, "nIsoTrack");
  TTreeReaderValue<int>   _nVert(reader, "nVert");
  TTreeReaderValue<float> _vertex_x(reader, "vertex_x");
  TTreeReaderValue<float> _vertex_y(reader, "vertex_y");
  TTreeReaderValue<float> _vertex_z(reader, "vertex_z");
  TTreeReaderValue<int>   _nJets(reader, "nJet");
  TTreeReaderValue<int>   _nJetsFwd(reader, "nJetFwd");
  TTreeReaderValue<int>   _nJet30(reader, "nJet30");
  TTreeReaderValue<int>   _nJet30a(reader, "nJet30a");
  TTreeReaderValue<int>   _nLepton(reader, "nLepGood");
  TTreeReaderValue<int>   _nTau(reader, "nTauGood");
  TTreeReaderValue<int>   _nGenChargino(reader, "nGenChargino");
  
  TTreeReaderValue<float> _xSec  (reader,(tree->GetBranchStatus("xsec")) ? "xsec" : "rho");
  TTreeReaderValue<float> _sumWgt(reader,(tree->GetBranchStatus("wgtsum")) ? "wgtsum" : "rho");
  TTreeReaderValue<float> _genWgt(reader,(tree->GetBranchStatus("genWeight")) ? "genWeight" : "rho");
  
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
  
  int iter=0;
  
  while(reader.Next()){
    if(maxNevents>0 && iter>maxNevents) break;
    
    shared_ptr<Event> newEvent = shared_ptr<Event>(new Event());
    
    tree->GetEntry(iter++);
    vector<shared_ptr<Track>> tracks = trackProcessor->GetTracksFromTree();
    
    for(int iTrack=0;iTrack<*_nTracks;iTrack++){
      auto track = tracks[iTrack];
      
      track->SetEventMetPt(*_metPt);
      track->SetEventMetEta(*_metEta);
      track->SetEventMetPhi(*_metPhi);
      track->SetEventMetMass(*_metMass);
      
      // Decay point is set only for random pion generation. It's not used in the fitter!!
      double minR = layerR[track->GetLastBarrelLayer()];
      double maxR = layerR[track->GetLastBarrelLayer()+1];
      double decayR = (maxR+minR)/2.;
      
      track->SetDecayPoint(make_unique<Point>(decayR*cos(track->GetPhi()) + *_vertex_x*10,
                                              decayR*sin(track->GetPhi()) + *_vertex_y*10,
                                              decayR/sin(track->GetTheta())*cos(track->GetTheta())+*_vertex_z*10)
                           );
      
      newEvent->AddTrack(track);
    }
    
    vector<shared_ptr<Helix>> helices = helixProcessor->GetHelicesFromTree();
    
    for(auto helix : helices){
      newEvent->AddHelix(helix);
    }

    vector<shared_ptr<Jet>> jets = jetProcessor->GetJetsFromTree();
    
    for(auto jet : jets){
      newEvent->AddJet(jet);
    }
    
    vector<shared_ptr<Lepton>> leptons = leptonProcessor->GetLeptonsFromTree();
    
    for(auto lepton : leptons){
      newEvent->AddLepton(lepton);
    }
    
    double lumi = config->totalLuminosity * 1000.; // transform from fb^-1 to pb^-1
    double weight = lumi * (*_genWgt) / (*_sumWgt);
    
    //    static map<string,set<double>> wgts;
    //
    //    if(*_genWgt != 1.0 && wgts[fileName].find(*_genWgt) == wgts[fileName].end() ){
    //      wgts[fileName].insert(*_genWgt);
    //      cout<<*_genWgt<<"\t"<<fileName<<endl;
    //    }
    
    if(dataType==xtracks::kBackground){
      weight *= (*_xSec);
    }
    if(dataType==xtracks::kSignal){
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
    else if(dataType==xtracks::kData){
      weight = 1;
    }
    
    newEvent->SetLumiSection(*_lumi);
    newEvent->SetRunNumber(*_run);
    newEvent->SetEventNumber(*_evt);
    newEvent->SetWeight(weight);
    
    newEvent->SetNvertices(*_nVert);
    newEvent->SetVertex(make_unique<Point>(*_vertex_x, *_vertex_y, *_vertex_z));
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
    
    newEvent->SetDataType(dataType);
    newEvent->SetSetIter(setIter);
    
    if(dataType == xtracks::kSignal){
      eventsSignal[setIter].push_back(newEvent);
    }
    else if(dataType == xtracks::kBackground){
      eventsBackground[setIter].push_back(newEvent);
    }
    else if(dataType == xtracks::kData){
      eventsData[setIter].push_back(newEvent);
    }
    else{
      throw out_of_range("Unknown data type provided");
    }
    
  }
}
