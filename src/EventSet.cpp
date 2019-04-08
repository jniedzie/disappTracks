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
jetProcessor(make_unique<JetProcessor>()),
helixProcessor(make_unique<HelixProcessor>()),
eventProcessor(make_unique<EventProcessor>()),
leptonProcessor(make_unique<LeptonProcessor>())
{
  AddEventsFromFile(fileName,dataType,maxNevents,iSig);
}

EventSet::EventSet(const EventSet &e) :
trackProcessor(make_unique<TrackProcessor>()),
jetProcessor(make_unique<JetProcessor>()),
helixProcessor(make_unique<HelixProcessor>()),
eventProcessor(make_unique<EventProcessor>()),
leptonProcessor(make_unique<LeptonProcessor>())
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

void EventSet::SaveToTree(string fileName, xtracks::EDataType dataType, int setIter)
{
  TFile outFile(fileName.c_str(),"RECREATE");
  outFile.cd();
  TTree *tree = new TTree("tree","tree");
	
  eventProcessor->SetupBranchesForWriting(tree);
	trackProcessor->SetupBranchesForWriting(tree);
  jetProcessor->SetupBranchesForWriting(tree);
  leptonProcessor->SetupBranchesForWriting(tree);
  helixProcessor->SetupBranchesForWriting(tree);
  
  function<void(shared_ptr<Event>, TTree*)> func = [&](shared_ptr<Event> event, TTree *tree) -> void {
    
    eventProcessor->SaveEventToTree(event);
		trackProcessor->SaveTracksToTree(event->GetTracks());
    leptonProcessor->SaveLeptonsToTree(event->GetLeptons());
    jetProcessor->SaveJetsToTree(event->GetJets());
    helixProcessor->SaveHelicesToTree(event->GetHelices());
    
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
        results.push_back(-inf);
        continue;
      }
      results.push_back((weightedSize(xtracks::kData,iData)-nBackgroundTotal)/sqrt(weightedSize(xtracks::kData,iData)));
    }
  }
  else{
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!config->runSignal[iSig]){
        results.push_back(-inf);
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
  hists["dedx"]->Draw(canvasTrack,2);
  hists["caloEm"]->Draw(canvasTrack,3);
  hists["caloHad"]->Draw(canvasTrack,4);
  hists["eta"]->Draw(canvasTrack,5);
  hists["phi"]->Draw(canvasTrack,6);
  hists["trackMetDphi"]->Draw(canvasTrack,7);
  hists["trackerLayers"]->Draw(canvasTrack,8);
  hists["dz"]->Draw(canvasTrack,9);
  hists["isolation"]->Draw(canvasTrack,10);
  hists["absIsolation"]->Draw(canvasTrack,11);
  
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
  
  eventProcessor->SetupBranchesForReading(tree);
  trackProcessor->SetupBranchesForReading(tree);
  jetProcessor->SetupBranchesForReading(tree);
  leptonProcessor->SetupBranchesForReading(tree);
  helixProcessor->SetupBranchesForReading(tree);
  
  for(int iEntry=0;iEntry<tree->GetEntries();iEntry++){
    if(maxNevents>0 && iEntry>maxNevents) break;
    tree->GetEntry(iEntry);
    
    auto event = eventProcessor->GetEventFromTree(dataType, setIter);
    
    vector<shared_ptr<Track>> tracks = trackProcessor->GetTracksFromTree();
    
    for(int iTrack=0;iTrack<tracks.size();iTrack++){
      auto track = tracks[iTrack];
      
      track->SetEventMetPt(event->GetMetPt());
      track->SetEventMetEta(event->GetMetEta());
      track->SetEventMetPhi(event->GetMetPhi());
      track->SetEventMetMass(event->GetMetMass());
      
      // Decay point is set only for random pion generation. It's not used in the fitter!!
      double minR = layerR[track->GetNtrackerLayers()-1];
      double maxR = layerR[track->GetNtrackerLayers()];
      double decayR = (maxR+minR)/2.;
      
      track->SetDecayPoint(Point(decayR*cos(track->GetPhi()) + 10*event->GetVertex()->GetX(),
                                 decayR*sin(track->GetPhi()) + 10*event->GetVertex()->GetY(),
                                 decayR/sin(track->GetTheta())*cos(track->GetTheta()) + 10*event->GetVertex()->GetZ())
                           );
      
      event->AddTrack(track);
    }
    
    vector<shared_ptr<Helix>> helices = helixProcessor->GetHelicesFromTree();
    
    for(auto helix : helices){
      event->AddHelix(helix);
    }

    vector<shared_ptr<Jet>> jets = jetProcessor->GetJetsFromTree();
    
    for(auto jet : jets){
      event->AddJet(jet);
    }
    
    vector<shared_ptr<Lepton>> leptons = leptonProcessor->GetLeptonsFromTree();
    
    for(auto lepton : leptons){
      event->AddLepton(lepton);
    }
    
   
    if(dataType == xtracks::kSignal)          eventsSignal[setIter].push_back(event);
    else if(dataType == xtracks::kBackground) eventsBackground[setIter].push_back(event);
    else if(dataType == xtracks::kData)       eventsData[setIter].push_back(event);
    else                                      throw out_of_range("Unknown data type provided");
  }
}
