//  EventSet.cpp
//
//  Created by Jeremi Niedziela on 08/11/2018.

#include "EventSet.hpp"
#include "HistSet.hpp"
#include "Logger.hpp"

EventSet::EventSet()
{
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(ESignal iSig : signals)         eventsSignal[iSig][year]      = vector<shared_ptr<Event>>();
    for(EBackground iBck : backgrounds) eventsBackground[iBck][year]  = vector<shared_ptr<Event>>();
    for(EData iData : datas)            eventsData[iData][year]       = vector<shared_ptr<Event>>();
  }
}

EventSet::EventSet(string fileName, xtracks::EDataType dataType, int year, int maxNevents, ESignal iSig)
{
  AddEventsFromFile(fileName, dataType, year, maxNevents, iSig);
}

EventSet::EventSet(const EventSet &e)
{
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(ESignal iSig : signals){
      eventsSignal[iSig][year] = vector<shared_ptr<Event>>();
      for(auto &event : e.eventsSignal.at(iSig).at(year)){
        eventsSignal[iSig][year].push_back(make_shared<Event>(*event));
      }
    }
    for(EBackground iBck : backgrounds){
      eventsBackground[iBck][year] = vector<shared_ptr<Event>>();
      for(auto &event : e.eventsBackground.at(iBck).at(year)){
        eventsBackground[iBck][year].push_back(make_shared<Event>(*event));
      }
    }
    for(EData iData : datas){
      eventsData[iData][year] = vector<shared_ptr<Event>>();
      for(auto &event : e.eventsData.at(iData).at(year)){
        eventsData[iData][year].push_back(make_shared<Event>(*event));
      }
    }
  }
}

EventSet EventSet::operator=(const EventSet &e)
{
  return EventSet(e);
}

EventSet::~EventSet()
{
  
}

void EventSet::SaveToTree(string fileName, xtracks::EDataType dataType, int setIter, int year) const
{
  if(find(years.begin(), years.end(), year) == years.end()){
    cout<<"ERROR -- privided year "<<year<<" not in years vector!!"<<endl;
  }
  
  TFile outFile(fileName.c_str(),"RECREATE");
  outFile.cd();
  TTree *tree = new TTree("tree","tree");
	
  eventProcessor.SetupBranchesForWriting(tree);
  trackProcessor.SetupBranchesForWriting(tree);
  leptonProcessor.SetupBranchesForWriting(tree);
  jetProcessor.SetupBranchesForWriting(tree);
  helixProcessor.SetupBranchesForWriting(tree);
  
  function<void(shared_ptr<Event>, TTree*)> func = [&](shared_ptr<Event> event, TTree *tree) -> void {
    eventProcessor.SaveEventToTree(event);
    trackProcessor.SaveTracksToTree(event->GetTracks());
    leptonProcessor.SaveLeptonsToTree(event->GetLeptons());
    jetProcessor.SaveJetsToTree(event->GetJets());
    helixProcessor.SaveHelicesToTree(event->GetHelices());
    
    tree->Fill();
  };
  
  if(dataType == xtracks::kSignal){
    for(auto &event : eventsSignal.at((ESignal)setIter).at(year)){func(event, tree);}
  }
  else if(dataType == xtracks::kBackground){
    for(auto &event : eventsBackground.at((EBackground)setIter).at(year)){func(event, tree);}
  }
  else if(dataType == xtracks::kData){
    for(auto &event : eventsData.at((EData)setIter).at(year)){func(event, tree);}
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  
  tree->Write();
}

void EventSet::LoadEventsFromFiles(string prefix)
{
  for(EBackground iBck : backgrounds){
    if(!config.runBackground[iBck]) continue;
    LoadEventsFromFiles(xtracks::kBackground, iBck, prefix);
  }
  for(ESignal iSig : signals){
    if(!config.runSignal[iSig]) continue;
    LoadEventsFromFiles(xtracks::kSignal, iSig, prefix);
  }
  for(EData iData : datas){
    if(!config.runData[iData]) continue;
    LoadEventsFromFiles(xtracks::kData, iData, prefix);
  }
}

void EventSet::LoadEventsFromFiles(xtracks::EDataType dataType, int setIter, string prefix, int iEvent)
{
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    pair<string, vector<string>> paths;
    int maxEvents = -1;
    
    if(dataType == xtracks::kBackground){
      paths = inFileNameBackground.at((EBackground)setIter).at(year);
      maxEvents = config.params["max_N_events_background"];
    }
    else if(dataType == xtracks::kSignal){
      paths = inFileNameSignal.at((ESignal)setIter).at(year);
      maxEvents = config.params["max_N_events_signal"];
    }
    else if(dataType == xtracks::kData){
      paths = inFileNameData.at((EData)setIter).at(year);
      maxEvents = config.params["max_N_events_data"];
    }
    
    
    if(config.params["load_single_subpath"]){
      if(prefix == ""){
        string fullPath = baseDataPath.at(year) + paths.first + paths.second[config.params["subpath_index"]] + commonDataSuffix + "tree.root";
        AddEventsFromFile(fullPath, dataType, year, maxEvents, setIter, iEvent);
      }
      else{
        for(string path : paths.second){
          string fullPath = baseDataPath.at(year) + paths.first + path + commonDataSuffix + prefix + "tree.root";
          AddEventsFromFile(fullPath, dataType, year, maxEvents, setIter, iEvent);
          
          if(prefix.find("after_L0/") == string::npos) break;
        }
      }
    }
    else{
      for(string path : paths.second){
        string fullPath = baseDataPath.at(year) + paths.first + path + commonDataSuffix + prefix + "tree.root";
        AddEventsFromFile(fullPath, dataType, year, maxEvents, setIter, iEvent);
        
        if(prefix != "") break;
      }
    }
  }
}

void EventSet::SaveEventsToFiles(string prefix) const
{
  int subpathIndex = 0;
  if(config.params["load_single_subpath"]) subpathIndex = config.params["subpath_index"];
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(EBackground iBck : backgrounds){
      if(!config.runBackground[iBck]) continue;
      
      auto &[basePath, paths] = inFileNameBackground.at(iBck).at(year);
      if(paths.size() == 0) continue;
      
      string path = baseDataPath.at(year) + basePath + paths[subpathIndex] + commonDataSuffix;
      system(("mkdir -p "+path+prefix).c_str());
      SaveToTree((path+prefix+"tree.root").c_str(), xtracks::kBackground, iBck, year);
    }
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      
      auto &[basePath, paths] = inFileNameSignal.at(iSig).at(year);
      if(paths.size() == 0) continue;
      string path = baseDataPath.at(year) + basePath + paths[subpathIndex] + commonDataSuffix;
      system(("mkdir -p "+path+prefix).c_str());
      SaveToTree((path+prefix+"tree.root").c_str(), xtracks::kSignal, iSig, year);
    }
    for(EData iData : datas){
      if(!config.runData[iData]) continue;
      
      auto &[basePath, paths] = inFileNameData.at(iData).at(year);
      if(paths.size() == 0) continue;
      string path = baseDataPath.at(year) + basePath + paths[subpathIndex] + commonDataSuffix;
      system(("mkdir -p "+path+prefix).c_str());
      SaveToTree((path+prefix+"tree.root").c_str(), xtracks::kData, iData, year);
    }
    
  }
}

void EventSet::PrintYields() const
{
  double nBackgroundTotal=0;
  int nBackgroundTotalRaw=0;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(EBackground iBck : backgrounds){
      if(!config.runBackground[iBck]) continue;

      nBackgroundTotal += weightedSize(xtracks::kBackground, iBck, year);
      nBackgroundTotalRaw += size(xtracks::kBackground, iBck, year);
      
      if(config.params["print_yields"]){
        cout<<backgroundTitle.at(iBck)<<"("<<year<<")\t";
        cout<<weightedSize(xtracks::kBackground, (int)iBck, year);
        cout<<"\t("<<size(xtracks::kBackground,(int)iBck, year)<<")"<<endl;
      }
    }
  }
  
  if(config.params["print_yields"]){
    cout<<"Background total:\t"<<nBackgroundTotal<<"\t("<<nBackgroundTotalRaw<<")"<<endl;
    for(int year : years){
      if(!config.params["load_"+to_string(year)]) continue;
      
      for(ESignal iSig : signals){
        if(!config.runSignal[iSig]) continue;
        cout<<signalTitle.at(iSig)<<"\tN events:\t";
        cout<<weightedSize(xtracks::kSignal, iSig, year);
        cout<<"\t("<<size(xtracks::kSignal, iSig, year)<<")"<<endl;
      }
      
      for(EData iData : datas){
        if(!config.runData[iData]) continue;
        cout<<dataTitle.at(iData)<<"\tsize:\t";
        cout<<weightedSize(xtracks::kData, iData, year)<<"\n";
      }
    }
  }
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      cout<<signalTitle.at(iSig)<<"\tS/sqrt(S+B):\t";
      cout<<weightedSize(xtracks::kSignal, iSig, year)/sqrt(nBackgroundTotal+weightedSize(xtracks::kSignal, iSig, year))<<endl;
    }
    
    for(EData iData : datas){
      if(!config.runData[iData]) continue;
      cout<<dataTitle.at(iData)<<"\t(M-B)/sqrt(M):\t";
      cout<<(weightedSize(xtracks::kData, iData, year)-nBackgroundTotal)/sqrt(weightedSize(xtracks::kData, iData, year))<<endl;
    }
    
  }
}

vector<double> EventSet::GetSignificance(bool inData) const
{
  double nBackgroundTotal=0;
  
  for(EBackground iBck : backgrounds){
    if(!config.runBackground[iBck]) continue;
    for(int year : years){
      if(!config.params["load_"+to_string(year)]) continue;
      nBackgroundTotal += weightedSize(xtracks::kBackground, iBck, year);
    }
  }
  
  vector<double> results;
  
  if(inData){
    for(EData iData : datas){
      if(!config.runData[iData]){
        results.push_back(-inf);
        continue;
      }
      for(int year : years){
        if(!config.params["load_"+to_string(year)]) continue;
        results.push_back((weightedSize(xtracks::kData, iData, year) - nBackgroundTotal)/sqrt(weightedSize(xtracks::kData, iData, year)));
      }
    }
  }
  else{
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]){
        results.push_back(-inf);
        continue;
      }
      for(int year : years){
        if(!config.params["load_"+to_string(year)]) continue;
        results.push_back(weightedSize(xtracks::kSignal, iSig, year)/sqrt(nBackgroundTotal+weightedSize(xtracks::kSignal, iSig, year)));
      }
    }
  }
  return results;
}

void EventSet::ApplyCuts(const EventCut   &eventCut,
                         const TrackCut   &trackCut,
                         const JetCut     &jetCut,
                         const LeptonCut  &leptonCut)
{
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      
      vector<int> *cutReasons = new vector<int>(100);
      
      for(int iEvent=0;iEvent<(int)eventsSignal[iSig][year].size();){
        
        eventProcessor.ApplyTrackCut(eventsSignal[iSig][year][iEvent], trackCut);
        eventProcessor.ApplyJetCut(eventsSignal[iSig][year][iEvent], jetCut);
        eventProcessor.ApplyLeptonCut(eventsSignal[iSig][year][iEvent], leptonCut);
        
        if(!eventProcessor.IsPassingCut(eventsSignal[iSig][year][iEvent], eventCut, cutReasons)){
          EraseFast(eventsSignal[iSig][year], iEvent);
        }
        else{
          iEvent++;
        }
      }
      
      Log(1)<<"Cut-through for signal sample: "<<signalTitle.at(iSig)<<"\n";
      int iter=0;
      for(int nEventsPassing : *cutReasons){
        if(nEventsPassing!=0) Log(1)<<nEventsPassing<<" events passing cut "<<iter++<<"\n";
      }
      
    }
    
    for(EBackground iBck : backgrounds){
      if(!config.runBackground[iBck]) continue;
      
      vector<int> *cutReasons = new vector<int>(100);
      
      for(int iEvent=0;iEvent<(int)eventsBackground[iBck][year].size();){
        
        eventProcessor.ApplyTrackCut(eventsBackground[iBck][year][iEvent], trackCut);
        eventProcessor.ApplyJetCut(eventsBackground[iBck][year][iEvent], jetCut);
        eventProcessor.ApplyLeptonCut(eventsBackground[iBck][year][iEvent], leptonCut);
        
        if(!eventProcessor.IsPassingCut(eventsBackground[iBck][year][iEvent], eventCut, cutReasons)){
          EraseFast(eventsBackground[iBck][year], iEvent);
        }
        else{
          iEvent++;
        }
      }
      
      Log(1)<<"Cut-through for background sample: "<<backgroundTitle.at(iBck)<<"\n";
      int iter=0;
      for(int nEventsPassing : *cutReasons){
        if(nEventsPassing!=0) Log(1)<<nEventsPassing<<" events passing cut "<<iter++<<"\n";
      }
      
      auto survivors = eventProcessor.survivingEvents;
      
      vector<pair<uint, unsigned long long>> lumi_event;
      
      for(auto event : survivors){
        lumi_event.push_back(make_pair(event->GetLumiSection(), event->GetEventNumber()));
      }
      
      sort(lumi_event.begin(), lumi_event.end());
      
      //    for(auto &[lumi, event] : lumi_event){
      //      cout<<lumi<<":"<<event<<endl;
      //    }
    }
    
    for(EData iData : datas){
      if(!config.runData[iData]) continue;
      
      vector<int> *cutReasons = new vector<int>(100);
      
      for(int iEvent=0;iEvent<(int)eventsData[iData][year].size();){
        
        eventProcessor.ApplyTrackCut(eventsData[iData][year][iEvent], trackCut);
        eventProcessor.ApplyJetCut(eventsData[iData][year][iEvent], jetCut);
        eventProcessor.ApplyLeptonCut(eventsData[iData][year][iEvent], leptonCut);
        
        if(!eventProcessor.IsPassingCut(eventsData[iData][year][iEvent], eventCut, cutReasons)){
          EraseFast(eventsData[iData][year], iEvent);
        }
        else{
          iEvent++;
        }
      }
      
      Log(1)<<"Cut-through for data sample: "<<dataTitle.at(iData)<<"\n";
      int iter=0;
      for(int nEventsPassing : *cutReasons){
        if(nEventsPassing!=0) Log(1)<<nEventsPassing<<" events passing cut "<<iter++<<"\n";
      }
      
    }
  }
}

void EventSet::DrawStandardPlots(string prefix) const
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
  hists["minDedx"]      = new HistSet(kTrackMinDedx);
  
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
  
  string outPath;
  if(config.params["cuts_level"]==0)       outPath = "results/plots_after_L0/";
  else if(config.params["cuts_level"]==1)  outPath = "results/plots_after_L1/";
  else{
    cout<<"ERROR -- unknown cuts level: "<<config.params["cuts_level"]<<endl;
    exit(0);
  }
  
    
  // Plot histograms
  TCanvas *canvasEvents = new TCanvas((prefix+"Events").c_str(),(prefix+"Events").c_str(),800,1200);
  canvasEvents->Divide(2,3);
  
  hists["nVertices"]->Draw(canvasEvents,1);
  hists["nJet"]->Draw(canvasEvents,2);
  hists["nIsoTrack"]->Draw(canvasEvents,3);
  hists["nMetPt"]->Draw(canvasEvents,4);
  hists["nMetJetDphi"]->Draw(canvasEvents,5);
  hists["trackMetDphi"]->Draw(canvasEvents,6);
  
  TCanvas *canvasTrack = new TCanvas((prefix+"Tracks").c_str(),(prefix+"Tracks").c_str(),800,1200);
  canvasTrack->Divide(2,3);
  
  hists["pt"]->Draw(canvasTrack,1);
  hists["trackerLayers"]->Draw(canvasTrack,2);
  hists["minDedx"]->Draw(canvasTrack,3);
  hists["isolation"]->Draw(canvasTrack,4);
  hists["caloEm"]->Draw(canvasTrack,5);
  hists["caloHad"]->Draw(canvasTrack,6);
  
//  hists["eta"]->Draw(canvasTrack,2);
//  hists["phi"]->Draw(canvasTrack,3);
//  hists["dz"]->Draw(canvasTrack,9);
//  hists["absIsolation"]->Draw(canvasTrack,11);
  
  TCanvas *canvasHelix = new TCanvas((prefix+"Helix").c_str(),(prefix+"Helix").c_str(),2880,1800);
  canvasHelix->Divide(3,3);
  hists["nHelices"]->Draw(canvasHelix,1);
  hists["helixX"]->Draw(canvasHelix,2);
  hists["helixY"]->Draw(canvasHelix,3);
  hists["helixZ"]->Draw(canvasHelix,4);
  hists["helixPx"]->Draw(canvasHelix,5);
  hists["helixPy"]->Draw(canvasHelix,6);
  hists["helixPz"]->Draw(canvasHelix,7);
  hists["helixCharge"]->Draw(canvasHelix,8);
  
  TCanvas *canvasJets = new TCanvas("Jets","Jets",800,1200);
  canvasJets->Divide(2,3);
  
  hists["jet_pt"]->Draw(canvasJets, 1);
  hists["jet_eta"]->Draw(canvasJets, 2);
  hists["jet_phi"]->Draw(canvasJets, 3);
  hists["jetTrackDr"]->Draw(canvasJets, 4);
  hists["jetCHF"]->Draw(canvasJets, 5);
  hists["jetNHF"]->Draw(canvasJets, 6);
  
  canvasEvents->SaveAs((outPath+"canvas_events.pdf").c_str());
  canvasTrack->SaveAs((outPath+"canvas_tracks.pdf").c_str());
  canvasJets->SaveAs((outPath+"canvas_jets.pdf").c_str());
  canvasHelix->SaveAs((outPath+"canvas_helices.pdf").c_str());
}

void EventSet::DrawPerLayerPlots() const
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


int EventSet::size(xtracks::EDataType dataType, int setIter, int year) const
{
  if(find(years.begin(), years.end(), year) == years.end()){
    cout<<"ERROR -- privided year "<<year<<" not in years vector!!"<<endl;
  }
  
  if(dataType == xtracks::kSignal)          return (int)eventsSignal.at((ESignal)setIter).at(year).size();
  else if(dataType == xtracks::kBackground) return (int)eventsBackground.at((EBackground)setIter).at(year).size();
  else if(dataType == xtracks::kData)       return (int)eventsData.at((EData)setIter).at(year).size();
  else throw out_of_range("Unknown data type provided");
  
  return 1;
}

double EventSet::weightedSize(xtracks::EDataType dataType, int setIter, int year) const
{
  if(find(years.begin(), years.end(), year) == years.end()){
    cout<<"ERROR -- privided year "<<year<<" not in years vector!!"<<endl;
  }
  
  double sum=0;
  
  if(dataType == xtracks::kSignal){
    for(auto &ev : eventsSignal.at((ESignal)setIter).at(year)){sum += ev->GetWeight();}
  }
  else if(dataType == xtracks::kBackground){
    for(auto &ev : eventsBackground.at((EBackground)setIter).at(year)){sum += ev->GetWeight();}
  }
  else if(dataType == xtracks::kData){
    for(auto &ev : eventsData.at((EData)setIter).at(year)){sum += ev->GetWeight();}
  }
  else{
    throw out_of_range("Unknown data type provided");
  }
  
  return sum;
}

shared_ptr<Event> EventSet::At(xtracks::EDataType dataType, int setIter, int year, int index) const
{
  if(find(years.begin(), years.end(), year) == years.end()){
    cout<<"ERROR -- privided year "<<year<<" not in years vector!!"<<endl;
  }
  
  if(dataType == xtracks::kSignal)          return eventsSignal.at((ESignal)setIter).at(year)[index];
  else if(dataType == xtracks::kBackground) return eventsBackground.at((EBackground)setIter).at(year)[index];
  else if(dataType == xtracks::kData)       return eventsData.at((EData)setIter).at(year)[index];
  else  throw out_of_range("Unknown data type provided");
}

shared_ptr<Event> EventSet::GetEvent(xtracks::EDataType dataType,
                                     uint run, uint lumi, unsigned long long eventNumber) const
{
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    if(dataType == xtracks::kSignal){
      for(ESignal iSig : signals){
        for(auto event : eventsSignal.at(iSig).at(year)){
          if(event->GetRunNumber() == run &&
             event->GetLumiSection() == lumi &&
             event->GetEventNumber() == eventNumber){
            return event;
          }
        }
      }
    }
    else if(dataType == xtracks::kBackground){
      for(EBackground iBck : backgrounds){
        for(auto event : eventsBackground.at(iBck).at(year)){
          if(event->GetRunNumber() == run &&
             event->GetLumiSection() == lumi &&
             event->GetEventNumber() == eventNumber){
            return event;
          }
        }
      }
    }
    else if(dataType == xtracks::kData){
      for(EData iData : datas){
        for(auto event : eventsData.at(iData).at(year)){
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
  }
  return nullptr;
}

void EventSet::AddEvent(shared_ptr<Event> event, xtracks::EDataType dataType, int setIter, int year)
{
  if(find(years.begin(), years.end(), year) == years.end()){
    cout<<"ERROR -- privided year "<<year<<" not in years vector!!"<<endl;
  }
  
  if(dataType == xtracks::kSignal)          eventsSignal.at((ESignal)setIter).at(year).push_back(event);
  else if(dataType == xtracks::kBackground) eventsBackground.at((EBackground)setIter).at(year).push_back(event);
  else if(dataType == xtracks::kData)       eventsData.at((EData)setIter).at(year).push_back(event);
  else throw out_of_range("Unknown data type provided");
}

void EventSet::AddEventsFromFile(string fileName, xtracks::EDataType dataType, int year,
                                 int maxNevents, int setIter, int iEvent)
{
  if(find(years.begin(), years.end(), year) == years.end()){
    cout<<"ERROR -- privided year "<<year<<" not in years vector!!"<<endl;
  }
  
  cout<<"Reading events from:"<<fileName<<endl;
  TFile *inFile = TFile::Open(fileName.c_str());
  TTree *tree = (TTree*)inFile->Get("tree");
  
  if(!tree){
    Log(0)<<"Error -- could not load tree from file: "<<fileName<<"\n";
    return;
  }
  
  TTreeReader reader("tree", inFile);
  
  string basePath;
  pair<string, vector<string>> paths;
  
  if(dataType == xtracks::kSignal)      paths = inFileNameSignal.at((ESignal)setIter).at(year);
  if(dataType == xtracks::kBackground)  paths = inFileNameBackground.at((EBackground)setIter).at(year);
  if(dataType == xtracks::kData)        paths = inFileNameData.at((EData)setIter).at(year);

//  basePath = paths.first + paths.second[0];
  basePath = baseDataPath.at(year) + paths.first + paths.second[config.params["subpath_index"]];
  
  TFile *inFileFriend = nullptr;
  TTree *treeFriend = nullptr;
  
  if(config.params["load_friend_tree"]){
    cout<<"Opening friend file"<<endl;

    inFileFriend = TFile::Open(Form("%s/%s/tree_friend.root",basePath.c_str(), commonDataSuffix.c_str()));
    
    if(!inFileFriend){
      cout<<"WARNING -- no friend file was found in path: "<<basePath<<"/tree_friend.root"<<endl;
      cout<<"WARNING -- some additional info will not be available"<<endl;
    }
    else{
      cout<<"Reading friend tree...";
      treeFriend = (TTree*)inFileFriend->Get("CharginoAnalyzer/tree");
      cout<<" done"<<endl;
    }
    if(!treeFriend) cout<<"ERROR -- Could not find friend tree in the friend file!"<<endl;
  }
  
  TTree *prefireTree = nullptr;
  if(year == 2017 && dataType==kSignal){
    cout<<"Opening prefire file"<<endl;
    string prefirePath = basePath + "/" + commonDataSuffix + "/tree_prefire.root";
    TFile *inFilePrefire = TFile::Open(prefirePath.c_str());
    
    if(!inFilePrefire){
      cout<<"WARNING -- no prefire file was found in path: "<<prefirePath<<endl;
    }
    else{
      cout<<"Reading prefire tree...";
      prefireTree = (TTree*)inFilePrefire->Get("tree");
      if(!prefireTree)  cout<<"ERROR -- Could not find prefire tree in the prefire file!"<<endl;
      else              cout<<" done"<<endl;
    }
  }
  TH1D *metWeights = nullptr;
  if(setIter == kChargino300_1  || setIter == kChargino300_10 ||
     setIter == kChargino400_1  ||
     setIter == kChargino500_1  || setIter == kChargino500_10 ||
     setIter == kChargino600_10 ||
     setIter == kChargino700_10 || setIter == kChargino700_30 ||
     setIter == kChargino800_10 || setIter == kChargino800_30 ||
     setIter == kChargino900_30){
    
    TFile *metFile = TFile::Open((basePath+"../metWeights.root").c_str());
    metWeights = (TH1D*)metFile->Get("metRatio");
  }
  
  eventProcessor.SetupBranchesForReading(tree, treeFriend, prefireTree);
  trackProcessor.SetupBranchesForReading(tree);
  jetProcessor.SetupBranchesForReading(tree);
  leptonProcessor.SetupBranchesForReading(tree);
  helixProcessor.SetupBranchesForReading(tree);
  
  cout<<"Loading events"<<endl;
  for(int iEntry=0;iEntry<tree->GetEntries();iEntry++){
    if(iEvent >= 0){
      tree->GetEntry(iEvent);
    }
    else{
      if(maxNevents>0 && iEntry>=maxNevents) break;
      tree->GetEntry(iEntry);
    }
    if(iEvent%100000 == 0)  Log(1)<<"Events loaded: "<<iEvent<<"\n";
    if(iEntry%100 == 0)     Log(2)<<"Events loaded: "<<iEntry<<"\n";

    auto event = eventProcessor.GetEventFromTree(dataType, setIter, year, treeFriend, prefireTree, metWeights);
    
    vector<shared_ptr<Track>> tracks = trackProcessor.GetTracksFromTree();
    
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
    
    Helices helices = helixProcessor.GetHelicesFromTree();
    for(auto helix : helices) event->AddHelix(helix);

    vector<shared_ptr<Jet>> jets = jetProcessor.GetJetsFromTree();
    for(auto jet : jets) event->AddJet(jet);
    
    vector<shared_ptr<Lepton>> leptons = leptonProcessor.GetLeptonsFromTree();
    for(auto lepton : leptons) event->AddLepton(lepton);
    
    if(dataType == xtracks::kSignal)          eventsSignal[(ESignal)setIter][year].push_back(event);
    else if(dataType == xtracks::kBackground) eventsBackground[(EBackground)setIter][year].push_back(event);
    else if(dataType == xtracks::kData)       eventsData[(EData)setIter][year].push_back(event);
    else                                      throw out_of_range("Unknown data type provided");
    
    if(iEvent >= 0) break;
  }
  if(inFileFriend) inFileFriend->Close();
}
