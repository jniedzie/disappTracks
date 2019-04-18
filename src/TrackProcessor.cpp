//  TrackProcessor.cpp
//
//  Created by Jeremi Niedziela on 30/01/2019.

#include "TrackProcessor.hpp"

TrackProcessor trackProcessor = TrackProcessor();

TrackProcessor::TrackProcessor()
{
  arrayNamesFloat = {
    "IsoTrack_eta",
    "IsoTrack_phi",
    "IsoTrack_caloEmEnergy",
    "IsoTrack_caloHadEnergy",
    "IsoTrack_dxy",
    "IsoTrack_edxy",
    "IsoTrack_dz",
    "IsoTrack_edz",
    "IsoTrack_mass",
    "IsoTrack_pt",
    "IsoTrack_relIso03"
  };
  
  arrayNamesInt = {
    "IsoTrack_pdgId",
    "IsoTrack_charge",
    "IsoTrack_trackerLayers",
    "IsoTrack_pixelLayers",
    "IsoTrack_trackerHits",
    "IsoTrack_pixelHits",
    "IsoTrack_missingInnerPixelHits",
    "IsoTrack_missingOuterPixelHits",
    "IsoTrack_missingInnerStripHits",
    "IsoTrack_missingOuterStripHits",
    "IsoTrack_missingInnerTrackerHits",
    "IsoTrack_missingOuterTrackerHits",
    "IsoTrack_missingMiddleTrackerHits",
    "IsoTrack_mcMatch"
  };
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    arrayNamesFloat.push_back(Form("IsoTrack_dedxByHit%i",iLayer));
    
    arrayNamesInt.push_back(Form("IsoTrack_subDetIdByHit%i",iLayer));
    arrayNamesInt.push_back(Form("IsoTrack_sizeXbyHit%i",iLayer));
    arrayNamesInt.push_back(Form("IsoTrack_sizeYbyHit%i",iLayer));
    arrayNamesInt.push_back(Form("IsoTrack_pixByHit%i",iLayer));
    arrayNamesInt.push_back(Form("IsoTrack_layerOrSideByHit%i",iLayer));
    arrayNamesInt.push_back(Form("IsoTrack_ladderOrBladeByHit%i",iLayer));
  }
}

TrackProcessor::~TrackProcessor()
{
  
}

shared_ptr<Track> TrackProcessor::GetRandomTrack(int nLayers, double maxEta)
{
  auto track = make_shared<Track>();
  
  track->eta = RandDouble(-maxEta, maxEta);
  track->phi = RandDouble(0, 2*TMath::Pi());
  track->nTrackerLayers = track->nPixelLayers = nLayers;
	
	for(int iHit=0;iHit<nLayers;iHit++){
		track->layer[iHit] = iHit+1;
		track->dedx[iHit] = 5.0;
	}
	
	double theta = track->GetTheta();
	
  double minR = layerR[nLayers-1];
  double maxR = layerR[nLayers];
  double decayR = RandDouble(minR, maxR);
	
  // This assumes vertex at 0,0,0!!
  double decayX = decayR*cos(track->phi);
  double decayY = decayR*sin(track->phi);
  double decayZ = decayR/sin(theta)*cos(theta);
  
  track->decayPoint = Point(decayX, decayY, decayZ);
  
  return track;
}

bool TrackProcessor::IsPassingCut(const shared_ptr<Track> track,
                                  const TrackCut &cut)
{
  if(cut.requiresMcMatch && track->mcMatch == 0) return false;
  
  
  if(cut.nPixelHits.IsOutside(track->nPixelHits))  return false;
  if(cut.nPixelLayers.IsOutside(track->nPixelLayers)) return false;
  if(cut.nLayers.IsOutside(track->GetNtrackerLayers())) return false;
  if(cut.nMissingInnerPixel.IsOutside(track->nMissingInnerPixelHits)) return false;
  if(cut.nMissingMiddleTracker.IsOutside(track->nMissingMiddleTrackerHits)) return false;
  if(cut.nMissingOuterTracker.IsOutside(track->nMissingOuterTrackerHits))  return false;
  
  // check number of dedx, number of detectors, number of clusters
  if(cut.nDedxClusters.IsOutside(track->nDedxClusters)) return false;
  if(cut.nDetIDs.IsOutside(track->nDetIDs))  return false;
  
  
  if(cut.totalDeDx.IsOutside(track->GetTotalDedx())) return false;
  
  for(int iCluster=0;iCluster<track->nDedxClusters;iCluster++){
    if(cut.dedxPerCluster.IsOutside(track->dedx[iCluster])) return false;
  }
  
  // check basic kinematical variables
  if(cut.pt.IsOutside(track->pt)) return false;
  if(cut.eta.IsOutside(track->eta)) return false;
  
  // check calo energy
  if(cut.caloEmEnergy.IsOutside(track->caloEmEnergy))  return false;
  if(cut.caloHadEnergy.IsOutside(track->caloHadEnergy))  return false;
  
  // check isolation
  if(cut.relativeIsolation.IsOutside(track->relativeIsolation))  return false;
  
  // Check track-met ΔΦ
  if(cut.trackMetDeltaPhi.GetMin() > -inf){
    TLorentzVector metVector, trackVector;
    metVector.SetPtEtaPhiM(track->eventMetPt, track->eventMetEta, track->eventMetPhi, track->eventMetMass);
    trackVector.SetPtEtaPhiM(track->pt, track->eta, track->phi, track->mass);
    
    if(cut.trackMetDeltaPhi.IsOutside(metVector.DeltaPhi(trackVector))) return false;
  }
  
  return true;
}

vector<shared_ptr<Track>> TrackProcessor::GetTracksFromTree()
{
  vector<shared_ptr<Track>> tracks = vector<shared_ptr<Track>>();
  
  for(int iTrack=0;iTrack<nTracks;iTrack++){
    auto track = make_shared<Track>();
    
    // float variables
    track->eta               = arrayValuesFloat["IsoTrack_eta"][iTrack];
    track->phi               = arrayValuesFloat["IsoTrack_phi"][iTrack];
    track->caloEmEnergy      = arrayValuesFloat["IsoTrack_caloEmEnergy"][iTrack];
    track->caloHadEnergy     = arrayValuesFloat["IsoTrack_caloHadEnergy"][iTrack];
    track->dxy               = arrayValuesFloat["IsoTrack_dxy"][iTrack];
    track->dxyErr            = arrayValuesFloat["IsoTrack_edxy"][iTrack];
    track->dz                = arrayValuesFloat["IsoTrack_dz"][iTrack];
    track->dzErr             = arrayValuesFloat["IsoTrack_edz"][iTrack];
    track->mass              = arrayValuesFloat["IsoTrack_mass"][iTrack];
    track->pt                = arrayValuesFloat["IsoTrack_pt"][iTrack];
    track->relativeIsolation = arrayValuesFloat["IsoTrack_relIso03"][iTrack];

    // int variables
    track->pid                       = arrayValuesInt["IsoTrack_pdgId"][iTrack];
    track->charge                    = arrayValuesInt["IsoTrack_charge"][iTrack];
    track->mcMatch                   = arrayValuesInt["IsoTrack_mcMatch"][iTrack];
    track->nTrackerLayers            = arrayValuesInt["IsoTrack_trackerLayers"][iTrack];
    track->nPixelLayers              = arrayValuesInt["IsoTrack_pixelLayers"][iTrack];
    track->nTrackerHits              = arrayValuesInt["IsoTrack_trackerHits"][iTrack];
    track->nPixelHits                = arrayValuesInt["IsoTrack_pixelHits"][iTrack];
    track->nMissingInnerPixelHits    = arrayValuesInt["IsoTrack_missingInnerPixelHits"][iTrack];
    track->nMissingOuterPixelHits    = arrayValuesInt["IsoTrack_missingOuterPixelHits"][iTrack];
    track->nMissingInnerStripHits    = arrayValuesInt["IsoTrack_missingInnerStripHits"][iTrack];
    track->nMissingOuterStripHits    = arrayValuesInt["IsoTrack_missingOuterStripHits"][iTrack];
    track->nMissingInnerTrackerHits  = arrayValuesInt["IsoTrack_missingInnerTrackerHits"][iTrack];
    track->nMissingOuterTrackerHits  = arrayValuesInt["IsoTrack_missingOuterTrackerHits"][iTrack];
    track->nMissingMiddleTrackerHits = arrayValuesInt["IsoTrack_missingMiddleTrackerHits"][iTrack];

    for(int iLayer=0;iLayer<nLayers;iLayer++){
      track->dedx[iLayer]     = arrayValuesFloat[Form("IsoTrack_dedxByHit%i",iLayer)][iTrack];
      track->subDetId[iLayer] = arrayValuesInt[Form("IsoTrack_subDetIdByHit%i",iLayer)][iTrack];
      track->sizeX[iLayer]    = arrayValuesInt[Form("IsoTrack_sizeXbyHit%i",iLayer)][iTrack];
      track->sizeY[iLayer]    = arrayValuesInt[Form("IsoTrack_sizeYbyHit%i",iLayer)][iTrack];
      track->detType[iLayer]  = arrayValuesInt[Form("IsoTrack_pixByHit%i",iLayer)][iTrack];
      track->layer[iLayer]    = arrayValuesInt[Form("IsoTrack_layerOrSideByHit%i",iLayer)][iTrack];
      track->ladder[iLayer]   = arrayValuesInt[Form("IsoTrack_ladderOrBladeByHit%i",iLayer)][iTrack];
    }
    track->CalculateInternals();
    tracks.push_back(track);
  }
  
  return tracks;
}

void TrackProcessor::SaveTracksToTree(vector<shared_ptr<Track>> tracks)
{
	nTracks = (int)tracks.size();
	
	for(int iTrack=0;iTrack<nTracks;iTrack++)
	{
		// float array variables
		arrayValuesFloat["IsoTrack_eta"][iTrack]           = tracks[iTrack]->GetEta();
		arrayValuesFloat["IsoTrack_phi"][iTrack]           = tracks[iTrack]->GetPhi();
		arrayValuesFloat["IsoTrack_caloEmEnergy"][iTrack]  = tracks[iTrack]->GetCaloEmEnergy();
		arrayValuesFloat["IsoTrack_caloHadEnergy"][iTrack] = tracks[iTrack]->GetCaloHadEnergy();
		arrayValuesFloat["IsoTrack_dxy"][iTrack]           = tracks[iTrack]->GetDxy();
		arrayValuesFloat["IsoTrack_edxy"][iTrack]          = tracks[iTrack]->GetDxyErr();
		arrayValuesFloat["IsoTrack_dz"][iTrack]            = tracks[iTrack]->GetDz();
		arrayValuesFloat["IsoTrack_edz"][iTrack]           = tracks[iTrack]->GetDzErr();
		arrayValuesFloat["IsoTrack_mass"][iTrack]          = tracks[iTrack]->GetMass();
		arrayValuesFloat["IsoTrack_pt"][iTrack]            = tracks[iTrack]->GetPt();
		arrayValuesFloat["IsoTrack_relIso03"][iTrack]      = tracks[iTrack]->GetRelativeIsolation();
		
		// int array variables
		arrayValuesInt["IsoTrack_pdgId"][iTrack]                    = tracks[iTrack]->GetPid();
		arrayValuesInt["IsoTrack_charge"][iTrack]                   = tracks[iTrack]->GetCharge();
		arrayValuesInt["IsoTrack_mcMatch"][iTrack]                  = tracks[iTrack]->GetMcMatch();
		arrayValuesInt["IsoTrack_trackerLayers"][iTrack]            = tracks[iTrack]->GetNtrackerLayers();
		arrayValuesInt["IsoTrack_pixelLayers"][iTrack]              = tracks[iTrack]->GetNpixelLayers();
		arrayValuesInt["IsoTrack_trackerHits"][iTrack]              = tracks[iTrack]->GetNtrackerHits();
		arrayValuesInt["IsoTrack_pixelHits"][iTrack]                = tracks[iTrack]->GetNpixelHits();
		arrayValuesInt["IsoTrack_missingInnerPixelHits"][iTrack]    = tracks[iTrack]->GetNmissingInnerPixelHits();
		arrayValuesInt["IsoTrack_missingOuterPixelHits"][iTrack]    = tracks[iTrack]->GetNmissingOuterPixelHits();
		arrayValuesInt["IsoTrack_missingInnerStripHits"][iTrack]    = tracks[iTrack]->GetNmissingInnerStripHits();
		arrayValuesInt["IsoTrack_missingOuterStripHits"][iTrack]    = tracks[iTrack]->GetNmissingOuterStripHits();
		arrayValuesInt["IsoTrack_missingInnerTrackerHits"][iTrack]  = tracks[iTrack]->GetNmissingInnerTrackerHits();
		arrayValuesInt["IsoTrack_missingOuterTrackerHits"][iTrack]  = tracks[iTrack]->GetNmissingOuterTrackerHits();
		arrayValuesInt["IsoTrack_missingMiddleTrackerHits"][iTrack] = tracks[iTrack]->GetNmissingMiddleTrackerHits();
		
		for(int iHit=0;iHit<nLayers;iHit++){
			arrayValuesFloat[Form("IsoTrack_dedxByHit%i",iHit)][iTrack] 			 = tracks[iTrack]->GetDeDxForHit(iHit);
			arrayValuesInt[Form("IsoTrack_subDetIdByHit%i",iHit)][iTrack]      = tracks[iTrack]->GetSubDetIdForHit(iHit);
			arrayValuesInt[Form("IsoTrack_sizeXbyHit%i",iHit)][iTrack]         = tracks[iTrack]->GetSizeXforHit(iHit);
			arrayValuesInt[Form("IsoTrack_sizeYbyHit%i",iHit)][iTrack]         = tracks[iTrack]->GetSizeYforHit(iHit);
			arrayValuesInt[Form("IsoTrack_pixByHit%i",iHit)][iTrack]           = tracks[iTrack]->GetDetTypeForHit(iHit);
			arrayValuesInt[Form("IsoTrack_layerOrSideByHit%i",iHit)][iTrack]   = tracks[iTrack]->GetLayerForHit(iHit);
			arrayValuesInt[Form("IsoTrack_ladderOrBladeByHit%i",iHit)][iTrack] = tracks[iTrack]->GetLadderForHit(iHit);
		}
	}
}

void TrackProcessor::SetupBranchesForReading(TTree *tree)
{
  tree->SetBranchAddress("nIsoTrack",&nTracks);
  
  for(string name : arrayNamesFloat){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesFloat[name]);
  }
  
  for(string name : arrayNamesInt){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    
    // special check for mcMatch branch, which may not exist
    if(name=="IsoTrack_mcMatch" && !tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- branch IsoTrack_mcMatch was not found. Will assume true for MC match"<<endl;
      for(int i=0; i<maxNtracks;i++){arrayValuesInt[name][i] = 1;}
      continue;
    }
    
    tree->SetBranchAddress(name.c_str(), &arrayValuesInt[name]);
  }
}

void TrackProcessor::SetupBranchesForWriting(TTree *tree)
{
	tree->Branch("nIsoTrack", &nTracks, "nIsoTrack/I");

  for(string name : arrayNamesFloat){
    tree->Branch(name.c_str(), &arrayValuesFloat[name], Form("%s[nIsoTrack]/F", name.c_str()));
  }
  
  for(string name : arrayNamesInt){
      tree->Branch(name.c_str(), &arrayValuesInt[name], Form("%s[nIsoTrack]/I", name.c_str()));
  }
}


