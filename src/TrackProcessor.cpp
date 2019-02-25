//
//  TrackProcessor.cpp
//
//  Created by Jeremi Niedziela on 30/01/2019.
//

#include "TrackProcessor.hpp"

TrackProcessor::TrackProcessor()
{
  
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
	
	for(int iLayer=0;iLayer<nLayers;iLayer++){
		track->layer[iLayer] = iLayer+1;
		track->dedx[iLayer] = 5.0;
	}
	
	double theta = track->GetTheta();
	
  double minR = layerR[nLayers-1];
  double maxR = layerR[nLayers];
  double decayR = RandDouble(minR, maxR);
	
  double decayX = decayR*cos(track->phi);
  double decayY = decayR*sin(track->phi);
  double decayZ = decayR/sin(theta)*cos(theta);
  
  track->decayPoint = make_unique<Point>(decayX, decayY, decayZ);
  
  return track;
}

bool TrackProcessor::IsPassingCut(const shared_ptr<Track> track,
                                  const unique_ptr<TrackCut> &cut)
{
  // check number of hits in pixel, stip and tracker in general
  if(cut->GetRequireSameNpixelHitsLayers()){
    if(track->nPixelHits != track->nPixelLayers) return false;
  }
  
  if(cut->GetRequireSameNtrackerHitsLayers()){
    if(track->nTrackerHits != track->nTrackerLayers) return false;
  }
  
  if(cut->GetRequireMcMatch()){
    if(track->mcMatch == 0) return false;
  }
  
  if(cut->GetNpixelHits().IsOutside(track->nPixelHits))  return false;
  if(cut->GetNpixelLayers().IsOutside(track->nPixelLayers)) return false;
  if(cut->GetNlayers().IsOutside(track->GetLastBarrelLayer())) return false;
  if(cut->GetNmissingInnerPixel().IsOutside(track->nMissingInnerPixelHits)) return false;
  if(cut->GetNmissingMiddleTracker().IsOutside(track->nMissingMiddleTrackerHits)) return false;
  if(cut->GetNmissingOuterTracker().IsOutside(track->nMissingOuterTrackerHits))  return false;
  
  // check number of dedx, number of detectors, number of clusters
  if(cut->GetNdedxClusters().IsOutside(track->nDedxClusters)) return false;
  if(cut->GetNdetIDs().IsOutside(track->nDetIDs))  return false;
  
  
  if(cut->GetTotalDedx().IsOutside(track->GetTotalDedx())) return false;
  
  for(int iCluster=0;iCluster<track->nDedxClusters;iCluster++){
    if(cut->GetDedxPerCluster().IsOutside(track->dedx[iCluster])) return false;
  }
  
  // check basic kinematical variables
  if(cut->GetPt().IsOutside(track->pt)) return false;
  if(cut->GetEta().IsOutside(track->eta)) return false;
  
  // check calo energy
  if(cut->GetCaloEmEnergy().IsOutside(track->caloEmEnergy))  return false;
  if(cut->GetCaloHadEnergy().IsOutside(track->caloHadEnergy))  return false;
  
  // check isolation
  if(cut->GetRelativeIsolation().IsOutside(track->relativeIsolation))  return false;
  
  // Check track-met ΔΦ
  if(cut->GetTrackMetDeltaPhi().GetMin() > -1000){
    
    // We use TLorentzVector to automatically deal with shifting the angle to [-π,π]
    TLorentzVector metVector, trackVector;
    metVector.SetPtEtaPhiM(track->eventMetPt, track->eventMetEta, track->eventMetPhi, track->eventMetMass);
    trackVector.SetPtEtaPhiM(track->pt, track->eta, track->phi, track->mass);
    
    if(cut->GetTrackMetDeltaPhi().IsOutside(metVector.DeltaPhi(trackVector))) return false;
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
      track->dedx[iLayer]     = arrayValuesFloat[Form("IsoTrack_dedxByLayer%i",iLayer)][iTrack];
      track->subDetId[iLayer] = arrayValuesInt[Form("IsoTrack_subDetIdByLayer%i",iLayer)][iTrack];
      track->sizeX[iLayer]    = arrayValuesInt[Form("IsoTrack_sizeXbyLayer%i",iLayer)][iTrack];
      track->sizeY[iLayer]    = arrayValuesInt[Form("IsoTrack_sizeYbyLayer%i",iLayer)][iTrack];
      track->detType[iLayer]  = arrayValuesInt[Form("IsoTrack_pixByLayer%i",iLayer)][iTrack];
      track->layer[iLayer]    = arrayValuesInt[Form("IsoTrack_layerOrSideByLayer%i",iLayer)][iTrack];
      track->ladder[iLayer]   = arrayValuesInt[Form("IsoTrack_ladderOrBladeByLayer%i",iLayer)][iTrack];
    }
    track->CalculateInternals();
		
		// Decay point is set only for random pion generation. It's not used in the fitter!!
		double minR = layerR[track->GetLastBarrelLayer()];
		double maxR = layerR[track->GetLastBarrelLayer()+1];
		double decayR = (maxR+minR)/2.;
		
		track->decayPoint = make_unique<Point>( decayR*cos(track->phi),
																						decayR*sin(track->phi),
																						decayR/sin(track->GetTheta())*cos(track->GetTheta()));
		
    tracks.push_back(track);
  }
  
  return tracks;
}


void TrackProcessor::SetupBranches(TTree *tree)
{
  // single int variables
  tree->SetBranchAddress("nIsoTrack",&nTracks);
  
  // single float variables
  // there are none at the moment
  
  // float array variables
  tree->SetBranchAddress("IsoTrack_eta",&arrayValuesFloat["IsoTrack_eta"]);
  tree->SetBranchAddress("IsoTrack_phi",&arrayValuesFloat["IsoTrack_phi"]);
  tree->SetBranchAddress("IsoTrack_caloEmEnergy",&arrayValuesFloat["IsoTrack_caloEmEnergy"]);
  tree->SetBranchAddress("IsoTrack_caloHadEnergy",&arrayValuesFloat["IsoTrack_caloHadEnergy"]);
  tree->SetBranchAddress("IsoTrack_dxy",&arrayValuesFloat["IsoTrack_dxy"]);
  tree->SetBranchAddress("IsoTrack_edxy",&arrayValuesFloat["IsoTrack_edxy"]);
  tree->SetBranchAddress("IsoTrack_dz",&arrayValuesFloat["IsoTrack_dz"]);
  tree->SetBranchAddress("IsoTrack_edz",&arrayValuesFloat["IsoTrack_edz"]);
  tree->SetBranchAddress("IsoTrack_mass",&arrayValuesFloat["IsoTrack_mass"]);
  tree->SetBranchAddress("IsoTrack_pt",&arrayValuesFloat["IsoTrack_pt"]);
  tree->SetBranchAddress("IsoTrack_relIso03",&arrayValuesFloat["IsoTrack_relIso03"]);
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    tree->SetBranchAddress(Form("IsoTrack_dedxByLayer%i",iLayer),
                           &arrayValuesFloat[Form("IsoTrack_dedxByLayer%i",iLayer)]);
  }
  
  // int array variables
  tree->SetBranchAddress("IsoTrack_pdgId",&arrayValuesInt["IsoTrack_pdgId"]);
  tree->SetBranchAddress("IsoTrack_charge",&arrayValuesInt["IsoTrack_charge"]);
  
  if(tree->GetBranchStatus("IsoTrack_mcMatch")){
    tree->SetBranchAddress("IsoTrack_mcMatch",&arrayValuesInt["IsoTrack_mcMatch"]);
  }
  else{
    cout<<"WARNING -- branch IsoTrack_mcMatch was not found. Will assume true for MC match"<<endl;
    for(int i=0; i<maxNtracks;i++){arrayValuesInt["IsoTrack_mcMatch"][i] = 1;}
  }
  tree->SetBranchAddress("IsoTrack_trackerLayers",&arrayValuesInt["IsoTrack_trackerLayers"]);
  tree->SetBranchAddress("IsoTrack_pixelLayers",&arrayValuesInt["IsoTrack_pixelLayers"]);
  tree->SetBranchAddress("IsoTrack_trackerHits",&arrayValuesInt["IsoTrack_trackerHits"]);
  tree->SetBranchAddress("IsoTrack_pixelHits",&arrayValuesInt["IsoTrack_pixelHits"]);
  tree->SetBranchAddress("IsoTrack_missingInnerPixelHits",&arrayValuesInt["IsoTrack_missingInnerPixelHits"]);
  tree->SetBranchAddress("IsoTrack_missingOuterPixelHits",&arrayValuesInt["IsoTrack_missingOuterPixelHits"]);
  tree->SetBranchAddress("IsoTrack_missingInnerStripHits",&arrayValuesInt["IsoTrack_missingInnerStripHits"]);
  tree->SetBranchAddress("IsoTrack_missingOuterStripHits",&arrayValuesInt["IsoTrack_missingOuterStripHits"]);
  tree->SetBranchAddress("IsoTrack_missingInnerTrackerHits",&arrayValuesInt["IsoTrack_missingInnerTrackerHits"]);
  tree->SetBranchAddress("IsoTrack_missingOuterTrackerHits",&arrayValuesInt["IsoTrack_missingOuterTrackerHits"]);
  tree->SetBranchAddress("IsoTrack_missingMiddleTrackerHits",&arrayValuesInt["IsoTrack_missingMiddleTrackerHits"]);
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    tree->SetBranchAddress(Form("IsoTrack_subDetIdByLayer%i",iLayer),
                           &arrayValuesInt[Form("IsoTrack_subDetIdByLayer%i",iLayer)]);
    tree->SetBranchAddress(Form("IsoTrack_sizeXbyLayer%i",iLayer),
                           &arrayValuesInt[Form("IsoTrack_sizeXbyLayer%i",iLayer)]);
    tree->SetBranchAddress(Form("IsoTrack_sizeYbyLayer%i",iLayer),
                           &arrayValuesInt[Form("IsoTrack_sizeYbyLayer%i",iLayer)]);
    tree->SetBranchAddress(Form("IsoTrack_pixByLayer%i",iLayer),
                           &arrayValuesInt[Form("IsoTrack_pixByLayer%i",iLayer)]);
    
    tree->SetBranchAddress(Form("IsoTrack_layerOrSideByLayer%i",iLayer),
                           &arrayValuesInt[Form("IsoTrack_layerOrSideByLayer%i",iLayer)]);
    tree->SetBranchAddress(Form("IsoTrack_ladderOrBladeByLayer%i",iLayer),
                           &arrayValuesInt[Form("IsoTrack_ladderOrBladeByLayer%i",iLayer)]);
    
  }
  
}
