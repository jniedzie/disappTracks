// -*- C++ -*-
//
// Package:    CharginoAnalysis/HitsExtractor
// Class:      HitsExtractor
//
/**\class HitsExtractor HitsExtractor.cc CharginoAnalysis/HitsExtractor/plugins/HitsExtractor.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Jeremi Niedziela
//         Created:  Thu, 29 Nov 2018 07:52:53 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <TH1D.h>
#include <TTree.h>

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

using reco::TrackCollection;

class HitsExtractor : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
public:
  explicit HitsExtractor(const edm::ParameterSet&);
  ~HitsExtractor();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  TTree *outputHits;
  
  std::vector<double> hitX;
  std::vector<double> hitY;
  std::vector<double> hitZ;
  std::vector<double> hitCharge;
  std::vector<double> hitSizeX;
  std::vector<double> hitSizeY;
  std::vector<double> hitXX;
  std::vector<double> hitYY;
  std::vector<double> hitZZ;
  
  std::vector<double> stripX;
  std::vector<double> stripY;
  std::vector<double> stripZ;
  std::vector<double> stripCharge;
  
  uint lumi;
  uint run;
  unsigned long long event;
  
  edm::EDGetTokenT<std::vector<reco::DeDxHitInfo>> dedxToken;
  edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelClusterToken;
  edm::EDGetTokenT<edmNew::DetSetVector<SiStripCluster>> stripClusterToken;
  edm::EDGetTokenT<std::vector<reco::Track>> generalTracksToken;
};

HitsExtractor::HitsExtractor(const edm::ParameterSet& iConfig) :
dedxToken(consumes<std::vector<reco::DeDxHitInfo>>(edm::InputTag("dedxHitInfo"))),
pixelClusterToken(consumes<edmNew::DetSetVector<SiPixelCluster>>(edm::InputTag("siPixelClusters"))),
stripClusterToken(consumes<edmNew::DetSetVector<SiStripCluster>>(edm::InputTag("siStripClusters"))),
generalTracksToken(consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks")))
{
  edm::Service<TFileService> fs;
  
  outputHits = fs->make<TTree>("hits","hits");
  outputHits->Branch("hitX", &hitX);
  outputHits->Branch("hitY", &hitY);
  outputHits->Branch("hitZ", &hitZ);
  outputHits->Branch("hitCharge", &hitCharge);
  outputHits->Branch("hitSizeX", &hitSizeX);
  outputHits->Branch("hitSizeY", &hitSizeY);
  outputHits->Branch("hitXX", &hitXX);
  outputHits->Branch("hitYY", &hitYY);
  outputHits->Branch("hitZZ", &hitZZ);
  
  outputHits->Branch("stripX", &stripX);
  outputHits->Branch("stripY", &stripY);
  outputHits->Branch("stripZ", &stripZ);
  outputHits->Branch("stripCharge", &stripCharge);
  
  outputHits->Branch("runNumber", &run);
  outputHits->Branch("lumiBlock", &lumi);
  outputHits->Branch("eventNumber", &event);
  
}

HitsExtractor::~HitsExtractor()
{
}

void HitsExtractor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ // method called for each event
  using namespace edm;
  
  bool filterTrackClusters = true;
  
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  event = iEvent.id().event();
  
//  uint searchRun = 297101;
//  if(searchRun != runNumber) return;
  
//  uint searchLumi = 113;
//  if(searchLumi != lumi) return;
  
//  uint searchEvent = 162801086;
//  if(searchEvent != eventNumber) return;
  
  ESHandle<TrackerGeometry> theTrackerGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  
  TrackingGeometry::DetContainer pixelBarrel = theTrackerGeometry->detsPXB();
  //  std::cout << "\n\nThere are " << pixelBarrel.size() << " detector elements in the PXB.\n\n" << std::endl;
  
  Handle<edmNew::DetSetVector<SiPixelCluster>> pixelClustersHandle;
  iEvent.getByToken(pixelClusterToken, pixelClustersHandle);
  
  Handle<edmNew::DetSetVector<SiStripCluster>> stripClustersHandle;
  iEvent.getByToken(stripClusterToken, stripClustersHandle);
  
  const edmNew::DetSetVector<SiPixelCluster> *pixelClusters = pixelClustersHandle.product();
  const edmNew::DetSetVector<SiStripCluster> *stripClusters = stripClustersHandle.product();
  
  
  // Store IDs of detectors which contain hits associated with tracks
  Handle<std::vector<reco::Track>> tracksHandle;
  iEvent.getByToken(generalTracksToken, tracksHandle);
  
  const std::vector<reco::Track> *tracks = tracksHandle.product();
  
  std::vector<uint> trackClustersIDs;
  
  std::cout<<"N tracks:"<<tracks->size()<<std::endl;
  
  for(uint iTrack=0;iTrack<tracks->size();iTrack++){
    reco::Track track = tracks->at(iTrack);
    
    for(uint iHit=0;iHit<track.recHitsSize();iHit++){
      const TrackingRecHit *hit = track.recHit(iHit).get();
      trackClustersIDs.push_back(hit->rawId());
    }
  }
  std::cout<<"N dets with tracks' hits:"<<trackClustersIDs.size()<<std::endl;
  
  hitX.clear();
  hitY.clear();
  hitZ.clear();
  hitCharge.clear();
  hitSizeX.clear();
  hitSizeY.clear();
  hitXX.clear();
  hitYY.clear();
  hitZZ.clear();
  stripX.clear();
  stripY.clear();
  stripZ.clear();
  stripCharge.clear();
  
  for(auto myDetSet=stripClusters->begin();myDetSet != stripClusters->end();myDetSet++){
    auto detSet = *myDetSet;
    
    // get detector geometry
    const GeomDet *geometry = nullptr;
    DetId detid(detSet.detId());
    TrackingGeometry::DetContainer stripGeom;
    
    if(detid.subdetId() == StripSubdetector::TIB)       stripGeom = theTrackerGeometry->detsTIB();
    else if(detid.subdetId() == StripSubdetector::TOB)  stripGeom = theTrackerGeometry->detsTOB();
    else                                                continue;
    
    for(auto i=stripGeom.begin();i != stripGeom.end(); i++){
      geometry = *i;
      if(geometry->geographicalId().rawId() == detSet.detId()) break;
    }
    
    if(!geometry){
      std::cout<<"Stips extraction -- geometry for detset not found..."<<std::endl;
      break;
    }
    
    // get detector topology
    auto detUnit = dynamic_cast<const StripGeomDetUnit*>(theTrackerGeometry->idToDet(detid));
    const StripTopology &stripTopology = detUnit->specificTopology();
    
//    std::cout<<"Adding "<<detSet.size()<<" strip clusters"<<std::endl;
    
    for(uint iCluster=0;iCluster<detSet.size();iCluster++){
      auto cluster = detSet[iCluster];
      
      LocalPoint l = stripTopology.localPosition(cluster.firstStrip());
      GlobalPoint p = detUnit->surface().toGlobal(l);

      stripX.push_back(p.x());
      stripY.push_back(p.y());
      stripZ.push_back(p.z());
      stripCharge.push_back(cluster.charge());
    }
  }
  std::cout<<"N pixel clusters:"<< pixelClusters->size()<<std::endl;
  int nRejected = 0;
  for(auto myDetSet=pixelClusters->begin();myDetSet != pixelClusters->end();myDetSet++){
    auto detSet = *myDetSet;
    
    // make sure that cluster belongs to te correct subdetector
    DetId detid(detSet.detId());
    
    if(detid.subdetId() != PixelSubdetector::PixelBarrel){
      continue;
    }
    
    // get detector geometry
    const GeomDet *geometry = nullptr;
    std::cout<<"Number of subdets in pixel barrel:"<<pixelBarrel.size()<<std::endl;
    std::cout<<"Hit det id:"<<detid.rawId()<<std::endl;
    
    for(auto i=pixelBarrel.begin();i != pixelBarrel.end(); i++){
      const GeomDet *det = *i;
      PXBDetId id = det->geographicalId();
      
      if(i==pixelBarrel.begin())  std::cout<<"\t\tFirst geom det id:"<<id.rawId()<<std::endl;
      
      if(id.rawId() == detid.rawId()){
        std::cout<<"Hit's det id matches with one of the subdet ids"<<std::endl;
        geometry = det;
        break;
      }
    }
    
    if(!geometry){
      std::cout<<"Pixels extraction -- geometry for detset not found..."<<std::endl;
      continue;
//      break;
    }
    
    if(filterTrackClusters){
      uint thisDetID = geometry->geographicalId().rawId();
      if(std::find(trackClustersIDs.begin(), trackClustersIDs.end(), thisDetID) != trackClustersIDs.end()) {
        nRejected++;
        continue;
      }
    }
    
    // get detector topology
    auto detUnit = dynamic_cast<const PixelGeomDetUnit*>(theTrackerGeometry->idToDet(detid));
    const PixelTopology &pixelTopology = detUnit->specificTopology();
  
    std::cout<<"Adding "<<detSet.size()<<" pixel clusters"<<std::endl;
    
    for(uint iCluster=0;iCluster<detSet.size();iCluster++){
      SiPixelCluster cluster = detSet[iCluster];
      
      
      
      LocalPoint l = pixelTopology.localPosition(MeasurementPoint(cluster.x(),cluster.y()));
      GlobalPoint p = detUnit->surface().toGlobal(l);
      LocalPoint ll = pixelTopology.localPosition(MeasurementPoint(cluster.minPixelRow(),
                                                                   cluster.minPixelCol()
                                                                   ));
      GlobalPoint pp = detUnit->surface().toGlobal(ll);

      hitX.push_back(p.x());
      hitY.push_back(p.y());
      hitZ.push_back(p.z());
      hitCharge.push_back(cluster.charge());
      hitSizeX.push_back(cluster.sizeX());
      hitSizeY.push_back(cluster.sizeY());
      hitXX.push_back(pp.x());
      hitYY.push_back(pp.y());
      hitZZ.push_back(pp.z());
    }
  }
  std::cout<<"N hits associated with tracks:"<<nRejected<<std::endl;
  outputHits->Fill();
}

void HitsExtractor::beginJob()
{
// method called once each job just before starting event loop
}

void HitsExtractor::endJob()
{
// method called once each job just after ending the event loop
}

void HitsExtractor::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("hits",edm::InputTag("dedxHitInfo"));
  desc.addUntracked<edm::InputTag>("pixelClusters",edm::InputTag("siPixelClusters"));
  desc.addUntracked<edm::InputTag>("stripClusters",edm::InputTag("siStripClusters"));
  desc.addUntracked<edm::InputTag>("generalTracks",edm::InputTag("generalTracks"));
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HitsExtractor);
