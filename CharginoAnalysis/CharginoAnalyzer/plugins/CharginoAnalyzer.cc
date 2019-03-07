// -*- C++ -*-
//
// Package:    CharginoAnalysis/CharginoAnalyzer
// Class:      CharginoAnalyzer
// 
/**\class CharginoAnalyzer CharginoAnalyzer.cc CharginoAnalysis/CharginoAnalyzer/plugins/CharginoAnalyzer.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Jeremi Niedziela
//         Created:  Tue, 19 Feb 2019 12:48:21 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <iomanip>

#include <TLorentzVector.h>
#include <TH1D.h>
#include <TTree.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include <DataFormats/TrackReco/interface/Track.h>

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetEnumerators.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using namespace reco;
using namespace std;
using namespace edm;

class CharginoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit CharginoAnalyzer(const edm::ParameterSet&);
  ~CharginoAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  EDGetTokenT<vector<GenParticle>> genParticlesToken;
  EDGetTokenT<vector<SimTrack>> simTracksToken;
  EDGetTokenT<vector<PSimHit>> simHitsPixelHighToken;
  EDGetTokenT<vector<PSimHit>> simHitsPixelLowToken;
  EDGetTokenT<OwnVector<TrackingRecHit,ClonePolicy<TrackingRecHit>>> trackingRecHitsToken;
  EDGetTokenT<vector<Track>> recTracksToken;
  EDGetTokenT<vector<Track>> trackingParticlesToken;
  
  Handle<vector<GenParticle>> genParticles;
  Handle<vector<SimTrack>> simTracks;
  Handle<vector<PSimHit>> simHitsPixelHigh;
  Handle<vector<PSimHit>> simHitsPixelLow;
  Handle<OwnVector<TrackingRecHit,ClonePolicy<TrackingRecHit>>> trackingRecHitsRef;
  Handle<vector<Track>> recTracks;
  
  ESHandle<TrackerGeometry> theTrackerGeometry;
  
  TH1D *pionPt;
  TH1D *pionPz;
  
  TTree *outputTree;
  
  std::vector<double> pion_vx;
  std::vector<double> pion_vy;
  std::vector<double> pion_vz;
  std::vector<double> pion_px;
  std::vector<double> pion_py;
  std::vector<double> pion_pz;
  std::vector<int> pion_charge;
  
  uint lumi;
  uint run;
  unsigned long long event;
  
  void printDeeply(const Candidate *p);
};

CharginoAnalyzer::CharginoAnalyzer(const edm::ParameterSet& iConfig) :
genParticlesToken(consumes<vector<GenParticle>>(InputTag("genParticles"))),
simTracksToken(consumes<vector<SimTrack>>(InputTag("g4SimHits"))),
simHitsPixelHighToken(consumes<vector<PSimHit>>(InputTag("g4SimHits", "TrackerHitsPixelBarrelHighTof"))),
simHitsPixelLowToken(consumes<vector<PSimHit>>(InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"))),
trackingRecHitsToken(consumes<OwnVector<TrackingRecHit,ClonePolicy<TrackingRecHit>>>(InputTag("generalTracks"))),
recTracksToken(consumes<vector<Track>>(InputTag("generalTracks")))
{
  edm::Service<TFileService> fs;
  
  pionPt = fs->make<TH1D>("pionPt","pionPtt",100,0.,1.);
  pionPz = fs->make<TH1D>("pionPz","pionPzz",100,0.,1.);
  
  outputTree = fs->make<TTree>("pions","pions");
  outputTree->Branch("pion_vx", &pion_vx);
  outputTree->Branch("pion_vy", &pion_vy);
  outputTree->Branch("pion_vz", &pion_vz);
  outputTree->Branch("pion_px", &pion_px);
  outputTree->Branch("pion_py", &pion_py);
  outputTree->Branch("pion_pz", &pion_pz);
  outputTree->Branch("pion_charge", &pion_charge);
  
  outputTree->Branch("runNumber", &run);
  outputTree->Branch("lumiBlock", &lumi);
  outputTree->Branch("eventNumber", &event);
}


CharginoAnalyzer::~CharginoAnalyzer()
{
}

void CharginoAnalyzer::printDeeply(const Candidate *p)
{
  auto mom = p->momentum();
  auto vertex = p->vertex();
  
  const std::vector<TrackingRecHit*> &trackingRecHits = trackingRecHitsRef.product()->data();
  
  string pdg = "unknown";
  if(p->pdgId() ==  1000022)       pdg = "Chi0";
  else if(p->pdgId() ==  1000024)  pdg = "Chi+";
  else if(p->pdgId() == -1000024)  pdg = "Chi-";
  else if(p->pdgId() ==  211)      pdg = "pi+";
  else if(p->pdgId() == -211)      pdg = "pi-";
  else                             pdg = to_string(p->pdgId());
  
  cout<<pdg<<"\t("<<p->status()<<")"<<endl;
  cout<<setprecision(4);
  cout<<"\tp = ("<<mom.x()<<","<<mom.y()<<","<<mom.z()<<")\tpt:"<<p->pt()<<endl;
  cout<<"\tv = ("<<vertex.x()<<","<<vertex.y()<<","<<vertex.z()<<")\tXY:"<<sqrt(p->vx()*p->vx()+p->vy()*p->vy())<<endl;
  cout<<"\teta:"<<mom.eta()<<"\tphi:"<<mom.phi()<<endl;
  
  if(abs(p->pdgId()) ==  211){
    pion_vx.push_back(vertex.x());
    pion_vy.push_back(vertex.y());
    pion_vz.push_back(vertex.z());
    pion_px.push_back(mom.x());
    pion_py.push_back(mom.y());
    pion_pz.push_back(mom.z());
    pion_charge.push_back(p->pdgId() < 0 ? -1 : 1);
  }
  
  // Print rec tracks and rec hits
  for(uint iTrack=0; iTrack<recTracks->size(); iTrack++){
    const Track &track = (*recTracks)[iTrack];
  
    if(fabs(track.eta()-mom.eta()) < 0.1 &&
       fabs(track.phi()-mom.phi()) < 0.1){
       
      cout<<"Rec track ("<<track.recHitsSize()<<" hits):"<<endl;
      cout<<"eta:"<<track.eta()<<"\tphi:"<<track.phi()<<endl;
      
      for(size_t iHit=0;iHit<track.recHitsSize();iHit++){
        const TrackingRecHit *recHit = track.recHit(iHit).get();
        if(!recHit){
          cout<<"No valid hit..."<<endl;
          continue;
        }
//        GlobalPoint recGlobal = recHit->globalPosition();
//        
//        cout<<"Rec hit:";
//        cout<<recGlobal.x()<<", ";
//        cout<<recGlobal.y()<<", ";
//        cout<<recGlobal.z()<<endl;
        
      }
      
    }
    
  }
  
  // Print sim tracks and sim hits
  for(uint iTrack=0; iTrack<simTracks->size(); iTrack++){
    const SimTrack &track = (*simTracks)[iTrack];
    
    auto trackMom = track.momentum();
    
    if(fabs(trackMom.eta()-mom.eta()) < 0.1 &&
       fabs(trackMom.phi()-mom.phi()) < 0.1 &&
       fabs(trackMom.x()-mom.x()) < 0.1 &&
       fabs(trackMom.y()-mom.y()) < 0.1 &&
       fabs(trackMom.z()-mom.z()) < 0.1){
    
      auto mom = track.momentum();
      
      cout<<"\tSim track:"<<endl;
      cout<<"\t\teta:"<<mom.eta()<<"\tphi:"<<mom.phi()<<endl;
      cout<<"\t\tmomentum:"<<mom.x()<<","<<mom.y()<<","<<mom.z()<<endl;
      cout<<"\t\tLow hits:"<<endl;
      
      for(uint iHit=0;iHit<simHitsPixelLow->size();iHit++){
        const PSimHit &hit = (*simHitsPixelLow)[iHit];
        
        if(hit.trackId() == track.trackId()){
        
          
          
          auto detUnit = dynamic_cast<const PixelGeomDetUnit*>(theTrackerGeometry->idToDet(hit.detUnitId()));
          const PixelTopology &pixelTopology = detUnit->specificTopology();
          
          LocalPoint l = pixelTopology.localPosition(MeasurementPoint(hit.localPosition().x(),
                                                                      hit.localPosition().y()));
          
          LocalPoint entryPointLocal = pixelTopology.localPosition(MeasurementPoint(hit.entryPoint().x(),
                                                                                    hit.entryPoint().y()));
          
          LocalPoint exitPointLocal = pixelTopology.localPosition(MeasurementPoint(hit.exitPoint().x(),
                                                                                   hit.exitPoint().y()));
          
          GlobalPoint p = detUnit->surface().toGlobal(l);
          GlobalPoint entryPointGlobal = detUnit->surface().toGlobal(entryPointLocal);
          GlobalPoint exitPointGlobal = detUnit->surface().toGlobal(exitPointLocal);
          
          cout<<"Center:";
          cout<<p.x()<<", ";
          cout<<p.y()<<", ";
          cout<<p.z()<<endl;
          
          cout<<"entry: ";
          cout<<entryPointGlobal.x()<<", ";
          cout<<entryPointGlobal.y()<<", ";
          cout<<entryPointGlobal.z()<<endl;
          
          cout<<"exit: ";
          cout<<exitPointGlobal.x()<<", ";
          cout<<exitPointGlobal.y()<<", ";
          cout<<exitPointGlobal.z()<<endl;
          
          for(uint iRecHit=0;iRecHit<trackingRecHits.size();iRecHit++){
            TrackingRecHit *recHit = trackingRecHits[iRecHit];
            
            if(hit.detUnitId() == recHit->geographicalId().rawId()){
            
              if(!recHit->isValid()){
                cout<<"Has rec hits, but not valid"<<endl;
                continue;
              }
              
              
              /*
              LocalPoint recLocal = pixelTopology.localPosition(MeasurementPoint(recHit->localPosition().x(),
                                                                                 recHit->localPosition().y()));
              
              GlobalPoint recGlobal = recHit->globalPosition();
              
              cout<<"Rec hit:";
              cout<<recGlobal.x()<<", ";
              cout<<recGlobal.y()<<", ";
              cout<<recGlobal.z()<<endl;
               */
            }
            
          }
          
        }
      }
      cout<<"\t\tHigh hits:"<<endl;
      
      for(uint iHit=0;iHit<simHitsPixelHigh->size();iHit++){
        const PSimHit &hit = (*simHitsPixelHigh)[iHit];
        if(hit.trackId() == track.trackId()){
          
          auto detUnit = dynamic_cast<const PixelGeomDetUnit*>(theTrackerGeometry->idToDet(hit.detUnitId()));
          const PixelTopology &pixelTopology = detUnit->specificTopology();
          
          LocalPoint l = pixelTopology.localPosition(MeasurementPoint(hit.localPosition().x(),
                                                                      hit.localPosition().y()));
          GlobalPoint p = detUnit->surface().toGlobal(l);
          
          cout<<p.x()<<", ";
          cout<<p.y()<<", ";
          cout<<p.z()<<endl;
          
        }
      }
      
      break;
    }
  }
  
  if(p->numberOfDaughters() != 0){
    cout<<"\n|"<<endl;
    cout<<"V\n"<<endl;
  }
  
  for(uint iDaughter=0; iDaughter<p->numberOfDaughters(); iDaughter++) {
    printDeeply(p->daughter(iDaughter));
  }
  
  
}

// ------------ method called for each event  ------------
void CharginoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool verbose = true;
  
  pion_vx.clear();
  pion_vy.clear();
  pion_vz.clear();
  pion_px.clear();
  pion_py.clear();
  pion_pz.clear();
  pion_charge.clear();
  
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  event = iEvent.id().event();
  
  if(verbose){
    cout<<"\n\n================================================================="<<endl;
    cout<<"Event:"<<event<<"\trun:"<<run<<"\tlumi:"<<lumi<<endl;
    cout<<"\n"<<endl;
  }
  
  iEvent.getByToken(genParticlesToken, genParticles);
  iEvent.getByToken(simTracksToken, simTracks);
  iEvent.getByToken(simHitsPixelHighToken, simHitsPixelHigh);
  iEvent.getByToken(simHitsPixelLowToken, simHitsPixelLow);
  iEvent.getByToken(trackingRecHitsToken, trackingRecHitsRef);
  iEvent.getByToken(recTracksToken, recTracks);
  
  iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  
  // pion id:       211
  // chargino id:   1000024
  // neutralino id: 1000022

  for(size_t i=0; i<genParticles->size(); i++){
    const GenParticle *p = &(*genParticles)[i];
    const Candidate *mom = p->mother();
    if(!mom) continue;
    
    int id = p->pdgId();
    int momId = mom->pdgId();
    
    // check if this is the first chargino/neutralino in the chain
    if(abs(id)!=1000022 && abs(id)!=1000024) continue;
    if(abs(momId) == 1000022 || abs(momId) == 1000024) continue;
      
    if(verbose){
      cout<<"\n\nParticle chain:\n\n"<<endl;
      printDeeply(p);
    }
    
  }
  
  outputTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
CharginoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CharginoAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CharginoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CharginoAnalyzer);
