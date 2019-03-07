import FWCore.ParameterSet.Config as cms

hitsExtractor = cms.EDAnalyzer('HitsExtractor',
                               hits = cms.untracked.InputTag('dedxHitInfo'),
                               pixelClusters = cms.untracked.InputTag('siPixelClusters'),
                               stripClusters = cms.untracked.InputTag('siStripClusters'),
                               generalTracks = cms.untracked.InputTag('generalTracks'),
)
