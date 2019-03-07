import FWCore.ParameterSet.Config as cms

process = cms.Process("HitsExtractor")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
fileNames = cms.untracked.vstring(
#        'root://eoscms.cern.ch///eos/cms/store/data/Run2016B/SingleMuon/MINIAOD/03Feb2017_ver2-v2/100000/AA09E663-A0EC-E611-9573-0025904C66EC.root'
#        'root://eoscms.cern.ch///eos/cms/store/data/Run2016B/SingleMuon/RAW-RECO/ZMu-PromptReco-v2/000/274/198/00000/501DC2B8-7B27-E611-90E7-02163E011EA7.root'
                                  
#'root://eoscms.cern.ch///eos/cms/store/data/Run2017B/MET/RAW-RECO/HighMET-17Nov2017-v1/70000/8A005F21-0EE4-E711-817D-0CC47A5FBE31.root'
#'root://eoscms.cern.ch///eos/cms/store/group/phys_diffraction/jniedzie/pickevents.root'
																		'root://eoscms.cern.ch///eos/cms/store/group/phys_diffraction/jniedzie/susy/chargino300GeV_ctau10cm_GEN-SIM-RAW-RECO.root'
#'root://xrootd-cms.infn.it///store/data/Run2017B/MET/RAW-RECO/HighMET-17Nov2017-v1/70000/8A005F21-0EE4-E711-817D-0CC47A5FBE31.root'
                                  
#'root://eoscms.cern.ch///eos/cms/store/group/phys_diffraction/jniedzie/8A005F21-0EE4-E711-817D-0CC47A5FBE31.root'
                                  )
                            )

process.TFileService = cms.Service("TFileService",fileName = cms.string("/afs/cern.ch/work/j/jniedzie/private/charginoAllHits.root"))

process.load("CharginoAnalysis.HitsExtractor.HitsExtractor_cfi")
process.hitsExtractor.hits="dedxHitInfo"
process.hitsExtractor.pixelClusters="siPixelClusters"
process.hitsExtractor.stripClusters="siStripClusters"
process.hitsExtractor.generalTracks="generalTracks"

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14', '')



process.p = cms.Path(process.hitsExtractor)
