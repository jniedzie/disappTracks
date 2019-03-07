import FWCore.ParameterSet.Config as cms

process = cms.Process("CharginoAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch///eos/cms/store/group/phys_diffraction/jniedzie/susy/chargino300GeV_ctau10cm_GEN-SIM-RAW-RECO.root'
    ),
    secondaryFileNames=cms.untracked.vstring(
        'root://eoscms.cern.ch///eos/cms/store/group/phys_diffraction/jniedzie/susy/chargino300GeV_ctau10cm_GEN-SIM-RAW.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("/afs/cern.ch/work/j/jniedzie/private/charginoAnalysisPions.root"))

process.load("CharginoAnalysis.CharginoAnalyzer.CfiFile_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v1', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14', '')

process.p = cms.Path(process.CharginoAnalyzer)
