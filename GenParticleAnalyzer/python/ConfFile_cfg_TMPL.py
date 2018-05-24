import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
XXXFILES
    )


)

process.demo = cms.EDAnalyzer('GenParticleAnalyzer'
)


process.TFileService = cms.Service("TFileService",
      fileName = cms.string("XXXOUTFILE"),
      closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.demo)
