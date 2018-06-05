import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.source = cms.Source("PoolSource",
    ## replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
    #    'file:myfile.root'
    #)
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2017E/DoubleMuon/AOD/17Nov2017-v1/20000/007A4914-57D3-E711-AD6E-5065F3819241.root')



)

process.daphne = cms.EDAnalyzer('DaphneAnalyzer'
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string('daphneTree.root'),
      closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.daphne)
