import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.MessageLogger = cms.Service("MessageLogger",
#                                              destinations   = cms.untracked.vstring('detailedInfo'),
#                                              categories      = cms.untracked.vstring('eventNumber'),
#                                              detailedInfo    = cms.untracked.PSet(
#                                                                        eventNumber = cms.untracked.PSet(
#                                                                                                  reportEvery = cms.untracked.int32(100)
#                                                                                                )
#                                                                                                                                 ),
#)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    ## replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
    #    'file:myfile.root'
    #)
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/SingleNeutrino/AODSIM/94X_mc2017_realistic_v10-v1/60000/00232899-6ED7-E711-9355-0025905A4964.root')
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseISpring17GS/SingleNeutrino/GEN-SIM/90X_upgrade2017_realistic_v20_ext1-v1/100000/E8051D36-B52F-E711-B317-4C79BA18103F.root')
    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseISpring17DR/SingleNeutrino/GEN-SIM-RECO/FlatPU0to75RECO_90X_upgrade2017_realistic_v20-v1/120000/001DD274-7240-E711-BB47-A0369FC5EEF4.root')





)

process.demo = cms.EDAnalyzer('GenParticleAnalyzer'
)

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('genTree.root')
                                   )

process.p = cms.Path(process.demo)
