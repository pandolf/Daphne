import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100


process.source = cms.Source("PoolSource",
    ## replace 'myfile.root' with the source file you want to use
    #fileNames = cms.untracked.vstring(
    #    'file:myfile.root'
    #)
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/SingleNeutrino/AODSIM/94X_mc2017_realistic_v10-v1/60000/00232899-6ED7-E711-9355-0025905A4964.root')
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseISpring17DR/SingleNeutrino/GEN-SIM-RECO/FlatPU0to75RECO_90X_upgrade2017_realistic_v20-v1/120000/001DD274-7240-E711-BB47-A0369FC5EEF4.root')
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/50000/0024A35F-86EB-E711-A5BC-A4BF0112BC06.root')
    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/50000/0024A35F-86EB-E711-A5BC-A4BF0112BC06.root')

    #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIWinter17DR/QCD_Pt-15to3000_TuneCP5_Flat_13TeV_pythia8/GEN-SIM-RAW/NZSPU0to70_94X_upgrade2018_realistic_v8-v1/00000/0003F19E-74FF-E711-B73B-FA163E697F46.root')

    fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIIFall17GS/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/GEN-SIM/93X_mc2017_realistic_v3-v1/150000/008099F0-81B2-E711-8F1A-0CC47A4C8EE8.root')



)

process.demo = cms.EDAnalyzer('GenParticleAnalyzer'
)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string('genTree.root'),
      closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(process.demo)
