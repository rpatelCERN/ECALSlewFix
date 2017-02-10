import FWCore.ParameterSet.Config as cms
import sys

index=int(sys.argv[2])
process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/fdata/hepx/store/user/rish/CombineCards/EGM/ReMiniAOD/CMSSW_8_0_25/src/PAT_%d.root' %index
    )
)
process.TFileService = cms.Service("TFileService", fileName = cms.string("histo%d.root" %index) )
process.PatJetFix = cms.EDProducer('PatJetFix',
electronsFixed=cms.InputTag("slimmedElectrons"),
photonsFixed=cms.InputTag("slimmedPhotons"),
electrons=cms.InputTag("slimmedElectronsBeforeGSFix"),
PackedPart=cms.InputTag("packedPFCandidates"),
photons=cms.InputTag("slimmedPhotonsBeforeGSFix"),
jets=cms.InputTag("slimmedJets"),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)

  
process.p = cms.Path(process.PatJetFix)
process.e = cms.EndPath(process.out)
