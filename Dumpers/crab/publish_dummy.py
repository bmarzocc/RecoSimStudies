import FWCore.ParameterSet.Config as cms

process = cms.Process("PUBLISH")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGunPt1-100_pythia8_PU_Flat55To75_14TeV_112X_mcRun3_2021_realistic_v15_GEN-SIM-RAW/cluster_job3235_step2.root'))

# Output definition
process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
                                         splitLevel = cms.untracked.int32(0),
                                         outputCommands = cms.untracked.vstring("keep *"),
                                         fileName = cms.untracked.string("cluster_job3235_step2.root")
)

process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.schedule = cms.Schedule(process.RECOSIMoutput_step)
