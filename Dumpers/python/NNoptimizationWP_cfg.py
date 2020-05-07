import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFile_EB = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v7_detseparated/working_points_data_EB_modelv7.root'),
    inputFile_EE = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v8_shrinkwindows/working_point_data_EE_modelv8.root'), 
    inputFileMustache_EB = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v7_detseparated/working_points_data_EB_mustache.root'),
    inputFileMustache_EE = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v8_shrinkwindows/working_point_data_EE_mustache.root'),
    
    ## maxEvents
    maxEvents = cms.untracked.int32(-1),
 
    ##fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'),   
 
    etBinning = cms.string('0 10 20 30 40 50 60 70 80 90 100'),
    etaBinning = cms.string('-3.0 -2.5 -1.75 -1.479 -0.75 0.0 0.75 1.479 1.75 2.5 3.0'),
    dnnBinning = cms.string('0.3 0.32413793 0.34827586 0.37241379 0.39655172 0.42068966 0.4448275 0.4689655 0.4931034 0.51724138 0.54137931 0.56551724 0.58965517 0.6137931 0.63793103 0.66206897 0.6862069 0.71034483 0.73448276 0.75862069 0.78275862 0.80689655 0.83103448 0.85517241 0.87931034 0.90344828 0.92758621 0.95172414 0.97586207' )

)
