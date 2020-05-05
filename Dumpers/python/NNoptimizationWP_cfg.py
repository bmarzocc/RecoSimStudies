import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.ioFilesOpt = cms.PSet(

    ##input file
    inputFile_EB = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v5_dynamicwindow/result_scan_model_v5_EB.root'),
    inputFile_EE = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v5_dynamicwindow/result_scan_model_v5_EE.root'), 
    inputFileMustache_EB = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v5_dynamicwindow/result_scan_mustache_EB.root'),
    inputFileMustache_EE = cms.string('/eos/user/d/dvalsecc/ECAL/EcalClustering/DeepCluster/models/v5_dynamicwindow/result_scan_mustache_EE.root'),
    
    ##fitFunction
    fitFunction = cms.string('cruijff'),
    #fitFunction = cms.string('doubleCB'),   
 
    dnnBinning = cms.string('0.5 0.51724138 0.53448276 0.55172414 0.56896552 0.5862069 0.60344828 0.62068966 0.63793103 0.65517241 0.67241379 0.68965517 0.70689655 0.72413793 0.74137931 0.75862069 0.77586207 0.79310345 0.81034483 0.82758621 0.84482759 0.86206897 0.87931034 0.89655172 0.9137931 0.93103448 0.94827586 0.96551724 0.98275862' )
)
