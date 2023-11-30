import FWCore.ParameterSet.Config as cms

def dumperStepData(process):
 
 #dumper step
 process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
 )
 process.load('RecoSimStudies.Dumpers.RecoSimDumperData_cfi')
 process.dumper_step = cms.Path(process.recosimdumper)
 
 process.schedule.remove(process.MINIAODoutput_step)
 process.schedule += cms.Schedule(process.dumper_step)

 return process

def dumperStepMC_training(process):

 #dumper step
 process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string(options.outputFile)
                                   fileName = cms.string('output.root')
 )
 process.load('RecoSimStudies.Dumpers.RecoSimDumperMC_training_cfi')
 process.dumper_step = cms.Path(process.recosimdumper)

 process.schedule.remove(process.MINIAODSIMoutput_step)
 process.schedule += cms.Schedule(process.dumper_step)

 return process

def idealICs_UL18(process):

 process.myICs = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalIntercalibConstantsRcd'),
             tag = cms.string('EcalIntercalibConstants_MC_Digi_2018')
         )
     )
 )
 process.es_prefer_icReco = cms.ESPrefer("PoolDBESSource","myICs")

 return process

def noise_UL18(process):

 process.myPedestals = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalPedestalsRcd'),
             tag = cms.string('EcalPedestals_UL_2018_mc')
         )
     )
 )
 process.es_prefer_pedestals = cms.ESPrefer("PoolDBESSource","myPedestals")

 return process

def noise_235fb(process):

 process.myICs = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalIntercalibConstantsRcd'),
             tag = cms.string('EcalIntercalibConstants_MC_Digi_2018')
         )
     )
 )
 process.es_prefer_icReco = cms.ESPrefer("PoolDBESSource","myICs")
 
 return process

def pfThres_UL18(process):

 process.myPFRechitThres = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalPFRecHitThresholdsRcd'),
             tag = cms.string('EcalPFRecHitThresholds_UL_2018_2e3sig')
         )
     )
 )
 process.es_prefer_pfRechitThres = cms.ESPrefer("PoolDBESSource","myPFRechitThres")

 return process

def pfThres_235fb(process):

 process.myPFRechitThres = cms.ESSource("PoolDBESSource",
     connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
     toGet = cms.VPSet(
         cms.PSet(
             record = cms.string('EcalPFRecHitThresholdsRcd'),
             tag = cms.string('EcalPFRecHitThresholds_34sigma_TL235')
         )
     )
 )
 process.es_prefer_pfRechitThres = cms.ESPrefer("PoolDBESSource","myPFRechitThres")

 return process


