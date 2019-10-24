
### Example for running

```
cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/test_5000_updated4_EB_doubleGamma.root inputFiles=file:test/test_5000_updated4_EB_douleGamma/step3.root
```
```
cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/photon_Et1to100GeV_closeEcal_EE_noPU_pfrh1.0_seed3.0_V01_v52_n15000.root inputFiles_load=data/samples/photon_Et1to100GeV_closeEcal_EE_noPU_pfrh1.0_seed3.0_V01_v52_n15000.txt
```

### EB no PU
cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/test_EB_noPU.root inputFiles=/store/user/mratti/EcalProd/photon_Et1to100GeV_closeEcal_EB_noPU_pfrh1_seed3_V01_v01_n10/step3.root

### EB w/ PU
cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/test_EB_wPU.root inputFiles=/store/user/mratti/EcalProd/photon_Et1to100GeV_closeEcal_EB_wPU_pfrh1_seed3_V01_v03_n10/step3.root


### EE no PU
cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/test_EE_noPU.root inputFiles=/store/user/mratti/EcalProd/photon_Et1to100GeV_closeEcal_EE_noPU_pfrh1_seed3_V01_v02_n10/step3.root

### EE w/ PU
cmsRun python/Cfg_RecoSimDumper_cfg.py outputFile=test/outputfiles/test_EE_wPU.root inputFiles=/store/user/mratti/EcalProd/photon_Et1to100GeV_closeEcal_EE_wPU_pfrh1_seed3_V01_v04_n10/step3.root

### Large sample EB no PU
cmsRun python/Cfg_RecoSimDumper_cfg.py  outputFile=test/outputfiles/photon_Et1to100GeV_closeEcal_EB_noPU_pfrh1_seed3_V01_v01_n15000.root inputFiles=/store/user/mratti/EcalProd/photon_Et1to100GeV_closeEcal_EB_noPU_pfrh1_seed3_V01_v01_n15000/step3.root

### Large sample EB w/PU
cmsRun python/Cfg_RecoSimDumper_cfg.py  outputFile=test/outputfiles/photon_Et1to100GeV_closeEcal_EB_wPU_pfrh1_seed3_V01_v03_n15000.root inputFiles=/store/user/mratti/EcalProd/photon_Et1to100GeV_closeEcal_EB_wPU_pfrh1_seed3_V01_v03_n15000/step3.root

