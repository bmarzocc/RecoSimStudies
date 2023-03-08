#ifndef RecoSimStudies_Dumpers_RecoSimDumper_H
#define RecoSimStudies_Dumpers_RecoSimDumper_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/Math/interface/libminifloat.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/DataRecord/interface/EcalLaserAlphasRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatiosRef.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRefRcd.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAPDPNRatios.h"
#include "CondFormats/DataRecord/interface/EcalLaserAPDPNRatiosRcd.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"

#include "DataFormats/Math/interface/deltaR.h"
//#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"

#include <memory>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <functional>
#include <set>
#include <assert.h>
#include <time.h>

#include <TMath.h>
#include <Math/VectorUtil.h>
//#include <boost/tokenizer.hpp>

class RecoSimDumper : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
      public:
         explicit RecoSimDumper(const edm::ParameterSet&);
	 ~RecoSimDumper();
  
  
      private:
	 void beginJob() override;
	 void analyze(const edm::Event&, const edm::EventSetup&) override;
         void endJob() override;
        
      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      void setTree(TTree* tree);
      void setVectors(int nGenParticles, int nCaloParticles, int nPFClusters, int nSuperClustersEB, int nSuperClustersEE, int nRetunedSuperClustersEB, int nRetunedSuperClustersEE, int nDeepSuperClustersEB, int nDeepSuperClustersEE); 
      double ptFast(const double energy, const math::XYZPoint& position, const math::XYZPoint& origin);
      int getGenStatusFlag(const reco::GenParticle* genParticle);
      int getGenMother(const reco::GenParticle* genParticle);
      int getGenParton(const std::vector<reco::GenParticle>* genParticles, const int genIndex);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<std::pair<DetId,std::pair<float,float>>>* getSharedHitsAndEnergies(const std::vector<std::pair<DetId, float> >* hitsAndEnergies1, const std::vector<std::pair<DetId, float> >* hitsAndEnergies2);
      std::pair<double,double> calculateCovariances(const reco::PFCluster* pfCluster, const EcalRecHitCollection* recHits, const CaloSubdetectorGeometry* geometry);
      std::vector<float> getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology);
      std::vector<double> getScores(const reco::PFCluster* pfCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<double> getScores(const reco::SuperCluster* superCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<double> getNoise(const reco::PFCluster* pfCluster, const std::vector<std::vector<std::pair<DetId, float> >> *hits_and_energies_CaloPart, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPartPU, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE, const EcalLaserAlphas* laserAlpha, const EcalLaserAPDPNRatios* laserRatio, const EcalIntercalibConstants* ical, const EcalIntercalibConstants* icalMC, const EcalPedestals* ped, const EcalADCToGeVConstant* adcToGeV, const EcalGainRatios* gr, bool useFractions);
      int getMatchedIndex(std::vector<std::vector<double>>* score, double selection, bool useMax, double scale, int iCl);
      void fillParticleMatchedIndex(std::vector<std::vector<int>>* particleMatchedIndex, std::vector<int>* clusterMatchedIndex);
      GlobalPoint calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES);
      bool passesJedID(const pat::Jet* jet, std::string version, std::string algo, std::string quality);  
       

      // ----------collection tokens-------------------
      edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
      edm::ESGetToken<EcalADCToGeVConstant, EcalADCToGeVConstantRcd> ADCtoGeVToken_;
      edm::ESGetToken<EcalLaserAlphas, EcalLaserAlphasRcd> alphaToken_;
      edm::ESGetToken<EcalLaserAPDPNRatios, EcalLaserAPDPNRatiosRcd> APDPNRatiosToken_;
      edm::ESGetToken<EcalIntercalibConstants, EcalIntercalibConstantsRcd> icalToken_;
      edm::ESGetToken<EcalIntercalibConstantsMC, EcalIntercalibConstantsMCRcd> icalMCToken_;
      edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> channelStatusToken_;
      edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> pedsToken_;
      edm::ESGetToken<EcalGainRatios, EcalGainRatiosRcd> ratioToken_;

      const CaloGeometry* geometry;  
      const CaloTopology* topology;
      const EcalADCToGeVConstant* adcToGeV;
      const EcalLaserAlphas* laserAlpha;
      const EcalLaserAPDPNRatios* laserRatio;
      const EcalIntercalibConstants* ical;
      const EcalIntercalibConstantsMC* icalMC;
      const EcalChannelStatus* chStatus;
      const EcalPedestals* ped; 
      const EcalGainRatios* gr;
      
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<std::vector<CaloParticle> > puCaloPartToken_;
      edm::EDGetTokenT<std::vector<CaloParticle> > ootpuCaloPartToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFRecHit>  > pfRecHitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCluster> > pfClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::GsfElectron> > gsfElectronToken_;
      edm::EDGetTokenT<std::vector<reco::Photon> > gedPhotonToken_;
      edm::EDGetTokenT<std::vector<pat::Electron> > patElectronToken_;
      edm::EDGetTokenT<std::vector<pat::Photon> > patPhotonToken_;
      edm::EDGetTokenT<std::vector<pat::Jet> > patJetToken_;
      edm::EDGetTokenT<std::vector<pat::MET> > patMETToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeSuperClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebRetunedSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeRetunedSuperClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebDeepSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeDeepSuperClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebDeepSuperClusterLWPToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeDeepSuperClusterLWPToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebDeepSuperClusterTWPToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeDeepSuperClusterTWPToken_;
      
      edm::Service<TFileService> iFile;
      const CaloSubdetectorGeometry* _ebGeom;
      const CaloSubdetectorGeometry* _eeGeom;
      const CaloSubdetectorGeometry* _esGeom;
      bool _esPlus;
      bool _esMinus;
      
      // ----------config inputs-------------------
      bool isMC_;
      bool doCompression_;
      int nBits_;
      bool saveGenParticles_; 
      std::vector<int> genParticlesToSave_;      
      bool saveCaloParticles_;  
      bool saveCaloParticlesPU_; 
      bool saveCaloParticlesOOTPU_;    
      bool subtractSignalCalo_;    
      bool saveSimhits_;       
      bool saveSimhitsPU_;          
      bool saveRechits_;          
      bool savePFRechits_;   
      bool savePFCluster_;    
      bool savePFClusterhits_;   
      bool saveShowerShapes_;  
      bool saveSuperCluster_;
      bool saveRetunedSC_;
      bool saveDeepSC_;  
      bool saveGsfElectrons_;  
      bool saveGedPhotons_;  
      bool savePatPhotons_; 
      bool savePatElectrons_; 
      bool savePatJets_;      
      
      std::string egmCutBasedElectronIDVeto_;
      std::string egmCutBasedElectronIDloose_;
      std::string egmCutBasedElectronIDmedium_;
      std::string egmCutBasedElectronIDtight_; 
      std::string egmMVAElectronIDloose_;
      std::string egmMVAElectronIDmedium_;
      std::string egmMVAElectronIDtight_; 
      std::string egmMVAElectronIDlooseNoIso_;
      std::string egmMVAElectronIDmediumNoIso_;
      std::string egmMVAElectronIDtightNoIso_; 
      std::string heepElectronID_; 
      std::string egmCutBasedPhotonIDloose_;
      std::string egmCutBasedPhotonIDmedium_;
      std::string egmCutBasedPhotonIDtight_;  
      std::string egmMVAPhotonIDmedium_;
      std::string egmMVAPhotonIDtight_; 
      std::vector<std::string> jmeCutBasedPFJetIDloose_;         
      std::vector<std::string> jmeCutBasedPFJetIDtight_;       
      std::vector<std::string> jmeCutBasedPFJetIDtightLepVeto_;
      std::map<std::string,std::map<std::string,std::map<std::string,std::map<std::string,bool> > > > jetIdMap;
      
      // ----------DNN inputs-------------------
      std::vector<double> HLF_VectorVar_;
      std::vector<std::vector<double>> PL_VectorVar_;
      std::vector<double> x_mean_, x_std_, list_mean_, list_std_;
     
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      std::vector<std::map<uint32_t,float> > caloParticleXtals_;
      std::pair<double,double> widths_;
      std::array<float,3> locCov_;
      std::array<float,3> full5x5_locCov_;
      std::vector<float> showerShapes_;
      std::map<int,std::vector<int>> genAllDaughters;
      
      long int eventId;
      int lumiId;
      int runId; 
      double truePU;
      double obsPU;
      int nVtx;
      float rho; 
      int genParticle_size; 
      int caloParticle_size;
      std::vector<int> genParticle_genMotherIndex;
      std::vector<std::vector<int> > genParticle_genDaughtersIndex;
      std::vector<int> genParticle_pdgId;
      std::vector<int> genParticle_status; 
      std::vector<int> genParticle_statusFlag; 
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;
      std::vector<std::vector<int> > genParticle_pfCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_superCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_retunedSuperCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_deepSuperCluster_dR_genScore_MatchedIndex; 
      int caloParticlePU_nHitsWithES; 
      int caloParticlePU_nHits; 
      float caloParticlePU_totEnergyWithES;
      float caloParticlePU_totEnergy;
      std::vector<float> caloParticlePU_xtalEnergy;
      std::vector<float> caloParticlePU_xtalEta;
      std::vector<float> caloParticlePU_xtalPhi;
      std::vector<int> caloParticlePU_xtalIeta;
      std::vector<int> caloParticlePU_xtalIphi; 
      std::vector<int> caloParticlePU_xtalIz;
      std::vector<int> caloParticlePU_xtalIplane;
      int caloParticleOOTPU_nHitsWithES; 
      int caloParticleOOTPU_nHits; 
      float caloParticleOOTPU_totEnergyWithES;
      float caloParticleOOTPU_totEnergy;
      std::vector<float> caloParticleOOTPU_xtalEnergy;
      std::vector<float> caloParticleOOTPU_xtalEta;
      std::vector<float> caloParticleOOTPU_xtalPhi;
      std::vector<int> caloParticleOOTPU_xtalIeta;
      std::vector<int> caloParticleOOTPU_xtalIphi; 
      std::vector<int> caloParticleOOTPU_xtalIz;
      std::vector<int> caloParticleOOTPU_xtalIplane;
      std::vector<int> caloParticle_index; 
      std::vector<int> caloParticle_nXtals;  
      std::vector<int> caloParticle_pdgId;
      std::vector<int> caloParticle_status;
      std::vector<int> caloParticle_charge;
      std::vector<float> caloParticle_genEnergy;
      std::vector<float> caloParticle_simEnergy;  
      std::vector<float> caloParticle_simEnergyGoodStatus;   
      std::vector<float> caloParticle_simEnergyWithES;
      std::vector<float> caloParticle_genPt;
      std::vector<float> caloParticle_simPt;
      std::vector<float> caloParticle_genEta;
      std::vector<float> caloParticle_simEta;
      std::vector<float> caloParticle_genPhi;
      std::vector<float> caloParticle_simPhi;
      std::vector<int> caloParticle_partonIndex;
      std::vector<int> caloParticle_partonPdgId; 
      std::vector<int> caloParticle_partonCharge;  
      std::vector<float> caloParticle_partonEnergy;
      std::vector<float> caloParticle_partonPt;
      std::vector<float> caloParticle_partonEta;
      std::vector<float> caloParticle_partonPhi; 
      std::vector<int> caloParticle_genMotherPdgId;
      std::vector<int> caloParticle_genMotherStatus; 
      std::vector<int> caloParticle_genMotherCharge;
      std::vector<float> caloParticle_genMotherEnergy;
      std::vector<float> caloParticle_genMotherPt;
      std::vector<float> caloParticle_genMotherEta;
      std::vector<float> caloParticle_genMotherPhi; 
      std::vector<int> caloParticle_simIeta;
      std::vector<int> caloParticle_simIphi;
      std::vector<int> caloParticle_simIz;
      std::vector<int> caloParticle_nSharedXtals;
      std::vector<int> caloParticle_sharedIndex1;
      std::vector<int> caloParticle_sharedIndex2;
      std::vector<float> caloParticle_sharedEnergyFrac1;
      std::vector<float> caloParticle_sharedEnergyFrac2;
      std::vector<std::vector<int> > caloParticle_pfCluster_dR_simScore_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_pfCluster_recoToSim_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<std::vector<int> > caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex;      
      std::vector<std::vector<int> > caloParticle_superCluster_dR_simScore_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_superCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_superCluster_recoToSim_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<std::vector<int> > caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex;      
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex;      
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_dR_simScore_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_sim_fraction_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<std::vector<int> > caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex;      
      std::vector<std::vector<float> > simHit_energy;
      std::vector<std::vector<float> > simHit_eta;
      std::vector<std::vector<float> > simHit_phi;
      std::vector<std::vector<int> > simHit_ieta;
      std::vector<std::vector<int> > simHit_iphi;
      std::vector<std::vector<int> > simHit_iz;
      std::vector<std::vector<int> > simHit_iplane; 
      std::vector<std::vector<int> > simHit_chStatus; 
      std::vector<float> recHit_noPF_energy;
      std::vector<float> recHit_noPF_eta;
      std::vector<float> recHit_noPF_phi;
      std::vector<int> recHit_noPF_ieta;
      std::vector<int> recHit_noPF_iphi;
      std::vector<int> recHit_noPF_iz;
      std::vector<float> pfRecHit_unClustered_energy;
      std::vector<float> pfRecHit_unClustered_eta;
      std::vector<float> pfRecHit_unClustered_phi;
      std::vector<int> pfRecHit_unClustered_ieta;
      std::vector<int> pfRecHit_unClustered_iphi;
      std::vector<int> pfRecHit_unClustered_iz;
      std::vector<std::vector<float> > pfClusterHit_fraction; 
      std::vector<std::vector<float> > pfClusterHit_rechitEnergy; 
      std::vector<std::vector<float> > pfClusterHit_eta;
      std::vector<std::vector<float> > pfClusterHit_phi;
      std::vector<std::vector<int> > pfClusterHit_ieta;
      std::vector<std::vector<int> > pfClusterHit_iphi;
      std::vector<std::vector<int> > pfClusterHit_iz;     
      std::vector<std::vector<float> > pfClusterHit_adcToGeV; 
      std::vector<std::vector<float> > pfClusterHit_laserCorr; 
      std::vector<std::vector<float> > pfClusterHit_ic; 
      std::vector<std::vector<float> > pfClusterHit_icMC; 
      std::vector<std::vector<int> > pfClusterHit_chStatus; 
      std::vector<float> pfCluster_rawEnergy;    
      std::vector<float> pfCluster_rawEnergyUncalib;
      std::vector<float> pfCluster_energy;
      std::vector<float> pfCluster_rawPt;
      std::vector<float> pfCluster_pt;
      std::vector<float> pfCluster_noise;
      std::vector<float> pfCluster_noiseUncalib;    
      std::vector<float> pfCluster_noiseNoFractions;
      std::vector<float> pfCluster_noiseUncalibNoFractions;
      std::vector<float> pfCluster_noiseDB;
      std::vector<float> pfCluster_noiseDBUncalib;    
      std::vector<float> pfCluster_noiseDBNoFractions;
      std::vector<float> pfCluster_noiseDBUncalibNoFractions;
      std::vector<float> pfCluster_eta;
      std::vector<float> pfCluster_phi; 
      std::vector<int> pfCluster_ieta;
      std::vector<int> pfCluster_iphi;
      std::vector<int> pfCluster_iz;
      std::vector<double> pfCluster_nXtals;
      std::vector<std::vector<int> > pfCluster_superClustersIndex;
      std::vector<std::vector<int> > pfCluster_retunedSuperClustersIndex;
      std::vector<std::vector<int> > pfCluster_deepSuperClustersIndex;
      std::vector<float> pfCluster_etaWidth; 
      std::vector<float> pfCluster_phiWidth; 
      std::vector<float> pfCluster_e5x5;
      std::vector<float> pfCluster_e2x2Ratio;
      std::vector<float> pfCluster_e3x3Ratio;
      std::vector<float> pfCluster_eMaxRatio;
      std::vector<float> pfCluster_e2ndRatio;
      std::vector<float> pfCluster_eTopRatio;
      std::vector<float> pfCluster_eRightRatio;
      std::vector<float> pfCluster_eBottomRatio;
      std::vector<float> pfCluster_eLeftRatio;
      std::vector<float> pfCluster_e2x5MaxRatio;
      std::vector<float> pfCluster_e2x5TopRatio;
      std::vector<float> pfCluster_e2x5RightRatio;
      std::vector<float> pfCluster_e2x5BottomRatio;
      std::vector<float> pfCluster_e2x5LeftRatio;
      std::vector<float> pfCluster_swissCross;
      std::vector<float> pfCluster_r9;
      std::vector<float> pfCluster_sigmaIetaIeta; 
      std::vector<float> pfCluster_sigmaIetaIphi; 
      std::vector<float> pfCluster_sigmaIphiIphi; 
      std::vector<float> pfCluster_full5x5_e5x5;
      std::vector<float> pfCluster_full5x5_e2x2Ratio;
      std::vector<float> pfCluster_full5x5_e3x3Ratio;
      std::vector<float> pfCluster_full5x5_eMaxRatio;
      std::vector<float> pfCluster_full5x5_e2ndRatio;
      std::vector<float> pfCluster_full5x5_eTopRatio;
      std::vector<float> pfCluster_full5x5_eRightRatio;
      std::vector<float> pfCluster_full5x5_eBottomRatio;
      std::vector<float> pfCluster_full5x5_eLeftRatio;
      std::vector<float> pfCluster_full5x5_e2x5MaxRatio;
      std::vector<float> pfCluster_full5x5_e2x5TopRatio;
      std::vector<float> pfCluster_full5x5_e2x5RightRatio;
      std::vector<float> pfCluster_full5x5_e2x5BottomRatio;
      std::vector<float> pfCluster_full5x5_e2x5LeftRatio;
      std::vector<float> pfCluster_full5x5_swissCross;
      std::vector<float> pfCluster_full5x5_r9;
      std::vector<float> pfCluster_full5x5_sigmaIetaIeta; 
      std::vector<float> pfCluster_full5x5_sigmaIetaIphi; 
      std::vector<float> pfCluster_full5x5_sigmaIphiIphi; 
      std::vector<int> pfCluster_dR_genScore_MatchedIndex;
      std::vector<int> pfCluster_dR_simScore_MatchedIndex;
      std::vector<int> pfCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_MatchedIndex;
      std::vector<int> pfCluster_recoToSim_fraction_MatchedIndex;
      std::vector<int> pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<int> pfCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<int> pfCluster_recoEnergy_sharedXtals_MatchedIndex;  
      std::vector<std::vector<double> > pfCluster_dR_genScore;
      std::vector<std::vector<double> > pfCluster_dR_simScore;
      std::vector<std::vector<double> > pfCluster_sim_nSharedXtals;
      std::vector<std::vector<double> > pfCluster_sim_fraction_noHitsFraction;
      std::vector<std::vector<double> > pfCluster_sim_fraction;
      std::vector<std::vector<double> > pfCluster_recoToSim_fraction;
      std::vector<std::vector<double> > pfCluster_recoToSim_fraction_sharedXtals;
      std::vector<std::vector<double> > pfCluster_simEnergy_sharedXtals; 
      std::vector<std::vector<double> > pfCluster_recoEnergy_sharedXtals;  
      std::vector<double> pfCluster_simPU_nSharedXtals;
      std::vector<double> pfCluster_recoEnergy_sharedXtalsPU;  
      std::vector<double> pfCluster_simEnergy_sharedXtalsPU; 
      std::vector<double> pfCluster_recoEnergy_noHitsFraction_sharedXtalsPU;  
      std::vector<double> pfCluster_simEnergy_noHitsFraction_sharedXtalsPU; 
      std::vector<double> pfCluster_simOOTPU_nSharedXtals;
      std::vector<double> pfCluster_recoEnergy_sharedXtalsOOTPU;  
      std::vector<double> pfCluster_simEnergy_sharedXtalsOOTPU; 
      std::vector<double> pfCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU;  
      std::vector<double> pfCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU; 
      std::vector<int> gsfElectron_index;
      std::vector<uint32_t> gsfElectron_seedRawId;
      std::vector<bool> gsfElectron_isEB;
      std::vector<bool> gsfElectron_isEE;
      std::vector<bool> gsfElectron_isEBEEGap;
      std::vector<bool> gsfElectron_isEBEtaGap;
      std::vector<bool> gsfElectron_isEBPhiGap;
      std::vector<bool> gsfElectron_isEEDeeGap;
      std::vector<bool> gsfElectron_isEERingGap;
      std::vector<bool> gsfElectron_isEcalDriven;
      std::vector<bool> gsfElectron_isTrackerDriven;
      std::vector<int> gsfElectron_classification;
      std::vector<int> gsfElectron_refinedSCNPFClusters;
      std::vector<int> gsfElectron_ecalSCNPFClusters;
      std::vector<float> gsfElectron_p;
      std::vector<float> gsfElectron_pt; 
      std::vector<float> gsfElectron_et;
      std::vector<float> gsfElectron_energy;
      std::vector<float> gsfElectron_energyErr;
      std::vector<float> gsfElectron_ecalEnergy;
      std::vector<float> gsfElectron_ecalEnergyErr;
      std::vector<float> gsfElectron_eta;
      std::vector<float> gsfElectron_phi;
      std::vector<float> gsfElectron_trkEtaMode;
      std::vector<float> gsfElectron_trkPhiMode;
      std::vector<float> gsfElectron_trkPMode;
      std::vector<float> gsfElectron_trkPModeErr;
      std::vector<float> gsfElectron_trkPInn;
      std::vector<float> gsfElectron_trkPtInn;
      std::vector<float> gsfElectron_trkPVtx;
      std::vector<float> gsfElectron_trkPOut;
      std::vector<float> gsfElectron_trkChi2;
      std::vector<float> gsfElectron_trkNDof;
      std::vector<float> gsfElectron_trackFbrem;
      std::vector<float> gsfElectron_superClusterFbrem;
      std::vector<float> gsfElectron_hademTow;
      std::vector<float> gsfElectron_hademCone;
      std::vector<float> gsfElectron_ecalDrivenSeed;
      std::vector<float> gsfElectron_nrSatCrys;
      std::vector<float> gsfElectron_refinedSCEnergy;
      std::vector<float> gsfElectron_refinedSCRawEnergy;
      std::vector<float> gsfElectron_refinedSCRawESEnergy;
      std::vector<float> gsfElectron_refinedSCEt;
      std::vector<float> gsfElectron_refinedSCPhiWidth;
      std::vector<float> gsfElectron_refinedSCEtaWidth;
      std::vector<float> gsfElectron_refinedSCEoP;
      std::vector<float> gsfElectron_ecalSCEnergy;
      std::vector<float> gsfElectron_ecalSCRawEnergy;
      std::vector<float> gsfElectron_ecalSCRawESEnergy;
      std::vector<float> gsfElectron_ecalSCEt;
      std::vector<float> gsfElectron_ecalSCPhiWidth;
      std::vector<float> gsfElectron_ecalSCEtaWidth;
      std::vector<float> gsfElectron_ecalSCEoP;
      std::vector<float> gsfElectron_scEta;
      std::vector<float> gsfElectron_scPhi;
      std::vector<float> gsfElectron_scSwissCross;
      std::vector<float> gsfElectron_scE2x2;
      std::vector<float> gsfElectron_scE3x3; 
      std::vector<float> gsfElectron_scE5x5; 
      std::vector<float> gsfElectron_scEMax;
      std::vector<float> gsfElectron_scR9; 
      std::vector<float> gsfElectron_scSigmaIEtaIEta;
      std::vector<float> gsfElectron_scSigmaIEtaIPhi;
      std::vector<float> gsfElectron_scSigmaIPhiIPhi;
      std::vector<float> gsfElectron_full5x5_scE2x2;
      std::vector<float> gsfElectron_full5x5_scE3x3;
      std::vector<float> gsfElectron_full5x5_scE5x5;
      std::vector<float> gsfElectron_full5x5_scEMax;
      std::vector<float> gsfElectron_full5x5_scR9;  
      std::vector<float> gsfElectron_full5x5_scSigmaIEtaIEta;
      std::vector<float> gsfElectron_full5x5_scSigmaIEtaIPhi;
      std::vector<float> gsfElectron_full5x5_scSigmaIPhiIPhi;
      std::vector<float> gsfElectron_HoE; 
      std::vector<float> gsfElectron_trkIso03; 
      std::vector<float> gsfElectron_ecalIso03;  
      std::vector<float> gsfElectron_hcalIso03;
      std::vector<float> gsfElectron_trkIso04; 
      std::vector<float> gsfElectron_ecalIso04;  
      std::vector<float> gsfElectron_hcalIso04;
      std::vector<float> gsfElectron_pfPhotonIso;     
      std::vector<float> gsfElectron_pfChargedHadronIso;     
      std::vector<float> gsfElectron_pfNeutralHadronIso;
      std::vector<float> gsfElectron_mva_Isolated; 
      std::vector<float> gsfElectron_mva_e_pi; 
      std::vector<float> gsfElectron_dnn_signal_Isolated;
      std::vector<float> gsfElectron_dnn_signal_nonIsolated;  
      std::vector<float> gsfElectron_dnn_bkg_nonIsolated; 
      std::vector<float> gsfElectron_dnn_bkg_Tau;
      std::vector<float> gsfElectron_dnn_bkg_Photon;
      std::vector<int> gedPhoton_index;
      std::vector<uint32_t> gedPhoton_seedRawId;
      std::vector<bool> gedPhoton_isEB;
      std::vector<bool> gedPhoton_isEE;
      std::vector<bool> gedPhoton_isEBEEGap;
      std::vector<bool> gedPhoton_isEBEtaGap;
      std::vector<bool> gedPhoton_isEBPhiGap;
      std::vector<bool> gedPhoton_isEEDeeGap;
      std::vector<bool> gedPhoton_isEERingGap;
      std::vector<int> gedPhoton_refinedSCNPFClusters;
      std::vector<int> gedPhoton_ecalSCNPFClusters;
      std::vector<bool> gedPhoton_hasConversionTracks;
      std::vector<int> gedPhoton_nConversions;
      std::vector<int> gedPhoton_nConversionsOneLeg;
      std::vector<float> gedPhoton_et; 
      std::vector<float> gedPhoton_energy;
      std::vector<float> gedPhoton_energyErr;
      std::vector<float> gedPhoton_ecalEnergy;
      std::vector<float> gedPhoton_ecalEnergyErr;
      std::vector<float> gedPhoton_eta;
      std::vector<float> gedPhoton_phi;
      std::vector<float> gedPhoton_hademTow;
      std::vector<float> gedPhoton_hademCone;
      std::vector<float> gedPhoton_nrSatCrys;    
      std::vector<float> gedPhoton_refinedSCEnergy;  
      std::vector<float> gedPhoton_refinedSCRawEnergy;  
      std::vector<float> gedPhoton_refinedSCRawESEnergy;
      std::vector<float> gedPhoton_refinedSCEt;
      std::vector<float> gedPhoton_refinedSCEtaWidth;
      std::vector<float> gedPhoton_refinedSCPhiWidth; 
      std::vector<float> gedPhoton_ecalSCEnergy;  
      std::vector<float> gedPhoton_ecalSCRawEnergy;  
      std::vector<float> gedPhoton_ecalSCRawESEnergy;
      std::vector<float> gedPhoton_ecalSCEt;
      std::vector<float> gedPhoton_ecalSCEtaWidth;
      std::vector<float> gedPhoton_ecalSCPhiWidth; 
      std::vector<float> gedPhoton_scEta;
      std::vector<float> gedPhoton_scPhi; 
      std::vector<float> gedPhoton_scSwissCross;
      std::vector<float> gedPhoton_scE2x2;
      std::vector<float> gedPhoton_scE3x3; 
      std::vector<float> gedPhoton_scE5x5; 
      std::vector<float> gedPhoton_scEMax;
      std::vector<float> gedPhoton_scR9;
      std::vector<float> gedPhoton_scSigmaIEtaIEta;
      std::vector<float> gedPhoton_scSigmaIEtaIPhi;
      std::vector<float> gedPhoton_scSigmaIPhiIPhi;
      std::vector<float> gedPhoton_full5x5_scE2x2;
      std::vector<float> gedPhoton_full5x5_scE3x3;
      std::vector<float> gedPhoton_full5x5_scE5x5;
      std::vector<float> gedPhoton_full5x5_scEMax;
      std::vector<float> gedPhoton_full5x5_scR9;
      std::vector<float> gedPhoton_full5x5_scSigmaIEtaIEta;
      std::vector<float> gedPhoton_full5x5_scSigmaIEtaIPhi;
      std::vector<float> gedPhoton_full5x5_scSigmaIPhiIPhi;
      std::vector<float> gedPhoton_HoE;
      std::vector<float> gedPhoton_trkIso03; 
      std::vector<float> gedPhoton_ecalIso03;  
      std::vector<float> gedPhoton_hcalIso03;
      std::vector<float> gedPhoton_trkIso04; 
      std::vector<float> gedPhoton_ecalIso04;  
      std::vector<float> gedPhoton_hcalIso04;
      std::vector<float> gedPhoton_pfPhotonIso;     
      std::vector<float> gedPhoton_pfChargedHadronIso;     
      std::vector<float> gedPhoton_pfNeutralHadronIso;
      std::vector<int> gedPhoton_nClusterOutsideMustache; 
      std::vector<float> gedPhoton_etOutsideMustache; 
      std::vector<float> gedPhoton_pfMVA; 
      std::vector<float> gedPhoton_pfDNN; 
      float patMET_sumEt;
      float patMET_et;
      float mll;
      std::vector<int> patElectron_index; 
      std::vector<uint32_t> patElectron_seedRawId;
      std::vector<int> patElectron_classification;
      std::vector<int> patElectron_refinedSCNPFClusters;
      std::vector<int> patElectron_ecalSCNPFClusters;
      std::vector<int> patElectron_charge;
      std::vector<bool> patElectron_isEB;
      std::vector<bool> patElectron_isEE;
      std::vector<bool> patElectron_isEBEEGap;
      std::vector<bool> patElectron_isEBEtaGap;
      std::vector<bool> patElectron_isEBPhiGap;
      std::vector<bool> patElectron_isEEDeeGap;
      std::vector<bool> patElectron_isEERingGap;
      std::vector<bool> patElectron_isEcalDriven;
      std::vector<bool> patElectron_isTrackerDriven;
      std::vector<bool> patElectron_passConversionVeto;
      std::vector<int> patElectron_nOverlapPhotons;
      std::vector<std::vector<int> > patElectron_overlapPhotonIndices;
      std::vector<bool> patElectron_hasOverlapJet;
      std::vector<int> patElectron_overlapJetIndex; 
      std::vector<float> patElectron_eta;
      std::vector<float> patElectron_phi;
      std::vector<float> patElectron_p;
      std::vector<float> patElectron_pt;
      std::vector<float> patElectron_pIn;
      std::vector<float> patElectron_pOut;
      std::vector<float> patElectron_pAtCalo;
      std::vector<float> patElectron_deltaEtaIn;
      std::vector<float> patElectron_deltaPhiIn;
      std::vector<float> patElectron_deltaEtaSeedClusterAtCalo;
      std::vector<float> patElectron_deltaEtaEleClusterAtCalo;
      std::vector<float> patElectron_deltaPhiEleClusterAtCalo;
      std::vector<float> patElectron_deltaPhiSeedClusterAtCalo;
      std::vector<int> patElectron_misHits;
      std::vector<int> patElectron_nAmbiguousGsfTracks; 
      std::vector<float> patElectron_trackFbrem;
      std::vector<float> patElectron_superClusterFbrem;
      std::vector<float> patElectron_dz;
      std::vector<float> patElectron_dzError;
      std::vector<float> patElectron_dxy;
      std::vector<float> patElectron_dxyError;
      std::vector<float> patElectron_energy;
      std::vector<float> patElectron_energyErr;
      std::vector<float> patElectron_ecalEnergy;
      std::vector<float> patElectron_ecalEnergyErr;
      std::vector<float> patElectron_et;
      std::vector<float> patElectron_mt;
      std::vector<float> patElectron_dphiMET;
      std::vector<float> patElectron_refinedSCEnergy;
      std::vector<float> patElectron_refinedSCRawEnergy;
      std::vector<float> patElectron_refinedSCRawESEnergy;
      std::vector<float> patElectron_refinedSCEt;
      std::vector<float> patElectron_refinedSCPhiWidth;
      std::vector<float> patElectron_refinedSCEtaWidth;
      std::vector<float> patElectron_refinedSCEoP;
      std::vector<float> patElectron_ecalSCEnergy;
      std::vector<float> patElectron_ecalSCRawEnergy;
      std::vector<float> patElectron_ecalSCRawESEnergy;
      std::vector<float> patElectron_ecalSCEt;
      std::vector<float> patElectron_ecalSCPhiWidth;
      std::vector<float> patElectron_ecalSCEtaWidth;
      std::vector<float> patElectron_ecalSCEoP;
      std::vector<float> patElectron_scEta;
      std::vector<float> patElectron_scPhi;
      std::vector<float> patElectron_scSwissCross;
      std::vector<float> patElectron_scE2x2;
      std::vector<float> patElectron_scE3x3; 
      std::vector<float> patElectron_scE5x5; 
      std::vector<float> patElectron_scEMax;
      std::vector<float> patElectron_scR9; 
      std::vector<float> patElectron_scSigmaIEtaIEta;
      std::vector<float> patElectron_scSigmaIEtaIPhi;
      std::vector<float> patElectron_scSigmaIPhiIPhi;
      std::vector<float> patElectron_full5x5_scE2x2;
      std::vector<float> patElectron_full5x5_scE3x3;
      std::vector<float> patElectron_full5x5_scE5x5;
      std::vector<float> patElectron_full5x5_scEMax;
      std::vector<float> patElectron_full5x5_scR9;  
      std::vector<float> patElectron_full5x5_scSigmaIEtaIEta;
      std::vector<float> patElectron_full5x5_scSigmaIEtaIPhi;
      std::vector<float> patElectron_full5x5_scSigmaIPhiIPhi;
      std::vector<float> patElectron_HoE;
      std::vector<float> patElectron_trkIso03;
      std::vector<float> patElectron_ecalIso03;
      std::vector<float> patElectron_hcalIso03;
      std::vector<float> patElectron_trkIso04;
      std::vector<float> patElectron_ecalIso04;
      std::vector<float> patElectron_hcalIso04;
      std::vector<float> patElectron_pfPhotonIso;
      std::vector<float> patElectron_pfNeutralHadronIso;
      std::vector<float> patElectron_pfChargedHadronIso;
      std::vector<float> patElectron_mva_Isolated; 
      std::vector<float> patElectron_mva_e_pi; 
      std::vector<float> patElectron_dnn_signal_Isolated;
      std::vector<float> patElectron_dnn_signal_nonIsolated;  
      std::vector<float> patElectron_dnn_bkg_nonIsolated; 
      std::vector<float> patElectron_dnn_bkg_Tau;
      std::vector<float> patElectron_dnn_bkg_Photon;
      std::vector<int> patElectron_egmCutBasedElectronIDVeto;
      std::vector<int> patElectron_egmCutBasedElectronIDloose;
      std::vector<int> patElectron_egmCutBasedElectronIDmedium;
      std::vector<int> patElectron_egmCutBasedElectronIDtight;
      std::vector<int> patElectron_egmMVAElectronIDloose;
      std::vector<int> patElectron_egmMVAElectronIDmedium;
      std::vector<int> patElectron_egmMVAElectronIDtight;
      std::vector<int> patElectron_egmMVAElectronIDlooseNoIso;
      std::vector<int> patElectron_egmMVAElectronIDmediumNoIso;
      std::vector<int> patElectron_egmMVAElectronIDtightNoIso;
      std::vector<int> patElectron_heepElectronID;
      std::vector<int> patPhoton_index; 
      std::vector<uint32_t> patPhoton_seedRawId;
      std::vector<int> patPhoton_refinedSCNPFClusters;
      std::vector<int> patPhoton_ecalSCNPFClusters;
      std::vector<bool> patPhoton_isEB;
      std::vector<bool> patPhoton_isEE;
      std::vector<bool> patPhoton_isEBEEGap;
      std::vector<bool> patPhoton_isEBEtaGap;
      std::vector<bool> patPhoton_isEBPhiGap;
      std::vector<bool> patPhoton_isEEDeeGap;
      std::vector<bool> patPhoton_isEERingGap;
      std::vector<bool> patPhoton_passElectronVeto;
      std::vector<bool> patPhoton_hasPixelSeed;
      std::vector<bool> patPhoton_hasConversionTracks;
      std::vector<int> patPhoton_nConversions;
      std::vector<int> patPhoton_nConversionsOneLeg; 
      std::vector<bool> patPhoton_hasOverlapElectron;
      std::vector<int> patPhoton_overlapElectronIndex;
      std::vector<bool> patPhoton_hasOverlapJet;
      std::vector<int> patPhoton_overlapJetIndex;  
      std::vector<float> patPhoton_eta;
      std::vector<float> patPhoton_phi;
      std::vector<float> patPhoton_energy; 
      std::vector<float> patPhoton_energyErr;
      std::vector<float> patPhoton_ecalEnergy;
      std::vector<float> patPhoton_ecalEnergyErr;
      std::vector<float> patPhoton_et;
      std::vector<float> patPhoton_mt;
      std::vector<float> patPhoton_dphiMET;     
      std::vector<float> patPhoton_refinedSCEnergy;  
      std::vector<float> patPhoton_refinedSCRawEnergy;  
      std::vector<float> patPhoton_refinedSCRawESEnergy;
      std::vector<float> patPhoton_refinedSCEt;
      std::vector<float> patPhoton_refinedSCEtaWidth;
      std::vector<float> patPhoton_refinedSCPhiWidth;    
      std::vector<float> patPhoton_ecalSCEnergy;  
      std::vector<float> patPhoton_ecalSCRawEnergy;  
      std::vector<float> patPhoton_ecalSCRawESEnergy;
      std::vector<float> patPhoton_ecalSCEt;
      std::vector<float> patPhoton_ecalSCEtaWidth;
      std::vector<float> patPhoton_ecalSCPhiWidth; 
      std::vector<float> patPhoton_scEta;
      std::vector<float> patPhoton_scPhi;
      std::vector<float> patPhoton_scSwissCross;
      std::vector<float> patPhoton_scE2x2;
      std::vector<float> patPhoton_scE3x3; 
      std::vector<float> patPhoton_scE5x5; 
      std::vector<float> patPhoton_scEMax;
      std::vector<float> patPhoton_scR9;
      std::vector<float> patPhoton_scSigmaIEtaIEta;
      std::vector<float> patPhoton_scSigmaIEtaIPhi;
      std::vector<float> patPhoton_scSigmaIPhiIPhi;
      std::vector<float> patPhoton_full5x5_scE2x2;
      std::vector<float> patPhoton_full5x5_scE3x3;
      std::vector<float> patPhoton_full5x5_scE5x5;
      std::vector<float> patPhoton_full5x5_scEMax;
      std::vector<float> patPhoton_full5x5_scR9;
      std::vector<float> patPhoton_full5x5_scSigmaIEtaIEta;
      std::vector<float> patPhoton_full5x5_scSigmaIEtaIPhi;
      std::vector<float> patPhoton_full5x5_scSigmaIPhiIPhi;
      std::vector<float> patPhoton_HoE;
      std::vector<float> patPhoton_trkIso03;
      std::vector<float> patPhoton_ecalIso03;
      std::vector<float> patPhoton_hcalIso03;
      std::vector<float> patPhoton_trkIso04;
      std::vector<float> patPhoton_ecalIso04;
      std::vector<float> patPhoton_hcalIso04;
      std::vector<float> patPhoton_patParticleIso;
      std::vector<float> patPhoton_pfChargedHadronIso;
      std::vector<float> patPhoton_pfNeutralHadronIso;
      std::vector<float> patPhoton_pfPhotonIso;
      std::vector<float> patPhoton_pfPuChargedHadronIso;
      std::vector<int> patPhoton_nClusterOutsideMustache; 
      std::vector<float> patPhoton_etOutsideMustache; 
      std::vector<float> patPhoton_pfMVA; 
      std::vector<float> patPhoton_pfDNN; 
      std::vector<int> patPhoton_egmCutBasedPhotonIDloose;
      std::vector<int> patPhoton_egmCutBasedPhotonIDmedium;
      std::vector<int> patPhoton_egmCutBasedPhotonIDtight;
      std::vector<int> patPhoton_egmMVAPhotonIDmedium;
      std::vector<int> patPhoton_egmMVAPhotonIDtight;
      std::vector<int> patJet_index;
      std::vector<bool> patJet_isCaloJet;
      std::vector<bool> patJet_isJPTJet;
      std::vector<bool> patJet_isPFJet;
      std::vector<bool> patJet_isBasicJet;
      std::vector<float> patJet_charge;
      std::vector<float> patJet_energy;  
      std::vector<float> patJet_eta;    
      std::vector<float> patJet_phi; 
      std::vector<float> patJet_area; 
      std::vector<float> patJet_pt;  
      std::vector<float> patJet_uncorrectedEnergy;  
      std::vector<float> patJet_uncorrectedPt;      
      std::vector<float> patJet_energyFractionHadronic;
      std::vector<float> patJet_hadEnergyInHB;
      std::vector<float> patJet_hadEnergyInHO;
      std::vector<float> patJet_hadEnergyInHE;
      std::vector<float> patJet_hadEnergyInHF;
      std::vector<float> patJet_emEnergyInEB;
      std::vector<float> patJet_emEnergyInEE;
      std::vector<float> patJet_emEnergyInHF;
      std::vector<float> patJet_chargedHadronEnergyFraction;
      std::vector<float> patJet_neutralHadronEnergyFraction;
      std::vector<float> patJet_chargedEmEnergyFraction;
      std::vector<float> patJet_neutralEmEnergyFraction;
      std::vector<int> patJet_nOverlapMuons;
      std::vector<int> patJet_nOverlapTaus;
      std::vector<int> patJet_nOverlapElectrons;
      std::vector<std::vector<int> > patJet_overlapElectronIndices;
      std::vector<int> patJet_nOverlapPhotons;
      std::vector<std::vector<int> > patJet_overlapPhotonIndices;
      std::vector<float> patJet_photonEnergy;
      std::vector<float> patJet_photonEnergyFraction;
      std::vector<float> patJet_electronEnergy;
      std::vector<float> patJet_electronEnergyFraction;
      std::vector<float> patJet_muonEnergy;
      std::vector<float> patJet_muonEnergyFraction;
      std::vector<float> patJet_HFHadronEnergy;       
      std::vector<float> patJet_HFHadronEnergyFraction;
      std::vector<float> patJet_HFEMEnergy;   
      std::vector<float> patJet_HFEMEnergyFraction;
      std::vector<float> patJet_chargedMuEnergy; 
      std::vector<float> patJet_chargedMuEnergyFraction;
      std::vector<float> patJet_hoEnergy;  
      std::vector<float> patJet_hoEnergyFraction;
      std::vector<int> patJet_nCandidates;
      std::vector<int> patJet_nCandInEcal;
      std::vector<int> patJet_nCandInEcalWithCharge;
      std::vector<std::vector<float> > patJet_candInEcal_charge;
      std::vector<std::vector<float> > patJet_candInEcal_ecalEnergy;
      std::vector<std::vector<float> > patJet_candInEcal_ecalEnergyFraction;
      std::vector<std::vector<float> > patJet_candInEcal_hcalEnergy;
      std::vector<std::vector<float> > patJet_candInEcal_hcalEnergyFraction;
      std::vector<std::vector<float> > patJet_candInEcal_eta;
      std::vector<std::vector<float> > patJet_candInEcal_phi;
      std::vector<float> patJet_bTagScore_pfDeepCSV;
      std::vector<float> patJet_puIDScore;
      std::vector<int> patJet_puID; 
      std::vector<float> patJet_qgL;
      std::vector<bool> patJet_jmeCutBasedPFJetIDloose; 
      std::vector<bool> patJet_jmeCutBasedPFJetIDtight;  
      std::vector<bool> patJet_jmeCutBasedPFJetIDtightLepVeto;   
      std::vector<uint32_t> superCluster_seedRawId;
      std::vector<float> superCluster_rawEnergy;
      std::vector<float> superCluster_rawESEnergy;
      std::vector<float> superCluster_energy;
      std::vector<float> superCluster_eta;
      std::vector<float> superCluster_phi;  
      std::vector<float> superCluster_etaWidth;  
      std::vector<float> superCluster_phiWidth;   
      std::vector<float> superCluster_R; 
      std::vector<int> superCluster_nPFClusters;    
      std::vector<int> superCluster_ieta;
      std::vector<int> superCluster_iphi;    
      std::vector<int> superCluster_iz;  
      std::vector<int> superCluster_seedIndex;
      std::vector<std::vector<int> > superCluster_pfClustersIndex;
      std::vector<int> superCluster_nXtals;
      std::vector<int> superCluster_dR_genScore_MatchedIndex;
      std::vector<int> superCluster_dR_simScore_MatchedIndex;
      std::vector<int> superCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<int> superCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<int> superCluster_sim_fraction_MatchedIndex;
      std::vector<int> superCluster_recoToSim_fraction_MatchedIndex;
      std::vector<int> superCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<int> superCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<int> superCluster_recoEnergy_sharedXtals_MatchedIndex;  
      std::vector<std::vector<double> > superCluster_dR_genScore;
      std::vector<std::vector<double> > superCluster_dR_simScore;
      std::vector<std::vector<double> > superCluster_sim_nSharedXtals;
      std::vector<std::vector<double> > superCluster_sim_fraction_noHitsFraction;
      std::vector<std::vector<double> > superCluster_sim_fraction;
      std::vector<std::vector<double> > superCluster_recoToSim_fraction;
      std::vector<std::vector<double> > superCluster_recoToSim_fraction_sharedXtals;
      std::vector<std::vector<double> > superCluster_simEnergy_sharedXtals; 
      std::vector<std::vector<double> > superCluster_recoEnergy_sharedXtals; 
      std::vector<double> superCluster_simPU_nSharedXtals;
      std::vector<double> superCluster_recoEnergy_sharedXtalsPU;  
      std::vector<double> superCluster_simEnergy_sharedXtalsPU; 
      std::vector<double> superCluster_recoEnergy_noHitsFraction_sharedXtalsPU;  
      std::vector<double> superCluster_simEnergy_noHitsFraction_sharedXtalsPU; 
      std::vector<double> superCluster_simOOTPU_nSharedXtals;
      std::vector<double> superCluster_recoEnergy_sharedXtalsOOTPU;  
      std::vector<double> superCluster_simEnergy_sharedXtalsOOTPU; 
      std::vector<double> superCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU;  
      std::vector<double> superCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU; 
      std::vector<float> superCluster_e5x5;
      std::vector<float> superCluster_e2x2Ratio;
      std::vector<float> superCluster_e3x3Ratio;
      std::vector<float> superCluster_eMaxRatio;
      std::vector<float> superCluster_e2ndRatio;
      std::vector<float> superCluster_eTopRatio;
      std::vector<float> superCluster_eRightRatio;
      std::vector<float> superCluster_eBottomRatio;
      std::vector<float> superCluster_eLeftRatio;
      std::vector<float> superCluster_e2x5MaxRatio;
      std::vector<float> superCluster_e2x5TopRatio;
      std::vector<float> superCluster_e2x5RightRatio;
      std::vector<float> superCluster_e2x5BottomRatio;
      std::vector<float> superCluster_e2x5LeftRatio;
      std::vector<float> superCluster_swissCross;
      std::vector<float> superCluster_r9;
      std::vector<float> superCluster_sigmaIetaIeta; 
      std::vector<float> superCluster_sigmaIetaIphi; 
      std::vector<float> superCluster_sigmaIphiIphi; 
      std::vector<float> superCluster_full5x5_e5x5;
      std::vector<float> superCluster_full5x5_e2x2Ratio;
      std::vector<float> superCluster_full5x5_e3x3Ratio;
      std::vector<float> superCluster_full5x5_eMaxRatio;
      std::vector<float> superCluster_full5x5_e2ndRatio;
      std::vector<float> superCluster_full5x5_eTopRatio;
      std::vector<float> superCluster_full5x5_eRightRatio;
      std::vector<float> superCluster_full5x5_eBottomRatio;
      std::vector<float> superCluster_full5x5_eLeftRatio;
      std::vector<float> superCluster_full5x5_e2x5MaxRatio;
      std::vector<float> superCluster_full5x5_e2x5TopRatio;
      std::vector<float> superCluster_full5x5_e2x5RightRatio;
      std::vector<float> superCluster_full5x5_e2x5BottomRatio;
      std::vector<float> superCluster_full5x5_e2x5LeftRatio;
      std::vector<float> superCluster_full5x5_swissCross;
      std::vector<float> superCluster_full5x5_r9;
      std::vector<float> superCluster_full5x5_sigmaIetaIeta; 
      std::vector<float> superCluster_full5x5_sigmaIetaIphi; 
      std::vector<float> superCluster_full5x5_sigmaIphiIphi;      
      std::vector<std::vector<float> > superCluster_psCluster_energy;
      std::vector<std::vector<float> > superCluster_psCluster_eta;
      std::vector<std::vector<float> > superCluster_psCluster_phi;
      std::vector<uint32_t> retunedSuperCluster_seedRawId;
      std::vector<float> retunedSuperCluster_rawEnergy;
      std::vector<float> retunedSuperCluster_rawESEnergy;
      std::vector<float> retunedSuperCluster_energy;
      std::vector<float> retunedSuperCluster_eta;
      std::vector<float> retunedSuperCluster_phi;  
      std::vector<float> retunedSuperCluster_etaWidth;  
      std::vector<float> retunedSuperCluster_phiWidth;   
      std::vector<float> retunedSuperCluster_R; 
      std::vector<int> retunedSuperCluster_nPFClusters;    
      std::vector<int> retunedSuperCluster_ieta;
      std::vector<int> retunedSuperCluster_iphi;    
      std::vector<int> retunedSuperCluster_iz;  
      std::vector<int> retunedSuperCluster_seedIndex;
      std::vector<std::vector<int> > retunedSuperCluster_pfClustersIndex;
      std::vector<int> retunedSuperCluster_nXtals; 
      std::vector<int> retunedSuperCluster_dR_genScore_MatchedIndex;
      std::vector<int> retunedSuperCluster_dR_simScore_MatchedIndex;
      std::vector<int> retunedSuperCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<int> retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<int> retunedSuperCluster_sim_fraction_MatchedIndex;
      std::vector<int> retunedSuperCluster_recoToSim_fraction_MatchedIndex;
      std::vector<int> retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<int> retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<int> retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex;  
      std::vector<std::vector<double> > retunedSuperCluster_dR_genScore;
      std::vector<std::vector<double> > retunedSuperCluster_dR_simScore;
      std::vector<std::vector<double> > retunedSuperCluster_sim_nSharedXtals;
      std::vector<std::vector<double> > retunedSuperCluster_sim_fraction_noHitsFraction;
      std::vector<std::vector<double> > retunedSuperCluster_sim_fraction;
      std::vector<std::vector<double> > retunedSuperCluster_recoToSim_fraction;
      std::vector<std::vector<double> > retunedSuperCluster_recoToSim_fraction_sharedXtals;
      std::vector<std::vector<double> > retunedSuperCluster_simEnergy_sharedXtals; 
      std::vector<std::vector<double> > retunedSuperCluster_recoEnergy_sharedXtals; 
      std::vector<double> retunedSuperCluster_simPU_nSharedXtals;
      std::vector<double> retunedSuperCluster_recoEnergy_sharedXtalsPU;  
      std::vector<double> retunedSuperCluster_simEnergy_sharedXtalsPU; 
      std::vector<double> retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU;  
      std::vector<double> retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU; 
      std::vector<double> retunedSuperCluster_simOOTPU_nSharedXtals;
      std::vector<double> retunedSuperCluster_recoEnergy_sharedXtalsOOTPU;  
      std::vector<double> retunedSuperCluster_simEnergy_sharedXtalsOOTPU; 
      std::vector<double> retunedSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU;  
      std::vector<double> retunedSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU; 
      std::vector<float> retunedSuperCluster_e5x5;
      std::vector<float> retunedSuperCluster_e2x2Ratio;
      std::vector<float> retunedSuperCluster_e3x3Ratio;
      std::vector<float> retunedSuperCluster_eMaxRatio;
      std::vector<float> retunedSuperCluster_e2ndRatio;
      std::vector<float> retunedSuperCluster_eTopRatio;
      std::vector<float> retunedSuperCluster_eRightRatio;
      std::vector<float> retunedSuperCluster_eBottomRatio;
      std::vector<float> retunedSuperCluster_eLeftRatio;
      std::vector<float> retunedSuperCluster_e2x5MaxRatio;
      std::vector<float> retunedSuperCluster_e2x5TopRatio;
      std::vector<float> retunedSuperCluster_e2x5RightRatio;
      std::vector<float> retunedSuperCluster_e2x5BottomRatio;
      std::vector<float> retunedSuperCluster_e2x5LeftRatio;
      std::vector<float> retunedSuperCluster_swissCross;
      std::vector<float> retunedSuperCluster_r9;
      std::vector<float> retunedSuperCluster_sigmaIetaIeta; 
      std::vector<float> retunedSuperCluster_sigmaIetaIphi; 
      std::vector<float> retunedSuperCluster_sigmaIphiIphi; 
      std::vector<float> retunedSuperCluster_full5x5_e5x5;
      std::vector<float> retunedSuperCluster_full5x5_e2x2Ratio;
      std::vector<float> retunedSuperCluster_full5x5_e3x3Ratio;
      std::vector<float> retunedSuperCluster_full5x5_eMaxRatio;
      std::vector<float> retunedSuperCluster_full5x5_e2ndRatio;
      std::vector<float> retunedSuperCluster_full5x5_eTopRatio;
      std::vector<float> retunedSuperCluster_full5x5_eRightRatio;
      std::vector<float> retunedSuperCluster_full5x5_eBottomRatio;
      std::vector<float> retunedSuperCluster_full5x5_eLeftRatio;
      std::vector<float> retunedSuperCluster_full5x5_e2x5MaxRatio;
      std::vector<float> retunedSuperCluster_full5x5_e2x5TopRatio;
      std::vector<float> retunedSuperCluster_full5x5_e2x5RightRatio;
      std::vector<float> retunedSuperCluster_full5x5_e2x5BottomRatio;
      std::vector<float> retunedSuperCluster_full5x5_e2x5LeftRatio;
      std::vector<float> retunedSuperCluster_full5x5_swissCross;
      std::vector<float> retunedSuperCluster_full5x5_r9;
      std::vector<float> retunedSuperCluster_full5x5_sigmaIetaIeta; 
      std::vector<float> retunedSuperCluster_full5x5_sigmaIetaIphi; 
      std::vector<float> retunedSuperCluster_full5x5_sigmaIphiIphi; 
      std::vector<std::vector<float> > retunedSuperCluster_psCluster_energy;
      std::vector<std::vector<float> > retunedSuperCluster_psCluster_eta;
      std::vector<std::vector<float> > retunedSuperCluster_psCluster_phi;
      std::vector<uint32_t> deepSuperCluster_seedRawId;
      std::vector<float> deepSuperCluster_rawEnergy;
      std::vector<float> deepSuperCluster_rawESEnergy;
      std::vector<float> deepSuperCluster_energy; 
      std::vector<float> deepSuperCluster_eta;
      std::vector<float> deepSuperCluster_phi;  
      std::vector<float> deepSuperCluster_etaWidth;  
      std::vector<float> deepSuperCluster_phiWidth;   
      std::vector<float> deepSuperCluster_R; 
      std::vector<int> deepSuperCluster_nPFClusters;    
      std::vector<int> deepSuperCluster_ieta;
      std::vector<int> deepSuperCluster_iphi;    
      std::vector<int> deepSuperCluster_iz;  
      std::vector<int> deepSuperCluster_seedIndex;
      std::vector<std::vector<int> > deepSuperCluster_pfClustersIndex;
      std::vector<int> deepSuperCluster_nXtals; 
      std::vector<int> deepSuperCluster_dR_genScore_MatchedIndex;
      std::vector<int> deepSuperCluster_dR_simScore_MatchedIndex;
      std::vector<int> deepSuperCluster_sim_nSharedXtals_MatchedIndex;
      std::vector<int> deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex;
      std::vector<int> deepSuperCluster_sim_fraction_MatchedIndex;
      std::vector<int> deepSuperCluster_recoToSim_fraction_MatchedIndex;
      std::vector<int> deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex;  
      std::vector<int> deepSuperCluster_simEnergy_sharedXtals_MatchedIndex; 
      std::vector<int> deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex;  
      std::vector<std::vector<double> > deepSuperCluster_dR_genScore;
      std::vector<std::vector<double> > deepSuperCluster_dR_simScore;
      std::vector<std::vector<double> > deepSuperCluster_sim_nSharedXtals;
      std::vector<std::vector<double> > deepSuperCluster_sim_fraction_noHitsFraction;
      std::vector<std::vector<double> > deepSuperCluster_sim_fraction;
      std::vector<std::vector<double> > deepSuperCluster_recoToSim_fraction;
      std::vector<std::vector<double> > deepSuperCluster_recoToSim_fraction_sharedXtals;
      std::vector<std::vector<double> > deepSuperCluster_simEnergy_sharedXtals; 
      std::vector<std::vector<double> > deepSuperCluster_recoEnergy_sharedXtals; 
      std::vector<double> deepSuperCluster_simPU_nSharedXtals;
      std::vector<double> deepSuperCluster_recoEnergy_sharedXtalsPU;  
      std::vector<double> deepSuperCluster_simEnergy_sharedXtalsPU; 
      std::vector<double> deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsPU;  
      std::vector<double> deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsPU; 
      std::vector<double> deepSuperCluster_simOOTPU_nSharedXtals;
      std::vector<double> deepSuperCluster_recoEnergy_sharedXtalsOOTPU;  
      std::vector<double> deepSuperCluster_simEnergy_sharedXtalsOOTPU; 
      std::vector<double> deepSuperCluster_recoEnergy_noHitsFraction_sharedXtalsOOTPU;  
      std::vector<double> deepSuperCluster_simEnergy_noHitsFraction_sharedXtalsOOTPU; 
      std::vector<float> deepSuperCluster_e5x5;
      std::vector<float> deepSuperCluster_e2x2Ratio;
      std::vector<float> deepSuperCluster_e3x3Ratio;
      std::vector<float> deepSuperCluster_eMaxRatio;
      std::vector<float> deepSuperCluster_e2ndRatio;
      std::vector<float> deepSuperCluster_eTopRatio;
      std::vector<float> deepSuperCluster_eRightRatio;
      std::vector<float> deepSuperCluster_eBottomRatio;
      std::vector<float> deepSuperCluster_eLeftRatio;
      std::vector<float> deepSuperCluster_e2x5MaxRatio;
      std::vector<float> deepSuperCluster_e2x5TopRatio;
      std::vector<float> deepSuperCluster_e2x5RightRatio;
      std::vector<float> deepSuperCluster_e2x5BottomRatio;
      std::vector<float> deepSuperCluster_e2x5LeftRatio;
      std::vector<float> deepSuperCluster_swissCross;
      std::vector<float> deepSuperCluster_r9;
      std::vector<float> deepSuperCluster_sigmaIetaIeta; 
      std::vector<float> deepSuperCluster_sigmaIetaIphi; 
      std::vector<float> deepSuperCluster_sigmaIphiIphi; 
      std::vector<float> deepSuperCluster_full5x5_e5x5;
      std::vector<float> deepSuperCluster_full5x5_e2x2Ratio;
      std::vector<float> deepSuperCluster_full5x5_e3x3Ratio;
      std::vector<float> deepSuperCluster_full5x5_eMaxRatio;
      std::vector<float> deepSuperCluster_full5x5_e2ndRatio;
      std::vector<float> deepSuperCluster_full5x5_eTopRatio;
      std::vector<float> deepSuperCluster_full5x5_eRightRatio;
      std::vector<float> deepSuperCluster_full5x5_eBottomRatio;
      std::vector<float> deepSuperCluster_full5x5_eLeftRatio;
      std::vector<float> deepSuperCluster_full5x5_e2x5MaxRatio;
      std::vector<float> deepSuperCluster_full5x5_e2x5TopRatio;
      std::vector<float> deepSuperCluster_full5x5_e2x5RightRatio;
      std::vector<float> deepSuperCluster_full5x5_e2x5BottomRatio;
      std::vector<float> deepSuperCluster_full5x5_e2x5LeftRatio;
      std::vector<float> deepSuperCluster_full5x5_swissCross;
      std::vector<float> deepSuperCluster_full5x5_r9;
      std::vector<float> deepSuperCluster_full5x5_sigmaIetaIeta; 
      std::vector<float> deepSuperCluster_full5x5_sigmaIetaIphi; 
      std::vector<float> deepSuperCluster_full5x5_sigmaIphiIphi; 
      std::vector<std::vector<float> > deepSuperCluster_psCluster_energy;
      std::vector<std::vector<float> > deepSuperCluster_psCluster_eta;
      std::vector<std::vector<float> > deepSuperCluster_psCluster_phi;
      std::vector<uint32_t> gsfgsfElectron_seedRawId;
      std::vector<uint32_t> gsfgedPhoton_seedRawId;
      std::vector<double> dR_genScore;
      std::vector<double> dR_simScore;
      std::vector<double> sim_nSharedXtals;
      std::vector<double> sim_fraction_noHitsFraction;
      std::vector<double> sim_fraction;
      std::vector<double> recoToSim_fraction;
      std::vector<double> recoToSim_fraction_sharedXtals;  
      std::vector<double> simEnergy_sharedXtals; 
      std::vector<double> recoEnergy_sharedXtals; 
      std::vector<double> simPU_nSharedXtals;
      std::vector<double> simEnergy_sharedXtalsPU; 
      std::vector<double> recoEnergy_sharedXtalsPU; 
      std::vector<double> simEnergy_noHitsFraction_sharedXtalsPU; 
      std::vector<double> recoEnergy_noHitsFraction_sharedXtalsPU; 
      std::vector<double> simOOTPU_nSharedXtals;
      std::vector<double> simEnergy_sharedXtalsOOTPU; 
      std::vector<double> recoEnergy_sharedXtalsOOTPU; 
      std::vector<double> simEnergy_noHitsFraction_sharedXtalsOOTPU; 
      std::vector<double> recoEnergy_noHitsFraction_sharedXtalsOOTPU; 

      std::vector<reco::GenParticle> genParts;  
      std::map<int,int> genMapMother;
      std::map<int,std::vector<int> > genMapDaughters; 
      std::vector<DetId> pfRechit_unClustered;
      std::vector<std::vector<DetId>> hits_PFCluster;
      std::vector<std::vector<DetId>> hits_CaloParticle;
      std::vector<std::vector<float>> energies_CaloParticle;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart; 
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPartPU; 
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPartOOTPU; 
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_PFCluster;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_SuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_SuperClusterEE;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_RetunedSuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_RetunedSuperClusterEE;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_DeepSuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_DeepSuperClusterEE;   
  
};

#endif
