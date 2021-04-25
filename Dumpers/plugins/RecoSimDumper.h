#ifndef RecoSimStudies_Dumpers_RecoSimDumper_H
#define RecoSimStudies_Dumpers_RecoSimDumper_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/LooperFactory.h"
#include "FWCore/Framework/interface/ESProducerLooper.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/Event.h"

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"

#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
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
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
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

class RecoSimDumper : public edm::EDAnalyzer
{
      public:
         explicit RecoSimDumper(const edm::ParameterSet&);
	 ~RecoSimDumper();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
         virtual void endJob() ;
        
      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      void setTree(TTree* tree);
      void setVectors(int nGenParticles, int nCaloParticles, int nPFClusters, int nSuperClustersEB, int nSuperClustersEE, int nRetunedSuperClustersEB, int nRetunedSuperClustersEE, int nDeepSuperClustersEB, int nDeepSuperClustersEE); 
      double ptFast(const double energy, const math::XYZPoint& position, const math::XYZPoint& origin);
      void addDaughters(const std::vector<reco::GenParticle>* genParticles, const int genIndex1, const int genIndex2);
      int getGenMother(const std::vector<reco::GenParticle>* genParticles, const int genIndex);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<std::pair<DetId,std::pair<float,float>>>* getSharedHitsAndEnergies(const std::vector<std::pair<DetId, float> >* hitsAndEnergies1, const std::vector<std::pair<DetId, float> >* hitsAndEnergies2);
      std::pair<double,double> calculateCovariances(const reco::PFCluster* pfCluster, const EcalRecHitCollection* recHits, const CaloSubdetectorGeometry* geometry);
      std::vector<float> getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology);
      std::vector<double> getScores(const reco::PFCluster* pfCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      std::vector<double> getScores(const reco::SuperCluster* superCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE);
      int getMatchedIndex(std::vector<std::vector<double>>* score, double selection, bool useMax, double scale, int iCl);
      void fillParticleMatchedIndex(std::vector<std::vector<int>>* particleMatchedIndex, std::vector<int>* clusterMatchedIndex);
      GlobalPoint calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<std::vector<SimVertex> > simVtxToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFRecHit>  > pfRecHitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCluster> > pfClusterToken_; 
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
      edm::EDGetTokenT<CaloTowerCollection> hcalTowersToken_;   

      edm::Service<TFileService> iFile;
      const CaloSubdetectorGeometry* _ebGeom;
      const CaloSubdetectorGeometry* _eeGeom;
      const CaloSubdetectorGeometry* _esGeom;
      bool _esPlus;
      bool _esMinus;
      
      // ----------config inputs-------------------
      bool useRetunedSC_;
      bool useDeepSC_;
      bool doCompression_;
      int nBits_;
      bool saveGenParticles_;       
      bool saveCaloParticles_;  
      bool saveCaloParticlesPU_; 
      bool saveCaloParticlesOOTPU_;    
      bool saveSimhits_;          
      bool saveRechits_;          
      bool savePFRechits_;   
      bool savePFCluster_;    
      bool savePFClusterhits_;   
      bool saveSuperCluster_;     
      bool saveShowerShapes_;   
      
      // ----------DNN inputs-------------------
      std::vector<double> HLF_VectorVar_;
      std::vector<std::vector<double>> PL_VectorVar_;
      std::vector<double> x_mean_, x_std_, list_mean_, list_std_;
     
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      std::vector<std::map<uint32_t,float> > caloParticleXtals_;
      std::pair<double,double> widths_;
      std::vector<float> locCov_;
      std::vector<float> full5x5_locCov_;
      std::vector<float> showerShapes_;
      std::map<int,std::vector<int>> genDaughters;
      
      long int eventId;
      int lumiId;
      int runId; 
      double truePU;
      double obsPU;
      int nVtx;
      float rho; 
      int genParticle_size; 
      int caloParticle_size;
      int caloParticlePU_size;
      int caloParticleOOTPU_size; 
      std::vector<int> genParticle_pdgId;
      std::vector<int> genParticle_status; 
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;
      std::vector<std::vector<int> > genParticle_pfCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_superCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_retunedSuperCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_deepSuperCluster_dR_genScore_MatchedIndex; 
      std::vector<int> caloParticle_index; 
      std::vector<int> caloParticle_nXtals;  
      std::vector<bool> caloParticle_isPU;
      std::vector<bool> caloParticle_isOOTPU;
      std::vector<int> caloParticle_pdgId;
      std::vector<int> caloParticle_status;
      std::vector<int> caloParticle_g4TracksEventID;
      std::vector<int> caloParticle_g4TracksBX;
      std::vector<float> caloParticle_genEnergy;
      std::vector<float> caloParticle_simEnergy;
      std::vector<float> caloParticle_genPt;
      std::vector<float> caloParticle_simPt;
      std::vector<float> caloParticle_genEta;
      std::vector<float> caloParticle_simEta;
      std::vector<float> caloParticle_genPhi;
      std::vector<float> caloParticle_simPhi;
      std::vector<int> caloParticle_genMotherIndex;
      std::vector<int> caloParticle_genMotherPdgId; 
      std::vector<float> caloParticle_genMotherEnergy;
      std::vector<float> caloParticle_genMotherPt;
      std::vector<float> caloParticle_genMotherEta;
      std::vector<float> caloParticle_genMotherPhi; 
      std::vector<float> caloParticle_genMotherDR; 
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
      std::vector<float> pfCluster_rawEnergy;
      std::vector<float> pfCluster_energy;
      std::vector<float> pfCluster_rawPt;
      std::vector<float> pfCluster_pt;
      std::vector<float> pfCluster_eta;
      std::vector<float> pfCluster_phi; 
      std::vector<int> pfCluster_ieta;
      std::vector<int> pfCluster_iphi;
      std::vector<int> pfCluster_iz;
      std::vector<int> pfCluster_nXtals;
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
      std::vector<float> superCluster_rawEnergy;
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
      std::vector<float> retunedSuperCluster_rawEnergy;
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
      std::vector<float> deepSuperCluster_rawEnergy;
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
      std::vector<double> dR_genScore;
      std::vector<double> dR_simScore;
      std::vector<double> sim_nSharedXtals;
      std::vector<double> sim_fraction_noHitsFraction;
      std::vector<double> sim_fraction;
      std::vector<double> recoToSim_fraction;
      std::vector<double> recoToSim_fraction_sharedXtals;  
      std::vector<double> simEnergy_sharedXtals; 
      std::vector<double> recoEnergy_sharedXtals; 
      std::vector<DetId> pfRechit_unClustered;
      std::vector<std::vector<DetId>> hits_PFCluster;
      std::vector<std::vector<DetId>> hits_CaloParticle;
      std::vector<std::vector<float>> energies_CaloParticle;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart; 
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_PFCluster;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_SuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_SuperClusterEE;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_RetunedSuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_RetunedSuperClusterEE;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_DeepSuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_DeepSuperClusterEE;   
  
};

#endif
