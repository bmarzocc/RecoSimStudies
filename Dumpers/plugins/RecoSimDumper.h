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
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle,  float simHitEnergy_cut);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE);
      std::vector<double> getScores(const std::vector<std::pair<DetId, float> >*hits_and_energies_Cluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE);
      GlobalPoint calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFRecHit>  > pfRecHitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCluster> > pfClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeSuperClusterToken_; 

      edm::Service<TFileService> iFile;
      const CaloSubdetectorGeometry* _ebGeom;
      const CaloSubdetectorGeometry* _eeGeom;
      const CaloSubdetectorGeometry* _esGeom;
      bool _esPlus;
      bool _esMinus;

      // ----------config inputs-------------------
      bool doCompression_;
      int nBits_;
      bool saveGenParticles_;       
      bool saveCaloParticles_;    
      bool saveSimhits_;          
      bool saveRechits_;          
      bool savePFRechits_;   
      bool savePFCluster_;    
      bool savePFClusterhits_;   
      bool saveSuperCluster_;     
      bool saveShowerShapes_;   
      bool saveScores_;      
      std::vector<int> genID_;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      std::vector<std::map<uint32_t,float> > caloParticleXtals_;
      std::vector<float> locCov;
      std::vector<float> full5x5_locCov;
      
      long int eventId;
      int lumiId;
      int runId; 
      int nVtx;
      std::vector<int> genParticle_id;
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;
      std::vector<std::vector<int> > genParticle_pfCluster_dR_genScore_MatchedIndex;
      std::vector<std::vector<int> > genParticle_superCluster_dR_genScore_MatchedIndex;
      std::vector<int> caloParticle_id;
      std::vector<float> caloParticle_genEnergy;
      std::vector<float> caloParticle_simEnergy;
      std::vector<float> caloParticle_genPt;
      std::vector<float> caloParticle_simPt;
      std::vector<float> caloParticle_genEta;
      std::vector<float> caloParticle_simEta;
      std::vector<float> caloParticle_genPhi;
      std::vector<float> caloParticle_simPhi;
      std::vector<int> caloParticle_simIeta;
      std::vector<int> caloParticle_simIphi;
      std::vector<int> caloParticle_simIz;
      std::vector<std::vector<int> > caloParticle_pfCluster_dR_simScore_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_n_shared_xtals_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_old_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min1_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min3_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min3_1MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min3_5MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min3_10MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min3_50MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_fraction_min3_100MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_rechit_diff_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_dR_simScore_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_n_shared_xtals_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_old_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min1_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min3_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min3_1MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min3_5MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min3_10MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min3_50MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_fraction_min3_100MeVCut_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_rechit_diff_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_sim_rechit_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex;
      std::vector<std::vector<int> > caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex;
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
      std::vector<std::vector<float> > pfClusterHit_energy; 
      std::vector<std::vector<float> > pfClusterHit_rechitEnergy; 
      std::vector<std::vector<float> > pfClusterHit_eta;
      std::vector<std::vector<float> > pfClusterHit_phi;
      std::vector<std::vector<int> > pfClusterHit_ieta;
      std::vector<std::vector<int> > pfClusterHit_iphi;
      std::vector<std::vector<int> > pfClusterHit_iz;
      std::vector<float> pfCluster_energy;
      std::vector<float> pfCluster_eta;
      std::vector<float> pfCluster_phi;
      std::vector<int> pfCluster_ieta;
      std::vector<int> pfCluster_iphi;
      std::vector<int> pfCluster_iz;
      std::vector<std::vector<int> > pfCluster_superClustersIndex;
      std::vector<int> pfCluster_dR_genScore_MatchedIndex;
      std::vector<int> pfCluster_dR_simScore_MatchedIndex;
      std::vector<int> pfCluster_n_shared_xtals_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_old_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min1_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min3_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min3_1MeVCut_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min3_5MeVCut_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min3_10MeVCut_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min3_50MeVCut_MatchedIndex;
      std::vector<int> pfCluster_sim_fraction_min3_100MeVCut_MatchedIndex; 
      std::vector<int> pfCluster_sim_rechit_diff_MatchedIndex;
      std::vector<int> pfCluster_sim_rechit_fraction_MatchedIndex;
      std::vector<int> pfCluster_global_sim_rechit_fraction_MatchedIndex;
      std::vector<int> pfCluster_hgcal_caloToCluster_MatchedIndex;
      std::vector<int> pfCluster_hgcal_clusterToCalo_MatchedIndex;   
      std::vector<std::vector<double> > pfCluster_dR_genScore;
      std::vector<std::vector<double> > pfCluster_dR_simScore;
      std::vector<std::vector<int> > pfCluster_n_shared_xtals;
      std::vector<std::vector<double> > pfCluster_sim_fraction_old;
      std::vector<std::vector<double> > pfCluster_sim_fraction;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min1;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min3;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min3_1MeVCut;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min3_5MeVCut;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min3_10MeVCut;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min3_50MeVCut;
      std::vector<std::vector<double> > pfCluster_sim_fraction_min3_100MeVCut; 
      std::vector<std::vector<double> > pfCluster_sim_rechit_diff;
      std::vector<std::vector<double> > pfCluster_sim_rechit_fraction;
      std::vector<std::vector<double> > pfCluster_global_sim_rechit_fraction;
      std::vector<std::vector<double> > pfCluster_hgcal_caloToCluster;
      std::vector<std::vector<double> > pfCluster_hgcal_clusterToCalo;
      std::vector<float> superCluster_energy;
      std::vector<float> superCluster_eta;
      std::vector<float> superCluster_phi;  
      std::vector<float> superCluster_etaWidth;  
      std::vector<float> superCluster_phiWidth;   
      std::vector<float> superCluster_R;   
      std::vector<int> superCluster_ieta;
      std::vector<int> superCluster_iphi;    
      std::vector<int> superCluster_iz;  
      std::vector<int> superCluster_seedIndex;
      std::vector<std::vector<int> > superCluster_pfClustersIndex;
      std::vector<int> superCluster_dR_genScore_MatchedIndex;
      std::vector<int> superCluster_dR_simScore_MatchedIndex;
      std::vector<int> superCluster_n_shared_xtals_MatchedIndex;
      std::vector<int> superCluster_sim_fraction_old_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_MatchedIndex;
      std::vector<int> superCluster_sim_fraction_min1_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_min3_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_min3_1MeVCut_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_min3_5MeVCut_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_min3_10MeVCut_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_min3_50MeVCut_MatchedIndex; 
      std::vector<int> superCluster_sim_fraction_min3_100MeVCut_MatchedIndex; 
      std::vector<int> superCluster_sim_rechit_diff_MatchedIndex;
      std::vector<int> superCluster_sim_rechit_fraction_MatchedIndex;
      std::vector<int> superCluster_global_sim_rechit_fraction_MatchedIndex;
      std::vector<int> superCluster_hgcal_caloToCluster_MatchedIndex;
      std::vector<int> superCluster_hgcal_clusterToCalo_MatchedIndex;
      std::vector<std::vector<double> > superCluster_dR_genScore;
      std::vector<std::vector<double> > superCluster_dR_simScore;
      std::vector<std::vector<int> > superCluster_n_shared_xtals;
      std::vector<std::vector<double> > superCluster_sim_fraction_old;
      std::vector<std::vector<double> > superCluster_sim_fraction;
      std::vector<std::vector<double> > superCluster_sim_fraction_min1;
      std::vector<std::vector<double> > superCluster_sim_fraction_min3;
      std::vector<std::vector<double> > superCluster_sim_fraction_min3_1MeVCut;
      std::vector<std::vector<double> > superCluster_sim_fraction_min3_5MeVCut;
      std::vector<std::vector<double> > superCluster_sim_fraction_min3_10MeVCut;
      std::vector<std::vector<double> > superCluster_sim_fraction_min3_50MeVCut;
      std::vector<std::vector<double> > superCluster_sim_fraction_min3_100MeVCut;
      std::vector<std::vector<double> > superCluster_sim_rechit_diff;
      std::vector<std::vector<double> > superCluster_sim_rechit_fraction;
      std::vector<std::vector<double> > superCluster_global_sim_rechit_fraction;
      std::vector<std::vector<double> > superCluster_hgcal_caloToCluster;
      std::vector<std::vector<double> > superCluster_hgcal_clusterToCalo;
      std::vector<float> superCluster_r9;
      std::vector<float> superCluster_sigmaIetaIeta; 
      std::vector<float> superCluster_sigmaIetaIphi; 
      std::vector<float> superCluster_sigmaIphiIphi; 
      std::vector<float> superCluster_full5x5_r9; 
      std::vector<float> superCluster_full5x5_sigmaIetaIeta;
      std::vector<float> superCluster_full5x5_sigmaIetaIphi;
      std::vector<float> superCluster_full5x5_sigmaIphiIphi; 
      std::vector<std::vector<float> > psCluster_energy;
      std::vector<std::vector<float> > psCluster_eta;
      std::vector<std::vector<float> > psCluster_phi;

      std::vector<double> dR_genScore;
      std::vector<double> dR_simScore;
      std::vector<int> n_shared_xtals;
      std::vector<double> sim_fraction_old;
      std::vector<double> sim_fraction;
      std::vector<double> sim_fraction_min1;
      std::vector<double> sim_fraction_min3;
      std::vector<double> sim_fraction_min3_1MeVCut;
      std::vector<double> sim_fraction_min3_5MeVCut;
      std::vector<double> sim_fraction_min3_10MeVCut;
      std::vector<double> sim_fraction_min3_50MeVCut;
      std::vector<double> sim_fraction_min3_100MeVCut;
      std::vector<double> sim_rechit_diff;
      std::vector<double> sim_rechit_fraction;
      std::vector<double> global_sim_rechit_fraction;
      std::vector<double> hgcal_caloToCluster;
      std::vector<double> hgcal_clusterToCalo;
      std::vector<DetId> pfRechit_unClustered;
      std::vector<std::vector<DetId>> hits_CaloPart;
      std::vector<std::vector<DetId>> hits_PFCluster;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart; 
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart_1MeVCut;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart_5MeVCut;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart_10MeVCut;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart_50MeVCut;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart_100MeVCut;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_PFCluster;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_SuperClusterEB;
      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_SuperClusterEE;
};

#endif
