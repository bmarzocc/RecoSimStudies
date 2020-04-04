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
#include "RecoSimStudies/Dumpers/plugins/RecoSimDumper.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

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
#include "TCanvas.h"
#include "TMath.h"
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

using namespace cms;
using namespace edm;
using namespace std;
using namespace reco;

//
// constructors and destructor
//
RecoSimDumper::RecoSimDumper(const edm::ParameterSet& iConfig)
{

   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_                 = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   pfRecHitToken_                 = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection")); 
   pfClusterToken_                = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   ebSuperClusterToken_           = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_           = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   useHcalTowers_                 = iConfig.getParameter<bool>("useHcalTowers");  
   hcalTowersToken_               = consumes<CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("hcalTowersCollection"));
   useRetunedSC_                  = iConfig.getParameter<bool>("useRetunedSC");  
   if(useRetunedSC_){
      ebRetunedSuperClusterToken_ = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebRetunedSuperClusterCollection"));
      eeRetunedSuperClusterToken_ = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeRetunedSuperClusterCollection"));
   } 
   useDeepSC_                     = iConfig.getParameter<bool>("useDeepSC");  
   if(useDeepSC_){
      ebDeepSuperClusterToken_    = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebDeepSuperClusterCollection"));
      eeDeepSuperClusterToken_    = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeDeepSuperClusterCollection"));
      ebDeepSuperClusterLWPToken_ = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebDeepSuperClusterLWPCollection"));
      eeDeepSuperClusterLWPToken_ = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeDeepSuperClusterLWPCollection"));
      ebDeepSuperClusterTWPToken_ = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebDeepSuperClusterTWPCollection"));
      eeDeepSuperClusterTWPToken_ = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeDeepSuperClusterTWPCollection"));
   }
   doCompression_                 = iConfig.getParameter<bool>("doCompression");
   nBits_                         = iConfig.getParameter<int>("nBits");
   saveGenParticles_              = iConfig.getParameter<bool>("saveGenParticles");
   saveCaloParticles_             = iConfig.getParameter<bool>("saveCaloParticles");
   saveSimhits_             	  = iConfig.getParameter<bool>("saveSimhits");
   saveRechits_                   = iConfig.getParameter<bool>("saveRechits");
   savePFRechits_                 = iConfig.getParameter<bool>("savePFRechits"); 
   savePFCluster_                 = iConfig.getParameter<bool>("savePFCluster");
   savePFClusterhits_             = iConfig.getParameter<bool>("savePFClusterhits");
   saveSuperCluster_              = iConfig.getParameter<bool>("saveSuperCluster");
   saveShowerShapes_              = iConfig.getParameter<bool>("saveShowerShapes");

   scoreType_                     = iConfig.getParameter<std::string>("scoreType");
   saveScores_                    = iConfig.getParameter<bool>("saveScores");
   genID_                         = iConfig.getParameter<std::vector<int>>("genID");
   
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   tree->Branch("rho", &rho, "rho/F"); 
   if(saveGenParticles_){
      tree->Branch("genParticle_id","std::vector<int>",&genParticle_id);
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
      if(savePFCluster_) tree->Branch("genParticle_pfCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_pfCluster_dR_genScore_MatchedIndex);
      if(saveSuperCluster_) tree->Branch("genParticle_superCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_superCluster_dR_genScore_MatchedIndex);
      if(saveSuperCluster_ && useRetunedSC_) tree->Branch("genParticle_retunedSuperCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_retunedSuperCluster_dR_genScore_MatchedIndex); 
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("genParticle_deepSuperCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_deepSuperCluster_dR_genScore_MatchedIndex); 
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("genParticle_deepSuperClusterLWP_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_deepSuperClusterLWP_dR_genScore_MatchedIndex);  
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("genParticle_deepSuperClusterTWP_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_deepSuperClusterTWP_dR_genScore_MatchedIndex);  
   }
   if(saveCaloParticles_){
      tree->Branch("caloParticle_id","std::vector<int>",&caloParticle_id); 
      tree->Branch("caloParticle_genEnergy","std::vector<float>",&caloParticle_genEnergy);
      tree->Branch("caloParticle_simEnergy","std::vector<float>",&caloParticle_simEnergy); 
      tree->Branch("caloParticle_genPt","std::vector<float>",&caloParticle_genPt);
      tree->Branch("caloParticle_simPt","std::vector<float>",&caloParticle_simPt);
      tree->Branch("caloParticle_genEta","std::vector<float>",&caloParticle_genEta);
      tree->Branch("caloParticle_simEta","std::vector<float>",&caloParticle_simEta);
      tree->Branch("caloParticle_genPhi","std::vector<float>",&caloParticle_genPhi);
      tree->Branch("caloParticle_simPhi","std::vector<float>",&caloParticle_simPhi);
      tree->Branch("caloParticle_simIeta","std::vector<int>",&caloParticle_simIeta);
      tree->Branch("caloParticle_simIphi","std::vector<int>",&caloParticle_simIphi);
      tree->Branch("caloParticle_simIz","std::vector<int>",&caloParticle_simIz);
      if(savePFCluster_){   
         tree->Branch("caloParticle_pfCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_old_MatchedIndex);
         if(!saveScores_) tree->Branch("caloParticle_pfCluster_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_simScore_MatchedIndex);
         if(saveScores_){
            tree->Branch("caloParticle_pfCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_n_shared_xtals_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_diff_MatchedIndex);
            tree->Branch("caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex);   
            tree->Branch("caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex); 
            tree->Branch("caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex);  
            tree->Branch("caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex); 
            tree->Branch("caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex);  
            tree->Branch("caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex);        
         }   
      }
      if(saveSuperCluster_){    
         tree->Branch("caloParticle_superCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_old_MatchedIndex);
         if(!saveScores_) tree->Branch("caloParticle_superCluster_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_simScore_MatchedIndex);
         if(useRetunedSC_){
            tree->Branch("caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex);
            if(!saveScores_) tree->Branch("caloParticle_retunedSuperCluster_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_simScore_MatchedIndex);
         } 
         if(saveScores_){
            tree->Branch("caloParticle_superCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_n_shared_xtals_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_1MeVCut_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_5MeVCut_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_10MeVCut_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_50MeVCut_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_100MeVCut_MatchedIndex);  
            tree->Branch("caloParticle_superCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_500MeVCut_MatchedIndex);  
            tree->Branch("caloParticle_superCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_1GeVCut_MatchedIndex);    
            tree->Branch("caloParticle_superCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_rechit_diff_MatchedIndex);
            tree->Branch("caloParticle_superCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_rechit_fraction_MatchedIndex);   
            tree->Branch("caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex); 
            tree->Branch("caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex); 
            tree->Branch("caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex); 
            tree->Branch("caloParticle_superCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_rechit_combined_fraction_MatchedIndex); 
            tree->Branch("caloParticle_superCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_rechit_sim_combined_fraction_MatchedIndex);  
            if(useRetunedSC_){
               tree->Branch("caloParticle_retunedSuperCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_n_shared_xtals_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex);    
               tree->Branch("caloParticle_retunedSuperCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_rechit_diff_MatchedIndex);
               tree->Branch("caloParticle_retunedSuperCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("caloParticle_retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("caloParticle_retunedSuperCluster_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("caloParticle_retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("caloParticle_retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("caloParticle_retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex); 
            }
         }   
         if(useDeepSC_){
            tree->Branch("caloParticle_deepSuperCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperCluster_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_old_MatchedIndex);
            if(!saveScores_) tree->Branch("caloParticle_deepSuperCluster_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperClusterLWP_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_old_MatchedIndex);
            if(!saveScores_) tree->Branch("caloParticle_deepSuperClusterLWP_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperClusterTWP_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_old_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_old_MatchedIndex);
            if(!saveScores_) tree->Branch("caloParticle_deepSuperClusterTWP_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_simScore_MatchedIndex); 
            if(saveScores_){
               tree->Branch("caloParticle_deepSuperCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_n_shared_xtals_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex);    
               tree->Branch("caloParticle_deepSuperCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_rechit_diff_MatchedIndex);
               tree->Branch("caloParticle_deepSuperCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("caloParticle_deepSuperCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperCluster_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_n_shared_xtals_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex);    
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_rechit_diff_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("caloParticle_deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_n_shared_xtals_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex);  
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex);    
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_rechit_diff_MatchedIndex);
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("caloParticle_deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("caloParticle_deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex); 
            }  
         } 
      }
   }
   if(saveSimhits_){
      tree->Branch("simHit_energy","std::vector<std::vector<float> >",&simHit_energy);
      tree->Branch("simHit_eta","std::vector<std::vector<float> >",&simHit_eta);
      tree->Branch("simHit_phi","std::vector<std::vector<float> >",&simHit_phi);
      tree->Branch("simHit_ieta","std::vector<std::vector<int> >",&simHit_ieta);
      tree->Branch("simHit_iphi","std::vector<std::vector<int> >",&simHit_iphi);
      tree->Branch("simHit_iz","std::vector<std::vector<int> >",&simHit_iz);
   }
   if(saveRechits_){
      tree->Branch("recHit_noPF_energy","std::vector<float>",&recHit_noPF_energy);
      tree->Branch("recHit_noPF_eta","std::vector<float>",&recHit_noPF_eta); 
      tree->Branch("recHit_noPF_phi","std::vector<float>",&recHit_noPF_phi);
      tree->Branch("recHit_noPF_ieta","std::vector<int>",&recHit_noPF_ieta); 
      tree->Branch("recHit_noPF_iphi","std::vector<int>",&recHit_noPF_iphi);
      tree->Branch("recHit_noPF_iz","std::vector<int>",&recHit_noPF_iz);     
   }
   if(savePFRechits_){ 
      tree->Branch("pfRecHit_unClustered_energy","std::vector<float>",&pfRecHit_unClustered_energy);
      tree->Branch("pfRecHit_unClustered_eta","std::vector<float>",&pfRecHit_unClustered_eta); 
      tree->Branch("pfRecHit_unClustered_phi","std::vector<float>",&pfRecHit_unClustered_phi);
      tree->Branch("pfRecHit_unClustered_ieta","std::vector<int>",&pfRecHit_unClustered_ieta); 
      tree->Branch("pfRecHit_unClustered_iphi","std::vector<int>",&pfRecHit_unClustered_iphi);
      tree->Branch("pfRecHit_unClustered_iz","std::vector<int>",&pfRecHit_unClustered_iz);     
   }
   if(savePFCluster_){
      tree->Branch("pfCluster_energy","std::vector<float>",&pfCluster_energy);
      tree->Branch("pfCluster_eta","std::vector<float>",&pfCluster_eta);
      tree->Branch("pfCluster_phi","std::vector<float>",&pfCluster_phi);   
      tree->Branch("pfCluster_ieta","std::vector<int>",&pfCluster_ieta);
      tree->Branch("pfCluster_iphi","std::vector<int>",&pfCluster_iphi);   
      tree->Branch("pfCluster_iz","std::vector<int>",&pfCluster_iz);
      tree->Branch("pfCluster_nXtals","std::vector<int>",&pfCluster_nXtals);  
      if(saveSuperCluster_) tree->Branch("pfCluster_superClustersIndex","std::vector<std::vector<int> >",&pfCluster_superClustersIndex); 
      if(saveSuperCluster_ && useRetunedSC_) tree->Branch("pfCluster_retunedSuperClustersIndex","std::vector<std::vector<int> >",&pfCluster_retunedSuperClustersIndex);  
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("pfCluster_deepSuperClustersIndex","std::vector<std::vector<int> >",&pfCluster_deepSuperClustersIndex); 
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("pfCluster_deepSuperClusterLWPsIndex","std::vector<std::vector<int> >",&pfCluster_deepSuperClusterLWPsIndex); 
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("pfCluster_deepSuperClusterTWPsIndex","std::vector<std::vector<int> >",&pfCluster_deepSuperClusterTWPsIndex); 
      if(saveCaloParticles_){ 
         tree->Branch("pfCluster_dR_genScore_MatchedIndex","std::vector<int>",&pfCluster_dR_genScore_MatchedIndex);
         tree->Branch("pfCluster_dR_simScore_MatchedIndex","std::vector<int>",&pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_old_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_old_MatchedIndex);
         if(!saveScores_) tree->Branch("pfCluster_simScore_MatchedIndex","std::vector<int>",&pfCluster_simScore_MatchedIndex);
         if(saveScores_){
            tree->Branch("pfCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&pfCluster_n_shared_xtals_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_1MeVCut_MatchedIndex); 
            tree->Branch("pfCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_5MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_10MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_50MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_100MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_500MeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_1GeVCut_MatchedIndex);
            tree->Branch("pfCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_diff_MatchedIndex);
            tree->Branch("pfCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_fraction_MatchedIndex);   
            tree->Branch("pfCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&pfCluster_global_sim_rechit_fraction_MatchedIndex); 
            tree->Branch("pfCluster_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&pfCluster_hgcal_caloToCluster_MatchedIndex); 
            tree->Branch("pfCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&pfCluster_hgcal_clusterToCalo_MatchedIndex);   
            tree->Branch("pfCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_combined_fraction_MatchedIndex); 
            tree->Branch("pfCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&pfCluster_rechit_sim_combined_fraction_MatchedIndex);     
         } 
      } 
      if(saveCaloParticles_){
         tree->Branch("pfCluster_dR_genScore","std::vector<std::vector<double> >",&pfCluster_dR_genScore);
         tree->Branch("pfCluster_dR_simScore","std::vector<std::vector<double> >",&pfCluster_dR_simScore);
         tree->Branch("pfCluster_sim_fraction_old","std::vector<std::vector<double> >",&pfCluster_sim_fraction_old);
         if(!saveScores_) tree->Branch("pfCluster_simScore","std::vector<std::vector<double> >",&pfCluster_simScore);
         if(saveScores_){
            tree->Branch("pfCluster_n_shared_xtals","std::vector<std::vector<double> >",&pfCluster_n_shared_xtals);
            tree->Branch("pfCluster_sim_fraction","std::vector<std::vector<double> >",&pfCluster_sim_fraction);
            tree->Branch("pfCluster_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_1MeVCut);
            tree->Branch("pfCluster_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_5MeVCut);
            tree->Branch("pfCluster_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_10MeVCut);
            tree->Branch("pfCluster_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_50MeVCut);
            tree->Branch("pfCluster_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_100MeVCut); 
            tree->Branch("pfCluster_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_500MeVCut); 
            tree->Branch("pfCluster_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&pfCluster_sim_fraction_1GeVCut);      
            tree->Branch("pfCluster_sim_rechit_diff","std::vector<std::vector<double> >",&pfCluster_sim_rechit_diff);
            tree->Branch("pfCluster_sim_rechit_fraction","std::vector<std::vector<double> >",&pfCluster_sim_rechit_fraction);   
            tree->Branch("pfCluster_global_sim_rechit_fraction","std::vector<std::vector<double> >",&pfCluster_global_sim_rechit_fraction); 
            tree->Branch("pfCluster_hgcal_caloToCluster","std::vector<std::vector<double> >",&pfCluster_hgcal_caloToCluster); 
            tree->Branch("pfCluster_hgcal_clusterToCalo","std::vector<std::vector<double> >",&pfCluster_hgcal_clusterToCalo); 
            tree->Branch("pfCluster_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&pfCluster_sim_rechit_combined_fraction); 
            tree->Branch("pfCluster_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&pfCluster_rechit_sim_combined_fraction);   
         }
      }
      if(savePFClusterhits_){ 
         tree->Branch("pfClusterHit_energy","std::vector<std::vector<float> >",&pfClusterHit_energy);
         tree->Branch("pfClusterHit_rechitEnergy","std::vector<std::vector<float> >",&pfClusterHit_rechitEnergy);
         tree->Branch("pfClusterHit_eta","std::vector<std::vector<float> >",&pfClusterHit_eta);
         tree->Branch("pfClusterHit_phi","std::vector<std::vector<float> >",&pfClusterHit_phi);
         tree->Branch("pfClusterHit_ieta","std::vector<std::vector<int> >",&pfClusterHit_ieta);
         tree->Branch("pfClusterHit_iphi","std::vector<std::vector<int> >",&pfClusterHit_iphi);
         tree->Branch("pfClusterHit_iz","std::vector<std::vector<int> >",&pfClusterHit_iz);
      }   
   }
   if(saveSuperCluster_){
      tree->Branch("superCluster_energy","std::vector<float> ",&superCluster_energy);
      tree->Branch("superCluster_eta","std::vector<float>",&superCluster_eta);
      tree->Branch("superCluster_phi","std::vector<float>",&superCluster_phi);  
      tree->Branch("superCluster_etaWidth","std::vector<float>",&superCluster_etaWidth);
      tree->Branch("superCluster_phiWidth","std::vector<float>",&superCluster_phiWidth);  
      tree->Branch("superCluster_R","std::vector<float>",&superCluster_R);   
      tree->Branch("superCluster_nPFClusters","std::vector<int>",&superCluster_nPFClusters);   
      tree->Branch("superCluster_ieta","std::vector<int>",&superCluster_ieta);
      tree->Branch("superCluster_iphi","std::vector<int>",&superCluster_iphi);  
      tree->Branch("superCluster_iz","std::vector<int>",&superCluster_iz);    
      if(savePFCluster_) tree->Branch("superCluster_seedIndex","std::vector<int>",&superCluster_seedIndex);     
      if(savePFCluster_) tree->Branch("superCluster_pfClustersIndex","std::vector<std::vector<int> >",&superCluster_pfClustersIndex); 
      tree->Branch("superCluster_psCluster_energy", "std::vector<std::vector<float> >", &superCluster_psCluster_energy);
      tree->Branch("superCluster_psCluster_eta", "std::vector<std::vector<float> >", &superCluster_psCluster_eta);
      tree->Branch("superCluster_psCluster_phi", "std::vector<std::vector<float> >", &superCluster_psCluster_phi); 
      if(useRetunedSC_){   
         tree->Branch("retunedSuperCluster_energy","std::vector<float> ",&retunedSuperCluster_energy);
         tree->Branch("retunedSuperCluster_eta","std::vector<float>",&retunedSuperCluster_eta);
         tree->Branch("retunedSuperCluster_phi","std::vector<float>",&retunedSuperCluster_phi);  
         tree->Branch("retunedSuperCluster_etaWidth","std::vector<float>",&retunedSuperCluster_etaWidth);
         tree->Branch("retunedSuperCluster_phiWidth","std::vector<float>",&retunedSuperCluster_phiWidth);  
         tree->Branch("retunedSuperCluster_R","std::vector<float>",&retunedSuperCluster_R);   
         tree->Branch("retunedSuperCluster_nPFClusters","std::vector<int>",&retunedSuperCluster_nPFClusters);   
         tree->Branch("retunedSuperCluster_ieta","std::vector<int>",&retunedSuperCluster_ieta);
         tree->Branch("retunedSuperCluster_iphi","std::vector<int>",&retunedSuperCluster_iphi);  
         tree->Branch("retunedSuperCluster_iz","std::vector<int>",&retunedSuperCluster_iz);    
         if(savePFCluster_) tree->Branch("retunedSuperCluster_seedIndex","std::vector<int>",&retunedSuperCluster_seedIndex);     
         if(savePFCluster_) tree->Branch("retunedSuperCluster_pfClustersIndex","std::vector<std::vector<int> >",&retunedSuperCluster_pfClustersIndex); 
         tree->Branch("retunedSuperCluster_psCluster_energy", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_energy);
         tree->Branch("retunedSuperCluster_psCluster_eta", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_eta);
         tree->Branch("retunedSuperCluster_psCluster_phi", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_phi); 
      } 
      if(useDeepSC_){
         tree->Branch("deepSuperCluster_energy","std::vector<float> ",&deepSuperCluster_energy);
         tree->Branch("deepSuperCluster_eta","std::vector<float>",&deepSuperCluster_eta);
         tree->Branch("deepSuperCluster_phi","std::vector<float>",&deepSuperCluster_phi);  
         tree->Branch("deepSuperCluster_etaWidth","std::vector<float>",&deepSuperCluster_etaWidth);
         tree->Branch("deepSuperCluster_phiWidth","std::vector<float>",&deepSuperCluster_phiWidth);  
         tree->Branch("deepSuperCluster_R","std::vector<float>",&deepSuperCluster_R);   
         tree->Branch("deepSuperCluster_nPFClusters","std::vector<int>",&deepSuperCluster_nPFClusters);   
         tree->Branch("deepSuperCluster_ieta","std::vector<int>",&deepSuperCluster_ieta);
         tree->Branch("deepSuperCluster_iphi","std::vector<int>",&deepSuperCluster_iphi);  
         tree->Branch("deepSuperCluster_iz","std::vector<int>",&deepSuperCluster_iz);    
         if(savePFCluster_) tree->Branch("deepSuperCluster_seedIndex","std::vector<int>",&deepSuperCluster_seedIndex);     
         if(savePFCluster_) tree->Branch("deepSuperCluster_pfClustersIndex","std::vector<std::vector<int> >",&deepSuperCluster_pfClustersIndex); 
         tree->Branch("deepSuperCluster_psCluster_energy", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_energy);
         tree->Branch("deepSuperCluster_psCluster_eta", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_eta);
         tree->Branch("deepSuperCluster_psCluster_phi", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_phi); 
         tree->Branch("deepSuperClusterLWP_energy","std::vector<float> ",&deepSuperClusterLWP_energy);
         tree->Branch("deepSuperClusterLWP_eta","std::vector<float>",&deepSuperClusterLWP_eta);
         tree->Branch("deepSuperClusterLWP_phi","std::vector<float>",&deepSuperClusterLWP_phi);  
         tree->Branch("deepSuperClusterLWP_etaWidth","std::vector<float>",&deepSuperClusterLWP_etaWidth);
         tree->Branch("deepSuperClusterLWP_phiWidth","std::vector<float>",&deepSuperClusterLWP_phiWidth);  
         tree->Branch("deepSuperClusterLWP_R","std::vector<float>",&deepSuperClusterLWP_R);   
         tree->Branch("deepSuperClusterLWP_nPFClusters","std::vector<int>",&deepSuperClusterLWP_nPFClusters);   
         tree->Branch("deepSuperClusterLWP_ieta","std::vector<int>",&deepSuperClusterLWP_ieta);
         tree->Branch("deepSuperClusterLWP_iphi","std::vector<int>",&deepSuperClusterLWP_iphi);  
         tree->Branch("deepSuperClusterLWP_iz","std::vector<int>",&deepSuperClusterLWP_iz);    
         if(savePFCluster_) tree->Branch("deepSuperClusterLWP_seedIndex","std::vector<int>",&deepSuperClusterLWP_seedIndex);     
         if(savePFCluster_) tree->Branch("deepSuperClusterLWP_pfClustersIndex","std::vector<std::vector<int> >",&deepSuperClusterLWP_pfClustersIndex); 
         tree->Branch("deepSuperClusterLWP_psCluster_energy", "std::vector<std::vector<float> >", &deepSuperClusterLWP_psCluster_energy);
         tree->Branch("deepSuperClusterLWP_psCluster_eta", "std::vector<std::vector<float> >", &deepSuperClusterLWP_psCluster_eta);
         tree->Branch("deepSuperClusterLWP_psCluster_phi", "std::vector<std::vector<float> >", &deepSuperClusterLWP_psCluster_phi);
         tree->Branch("deepSuperClusterTWP_energy","std::vector<float> ",&deepSuperClusterTWP_energy);
         tree->Branch("deepSuperClusterTWP_eta","std::vector<float>",&deepSuperClusterTWP_eta);
         tree->Branch("deepSuperClusterTWP_phi","std::vector<float>",&deepSuperClusterTWP_phi);  
         tree->Branch("deepSuperClusterTWP_etaWidth","std::vector<float>",&deepSuperClusterTWP_etaWidth);
         tree->Branch("deepSuperClusterTWP_phiWidth","std::vector<float>",&deepSuperClusterTWP_phiWidth);  
         tree->Branch("deepSuperClusterTWP_R","std::vector<float>",&deepSuperClusterTWP_R);   
         tree->Branch("deepSuperClusterTWP_nPFClusters","std::vector<int>",&deepSuperClusterTWP_nPFClusters);   
         tree->Branch("deepSuperClusterTWP_ieta","std::vector<int>",&deepSuperClusterTWP_ieta);
         tree->Branch("deepSuperClusterTWP_iphi","std::vector<int>",&deepSuperClusterTWP_iphi);  
         tree->Branch("deepSuperClusterTWP_iz","std::vector<int>",&deepSuperClusterTWP_iz);    
         if(savePFCluster_) tree->Branch("deepSuperClusterTWP_seedIndex","std::vector<int>",&deepSuperClusterTWP_seedIndex);     
         if(savePFCluster_) tree->Branch("deepSuperClusterTWP_pfClustersIndex","std::vector<std::vector<int> >",&deepSuperClusterTWP_pfClustersIndex); 
         tree->Branch("deepSuperClusterTWP_psCluster_energy", "std::vector<std::vector<float> >", &deepSuperClusterTWP_psCluster_energy);
         tree->Branch("deepSuperClusterTWP_psCluster_eta", "std::vector<std::vector<float> >", &deepSuperClusterTWP_psCluster_eta);
         tree->Branch("deepSuperClusterTWP_psCluster_phi", "std::vector<std::vector<float> >", &deepSuperClusterTWP_psCluster_phi);
      }   
      if(saveCaloParticles_){
         tree->Branch("superCluster_dR_genScore_MatchedIndex","std::vector<int>",&superCluster_dR_genScore_MatchedIndex);
         tree->Branch("superCluster_dR_simScore_MatchedIndex","std::vector<int>",&superCluster_dR_simScore_MatchedIndex);
         tree->Branch("superCluster_sim_fraction_old_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_old_MatchedIndex); 
         if(!saveScores_) tree->Branch("superCluster_simScore_MatchedIndex","std::vector<int>",&superCluster_simScore_MatchedIndex);
         if(useRetunedSC_){  
            tree->Branch("retunedSuperCluster_dR_genScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_dR_genScore_MatchedIndex);
            tree->Branch("retunedSuperCluster_dR_simScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("retunedSuperCluster_sim_fraction_old_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_old_MatchedIndex); 
            if(!saveScores_) tree->Branch("retunedSuperCluster_simScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_simScore_MatchedIndex);
         } 
         if(saveScores_){
            tree->Branch("superCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&superCluster_n_shared_xtals_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_1MeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_5MeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_10MeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_50MeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_100MeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_500MeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_1GeVCut_MatchedIndex);
            tree->Branch("superCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&superCluster_sim_rechit_diff_MatchedIndex);
            tree->Branch("superCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_rechit_fraction_MatchedIndex);   
            tree->Branch("superCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&superCluster_global_sim_rechit_fraction_MatchedIndex); 
            tree->Branch("superCluster_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&superCluster_hgcal_caloToCluster_MatchedIndex); 
            tree->Branch("superCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&superCluster_hgcal_clusterToCalo_MatchedIndex); 
            tree->Branch("superCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_rechit_combined_fraction_MatchedIndex); 
            tree->Branch("superCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&superCluster_rechit_sim_combined_fraction_MatchedIndex);
            if(useRetunedSC_){  
               tree->Branch("retunedSuperCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_n_shared_xtals_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_rechit_diff_MatchedIndex);
               tree->Branch("retunedSuperCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("retunedSuperCluster_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&retunedSuperCluster_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex);   
            } 
         } 
         if(useDeepSC_){
            tree->Branch("deepSuperCluster_dR_genScore_MatchedIndex","std::vector<int>",&deepSuperCluster_dR_genScore_MatchedIndex);
            tree->Branch("deepSuperCluster_dR_simScore_MatchedIndex","std::vector<int>",&deepSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("deepSuperCluster_sim_fraction_old_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_old_MatchedIndex); 
            if(!saveScores_) tree->Branch("deepSuperCluster_simScore_MatchedIndex","std::vector<int>",&deepSuperCluster_simScore_MatchedIndex);
            tree->Branch("deepSuperClusterLWP_dR_genScore_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_dR_genScore_MatchedIndex);
            tree->Branch("deepSuperClusterLWP_dR_simScore_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_dR_simScore_MatchedIndex);
            tree->Branch("deepSuperClusterLWP_sim_fraction_old_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_old_MatchedIndex); 
            if(!saveScores_) tree->Branch("deepSuperClusterLWP_simScore_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_simScore_MatchedIndex);
            tree->Branch("deepSuperClusterTWP_dR_genScore_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_dR_genScore_MatchedIndex);
            tree->Branch("deepSuperClusterTWP_dR_simScore_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_dR_simScore_MatchedIndex);
            tree->Branch("deepSuperClusterTWP_sim_fraction_old_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_old_MatchedIndex); 
            if(!saveScores_) tree->Branch("deepSuperClusterTWP_simScore_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_simScore_MatchedIndex);
            if(saveScores_){
               tree->Branch("deepSuperCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&deepSuperCluster_n_shared_xtals_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_rechit_diff_MatchedIndex);
               tree->Branch("deepSuperCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("deepSuperCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("deepSuperCluster_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&deepSuperCluster_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("deepSuperCluster_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&deepSuperCluster_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex); 
               tree->Branch("deepSuperClusterLWP_n_shared_xtals_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_n_shared_xtals_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_rechit_diff_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_rechit_diff_MatchedIndex);
               tree->Branch("deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex); 
               tree->Branch("deepSuperClusterTWP_n_shared_xtals_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_n_shared_xtals_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_rechit_diff_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_rechit_diff_MatchedIndex);
               tree->Branch("deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex);   
               tree->Branch("deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex); 
               tree->Branch("deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex); 
               tree->Branch("deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex); 
               tree->Branch("deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex); 
               tree->Branch("deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex","std::vector<int>",&deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex);  
            } 
         } 
      } 
      if(saveCaloParticles_){
         tree->Branch("superCluster_dR_genScore","std::vector<std::vector<double> >",&superCluster_dR_genScore);
         tree->Branch("superCluster_dR_simScore","std::vector<std::vector<double> >",&superCluster_dR_simScore);
         tree->Branch("superCluster_sim_fraction_old","std::vector<std::vector<double> >",&superCluster_sim_fraction_old);
         if(!saveScores_) tree->Branch("superCluster_simScore","std::vector<std::vector<double> >",&superCluster_simScore);
         if(useRetunedSC_){ 
            tree->Branch("retunedSuperCluster_dR_genScore","std::vector<std::vector<double> >",&retunedSuperCluster_dR_genScore);
            tree->Branch("retunedSuperCluster_dR_simScore","std::vector<std::vector<double> >",&retunedSuperCluster_dR_simScore);
            tree->Branch("retunedSuperCluster_sim_fraction_old","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_old);
            if(!saveScores_) tree->Branch("retunedSuperCluster_simScore","std::vector<std::vector<double> >",&retunedSuperCluster_simScore);
         } 
         if(saveScores_){
            tree->Branch("superCluster_n_shared_xtals","std::vector<std::vector<double> >",&superCluster_n_shared_xtals);
            tree->Branch("superCluster_sim_fraction","std::vector<std::vector<double> >",&superCluster_sim_fraction);
            tree->Branch("superCluster_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_1MeVCut); 
            tree->Branch("superCluster_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_5MeVCut);
            tree->Branch("superCluster_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_10MeVCut);
            tree->Branch("superCluster_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_50MeVCut);
            tree->Branch("superCluster_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_100MeVCut);  
            tree->Branch("superCluster_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_500MeVCut);  
            tree->Branch("superCluster_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&superCluster_sim_fraction_1GeVCut);     
            tree->Branch("superCluster_sim_rechit_diff","std::vector<std::vector<double> >",&superCluster_sim_rechit_diff);
            tree->Branch("superCluster_sim_rechit_fraction","std::vector<std::vector<double> >",&superCluster_sim_rechit_fraction);   
            tree->Branch("superCluster_global_sim_rechit_fraction","std::vector<std::vector<double> >",&superCluster_global_sim_rechit_fraction); 
            tree->Branch("superCluster_hgcal_caloToCluster","std::vector<std::vector<double> >",&superCluster_hgcal_caloToCluster); 
            tree->Branch("superCluster_hgcal_clusterToCalo","std::vector<std::vector<double> >",&superCluster_hgcal_clusterToCalo); 
            tree->Branch("superCluster_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&superCluster_sim_rechit_combined_fraction); 
            tree->Branch("superCluster_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&superCluster_rechit_sim_combined_fraction);  
            if(useRetunedSC_){
               tree->Branch("retunedSuperCluster_n_shared_xtals","std::vector<std::vector<double> >",&retunedSuperCluster_n_shared_xtals);
               tree->Branch("retunedSuperCluster_sim_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction);
               tree->Branch("retunedSuperCluster_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_1MeVCut); 
               tree->Branch("retunedSuperCluster_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_5MeVCut);
               tree->Branch("retunedSuperCluster_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_10MeVCut);
               tree->Branch("retunedSuperCluster_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_50MeVCut);
               tree->Branch("retunedSuperCluster_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_100MeVCut);  
               tree->Branch("retunedSuperCluster_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_500MeVCut);  
               tree->Branch("retunedSuperCluster_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_1GeVCut);     
               tree->Branch("retunedSuperCluster_sim_rechit_diff","std::vector<std::vector<double> >",&retunedSuperCluster_sim_rechit_diff);
               tree->Branch("retunedSuperCluster_sim_rechit_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_rechit_fraction);   
               tree->Branch("retunedSuperCluster_global_sim_rechit_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_global_sim_rechit_fraction); 
               tree->Branch("retunedSuperCluster_hgcal_caloToCluster","std::vector<std::vector<double> >",&retunedSuperCluster_hgcal_caloToCluster); 
               tree->Branch("retunedSuperCluster_hgcal_clusterToCalo","std::vector<std::vector<double> >",&retunedSuperCluster_hgcal_clusterToCalo); 
               tree->Branch("retunedSuperCluster_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_rechit_combined_fraction); 
               tree->Branch("retunedSuperCluster_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_rechit_sim_combined_fraction);    
            }  
         }
         if(useDeepSC_){
            tree->Branch("deepSuperCluster_dR_genScore","std::vector<std::vector<double> >",&deepSuperCluster_dR_genScore);
            tree->Branch("deepSuperCluster_dR_simScore","std::vector<std::vector<double> >",&deepSuperCluster_dR_simScore);
            tree->Branch("deepSuperCluster_sim_fraction_old","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_old);
            if(!saveScores_) tree->Branch("deepSuperCluster_simScore","std::vector<std::vector<double> >",&deepSuperCluster_simScore);
            tree->Branch("deepSuperClusterLWP_dR_genScore","std::vector<std::vector<double> >",&deepSuperClusterLWP_dR_genScore);
            tree->Branch("deepSuperClusterLWP_dR_simScore","std::vector<std::vector<double> >",&deepSuperClusterLWP_dR_simScore);
            tree->Branch("deepSuperClusterLWP_sim_fraction_old","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_old);
            if(!saveScores_) tree->Branch("deepSuperClusterLWP_simScore","std::vector<std::vector<double> >",&deepSuperClusterLWP_simScore);
            tree->Branch("deepSuperClusterTWP_dR_genScore","std::vector<std::vector<double> >",&deepSuperClusterTWP_dR_genScore);
            tree->Branch("deepSuperClusterTWP_dR_simScore","std::vector<std::vector<double> >",&deepSuperClusterTWP_dR_simScore);
            tree->Branch("deepSuperClusterTWP_sim_fraction_old","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_old);
            if(!saveScores_) tree->Branch("deepSuperClusterTWP_simScore","std::vector<std::vector<double> >",&deepSuperClusterTWP_simScore);
            if(saveScores_){
               tree->Branch("deepSuperCluster_n_shared_xtals","std::vector<std::vector<double> >",&deepSuperCluster_n_shared_xtals);
               tree->Branch("deepSuperCluster_sim_fraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction);
               tree->Branch("deepSuperCluster_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_1MeVCut); 
               tree->Branch("deepSuperCluster_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_5MeVCut);
               tree->Branch("deepSuperCluster_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_10MeVCut);
               tree->Branch("deepSuperCluster_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_50MeVCut);
               tree->Branch("deepSuperCluster_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_100MeVCut);  
               tree->Branch("deepSuperCluster_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_500MeVCut);  
               tree->Branch("deepSuperCluster_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_1GeVCut);      
               tree->Branch("deepSuperCluster_sim_rechit_diff","std::vector<std::vector<double> >",&deepSuperCluster_sim_rechit_diff);
               tree->Branch("deepSuperCluster_sim_rechit_fraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_rechit_fraction);   
               tree->Branch("deepSuperCluster_global_sim_rechit_fraction","std::vector<std::vector<double> >",&deepSuperCluster_global_sim_rechit_fraction); 
               tree->Branch("deepSuperCluster_hgcal_caloToCluster","std::vector<std::vector<double> >",&deepSuperCluster_hgcal_caloToCluster); 
               tree->Branch("deepSuperCluster_hgcal_clusterToCalo","std::vector<std::vector<double> >",&deepSuperCluster_hgcal_clusterToCalo); 
               tree->Branch("deepSuperCluster_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_rechit_combined_fraction); 
               tree->Branch("deepSuperCluster_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&deepSuperCluster_rechit_sim_combined_fraction); 
               tree->Branch("deepSuperClusterLWP_n_shared_xtals","std::vector<std::vector<double> >",&deepSuperClusterLWP_n_shared_xtals);
               tree->Branch("deepSuperClusterLWP_sim_fraction","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction);
               tree->Branch("deepSuperClusterLWP_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_1MeVCut); 
               tree->Branch("deepSuperClusterLWP_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_5MeVCut);
               tree->Branch("deepSuperClusterLWP_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_10MeVCut);
               tree->Branch("deepSuperClusterLWP_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_50MeVCut);
               tree->Branch("deepSuperClusterLWP_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_100MeVCut);  
               tree->Branch("deepSuperClusterLWP_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_500MeVCut);  
               tree->Branch("deepSuperClusterLWP_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_fraction_1GeVCut);      
               tree->Branch("deepSuperClusterLWP_sim_rechit_diff","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_rechit_diff);
               tree->Branch("deepSuperClusterLWP_sim_rechit_fraction","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_rechit_fraction);   
               tree->Branch("deepSuperClusterLWP_global_sim_rechit_fraction","std::vector<std::vector<double> >",&deepSuperClusterLWP_global_sim_rechit_fraction); 
               tree->Branch("deepSuperClusterLWP_hgcal_caloToCluster","std::vector<std::vector<double> >",&deepSuperClusterLWP_hgcal_caloToCluster); 
               tree->Branch("deepSuperClusterLWP_hgcal_clusterToCalo","std::vector<std::vector<double> >",&deepSuperClusterLWP_hgcal_clusterToCalo); 
               tree->Branch("deepSuperClusterLWP_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&deepSuperClusterLWP_sim_rechit_combined_fraction); 
               tree->Branch("deepSuperClusterLWP_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&deepSuperClusterLWP_rechit_sim_combined_fraction); 
               tree->Branch("deepSuperClusterTWP_n_shared_xtals","std::vector<std::vector<double> >",&deepSuperClusterTWP_n_shared_xtals);
               tree->Branch("deepSuperClusterTWP_sim_fraction","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction);
               tree->Branch("deepSuperClusterTWP_sim_fraction_1MeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_1MeVCut); 
               tree->Branch("deepSuperClusterTWP_sim_fraction_5MeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_5MeVCut);
               tree->Branch("deepSuperClusterTWP_sim_fraction_10MeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_10MeVCut);
               tree->Branch("deepSuperClusterTWP_sim_fraction_50MeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_50MeVCut);
               tree->Branch("deepSuperClusterTWP_sim_fraction_100MeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_100MeVCut);  
               tree->Branch("deepSuperClusterTWP_sim_fraction_500MeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_500MeVCut);  
               tree->Branch("deepSuperClusterTWP_sim_fraction_1GeVCut","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_fraction_1GeVCut);      
               tree->Branch("deepSuperClusterTWP_sim_rechit_diff","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_rechit_diff);
               tree->Branch("deepSuperClusterTWP_sim_rechit_fraction","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_rechit_fraction);   
               tree->Branch("deepSuperClusterTWP_global_sim_rechit_fraction","std::vector<std::vector<double> >",&deepSuperClusterTWP_global_sim_rechit_fraction); 
               tree->Branch("deepSuperClusterTWP_hgcal_caloToCluster","std::vector<std::vector<double> >",&deepSuperClusterTWP_hgcal_caloToCluster); 
               tree->Branch("deepSuperClusterTWP_hgcal_clusterToCalo","std::vector<std::vector<double> >",&deepSuperClusterTWP_hgcal_clusterToCalo); 
               tree->Branch("deepSuperClusterTWP_sim_rechit_combined_fraction","std::vector<std::vector<double> >",&deepSuperClusterTWP_sim_rechit_combined_fraction); 
               tree->Branch("deepSuperClusterTWP_rechit_sim_combined_fraction","std::vector<std::vector<double> >",&deepSuperClusterTWP_rechit_sim_combined_fraction);
            }
         }  
      }  
   }
   if(savePFCluster_ && saveShowerShapes_){  
      tree->Branch("pfCluster_e5x5","std::vector<float>",&pfCluster_e5x5);
      tree->Branch("pfCluster_e2x2Ratio","std::vector<float>",&pfCluster_e2x2Ratio);
      tree->Branch("pfCluster_e3x3Ratio","std::vector<float>",&pfCluster_e3x3Ratio);
      tree->Branch("pfCluster_eMaxRatio","std::vector<float>",&pfCluster_eMaxRatio);
      tree->Branch("pfCluster_e2ndRatio","std::vector<float>",&pfCluster_e2ndRatio);
      tree->Branch("pfCluster_eTopRatio","std::vector<float>",&pfCluster_eTopRatio);
      tree->Branch("pfCluster_eRightRatio","std::vector<float>",&pfCluster_eRightRatio);
      tree->Branch("pfCluster_eBottomRatio","std::vector<float>",&pfCluster_eBottomRatio);
      tree->Branch("pfCluster_eLeftRatio","std::vector<float>",&pfCluster_eLeftRatio);
      tree->Branch("pfCluster_e2x5MaxRatio","std::vector<float>",&pfCluster_e2x5MaxRatio);
      tree->Branch("pfCluster_e2x5TopRatio","std::vector<float>",&pfCluster_e2x5TopRatio);
      tree->Branch("pfCluster_e2x5RightRatio","std::vector<float>",&pfCluster_e2x5RightRatio);
      tree->Branch("pfCluster_e2x5BottomRatio","std::vector<float>",&pfCluster_e2x5BottomRatio); 
      tree->Branch("pfCluster_e2x5LeftRatio","std::vector<float>",&pfCluster_e2x5LeftRatio); 
      tree->Branch("pfCluster_swissCross","std::vector<float>",&pfCluster_swissCross); 
      tree->Branch("pfCluster_r9","std::vector<float>",&pfCluster_r9);
      tree->Branch("pfCluster_sigmaIetaIeta","std::vector<float>",&pfCluster_sigmaIetaIeta);
      tree->Branch("pfCluster_sigmaIetaIphi","std::vector<float>",&pfCluster_sigmaIetaIphi);
      tree->Branch("pfCluster_sigmaIphiIphi","std::vector<float>",&pfCluster_sigmaIphiIphi);
      tree->Branch("pfCluster_full5x5_e5x5","std::vector<float>",&pfCluster_full5x5_e5x5);
      tree->Branch("pfCluster_full5x5_e2x2Ratio","std::vector<float>",&pfCluster_full5x5_e2x2Ratio);
      tree->Branch("pfCluster_full5x5_e3x3Ratio","std::vector<float>",&pfCluster_full5x5_e3x3Ratio);
      tree->Branch("pfCluster_full5x5_eMaxRatio","std::vector<float>",&pfCluster_full5x5_eMaxRatio);
      tree->Branch("pfCluster_full5x5_e2ndRatio","std::vector<float>",&pfCluster_full5x5_e2ndRatio);
      tree->Branch("pfCluster_full5x5_eTopRatio","std::vector<float>",&pfCluster_full5x5_eTopRatio);
      tree->Branch("pfCluster_full5x5_eRightRatio","std::vector<float>",&pfCluster_full5x5_eRightRatio);
      tree->Branch("pfCluster_full5x5_eBottomRatio","std::vector<float>",&pfCluster_full5x5_eBottomRatio);
      tree->Branch("pfCluster_full5x5_eLeftRatio","std::vector<float>",&pfCluster_full5x5_eLeftRatio);
      tree->Branch("pfCluster_full5x5_e2x5MaxRatio","std::vector<float>",&pfCluster_full5x5_e2x5MaxRatio);
      tree->Branch("pfCluster_full5x5_e2x5TopRatio","std::vector<float>",&pfCluster_full5x5_e2x5TopRatio);
      tree->Branch("pfCluster_full5x5_e2x5RightRatio","std::vector<float>",&pfCluster_full5x5_e2x5RightRatio);
      tree->Branch("pfCluster_full5x5_e2x5BottomRatio","std::vector<float>",&pfCluster_full5x5_e2x5BottomRatio); 
      tree->Branch("pfCluster_full5x5_e2x5LeftRatio","std::vector<float>",&pfCluster_full5x5_e2x5LeftRatio); 
      tree->Branch("pfCluster_full5x5_swissCross","std::vector<float>",&pfCluster_full5x5_swissCross); 
      tree->Branch("pfCluster_full5x5_r9","std::vector<float>",&pfCluster_full5x5_r9);
      tree->Branch("pfCluster_full5x5_sigmaIetaIeta","std::vector<float>",&pfCluster_full5x5_sigmaIetaIeta);
      tree->Branch("pfCluster_full5x5_sigmaIetaIphi","std::vector<float>",&pfCluster_full5x5_sigmaIetaIphi);
      tree->Branch("pfCluster_full5x5_sigmaIphiIphi","std::vector<float>",&pfCluster_full5x5_sigmaIphiIphi);
   }
   if(saveSuperCluster_ && saveShowerShapes_){  
      tree->Branch("superCluster_e5x5","std::vector<float>",&superCluster_e5x5);
      tree->Branch("superCluster_e2x2Ratio","std::vector<float>",&superCluster_e2x2Ratio);
      tree->Branch("superCluster_e3x3Ratio","std::vector<float>",&superCluster_e3x3Ratio);
      tree->Branch("superCluster_eMaxRatio","std::vector<float>",&superCluster_eMaxRatio);
      tree->Branch("superCluster_e2ndRatio","std::vector<float>",&superCluster_e2ndRatio);
      tree->Branch("superCluster_eTopRatio","std::vector<float>",&superCluster_eTopRatio);
      tree->Branch("superCluster_eRightRatio","std::vector<float>",&superCluster_eRightRatio);
      tree->Branch("superCluster_eBottomRatio","std::vector<float>",&superCluster_eBottomRatio);
      tree->Branch("superCluster_eLeftRatio","std::vector<float>",&superCluster_eLeftRatio);
      tree->Branch("superCluster_e2x5MaxRatio","std::vector<float>",&superCluster_e2x5MaxRatio);
      tree->Branch("superCluster_e2x5TopRatio","std::vector<float>",&superCluster_e2x5TopRatio);
      tree->Branch("superCluster_e2x5RightRatio","std::vector<float>",&superCluster_e2x5RightRatio);
      tree->Branch("superCluster_e2x5BottomRatio","std::vector<float>",&superCluster_e2x5BottomRatio); 
      tree->Branch("superCluster_e2x5LeftRatio","std::vector<float>",&superCluster_e2x5LeftRatio); 
      tree->Branch("superCluster_swissCross","std::vector<float>",&superCluster_swissCross); 
      tree->Branch("superCluster_r9","std::vector<float>",&superCluster_r9);
      tree->Branch("superCluster_sigmaIetaIeta","std::vector<float>",&superCluster_sigmaIetaIeta);
      tree->Branch("superCluster_sigmaIetaIphi","std::vector<float>",&superCluster_sigmaIetaIphi);
      tree->Branch("superCluster_sigmaIphiIphi","std::vector<float>",&superCluster_sigmaIphiIphi);
      tree->Branch("superCluster_full5x5_e5x5","std::vector<float>",&superCluster_full5x5_e5x5);
      tree->Branch("superCluster_full5x5_e2x2Ratio","std::vector<float>",&superCluster_full5x5_e2x2Ratio);
      tree->Branch("superCluster_full5x5_e3x3Ratio","std::vector<float>",&superCluster_full5x5_e3x3Ratio);
      tree->Branch("superCluster_full5x5_eMaxRatio","std::vector<float>",&superCluster_full5x5_eMaxRatio);
      tree->Branch("superCluster_full5x5_e2ndRatio","std::vector<float>",&superCluster_full5x5_e2ndRatio);
      tree->Branch("superCluster_full5x5_eTopRatio","std::vector<float>",&superCluster_full5x5_eTopRatio);
      tree->Branch("superCluster_full5x5_eRightRatio","std::vector<float>",&superCluster_full5x5_eRightRatio);
      tree->Branch("superCluster_full5x5_eBottomRatio","std::vector<float>",&superCluster_full5x5_eBottomRatio);
      tree->Branch("superCluster_full5x5_eLeftRatio","std::vector<float>",&superCluster_full5x5_eLeftRatio);
      tree->Branch("superCluster_full5x5_e2x5MaxRatio","std::vector<float>",&superCluster_full5x5_e2x5MaxRatio);
      tree->Branch("superCluster_full5x5_e2x5TopRatio","std::vector<float>",&superCluster_full5x5_e2x5TopRatio);
      tree->Branch("superCluster_full5x5_e2x5RightRatio","std::vector<float>",&superCluster_full5x5_e2x5RightRatio);
      tree->Branch("superCluster_full5x5_e2x5BottomRatio","std::vector<float>",&superCluster_full5x5_e2x5BottomRatio); 
      tree->Branch("superCluster_full5x5_e2x5LeftRatio","std::vector<float>",&superCluster_full5x5_e2x5LeftRatio); 
      tree->Branch("superCluster_full5x5_swissCross","std::vector<float>",&superCluster_full5x5_swissCross); 
      tree->Branch("superCluster_full5x5_r9","std::vector<float>",&superCluster_full5x5_r9);
      tree->Branch("superCluster_full5x5_sigmaIetaIeta","std::vector<float>",&superCluster_full5x5_sigmaIetaIeta);
      tree->Branch("superCluster_full5x5_sigmaIetaIphi","std::vector<float>",&superCluster_full5x5_sigmaIetaIphi);
      tree->Branch("superCluster_full5x5_sigmaIphiIphi","std::vector<float>",&superCluster_full5x5_sigmaIphiIphi);
      if(useHcalTowers_ ){
         tree->Branch("superCluster_HoEraw","std::vector<float>",&superCluster_HoEraw);
         tree->Branch("superCluster_HoErawBC","std::vector<float>",&superCluster_HoErawBC);
      }   
      if(useRetunedSC_){
         tree->Branch("retunedSuperCluster_e5x5","std::vector<float>",&retunedSuperCluster_e5x5);
         tree->Branch("retunedSuperCluster_e2x2Ratio","std::vector<float>",&retunedSuperCluster_e2x2Ratio);
         tree->Branch("retunedSuperCluster_e3x3Ratio","std::vector<float>",&retunedSuperCluster_e3x3Ratio);
         tree->Branch("retunedSuperCluster_eMaxRatio","std::vector<float>",&retunedSuperCluster_eMaxRatio);
         tree->Branch("retunedSuperCluster_e2ndRatio","std::vector<float>",&retunedSuperCluster_e2ndRatio);
         tree->Branch("retunedSuperCluster_eTopRatio","std::vector<float>",&retunedSuperCluster_eTopRatio);
         tree->Branch("retunedSuperCluster_eRightRatio","std::vector<float>",&retunedSuperCluster_eRightRatio);
         tree->Branch("retunedSuperCluster_eBottomRatio","std::vector<float>",&retunedSuperCluster_eBottomRatio);
         tree->Branch("retunedSuperCluster_eLeftRatio","std::vector<float>",&retunedSuperCluster_eLeftRatio);
         tree->Branch("retunedSuperCluster_e2x5MaxRatio","std::vector<float>",&retunedSuperCluster_e2x5MaxRatio);
         tree->Branch("retunedSuperCluster_e2x5TopRatio","std::vector<float>",&retunedSuperCluster_e2x5TopRatio);
         tree->Branch("retunedSuperCluster_e2x5RightRatio","std::vector<float>",&retunedSuperCluster_e2x5RightRatio);
         tree->Branch("retunedSuperCluster_e2x5BottomRatio","std::vector<float>",&retunedSuperCluster_e2x5BottomRatio); 
         tree->Branch("retunedSuperCluster_e2x5LeftRatio","std::vector<float>",&retunedSuperCluster_e2x5LeftRatio); 
         tree->Branch("retunedSuperCluster_swissCross","std::vector<float>",&retunedSuperCluster_swissCross); 
         tree->Branch("retunedSuperCluster_r9","std::vector<float>",&retunedSuperCluster_r9);
         tree->Branch("retunedSuperCluster_sigmaIetaIeta","std::vector<float>",&retunedSuperCluster_sigmaIetaIeta);
         tree->Branch("retunedSuperCluster_sigmaIetaIphi","std::vector<float>",&retunedSuperCluster_sigmaIetaIphi);
         tree->Branch("retunedSuperCluster_sigmaIphiIphi","std::vector<float>",&retunedSuperCluster_sigmaIphiIphi);
         tree->Branch("retunedSuperCluster_full5x5_e5x5","std::vector<float>",&retunedSuperCluster_full5x5_e5x5);
         tree->Branch("retunedSuperCluster_full5x5_e2x2Ratio","std::vector<float>",&retunedSuperCluster_full5x5_e2x2Ratio);
         tree->Branch("retunedSuperCluster_full5x5_e3x3Ratio","std::vector<float>",&retunedSuperCluster_full5x5_e3x3Ratio);
         tree->Branch("retunedSuperCluster_full5x5_eMaxRatio","std::vector<float>",&retunedSuperCluster_full5x5_eMaxRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2ndRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2ndRatio);
         tree->Branch("retunedSuperCluster_full5x5_eTopRatio","std::vector<float>",&retunedSuperCluster_full5x5_eTopRatio);
         tree->Branch("retunedSuperCluster_full5x5_eRightRatio","std::vector<float>",&retunedSuperCluster_full5x5_eRightRatio);
         tree->Branch("retunedSuperCluster_full5x5_eBottomRatio","std::vector<float>",&retunedSuperCluster_full5x5_eBottomRatio);
         tree->Branch("retunedSuperCluster_full5x5_eLeftRatio","std::vector<float>",&retunedSuperCluster_full5x5_eLeftRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5MaxRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5MaxRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5TopRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5TopRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5RightRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5RightRatio);
         tree->Branch("retunedSuperCluster_full5x5_e2x5BottomRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5BottomRatio); 
         tree->Branch("retunedSuperCluster_full5x5_e2x5LeftRatio","std::vector<float>",&retunedSuperCluster_full5x5_e2x5LeftRatio); 
         tree->Branch("retunedSuperCluster_full5x5_swissCross","std::vector<float>",&retunedSuperCluster_full5x5_swissCross); 
         tree->Branch("retunedSuperCluster_full5x5_r9","std::vector<float>",&retunedSuperCluster_full5x5_r9);
         tree->Branch("retunedSuperCluster_full5x5_sigmaIetaIeta","std::vector<float>",&retunedSuperCluster_full5x5_sigmaIetaIeta);
         tree->Branch("retunedSuperCluster_full5x5_sigmaIetaIphi","std::vector<float>",&retunedSuperCluster_full5x5_sigmaIetaIphi);
         tree->Branch("retunedSuperCluster_full5x5_sigmaIphiIphi","std::vector<float>",&retunedSuperCluster_full5x5_sigmaIphiIphi);
         if(useHcalTowers_ ){
            tree->Branch("retunedSuperCluster_HoEraw","std::vector<float>",&retunedSuperCluster_HoEraw);
            tree->Branch("retunedSuperCluster_HoErawBC","std::vector<float>",&retunedSuperCluster_HoErawBC);
         }   
      }
      if(useDeepSC_){
         tree->Branch("deepSuperCluster_e5x5","std::vector<float>",&deepSuperCluster_e5x5);
         tree->Branch("deepSuperCluster_e2x2Ratio","std::vector<float>",&deepSuperCluster_e2x2Ratio);
         tree->Branch("deepSuperCluster_e3x3Ratio","std::vector<float>",&deepSuperCluster_e3x3Ratio);
         tree->Branch("deepSuperCluster_eMaxRatio","std::vector<float>",&deepSuperCluster_eMaxRatio);
         tree->Branch("deepSuperCluster_e2ndRatio","std::vector<float>",&deepSuperCluster_e2ndRatio);
         tree->Branch("deepSuperCluster_eTopRatio","std::vector<float>",&deepSuperCluster_eTopRatio);
         tree->Branch("deepSuperCluster_eRightRatio","std::vector<float>",&deepSuperCluster_eRightRatio);
         tree->Branch("deepSuperCluster_eBottomRatio","std::vector<float>",&deepSuperCluster_eBottomRatio);
         tree->Branch("deepSuperCluster_eLeftRatio","std::vector<float>",&deepSuperCluster_eLeftRatio);
         tree->Branch("deepSuperCluster_e2x5MaxRatio","std::vector<float>",&deepSuperCluster_e2x5MaxRatio);
         tree->Branch("deepSuperCluster_e2x5TopRatio","std::vector<float>",&deepSuperCluster_e2x5TopRatio);
         tree->Branch("deepSuperCluster_e2x5RightRatio","std::vector<float>",&deepSuperCluster_e2x5RightRatio);
         tree->Branch("deepSuperCluster_e2x5BottomRatio","std::vector<float>",&deepSuperCluster_e2x5BottomRatio); 
         tree->Branch("deepSuperCluster_e2x5LeftRatio","std::vector<float>",&deepSuperCluster_e2x5LeftRatio); 
         tree->Branch("deepSuperCluster_swissCross","std::vector<float>",&deepSuperCluster_swissCross); 
         tree->Branch("deepSuperCluster_r9","std::vector<float>",&deepSuperCluster_r9);
         tree->Branch("deepSuperCluster_sigmaIetaIeta","std::vector<float>",&deepSuperCluster_sigmaIetaIeta);
         tree->Branch("deepSuperCluster_sigmaIetaIphi","std::vector<float>",&deepSuperCluster_sigmaIetaIphi);
         tree->Branch("deepSuperCluster_sigmaIphiIphi","std::vector<float>",&deepSuperCluster_sigmaIphiIphi);
         tree->Branch("deepSuperCluster_full5x5_e5x5","std::vector<float>",&deepSuperCluster_full5x5_e5x5);
         tree->Branch("deepSuperCluster_full5x5_e2x2Ratio","std::vector<float>",&deepSuperCluster_full5x5_e2x2Ratio);
         tree->Branch("deepSuperCluster_full5x5_e3x3Ratio","std::vector<float>",&deepSuperCluster_full5x5_e3x3Ratio);
         tree->Branch("deepSuperCluster_full5x5_eMaxRatio","std::vector<float>",&deepSuperCluster_full5x5_eMaxRatio);
         tree->Branch("deepSuperCluster_full5x5_e2ndRatio","std::vector<float>",&deepSuperCluster_full5x5_e2ndRatio);
         tree->Branch("deepSuperCluster_full5x5_eTopRatio","std::vector<float>",&deepSuperCluster_full5x5_eTopRatio);
         tree->Branch("deepSuperCluster_full5x5_eRightRatio","std::vector<float>",&deepSuperCluster_full5x5_eRightRatio);
         tree->Branch("deepSuperCluster_full5x5_eBottomRatio","std::vector<float>",&deepSuperCluster_full5x5_eBottomRatio);
         tree->Branch("deepSuperCluster_full5x5_eLeftRatio","std::vector<float>",&deepSuperCluster_full5x5_eLeftRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5MaxRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5MaxRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5TopRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5TopRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5RightRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5RightRatio);
         tree->Branch("deepSuperCluster_full5x5_e2x5BottomRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5BottomRatio); 
         tree->Branch("deepSuperCluster_full5x5_e2x5LeftRatio","std::vector<float>",&deepSuperCluster_full5x5_e2x5LeftRatio); 
         tree->Branch("deepSuperCluster_full5x5_swissCross","std::vector<float>",&deepSuperCluster_full5x5_swissCross); 
         tree->Branch("deepSuperCluster_full5x5_r9","std::vector<float>",&deepSuperCluster_full5x5_r9);
         tree->Branch("deepSuperCluster_full5x5_sigmaIetaIeta","std::vector<float>",&deepSuperCluster_full5x5_sigmaIetaIeta);
         tree->Branch("deepSuperCluster_full5x5_sigmaIetaIphi","std::vector<float>",&deepSuperCluster_full5x5_sigmaIetaIphi);
         tree->Branch("deepSuperCluster_full5x5_sigmaIphiIphi","std::vector<float>",&deepSuperCluster_full5x5_sigmaIphiIphi); 
         tree->Branch("deepSuperClusterLWP_e5x5","std::vector<float>",&deepSuperClusterLWP_e5x5);
         tree->Branch("deepSuperClusterLWP_e2x2Ratio","std::vector<float>",&deepSuperClusterLWP_e2x2Ratio);
         tree->Branch("deepSuperClusterLWP_e3x3Ratio","std::vector<float>",&deepSuperClusterLWP_e3x3Ratio);
         tree->Branch("deepSuperClusterLWP_eMaxRatio","std::vector<float>",&deepSuperClusterLWP_eMaxRatio);
         tree->Branch("deepSuperClusterLWP_e2ndRatio","std::vector<float>",&deepSuperClusterLWP_e2ndRatio);
         tree->Branch("deepSuperClusterLWP_eTopRatio","std::vector<float>",&deepSuperClusterLWP_eTopRatio);
         tree->Branch("deepSuperClusterLWP_eRightRatio","std::vector<float>",&deepSuperClusterLWP_eRightRatio);
         tree->Branch("deepSuperClusterLWP_eBottomRatio","std::vector<float>",&deepSuperClusterLWP_eBottomRatio);
         tree->Branch("deepSuperClusterLWP_eLeftRatio","std::vector<float>",&deepSuperClusterLWP_eLeftRatio);
         tree->Branch("deepSuperClusterLWP_e2x5MaxRatio","std::vector<float>",&deepSuperClusterLWP_e2x5MaxRatio);
         tree->Branch("deepSuperClusterLWP_e2x5TopRatio","std::vector<float>",&deepSuperClusterLWP_e2x5TopRatio);
         tree->Branch("deepSuperClusterLWP_e2x5RightRatio","std::vector<float>",&deepSuperClusterLWP_e2x5RightRatio);
         tree->Branch("deepSuperClusterLWP_e2x5BottomRatio","std::vector<float>",&deepSuperClusterLWP_e2x5BottomRatio); 
         tree->Branch("deepSuperClusterLWP_e2x5LeftRatio","std::vector<float>",&deepSuperClusterLWP_e2x5LeftRatio); 
         tree->Branch("deepSuperClusterLWP_swissCross","std::vector<float>",&deepSuperClusterLWP_swissCross); 
         tree->Branch("deepSuperClusterLWP_r9","std::vector<float>",&deepSuperClusterLWP_r9);
         tree->Branch("deepSuperClusterLWP_sigmaIetaIeta","std::vector<float>",&deepSuperClusterLWP_sigmaIetaIeta);
         tree->Branch("deepSuperClusterLWP_sigmaIetaIphi","std::vector<float>",&deepSuperClusterLWP_sigmaIetaIphi);
         tree->Branch("deepSuperClusterLWP_sigmaIphiIphi","std::vector<float>",&deepSuperClusterLWP_sigmaIphiIphi);
         tree->Branch("deepSuperClusterLWP_full5x5_e5x5","std::vector<float>",&deepSuperClusterLWP_full5x5_e5x5);
         tree->Branch("deepSuperClusterLWP_full5x5_e2x2Ratio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2x2Ratio);
         tree->Branch("deepSuperClusterLWP_full5x5_e3x3Ratio","std::vector<float>",&deepSuperClusterLWP_full5x5_e3x3Ratio);
         tree->Branch("deepSuperClusterLWP_full5x5_eMaxRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_eMaxRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_e2ndRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2ndRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_eTopRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_eTopRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_eRightRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_eRightRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_eBottomRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_eBottomRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_eLeftRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_eLeftRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_e2x5MaxRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2x5MaxRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_e2x5TopRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2x5TopRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_e2x5RightRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2x5RightRatio);
         tree->Branch("deepSuperClusterLWP_full5x5_e2x5BottomRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2x5BottomRatio); 
         tree->Branch("deepSuperClusterLWP_full5x5_e2x5LeftRatio","std::vector<float>",&deepSuperClusterLWP_full5x5_e2x5LeftRatio); 
         tree->Branch("deepSuperClusterLWP_full5x5_swissCross","std::vector<float>",&deepSuperClusterLWP_full5x5_swissCross); 
         tree->Branch("deepSuperClusterLWP_full5x5_r9","std::vector<float>",&deepSuperClusterLWP_full5x5_r9);
         tree->Branch("deepSuperClusterLWP_full5x5_sigmaIetaIeta","std::vector<float>",&deepSuperClusterLWP_full5x5_sigmaIetaIeta);
         tree->Branch("deepSuperClusterLWP_full5x5_sigmaIetaIphi","std::vector<float>",&deepSuperClusterLWP_full5x5_sigmaIetaIphi);
         tree->Branch("deepSuperClusterLWP_full5x5_sigmaIphiIphi","std::vector<float>",&deepSuperClusterLWP_full5x5_sigmaIphiIphi);
         tree->Branch("deepSuperClusterTWP_e5x5","std::vector<float>",&deepSuperClusterTWP_e5x5);
         tree->Branch("deepSuperClusterTWP_e2x2Ratio","std::vector<float>",&deepSuperClusterTWP_e2x2Ratio);
         tree->Branch("deepSuperClusterTWP_e3x3Ratio","std::vector<float>",&deepSuperClusterTWP_e3x3Ratio);
         tree->Branch("deepSuperClusterTWP_eMaxRatio","std::vector<float>",&deepSuperClusterTWP_eMaxRatio);
         tree->Branch("deepSuperClusterTWP_e2ndRatio","std::vector<float>",&deepSuperClusterTWP_e2ndRatio);
         tree->Branch("deepSuperClusterTWP_eTopRatio","std::vector<float>",&deepSuperClusterTWP_eTopRatio);
         tree->Branch("deepSuperClusterTWP_eRightRatio","std::vector<float>",&deepSuperClusterTWP_eRightRatio);
         tree->Branch("deepSuperClusterTWP_eBottomRatio","std::vector<float>",&deepSuperClusterTWP_eBottomRatio);
         tree->Branch("deepSuperClusterTWP_eLeftRatio","std::vector<float>",&deepSuperClusterTWP_eLeftRatio);
         tree->Branch("deepSuperClusterTWP_e2x5MaxRatio","std::vector<float>",&deepSuperClusterTWP_e2x5MaxRatio);
         tree->Branch("deepSuperClusterTWP_e2x5TopRatio","std::vector<float>",&deepSuperClusterTWP_e2x5TopRatio);
         tree->Branch("deepSuperClusterTWP_e2x5RightRatio","std::vector<float>",&deepSuperClusterTWP_e2x5RightRatio);
         tree->Branch("deepSuperClusterTWP_e2x5BottomRatio","std::vector<float>",&deepSuperClusterTWP_e2x5BottomRatio); 
         tree->Branch("deepSuperClusterTWP_e2x5LeftRatio","std::vector<float>",&deepSuperClusterTWP_e2x5LeftRatio); 
         tree->Branch("deepSuperClusterTWP_swissCross","std::vector<float>",&deepSuperClusterTWP_swissCross); 
         tree->Branch("deepSuperClusterTWP_r9","std::vector<float>",&deepSuperClusterTWP_r9);
         tree->Branch("deepSuperClusterTWP_sigmaIetaIeta","std::vector<float>",&deepSuperClusterTWP_sigmaIetaIeta);
         tree->Branch("deepSuperClusterTWP_sigmaIetaIphi","std::vector<float>",&deepSuperClusterTWP_sigmaIetaIphi);
         tree->Branch("deepSuperClusterTWP_sigmaIphiIphi","std::vector<float>",&deepSuperClusterTWP_sigmaIphiIphi);
         tree->Branch("deepSuperClusterTWP_full5x5_e5x5","std::vector<float>",&deepSuperClusterTWP_full5x5_e5x5);
         tree->Branch("deepSuperClusterTWP_full5x5_e2x2Ratio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2x2Ratio);
         tree->Branch("deepSuperClusterTWP_full5x5_e3x3Ratio","std::vector<float>",&deepSuperClusterTWP_full5x5_e3x3Ratio);
         tree->Branch("deepSuperClusterTWP_full5x5_eMaxRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_eMaxRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_e2ndRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2ndRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_eTopRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_eTopRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_eRightRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_eRightRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_eBottomRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_eBottomRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_eLeftRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_eLeftRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_e2x5MaxRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2x5MaxRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_e2x5TopRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2x5TopRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_e2x5RightRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2x5RightRatio);
         tree->Branch("deepSuperClusterTWP_full5x5_e2x5BottomRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2x5BottomRatio); 
         tree->Branch("deepSuperClusterTWP_full5x5_e2x5LeftRatio","std::vector<float>",&deepSuperClusterTWP_full5x5_e2x5LeftRatio); 
         tree->Branch("deepSuperClusterTWP_full5x5_swissCross","std::vector<float>",&deepSuperClusterTWP_full5x5_swissCross); 
         tree->Branch("deepSuperClusterTWP_full5x5_r9","std::vector<float>",&deepSuperClusterTWP_full5x5_r9);
         tree->Branch("deepSuperClusterTWP_full5x5_sigmaIetaIeta","std::vector<float>",&deepSuperClusterTWP_full5x5_sigmaIetaIeta);
         tree->Branch("deepSuperClusterTWP_full5x5_sigmaIetaIphi","std::vector<float>",&deepSuperClusterTWP_full5x5_sigmaIetaIphi);
         tree->Branch("deepSuperClusterTWP_full5x5_sigmaIphiIphi","std::vector<float>",&deepSuperClusterTWP_full5x5_sigmaIphiIphi);
         if(useHcalTowers_ ){
            tree->Branch("deepSuperCluster_HoEraw","std::vector<float>",&deepSuperCluster_HoEraw);
            tree->Branch("deepSuperCluster_HoErawBC","std::vector<float>",&deepSuperCluster_HoErawBC);
            tree->Branch("deepSuperClusterLWP_HoEraw","std::vector<float>",&deepSuperClusterLWP_HoEraw);
            tree->Branch("deepSuperClusterLWP_HoErawBC","std::vector<float>",&deepSuperClusterLWP_HoErawBC);
            tree->Branch("deepSuperClusterTWP_HoEraw","std::vector<float>",&deepSuperClusterTWP_HoEraw);
            tree->Branch("deepSuperClusterTWP_HoErawBC","std::vector<float>",&deepSuperClusterTWP_HoErawBC);
         }  
      } 
   }
}

RecoSimDumper::~RecoSimDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void RecoSimDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

   //calo geometry
   edm::ESHandle<CaloGeometry> caloGeometry;
   iSetup.get<CaloGeometryRecord>().get(caloGeometry);
   const CaloGeometry *geometry = caloGeometry.product();
   _ebGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
   _eeGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
   _esGeom = caloGeometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
   if (_esGeom) {
    for (uint32_t ic = 0; ic < _esGeom->getValidDetIds().size() && (!_esPlus || !_esMinus); ++ic) {
      const double z = _esGeom->getGeometry(_esGeom->getValidDetIds()[ic])->getPosition().z();
      _esPlus = _esPlus || (0 < z);
      _esMinus = _esMinus || (0 > z);
    }
   }

   edm::ESHandle<CaloTopology> caloTopology;
   iSetup.get<CaloTopologyRecord>().get(caloTopology);
   const CaloTopology* topology = caloTopology.product();
   
   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   ev.getByToken(genToken_,genParticles);
   if (!genParticles.isValid()) {
       std::cerr << "Analyze --> genParticles not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<CaloParticle> > caloParticles;
   ev.getByToken(caloPartToken_,caloParticles);
   if (!caloParticles.isValid()) {
       std::cerr << "Analyze --> caloParticles not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if(saveRechits_) {
      if (!recHitsEB.isValid()) {
          std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
          return;
      }
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if(saveRechits_) {
      if (!recHitsEE.isValid()) {
          std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFRecHit> > pfRecHits;
   ev.getByToken(pfRecHitToken_, pfRecHits);
   if(savePFRechits_) {
      if (!pfRecHits.isValid()) {
          std::cerr << "Analyze --> pfRecHits not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFCluster> > pfClusters;
   ev.getByToken(pfClusterToken_, pfClusters);
   if(savePFCluster_) {
      if (!pfClusters.isValid()) {
          std::cerr << "Analyze --> pfClusters not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   ev.getByToken(ebSuperClusterToken_, superClusterEB);
   if(saveSuperCluster_) {
      if (!superClusterEB.isValid()) {
          std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   ev.getByToken(eeSuperClusterToken_, superClusterEE);
   if(saveSuperCluster_) {
      if (!superClusterEE.isValid()) {
          std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<reco::SuperCluster> > retunedSuperClusterEB;
   ev.getByToken(ebRetunedSuperClusterToken_, retunedSuperClusterEB);
   if(saveSuperCluster_ && useRetunedSC_) {
      if (!retunedSuperClusterEB.isValid()) {
          std::cerr << "Analyze --> retunedSuperClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > retunedSuperClusterEE;
   ev.getByToken(eeRetunedSuperClusterToken_, retunedSuperClusterEE);
   if(saveSuperCluster_ && useRetunedSC_) {
      if (!retunedSuperClusterEE.isValid()) {
          std::cerr << "Analyze --> retunedSuperClusterEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterEB;
   ev.getByToken(ebDeepSuperClusterToken_, deepSuperClusterEB);
   if(saveSuperCluster_) {
      if (!deepSuperClusterEB.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterEE;
   ev.getByToken(eeDeepSuperClusterToken_, deepSuperClusterEE);
   if(saveSuperCluster_) {
      if (!deepSuperClusterEE.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterLWPEB;
   ev.getByToken(ebDeepSuperClusterLWPToken_, deepSuperClusterLWPEB);
   if(saveSuperCluster_) {
      if (!deepSuperClusterLWPEB.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterLWPEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterLWPEE;
   ev.getByToken(eeDeepSuperClusterLWPToken_, deepSuperClusterLWPEE);
   if(saveSuperCluster_) {
      if (!deepSuperClusterLWPEE.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterLWPEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterTWPEB;
   ev.getByToken(ebDeepSuperClusterTWPToken_, deepSuperClusterTWPEB);
   if(saveSuperCluster_) {
      if (!deepSuperClusterTWPEB.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterTWPEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterTWPEE;
   ev.getByToken(eeDeepSuperClusterTWPToken_, deepSuperClusterTWPEE);
   if(saveSuperCluster_) {
      if (!deepSuperClusterTWPEE.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterTWPEE not found" << std::endl; 
          return;
      }
   }  

   //compute EgammaTowers;
   Handle<CaloTowerCollection> hcalTowers;
   ev.getByToken(hcalTowersToken_, hcalTowers);
   if(useHcalTowers_){
      if (!hcalTowers.isValid()) {
          std::cerr << "Analyze --> hcalTowers not found" << std::endl; 
          return;
      } 
      //hcalTowersColl = hcalTowers.product();
      towerIso1_ = new EgammaTowerIsolation(0.15, 0., 0., 1, hcalTowers.product());
      towerIso2_ = new EgammaTowerIsolation(0.15, 0., 0., 2, hcalTowers.product());
      //egammaHadTower_ = new EgammaHadTower(iSetup);
      //egammaHadTower_->setTowerCollection(hcalTowers.product()); 
   } 

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();

   nVtx = vertices->size();
   rho = *(rhos.product());

   genParticle_id.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();
   std::vector<GenParticle> genParts;
   for(const auto& iGen : *(genParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++)
           if((iGen.pdgId()==genID_.at(id) || genID_.at(id)==0) && iGen.status()==1) isGoodParticle=true;
      
       if(!isGoodParticle) continue; 
       genParticle_id.push_back(iGen.pdgId()); 
       genParticle_energy.push_back(iGen.energy()); 
       genParticle_pt.push_back(iGen.pt());
       genParticle_eta.push_back(iGen.eta());
       genParticle_phi.push_back(iGen.phi());
       
       genParts.push_back(iGen); 
   } 
   
   int nGenParticles = genParts.size(); 
   //std::cout << "GenParticles size  : " << nGenParticles << std::endl;
   
   std::vector<CaloParticle> caloParts;
   for(const auto& iCalo : *(caloParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++) 
           if(iCalo.pdgId()==genID_.at(id) || genID_.at(id)==0) isGoodParticle=true;

       if(!isGoodParticle) continue;     

       caloParts.push_back(iCalo); 
   }

   int nCaloParticles = caloParts.size(); 
   //std::cout << "CaloParticles size  : " << nCaloParticles << std::endl;
  
   genParticle_pfCluster_dR_genScore_MatchedIndex.clear();
   genParticle_superCluster_dR_genScore_MatchedIndex.clear();
   genParticle_retunedSuperCluster_dR_genScore_MatchedIndex.clear();
   genParticle_deepSuperCluster_dR_genScore_MatchedIndex.clear();
   genParticle_deepSuperClusterLWP_dR_genScore_MatchedIndex.clear();
   genParticle_deepSuperClusterTWP_dR_genScore_MatchedIndex.clear();
   genParticle_pfCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_superCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_retunedSuperCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_deepSuperCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_deepSuperClusterLWP_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_deepSuperClusterTWP_dR_genScore_MatchedIndex.resize(nGenParticles);

   caloParticle_id.clear();
   caloParticle_genEnergy.clear();
   caloParticle_simEnergy.clear();
   caloParticle_genPt.clear();
   caloParticle_simPt.clear();
   caloParticle_genEta.clear();
   caloParticle_simEta.clear();
   caloParticle_genPhi.clear();
   caloParticle_simPhi.clear();
   caloParticle_simIeta.clear();
   caloParticle_simIphi.clear();
   caloParticle_simIz.clear();
   caloParticle_pfCluster_dR_simScore_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_old_MatchedIndex.clear();
   caloParticle_pfCluster_simScore_MatchedIndex.clear();
   caloParticle_pfCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex.clear();  
   caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex.clear();  
   caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex.clear();  
   caloParticle_pfCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex.clear();
   caloParticle_superCluster_dR_simScore_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_old_MatchedIndex.clear();
   caloParticle_superCluster_simScore_MatchedIndex.clear();
   caloParticle_superCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_1MeVCut_MatchedIndex.clear(); 
   caloParticle_superCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_100MeVCut_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_500MeVCut_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_1GeVCut_MatchedIndex.clear();
   caloParticle_superCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_superCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_superCluster_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_superCluster_rechit_sim_combined_fraction_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_simScore_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex.clear();  
   caloParticle_deepSuperCluster_dR_simScore_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_old_MatchedIndex.clear();
   caloParticle_deepSuperCluster_simScore_MatchedIndex.clear();
   caloParticle_deepSuperCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_deepSuperCluster_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_dR_simScore_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_old_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_simScore_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_n_shared_xtals_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex.clear(); 
   caloParticle_deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_dR_simScore_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_old_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_simScore_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_n_shared_xtals_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex.clear(); 
   caloParticle_deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex.clear();
   caloParticle_deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles);  
   caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles);  
   caloParticle_superCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles); 
   caloParticle_superCluster_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles); 
   caloParticle_retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperClusterLWP_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperClusterTWP_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_old_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex.resize(nCaloParticles);
   
   simHit_energy.clear();
   simHit_eta.clear();
   simHit_phi.clear();
   simHit_ieta.clear();
   simHit_iphi.clear();
   simHit_iz.clear();
   simHit_energy.resize(nCaloParticles);
   simHit_eta.resize(nCaloParticles);
   simHit_phi.resize(nCaloParticles);
   simHit_ieta.resize(nCaloParticles);
   simHit_iphi.resize(nCaloParticles);
   simHit_iz.resize(nCaloParticles);

   recHit_noPF_energy.clear();
   recHit_noPF_eta.clear();
   recHit_noPF_phi.clear();
   recHit_noPF_ieta.clear();
   recHit_noPF_iphi.clear();
   recHit_noPF_iz.clear();  

   pfRecHit_unClustered_energy.clear();
   pfRecHit_unClustered_eta.clear();
   pfRecHit_unClustered_phi.clear();
   pfRecHit_unClustered_ieta.clear();
   pfRecHit_unClustered_iphi.clear();
   pfRecHit_unClustered_iz.clear();

   int nPFClusters = (pfClusters.product())->size();
   pfCluster_energy.clear();
   pfCluster_eta.clear();
   pfCluster_phi.clear();
   pfCluster_ieta.clear();
   pfCluster_iphi.clear();
   pfCluster_iz.clear();
   pfCluster_superClustersIndex.clear();
   pfCluster_retunedSuperClustersIndex.clear();
   pfCluster_deepSuperClustersIndex.clear();
   pfCluster_deepSuperClusterLWPsIndex.clear();
   pfCluster_deepSuperClusterTWPsIndex.clear(); 
   pfCluster_e5x5.clear();
   pfCluster_e2x2Ratio.clear();
   pfCluster_e3x3Ratio.clear();
   pfCluster_eMaxRatio.clear();
   pfCluster_e2ndRatio.clear();
   pfCluster_eTopRatio.clear();
   pfCluster_eRightRatio.clear();
   pfCluster_eBottomRatio.clear();
   pfCluster_eLeftRatio.clear();
   pfCluster_e2x5MaxRatio.clear();
   pfCluster_e2x5TopRatio.clear();
   pfCluster_e2x5RightRatio.clear();
   pfCluster_e2x5BottomRatio.clear();
   pfCluster_e2x5LeftRatio.clear();
   pfCluster_swissCross.clear();
   pfCluster_r9.clear();
   pfCluster_sigmaIetaIeta.clear(); 
   pfCluster_sigmaIetaIphi.clear(); 
   pfCluster_sigmaIphiIphi.clear(); 
   pfCluster_full5x5_e5x5.clear();
   pfCluster_full5x5_e2x2Ratio.clear();
   pfCluster_full5x5_e3x3Ratio.clear();
   pfCluster_full5x5_eMaxRatio.clear();
   pfCluster_full5x5_e2ndRatio.clear();
   pfCluster_full5x5_eTopRatio.clear();
   pfCluster_full5x5_eRightRatio.clear();
   pfCluster_full5x5_eBottomRatio.clear();
   pfCluster_full5x5_eLeftRatio.clear();
   pfCluster_full5x5_e2x5MaxRatio.clear();
   pfCluster_full5x5_e2x5TopRatio.clear();
   pfCluster_full5x5_e2x5RightRatio.clear();
   pfCluster_full5x5_e2x5BottomRatio.clear();
   pfCluster_full5x5_e2x5LeftRatio.clear();
   pfCluster_full5x5_swissCross.clear();
   pfCluster_full5x5_r9.clear();
   pfCluster_full5x5_sigmaIetaIeta.clear(); 
   pfCluster_full5x5_sigmaIetaIphi.clear(); 
   pfCluster_full5x5_sigmaIphiIphi.clear();    
   pfCluster_dR_genScore_MatchedIndex.clear();
   pfCluster_dR_simScore_MatchedIndex.clear();
   pfCluster_sim_fraction_old_MatchedIndex.clear();
   pfCluster_simScore_MatchedIndex.clear();
   pfCluster_n_shared_xtals_MatchedIndex.clear();
   pfCluster_sim_fraction_MatchedIndex.clear();
   pfCluster_sim_fraction_1MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_5MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_10MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_50MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_100MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_500MeVCut_MatchedIndex.clear();
   pfCluster_sim_fraction_1GeVCut_MatchedIndex.clear();
   pfCluster_sim_rechit_diff_MatchedIndex.clear();
   pfCluster_sim_rechit_fraction_MatchedIndex.clear();
   pfCluster_global_sim_rechit_fraction_MatchedIndex.clear();  
   pfCluster_hgcal_caloToCluster_MatchedIndex.clear();  
   pfCluster_hgcal_clusterToCalo_MatchedIndex.clear();   
   pfCluster_sim_rechit_combined_fraction_MatchedIndex.clear();  
   pfCluster_rechit_sim_combined_fraction_MatchedIndex.clear();   
   pfCluster_dR_genScore.clear();
   pfCluster_dR_simScore.clear();
   pfCluster_sim_fraction_old.clear();
   pfCluster_simScore.clear();
   pfCluster_n_shared_xtals.clear();
   pfCluster_sim_fraction.clear();
   pfCluster_sim_fraction_1MeVCut.clear();
   pfCluster_sim_fraction_5MeVCut.clear();
   pfCluster_sim_fraction_10MeVCut.clear();
   pfCluster_sim_fraction_50MeVCut.clear();
   pfCluster_sim_fraction_100MeVCut.clear();
   pfCluster_sim_fraction_500MeVCut.clear();
   pfCluster_sim_fraction_1GeVCut.clear();
   pfCluster_sim_rechit_diff.clear();
   pfCluster_sim_rechit_fraction.clear();
   pfCluster_global_sim_rechit_fraction.clear();  
   pfCluster_hgcal_caloToCluster.clear();  
   pfCluster_hgcal_clusterToCalo.clear();    
   pfCluster_sim_rechit_combined_fraction.clear();  
   pfCluster_rechit_sim_combined_fraction.clear();    
   pfCluster_dR_genScore.resize(nPFClusters);
   pfCluster_dR_simScore.resize(nPFClusters);
   pfCluster_sim_fraction_old.resize(nPFClusters);
   pfCluster_simScore.resize(nPFClusters);
   pfCluster_n_shared_xtals.resize(nPFClusters);
   pfCluster_sim_fraction.resize(nPFClusters);
   pfCluster_sim_fraction_1MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_5MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_10MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_50MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_100MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_500MeVCut.resize(nPFClusters);
   pfCluster_sim_fraction_1GeVCut.resize(nPFClusters);
   pfCluster_sim_rechit_diff.resize(nPFClusters);
   pfCluster_sim_rechit_fraction.resize(nPFClusters);
   pfCluster_global_sim_rechit_fraction.resize(nPFClusters); 
   pfCluster_hgcal_caloToCluster.resize(nPFClusters);   
   pfCluster_hgcal_clusterToCalo.resize(nPFClusters);  
   pfCluster_sim_rechit_combined_fraction.resize(nPFClusters);   
   pfCluster_rechit_sim_combined_fraction.resize(nPFClusters);  
   pfCluster_superClustersIndex.resize(nPFClusters); 
   pfCluster_retunedSuperClustersIndex.resize(nPFClusters);  
   pfCluster_deepSuperClustersIndex.resize(nPFClusters); 
   pfCluster_deepSuperClusterLWPsIndex.resize(nPFClusters); 
   pfCluster_deepSuperClusterTWPsIndex.resize(nPFClusters); 

   pfClusterHit_energy.clear();
   pfClusterHit_rechitEnergy.clear();
   pfClusterHit_eta.clear();
   pfClusterHit_phi.clear();   
   pfClusterHit_ieta.clear();
   pfClusterHit_iphi.clear(); 
   pfClusterHit_iz.clear();           
   pfClusterHit_energy.resize(nPFClusters);  
   pfClusterHit_rechitEnergy.resize(nPFClusters);  
   pfClusterHit_eta.resize(nPFClusters);  
   pfClusterHit_phi.resize(nPFClusters);     
   pfClusterHit_ieta.resize(nPFClusters);  
   pfClusterHit_iphi.resize(nPFClusters);   
   pfClusterHit_iz.resize(nPFClusters);   
  
   superCluster_energy.clear(); 
   superCluster_eta.clear(); 
   superCluster_phi.clear();  
   superCluster_etaWidth.clear();     
   superCluster_phiWidth.clear(); 
   superCluster_R.clear(); 
   superCluster_nPFClusters.clear(); 
   superCluster_ieta.clear(); 
   superCluster_iphi.clear();
   superCluster_iz.clear(); 
   superCluster_seedIndex.clear(); 
   superCluster_pfClustersIndex.clear(); 
   superCluster_e5x5.clear();
   superCluster_e2x2Ratio.clear();
   superCluster_e3x3Ratio.clear();
   superCluster_eMaxRatio.clear();
   superCluster_e2ndRatio.clear();
   superCluster_eTopRatio.clear();
   superCluster_eRightRatio.clear();
   superCluster_eBottomRatio.clear();
   superCluster_eLeftRatio.clear();
   superCluster_e2x5MaxRatio.clear();
   superCluster_e2x5TopRatio.clear();
   superCluster_e2x5RightRatio.clear();
   superCluster_e2x5BottomRatio.clear();
   superCluster_e2x5LeftRatio.clear();
   superCluster_swissCross.clear();
   superCluster_r9.clear();
   superCluster_sigmaIetaIeta.clear(); 
   superCluster_sigmaIetaIphi.clear(); 
   superCluster_sigmaIphiIphi.clear(); 
   superCluster_full5x5_e5x5.clear();
   superCluster_full5x5_e2x2Ratio.clear();
   superCluster_full5x5_e3x3Ratio.clear();
   superCluster_full5x5_eMaxRatio.clear();
   superCluster_full5x5_e2ndRatio.clear();
   superCluster_full5x5_eTopRatio.clear();
   superCluster_full5x5_eRightRatio.clear();
   superCluster_full5x5_eBottomRatio.clear();
   superCluster_full5x5_eLeftRatio.clear();
   superCluster_full5x5_e2x5MaxRatio.clear();
   superCluster_full5x5_e2x5TopRatio.clear();
   superCluster_full5x5_e2x5RightRatio.clear();
   superCluster_full5x5_e2x5BottomRatio.clear();
   superCluster_full5x5_e2x5LeftRatio.clear();
   superCluster_full5x5_swissCross.clear();
   superCluster_full5x5_r9.clear();
   superCluster_full5x5_sigmaIetaIeta.clear(); 
   superCluster_full5x5_sigmaIetaIphi.clear(); 
   superCluster_full5x5_sigmaIphiIphi.clear();    
   superCluster_dR_genScore_MatchedIndex.clear();  
   superCluster_dR_simScore_MatchedIndex.clear();  
   superCluster_sim_fraction_old_MatchedIndex.clear();  
   superCluster_simScore_MatchedIndex.clear();  
   superCluster_n_shared_xtals_MatchedIndex.clear();   
   superCluster_sim_fraction_MatchedIndex.clear();  
   superCluster_sim_fraction_1MeVCut_MatchedIndex.clear();  
   superCluster_sim_fraction_5MeVCut_MatchedIndex.clear();  
   superCluster_sim_fraction_10MeVCut_MatchedIndex.clear();  
   superCluster_sim_fraction_50MeVCut_MatchedIndex.clear();  
   superCluster_sim_fraction_100MeVCut_MatchedIndex.clear();  
   superCluster_sim_fraction_500MeVCut_MatchedIndex.clear();  
   superCluster_sim_fraction_1GeVCut_MatchedIndex.clear();  
   superCluster_sim_rechit_diff_MatchedIndex.clear();  
   superCluster_sim_rechit_fraction_MatchedIndex.clear();  
   superCluster_global_sim_rechit_fraction_MatchedIndex.clear();  
   superCluster_hgcal_caloToCluster_MatchedIndex.clear();  
   superCluster_hgcal_clusterToCalo_MatchedIndex.clear();   
   superCluster_sim_rechit_combined_fraction_MatchedIndex.clear();  
   superCluster_rechit_sim_combined_fraction_MatchedIndex.clear();   
   superCluster_dR_genScore.clear();
   superCluster_dR_simScore.clear();
   superCluster_sim_fraction_old.clear();
   superCluster_simScore.clear();
   superCluster_n_shared_xtals.clear();
   superCluster_sim_fraction.clear();
   superCluster_sim_fraction_1MeVCut.clear();
   superCluster_sim_fraction_5MeVCut.clear();
   superCluster_sim_fraction_10MeVCut.clear();
   superCluster_sim_fraction_50MeVCut.clear();
   superCluster_sim_fraction_100MeVCut.clear();
   superCluster_sim_fraction_500MeVCut.clear();
   superCluster_sim_fraction_1GeVCut.clear();
   superCluster_sim_rechit_diff.clear();
   superCluster_sim_rechit_fraction.clear();
   superCluster_global_sim_rechit_fraction.clear();  
   superCluster_hgcal_caloToCluster.clear();  
   superCluster_hgcal_clusterToCalo.clear();   
   superCluster_sim_rechit_combined_fraction.clear();  
   superCluster_rechit_sim_combined_fraction.clear();   
   superCluster_psCluster_energy.clear();
   superCluster_psCluster_eta.clear();
   superCluster_psCluster_phi.clear();
   superCluster_HoEraw.clear(); 
   superCluster_HoErawBC.clear(); 
   retunedSuperCluster_energy.clear(); 
   retunedSuperCluster_eta.clear(); 
   retunedSuperCluster_phi.clear();  
   retunedSuperCluster_etaWidth.clear();     
   retunedSuperCluster_phiWidth.clear(); 
   retunedSuperCluster_R.clear(); 
   retunedSuperCluster_nPFClusters.clear(); 
   retunedSuperCluster_ieta.clear(); 
   retunedSuperCluster_iphi.clear();
   retunedSuperCluster_iz.clear(); 
   retunedSuperCluster_seedIndex.clear(); 
   retunedSuperCluster_pfClustersIndex.clear(); 
   retunedSuperCluster_e5x5.clear();
   retunedSuperCluster_e2x2Ratio.clear();
   retunedSuperCluster_e3x3Ratio.clear();
   retunedSuperCluster_eMaxRatio.clear();
   retunedSuperCluster_e2ndRatio.clear();
   retunedSuperCluster_eTopRatio.clear();
   retunedSuperCluster_eRightRatio.clear();
   retunedSuperCluster_eBottomRatio.clear();
   retunedSuperCluster_eLeftRatio.clear();
   retunedSuperCluster_e2x5MaxRatio.clear();
   retunedSuperCluster_e2x5TopRatio.clear();
   retunedSuperCluster_e2x5RightRatio.clear();
   retunedSuperCluster_e2x5BottomRatio.clear();
   retunedSuperCluster_e2x5LeftRatio.clear();
   retunedSuperCluster_swissCross.clear();
   retunedSuperCluster_r9.clear();
   retunedSuperCluster_sigmaIetaIeta.clear(); 
   retunedSuperCluster_sigmaIetaIphi.clear(); 
   retunedSuperCluster_sigmaIphiIphi.clear(); 
   retunedSuperCluster_full5x5_e5x5.clear();
   retunedSuperCluster_full5x5_e2x2Ratio.clear();
   retunedSuperCluster_full5x5_e3x3Ratio.clear();
   retunedSuperCluster_full5x5_eMaxRatio.clear();
   retunedSuperCluster_full5x5_e2ndRatio.clear();
   retunedSuperCluster_full5x5_eTopRatio.clear();
   retunedSuperCluster_full5x5_eRightRatio.clear();
   retunedSuperCluster_full5x5_eBottomRatio.clear();
   retunedSuperCluster_full5x5_eLeftRatio.clear();
   retunedSuperCluster_full5x5_e2x5MaxRatio.clear();
   retunedSuperCluster_full5x5_e2x5TopRatio.clear();
   retunedSuperCluster_full5x5_e2x5RightRatio.clear();
   retunedSuperCluster_full5x5_e2x5BottomRatio.clear();
   retunedSuperCluster_full5x5_e2x5LeftRatio.clear();
   retunedSuperCluster_full5x5_swissCross.clear();
   retunedSuperCluster_full5x5_r9.clear();
   retunedSuperCluster_full5x5_sigmaIetaIeta.clear(); 
   retunedSuperCluster_full5x5_sigmaIetaIphi.clear(); 
   retunedSuperCluster_full5x5_sigmaIphiIphi.clear();    
   retunedSuperCluster_dR_genScore_MatchedIndex.clear();  
   retunedSuperCluster_dR_simScore_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_old_MatchedIndex.clear();  
   retunedSuperCluster_simScore_MatchedIndex.clear();  
   retunedSuperCluster_n_shared_xtals_MatchedIndex.clear();   
   retunedSuperCluster_sim_fraction_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex.clear();  
   retunedSuperCluster_sim_rechit_diff_MatchedIndex.clear();  
   retunedSuperCluster_sim_rechit_fraction_MatchedIndex.clear();  
   retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex.clear();  
   retunedSuperCluster_hgcal_caloToCluster_MatchedIndex.clear();  
   retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex.clear();   
   retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex.clear();  
   retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex.clear();   
   retunedSuperCluster_dR_genScore.clear();
   retunedSuperCluster_dR_simScore.clear();
   retunedSuperCluster_sim_fraction_old.clear();
   retunedSuperCluster_simScore.clear();
   retunedSuperCluster_n_shared_xtals.clear();
   retunedSuperCluster_sim_fraction.clear();
   retunedSuperCluster_sim_fraction_1MeVCut.clear();
   retunedSuperCluster_sim_fraction_5MeVCut.clear();
   retunedSuperCluster_sim_fraction_10MeVCut.clear();
   retunedSuperCluster_sim_fraction_50MeVCut.clear();
   retunedSuperCluster_sim_fraction_100MeVCut.clear();
   retunedSuperCluster_sim_fraction_500MeVCut.clear();
   retunedSuperCluster_sim_fraction_1GeVCut.clear();
   retunedSuperCluster_sim_rechit_diff.clear();
   retunedSuperCluster_sim_rechit_fraction.clear();
   retunedSuperCluster_global_sim_rechit_fraction.clear();  
   retunedSuperCluster_hgcal_caloToCluster.clear();  
   retunedSuperCluster_hgcal_clusterToCalo.clear();   
   retunedSuperCluster_sim_rechit_combined_fraction.clear();  
   retunedSuperCluster_rechit_sim_combined_fraction.clear();   
   retunedSuperCluster_psCluster_energy.clear();
   retunedSuperCluster_psCluster_eta.clear();
   retunedSuperCluster_psCluster_phi.clear();
   retunedSuperCluster_HoEraw.clear(); 
   retunedSuperCluster_HoErawBC.clear(); 
   deepSuperCluster_energy.clear(); 
   deepSuperCluster_eta.clear(); 
   deepSuperCluster_phi.clear();  
   deepSuperCluster_etaWidth.clear();     
   deepSuperCluster_phiWidth.clear(); 
   deepSuperCluster_R.clear(); 
   deepSuperCluster_nPFClusters.clear(); 
   deepSuperCluster_ieta.clear(); 
   deepSuperCluster_iphi.clear();
   deepSuperCluster_iz.clear(); 
   deepSuperCluster_seedIndex.clear(); 
   deepSuperCluster_pfClustersIndex.clear(); 
   deepSuperCluster_e5x5.clear();
   deepSuperCluster_e2x2Ratio.clear();
   deepSuperCluster_e3x3Ratio.clear();
   deepSuperCluster_eMaxRatio.clear();
   deepSuperCluster_e2ndRatio.clear();
   deepSuperCluster_eTopRatio.clear();
   deepSuperCluster_eRightRatio.clear();
   deepSuperCluster_eBottomRatio.clear();
   deepSuperCluster_eLeftRatio.clear();
   deepSuperCluster_e2x5MaxRatio.clear();
   deepSuperCluster_e2x5TopRatio.clear();
   deepSuperCluster_e2x5RightRatio.clear();
   deepSuperCluster_e2x5BottomRatio.clear();
   deepSuperCluster_e2x5LeftRatio.clear();
   deepSuperCluster_swissCross.clear();
   deepSuperCluster_r9.clear();
   deepSuperCluster_sigmaIetaIeta.clear(); 
   deepSuperCluster_sigmaIetaIphi.clear(); 
   deepSuperCluster_sigmaIphiIphi.clear(); 
   deepSuperCluster_full5x5_e5x5.clear();
   deepSuperCluster_full5x5_e2x2Ratio.clear();
   deepSuperCluster_full5x5_e3x3Ratio.clear();
   deepSuperCluster_full5x5_eMaxRatio.clear();
   deepSuperCluster_full5x5_e2ndRatio.clear();
   deepSuperCluster_full5x5_eTopRatio.clear();
   deepSuperCluster_full5x5_eRightRatio.clear();
   deepSuperCluster_full5x5_eBottomRatio.clear();
   deepSuperCluster_full5x5_eLeftRatio.clear();
   deepSuperCluster_full5x5_e2x5MaxRatio.clear();
   deepSuperCluster_full5x5_e2x5TopRatio.clear();
   deepSuperCluster_full5x5_e2x5RightRatio.clear();
   deepSuperCluster_full5x5_e2x5BottomRatio.clear();
   deepSuperCluster_full5x5_e2x5LeftRatio.clear();
   deepSuperCluster_full5x5_swissCross.clear();
   deepSuperCluster_full5x5_r9.clear();
   deepSuperCluster_full5x5_sigmaIetaIeta.clear(); 
   deepSuperCluster_full5x5_sigmaIetaIphi.clear(); 
   deepSuperCluster_full5x5_sigmaIphiIphi.clear();    
   deepSuperCluster_dR_genScore_MatchedIndex.clear();  
   deepSuperCluster_dR_simScore_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_old_MatchedIndex.clear();  
   deepSuperCluster_simScore_MatchedIndex.clear();  
   deepSuperCluster_n_shared_xtals_MatchedIndex.clear();   
   deepSuperCluster_sim_fraction_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex.clear();  
   deepSuperCluster_sim_rechit_diff_MatchedIndex.clear();  
   deepSuperCluster_sim_rechit_fraction_MatchedIndex.clear();  
   deepSuperCluster_global_sim_rechit_fraction_MatchedIndex.clear();  
   deepSuperCluster_hgcal_caloToCluster_MatchedIndex.clear();  
   deepSuperCluster_hgcal_clusterToCalo_MatchedIndex.clear();   
   deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex.clear();  
   deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex.clear();   
   deepSuperCluster_dR_genScore.clear();
   deepSuperCluster_dR_simScore.clear();
   deepSuperCluster_sim_fraction_old.clear();
   deepSuperCluster_simScore.clear();
   deepSuperCluster_n_shared_xtals.clear();
   deepSuperCluster_sim_fraction.clear();
   deepSuperCluster_sim_fraction_1MeVCut.clear();
   deepSuperCluster_sim_fraction_5MeVCut.clear();
   deepSuperCluster_sim_fraction_10MeVCut.clear();
   deepSuperCluster_sim_fraction_50MeVCut.clear();
   deepSuperCluster_sim_fraction_100MeVCut.clear();
   deepSuperCluster_sim_fraction_500MeVCut.clear();
   deepSuperCluster_sim_fraction_1GeVCut.clear();
   deepSuperCluster_sim_rechit_diff.clear();
   deepSuperCluster_sim_rechit_fraction.clear();
   deepSuperCluster_global_sim_rechit_fraction.clear();  
   deepSuperCluster_hgcal_caloToCluster.clear();  
   deepSuperCluster_hgcal_clusterToCalo.clear(); 
   deepSuperCluster_sim_rechit_combined_fraction.clear();  
   deepSuperCluster_rechit_sim_combined_fraction.clear();   
   deepSuperCluster_psCluster_energy.clear();
   deepSuperCluster_psCluster_eta.clear();
   deepSuperCluster_psCluster_phi.clear();
   deepSuperCluster_HoEraw.clear(); 
   deepSuperCluster_HoErawBC.clear(); 
   deepSuperClusterLWP_energy.clear(); 
   deepSuperClusterLWP_eta.clear(); 
   deepSuperClusterLWP_phi.clear();  
   deepSuperClusterLWP_etaWidth.clear();     
   deepSuperClusterLWP_phiWidth.clear(); 
   deepSuperClusterLWP_R.clear(); 
   deepSuperClusterLWP_nPFClusters.clear(); 
   deepSuperClusterLWP_ieta.clear(); 
   deepSuperClusterLWP_iphi.clear();
   deepSuperClusterLWP_iz.clear(); 
   deepSuperClusterLWP_seedIndex.clear(); 
   deepSuperClusterLWP_pfClustersIndex.clear(); 
   deepSuperClusterLWP_e5x5.clear();
   deepSuperClusterLWP_e2x2Ratio.clear();
   deepSuperClusterLWP_e3x3Ratio.clear();
   deepSuperClusterLWP_eMaxRatio.clear();
   deepSuperClusterLWP_e2ndRatio.clear();
   deepSuperClusterLWP_eTopRatio.clear();
   deepSuperClusterLWP_eRightRatio.clear();
   deepSuperClusterLWP_eBottomRatio.clear();
   deepSuperClusterLWP_eLeftRatio.clear();
   deepSuperClusterLWP_e2x5MaxRatio.clear();
   deepSuperClusterLWP_e2x5TopRatio.clear();
   deepSuperClusterLWP_e2x5RightRatio.clear();
   deepSuperClusterLWP_e2x5BottomRatio.clear();
   deepSuperClusterLWP_e2x5LeftRatio.clear();
   deepSuperClusterLWP_swissCross.clear();
   deepSuperClusterLWP_r9.clear();
   deepSuperClusterLWP_sigmaIetaIeta.clear(); 
   deepSuperClusterLWP_sigmaIetaIphi.clear(); 
   deepSuperClusterLWP_sigmaIphiIphi.clear(); 
   deepSuperClusterLWP_full5x5_e5x5.clear();
   deepSuperClusterLWP_full5x5_e2x2Ratio.clear();
   deepSuperClusterLWP_full5x5_e3x3Ratio.clear();
   deepSuperClusterLWP_full5x5_eMaxRatio.clear();
   deepSuperClusterLWP_full5x5_e2ndRatio.clear();
   deepSuperClusterLWP_full5x5_eTopRatio.clear();
   deepSuperClusterLWP_full5x5_eRightRatio.clear();
   deepSuperClusterLWP_full5x5_eBottomRatio.clear();
   deepSuperClusterLWP_full5x5_eLeftRatio.clear();
   deepSuperClusterLWP_full5x5_e2x5MaxRatio.clear();
   deepSuperClusterLWP_full5x5_e2x5TopRatio.clear();
   deepSuperClusterLWP_full5x5_e2x5RightRatio.clear();
   deepSuperClusterLWP_full5x5_e2x5BottomRatio.clear();
   deepSuperClusterLWP_full5x5_e2x5LeftRatio.clear();
   deepSuperClusterLWP_full5x5_swissCross.clear();
   deepSuperClusterLWP_full5x5_r9.clear();
   deepSuperClusterLWP_full5x5_sigmaIetaIeta.clear(); 
   deepSuperClusterLWP_full5x5_sigmaIetaIphi.clear(); 
   deepSuperClusterLWP_full5x5_sigmaIphiIphi.clear();    
   deepSuperClusterLWP_dR_genScore_MatchedIndex.clear();  
   deepSuperClusterLWP_dR_simScore_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_old_MatchedIndex.clear();  
   deepSuperClusterLWP_simScore_MatchedIndex.clear();  
   deepSuperClusterLWP_n_shared_xtals_MatchedIndex.clear();   
   deepSuperClusterLWP_sim_fraction_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_rechit_diff_MatchedIndex.clear();  
   deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex.clear();  
   deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex.clear();  
   deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex.clear();  
   deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex.clear();   
   deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex.clear();  
   deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex.clear();   
   deepSuperClusterLWP_dR_genScore.clear();
   deepSuperClusterLWP_dR_simScore.clear();
   deepSuperClusterLWP_sim_fraction_old.clear();
   deepSuperClusterLWP_simScore.clear();
   deepSuperClusterLWP_n_shared_xtals.clear();
   deepSuperClusterLWP_sim_fraction.clear();
   deepSuperClusterLWP_sim_fraction_1MeVCut.clear();
   deepSuperClusterLWP_sim_fraction_5MeVCut.clear();
   deepSuperClusterLWP_sim_fraction_10MeVCut.clear();
   deepSuperClusterLWP_sim_fraction_50MeVCut.clear();
   deepSuperClusterLWP_sim_fraction_100MeVCut.clear();
   deepSuperClusterLWP_sim_fraction_500MeVCut.clear();
   deepSuperClusterLWP_sim_fraction_1GeVCut.clear();
   deepSuperClusterLWP_sim_rechit_diff.clear();
   deepSuperClusterLWP_sim_rechit_fraction.clear();
   deepSuperClusterLWP_global_sim_rechit_fraction.clear();  
   deepSuperClusterLWP_hgcal_caloToCluster.clear();  
   deepSuperClusterLWP_hgcal_clusterToCalo.clear(); 
   deepSuperClusterLWP_sim_rechit_combined_fraction.clear();  
   deepSuperClusterLWP_rechit_sim_combined_fraction.clear();   
   deepSuperClusterLWP_psCluster_energy.clear();
   deepSuperClusterLWP_psCluster_eta.clear();
   deepSuperClusterLWP_psCluster_phi.clear();
   deepSuperClusterLWP_HoEraw.clear(); 
   deepSuperClusterLWP_HoErawBC.clear(); 
   deepSuperClusterTWP_energy.clear(); 
   deepSuperClusterTWP_eta.clear(); 
   deepSuperClusterTWP_phi.clear();  
   deepSuperClusterTWP_etaWidth.clear();     
   deepSuperClusterTWP_phiWidth.clear(); 
   deepSuperClusterTWP_R.clear(); 
   deepSuperClusterTWP_nPFClusters.clear(); 
   deepSuperClusterTWP_ieta.clear(); 
   deepSuperClusterTWP_iphi.clear();
   deepSuperClusterTWP_iz.clear(); 
   deepSuperClusterTWP_seedIndex.clear(); 
   deepSuperClusterTWP_pfClustersIndex.clear(); 
   deepSuperClusterTWP_e5x5.clear();
   deepSuperClusterTWP_e2x2Ratio.clear();
   deepSuperClusterTWP_e3x3Ratio.clear();
   deepSuperClusterTWP_eMaxRatio.clear();
   deepSuperClusterTWP_e2ndRatio.clear();
   deepSuperClusterTWP_eTopRatio.clear();
   deepSuperClusterTWP_eRightRatio.clear();
   deepSuperClusterTWP_eBottomRatio.clear();
   deepSuperClusterTWP_eLeftRatio.clear();
   deepSuperClusterTWP_e2x5MaxRatio.clear();
   deepSuperClusterTWP_e2x5TopRatio.clear();
   deepSuperClusterTWP_e2x5RightRatio.clear();
   deepSuperClusterTWP_e2x5BottomRatio.clear();
   deepSuperClusterTWP_e2x5LeftRatio.clear();
   deepSuperClusterTWP_swissCross.clear();
   deepSuperClusterTWP_r9.clear();
   deepSuperClusterTWP_sigmaIetaIeta.clear(); 
   deepSuperClusterTWP_sigmaIetaIphi.clear(); 
   deepSuperClusterTWP_sigmaIphiIphi.clear(); 
   deepSuperClusterTWP_full5x5_e5x5.clear();
   deepSuperClusterTWP_full5x5_e2x2Ratio.clear();
   deepSuperClusterTWP_full5x5_e3x3Ratio.clear();
   deepSuperClusterTWP_full5x5_eMaxRatio.clear();
   deepSuperClusterTWP_full5x5_e2ndRatio.clear();
   deepSuperClusterTWP_full5x5_eTopRatio.clear();
   deepSuperClusterTWP_full5x5_eRightRatio.clear();
   deepSuperClusterTWP_full5x5_eBottomRatio.clear();
   deepSuperClusterTWP_full5x5_eLeftRatio.clear();
   deepSuperClusterTWP_full5x5_e2x5MaxRatio.clear();
   deepSuperClusterTWP_full5x5_e2x5TopRatio.clear();
   deepSuperClusterTWP_full5x5_e2x5RightRatio.clear();
   deepSuperClusterTWP_full5x5_e2x5BottomRatio.clear();
   deepSuperClusterTWP_full5x5_e2x5LeftRatio.clear();
   deepSuperClusterTWP_full5x5_swissCross.clear();
   deepSuperClusterTWP_full5x5_r9.clear();
   deepSuperClusterTWP_full5x5_sigmaIetaIeta.clear(); 
   deepSuperClusterTWP_full5x5_sigmaIetaIphi.clear(); 
   deepSuperClusterTWP_full5x5_sigmaIphiIphi.clear();    
   deepSuperClusterTWP_dR_genScore_MatchedIndex.clear();  
   deepSuperClusterTWP_dR_simScore_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_old_MatchedIndex.clear();  
   deepSuperClusterTWP_simScore_MatchedIndex.clear();  
   deepSuperClusterTWP_n_shared_xtals_MatchedIndex.clear();   
   deepSuperClusterTWP_sim_fraction_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_rechit_diff_MatchedIndex.clear();  
   deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex.clear();  
   deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex.clear();  
   deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex.clear();  
   deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex.clear();   
   deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex.clear();  
   deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex.clear();   
   deepSuperClusterTWP_dR_genScore.clear();
   deepSuperClusterTWP_dR_simScore.clear();
   deepSuperClusterTWP_sim_fraction_old.clear();
   deepSuperClusterTWP_simScore.clear();
   deepSuperClusterTWP_n_shared_xtals.clear();
   deepSuperClusterTWP_sim_fraction.clear();
   deepSuperClusterTWP_sim_fraction_1MeVCut.clear();
   deepSuperClusterTWP_sim_fraction_5MeVCut.clear();
   deepSuperClusterTWP_sim_fraction_10MeVCut.clear();
   deepSuperClusterTWP_sim_fraction_50MeVCut.clear();
   deepSuperClusterTWP_sim_fraction_100MeVCut.clear();
   deepSuperClusterTWP_sim_fraction_500MeVCut.clear();
   deepSuperClusterTWP_sim_fraction_1GeVCut.clear();
   deepSuperClusterTWP_sim_rechit_diff.clear();
   deepSuperClusterTWP_sim_rechit_fraction.clear();
   deepSuperClusterTWP_global_sim_rechit_fraction.clear();  
   deepSuperClusterTWP_hgcal_caloToCluster.clear();  
   deepSuperClusterTWP_hgcal_clusterToCalo.clear(); 
   deepSuperClusterTWP_sim_rechit_combined_fraction.clear();  
   deepSuperClusterTWP_rechit_sim_combined_fraction.clear();   
   deepSuperClusterTWP_psCluster_energy.clear();
   deepSuperClusterTWP_psCluster_eta.clear();
   deepSuperClusterTWP_psCluster_phi.clear(); 
   deepSuperClusterTWP_HoEraw.clear(); 
   deepSuperClusterTWP_HoErawBC.clear(); 

   int nSuperClusters = (superClusterEB.product())->size() + (superClusterEE.product())->size();
   superCluster_seedIndex.resize(nSuperClusters); 
   superCluster_dR_genScore.resize(nSuperClusters);
   superCluster_dR_simScore.resize(nSuperClusters);
   superCluster_sim_fraction_old.resize(nSuperClusters);
   superCluster_simScore.resize(nSuperClusters);
   superCluster_n_shared_xtals.resize(nSuperClusters);
   superCluster_sim_fraction.resize(nSuperClusters);
   superCluster_sim_fraction_1MeVCut.resize(nSuperClusters);
   superCluster_sim_fraction_5MeVCut.resize(nSuperClusters);
   superCluster_sim_fraction_10MeVCut.resize(nSuperClusters);
   superCluster_sim_fraction_50MeVCut.resize(nSuperClusters);
   superCluster_sim_fraction_100MeVCut.resize(nSuperClusters);
   superCluster_sim_fraction_500MeVCut.resize(nSuperClusters);
   superCluster_sim_fraction_1GeVCut.resize(nSuperClusters);
   superCluster_sim_rechit_diff.resize(nSuperClusters);
   superCluster_sim_rechit_fraction.resize(nSuperClusters);
   superCluster_global_sim_rechit_fraction.resize(nSuperClusters);  
   superCluster_hgcal_caloToCluster.resize(nSuperClusters);  
   superCluster_hgcal_clusterToCalo.resize(nSuperClusters);   
   superCluster_sim_rechit_combined_fraction.resize(nSuperClusters);  
   superCluster_rechit_sim_combined_fraction.resize(nSuperClusters);    
   superCluster_pfClustersIndex.resize(nSuperClusters);
   superCluster_psCluster_energy.resize((int)(superClusterEE.product())->size());
   superCluster_psCluster_eta.resize((int)(superClusterEE.product())->size());
   superCluster_psCluster_phi.resize((int)(superClusterEE.product())->size());

   int nRetunedSuperClusters = (retunedSuperClusterEB.product())->size() + (retunedSuperClusterEE.product())->size();
   retunedSuperCluster_seedIndex.resize(nRetunedSuperClusters); 
   retunedSuperCluster_dR_genScore.resize(nRetunedSuperClusters);
   retunedSuperCluster_dR_simScore.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_old.resize(nRetunedSuperClusters);
   retunedSuperCluster_simScore.resize(nRetunedSuperClusters);
   retunedSuperCluster_n_shared_xtals.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_1MeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_5MeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_10MeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_50MeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_100MeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_500MeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_fraction_1GeVCut.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_rechit_diff.resize(nRetunedSuperClusters);
   retunedSuperCluster_sim_rechit_fraction.resize(nRetunedSuperClusters);
   retunedSuperCluster_global_sim_rechit_fraction.resize(nRetunedSuperClusters);  
   retunedSuperCluster_hgcal_caloToCluster.resize(nRetunedSuperClusters);  
   retunedSuperCluster_hgcal_clusterToCalo.resize(nRetunedSuperClusters);   
   retunedSuperCluster_sim_rechit_combined_fraction.resize(nRetunedSuperClusters);  
   retunedSuperCluster_rechit_sim_combined_fraction.resize(nRetunedSuperClusters);    
   retunedSuperCluster_pfClustersIndex.resize(nRetunedSuperClusters);
   retunedSuperCluster_psCluster_energy.resize((int)(retunedSuperClusterEE.product())->size());
   retunedSuperCluster_psCluster_eta.resize((int)(retunedSuperClusterEE.product())->size());
   retunedSuperCluster_psCluster_phi.resize((int)(retunedSuperClusterEE.product())->size());

   int nDeepSuperClusters = (deepSuperClusterEB.product())->size() + (deepSuperClusterEE.product())->size();
   deepSuperCluster_seedIndex.resize(nDeepSuperClusters); 
   deepSuperCluster_dR_genScore.resize(nDeepSuperClusters);
   deepSuperCluster_dR_simScore.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_old.resize(nDeepSuperClusters);
   deepSuperCluster_simScore.resize(nDeepSuperClusters);
   deepSuperCluster_n_shared_xtals.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_1MeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_5MeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_10MeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_50MeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_100MeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_500MeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_fraction_1GeVCut.resize(nDeepSuperClusters);
   deepSuperCluster_sim_rechit_diff.resize(nDeepSuperClusters);
   deepSuperCluster_sim_rechit_fraction.resize(nDeepSuperClusters);
   deepSuperCluster_global_sim_rechit_fraction.resize(nDeepSuperClusters);  
   deepSuperCluster_hgcal_caloToCluster.resize(nDeepSuperClusters);  
   deepSuperCluster_hgcal_clusterToCalo.resize(nDeepSuperClusters);
   deepSuperCluster_sim_rechit_combined_fraction.resize(nDeepSuperClusters);  
   deepSuperCluster_rechit_sim_combined_fraction.resize(nDeepSuperClusters);
   deepSuperCluster_pfClustersIndex.resize(nDeepSuperClusters);
   deepSuperCluster_psCluster_energy.resize((int)(deepSuperClusterEE.product())->size());
   deepSuperCluster_psCluster_eta.resize((int)(deepSuperClusterEE.product())->size());
   deepSuperCluster_psCluster_phi.resize((int)(deepSuperClusterEE.product())->size());

   int nDeepSuperClusterLWPs = (deepSuperClusterLWPEB.product())->size() + (deepSuperClusterLWPEE.product())->size();
   deepSuperClusterLWP_seedIndex.resize(nDeepSuperClusterLWPs); 
   deepSuperClusterLWP_dR_genScore.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_dR_simScore.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_old.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_simScore.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_n_shared_xtals.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_1MeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_5MeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_10MeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_50MeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_100MeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_500MeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_fraction_1GeVCut.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_rechit_diff.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_rechit_fraction.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_global_sim_rechit_fraction.resize(nDeepSuperClusterLWPs);  
   deepSuperClusterLWP_hgcal_caloToCluster.resize(nDeepSuperClusterLWPs);  
   deepSuperClusterLWP_hgcal_clusterToCalo.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_sim_rechit_combined_fraction.resize(nDeepSuperClusterLWPs);  
   deepSuperClusterLWP_rechit_sim_combined_fraction.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_pfClustersIndex.resize(nDeepSuperClusterLWPs);
   deepSuperClusterLWP_psCluster_energy.resize((int)(deepSuperClusterLWPEE.product())->size());
   deepSuperClusterLWP_psCluster_eta.resize((int)(deepSuperClusterLWPEE.product())->size());
   deepSuperClusterLWP_psCluster_phi.resize((int)(deepSuperClusterLWPEE.product())->size());

   int nDeepSuperClusterTWPs = (deepSuperClusterTWPEB.product())->size() + (deepSuperClusterTWPEE.product())->size();
   deepSuperClusterTWP_seedIndex.resize(nDeepSuperClusterTWPs); 
   deepSuperClusterTWP_dR_genScore.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_dR_simScore.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_old.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_simScore.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_n_shared_xtals.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_1MeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_5MeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_10MeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_50MeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_100MeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_500MeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_fraction_1GeVCut.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_rechit_diff.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_rechit_fraction.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_global_sim_rechit_fraction.resize(nDeepSuperClusterTWPs);  
   deepSuperClusterTWP_hgcal_caloToCluster.resize(nDeepSuperClusterTWPs);  
   deepSuperClusterTWP_hgcal_clusterToCalo.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_sim_rechit_combined_fraction.resize(nDeepSuperClusterTWPs);  
   deepSuperClusterTWP_rechit_sim_combined_fraction.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_pfClustersIndex.resize(nDeepSuperClusterTWPs);
   deepSuperClusterTWP_psCluster_energy.resize((int)(deepSuperClusterTWPEE.product())->size());
   deepSuperClusterTWP_psCluster_eta.resize((int)(deepSuperClusterTWPEE.product())->size());
   deepSuperClusterTWP_psCluster_phi.resize((int)(deepSuperClusterTWPEE.product())->size());
  
   hitsAndEnergies_CaloPart.clear();
   hitsAndEnergies_CaloPart_1MeVCut.clear();
   hitsAndEnergies_CaloPart_5MeVCut.clear();
   hitsAndEnergies_CaloPart_10MeVCut.clear();
   hitsAndEnergies_CaloPart_50MeVCut.clear();
   hitsAndEnergies_CaloPart_100MeVCut.clear();
   hitsAndEnergies_PFCluster.clear();
   hitsAndEnergies_SuperClusterEB.clear();
   hitsAndEnergies_SuperClusterEE.clear();
   hitsAndEnergies_RetunedSuperClusterEB.clear();
   hitsAndEnergies_RetunedSuperClusterEE.clear();
   hitsAndEnergies_DeepSuperClusterEB.clear();
   hitsAndEnergies_DeepSuperClusterEE.clear();
   hitsAndEnergies_DeepSuperClusterLWPEB.clear();
   hitsAndEnergies_DeepSuperClusterLWPEE.clear();
   hitsAndEnergies_DeepSuperClusterTWPEB.clear();
   hitsAndEnergies_DeepSuperClusterTWPEE.clear(); 

   GlobalPoint caloParticle_position;
   GlobalPoint cell;

   for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
   
       const auto& genParticles_caloPart = caloParts.at(iCalo).genParticles();
       caloParticle_id.push_back(caloParts.at(iCalo).pdgId());
       if(genParticles_caloPart.empty()){
          cout << "WARNING: no associated genParticle found, making standard dR matching" << endl;
          float dR=999.;
          int igen_tmp=-1; 
          int igen=0; 
          for(const auto& iGen : *(genParticles.product()))
          {
              float dR_tmp = deltaR(caloParts.at(iCalo).eta(),caloParts.at(iCalo).phi(),iGen.eta(),iGen.phi());  
              if(dR_tmp<dR && iGen.status()==1){
                 dR=dR_tmp;
                 igen_tmp=igen;
              }  
              igen++;
          } 
          const auto& genParticles_tmp = *(genParticles.product());
          auto genParticle = genParticles_tmp[igen_tmp]; 
          caloParticle_genEnergy.push_back(reduceFloat(genParticle.energy(),nBits_));
          caloParticle_genPt.push_back(reduceFloat(genParticle.pt(),nBits_));
          caloParticle_genEta.push_back(reduceFloat(genParticle.eta(),nBits_));
          caloParticle_genPhi.push_back(reduceFloat(genParticle.phi(),nBits_));
       }else{
          caloParticle_genEnergy.push_back(reduceFloat((*genParticles_caloPart.begin())->energy(),nBits_));
          caloParticle_genPt.push_back(reduceFloat((*genParticles_caloPart.begin())->pt(),nBits_));
          caloParticle_genEta.push_back(reduceFloat((*genParticles_caloPart.begin())->eta(),nBits_));
          caloParticle_genPhi.push_back(reduceFloat((*genParticles_caloPart.begin())->phi(),nBits_));
       }
      
       hitsAndEnergies_CaloPart.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),-1.));
       hitsAndEnergies_CaloPart_1MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.001)); 
       hitsAndEnergies_CaloPart_5MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.005));    
       hitsAndEnergies_CaloPart_10MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.01)); 
       hitsAndEnergies_CaloPart_50MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.05)); 
       hitsAndEnergies_CaloPart_100MeVCut.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)),0.1));     
       GlobalPoint caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
       caloParticle_simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
       caloParticle_simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
       caloParticle_simPt.push_back(reduceFloat(sqrt(caloParts.at(iCalo).px()*caloParts.at(iCalo).px() + caloParts.at(iCalo).py()*caloParts.at(iCalo).py() + caloParts.at(iCalo).pz()*caloParts.at(iCalo).pz())/TMath::CosH(caloParticle_position.eta()),nBits_));   
       if(std::abs(caloParticle_position.eta()) < 1.479){  
          EBDetId eb_id(_ebGeom->getClosestCell(caloParticle_position));  
          caloParticle_simIeta.push_back(eb_id.ieta());
          caloParticle_simIphi.push_back(eb_id.iphi());
          caloParticle_simIz.push_back(0); 
       }else{            
          int iz=-99;
          EEDetId ee_id(_eeGeom->getClosestCell(caloParticle_position));   
          caloParticle_simIeta.push_back(ee_id.ix());
          caloParticle_simIphi.push_back(ee_id.iy());
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1;  
          caloParticle_simIz.push_back(iz); 
       }   
   }
   
   //save hitsAndEnergies for each CaloParticle, PFcluster and SuperCluster
   if(saveCaloParticles_){
      hits_CaloPart.resize(nCaloParticles);
      for(unsigned int iCalo = 0; iCalo < hitsAndEnergies_CaloPart.size(); iCalo++){
          for(unsigned int i = 0; i < hitsAndEnergies_CaloPart.at(iCalo).size(); i++)
              hits_CaloPart[iCalo].push_back(hitsAndEnergies_CaloPart.at(iCalo).at(i).first);
      }
   } 

   if(savePFCluster_){
      hits_PFCluster.resize(nPFClusters);
      for(const auto& iPFCluster : *(pfClusters.product())){  
          reco::CaloCluster caloBC(iPFCluster);
          hitsAndEnergies_PFCluster.push_back(*getHitsAndEnergiesBC(&caloBC,recHitsEB,recHitsEE));
      }
      for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){
          for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++)
              hits_PFCluster[iPFCl].push_back(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);
      }
   }
        
   if(saveSuperCluster_){
     for(const auto& iSuperCluster : *(superClusterEB.product())) 
         hitsAndEnergies_SuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE));
     for(const auto& iSuperCluster : *(superClusterEE.product())) 
         hitsAndEnergies_SuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE));
   }

   if(saveSuperCluster_ && useRetunedSC_){
     for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEB.product())) 
         hitsAndEnergies_RetunedSuperClusterEB.push_back(*getHitsAndEnergiesSC(&iRetunedSuperCluster,recHitsEB,recHitsEE));
     for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEE.product())) 
         hitsAndEnergies_RetunedSuperClusterEE.push_back(*getHitsAndEnergiesSC(&iRetunedSuperCluster,recHitsEB,recHitsEE));
   }

   if(saveSuperCluster_ && useDeepSC_){
     for(const auto& iSuperCluster : *(deepSuperClusterEB.product())) 
         hitsAndEnergies_DeepSuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE));
     for(const auto& iSuperCluster : *(deepSuperClusterEE.product())) 
         hitsAndEnergies_DeepSuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE));
     for(const auto& iSuperClusterLWP : *(deepSuperClusterLWPEB.product())) 
         hitsAndEnergies_DeepSuperClusterLWPEB.push_back(*getHitsAndEnergiesSC(&iSuperClusterLWP,recHitsEB,recHitsEE));
     for(const auto& iSuperClusterLWP : *(deepSuperClusterLWPEE.product())) 
         hitsAndEnergies_DeepSuperClusterLWPEE.push_back(*getHitsAndEnergiesSC(&iSuperClusterLWP,recHitsEB,recHitsEE));
     for(const auto& iSuperClusterTWP : *(deepSuperClusterTWPEB.product())) 
         hitsAndEnergies_DeepSuperClusterTWPEB.push_back(*getHitsAndEnergiesSC(&iSuperClusterTWP,recHitsEB,recHitsEE));
     for(const auto& iSuperClusterTWP : *(deepSuperClusterTWPEE.product())) 
         hitsAndEnergies_DeepSuperClusterTWPEE.push_back(*getHitsAndEnergiesSC(&iSuperClusterTWP,recHitsEB,recHitsEE));
   }   

   //save simhits information
   for(unsigned int iCaloCount=0; iCaloCount<hitsAndEnergies_CaloPart.size(); iCaloCount++) 
   {
       float calo_simEnergy=0.;

       for(auto const& hit: hitsAndEnergies_CaloPart[iCaloCount])
       {
           DetId id(hit.first);
           if(id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap) continue;
               
           calo_simEnergy += hit.second; 

           cell = geometry->getPosition(id);
           float eta = cell.eta();  
           float phi = cell.phi();  
           int ieta = -99; 
           int iphi = -99;
           int iz = -99;  
           if(id.subdetId()==EcalBarrel){
              EBDetId eb_id(id);
              ieta = eb_id.ieta(); 
              iphi = eb_id.iphi();  
              iz = 0;   
           }
           if(id.subdetId()==EcalEndcap){
              EEDetId ee_id(id);
              ieta = ee_id.ix(); 
              iphi = ee_id.iy();
              if(ee_id.zside()<0) iz=-1;
              if(ee_id.zside()>0) iz=1;   
           }

           if(saveSimhits_){
              simHit_energy[iCaloCount].push_back(reduceFloat(hit.second,nBits_));
              simHit_eta[iCaloCount].push_back(reduceFloat(eta,nBits_));
              simHit_phi[iCaloCount].push_back(reduceFloat(phi,nBits_));
              simHit_ieta[iCaloCount].push_back(ieta);
              simHit_iphi[iCaloCount].push_back(iphi);
              simHit_iz[iCaloCount].push_back(iz); 
           }
       } 
       caloParticle_simEnergy.push_back(reduceFloat(calo_simEnergy,nBits_));
   }
   
   //Save PFClusters 
   if(savePFCluster_){
      
      int iPFCl=0;
      //std::cout << "PFClusters size     : " << (pfClusters.product())->size() << std::endl;
      for(const auto& iPFCluster : *(pfClusters.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();    
          sim_rechit_combined_fraction.clear();
          rechit_sim_combined_fraction.clear();

          pfCluster_energy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_eta.push_back(reduceFloat(iPFCluster.eta(),nBits_));
          pfCluster_phi.push_back(reduceFloat(iPFCluster.phi(),nBits_));
          reco::CaloCluster caloBC(iPFCluster);

          math::XYZPoint caloPos = caloBC.position();
          if(iPFCluster.layer() == PFLayer::ECAL_BARREL){  
             EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             pfCluster_ieta.push_back(eb_id.ieta());
             pfCluster_iphi.push_back(eb_id.iphi());
             pfCluster_iz.push_back(0); 
          }else if(iPFCluster.layer() == PFLayer::ECAL_ENDCAP){ 
             int iz=-99;
             EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1;     
             pfCluster_ieta.push_back(ee_id.ix());
             pfCluster_iphi.push_back(ee_id.iy());
             pfCluster_iz.push_back(iz); 
          } 
           
          if(saveShowerShapes_ && iPFCluster.layer() == PFLayer::ECAL_BARREL){
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), &(*topology));  
             pfCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             pfCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             pfCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     pfCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             pfCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             pfCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             pfCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             pfCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             pfCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             pfCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             pfCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             pfCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             pfCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             pfCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             pfCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             pfCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             pfCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             pfCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             pfCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             pfCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             pfCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             pfCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             pfCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             pfCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             pfCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             pfCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             pfCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             pfCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             pfCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             pfCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             pfCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             pfCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             pfCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             pfCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             pfCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             pfCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             pfCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             pfCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }else if(saveShowerShapes_ && iPFCluster.layer() == PFLayer::ECAL_ENDCAP){ 
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), &(*topology));  
             pfCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             pfCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             pfCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     pfCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             pfCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             pfCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             pfCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             pfCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             pfCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             pfCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             pfCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             pfCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             pfCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             pfCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             pfCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             pfCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             pfCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             pfCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             pfCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             pfCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             pfCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             pfCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             pfCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             pfCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             pfCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             pfCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             pfCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             pfCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             pfCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             pfCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             pfCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             pfCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             pfCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             pfCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             pfCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             pfCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             pfCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             pfCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 
          }  
          
          if(savePFClusterhits_){
             //for save PFClusterHit      
             for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){      
                 cell = geometry->getPosition(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);
                 pfClusterHit_energy[iPFCl].push_back(reduceFloat(hitsAndEnergies_PFCluster.at(iPFCl).at(i).second,nBits_));
                 pfClusterHit_eta[iPFCl].push_back(reduceFloat(cell.eta(),nBits_));
                 pfClusterHit_phi[iPFCl].push_back(reduceFloat(cell.phi(),nBits_));
                 if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalBarrel){ 
                    EBDetId eb_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first); 
                    pfClusterHit_rechitEnergy[iPFCl].push_back(reduceFloat((*(recHitsEB.product())->find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first)).energy(),nBits_)); 
                    pfClusterHit_ieta[iPFCl].push_back(eb_id.ieta());
                    pfClusterHit_iphi[iPFCl].push_back(eb_id.iphi());
                    pfClusterHit_iz[iPFCl].push_back(0); 
                 }else if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalEndcap){  
                    int iz=-99;
                    EEDetId ee_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);  
                    pfClusterHit_rechitEnergy[iPFCl].push_back(reduceFloat((*(recHitsEE.product())->find(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first)).energy(),nBits_)); 
                    pfClusterHit_ieta[iPFCl].push_back(ee_id.ix());
                    pfClusterHit_iphi[iPFCl].push_back(ee_id.iy());
                    if(ee_id.zside()<0) iz=-1;
                    if(ee_id.zside()>0) iz=1;   
                    pfClusterHit_iz[iPFCl].push_back(iz); 
                 } 
             }
          }
   
          //compute scores     
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             pfCluster_dR_genScore[iPFCl] = dR_genScore;        
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) pfCluster_dR_genScore_MatchedIndex.push_back(-1);
             else pfCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_PFCluster.at(iPFCl),&hitsAndEnergies_CaloPart.at(iCalo),recHitsEB,recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE);    
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_PFCluster.at(iPFCl), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);              
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_simScore.push_back(999.);  
 
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);  
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);  
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 if(scoreType_=="sim_rechit_combined_fraction") simScore.push_back(scores[17]);  
                 if(scoreType_=="rechit_sim_combined_fraction") simScore.push_back(scores[18]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);
                 sim_rechit_combined_fraction.push_back(scores[17]);  
                 rechit_sim_combined_fraction.push_back(scores[18]);        
                     
             } 
             pfCluster_nXtals.push_back((iPFCluster.hitsAndFractions()).size());   
             pfCluster_dR_simScore[iPFCl] = dR_simScore;  
             pfCluster_sim_fraction_old[iPFCl] = sim_fraction_old;        
             pfCluster_simScore[iPFCl] = simScore; 
             pfCluster_n_shared_xtals[iPFCl] = n_shared_xtals;
             pfCluster_sim_fraction[iPFCl] = sim_fraction; 
             pfCluster_sim_fraction_1MeVCut[iPFCl] = sim_fraction_1MeVCut; 
             pfCluster_sim_fraction_5MeVCut[iPFCl] = sim_fraction_5MeVCut;  
             pfCluster_sim_fraction_10MeVCut[iPFCl] = sim_fraction_10MeVCut;  
             pfCluster_sim_fraction_50MeVCut[iPFCl] = sim_fraction_50MeVCut;  
             pfCluster_sim_fraction_100MeVCut[iPFCl] = sim_fraction_100MeVCut;       
             pfCluster_sim_fraction_500MeVCut[iPFCl] = sim_fraction_500MeVCut;       
             pfCluster_sim_fraction_1GeVCut[iPFCl] = sim_fraction_1GeVCut;            
             pfCluster_sim_rechit_diff[iPFCl] = sim_rechit_diff; 
             pfCluster_sim_rechit_fraction[iPFCl] = sim_rechit_fraction;           
             pfCluster_global_sim_rechit_fraction[iPFCl] = global_sim_rechit_fraction;
             pfCluster_hgcal_caloToCluster[iPFCl] = hgcal_caloToCluster; 
             pfCluster_hgcal_clusterToCalo[iPFCl] = hgcal_clusterToCalo;  
             pfCluster_sim_rechit_combined_fraction[iPFCl] = sim_rechit_combined_fraction; 
             pfCluster_rechit_sim_combined_fraction[iPFCl] = rechit_sim_combined_fraction;     
           
             pfCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));
             pfCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                   else pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                }else{  
                   if(iPFCluster.layer() == PFLayer::ECAL_BARREL) pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_old, pfCluster_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iPFCl));
                   else if(iPFCluster.layer() == PFLayer::ECAL_ENDCAP) pfCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({pfCluster_sim_fraction_old, pfCluster_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iPFCl));                
                } 
             }else{
                pfCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));     
                pfCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));    
                pfCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));   
                pfCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));   
                pfCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));   
                pfCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));    
                pfCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));      
                pfCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));          
                pfCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));          
                pfCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl)); 
                pfCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&pfCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));      
                pfCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&pfCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));  
                pfCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPFCl));          
             }   
          }    
    
          iPFCl++;        
      } 
   }
  
   //save inverse of matchings
   if(saveCaloParticles_ && savePFCluster_){ 
      fillParticleMatchedIndex(&genParticle_pfCluster_dR_genScore_MatchedIndex,&pfCluster_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && savePFCluster_){ 
      fillParticleMatchedIndex(&caloParticle_pfCluster_dR_simScore_MatchedIndex,&pfCluster_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_old_MatchedIndex,&pfCluster_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_pfCluster_simScore_MatchedIndex,&pfCluster_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_pfCluster_n_shared_xtals_MatchedIndex,&pfCluster_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_MatchedIndex,&pfCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_1MeVCut_MatchedIndex,&pfCluster_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_5MeVCut_MatchedIndex,&pfCluster_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_10MeVCut_MatchedIndex,&pfCluster_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_50MeVCut_MatchedIndex,&pfCluster_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_100MeVCut_MatchedIndex,&pfCluster_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_500MeVCut_MatchedIndex,&pfCluster_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_1GeVCut_MatchedIndex,&pfCluster_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_rechit_diff_MatchedIndex,&pfCluster_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex,&pfCluster_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex,&pfCluster_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_pfCluster_hgcal_caloToCluster_MatchedIndex,&pfCluster_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_hgcal_clusterToCalo_MatchedIndex,&pfCluster_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_rechit_combined_fraction_MatchedIndex,&pfCluster_sim_rechit_combined_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_pfCluster_rechit_sim_combined_fraction_MatchedIndex,&pfCluster_rechit_sim_combined_fraction_MatchedIndex);            
      }      
   } 
   
   //Save SuperClusters 
   locCov_.clear();
   full5x5_locCov_.clear();
   if(saveSuperCluster_){
      int iSC=0;
      //std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
      for(const auto& iSuperCluster : *(superClusterEB.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();   

          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          superCluster_etaWidth.push_back(reduceFloat(iSuperCluster.etaWidth(),nBits_));
          superCluster_phiWidth.push_back(reduceFloat(iSuperCluster.phiWidth(),nBits_));
          superCluster_R.push_back(reduceFloat(iSuperCluster.position().R(),nBits_));
          superCluster_nPFClusters.push_back(iSuperCluster.clusters().size());
          math::XYZPoint caloPos = iSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          superCluster_ieta.push_back(eb_id.ieta());
          superCluster_iphi.push_back(eb_id.iphi());
          superCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             superCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             superCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             superCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     superCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             superCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             superCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             superCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             superCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             superCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             superCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             superCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             superCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             superCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             superCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             superCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             superCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             superCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             superCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             superCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             superCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             superCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             superCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             superCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             superCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             superCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             superCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             superCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             superCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             superCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             superCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             superCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             superCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             superCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             superCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             superCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             superCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) superCluster_dR_genScore_MatchedIndex.push_back(-1);
             else superCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);        
                   
             } 

             superCluster_dR_simScore[iSC] = dR_simScore;  
             superCluster_sim_fraction_old[iSC] = sim_fraction_old;     
             superCluster_simScore[iSC] = simScore;  
             superCluster_n_shared_xtals[iSC] = n_shared_xtals;   
             superCluster_sim_fraction[iSC] = sim_fraction;  
             superCluster_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             superCluster_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             superCluster_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             superCluster_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             superCluster_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             superCluster_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             superCluster_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             superCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
             superCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             superCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             superCluster_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             superCluster_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 
            
             superCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             superCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") superCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else superCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   superCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({superCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({superCluster_sim_fraction_old, superCluster_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iSC));    
                } 
             }else{
                superCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                superCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                superCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                superCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                superCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                superCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                superCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                superCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                superCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                superCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                superCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&superCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                superCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&superCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));        
             }
          }
          
          if(savePFCluster_){   
             //save clusters and superClusters mutual info
             reco::CaloCluster caloSeed(*iSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) superCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) superCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global SuperCluster indexing for EE has an offset = nSuperClusterEB
      iSC = (superClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;

      for(const auto& iSuperCluster : *(superClusterEE.product())){    

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();   
          iSC_tmp++;
        
          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          superCluster_etaWidth.push_back(reduceFloat(iSuperCluster.etaWidth(),nBits_));
          superCluster_phiWidth.push_back(reduceFloat(iSuperCluster.phiWidth(),nBits_));
          superCluster_R.push_back(reduceFloat(iSuperCluster.position().R(),nBits_)); 
          superCluster_nPFClusters.push_back(iSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          superCluster_ieta.push_back(ee_id.ix());
          superCluster_iphi.push_back(ee_id.iy());
          superCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             superCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             superCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             superCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     superCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             superCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             superCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             superCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             superCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             superCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             superCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             superCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             superCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             superCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             superCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             superCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             superCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             superCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             superCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             superCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             superCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             superCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             superCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             superCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             superCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             superCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             superCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             superCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             superCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             superCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             superCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             superCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             superCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             superCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             superCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             superCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             superCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) superCluster_dR_genScore_MatchedIndex.push_back(-1);
             else superCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 

          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);
                 
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
   
             } 

             superCluster_dR_simScore[iSC] = dR_simScore;  
             superCluster_sim_fraction_old[iSC] = sim_fraction_old;     
             superCluster_simScore[iSC] = simScore;  
             superCluster_n_shared_xtals[iSC] = n_shared_xtals;   
             superCluster_sim_fraction[iSC] = sim_fraction;  
             superCluster_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             superCluster_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             superCluster_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             superCluster_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             superCluster_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             superCluster_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             superCluster_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             superCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
             superCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             superCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             superCluster_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             superCluster_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 
             
             superCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             superCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") superCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else superCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   superCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({superCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({superCluster_sim_fraction_old, superCluster_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iSC));                     
                } 
             }else{
                superCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                superCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                superCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                superCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                superCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                superCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                superCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                superCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                superCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                superCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                superCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&superCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                superCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&superCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                superCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));         
             } 
          }
          
          if(iSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  superCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  superCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  superCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_){   
             //save clusters and superClusters mutual info
             reco::CaloCluster caloSeed(*iSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) superCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) superCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }
   }

   //save pfCluster_superClustersIndex
   if(savePFCluster_ && saveSuperCluster_){
      for(unsigned int iSC=0; iSC<superCluster_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex.at(iSC).size(); iPF++)
              if(superCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_superClustersIndex[superCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
    }

   //save inverse of matchings
   if(saveCaloParticles_ && saveSuperCluster_){ 
      fillParticleMatchedIndex(&genParticle_superCluster_dR_genScore_MatchedIndex,&superCluster_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && saveSuperCluster_){ 
      fillParticleMatchedIndex(&caloParticle_superCluster_dR_simScore_MatchedIndex,&superCluster_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_old_MatchedIndex,&superCluster_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_superCluster_simScore_MatchedIndex,&superCluster_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_superCluster_n_shared_xtals_MatchedIndex,&superCluster_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_MatchedIndex,&superCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_1MeVCut_MatchedIndex,&superCluster_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_5MeVCut_MatchedIndex,&superCluster_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_10MeVCut_MatchedIndex,&superCluster_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_50MeVCut_MatchedIndex,&superCluster_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_100MeVCut_MatchedIndex,&superCluster_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_500MeVCut_MatchedIndex,&superCluster_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_1GeVCut_MatchedIndex,&superCluster_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_rechit_diff_MatchedIndex,&superCluster_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_rechit_fraction_MatchedIndex,&superCluster_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex,&superCluster_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_superCluster_hgcal_caloToCluster_MatchedIndex,&superCluster_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_hgcal_clusterToCalo_MatchedIndex,&superCluster_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_rechit_combined_fraction_MatchedIndex,&superCluster_sim_rechit_combined_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_rechit_sim_combined_fraction_MatchedIndex,&superCluster_rechit_sim_combined_fraction_MatchedIndex);       
      }      
   }

   //Save retunedSuperClusters 
   //std::cout << "-----> retunedSuperClusters <-----" << std::endl;  
   locCov_.clear();
   full5x5_locCov_.clear();
   if(saveSuperCluster_ && useRetunedSC_){
      int iSC=0;
      //std::cout << "retunedSuperClustersEB size: " << (retunedSuperClusterEB.product())->size() << std::endl;
      for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEB.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();   

          retunedSuperCluster_energy.push_back(reduceFloat(iRetunedSuperCluster.energy(),nBits_));
          retunedSuperCluster_eta.push_back(reduceFloat(iRetunedSuperCluster.eta(),nBits_));
          retunedSuperCluster_phi.push_back(reduceFloat(iRetunedSuperCluster.phi(),nBits_));
          retunedSuperCluster_etaWidth.push_back(reduceFloat(iRetunedSuperCluster.etaWidth(),nBits_));
          retunedSuperCluster_phiWidth.push_back(reduceFloat(iRetunedSuperCluster.phiWidth(),nBits_));
          retunedSuperCluster_R.push_back(reduceFloat(iRetunedSuperCluster.position().R(),nBits_));
          retunedSuperCluster_nPFClusters.push_back(iRetunedSuperCluster.clusters().size());
          math::XYZPoint caloPos = iRetunedSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          retunedSuperCluster_ieta.push_back(eb_id.ieta());
          retunedSuperCluster_iphi.push_back(eb_id.iphi());
          retunedSuperCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iRetunedSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             retunedSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             retunedSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             retunedSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     retunedSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             retunedSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             retunedSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             retunedSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             retunedSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             retunedSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             retunedSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             retunedSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             retunedSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             retunedSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             retunedSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             retunedSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             retunedSuperCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             retunedSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             retunedSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             retunedSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             retunedSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             retunedSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             retunedSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             retunedSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             retunedSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             retunedSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             retunedSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             retunedSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             retunedSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             retunedSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             retunedSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             retunedSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             retunedSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             retunedSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             retunedSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             retunedSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             retunedSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iRetunedSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             retunedSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             retunedSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             retunedSuperCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) retunedSuperCluster_dR_genScore_MatchedIndex.push_back(-1);
             else retunedSuperCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_RetunedSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 

             retunedSuperCluster_dR_simScore[iSC] = dR_simScore;  
             retunedSuperCluster_sim_fraction_old[iSC] = sim_fraction_old;   
             retunedSuperCluster_simScore[iSC] = simScore;  
             retunedSuperCluster_n_shared_xtals[iSC] = n_shared_xtals;  
             retunedSuperCluster_sim_fraction[iSC] = sim_fraction;  
             retunedSuperCluster_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             retunedSuperCluster_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             retunedSuperCluster_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             retunedSuperCluster_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             retunedSuperCluster_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             retunedSuperCluster_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             retunedSuperCluster_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             retunedSuperCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
             retunedSuperCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             retunedSuperCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             retunedSuperCluster_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             retunedSuperCluster_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 
            
             retunedSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             retunedSuperCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") retunedSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else retunedSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   retunedSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({retunedSuperCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({retunedSuperCluster_sim_fraction_old, retunedSuperCluster_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iSC));      
                } 
             }else{
                retunedSuperCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                retunedSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                retunedSuperCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                retunedSuperCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));        
             }             
          }

          if(savePFCluster_){   
             //save clusters and retunedSuperClusters mutual info
             reco::CaloCluster caloSeed(*iRetunedSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iRetunedSuperCluster.clustersBegin(); iBC != iRetunedSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) retunedSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) retunedSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global retunedSuperCluster indexing for EE has an offset = nretunedSuperClusterEB
      iSC = (retunedSuperClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "retunedSuperClustersEE size: " << (retunedSuperClusterEE.product())->size() << std::endl;
      for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEE.product())){    

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();  
          iSC_tmp++;
        
          retunedSuperCluster_energy.push_back(reduceFloat(iRetunedSuperCluster.energy(),nBits_));
          retunedSuperCluster_eta.push_back(reduceFloat(iRetunedSuperCluster.eta(),nBits_));
          retunedSuperCluster_phi.push_back(reduceFloat(iRetunedSuperCluster.phi(),nBits_));
          retunedSuperCluster_etaWidth.push_back(reduceFloat(iRetunedSuperCluster.etaWidth(),nBits_));
          retunedSuperCluster_phiWidth.push_back(reduceFloat(iRetunedSuperCluster.phiWidth(),nBits_));
          retunedSuperCluster_R.push_back(reduceFloat(iRetunedSuperCluster.position().R(),nBits_)); 
          retunedSuperCluster_nPFClusters.push_back(iRetunedSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iRetunedSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          retunedSuperCluster_ieta.push_back(ee_id.ix());
          retunedSuperCluster_iphi.push_back(ee_id.iy());
          retunedSuperCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iRetunedSuperCluster.seed()); 
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             retunedSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             retunedSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             retunedSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     retunedSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             retunedSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             retunedSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             retunedSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             retunedSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             retunedSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             retunedSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             retunedSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             retunedSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             retunedSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             retunedSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             retunedSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             retunedSuperCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             retunedSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             retunedSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             retunedSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             retunedSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             retunedSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             retunedSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             retunedSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             retunedSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             retunedSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             retunedSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             retunedSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             retunedSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             retunedSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             retunedSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             retunedSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             retunedSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             retunedSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             retunedSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             retunedSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             retunedSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iRetunedSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             retunedSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             retunedSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             retunedSuperCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) retunedSuperCluster_dR_genScore_MatchedIndex.push_back(-1);
             else retunedSuperCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_RetunedSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_RetunedSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 
             
             retunedSuperCluster_dR_simScore[iSC] = dR_simScore;  
             retunedSuperCluster_sim_fraction_old[iSC] = sim_fraction_old;   
             retunedSuperCluster_simScore[iSC] = simScore;  
             retunedSuperCluster_n_shared_xtals[iSC] = n_shared_xtals;  
             retunedSuperCluster_sim_fraction[iSC] = sim_fraction;  
             retunedSuperCluster_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             retunedSuperCluster_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             retunedSuperCluster_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             retunedSuperCluster_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             retunedSuperCluster_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             retunedSuperCluster_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             retunedSuperCluster_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             retunedSuperCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
             retunedSuperCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             retunedSuperCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             retunedSuperCluster_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             retunedSuperCluster_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 

             retunedSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             retunedSuperCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") retunedSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else retunedSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   retunedSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({retunedSuperCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({retunedSuperCluster_sim_fraction_old, retunedSuperCluster_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iSC));                     
                } 
             }else{
                retunedSuperCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                retunedSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                retunedSuperCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                retunedSuperCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));           
             }
          }

          if(iRetunedSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iRetunedSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iRetunedSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  retunedSuperCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iRetunedSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  retunedSuperCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iRetunedSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  retunedSuperCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iRetunedSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_){   
             //save clusters and retunedSuperClusters mutual info
             reco::CaloCluster caloSeed(*iRetunedSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iRetunedSuperCluster.clustersBegin(); iBC != iRetunedSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) retunedSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) retunedSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }
   }

   //save pfCluster_retunedSuperClustersIndex
   if(savePFCluster_ && saveSuperCluster_ && useDeepSC_){
      for(unsigned int iSC=0; iSC<retunedSuperCluster_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<retunedSuperCluster_pfClustersIndex.at(iSC).size(); iPF++)
              if(retunedSuperCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_retunedSuperClustersIndex[retunedSuperCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
    }
   
   //save inverse of matchings
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&genParticle_retunedSuperCluster_dR_genScore_MatchedIndex,&retunedSuperCluster_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex,&retunedSuperCluster_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_old_MatchedIndex,&retunedSuperCluster_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_simScore_MatchedIndex,&retunedSuperCluster_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_n_shared_xtals_MatchedIndex,&retunedSuperCluster_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex,&retunedSuperCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex,&retunedSuperCluster_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_rechit_diff_MatchedIndex,&retunedSuperCluster_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_rechit_fraction_MatchedIndex,&retunedSuperCluster_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex,&retunedSuperCluster_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_hgcal_caloToCluster_MatchedIndex,&retunedSuperCluster_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex,&retunedSuperCluster_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex,&retunedSuperCluster_sim_rechit_combined_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex,&retunedSuperCluster_rechit_sim_combined_fraction_MatchedIndex);
      }      
   }
   
   //Save deepSuperClusters 
   //std::cout << "-----> deepSuperClusters <-----" << std::endl;  
   locCov_.clear();
   full5x5_locCov_.clear();
   if(saveSuperCluster_ && useDeepSC_){
      int iSC=0;
      //std::cout << "deepSuperClustersEB size: " << (deepSuperClusterEB.product())->size() << std::endl;
      for(const auto& iDeepSuperCluster : *(deepSuperClusterEB.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();   

          deepSuperCluster_energy.push_back(reduceFloat(iDeepSuperCluster.energy(),nBits_));
          deepSuperCluster_eta.push_back(reduceFloat(iDeepSuperCluster.eta(),nBits_));
          deepSuperCluster_phi.push_back(reduceFloat(iDeepSuperCluster.phi(),nBits_));
          deepSuperCluster_etaWidth.push_back(reduceFloat(iDeepSuperCluster.etaWidth(),nBits_));
          deepSuperCluster_phiWidth.push_back(reduceFloat(iDeepSuperCluster.phiWidth(),nBits_));
          deepSuperCluster_R.push_back(reduceFloat(iDeepSuperCluster.position().R(),nBits_));
          deepSuperCluster_nPFClusters.push_back(iDeepSuperCluster.clusters().size());
          math::XYZPoint caloPos = iDeepSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          deepSuperCluster_ieta.push_back(eb_id.ieta());
          deepSuperCluster_iphi.push_back(eb_id.iphi());
          deepSuperCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iDeepSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             deepSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             deepSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             deepSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             deepSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             deepSuperCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) deepSuperCluster_dR_genScore_MatchedIndex.push_back(-1);
             else deepSuperCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_DeepSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEB.at(iSC), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 

             deepSuperCluster_dR_simScore[iSC] = dR_simScore;  
             deepSuperCluster_sim_fraction_old[iSC] = sim_fraction_old;   
             deepSuperCluster_simScore[iSC] = simScore;  
             deepSuperCluster_n_shared_xtals[iSC] = n_shared_xtals;  
             deepSuperCluster_sim_fraction[iSC] = sim_fraction;  
             deepSuperCluster_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             deepSuperCluster_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             deepSuperCluster_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             deepSuperCluster_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             deepSuperCluster_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             deepSuperCluster_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             deepSuperCluster_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             deepSuperCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
             deepSuperCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             deepSuperCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             deepSuperCluster_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             deepSuperCluster_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 
            
             deepSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             deepSuperCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") deepSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else deepSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   deepSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({deepSuperCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({deepSuperCluster_sim_fraction_old, deepSuperCluster_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iSC));      
                } 
             }else{
                deepSuperCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                deepSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                deepSuperCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));        
             }             
          }

          if(savePFCluster_){   
             //save clusters and deepSuperClusters mutual info
             reco::CaloCluster caloSeed(*iDeepSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperCluster.clustersBegin(); iBC != iDeepSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global deepSuperCluster indexing for EE has an offset = ndeepSuperClusterEB
      iSC = (deepSuperClusterEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "deepSuperClustersEE size: " << (deepSuperClusterEE.product())->size() << std::endl;
      for(const auto& iDeepSuperCluster : *(deepSuperClusterEE.product())){    

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();  
          iSC_tmp++;
        
          deepSuperCluster_energy.push_back(reduceFloat(iDeepSuperCluster.energy(),nBits_));
          deepSuperCluster_eta.push_back(reduceFloat(iDeepSuperCluster.eta(),nBits_));
          deepSuperCluster_phi.push_back(reduceFloat(iDeepSuperCluster.phi(),nBits_));
          deepSuperCluster_etaWidth.push_back(reduceFloat(iDeepSuperCluster.etaWidth(),nBits_));
          deepSuperCluster_phiWidth.push_back(reduceFloat(iDeepSuperCluster.phiWidth(),nBits_));
          deepSuperCluster_R.push_back(reduceFloat(iDeepSuperCluster.position().R(),nBits_)); 
          deepSuperCluster_nPFClusters.push_back(iDeepSuperCluster.clusters().size());  
          math::XYZPoint caloPos = iDeepSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          deepSuperCluster_ieta.push_back(ee_id.ix());
          deepSuperCluster_iphi.push_back(ee_id.iy());
          deepSuperCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iDeepSuperCluster.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             deepSuperCluster_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperCluster_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperCluster_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperCluster_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperCluster_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperCluster_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperCluster_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperCluster_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperCluster_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperCluster_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperCluster_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperCluster_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperCluster_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperCluster_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperCluster_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperCluster_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             deepSuperCluster_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperCluster_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperCluster_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperCluster_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperCluster_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperCluster_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperCluster_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperCluster_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperCluster_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperCluster_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperCluster_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperCluster_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperCluster_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperCluster_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperCluster_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperCluster_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperCluster_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperCluster_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             deepSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_));

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperCluster, towerIso1_, towerIso2_, egammaHadTower_);
             deepSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));  
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             deepSuperCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) deepSuperCluster_dR_genScore_MatchedIndex.push_back(-1);
             else deepSuperCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_DeepSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 
             
             deepSuperCluster_dR_simScore[iSC] = dR_simScore;  
             deepSuperCluster_sim_fraction_old[iSC] = sim_fraction_old;   
             deepSuperCluster_simScore[iSC] = simScore;  
             deepSuperCluster_n_shared_xtals[iSC] = n_shared_xtals;  
             deepSuperCluster_sim_fraction[iSC] = sim_fraction;  
             deepSuperCluster_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             deepSuperCluster_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             deepSuperCluster_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             deepSuperCluster_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             deepSuperCluster_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             deepSuperCluster_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             deepSuperCluster_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             deepSuperCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
             deepSuperCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             deepSuperCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             deepSuperCluster_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             deepSuperCluster_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 

             deepSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             deepSuperCluster_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") deepSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else deepSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   deepSuperCluster_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({deepSuperCluster_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({deepSuperCluster_sim_fraction_old, deepSuperCluster_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iSC));                     
                } 
             }else{
                deepSuperCluster_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                deepSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperCluster_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                deepSuperCluster_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperCluster_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));           
             }
          }

          if(iDeepSuperCluster.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iDeepSuperCluster.preshowerClusters().size(); iPC++){
                  if(!iDeepSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  deepSuperCluster_psCluster_energy[iSC_tmp].push_back(reduceFloat(iDeepSuperCluster.preshowerClusters()[iPC]->energy(),nBits_));
                  deepSuperCluster_psCluster_eta[iSC_tmp].push_back(reduceFloat(iDeepSuperCluster.preshowerClusters()[iPC]->eta(),nBits_));
                  deepSuperCluster_psCluster_phi[iSC_tmp].push_back(reduceFloat(iDeepSuperCluster.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_){   
             //save clusters and deepSuperClusters mutual info
             reco::CaloCluster caloSeed(*iDeepSuperCluster.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperCluster.clustersBegin(); iBC != iDeepSuperCluster.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperCluster_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperCluster_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }
   }

   //save pfCluster_deepSuperClustersIndex
   if(savePFCluster_ && saveSuperCluster_ && useDeepSC_){
      for(unsigned int iSC=0; iSC<deepSuperCluster_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<deepSuperCluster_pfClustersIndex.at(iSC).size(); iPF++)
              if(deepSuperCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_deepSuperClustersIndex[deepSuperCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
    }

   //save inverse of matchings
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&genParticle_deepSuperCluster_dR_genScore_MatchedIndex,&deepSuperCluster_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&caloParticle_deepSuperCluster_dR_simScore_MatchedIndex,&deepSuperCluster_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_old_MatchedIndex,&deepSuperCluster_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_simScore_MatchedIndex,&deepSuperCluster_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_n_shared_xtals_MatchedIndex,&deepSuperCluster_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_MatchedIndex,&deepSuperCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex,&deepSuperCluster_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_rechit_diff_MatchedIndex,&deepSuperCluster_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_rechit_fraction_MatchedIndex,&deepSuperCluster_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_global_sim_rechit_fraction_MatchedIndex,&deepSuperCluster_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_hgcal_caloToCluster_MatchedIndex,&deepSuperCluster_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_hgcal_clusterToCalo_MatchedIndex,&deepSuperCluster_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex,&deepSuperCluster_sim_rechit_combined_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex,&deepSuperCluster_rechit_sim_combined_fraction_MatchedIndex);
      }      
   }

   //Save deepSuperClusterLWPs 
   //std::cout << "-----> deepSuperClusterLWPs <-----" << std::endl;
   locCov_.clear();
   full5x5_locCov_.clear();
   if(saveSuperCluster_ && useDeepSC_){
      int iSC=0;
      //std::cout << "deepSuperClusterLWPsEB size: " << (deepSuperClusterLWPEB.product())->size() << std::endl;
      for(const auto& iDeepSuperClusterLWP : *(deepSuperClusterLWPEB.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();   

          deepSuperClusterLWP_energy.push_back(reduceFloat(iDeepSuperClusterLWP.energy(),nBits_));
          deepSuperClusterLWP_eta.push_back(reduceFloat(iDeepSuperClusterLWP.eta(),nBits_));
          deepSuperClusterLWP_phi.push_back(reduceFloat(iDeepSuperClusterLWP.phi(),nBits_));
          deepSuperClusterLWP_etaWidth.push_back(reduceFloat(iDeepSuperClusterLWP.etaWidth(),nBits_));
          deepSuperClusterLWP_phiWidth.push_back(reduceFloat(iDeepSuperClusterLWP.phiWidth(),nBits_));
          deepSuperClusterLWP_R.push_back(reduceFloat(iDeepSuperClusterLWP.position().R(),nBits_));
          deepSuperClusterLWP_nPFClusters.push_back(iDeepSuperClusterLWP.clusters().size());
          math::XYZPoint caloPos = iDeepSuperClusterLWP.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          deepSuperClusterLWP_ieta.push_back(eb_id.ieta());
          deepSuperClusterLWP_iphi.push_back(eb_id.iphi());
          deepSuperClusterLWP_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iDeepSuperClusterLWP.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             deepSuperClusterLWP_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperClusterLWP_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperClusterLWP_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperClusterLWP_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperClusterLWP_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperClusterLWP_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperClusterLWP_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperClusterLWP_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperClusterLWP_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperClusterLWP_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperClusterLWP_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperClusterLWP_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperClusterLWP_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperClusterLWP_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperClusterLWP_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperClusterLWP_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             deepSuperClusterLWP_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperClusterLWP_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperClusterLWP_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperClusterLWP_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperClusterLWP_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperClusterLWP_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperClusterLWP_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperClusterLWP_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperClusterLWP_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperClusterLWP_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperClusterLWP_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperClusterLWP_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperClusterLWP_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperClusterLWP_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperClusterLWP_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperClusterLWP_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperClusterLWP_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperClusterLWP_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperClusterLWP_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             deepSuperClusterLWP_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperClusterLWP_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperClusterLWP_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperClusterLWP, towerIso1_, towerIso2_, egammaHadTower_);
             deepSuperClusterLWP_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperClusterLWP_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));  
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             deepSuperClusterLWP_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) deepSuperClusterLWP_dR_genScore_MatchedIndex.push_back(-1);
             else deepSuperClusterLWP_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_DeepSuperClusterLWPEB.at(iSC), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEB.at(iSC), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEB.at(iSC), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEB.at(iSC), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEB.at(iSC), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEB.at(iSC), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 

             deepSuperClusterLWP_dR_simScore[iSC] = dR_simScore;  
             deepSuperClusterLWP_sim_fraction_old[iSC] = sim_fraction_old;   
             deepSuperClusterLWP_simScore[iSC] = simScore;  
             deepSuperClusterLWP_n_shared_xtals[iSC] = n_shared_xtals;  
             deepSuperClusterLWP_sim_fraction[iSC] = sim_fraction;  
             deepSuperClusterLWP_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             deepSuperClusterLWP_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             deepSuperClusterLWP_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             deepSuperClusterLWP_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             deepSuperClusterLWP_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             deepSuperClusterLWP_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             deepSuperClusterLWP_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             deepSuperClusterLWP_sim_rechit_diff[iSC] = sim_rechit_diff; 
             deepSuperClusterLWP_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             deepSuperClusterLWP_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             deepSuperClusterLWP_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             deepSuperClusterLWP_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 
            
             deepSuperClusterLWP_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             deepSuperClusterLWP_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") deepSuperClusterLWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else deepSuperClusterLWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   deepSuperClusterLWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({deepSuperClusterLWP_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({deepSuperClusterLWP_sim_fraction_old, deepSuperClusterLWP_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iSC));      
                } 
             }else{
                deepSuperClusterLWP_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                deepSuperClusterLWP_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterLWP_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));        
             }             
          }

          if(savePFCluster_){   
             //save clusters and deepSuperClusterLWPs mutual info
             reco::CaloCluster caloSeed(*iDeepSuperClusterLWP.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperClusterLWP.clustersBegin(); iBC != iDeepSuperClusterLWP.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperClusterLWP_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperClusterLWP_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global deepSuperClusterLWP indexing for EE has an offset = ndeepSuperClusterLWPEB
      iSC = (deepSuperClusterLWPEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "deepSuperClusterLWPsEE size: " << (deepSuperClusterLWPEE.product())->size() << std::endl;
      for(const auto& iDeepSuperClusterLWP : *(deepSuperClusterLWPEE.product())){    

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();  
          iSC_tmp++;
        
          deepSuperClusterLWP_energy.push_back(reduceFloat(iDeepSuperClusterLWP.energy(),nBits_));
          deepSuperClusterLWP_eta.push_back(reduceFloat(iDeepSuperClusterLWP.eta(),nBits_));
          deepSuperClusterLWP_phi.push_back(reduceFloat(iDeepSuperClusterLWP.phi(),nBits_));
          deepSuperClusterLWP_etaWidth.push_back(reduceFloat(iDeepSuperClusterLWP.etaWidth(),nBits_));
          deepSuperClusterLWP_phiWidth.push_back(reduceFloat(iDeepSuperClusterLWP.phiWidth(),nBits_));
          deepSuperClusterLWP_R.push_back(reduceFloat(iDeepSuperClusterLWP.position().R(),nBits_)); 
          deepSuperClusterLWP_nPFClusters.push_back(iDeepSuperClusterLWP.clusters().size());  
          math::XYZPoint caloPos = iDeepSuperClusterLWP.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          deepSuperClusterLWP_ieta.push_back(ee_id.ix());
          deepSuperClusterLWP_iphi.push_back(ee_id.iy());
          deepSuperClusterLWP_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iDeepSuperClusterLWP.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             deepSuperClusterLWP_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperClusterLWP_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperClusterLWP_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperClusterLWP_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperClusterLWP_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperClusterLWP_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperClusterLWP_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperClusterLWP_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperClusterLWP_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperClusterLWP_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperClusterLWP_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperClusterLWP_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperClusterLWP_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperClusterLWP_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperClusterLWP_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperClusterLWP_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             deepSuperClusterLWP_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperClusterLWP_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperClusterLWP_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperClusterLWP_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperClusterLWP_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperClusterLWP_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperClusterLWP_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperClusterLWP_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperClusterLWP_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperClusterLWP_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperClusterLWP_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperClusterLWP_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperClusterLWP_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperClusterLWP_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperClusterLWP_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperClusterLWP_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperClusterLWP_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperClusterLWP_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperClusterLWP_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             deepSuperClusterLWP_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperClusterLWP_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperClusterLWP_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_));

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperClusterLWP, towerIso1_, towerIso2_, egammaHadTower_);
             deepSuperClusterLWP_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperClusterLWP_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));   
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             deepSuperClusterLWP_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) deepSuperClusterLWP_dR_genScore_MatchedIndex.push_back(-1);
             else deepSuperClusterLWP_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_DeepSuperClusterLWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterLWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterLWP.eta(),iDeepSuperClusterLWP.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 
             
             deepSuperClusterLWP_dR_simScore[iSC] = dR_simScore;  
             deepSuperClusterLWP_sim_fraction_old[iSC] = sim_fraction_old;   
             deepSuperClusterLWP_simScore[iSC] = simScore;  
             deepSuperClusterLWP_n_shared_xtals[iSC] = n_shared_xtals;  
             deepSuperClusterLWP_sim_fraction[iSC] = sim_fraction;  
             deepSuperClusterLWP_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             deepSuperClusterLWP_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             deepSuperClusterLWP_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             deepSuperClusterLWP_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             deepSuperClusterLWP_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             deepSuperClusterLWP_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             deepSuperClusterLWP_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             deepSuperClusterLWP_sim_rechit_diff[iSC] = sim_rechit_diff; 
             deepSuperClusterLWP_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             deepSuperClusterLWP_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             deepSuperClusterLWP_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             deepSuperClusterLWP_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 

             deepSuperClusterLWP_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             deepSuperClusterLWP_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") deepSuperClusterLWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else deepSuperClusterLWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   deepSuperClusterLWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({deepSuperClusterLWP_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({deepSuperClusterLWP_sim_fraction_old, deepSuperClusterLWP_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iSC));                     
                } 
             }else{
                deepSuperClusterLWP_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                deepSuperClusterLWP_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterLWP_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterLWP_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));           
             }
          }

          if(iDeepSuperClusterLWP.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iDeepSuperClusterLWP.preshowerClusters().size(); iPC++){
                  if(!iDeepSuperClusterLWP.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  deepSuperClusterLWP_psCluster_energy[iSC_tmp].push_back(reduceFloat(iDeepSuperClusterLWP.preshowerClusters()[iPC]->energy(),nBits_));
                  deepSuperClusterLWP_psCluster_eta[iSC_tmp].push_back(reduceFloat(iDeepSuperClusterLWP.preshowerClusters()[iPC]->eta(),nBits_));
                  deepSuperClusterLWP_psCluster_phi[iSC_tmp].push_back(reduceFloat(iDeepSuperClusterLWP.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_){   
             //save clusters and deepSuperClusterLWPs mutual info
             reco::CaloCluster caloSeed(*iDeepSuperClusterLWP.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperClusterLWP.clustersBegin(); iBC != iDeepSuperClusterLWP.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperClusterLWP_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperClusterLWP_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }
   }

   //save pfCluster_deepSuperClusterLWPsIndex
   if(savePFCluster_ && saveSuperCluster_ && useDeepSC_){
      for(unsigned int iSC=0; iSC<deepSuperClusterLWP_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<deepSuperClusterLWP_pfClustersIndex.at(iSC).size(); iPF++)
              if(deepSuperClusterLWP_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_deepSuperClusterLWPsIndex[deepSuperClusterLWP_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
    }

   //save inverse of matchings
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&genParticle_deepSuperClusterLWP_dR_genScore_MatchedIndex,&deepSuperClusterLWP_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_dR_simScore_MatchedIndex,&deepSuperClusterLWP_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_old_MatchedIndex,&deepSuperClusterLWP_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_simScore_MatchedIndex,&deepSuperClusterLWP_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_n_shared_xtals_MatchedIndex,&deepSuperClusterLWP_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_MatchedIndex,&deepSuperClusterLWP_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex,&deepSuperClusterLWP_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_rechit_diff_MatchedIndex,&deepSuperClusterLWP_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex,&deepSuperClusterLWP_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex,&deepSuperClusterLWP_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex,&deepSuperClusterLWP_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex,&deepSuperClusterLWP_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex,&deepSuperClusterLWP_sim_rechit_combined_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex,&deepSuperClusterLWP_rechit_sim_combined_fraction_MatchedIndex);
      }      
   }

   //Save deepSuperClusterTWPs 
   //std::cout << "-----> deepSuperClusterTWPs <-----" << std::endl; 
   locCov_.clear();
   full5x5_locCov_.clear();
   if(saveSuperCluster_ && useDeepSC_){
      int iSC=0;
      //std::cout << "deepSuperClusterTWPsEB size: " << (deepSuperClusterTWPEB.product())->size() << std::endl;
      for(const auto& iDeepSuperClusterTWP : *(deepSuperClusterTWPEB.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();   

          deepSuperClusterTWP_energy.push_back(reduceFloat(iDeepSuperClusterTWP.energy(),nBits_));
          deepSuperClusterTWP_eta.push_back(reduceFloat(iDeepSuperClusterTWP.eta(),nBits_));
          deepSuperClusterTWP_phi.push_back(reduceFloat(iDeepSuperClusterTWP.phi(),nBits_));
          deepSuperClusterTWP_etaWidth.push_back(reduceFloat(iDeepSuperClusterTWP.etaWidth(),nBits_));
          deepSuperClusterTWP_phiWidth.push_back(reduceFloat(iDeepSuperClusterTWP.phiWidth(),nBits_));
          deepSuperClusterTWP_R.push_back(reduceFloat(iDeepSuperClusterTWP.position().R(),nBits_));
          deepSuperClusterTWP_nPFClusters.push_back(iDeepSuperClusterTWP.clusters().size());
          math::XYZPoint caloPos = iDeepSuperClusterTWP.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          deepSuperClusterTWP_ieta.push_back(eb_id.ieta());
          deepSuperClusterTWP_iphi.push_back(eb_id.iphi());
          deepSuperClusterTWP_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iDeepSuperClusterTWP.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEB.product())), topology);  
             deepSuperClusterTWP_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperClusterTWP_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperClusterTWP_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperClusterTWP_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperClusterTWP_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperClusterTWP_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperClusterTWP_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperClusterTWP_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperClusterTWP_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperClusterTWP_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperClusterTWP_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperClusterTWP_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperClusterTWP_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperClusterTWP_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperClusterTWP_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperClusterTWP_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             deepSuperClusterTWP_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperClusterTWP_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperClusterTWP_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperClusterTWP_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperClusterTWP_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperClusterTWP_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperClusterTWP_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperClusterTWP_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperClusterTWP_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperClusterTWP_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperClusterTWP_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperClusterTWP_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperClusterTWP_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperClusterTWP_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperClusterTWP_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperClusterTWP_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperClusterTWP_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperClusterTWP_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperClusterTWP_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             deepSuperClusterTWP_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperClusterTWP_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperClusterTWP_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperClusterTWP, towerIso1_, towerIso2_, egammaHadTower_);
             deepSuperClusterTWP_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperClusterTWP_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));  
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             deepSuperClusterTWP_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) deepSuperClusterTWP_dR_genScore_MatchedIndex.push_back(-1);
             else deepSuperClusterTWP_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_DeepSuperClusterTWPEB.at(iSC), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEB.at(iSC), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEB.at(iSC), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEB.at(iSC), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEB.at(iSC), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEB.at(iSC), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 

             deepSuperClusterTWP_dR_simScore[iSC] = dR_simScore;  
             deepSuperClusterTWP_sim_fraction_old[iSC] = sim_fraction_old;   
             deepSuperClusterTWP_simScore[iSC] = simScore;  
             deepSuperClusterTWP_n_shared_xtals[iSC] = n_shared_xtals;  
             deepSuperClusterTWP_sim_fraction[iSC] = sim_fraction;  
             deepSuperClusterTWP_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             deepSuperClusterTWP_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             deepSuperClusterTWP_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             deepSuperClusterTWP_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             deepSuperClusterTWP_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             deepSuperClusterTWP_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             deepSuperClusterTWP_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             deepSuperClusterTWP_sim_rechit_diff[iSC] = sim_rechit_diff; 
             deepSuperClusterTWP_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             deepSuperClusterTWP_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             deepSuperClusterTWP_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             deepSuperClusterTWP_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 
            
             deepSuperClusterTWP_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             deepSuperClusterTWP_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") deepSuperClusterTWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else deepSuperClusterTWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   deepSuperClusterTWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({deepSuperClusterTWP_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({deepSuperClusterTWP_sim_fraction_old, deepSuperClusterTWP_global_sim_rechit_fraction}), std::vector<double>({0.8,0.5}), iSC));      
                } 
             }else{
                deepSuperClusterTWP_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                deepSuperClusterTWP_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterTWP_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));        
             }             
          }

          if(savePFCluster_){   
             //save clusters and deepSuperClusterTWPs mutual info
             reco::CaloCluster caloSeed(*iDeepSuperClusterTWP.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperClusterTWP.clustersBegin(); iBC != iDeepSuperClusterTWP.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC)); 
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperClusterTWP_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperClusterTWP_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      } 

      // The global deepSuperClusterTWP indexing for EE has an offset = ndeepSuperClusterTWPEB
      iSC = (deepSuperClusterTWPEB.product())->size();
      int iSC_tmp=-1;
      //std::cout << "deepSuperClusterTWPsEE size: " << (deepSuperClusterTWPEE.product())->size() << std::endl;
      for(const auto& iDeepSuperClusterTWP : *(deepSuperClusterTWPEE.product())){    

          dR_genScore.clear();
          dR_simScore.clear();
          sim_fraction_old.clear();
          simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_1MeVCut.clear();
          sim_fraction_5MeVCut.clear();
          sim_fraction_10MeVCut.clear();
          sim_fraction_50MeVCut.clear();
          sim_fraction_100MeVCut.clear();  
          sim_fraction_500MeVCut.clear(); 
          sim_fraction_1GeVCut.clear();    
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          hgcal_caloToCluster.clear();
          hgcal_clusterToCalo.clear();  
          iSC_tmp++;
        
          deepSuperClusterTWP_energy.push_back(reduceFloat(iDeepSuperClusterTWP.energy(),nBits_));
          deepSuperClusterTWP_eta.push_back(reduceFloat(iDeepSuperClusterTWP.eta(),nBits_));
          deepSuperClusterTWP_phi.push_back(reduceFloat(iDeepSuperClusterTWP.phi(),nBits_));
          deepSuperClusterTWP_etaWidth.push_back(reduceFloat(iDeepSuperClusterTWP.etaWidth(),nBits_));
          deepSuperClusterTWP_phiWidth.push_back(reduceFloat(iDeepSuperClusterTWP.phiWidth(),nBits_));
          deepSuperClusterTWP_R.push_back(reduceFloat(iDeepSuperClusterTWP.position().R(),nBits_)); 
          deepSuperClusterTWP_nPFClusters.push_back(iDeepSuperClusterTWP.clusters().size());  
          math::XYZPoint caloPos = iDeepSuperClusterTWP.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          deepSuperClusterTWP_ieta.push_back(ee_id.ix());
          deepSuperClusterTWP_iphi.push_back(ee_id.iy());
          deepSuperClusterTWP_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iDeepSuperClusterTWP.seed());  
             showerShapes_.clear();
             showerShapes_ = getShowerShapes(&caloBC, &(*(recHitsEE.product())), topology);  
             deepSuperClusterTWP_e5x5.push_back(reduceFloat(showerShapes_[0],nBits_));
             deepSuperClusterTWP_e2x2Ratio.push_back(reduceFloat(showerShapes_[1],nBits_));
             deepSuperClusterTWP_e3x3Ratio.push_back(reduceFloat(showerShapes_[2],nBits_));
      	     deepSuperClusterTWP_eMaxRatio.push_back(reduceFloat(showerShapes_[3],nBits_));
             deepSuperClusterTWP_e2ndRatio.push_back(reduceFloat(showerShapes_[4],nBits_));
             deepSuperClusterTWP_eTopRatio.push_back(reduceFloat(showerShapes_[5],nBits_));
             deepSuperClusterTWP_eRightRatio.push_back(reduceFloat(showerShapes_[6],nBits_));
             deepSuperClusterTWP_eBottomRatio.push_back(reduceFloat(showerShapes_[7],nBits_));
             deepSuperClusterTWP_eLeftRatio.push_back(reduceFloat(showerShapes_[8],nBits_));
             deepSuperClusterTWP_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[9],nBits_));
             deepSuperClusterTWP_e2x5TopRatio.push_back(reduceFloat(showerShapes_[10],nBits_));
             deepSuperClusterTWP_e2x5RightRatio.push_back(reduceFloat(showerShapes_[11],nBits_));
             deepSuperClusterTWP_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[12],nBits_));
             deepSuperClusterTWP_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[13],nBits_));
             deepSuperClusterTWP_swissCross.push_back(reduceFloat(showerShapes_[14],nBits_));
             deepSuperClusterTWP_r9.push_back(reduceFloat(showerShapes_[15],nBits_));
             deepSuperClusterTWP_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[16],nBits_)); 
             deepSuperClusterTWP_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[17],nBits_)); 
             deepSuperClusterTWP_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[18],nBits_)); 
             deepSuperClusterTWP_full5x5_e5x5.push_back(reduceFloat(showerShapes_[19],nBits_));
             deepSuperClusterTWP_full5x5_e2x2Ratio.push_back(reduceFloat(showerShapes_[20],nBits_));
             deepSuperClusterTWP_full5x5_e3x3Ratio.push_back(reduceFloat(showerShapes_[21],nBits_));
             deepSuperClusterTWP_full5x5_eMaxRatio.push_back(reduceFloat(showerShapes_[22],nBits_));
             deepSuperClusterTWP_full5x5_e2ndRatio.push_back(reduceFloat(showerShapes_[23],nBits_));
             deepSuperClusterTWP_full5x5_eTopRatio.push_back(reduceFloat(showerShapes_[24],nBits_));
             deepSuperClusterTWP_full5x5_eRightRatio.push_back(reduceFloat(showerShapes_[25],nBits_));
             deepSuperClusterTWP_full5x5_eBottomRatio.push_back(reduceFloat(showerShapes_[26],nBits_));
             deepSuperClusterTWP_full5x5_eLeftRatio.push_back(reduceFloat(showerShapes_[27],nBits_));
             deepSuperClusterTWP_full5x5_e2x5MaxRatio.push_back(reduceFloat(showerShapes_[28],nBits_));
             deepSuperClusterTWP_full5x5_e2x5TopRatio.push_back(reduceFloat(showerShapes_[29],nBits_));
             deepSuperClusterTWP_full5x5_e2x5RightRatio.push_back(reduceFloat(showerShapes_[30],nBits_));
             deepSuperClusterTWP_full5x5_e2x5BottomRatio.push_back(reduceFloat(showerShapes_[31],nBits_));
             deepSuperClusterTWP_full5x5_e2x5LeftRatio.push_back(reduceFloat(showerShapes_[32],nBits_));
             deepSuperClusterTWP_full5x5_swissCross.push_back(reduceFloat(showerShapes_[33],nBits_));
             deepSuperClusterTWP_full5x5_r9.push_back(reduceFloat(showerShapes_[34],nBits_));
             deepSuperClusterTWP_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperClusterTWP_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperClusterTWP_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperClusterTWP, towerIso1_, towerIso2_, egammaHadTower_);
             deepSuperClusterTWP_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperClusterTWP_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())<0.1) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())); 
                 else dR_genScore.push_back(999.);  
             }  
             deepSuperClusterTWP_dR_genScore[iSC] = dR_genScore; 
             if(std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==-999.;}) || std::all_of(dR_genScore.begin(),dR_genScore.end(),[](double i){return i==999.;})) deepSuperClusterTWP_dR_genScore_MatchedIndex.push_back(-1);
             else deepSuperClusterTWP_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin()); 
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&hitsAndEnergies_DeepSuperClusterTWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_1MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_1MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_5MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_5MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_10MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_10MeVCut.at(iCalo), recHitsEB, recHitsEE); 
                 std::vector<double> scores_50MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_50MeVCut.at(iCalo), recHitsEB,recHitsEE);  
                 std::vector<double> scores_100MeVCut = getScores(&hitsAndEnergies_DeepSuperClusterTWPEE.at(iSC_tmp), &hitsAndEnergies_CaloPart_100MeVCut.at(iCalo), recHitsEB,recHitsEE);                
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())<0.1) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperClusterTWP.eta(),iDeepSuperClusterTWP.phi())); 
                 else dR_simScore.push_back(999.);
          
                 if(scoreType_=="n_shared_xtals") simScore.push_back(scores[0]);  
                 if(scoreType_=="sim_fraction") simScore.push_back(scores[1]);
                 if(scoreType_=="simScore_final_combination") simScore.push_back(scores[1]);   
                 if(scoreType_=="sim_fraction_1MeVCut") simScore.push_back(scores[10]);  
                 if(scoreType_=="sim_fraction_5MeVCut") simScore.push_back(scores[11]);  
                 if(scoreType_=="sim_fraction_10MeVCut") simScore.push_back(scores[12]);  
                 if(scoreType_=="sim_fraction_50MeVCut") simScore.push_back(scores[13]);
                 if(scoreType_=="sim_fraction_100MeVCut") simScore.push_back(scores[14]);    
                 if(scoreType_=="sim_fraction_500MeVCut") simScore.push_back(scores[15]);    
                 if(scoreType_=="sim_fraction_1GeVCut") simScore.push_back(scores[16]);      
                 if(scoreType_=="sim_rechit_diff") simScore.push_back(scores[2]); 
                 if(scoreType_=="sim_rechit_fraction") simScore.push_back(scores[3]);           
                 if(scoreType_=="global_sim_rechit_fraction") simScore.push_back(scores[4]);
                 if(scoreType_=="hgcal_caloToCluster") simScore.push_back(scores[7]);  
                 if(scoreType_=="hgcal_clusterToCalo") simScore.push_back(scores[8]);  
                 
                 sim_fraction_old.push_back(scores[9]);  
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_fraction_1MeVCut.push_back(scores[10]);  
                 sim_fraction_5MeVCut.push_back(scores[11]);  
                 sim_fraction_10MeVCut.push_back(scores[12]);  
                 sim_fraction_50MeVCut.push_back(scores[13]);
                 sim_fraction_100MeVCut.push_back(scores[14]);    
                 sim_fraction_500MeVCut.push_back(scores[15]);    
                 sim_fraction_1GeVCut.push_back(scores[16]);      
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 hgcal_caloToCluster.push_back(scores[7]);  
                 hgcal_clusterToCalo.push_back(scores[8]);       
             } 
             
             deepSuperClusterTWP_dR_simScore[iSC] = dR_simScore;  
             deepSuperClusterTWP_sim_fraction_old[iSC] = sim_fraction_old;   
             deepSuperClusterTWP_simScore[iSC] = simScore;  
             deepSuperClusterTWP_n_shared_xtals[iSC] = n_shared_xtals;  
             deepSuperClusterTWP_sim_fraction[iSC] = sim_fraction;  
             deepSuperClusterTWP_sim_fraction_1MeVCut[iSC] = sim_fraction_1MeVCut; 
             deepSuperClusterTWP_sim_fraction_5MeVCut[iSC] = sim_fraction_5MeVCut;  
             deepSuperClusterTWP_sim_fraction_10MeVCut[iSC] = sim_fraction_10MeVCut;  
             deepSuperClusterTWP_sim_fraction_50MeVCut[iSC] = sim_fraction_50MeVCut;  
             deepSuperClusterTWP_sim_fraction_100MeVCut[iSC] = sim_fraction_100MeVCut;       
             deepSuperClusterTWP_sim_fraction_500MeVCut[iSC] = sim_fraction_500MeVCut;       
             deepSuperClusterTWP_sim_fraction_1GeVCut[iSC] = sim_fraction_1GeVCut;            
             deepSuperClusterTWP_sim_rechit_diff[iSC] = sim_rechit_diff; 
             deepSuperClusterTWP_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
             deepSuperClusterTWP_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
             deepSuperClusterTWP_hgcal_caloToCluster[iSC] = hgcal_caloToCluster; 
             deepSuperClusterTWP_hgcal_clusterToCalo[iSC] = hgcal_clusterToCalo; 

             deepSuperClusterTWP_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_dR_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             deepSuperClusterTWP_sim_fraction_old_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_old, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));
             if(!saveScores_){
                if(scoreType_!="simScore_final_combination"){ 
                   if(scoreType_=="sim_rechit_diff" || scoreType_=="sim_rechit_fraction" || scoreType_=="global_sim_rechit_fraction" || scoreType_=="hgcal_caloToCluster" || scoreType_=="hgcal_clusterToCalo" || scoreType_=="rechit_sim_combined_fraction" || scoreType_=="sim_rechit_combined_fraction") deepSuperClusterTWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_simScore, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                   else deepSuperClusterTWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_simScore, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                }else{
                   deepSuperClusterTWP_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_simScore, 0.04, true, std::vector<std::vector<std::vector<double>>>({deepSuperClusterTWP_sim_fraction_100MeVCut}), std::vector<double>({0.01}), std::vector<std::vector<std::vector<double>>>({deepSuperClusterTWP_sim_fraction_old, deepSuperClusterTWP_global_sim_rechit_fraction}), std::vector<double>({0.1,0.5}), iSC));                     
                } 
             }else{
                deepSuperClusterTWP_n_shared_xtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_n_shared_xtals, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));     
                deepSuperClusterTWP_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_1MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_5MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_10MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));   
                deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_50MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));    
                deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_100MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_500MeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_fraction_1GeVCut, -999., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));          
                deepSuperClusterTWP_sim_rechit_diff_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_rechit_diff, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_global_sim_rechit_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC)); 
                deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_hgcal_caloToCluster, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));      
                deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_hgcal_clusterToCalo, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_rechit_sim_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));  
                deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperClusterTWP_sim_rechit_combined_fraction, 999., false, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iSC));           
             }
          }

          if(iDeepSuperClusterTWP.preshowerClusters().isAvailable()){
              for(unsigned int iPC=0; iPC<iDeepSuperClusterTWP.preshowerClusters().size(); iPC++){
                  if(!iDeepSuperClusterTWP.preshowerClusters()[iPC].isAvailable()) { continue; } 
                  deepSuperClusterTWP_psCluster_energy[iSC_tmp].push_back(reduceFloat(iDeepSuperClusterTWP.preshowerClusters()[iPC]->energy(),nBits_));
                  deepSuperClusterTWP_psCluster_eta[iSC_tmp].push_back(reduceFloat(iDeepSuperClusterTWP.preshowerClusters()[iPC]->eta(),nBits_));
                  deepSuperClusterTWP_psCluster_phi[iSC_tmp].push_back(reduceFloat(iDeepSuperClusterTWP.preshowerClusters()[iPC]->phi(),nBits_));   
              }
          } 

          if(savePFCluster_){   
             //save clusters and deepSuperClusterTWPs mutual info
             reco::CaloCluster caloSeed(*iDeepSuperClusterTWP.seed());  
             for(reco::CaloCluster_iterator iBC = iDeepSuperClusterTWP.clustersBegin(); iBC != iDeepSuperClusterTWP.clustersEnd(); ++iBC){
                 reco::CaloCluster caloSCluster(*(*iBC));  
                 int iPF=0;   
                 for(const auto& iPFCluster : *(pfClusters.product())){
                     reco::CaloCluster caloPFCluster(iPFCluster);
                     if(caloPFCluster == caloSCluster) deepSuperClusterTWP_pfClustersIndex[iSC].push_back(iPF); 
                     if(caloPFCluster == caloSCluster && caloSCluster == caloSeed) deepSuperClusterTWP_seedIndex[iSC]=iPF;   
                     iPF++;   
                 }     
             }      
          }
          iSC++;  
      }
   }

   //save pfCluster_deepSuperClusterTWPsIndex
   if(savePFCluster_ && saveSuperCluster_ && useDeepSC_){
      for(unsigned int iSC=0; iSC<deepSuperClusterTWP_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<deepSuperClusterTWP_pfClustersIndex.at(iSC).size(); iPF++)
              if(deepSuperClusterTWP_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_deepSuperClusterTWPsIndex[deepSuperClusterTWP_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
    }

   //save inverse of matchings
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&genParticle_deepSuperClusterTWP_dR_genScore_MatchedIndex,&deepSuperClusterTWP_dR_genScore_MatchedIndex);
   } 
   if(saveCaloParticles_ && saveSuperCluster_ && useDeepSC_){ 
      fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_dR_simScore_MatchedIndex,&deepSuperClusterTWP_dR_simScore_MatchedIndex);
      fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_old_MatchedIndex,&deepSuperClusterTWP_sim_fraction_old_MatchedIndex);
      if(!saveScores_){
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_simScore_MatchedIndex,&deepSuperClusterTWP_simScore_MatchedIndex);
      }else{
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_n_shared_xtals_MatchedIndex,&deepSuperClusterTWP_n_shared_xtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_MatchedIndex,&deepSuperClusterTWP_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_1MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_5MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_10MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_50MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_100MeVCut_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_500MeVCut_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex,&deepSuperClusterTWP_sim_fraction_1GeVCut_MatchedIndex);      
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_rechit_diff_MatchedIndex,&deepSuperClusterTWP_sim_rechit_diff_MatchedIndex);     
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex,&deepSuperClusterTWP_sim_rechit_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex,&deepSuperClusterTWP_global_sim_rechit_fraction_MatchedIndex);   
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex,&deepSuperClusterTWP_hgcal_caloToCluster_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex,&deepSuperClusterTWP_hgcal_clusterToCalo_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex,&deepSuperClusterTWP_sim_rechit_combined_fraction_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex,&deepSuperClusterTWP_rechit_sim_combined_fraction_MatchedIndex);
      }      
   }
  
   //Save unClustered pfRechits 
   pfRechit_unClustered.clear();
   if(savePFRechits_ || saveRechits_){
      for(const auto& iPFRechit : *(pfRecHits.product())){

          DetId pf_id(iPFRechit.detId());
          bool pfRecHit_isMatched_ = false;
          
          for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){ 
              for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){ 
                  if(iPFRechit.detId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()) pfRecHit_isMatched_ = true;
                  break;   
              }
          }

          if(pf_id.subdetId()==EcalBarrel){
             for(unsigned int iSC = 0; iSC < hitsAndEnergies_SuperClusterEB.size(); iSC++){
                 for(unsigned int i = 0; i < hitsAndEnergies_SuperClusterEB.at(iSC).size(); i++){ 
                     if(iPFRechit.detId() == hitsAndEnergies_SuperClusterEB.at(iSC).at(i).first.rawId()) pfRecHit_isMatched_ = true;
                     break;                   
                 }
             }       
          }else if(pf_id.subdetId()==EcalEndcap){
             for(unsigned int iSC = 0; iSC < hitsAndEnergies_SuperClusterEE.size(); iSC++){
                 for(unsigned int i = 0; i < hitsAndEnergies_SuperClusterEE.at(iSC).size(); i++){ 
                     if(iPFRechit.detId() == hitsAndEnergies_SuperClusterEE.at(iSC).at(i).first.rawId()) pfRecHit_isMatched_ = true;
                     break;                   
                 }
             }        
          }

          if(pfRecHit_isMatched_) continue;

          pfRechit_unClustered.push_back(pf_id); 

          cell = geometry->getPosition(pf_id); 
          pfRecHit_unClustered_energy.push_back(reduceFloat(iPFRechit.energy(),nBits_));    
          pfRecHit_unClustered_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          pfRecHit_unClustered_phi.push_back(reduceFloat(cell.phi(),nBits_)); 
          if(pf_id.subdetId()==EcalBarrel){ 
             EBDetId eb_id(pf_id);  
             pfRecHit_unClustered_ieta.push_back(eb_id.ieta());  
             pfRecHit_unClustered_iphi.push_back(eb_id.iphi());  
             pfRecHit_unClustered_iz.push_back(0);     
          }else if(pf_id.subdetId()==EcalEndcap){
             int iz=-99;
             EEDetId ee_id(pf_id);  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1; 
             pfRecHit_unClustered_ieta.push_back(ee_id.ix());  
             pfRecHit_unClustered_iphi.push_back(ee_id.iy());  
             pfRecHit_unClustered_iz.push_back(iz);    
          } 
      }   
   }  
   
   //Save noPF rechits 
   if(saveRechits_){
      for(const auto& iRechit : *(recHitsEB.product())){
          
          DetId rechit_id(iRechit.detid());

          bool rechit_isMatched_;
          for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){ 
              for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){ 
                  if(rechit_id.rawId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()) rechit_isMatched_ = true;
                  break;   
              }
          }
          if(rechit_isMatched_) continue;  
           
          std::vector<DetId>::iterator it = std::find(pfRechit_unClustered.begin(), pfRechit_unClustered.end(), rechit_id);   
          if (it != pfRechit_unClustered.end()) continue;  
          
          cell = geometry->getPosition(rechit_id); 
          recHit_noPF_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          recHit_noPF_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          recHit_noPF_phi.push_back(reduceFloat(cell.phi(),nBits_)); 

          EBDetId eb_id(rechit_id);  
          recHit_noPF_ieta.push_back(eb_id.ieta());  
          recHit_noPF_iphi.push_back(eb_id.iphi());  
          recHit_noPF_iz.push_back(0);     
      }
      for(const auto& iRechit : *(recHitsEE.product())){

          DetId rechit_id(iRechit.detid());

          bool rechit_isMatched_;
          for(unsigned int iPFCl = 0; iPFCl < hitsAndEnergies_PFCluster.size(); iPFCl++){ 
              for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){ 
                  if(rechit_id.rawId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId()) rechit_isMatched_ = true;
                  break;   
              }
          }
          if(rechit_isMatched_) continue; 
          
          std::vector<DetId>::iterator it = std::find(pfRechit_unClustered.begin(), pfRechit_unClustered.end(), rechit_id);   
          if (it != pfRechit_unClustered.end()) continue;  
          
          cell = geometry->getPosition(rechit_id); 
          recHit_noPF_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          recHit_noPF_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          recHit_noPF_phi.push_back(reduceFloat(cell.phi(),nBits_)); 

          int iz=-99;
          EEDetId ee_id(rechit_id);  
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1; 
          recHit_noPF_ieta.push_back(ee_id.ix());  
          recHit_noPF_iphi.push_back(ee_id.iy());  
          recHit_noPF_iz.push_back(iz);    
      }   
   }  
   //fill tree for each event
   tree->Fill();
}

void RecoSimDumper::beginJob()
{

}

void RecoSimDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<float> RecoSimDumper::getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology)
{
    std::vector<float> shapes;
    shapes.resize(38); 
    locCov_.clear();
    full5x5_locCov_.clear();
    locCov_ = EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    full5x5_locCov_ = noZS::EcalClusterTools::localCovariances(*caloBC, recHits, topology);
    
    float e5x5 = EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
    float e3x3 = EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
    float eMax = EcalClusterTools::eMax(*caloBC, recHits); // eMax
    float eTop = EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
    float eRight = EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
    float eBottom = EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
    float eLeft = EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
    float e4 = eTop + eRight + eBottom + eLeft;

    shapes[0] = e5x5;
    shapes[1] = EcalClusterTools::e2x2(*caloBC, recHits, topology)/e5x5; // e2x2/e5x5
    shapes[2] = EcalClusterTools::e3x3(*caloBC, recHits, topology)/e5x5; // e3x3/e5x5
    shapes[3] = EcalClusterTools::eMax(*caloBC,  recHits)/e5x5; // eMax/e5x5
    shapes[4] = EcalClusterTools::e2nd(*caloBC, recHits)/e5x5; // e2nd/e5x5
    shapes[5] = EcalClusterTools::eTop(*caloBC, recHits, topology)/e5x5; // eTop/e5x5 
    shapes[6] = EcalClusterTools::eRight(*caloBC, recHits, topology)/e5x5; // eRight/e5x5
    shapes[7] = EcalClusterTools::eBottom(*caloBC, recHits, topology)/e5x5; // eBottom/e5x5
    shapes[8] = EcalClusterTools::eLeft(*caloBC, recHits, topology)/e5x5; // eLeft/e5x5
    shapes[9] = EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/e5x5; // e2x5Max/e5x5
    shapes[10] = EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/e5x5; // e2x5Top/e5x5  
    shapes[11] = EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/e5x5; // e2x5Bottom/e5x5  
    shapes[12] = EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/e5x5; // e2x5Left/e5x5  
    shapes[13] = EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/e5x5; // e2x5Right/e5x5   
    shapes[14] = 1.-e4/eMax; // swissCross 
    shapes[15] = e3x3/caloBC->energy(); // r9     
    shapes[16] = sqrt(locCov_[0]); // sigmaIetaIeta 
    shapes[17] = sqrt(locCov_[1]); // sigmaIetaIphi
    shapes[18] = !edm::isFinite(locCov_[2]) ? 0. : sqrt(locCov_[2]); // sigmaIphiIphi 

    // full_5x5 variables
    float full5x5_e5x5 = noZS::EcalClusterTools::e5x5(*caloBC, recHits, topology); // e5x5
    float full5x5_e3x3 = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology); // e3x3
    float full5x5_eMax = noZS::EcalClusterTools::eMax(*caloBC, recHits); // eMax
    float full5x5_eTop = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology); // eTop 
    float full5x5_eRight = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology); // eRight
    float full5x5_eBottom = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology); // eBottom
    float full5x5_eLeft = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology); // eLeft
    float full5x5_e4 = full5x5_eTop + full5x5_eRight + full5x5_eBottom + full5x5_eLeft;

    shapes[19] = full5x5_e5x5;
    shapes[20] = noZS::EcalClusterTools::e2x2(*caloBC, recHits, topology)/full5x5_e5x5; // e2x2/e5x5
    shapes[21] = noZS::EcalClusterTools::e3x3(*caloBC, recHits, topology)/full5x5_e5x5; // e3x3/e5x5
    shapes[22] = noZS::EcalClusterTools::eMax(*caloBC, recHits)/full5x5_e5x5; // eMax/e5x5
    shapes[23] = noZS::EcalClusterTools::e2nd(*caloBC, recHits)/full5x5_e5x5; // e2nd/e5x5
    shapes[24] = noZS::EcalClusterTools::eTop(*caloBC, recHits, topology)/full5x5_e5x5; // eTop/e5x5 
    shapes[25] = noZS::EcalClusterTools::eRight(*caloBC, recHits, topology)/full5x5_e5x5; // eRight/e5x5
    shapes[26] = noZS::EcalClusterTools::eBottom(*caloBC, recHits, topology)/full5x5_e5x5; // eBottom/e5x5
    shapes[27] = noZS::EcalClusterTools::eLeft(*caloBC, recHits, topology)/full5x5_e5x5; // eLeft/e5x5
    shapes[28] = noZS::EcalClusterTools::e2x5Max(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Max/e5x5
    shapes[29] = noZS::EcalClusterTools::e2x5Top(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Top/e5x5  
    shapes[30] = noZS::EcalClusterTools::e2x5Right(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Bottom/e5x5  
    shapes[31] = noZS::EcalClusterTools::e2x5Bottom(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Left/e5x5  
    shapes[32] = noZS::EcalClusterTools::e2x5Left(*caloBC, recHits, topology)/full5x5_e5x5; // e2x5Right/e5x5   
    shapes[33] = 1.-full5x5_e4/full5x5_eMax; // swissCross 
    shapes[34] = full5x5_e3x3/caloBC->energy(); // r9
    shapes[35] = sqrt(full5x5_locCov_[0]); // sigmaIetaIeta        
    shapes[36] = sqrt(full5x5_locCov_[1]); // sigmaIetaIphi          
    shapes[37] = !edm::isFinite(full5x5_locCov_[2]) ? 0. : sqrt(full5x5_locCov_[2]); // sigmaIphiIphi 

    return shapes; 
}

std::vector<float> RecoSimDumper::getHoE(const reco::SuperCluster* iSuperCluster, EgammaTowerIsolation* towerIso1, EgammaTowerIsolation* towerIso2, const EgammaHadTower* egammaHadTower)
{
     std::vector<float> HoEs;
     HoEs.resize(2);
  
     //std::vector<CaloTowerDetId> towersBehindCluster = egammaHadTower->towersOf(*iSuperCluster);
     double HoEraw1 = towerIso1->getTowerESum(iSuperCluster)/iSuperCluster->rawEnergy();
     double HoEraw2 = towerIso2->getTowerESum(iSuperCluster)/iSuperCluster->rawEnergy();        
     //float HoEraw1bc = egammaHadTower->getDepth1HcalESum(towersBehindCluster)/iSuperCluster->energy();
     //float HoEraw2bc = egammaHadTower->getDepth2HcalESum(towersBehindCluster)/iSuperCluster->energy(); 
     HoEs[0] = HoEraw1 + HoEraw2;
     //HoEs[1] = HoEraw1bc + HoEraw2bc;
     HoEs[1] = -1.;

     return HoEs;
}

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle, float simHitEnergy_cut)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CaloPart_tmp = new std::vector<std::pair<DetId, float> >;
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;
    
    const auto& simClusters = iCaloParticle->simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int i = 0; i < hits_and_energies.size(); i++){ 
            if(hits_and_energies[i].second < simHitEnergy_cut) continue; 
            HitsAndEnergies_tmp->push_back(make_pair(DetId(hits_and_energies[i].first),hits_and_energies[i].second));  
        }  
    }

    for(unsigned int i = 0; i < HitsAndEnergies_tmp->size(); i++){  
        if (HitsAndEnergies_map.find(HitsAndEnergies_tmp->at(i).first) == HitsAndEnergies_map.end()) {
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_tmp->at(i).second;      
        }else{
            HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]=HitsAndEnergies_map[HitsAndEnergies_tmp->at(i).first]+HitsAndEnergies_tmp->at(i).second; 
        }
    }

    for(auto const& hit : HitsAndEnergies_map) 
         HitsAndEnergies_CaloPart_tmp->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_CaloPart_tmp;
}

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    
    const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster->hitsAndFractions();
    for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
        if(hitsAndFractions.at(i).first.subdetId()==EcalBarrel){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*(recHitsEB.product())->find(hitsAndFractions[i].first)).energy()));
        }else if(hitsAndFractions.at(i).first.subdetId()==EcalEndcap){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*(recHitsEE.product())->find(hitsAndFractions[i].first)).energy()));
        }
    }

    return HitsAndEnergies_tmp;
}


std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_SuperCluster_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;

    for(reco::CaloCluster_iterator iBC = iSuperCluster->clustersBegin(); iBC != iSuperCluster->clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < seedrechits.size(); i++){  
            if(seedrechits.at(i).first.subdetId()==EcalBarrel){   
                if (HitsAndEnergies_map.find(seedrechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[seedrechits.at(i).first]=seedrechits.at(i).second * (*(recHitsEB.product())->find(seedrechits.at(i).first)).energy();    
                }else{
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*(recHitsEB.product())->find(seedrechits.at(i).first)).energy();
                } 
            }else if(seedrechits.at(i).first.subdetId()==EcalEndcap){   
                if (HitsAndEnergies_map.find(seedrechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[seedrechits.at(i).first]=seedrechits.at(i).second * (*(recHitsEE.product())->find(seedrechits.at(i).first)).energy();   
                }else{
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*(recHitsEE.product())->find(seedrechits.at(i).first)).energy();
                } 
            }
        }                      
    } 

    for(auto const& hit : HitsAndEnergies_map) 
        HitsAndEnergies_SuperCluster_tmp->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_SuperCluster_tmp;
}

std::vector<double> RecoSimDumper::getScores(const std::vector<std::pair<DetId, float> >*hits_and_energies_Cluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<double> scores;
    scores.resize(19);

    double nSharedXtals=0;
    double simFraction=0.;
    double simFraction_old=0.;
    double sim_rechit_diff=0.;
    double sim_rechit_fraction=0.;     
    double global_sim_rechit_fraction=0.;  
    double hgcal_caloToCluster=0.;      
    double hgcal_clusterToCalo=0.;   
    double sim_rechit_combined_fraction=0.;
    double rechit_sim_combined_fraction=0.;    
   
    double rechits_tot_CaloPart = 0.;
    double rechits_tot_CaloPart_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart) {
        rechits_tot_CaloPart+=hit_CaloPart.second;
        rechits_tot_CaloPart_noEnergy+=1.;
    }

    double rechits_tot_Cluster = 0.;
    double rechits_tot_Cluster_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster) {
        rechits_tot_Cluster+=hit_Cluster.second;
        rechits_tot_Cluster_noEnergy+=1.;
    }
   
    double rechits_match_Cluster = 0.;
    double rechits_match_CaloPart = 0.;
    double rechits_match_CaloPart_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){     
        for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){  
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){

               rechits_match_Cluster += hit_Cluster.second;
               rechits_match_CaloPart += hit_CaloPart.second;    
               rechits_match_CaloPart_noEnergy += 1.0;

               sim_rechit_diff += fabs(hit_CaloPart.second-hit_Cluster.second);

               double reco_ratio=0.; 
               double sim_ratio = 0.;  
               if(rechits_tot_Cluster!=0.) reco_ratio = (double)hit_Cluster.second/(double)rechits_tot_Cluster;     
               if(rechits_tot_CaloPart!=0.) sim_ratio = (double)hit_CaloPart.second/(double)rechits_tot_CaloPart; 
               sim_rechit_fraction += fabs(sim_ratio - reco_ratio);    
            }         
        }
    }

    double hgcal_caloToCluster_Num = 0.;
    double hgcal_caloToCluster_Denum = 0.;
    double hgcal_clusterToCalo_Num = 0.;
    double hgcal_clusterToCalo_Denum = 0.;
    double reco_fraction = 0.;
    double sim_fraction = 0.;
    
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){  

        double rechitE=0.;
        if(hit_CaloPart.first.subdetId()==EcalBarrel) rechitE = (*(recHitsEB.product())->find(hit_CaloPart.first)).energy();
        else if(hit_CaloPart.first.subdetId()==EcalEndcap) rechitE = (*(recHitsEE.product())->find(hit_CaloPart.first)).energy(); 
        sim_fraction = (double)hit_CaloPart.second/(double)rechits_tot_CaloPart;   

        reco_fraction = 0.;   
        for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){   
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               reco_fraction = (double)rechitE/(double)rechits_tot_Cluster; 
            }
        }

        hgcal_caloToCluster_Num += (reco_fraction-sim_fraction)*(reco_fraction-sim_fraction)*rechitE*rechitE; 
        hgcal_caloToCluster_Denum += sim_fraction*sim_fraction*rechitE*rechitE;      
    }

    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){  

        double rechitE=0.;
        if(hit_Cluster.first.subdetId()==EcalBarrel) rechitE = (*(recHitsEB.product())->find(hit_Cluster.first)).energy();
        else if(hit_Cluster.first.subdetId()==EcalEndcap) rechitE = (*(recHitsEE.product())->find(hit_Cluster.first)).energy(); 
        reco_fraction = (double)rechitE/(double)rechits_tot_Cluster;
        hgcal_clusterToCalo_Denum += reco_fraction*reco_fraction*rechitE*rechitE;  

        for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){
            reco_fraction = 0.;
            sim_fraction = 0.; 
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               reco_fraction = (double)rechitE/(double)rechits_tot_Cluster; 
               sim_fraction = (double)hit_CaloPart.second/(double)rechits_tot_CaloPart;      
               hgcal_clusterToCalo_Num += (reco_fraction-sim_fraction)*(reco_fraction-sim_fraction)*rechitE*rechitE;     
            }
        }   
    }

    nSharedXtals = (int)rechits_match_CaloPart_noEnergy;
    if(nSharedXtals==0) nSharedXtals=-999;

    if(rechits_tot_CaloPart!=0.) simFraction = (double)rechits_match_CaloPart/(double)rechits_tot_CaloPart;
    else simFraction=-999.; 
    if(simFraction==0.) simFraction=-999.;

    if(rechits_match_CaloPart!=0.) simFraction_old = fabs(1.-(double)rechits_match_Cluster/(double)rechits_match_CaloPart);
    else simFraction_old=999.; 
    if(simFraction_old==1.) simFraction_old=999.;

    if(rechits_match_CaloPart_noEnergy!=0.) sim_rechit_diff = (1./(double)rechits_match_CaloPart_noEnergy)*(double)sim_rechit_diff;
    else sim_rechit_diff=999.; 

    if(sim_rechit_fraction!=0. && rechits_match_CaloPart_noEnergy!=0.) sim_rechit_fraction = (1./(double)rechits_match_CaloPart_noEnergy)*sim_rechit_fraction;
    else sim_rechit_fraction=999.;
    
    if(rechits_tot_CaloPart!=0. && rechits_tot_Cluster!=0. && rechits_match_CaloPart!=0. && rechits_match_Cluster!=0. && rechits_match_CaloPart_noEnergy!=0.) global_sim_rechit_fraction = (1./(double)rechits_match_CaloPart_noEnergy)*fabs((double)rechits_match_CaloPart/(double)rechits_tot_CaloPart - (double)rechits_match_Cluster/(double)rechits_tot_Cluster);
    else global_sim_rechit_fraction=999.;  

    if(rechits_tot_Cluster!=0.) sim_rechit_combined_fraction = (double)rechits_match_CaloPart/(double)rechits_tot_Cluster;
    else sim_rechit_combined_fraction = 999.;
    if(rechits_tot_CaloPart!=0.) rechit_sim_combined_fraction = (double)rechits_match_Cluster/(double)rechits_tot_CaloPart;
    else rechit_sim_combined_fraction = 999.;   
 
    if(hgcal_caloToCluster_Denum!=0.) hgcal_caloToCluster = (double)hgcal_caloToCluster_Num/(double)hgcal_caloToCluster_Denum; 
    else hgcal_caloToCluster = 999.;
    if(hgcal_caloToCluster_Denum!=0.) hgcal_clusterToCalo = (double)hgcal_clusterToCalo_Num/(double)hgcal_caloToCluster_Denum;  
    else hgcal_caloToCluster = 999.;
    
    scores[0] = (double)nSharedXtals;
    scores[1] = simFraction;
    scores[2] = sim_rechit_diff;
    scores[3] = sim_rechit_fraction;     
    scores[4] = global_sim_rechit_fraction;  
    if(simFraction>0.01) scores[5] = simFraction; 
    else scores[5] = -999.; 
    if(simFraction>0.03) scores[6] = simFraction; 
    else scores[6] = -999.;
    scores[7] = hgcal_caloToCluster; 
    scores[8] = hgcal_clusterToCalo; 
    scores[9] = simFraction_old;  
    if((double)rechits_match_CaloPart>0.001) scores[10] = simFraction; 
    else scores[10] = -999.; 
    if((double)rechits_match_CaloPart>0.005) scores[11] = simFraction; 
    else scores[11] = -999.;  
    if((double)rechits_match_CaloPart>0.01) scores[12] = simFraction; 
    else scores[12] = -999.;  
    if((double)rechits_match_CaloPart>0.05) scores[13] = simFraction; 
    else scores[13] = -999.;   
    if((double)rechits_match_CaloPart>0.1) scores[14] = simFraction; 
    else scores[14] = -999.;  
    if((double)rechits_match_CaloPart>0.5) scores[15] = simFraction; 
    else scores[15] = -999.;  
    if((double)rechits_match_CaloPart>1.) scores[16] = simFraction; 
    else scores[16] = -999.;  
    scores[17] = sim_rechit_combined_fraction;
    scores[18] = rechit_sim_combined_fraction;

    return scores;
}

int RecoSimDumper::getMatchedIndex(std::vector<std::vector<double>>* score, double selection, bool useMax, std::vector<std::vector<std::vector<double>>> scoreSelMax, std::vector<double> selectionMax, std::vector<std::vector<std::vector<double>>> scoreSelMin, std::vector<double> selectionMin, int iCl)
{
   int matchedIndex = -1; 
   if(!useMax){ 
      std::replace(score->at(iCl).begin(),score->at(iCl).end(), -999., 999.);
      if(std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::min_element(score->at(iCl).begin(),score->at(iCl).end()) - score->at(iCl).begin();
   }else{
      std::replace(score->at(iCl).begin(),score->at(iCl).end(), 999., -999.); 
      if(std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::max_element(score->at(iCl).begin(),score->at(iCl).end()) - score->at(iCl).begin();
   }

   if(matchedIndex==-1) return -1;
   
   //std::cout << "getMatchedIndex - Max - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMax.size() << " - " << selectionMax.size() << std::endl;
   //std::cout << "getMatchedIndex - Min - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMin.size() << " - " << selectionMin.size() << std::endl;
    
   bool passSelection = true;
   for(unsigned int iSelMax=0; iSelMax < scoreSelMax.size(); iSelMax++)
       if(scoreSelMax.at(iSelMax).at(iCl).at(matchedIndex) < selectionMax.at(iSelMax)) passSelection = false;
      
   for(unsigned int iSelMin=0; iSelMin < scoreSelMin.size(); iSelMin++)
       if(scoreSelMin.at(iSelMin).at(iCl).at(matchedIndex) > selectionMin.at(iSelMin)) passSelection = false;
 
   if(useMax && score->at(iCl).at(matchedIndex) > selection && passSelection) return matchedIndex;
   if(!useMax && score->at(iCl).at(matchedIndex) < selection && passSelection) return matchedIndex; 
   
   return -1; 
}

void RecoSimDumper::fillParticleMatchedIndex(std::vector<std::vector<int>>* particleMatchedIndex, std::vector<int>* clusterMatchedIndex)
{
   for(unsigned int iCl=0; iCl<clusterMatchedIndex->size(); iCl++)
          if(clusterMatchedIndex->at(iCl)>=0) particleMatchedIndex->at(clusterMatchedIndex->at(iCl)).push_back(iCl);
}

GlobalPoint RecoSimDumper::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES)
{
  double preshowerStartEta = 1.653;
  double preshowerEndEta = 2.6;
  double cl_energy_float = 0;
  double max_e = 0.0;
  double clusterT0 = 0.0;
  DetId id_max; 
 
  // find the seed and max layer
  for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP){
    const double rh_energyf = hit_CP.second;
    cl_energy_float += rh_energyf;
    if (rh_energyf > max_e) {
      max_e = rh_energyf;
      id_max = hit_CP.first;
    }
  }

  const CaloSubdetectorGeometry* ecal_geom = nullptr;
  // get seed geometry information
  if(id_max.subdetId()==EcalBarrel){
      ecal_geom = _ebGeom;
      clusterT0 = _param_T0_EB;
  }else if(id_max.subdetId()==EcalEndcap){  
      ecal_geom = _eeGeom;
      clusterT0 = _param_T0_EE;
  }else{
      throw cms::Exception("InvalidLayer") << "ECAL Position Calc only accepts ECAL_BARREL or ECAL_ENDCAP";
  }

  auto center_cell = ecal_geom->getGeometry(id_max);
  const double ctreta = center_cell->etaPos();
  const double actreta = std::abs(ctreta);
  // need to change T0 if in ES
  if (actreta > preshowerStartEta && actreta < preshowerEndEta && useES) {
    if (ctreta > 0 && _esPlus)
      clusterT0 = _param_T0_ES;
    if (ctreta < 0 && _esMinus)
      clusterT0 = _param_T0_ES;
  }

  // floats to reproduce exactly the EGM code
  const float maxDepth = _param_X0 * (clusterT0 + log(cl_energy_float));
  const float maxToFront = center_cell->getPosition().mag();
  // calculate the position
  const double logETot_inv = -log(cl_energy_float);
  double position_norm = 0.0;
  double x(0.0), y(0.0), z(0.0);
 
  for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
    if(hit_CP.first.subdetId()!=EcalBarrel && hit_CP.first.subdetId()!=EcalEndcap) continue;
    auto cell = ecal_geom->getGeometry(hit_CP.first);
    if(!cell.get()) continue;
    double weight = 0.0;
    const double rh_energy = hit_CP.second;
    if (rh_energy > 0.0)
      weight = std::max(0.0, (_param_W0 + log(rh_energy) + logETot_inv));
    const float depth = maxDepth + maxToFront - cell->getPosition().mag();
    const GlobalPoint pos = static_cast<const TruncatedPyramid*>(cell.get())->getPosition(depth);
    
    x += weight * pos.x();
    y += weight * pos.y();
    z += weight * pos.z();

    position_norm += weight;
  }

  // FALL BACK to LINEAR WEIGHTS
  if (position_norm == 0.) {
    for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
      if(hit_CP.first.subdetId()!=EcalBarrel && hit_CP.first.subdetId()!=EcalEndcap) continue; 
      auto cell = ecal_geom->getGeometry(hit_CP.first);
      if(!cell.get()) continue;
      double weight = 0.0;
      const double rh_energy = hit_CP.second;
      if (rh_energy > 0.0)
        weight = rh_energy / cl_energy_float;

      const float depth = maxDepth + maxToFront - cell->getPosition().mag();
      const GlobalPoint pos = cell->getPosition(depth);

      x += weight * pos.x();
      y += weight * pos.y();
      z += weight * pos.z();

      position_norm += weight;
    }
  }

  if (position_norm < _minAllowedNorm) {
    edm::LogError("WeirdClusterNormalization") << "Cluster too far from seeding cell: set position to (0,0,0).";
    return GlobalPoint(0., 0., 0.);
  } else {
    const double norm_inverse = 1.0 / position_norm;
    x *= norm_inverse;
    y *= norm_inverse;
    z *= norm_inverse;
  }

  return GlobalPoint(x, y, z); 
}

float RecoSimDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumper);

