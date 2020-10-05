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
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
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
         tree->Branch("caloParticle_pfCluster_dR_simScore_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("caloParticle_pfCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_MatchedIndex);      
         tree->Branch("caloParticle_pfCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_recoToSim_fraction_MatchedIndex);   
         tree->Branch("caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
         tree->Branch("caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex);  
         tree->Branch("caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex);         
      }
      if(saveSuperCluster_){    
         tree->Branch("caloParticle_superCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >", &caloParticle_superCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("caloParticle_superCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_sim_fraction_MatchedIndex);      
         tree->Branch("caloParticle_superCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_recoToSim_fraction_MatchedIndex);   
         tree->Branch("caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
         tree->Branch("caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex);  
         tree->Branch("caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex);
         if(useRetunedSC_){
            tree->Branch("caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex);
            tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
            tree->Branch("caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex);      
            tree->Branch("caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex);   
            tree->Branch("caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
            tree->Branch("caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
            tree->Branch("caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex);
         } 
         if(useDeepSC_){
            tree->Branch("caloParticle_deepSuperCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex);
            tree->Branch("caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
            tree->Branch("caloParticle_deepSuperCluster_sim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_sim_fraction_MatchedIndex);      
            tree->Branch("caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex);   
            tree->Branch("caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);   
            tree->Branch("caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
            tree->Branch("caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex", "std::vector<std::vector<int> >", &caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex); 
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
      tree->Branch("pfCluster_rawEnergy","std::vector<float>",&pfCluster_rawEnergy);
      tree->Branch("pfCluster_energy","std::vector<float>",&pfCluster_energy);
      tree->Branch("pfCluster_rawPt","std::vector<float>",&pfCluster_rawPt); 
      tree->Branch("pfCluster_pt","std::vector<float>",&pfCluster_pt);  
      tree->Branch("pfCluster_eta","std::vector<float>",&pfCluster_eta);
      tree->Branch("pfCluster_phi","std::vector<float>",&pfCluster_phi); 
      tree->Branch("pfCluster_ieta","std::vector<int>",&pfCluster_ieta);
      tree->Branch("pfCluster_iphi","std::vector<int>",&pfCluster_iphi);   
      tree->Branch("pfCluster_iz","std::vector<int>",&pfCluster_iz);
      tree->Branch("pfCluster_nXtals","std::vector<int>",&pfCluster_nXtals);  
      if(saveSuperCluster_) tree->Branch("pfCluster_superClustersIndex","std::vector<std::vector<int> >",&pfCluster_superClustersIndex); 
      if(saveSuperCluster_ && useRetunedSC_) tree->Branch("pfCluster_retunedSuperClustersIndex","std::vector<std::vector<int> >",&pfCluster_retunedSuperClustersIndex);  
      if(saveSuperCluster_ && useDeepSC_) tree->Branch("pfCluster_deepSuperClustersIndex","std::vector<std::vector<int> >",&pfCluster_deepSuperClustersIndex); 
      if(saveCaloParticles_){ 
         tree->Branch("pfCluster_dR_genScore","std::vector<std::vector<double> >",&pfCluster_dR_genScore);
         tree->Branch("pfCluster_dR_simScore","std::vector<std::vector<double> >",&pfCluster_dR_simScore);
         tree->Branch("pfCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&pfCluster_sim_nSharedXtals);
         tree->Branch("pfCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&pfCluster_sim_fraction_noHitsFraction); 
         tree->Branch("pfCluster_sim_fraction","std::vector<std::vector<double> >",&pfCluster_sim_fraction);
         tree->Branch("pfCluster_recoToSim_fraction","std::vector<std::vector<double> >",&pfCluster_recoToSim_fraction);
         tree->Branch("pfCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&pfCluster_recoToSim_fraction_sharedXtals);
         tree->Branch("pfCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&pfCluster_simEnergy_sharedXtals);
         tree->Branch("pfCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&pfCluster_recoEnergy_sharedXtals);  
         tree->Branch("pfCluster_dR_genScore_MatchedIndex","std::vector<int>",&pfCluster_dR_genScore_MatchedIndex);
         tree->Branch("pfCluster_dR_simScore_MatchedIndex","std::vector<int>",&pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("pfCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&pfCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_noHitsFraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("pfCluster_sim_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_MatchedIndex);
         tree->Branch("pfCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&pfCluster_recoToSim_fraction_MatchedIndex);
         tree->Branch("pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex","std::vector<int>",&pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
         tree->Branch("pfCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&pfCluster_simEnergy_sharedXtals_MatchedIndex);
         tree->Branch("pfCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&pfCluster_recoEnergy_sharedXtals_MatchedIndex);
      }
      if(savePFClusterhits_){ 
         tree->Branch("pfClusterHit_fraction","std::vector<std::vector<float> >",&pfClusterHit_fraction);
         tree->Branch("pfClusterHit_rechitEnergy","std::vector<std::vector<float> >",&pfClusterHit_rechitEnergy);
         tree->Branch("pfClusterHit_eta","std::vector<std::vector<float> >",&pfClusterHit_eta);
         tree->Branch("pfClusterHit_phi","std::vector<std::vector<float> >",&pfClusterHit_phi);
         tree->Branch("pfClusterHit_ieta","std::vector<std::vector<int> >",&pfClusterHit_ieta);
         tree->Branch("pfClusterHit_iphi","std::vector<std::vector<int> >",&pfClusterHit_iphi);
         tree->Branch("pfClusterHit_iz","std::vector<std::vector<int> >",&pfClusterHit_iz);
      }   
   }
   if(saveSuperCluster_){
      tree->Branch("superCluster_rawEnergy","std::vector<float> ",&superCluster_rawEnergy);     
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
      tree->Branch("superCluster_nXtals","std::vector<int>",&superCluster_nXtals);   
      if(savePFCluster_) tree->Branch("superCluster_seedIndex","std::vector<int>",&superCluster_seedIndex);     
      if(savePFCluster_) tree->Branch("superCluster_pfClustersIndex","std::vector<std::vector<int> >",&superCluster_pfClustersIndex); 
      tree->Branch("superCluster_psCluster_energy", "std::vector<std::vector<float> >", &superCluster_psCluster_energy);
      tree->Branch("superCluster_psCluster_eta", "std::vector<std::vector<float> >", &superCluster_psCluster_eta);
      tree->Branch("superCluster_psCluster_phi", "std::vector<std::vector<float> >", &superCluster_psCluster_phi); 
      if(saveCaloParticles_){
         tree->Branch("superCluster_dR_genScore","std::vector<std::vector<double> >",&superCluster_dR_genScore);
         tree->Branch("superCluster_dR_simScore","std::vector<std::vector<double> >",&superCluster_dR_simScore);
         tree->Branch("superCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&superCluster_sim_nSharedXtals);
         tree->Branch("superCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&superCluster_sim_fraction_noHitsFraction); 
         tree->Branch("superCluster_sim_fraction","std::vector<std::vector<double> >",&superCluster_sim_fraction);
         tree->Branch("superCluster_recoToSim_fraction","std::vector<std::vector<double> >",&superCluster_recoToSim_fraction);
         tree->Branch("superCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&superCluster_recoToSim_fraction_sharedXtals);
         tree->Branch("superCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&superCluster_simEnergy_sharedXtals);
         tree->Branch("superCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&superCluster_recoEnergy_sharedXtals);  
         tree->Branch("superCluster_dR_genScore_MatchedIndex","std::vector<int>",&superCluster_dR_genScore_MatchedIndex);
         tree->Branch("superCluster_dR_simScore_MatchedIndex","std::vector<int>",&superCluster_dR_simScore_MatchedIndex);
         tree->Branch("superCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&superCluster_sim_nSharedXtals_MatchedIndex);
         tree->Branch("superCluster_sim_fraction_noHitsFraction_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_noHitsFraction_MatchedIndex); 
         tree->Branch("superCluster_sim_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_MatchedIndex);
         tree->Branch("superCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&superCluster_recoToSim_fraction_MatchedIndex);
         tree->Branch("superCluster_recoToSim_fraction_sharedXtals_MatchedIndex","std::vector<int>",&superCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
         tree->Branch("superCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&superCluster_simEnergy_sharedXtals_MatchedIndex);
         tree->Branch("superCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&superCluster_recoEnergy_sharedXtals_MatchedIndex);
      }    
      if(useRetunedSC_){   
         tree->Branch("retunedSuperCluster_rawEnergy","std::vector<float> ",&retunedSuperCluster_rawEnergy);
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
         tree->Branch("retunedSuperCluster_nXtals","std::vector<int>",&retunedSuperCluster_nXtals);  
         if(savePFCluster_) tree->Branch("retunedSuperCluster_seedIndex","std::vector<int>",&retunedSuperCluster_seedIndex);     
         if(savePFCluster_) tree->Branch("retunedSuperCluster_pfClustersIndex","std::vector<std::vector<int> >",&retunedSuperCluster_pfClustersIndex); 
         tree->Branch("retunedSuperCluster_psCluster_energy", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_energy);
         tree->Branch("retunedSuperCluster_psCluster_eta", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_eta);
         tree->Branch("retunedSuperCluster_psCluster_phi", "std::vector<std::vector<float> >", &retunedSuperCluster_psCluster_phi); 
         if(saveCaloParticles_){
            tree->Branch("retunedSuperCluster_dR_genScore","std::vector<std::vector<double> >",&retunedSuperCluster_dR_genScore);
            tree->Branch("retunedSuperCluster_dR_simScore","std::vector<std::vector<double> >",&retunedSuperCluster_dR_simScore);
            tree->Branch("retunedSuperCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_sim_nSharedXtals);
            tree->Branch("retunedSuperCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction_noHitsFraction); 
            tree->Branch("retunedSuperCluster_sim_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_sim_fraction);
            tree->Branch("retunedSuperCluster_recoToSim_fraction","std::vector<std::vector<double> >",&retunedSuperCluster_recoToSim_fraction);
            tree->Branch("retunedSuperCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_recoToSim_fraction_sharedXtals);
            tree->Branch("retunedSuperCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_simEnergy_sharedXtals);
            tree->Branch("retunedSuperCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&retunedSuperCluster_recoEnergy_sharedXtals);  
            tree->Branch("retunedSuperCluster_dR_genScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_dR_genScore_MatchedIndex);
            tree->Branch("retunedSuperCluster_dR_simScore_MatchedIndex","std::vector<int>",&retunedSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("retunedSuperCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_nSharedXtals_MatchedIndex);
            tree->Branch("retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<int>", &retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
            tree->Branch("retunedSuperCluster_sim_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_sim_fraction_MatchedIndex);
            tree->Branch("retunedSuperCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&retunedSuperCluster_recoToSim_fraction_MatchedIndex);
            tree->Branch("retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<int>", &retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
            tree->Branch("retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex);
            tree->Branch("retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex);
         } 
      } 
      if(useDeepSC_){
         tree->Branch("deepSuperCluster_rawEnergy","std::vector<float> ",&deepSuperCluster_rawEnergy);
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
         tree->Branch("deepSuperCluster_nXtals","std::vector<int>",&deepSuperCluster_nXtals);    
         if(savePFCluster_) tree->Branch("deepSuperCluster_seedIndex","std::vector<int>",&deepSuperCluster_seedIndex);     
         if(savePFCluster_) tree->Branch("deepSuperCluster_pfClustersIndex","std::vector<std::vector<int> >",&deepSuperCluster_pfClustersIndex); 
         tree->Branch("deepSuperCluster_psCluster_energy", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_energy);
         tree->Branch("deepSuperCluster_psCluster_eta", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_eta);
         tree->Branch("deepSuperCluster_psCluster_phi", "std::vector<std::vector<float> >", &deepSuperCluster_psCluster_phi);
         if(saveCaloParticles_){
            tree->Branch("deepSuperCluster_dR_genScore","std::vector<std::vector<double> >",&deepSuperCluster_dR_genScore);
            tree->Branch("deepSuperCluster_dR_simScore","std::vector<std::vector<double> >",&deepSuperCluster_dR_simScore);
            tree->Branch("deepSuperCluster_sim_nSharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_sim_nSharedXtals);
            tree->Branch("deepSuperCluster_sim_fraction_noHitsFraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction_noHitsFraction); 
            tree->Branch("deepSuperCluster_sim_fraction","std::vector<std::vector<double> >",&deepSuperCluster_sim_fraction);
            tree->Branch("deepSuperCluster_recoToSim_fraction","std::vector<std::vector<double> >",&deepSuperCluster_recoToSim_fraction);
            tree->Branch("deepSuperCluster_recoToSim_fraction_sharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_recoToSim_fraction_sharedXtals);
            tree->Branch("deepSuperCluster_simEnergy_sharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_simEnergy_sharedXtals);
            tree->Branch("deepSuperCluster_recoEnergy_sharedXtals","std::vector<std::vector<double> >",&deepSuperCluster_recoEnergy_sharedXtals);  
            tree->Branch("deepSuperCluster_dR_genScore_MatchedIndex","std::vector<int>",&deepSuperCluster_dR_genScore_MatchedIndex);
            tree->Branch("deepSuperCluster_dR_simScore_MatchedIndex","std::vector<int>",&deepSuperCluster_dR_simScore_MatchedIndex);
            tree->Branch("deepSuperCluster_sim_nSharedXtals_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_nSharedXtals_MatchedIndex);
            tree->Branch("deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex", "std::vector<int>", &deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex); 
            tree->Branch("deepSuperCluster_sim_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_sim_fraction_MatchedIndex);
            tree->Branch("deepSuperCluster_recoToSim_fraction_MatchedIndex","std::vector<int>",&deepSuperCluster_recoToSim_fraction_MatchedIndex);
            tree->Branch("deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex", "std::vector<int>", &deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);
            tree->Branch("deepSuperCluster_simEnergy_sharedXtals_MatchedIndex","std::vector<int>",&deepSuperCluster_simEnergy_sharedXtals_MatchedIndex);
            tree->Branch("deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex","std::vector<int>",&deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex);
         }   
      }   
   }
   if(savePFCluster_ && saveShowerShapes_){  
      tree->Branch("pfCluster_etaWidth","std::vector<float>",&pfCluster_etaWidth); 
      tree->Branch("pfCluster_phiWidth","std::vector<float>",&pfCluster_phiWidth); 
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
         if(useHcalTowers_ ){
            tree->Branch("deepSuperCluster_HoEraw","std::vector<float>",&deepSuperCluster_HoEraw);
            tree->Branch("deepSuperCluster_HoErawBC","std::vector<float>",&deepSuperCluster_HoErawBC);
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
   if(saveRechits_) {  
      ev.getByToken(ebRechitToken_, recHitsEB);
      if (!recHitsEB.isValid()) {
          std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
          return;
      }
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   if(saveRechits_) {
      ev.getByToken(eeRechitToken_, recHitsEE);
      if (!recHitsEE.isValid()) {
          std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFRecHit> > pfRecHits;
   if(savePFRechits_) {
      ev.getByToken(pfRecHitToken_, pfRecHits);
      if (!pfRecHits.isValid()) {
          std::cerr << "Analyze --> pfRecHits not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::PFCluster> > pfClusters;
   if(savePFCluster_) {
      ev.getByToken(pfClusterToken_, pfClusters);
      if (!pfClusters.isValid()) {
          std::cerr << "Analyze --> pfClusters not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   if(saveSuperCluster_) {
      ev.getByToken(ebSuperClusterToken_, superClusterEB);
      if (!superClusterEB.isValid()) {
          std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   if(saveSuperCluster_) {
      ev.getByToken(eeSuperClusterToken_, superClusterEE);
      if (!superClusterEE.isValid()) {
          std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
          return;
      }
   }

   edm::Handle<std::vector<reco::SuperCluster> > retunedSuperClusterEB;
   if(saveSuperCluster_ && useRetunedSC_) {
      ev.getByToken(ebRetunedSuperClusterToken_, retunedSuperClusterEB);
      if (!retunedSuperClusterEB.isValid()) {
          std::cerr << "Analyze --> retunedSuperClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > retunedSuperClusterEE;
   if(saveSuperCluster_ && useRetunedSC_) {
      ev.getByToken(eeRetunedSuperClusterToken_, retunedSuperClusterEE);
      if (!retunedSuperClusterEE.isValid()) {
          std::cerr << "Analyze --> retunedSuperClusterEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterEB;
   if(saveSuperCluster_ && useDeepSC_) {
      ev.getByToken(ebDeepSuperClusterToken_, deepSuperClusterEB);
      if (!deepSuperClusterEB.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > deepSuperClusterEE;
   if(saveSuperCluster_ && useDeepSC_) {
      ev.getByToken(eeDeepSuperClusterToken_, deepSuperClusterEE);
      if (!deepSuperClusterEE.isValid()) {
          std::cerr << "Analyze --> deepSuperClusterEE not found" << std::endl; 
          return;
      }
   }

   //compute EgammaTowers;
   Handle<CaloTowerCollection> hcalTowers;
   if(useHcalTowers_){
      ev.getByToken(hcalTowersToken_, hcalTowers);
      if (!hcalTowers.isValid()) {
          std::cerr << "Analyze --> hcalTowers not found" << std::endl; 
          return;
      } 
      towerIso1_ = new EgammaTowerIsolation(0.15, 0., 0., 1, hcalTowers.product());
      towerIso2_ = new EgammaTowerIsolation(0.15, 0., 0., 2, hcalTowers.product());
      egammaHadTower_ = new EgammaHadTower(iSetup);
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
   
   hitsAndEnergies_CaloPart.clear();
   GlobalPoint caloParticle_position;

   std::vector<CaloParticle> caloParts;
   for(const auto& iCalo : *(caloParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++) 
           if(iCalo.pdgId()==genID_.at(id) || genID_.at(id)==0) isGoodParticle=true;

       if(!isGoodParticle) continue; 
    
       GlobalPoint caloParticle_position = calculateAndSetPositionActual(getHitsAndEnergiesCaloPart(&iCalo,-1.), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
       if(caloParticle_position == GlobalPoint(-999999., -999999., -999999.)){
          std::cout << "Invalid position for caloParticle, skipping caloParticle!" << std::endl;
          continue;
       }    

       hitsAndEnergies_CaloPart.push_back(*getHitsAndEnergiesCaloPart(&iCalo,-1.));
       caloParts.push_back(iCalo); 
   }

   int nCaloParticles = caloParts.size(); 
   //std::cout << "CaloParticles size  : " << nCaloParticles << std::endl;
  
   genParticle_pfCluster_dR_genScore_MatchedIndex.clear();
   genParticle_superCluster_dR_genScore_MatchedIndex.clear();
   genParticle_retunedSuperCluster_dR_genScore_MatchedIndex.clear();
   genParticle_deepSuperCluster_dR_genScore_MatchedIndex.clear();
   genParticle_pfCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_superCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_retunedSuperCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_deepSuperCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
  
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
   caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_pfCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex.clear();     
   caloParticle_superCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_superCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_superCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();       
   caloParticle_deepSuperCluster_dR_simScore_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   caloParticle_pfCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
   caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_pfCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
   caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);     
   caloParticle_superCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
   caloParticle_superCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_superCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
   caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
   caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
   caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);       
   caloParticle_deepSuperCluster_dR_simScore_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_sim_fraction_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles); 
   caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.resize(nCaloParticles);  
   
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
   pfCluster_rawEnergy.clear();
   pfCluster_energy.clear();
   pfCluster_rawPt.clear();
   pfCluster_pt.clear();
   pfCluster_eta.clear();
   pfCluster_phi.clear();
   pfCluster_ieta.clear();
   pfCluster_iphi.clear();
   pfCluster_iz.clear();
   pfCluster_superClustersIndex.clear();
   pfCluster_retunedSuperClustersIndex.clear();
   pfCluster_deepSuperClustersIndex.clear(); 
   pfCluster_etaWidth.clear();
   pfCluster_phiWidth.clear();
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
   pfCluster_sim_nSharedXtals_MatchedIndex.clear();
   pfCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   pfCluster_sim_fraction_MatchedIndex.clear(); 
   pfCluster_recoToSim_fraction_MatchedIndex.clear();
   pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   pfCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   pfCluster_recoEnergy_sharedXtals_MatchedIndex.clear();     
   pfCluster_dR_genScore.resize(nPFClusters);
   pfCluster_dR_simScore.resize(nPFClusters);
   pfCluster_sim_nSharedXtals.resize(nPFClusters);
   pfCluster_sim_fraction_noHitsFraction.resize(nPFClusters);
   pfCluster_sim_fraction.resize(nPFClusters);
   pfCluster_recoToSim_fraction.resize(nPFClusters);
   pfCluster_recoToSim_fraction_sharedXtals.resize(nPFClusters);
   pfCluster_simEnergy_sharedXtals.resize(nPFClusters);
   pfCluster_recoEnergy_sharedXtals.resize(nPFClusters);   
   pfCluster_superClustersIndex.resize(nPFClusters); 
   pfCluster_retunedSuperClustersIndex.resize(nPFClusters);  
   pfCluster_deepSuperClustersIndex.resize(nPFClusters); 

   pfClusterHit_fraction.clear();
   pfClusterHit_rechitEnergy.clear();
   pfClusterHit_eta.clear();
   pfClusterHit_phi.clear();   
   pfClusterHit_ieta.clear();
   pfClusterHit_iphi.clear(); 
   pfClusterHit_iz.clear();           
   pfClusterHit_fraction.resize(nPFClusters);  
   pfClusterHit_rechitEnergy.resize(nPFClusters);  
   pfClusterHit_eta.resize(nPFClusters);  
   pfClusterHit_phi.resize(nPFClusters);     
   pfClusterHit_ieta.resize(nPFClusters);  
   pfClusterHit_iphi.resize(nPFClusters);   
   pfClusterHit_iz.resize(nPFClusters);   
  
   superCluster_rawEnergy.clear(); 
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
   superCluster_HoEraw.clear(); 
   superCluster_HoErawBC.clear();   
   superCluster_psCluster_energy.clear();
   superCluster_psCluster_eta.clear();
   superCluster_psCluster_phi.clear();
   superCluster_nXtals.clear();
   superCluster_dR_genScore_MatchedIndex.clear();
   superCluster_dR_simScore_MatchedIndex.clear();
   superCluster_sim_nSharedXtals_MatchedIndex.clear();
   superCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   superCluster_sim_fraction_MatchedIndex.clear(); 
   superCluster_recoToSim_fraction_MatchedIndex.clear();
   superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   superCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   superCluster_recoEnergy_sharedXtals_MatchedIndex.clear();    
   if(saveSuperCluster_){
      int nSuperClusters = (superClusterEB.product())->size() + (superClusterEE.product())->size();
      superCluster_seedIndex.resize(nSuperClusters);     
      superCluster_pfClustersIndex.resize(nSuperClusters);
      superCluster_psCluster_energy.resize((int)(superClusterEE.product())->size());
      superCluster_psCluster_eta.resize((int)(superClusterEE.product())->size());
      superCluster_psCluster_phi.resize((int)(superClusterEE.product())->size());
      superCluster_dR_genScore.resize(nSuperClusters);
      superCluster_dR_simScore.resize(nSuperClusters);
      superCluster_sim_nSharedXtals.resize(nSuperClusters);
      superCluster_sim_fraction_noHitsFraction.resize(nSuperClusters);
      superCluster_sim_fraction.resize(nSuperClusters);
      superCluster_recoToSim_fraction.resize(nSuperClusters);
      superCluster_recoToSim_fraction_sharedXtals.resize(nSuperClusters);
      superCluster_simEnergy_sharedXtals.resize(nSuperClusters);
      superCluster_recoEnergy_sharedXtals.resize(nSuperClusters);  
   }
   
   retunedSuperCluster_rawEnergy.clear(); 
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
   retunedSuperCluster_HoEraw.clear(); 
   retunedSuperCluster_HoErawBC.clear();    
   retunedSuperCluster_psCluster_energy.clear();
   retunedSuperCluster_psCluster_eta.clear();
   retunedSuperCluster_psCluster_phi.clear();
   retunedSuperCluster_nXtals.clear(); 
   retunedSuperCluster_dR_genScore_MatchedIndex.clear();
   retunedSuperCluster_dR_simScore_MatchedIndex.clear();
   retunedSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   retunedSuperCluster_sim_fraction_MatchedIndex.clear(); 
   retunedSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   if(useRetunedSC_ && saveSuperCluster_){
      int nRetunedSuperClusters = (retunedSuperClusterEB.product())->size() + (retunedSuperClusterEE.product())->size();
      retunedSuperCluster_seedIndex.resize(nRetunedSuperClusters);     
      retunedSuperCluster_pfClustersIndex.resize(nRetunedSuperClusters);
      retunedSuperCluster_psCluster_energy.resize((int)(retunedSuperClusterEE.product())->size());
      retunedSuperCluster_psCluster_eta.resize((int)(retunedSuperClusterEE.product())->size());
      retunedSuperCluster_psCluster_phi.resize((int)(retunedSuperClusterEE.product())->size());
      retunedSuperCluster_dR_genScore.resize(nRetunedSuperClusters);
      retunedSuperCluster_dR_simScore.resize(nRetunedSuperClusters);
      retunedSuperCluster_sim_nSharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_sim_fraction_noHitsFraction.resize(nRetunedSuperClusters);
      retunedSuperCluster_sim_fraction.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoToSim_fraction.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoToSim_fraction_sharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_simEnergy_sharedXtals.resize(nRetunedSuperClusters);
      retunedSuperCluster_recoEnergy_sharedXtals.resize(nRetunedSuperClusters);  
   }

   deepSuperCluster_rawEnergy.clear(); 
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
   deepSuperCluster_HoEraw.clear(); 
   deepSuperCluster_HoErawBC.clear();      
   deepSuperCluster_psCluster_energy.clear();
   deepSuperCluster_psCluster_eta.clear();
   deepSuperCluster_psCluster_phi.clear();
   deepSuperCluster_nXtals.clear();
   deepSuperCluster_dR_genScore_MatchedIndex.clear();
   deepSuperCluster_dR_simScore_MatchedIndex.clear();
   deepSuperCluster_sim_nSharedXtals_MatchedIndex.clear();
   deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.clear();
   deepSuperCluster_sim_fraction_MatchedIndex.clear(); 
   deepSuperCluster_recoToSim_fraction_MatchedIndex.clear();
   deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.clear();  
   deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.clear(); 
   deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.clear();  
   if(useDeepSC_){ 
      int nDeepSuperClusters = (deepSuperClusterEB.product())->size() + (deepSuperClusterEE.product())->size();
      deepSuperCluster_seedIndex.resize(nDeepSuperClusters); 
      deepSuperCluster_pfClustersIndex.resize(nDeepSuperClusters);
      deepSuperCluster_psCluster_energy.resize((int)(deepSuperClusterEE.product())->size());
      deepSuperCluster_psCluster_eta.resize((int)(deepSuperClusterEE.product())->size());
      deepSuperCluster_psCluster_phi.resize((int)(deepSuperClusterEE.product())->size());
      deepSuperCluster_dR_genScore.resize(nDeepSuperClusters);
      deepSuperCluster_dR_simScore.resize(nDeepSuperClusters);
      deepSuperCluster_sim_nSharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_sim_fraction_noHitsFraction.resize(nDeepSuperClusters);
      deepSuperCluster_sim_fraction.resize(nDeepSuperClusters);
      deepSuperCluster_recoToSim_fraction.resize(nDeepSuperClusters);
      deepSuperCluster_recoToSim_fraction_sharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_simEnergy_sharedXtals.resize(nDeepSuperClusters);
      deepSuperCluster_recoEnergy_sharedXtals.resize(nDeepSuperClusters);   
   }

   hitsAndEnergies_PFCluster.clear();
   hitsAndEnergies_SuperClusterEB.clear();
   hitsAndEnergies_SuperClusterEE.clear();
   hitsAndEnergies_RetunedSuperClusterEB.clear();
   hitsAndEnergies_RetunedSuperClusterEE.clear();
   hitsAndEnergies_DeepSuperClusterEB.clear();
   hitsAndEnergies_DeepSuperClusterEE.clear();

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
   if(savePFCluster_){
      for(const auto& iPFCluster : *(pfClusters.product())){  
          reco::CaloCluster caloBC(iPFCluster);
          hitsAndEnergies_PFCluster.push_back(*getHitsAndEnergiesBC(&caloBC,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
      }
   }
        
   if(saveSuperCluster_){
     for(const auto& iSuperCluster : *(superClusterEB.product())) 
         hitsAndEnergies_SuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
     for(const auto& iSuperCluster : *(superClusterEE.product())) 
         hitsAndEnergies_SuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
   }

   if(saveSuperCluster_ && useRetunedSC_){
     for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEB.product())) 
         hitsAndEnergies_RetunedSuperClusterEB.push_back(*getHitsAndEnergiesSC(&iRetunedSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
     for(const auto& iRetunedSuperCluster : *(retunedSuperClusterEE.product())) 
         hitsAndEnergies_RetunedSuperClusterEE.push_back(*getHitsAndEnergiesSC(&iRetunedSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
   }

   if(saveSuperCluster_ && useDeepSC_){
     for(const auto& iSuperCluster : *(deepSuperClusterEB.product())) 
         hitsAndEnergies_DeepSuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
     for(const auto& iSuperCluster : *(deepSuperClusterEE.product())) 
         hitsAndEnergies_DeepSuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,&(*(recHitsEB.product())), &(*(recHitsEE.product()))));
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear(); 

          pfCluster_rawEnergy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_energy.push_back(reduceFloat(iPFCluster.correctedEnergy(),nBits_));
          pfCluster_rawPt.push_back(reduceFloat(iPFCluster.energy()/TMath::CosH(iPFCluster.eta()),nBits_));
          pfCluster_pt.push_back(reduceFloat(iPFCluster.correctedEnergy()/TMath::CosH(iPFCluster.eta()),nBits_));
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
             widths_ = calculateCovariances(&iPFCluster, &(*(recHitsEB.product())), &(*_ebGeom));
             pfCluster_etaWidth.push_back(reduceFloat(widths_.first,nBits_));
             pfCluster_phiWidth.push_back(reduceFloat(widths_.second,nBits_));
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
             widths_ = calculateCovariances(&iPFCluster, &(*(recHitsEE.product())), &(*_eeGeom));
             pfCluster_etaWidth.push_back(reduceFloat(widths_.first,nBits_));
             pfCluster_phiWidth.push_back(reduceFloat(widths_.second,nBits_)); 
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
             const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster.hitsAndFractions();  
             for(unsigned int i = 0; i < hitsAndEnergies_PFCluster.at(iPFCl).size(); i++){      
                 cell = geometry->getPosition(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);
                 for(unsigned int hits=0; hits<hitsAndFractions.size(); hits++){
                      if(hitsAndFractions.at(hits).first.rawId() == hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.rawId())
                         pfClusterHit_fraction[iPFCl].push_back(hitsAndFractions.at(hits).second);
                 }
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
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             pfCluster_dR_genScore[iPFCl] = dR_genScore;        
             pfCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_dR_genScore, 999., false, 0., iPFCl));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iPFCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             pfCluster_nXtals.push_back((iPFCluster.hitsAndFractions()).size());   
             pfCluster_dR_simScore[iPFCl] = dR_simScore;  
             pfCluster_sim_nSharedXtals[iPFCl] = sim_nSharedXtals;   
             pfCluster_sim_fraction_noHitsFraction[iPFCl] = sim_fraction_noHitsFraction;    
             pfCluster_sim_fraction[iPFCl] = sim_fraction;   
             pfCluster_recoToSim_fraction[iPFCl] = recoToSim_fraction;   
             pfCluster_recoToSim_fraction_sharedXtals[iPFCl] = recoToSim_fraction_sharedXtals;        
             pfCluster_simEnergy_sharedXtals[iPFCl] = simEnergy_sharedXtals;   
             pfCluster_recoEnergy_sharedXtals[iPFCl] = recoEnergy_sharedXtals;        

             pfCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&pfCluster_dR_simScore, 999., false, 0., iPFCl));
             pfCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_nSharedXtals, -999., true, 0., iPFCl));
             pfCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction_noHitsFraction, -999., true, 0., iPFCl));
             pfCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_sim_fraction, -999., true, 0., iPFCl));
             pfCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&pfCluster_recoToSim_fraction, 999., false, 1., iPFCl)); 
             pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iPFCl)); 
             pfCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_simEnergy_sharedXtals, -999., true, 0., iPFCl)); 
             pfCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&pfCluster_recoEnergy_sharedXtals, -999., true, 0., iPFCl));
             
          }    
          iPFCl++;        
      }

      //save inverse of matchings
      if(saveGenParticles_){ 
         fillParticleMatchedIndex(&genParticle_pfCluster_dR_genScore_MatchedIndex,&pfCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_){ 
         fillParticleMatchedIndex(&caloParticle_pfCluster_dR_simScore_MatchedIndex,&pfCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_nSharedXtals_MatchedIndex,&pfCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_noHitsFraction_MatchedIndex,&pfCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_sim_fraction_MatchedIndex,&pfCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_recoToSim_fraction_MatchedIndex,&pfCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&pfCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_simEnergy_sharedXtals_MatchedIndex,&pfCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_pfCluster_recoEnergy_sharedXtals_MatchedIndex,&pfCluster_recoEnergy_sharedXtals_MatchedIndex);  
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear(); 
 
          superCluster_rawEnergy.push_back(reduceFloat(iSuperCluster.rawEnergy(),nBits_));
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
             superCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iSuperCluster.rawEnergy(),nBits_));
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
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iSuperCluster, towerIso1_, towerIso2_, egammaHadTower_, hcalTowers.product());
             superCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             superCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             superCluster_dR_genScore[iSC] = dR_genScore;        
             superCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_genScore, 999., false, 0., iSC));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             superCluster_nXtals.push_back((iSuperCluster.hitsAndFractions()).size());   
             superCluster_dR_simScore[iSC] = dR_simScore;  
             superCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             superCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             superCluster_sim_fraction[iSC] = sim_fraction;   
             superCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             superCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             superCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             superCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             superCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_simScore, 999., false, 0., iSC));
             superCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_nSharedXtals, -999., true, 0., iSC));
             superCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             superCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction, -999., true, 0., iSC));
             superCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             superCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             superCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear();    
          iSC_tmp++;
        
          superCluster_rawEnergy.push_back(reduceFloat(iSuperCluster.rawEnergy(),nBits_));
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
             superCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iSuperCluster.rawEnergy(),nBits_));
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
             superCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iSuperCluster.rawEnergy(),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iSuperCluster, towerIso1_, towerIso2_, egammaHadTower_, hcalTowers.product());
             superCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             superCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             superCluster_dR_genScore[iSC] = dR_genScore;        
             superCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_genScore, 999., false, 0., iSC));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             superCluster_nXtals.push_back((iSuperCluster.hitsAndFractions()).size());   
             superCluster_dR_simScore[iSC] = dR_simScore;  
             superCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             superCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             superCluster_sim_fraction[iSC] = sim_fraction;   
             superCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             superCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             superCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             superCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             superCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&superCluster_dR_simScore, 999., false, 0., iSC));
             superCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_nSharedXtals, -999., true, 0., iSC));
             superCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             superCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_sim_fraction, -999., true, 0., iSC));
             superCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             superCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             superCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             superCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&superCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
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

      //save pfCluster_superClustersIndex
      if(savePFCluster_ && saveSuperCluster_){
         for(unsigned int iSC=0; iSC<superCluster_pfClustersIndex.size(); iSC++)
             for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex.at(iSC).size(); iPF++)
                 if(superCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_superClustersIndex[superCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
      }

     //save inverse of matchings
     if(saveGenParticles_){ 
        fillParticleMatchedIndex(&genParticle_superCluster_dR_genScore_MatchedIndex,&superCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_){ 
         fillParticleMatchedIndex(&caloParticle_superCluster_dR_simScore_MatchedIndex,&superCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_nSharedXtals_MatchedIndex,&superCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_noHitsFraction_MatchedIndex,&superCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_sim_fraction_MatchedIndex,&superCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_recoToSim_fraction_MatchedIndex,&superCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&superCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_simEnergy_sharedXtals_MatchedIndex,&superCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_superCluster_recoEnergy_sharedXtals_MatchedIndex,&superCluster_recoEnergy_sharedXtals_MatchedIndex);  
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear(); 

          retunedSuperCluster_rawEnergy.push_back(reduceFloat(iRetunedSuperCluster.rawEnergy(),nBits_));
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
             retunedSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iRetunedSuperCluster.rawEnergy(),nBits_));
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
             retunedSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iRetunedSuperCluster.rawEnergy(),nBits_));
             retunedSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iRetunedSuperCluster, towerIso1_, towerIso2_, egammaHadTower_, hcalTowers.product());
             retunedSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             retunedSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             retunedSuperCluster_dR_genScore[iSC] = dR_genScore;        
             retunedSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             retunedSuperCluster_nXtals.push_back((iRetunedSuperCluster.hitsAndFractions()).size());   
             retunedSuperCluster_dR_simScore[iSC] = dR_simScore;  
             retunedSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             retunedSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             retunedSuperCluster_sim_fraction[iSC] = sim_fraction;   
             retunedSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             retunedSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             retunedSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             retunedSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             retunedSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_simScore, 999., false, 0., iSC));
             retunedSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction, -999., true, 0., iSC));
             retunedSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC));                
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear(); 
          iSC_tmp++;
        
          retunedSuperCluster_rawEnergy.push_back(reduceFloat(iRetunedSuperCluster.rawEnergy(),nBits_));
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
             retunedSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iRetunedSuperCluster.rawEnergy(),nBits_));
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
             retunedSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iRetunedSuperCluster.rawEnergy(),nBits_));
             retunedSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             retunedSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iRetunedSuperCluster, towerIso1_, towerIso2_, egammaHadTower_, hcalTowers.product());
             retunedSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             retunedSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             retunedSuperCluster_dR_genScore[iSC] = dR_genScore;        
             retunedSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iRetunedSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iRetunedSuperCluster.eta(),iRetunedSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             retunedSuperCluster_nXtals.push_back((iRetunedSuperCluster.hitsAndFractions()).size());   
             retunedSuperCluster_dR_simScore[iSC] = dR_simScore;  
             retunedSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             retunedSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             retunedSuperCluster_sim_fraction[iSC] = sim_fraction;   
             retunedSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             retunedSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             retunedSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             retunedSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             retunedSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_dR_simScore, 999., false, 0., iSC));
             retunedSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             retunedSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_sim_fraction, -999., true, 0., iSC));
             retunedSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&retunedSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
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

      //save pfCluster_retunedSuperClustersIndex
      if(savePFCluster_ && saveSuperCluster_ && useDeepSC_){
         for(unsigned int iSC=0; iSC<retunedSuperCluster_pfClustersIndex.size(); iSC++)
             for(unsigned int iPF=0; iPF<retunedSuperCluster_pfClustersIndex.at(iSC).size(); iPF++)
                 if(retunedSuperCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_retunedSuperClustersIndex[retunedSuperCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
      } 
      
      //save inverse of matchings
      if(saveGenParticles_){ 
         fillParticleMatchedIndex(&genParticle_retunedSuperCluster_dR_genScore_MatchedIndex,&retunedSuperCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_){ 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_dR_simScore_MatchedIndex,&retunedSuperCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_nSharedXtals_MatchedIndex,&retunedSuperCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex,&retunedSuperCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_sim_fraction_MatchedIndex,&retunedSuperCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_recoToSim_fraction_MatchedIndex,&retunedSuperCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&retunedSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex,&retunedSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex,&retunedSuperCluster_recoEnergy_sharedXtals_MatchedIndex);  
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear(); 

          deepSuperCluster_rawEnergy.push_back(reduceFloat(iDeepSuperCluster.rawEnergy(),nBits_));
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
             deepSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iDeepSuperCluster.rawEnergy(),nBits_));
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
             deepSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iDeepSuperCluster.rawEnergy(),nBits_));
             deepSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_)); 

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperCluster, towerIso1_, towerIso2_, egammaHadTower_, hcalTowers.product());
             deepSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_)); 
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             deepSuperCluster_dR_genScore[iSC] = dR_genScore;        
             deepSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             deepSuperCluster_nXtals.push_back((iDeepSuperCluster.hitsAndFractions()).size());   
             deepSuperCluster_dR_simScore[iSC] = dR_simScore;  
             deepSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             deepSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             deepSuperCluster_sim_fraction[iSC] = sim_fraction;   
             deepSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             deepSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             deepSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             deepSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             deepSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_simScore, 999., false, 0., iSC));
             deepSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction, -999., true, 0., iSC));
             deepSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC));         
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
          sim_nSharedXtals.clear();
          sim_fraction_noHitsFraction.clear();
          sim_fraction.clear();
          recoToSim_fraction.clear();
          recoToSim_fraction_sharedXtals.clear();  
          simEnergy_sharedXtals.clear(); 
          recoEnergy_sharedXtals.clear(); 
          iSC_tmp++;
        
          deepSuperCluster_rawEnergy.push_back(reduceFloat(iDeepSuperCluster.rawEnergy(),nBits_));     
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
             deepSuperCluster_r9.push_back(reduceFloat(showerShapes_[2]*showerShapes_[0]/iDeepSuperCluster.rawEnergy(),nBits_));
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
             deepSuperCluster_full5x5_r9.push_back(reduceFloat(showerShapes_[21]*showerShapes_[19]/iDeepSuperCluster.rawEnergy(),nBits_));
             deepSuperCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(showerShapes_[35],nBits_)); 
             deepSuperCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(showerShapes_[36],nBits_)); 
             deepSuperCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat(showerShapes_[37],nBits_));

             HoEs_.clear();
             HoEs_ = getHoE(&iDeepSuperCluster, towerIso1_, towerIso2_, egammaHadTower_, hcalTowers.product());
             deepSuperCluster_HoEraw.push_back(reduceFloat(HoEs_[0],nBits_)); 
             deepSuperCluster_HoErawBC.push_back(reduceFloat(HoEs_[1],nBits_));  
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_genScore.push_back(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_genScore.push_back(999.);     
             }    
             deepSuperCluster_dR_genScore[iSC] = dR_genScore;        
             deepSuperCluster_dR_genScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_genScore, 999., false, 0., iSC));
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<double> scores = getScores(&iDeepSuperCluster,&hitsAndEnergies_CaloPart.at(iCalo), &(*(recHitsEB.product())), &(*(recHitsEE.product())));
                 
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())<999.) dR_simScore.push_back(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iDeepSuperCluster.eta(),iDeepSuperCluster.phi())); 
                 else dR_simScore.push_back(999.);  

                 sim_nSharedXtals.push_back(scores[0]);  
                 sim_fraction_noHitsFraction.push_back(scores[1]);  
                 sim_fraction.push_back(scores[2]);  
                 recoToSim_fraction.push_back(scores[3]);  
                 recoToSim_fraction_sharedXtals.push_back(scores[4]);  
                 simEnergy_sharedXtals.push_back(scores[5]);  
                 recoEnergy_sharedXtals.push_back(scores[6]);  
             } 

             deepSuperCluster_nXtals.push_back((iDeepSuperCluster.hitsAndFractions()).size());   
             deepSuperCluster_dR_simScore[iSC] = dR_simScore;  
             deepSuperCluster_sim_nSharedXtals[iSC] = sim_nSharedXtals;   
             deepSuperCluster_sim_fraction_noHitsFraction[iSC] = sim_fraction_noHitsFraction;    
             deepSuperCluster_sim_fraction[iSC] = sim_fraction;   
             deepSuperCluster_recoToSim_fraction[iSC] = recoToSim_fraction;   
             deepSuperCluster_recoToSim_fraction_sharedXtals[iSC] = recoToSim_fraction_sharedXtals;        
             deepSuperCluster_simEnergy_sharedXtals[iSC] = simEnergy_sharedXtals;   
             deepSuperCluster_recoEnergy_sharedXtals[iSC] = recoEnergy_sharedXtals;        

             deepSuperCluster_dR_simScore_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_dR_simScore, 999., false, 0., iSC));
             deepSuperCluster_sim_nSharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_nSharedXtals, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction_noHitsFraction, -999., true, 0., iSC));
             deepSuperCluster_sim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_sim_fraction, -999., true, 0., iSC));
             deepSuperCluster_recoToSim_fraction_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction, 999., false, 1., iSC)); 
             deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoToSim_fraction_sharedXtals, 999., false, 1., iSC)); 
             deepSuperCluster_simEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_simEnergy_sharedXtals, -999., true, 0., iSC)); 
             deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex.push_back(getMatchedIndex(&deepSuperCluster_recoEnergy_sharedXtals, -999., true, 0., iSC)); 
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

      //save pfCluster_deepSuperClustersIndex
      if(savePFCluster_ && saveSuperCluster_ && useDeepSC_){
         for(unsigned int iSC=0; iSC<deepSuperCluster_pfClustersIndex.size(); iSC++)
             for(unsigned int iPF=0; iPF<deepSuperCluster_pfClustersIndex.at(iSC).size(); iPF++)
                 if(deepSuperCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_deepSuperClustersIndex[deepSuperCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);   
      } 

      //save inverse of matchings
      if(saveGenParticles_){ 
         fillParticleMatchedIndex(&genParticle_deepSuperCluster_dR_genScore_MatchedIndex,&deepSuperCluster_dR_genScore_MatchedIndex);
      } 
      if(saveCaloParticles_){ 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_dR_simScore_MatchedIndex,&deepSuperCluster_dR_simScore_MatchedIndex);
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_nSharedXtals_MatchedIndex,&deepSuperCluster_sim_nSharedXtals_MatchedIndex); 
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex,&deepSuperCluster_sim_fraction_noHitsFraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_sim_fraction_MatchedIndex,&deepSuperCluster_sim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_recoToSim_fraction_MatchedIndex,&deepSuperCluster_recoToSim_fraction_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex,&deepSuperCluster_recoToSim_fraction_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_simEnergy_sharedXtals_MatchedIndex,&deepSuperCluster_simEnergy_sharedXtals_MatchedIndex);  
         fillParticleMatchedIndex(&caloParticle_deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex,&deepSuperCluster_recoEnergy_sharedXtals_MatchedIndex);  
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
double RecoSimDumper::ptFast(const double energy, const math::XYZPoint& position, const math::XYZPoint& origin)
{
   const auto v = position - origin;
   return energy * std::sqrt(v.perp2() / v.mag2()); 
}

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesCaloPart(const CaloParticle* iCaloParticle, float simHitEnergy_cut)
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

std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster->hitsAndFractions();
    for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
        if(hitsAndFractions.at(i).first.subdetId()==EcalBarrel){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*recHitsEB->find(hitsAndFractions.at(i).first)).energy()));
        }else if(hitsAndFractions.at(i).first.subdetId()==EcalEndcap){
           HitsAndEnergies_tmp->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*recHitsEE->find(hitsAndFractions.at(i).first)).energy()));
        }
    }

    return HitsAndEnergies_tmp;
}


std::vector<std::pair<DetId, float> >* RecoSimDumper::getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_SuperCluster_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;

    for(reco::CaloCluster_iterator iBC = iSuperCluster->clustersBegin(); iBC != iSuperCluster->clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < seedrechits.size(); i++){  
            if(seedrechits.at(i).first.subdetId()==EcalBarrel){   
                if (HitsAndEnergies_map.find(seedrechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[seedrechits.at(i).first]=seedrechits.at(i).second * (*recHitsEB->find(seedrechits.at(i).first)).energy();    
                }else{
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*recHitsEB->find(seedrechits.at(i).first)).energy();
                } 
            }else if(seedrechits.at(i).first.subdetId()==EcalEndcap){   
                if (HitsAndEnergies_map.find(seedrechits.at(i).first) == HitsAndEnergies_map.end()) {
                    HitsAndEnergies_map[seedrechits.at(i).first]=seedrechits.at(i).second * (*recHitsEE->find(seedrechits.at(i).first)).energy();   
                }else{
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*recHitsEE->find(seedrechits.at(i).first)).energy();
                } 
            }
        }                      
    } 

    for(auto const& hit : HitsAndEnergies_map) 
        HitsAndEnergies_SuperCluster_tmp->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_SuperCluster_tmp;
}

std::pair<double,double> RecoSimDumper::calculateCovariances(const reco::PFCluster* pfCluster, const EcalRecHitCollection* recHits, const CaloSubdetectorGeometry* geometry) {

  double etaWidth = 0.;
  double phiWidth = 0.;
  double numeratorEtaWidth = 0;
  double numeratorPhiWidth = 0;

  double clEnergy = pfCluster->energy();
  double denominator = clEnergy;

  double clEta = pfCluster->position().eta();
  double clPhi = pfCluster->position().phi();

  const std::vector<std::pair<DetId, float> >& detId = pfCluster->hitsAndFractions();
  // Loop over recHits associated with the given SuperCluster
  for (std::vector<std::pair<DetId, float> >::const_iterator hit = detId.begin(); hit != detId.end(); ++hit) {
    EcalRecHitCollection::const_iterator rHit = recHits->find((*hit).first);
    //FIXME: THIS IS JUST A WORKAROUND A FIX SHOULD BE APPLIED
    if (rHit == recHits->end()) {
      continue;
    }
    auto this_cell = geometry->getGeometry(rHit->id());
    if (this_cell == nullptr) {
      //edm::LogInfo("SuperClusterShapeAlgo") << "pointer to the cell in Calculate_Covariances is NULL!";
      continue;
    }
    GlobalPoint position = this_cell->getPosition();
    //take into account energy fractions
    double energyHit = rHit->energy() * hit->second; 

    //form differences
    double dPhi = position.phi() - clPhi;
    if (dPhi > +Geom::pi()) {
      dPhi = Geom::twoPi() - dPhi;
    }
    if (dPhi < -Geom::pi()) {
      dPhi = Geom::twoPi() + dPhi;
    }

    double dEta = position.eta() - clEta;

    if (energyHit > 0) {
      numeratorEtaWidth += energyHit * dEta * dEta;
      numeratorPhiWidth += energyHit * dPhi * dPhi;
    }

    etaWidth = sqrt(numeratorEtaWidth / denominator);
    phiWidth = sqrt(numeratorPhiWidth / denominator);
  }

  return std::make_pair(etaWidth,phiWidth);
}

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
    shapes[17] = locCov_[1]; // sigmaIetaIphi
    shapes[18] = !edm::isFinite(locCov_[2]) ? 0. : sqrt(locCov_[2]); // sigmaIphiIphi 

    //if(std::isnan(shapes.at(17))) std::cout << shapes[0] << " - " << shapes[1] << " - " << shapes[2] << " - " << shapes[3] << " - " << shapes[4] << " - " << shapes[5]  << " - " << shapes[6] << " - " << shapes[7] << " - " << shapes[8] << " - " << shapes[9] << " - " << shapes[10] << " - " << shapes[11] << " - " << shapes[12] << " - " << shapes[13] << " - " << shapes[14] << " - " << shapes[15]  << " - " << shapes[16] << " - " << shapes[17] << " - " << shapes[18] << std::endl;

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
    shapes[36] = full5x5_locCov_[1]; // sigmaIetaIphi          
    shapes[37] = !edm::isFinite(full5x5_locCov_[2]) ? 0. : sqrt(full5x5_locCov_[2]); // sigmaIphiIphi

    for(unsigned iVar=0; iVar<shapes.size(); iVar++)
        if(std::isnan(shapes.at(iVar))) std::cout << "showerShape = " << iVar << " ---> NAN " << std::endl;  

    return shapes; 
}

std::vector<float> RecoSimDumper::getHoE(const reco::SuperCluster* iSuperCluster, EgammaTowerIsolation* towerIso1, EgammaTowerIsolation* towerIso2, const EgammaHadTower* egammaHadTower, const CaloTowerCollection* caloTower)
{
     std::vector<float> HoEs;
     HoEs.resize(2);
  
     std::vector<CaloTowerDetId> towersBehindCluster = egammaHadTower->towersOf(*iSuperCluster);
     double HoEraw1 = towerIso1->getTowerESum(iSuperCluster)/iSuperCluster->rawEnergy();
     double HoEraw2 = towerIso2->getTowerESum(iSuperCluster)/iSuperCluster->rawEnergy();        
     float HoEraw1bc = egammaHadTower->getDepth1HcalESum(towersBehindCluster, *caloTower)/iSuperCluster->energy();
     float HoEraw2bc = egammaHadTower->getDepth2HcalESum(towersBehindCluster, *caloTower)/iSuperCluster->energy(); 
     HoEs[0] = HoEraw1 + HoEraw2;
     HoEs[1] = HoEraw1bc + HoEraw2bc;
     
     return HoEs;
}

std::vector<double> RecoSimDumper::getScores(const reco::PFCluster* pfCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<double> scores;
    scores.resize(7);

    double nSharedXtals=0;
    double simFraction_noHitsFraction=0.;
    double simFraction=0.; 
    double recoToSim=0.; 
    double recoToSim_shared=0.; 

    double simEnergy=0.;
    double simEnergy_shared=0.;
    double simEnergy_shared_noHitsFraction=0.;
    double recoEnergy=pfCluster->energy();
    double recoEnergy_shared=0.;
   
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart)
        simEnergy+=hit_CaloPart.second;
   
    const std::vector<std::pair<DetId,float> >* hitsAndFractions = &pfCluster->hitsAndFractions();
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){
        for(const std::pair<DetId, float>& hit_Cluster : *hitsAndFractions){     
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
               nSharedXtals+=1.;
               simEnergy_shared+=hit_CaloPart.second*hit_Cluster.second;
               simEnergy_shared_noHitsFraction+=hit_CaloPart.second;
               if(hit_Cluster.first.subdetId()==EcalBarrel) recoEnergy_shared+=hit_Cluster.second*(*recHitsEB->find(hit_Cluster.first)).energy();
               if(hit_Cluster.first.subdetId()==EcalEndcap) recoEnergy_shared+=hit_Cluster.second*(*recHitsEE->find(hit_Cluster.first)).energy(); 
            } 
        }
    }

    if(nSharedXtals<=0.) nSharedXtals = -1.; 

    if(simEnergy_shared_noHitsFraction>0. && simEnergy>0.) simFraction_noHitsFraction = simEnergy_shared_noHitsFraction/simEnergy;
    else simFraction_noHitsFraction = -1.; 

    if(simEnergy_shared>0. && simEnergy>0.) simFraction = simEnergy_shared/simEnergy;
    else simFraction = -1.;
 
    if(recoEnergy>0. && simEnergy_shared>0.) recoToSim = recoEnergy/simEnergy_shared;
    else recoToSim = -1.; 

    if(recoEnergy_shared>0. && simEnergy_shared>0.) recoToSim_shared = recoEnergy_shared/simEnergy_shared;
    else recoToSim_shared = -1.;  

    if(simEnergy_shared<=0) simEnergy_shared = -1.;

    if(recoEnergy_shared<=0) recoEnergy_shared = -1.;
    
    scores[0] = (double)nSharedXtals;
    scores[1] = simFraction_noHitsFraction;
    scores[2] = simFraction;
    scores[3] = recoToSim;
    scores[4] = recoToSim_shared;
    scores[5] = simEnergy_shared;
    scores[6] = recoEnergy_shared;

    for(unsigned iVar=0; iVar<scores.size(); iVar++)
        if(std::isnan(scores.at(iVar))) std::cout << "score = " << iVar << " ---> NAN " << std::endl; 

    return scores;
}

std::vector<double> RecoSimDumper::getScores(const reco::SuperCluster* superCluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart, const EcalRecHitCollection* recHitsEB, const EcalRecHitCollection* recHitsEE)
{
    std::vector<double> scores;
    scores.resize(7);

    double nSharedXtals=0;
    double simFraction_noHitsFraction=0.;
    double simFraction=0.; 
    double recoToSim=0.; 
    double recoToSim_shared=0.; 

    double simEnergy=0.;
    double simEnergy_shared=0.;
    double simEnergy_shared_noHitsFraction=0.;
    double recoEnergy=superCluster->energy();
    double recoEnergy_shared=0.;
   
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart)
        simEnergy+=hit_CaloPart.second;
   
    std::vector<DetId> superCluster_IDs;
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){      
        for(reco::CaloCluster_iterator iBC = superCluster->clustersBegin(); iBC != superCluster->clustersEnd(); ++iBC){
            const std::vector<std::pair<DetId,float> >* hitsAndFractions = &( *iBC )->hitsAndFractions();
            for(const std::pair<DetId, float>& hit_Cluster : *hitsAndFractions){      
                if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){
                   simEnergy_shared+=hit_CaloPart.second*hit_Cluster.second;
                   simEnergy_shared_noHitsFraction+=hit_CaloPart.second;
                   if(hit_Cluster.first.subdetId()==EcalBarrel) recoEnergy_shared+=hit_Cluster.second*(*recHitsEB->find(hit_Cluster.first)).energy();
                   if(hit_Cluster.first.subdetId()==EcalEndcap) recoEnergy_shared+=hit_Cluster.second*(*recHitsEE->find(hit_Cluster.first)).energy();
                   if(std::find(superCluster_IDs.begin(),superCluster_IDs.end(),hit_Cluster.first) == superCluster_IDs.end()) superCluster_IDs.push_back(hit_Cluster.first);
                } 
            } 
        }
    }

    nSharedXtals = (int)superCluster_IDs.size();
    if(nSharedXtals<=0.) nSharedXtals = -1.; 

    if(simEnergy_shared_noHitsFraction>0. && simEnergy>0.) simFraction_noHitsFraction = simEnergy_shared_noHitsFraction/simEnergy;
    else simFraction_noHitsFraction = -1.; 

    if(simEnergy_shared>0. && simEnergy>0.) simFraction = simEnergy_shared/simEnergy;
    else simFraction = -1.;
 
    if(recoEnergy>0. && simEnergy_shared>0.) recoToSim = recoEnergy/simEnergy_shared;
    else recoToSim = -1.; 

    if(recoEnergy_shared>0. && simEnergy_shared>0.) recoToSim_shared = recoEnergy_shared/simEnergy_shared;
    else recoToSim_shared = -1.;  

    if(simEnergy_shared<=0) simEnergy_shared = -1.;

    if(recoEnergy_shared<=0) recoEnergy_shared = -1.;
    
    scores[0] = (double)nSharedXtals;
    scores[1] = simFraction_noHitsFraction;
    scores[2] = simFraction;
    scores[3] = recoToSim;
    scores[4] = recoToSim_shared;
    scores[5] = simEnergy_shared;
    scores[6] = recoEnergy_shared;

    for(unsigned iVar=0; iVar<scores.size(); iVar++)
        if(std::isnan(scores.at(iVar))) std::cout << "score = " << iVar << " ---> NAN " << std::endl; 

    return scores;
}

int RecoSimDumper::getMatchedIndex(std::vector<std::vector<double>>* scores, double selection, bool useMax, double scale, int iCl)
{
   int matchedIndex = -1; 

   std::vector<double> score;
   for(unsigned int iCalo=0; iCalo<scores->at(iCl).size(); iCalo++){
       if(scores->at(iCl).at(iCalo)>0.) score.push_back(fabs(scale-scores->at(iCl).at(iCalo))); 
       else score.push_back(-1.); 
   }

   if(!useMax){ 
      std::replace(score.begin(),score.end(), -1., 999.);
      if(std::all_of(score.begin(),score.end(),[](double i){return i==-1.;}) || std::all_of(score.begin(),score.end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::min_element(score.begin(),score.end()) - score.begin();
   }else{
      std::replace(score.begin(),score.end(), -1., -999.); 
      if(std::all_of(score.begin(),score.end(),[](double i){return i==-1.;}) || std::all_of(score.begin(),score.end(),[](double i){return i==-999.;})) matchedIndex=-1;
      else matchedIndex = std::max_element(score.begin(),score.end()) - score.begin();
   }

   if(matchedIndex==-1) return -1;
   
   if(useMax && score.at(matchedIndex) > selection) return matchedIndex;
   if(!useMax && score.at(matchedIndex) < selection) return matchedIndex; 
   
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
      //throw cms::Exception("InvalidLayer") << "ECAL Position Calc only accepts ECAL_BARREL or ECAL_ENDCAP";
      std::cout << "WARNING: wrong layer, ECAL Position Calc only accepts ECAL_BARREL or ECAL_ENDCAP, returning invalid position" << std::endl;
  }

  if (ecal_geom == nullptr)
     return GlobalPoint(-999999., -999999., -999999.);

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

