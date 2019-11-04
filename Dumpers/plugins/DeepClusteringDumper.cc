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
#include "RecoSimStudies/Dumpers/plugins/DeepClusteringDumper.h"

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
DeepClusteringDumper::DeepClusteringDumper(const edm::ParameterSet& iConfig)
{

   genToken_                = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_           = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   pfRecHitToken_           = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection")); 
   pfClusterToken_          = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   ebSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   
   doCompression_           = iConfig.getParameter<bool>("doCompression");
   nBits_                   = iConfig.getParameter<int>("nBits");
   saveGenParticles_        = iConfig.getParameter<bool>("saveGenParticles");
   saveCaloParticles_       = iConfig.getParameter<bool>("saveCaloParticles");
   saveSimhits_             = iConfig.getParameter<bool>("saveSimhits");
   saveRechits_             = iConfig.getParameter<bool>("saveRechits");
   savePFRechits_           = iConfig.getParameter<bool>("savePFRechits"); 
   savePFCluster_           = iConfig.getParameter<bool>("savePFCluster");
   savePFClusterhits_       = iConfig.getParameter<bool>("savePFClusterhits");
   saveSuperCluster_        = iConfig.getParameter<bool>("saveSuperCluster");
   saveShowerShapes_        = iConfig.getParameter<bool>("saveShowerShapes");
   saveScores_              = iConfig.getParameter<bool>("saveScores");
   genID_                   = iConfig.getParameter<std::vector<int>>("genID");
   
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   if(saveGenParticles_){
      tree->Branch("genParticle_id","std::vector<int>",&genParticle_id);
      tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
      if(savePFCluster_) tree->Branch("genParticle_pfCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_pfCluster_dR_genScore_MatchedIndex);
      if(saveSuperCluster_) tree->Branch("genParticle_superCluster_dR_genScore_MatchedIndex","std::vector<std::vector<int> >",&genParticle_superCluster_dR_genScore_MatchedIndex);
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
         tree->Branch("caloParticle_pfCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_n_shared_xtals_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_min1_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_min1_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_fraction_min3_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_fraction_min3_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_diff_MatchedIndex);
         tree->Branch("caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex);   
         tree->Branch("caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex); 
      }
      if(saveSuperCluster_){
         tree->Branch("caloParticle_superCluster_dR_simScore_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_dR_simScore_MatchedIndex);
         tree->Branch("caloParticle_superCluster_n_shared_xtals_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_n_shared_xtals_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_fraction_min1_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_min1_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_fraction_min3_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_fraction_min3_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_rechit_diff_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_rechit_diff_MatchedIndex);
         tree->Branch("caloParticle_superCluster_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_sim_rechit_fraction_MatchedIndex);   
         tree->Branch("caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<std::vector<int> >",&caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex); 
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
      if(saveSuperCluster_) tree->Branch("pfCluster_superClustersIndex","std::vector<std::vector<int> >",&pfCluster_superClustersIndex); 
      if(saveCaloParticles_){ 
         tree->Branch("pfCluster_dR_genScore_MatchedIndex","std::vector<int>",&pfCluster_dR_genScore_MatchedIndex);
         tree->Branch("pfCluster_dR_simScore_MatchedIndex","std::vector<int>",&pfCluster_dR_simScore_MatchedIndex);
         tree->Branch("pfCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&pfCluster_n_shared_xtals_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_min1_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_min1_MatchedIndex);
         tree->Branch("pfCluster_sim_fraction_min3_MatchedIndex","std::vector<int>",&pfCluster_sim_fraction_min3_MatchedIndex);
         tree->Branch("pfCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_diff_MatchedIndex);
         tree->Branch("pfCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&pfCluster_sim_rechit_fraction_MatchedIndex);   
         tree->Branch("pfCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&pfCluster_global_sim_rechit_fraction_MatchedIndex); 
      } 
      if(saveCaloParticles_ && saveScores_){
         tree->Branch("pfCluster_dR_genScore","std::vector<std::vector<float> >",&pfCluster_dR_genScore);
         tree->Branch("pfCluster_dR_simScore","std::vector<std::vector<float> >",&pfCluster_dR_simScore);
         tree->Branch("pfCluster_n_shared_xtals","std::vector<std::vector<int> >",&pfCluster_n_shared_xtals);
         tree->Branch("pfCluster_sim_fraction","std::vector<std::vector<float> >",&pfCluster_sim_fraction);
         tree->Branch("pfCluster_sim_fraction_min1","std::vector<std::vector<float> >",&pfCluster_sim_fraction_min1);
         tree->Branch("pfCluster_sim_fraction_min3","std::vector<std::vector<float> >",&pfCluster_sim_fraction_min3);
         tree->Branch("pfCluster_sim_rechit_diff","std::vector<std::vector<float> >",&pfCluster_sim_rechit_diff);
         tree->Branch("pfCluster_sim_rechit_fraction","std::vector<std::vector<float> >",&pfCluster_sim_rechit_fraction);   
         tree->Branch("pfCluster_global_sim_rechit_fraction","std::vector<std::vector<float> >",&pfCluster_global_sim_rechit_fraction); 
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
      tree->Branch("superCluster_ieta","std::vector<int>",&superCluster_ieta);
      tree->Branch("superCluster_iphi","std::vector<int>",&superCluster_iphi);  
      tree->Branch("superCluster_iz","std::vector<int>",&superCluster_iz);  
      if(savePFCluster_) tree->Branch("superCluster_pfClustersIndex","std::vector<std::vector<int> >",&superCluster_pfClustersIndex);   
      if(saveCaloParticles_){
         tree->Branch("superCluster_dR_genScore_MatchedIndex","std::vector<int>",&superCluster_dR_genScore_MatchedIndex);
         tree->Branch("superCluster_dR_simScore_MatchedIndex","std::vector<int>",&superCluster_dR_simScore_MatchedIndex);
         tree->Branch("superCluster_n_shared_xtals_MatchedIndex","std::vector<int>",&superCluster_n_shared_xtals_MatchedIndex);
         tree->Branch("superCluster_sim_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_MatchedIndex);
         tree->Branch("superCluster_sim_fraction_min1_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_min1_MatchedIndex);
         tree->Branch("superCluster_sim_fraction_min3_MatchedIndex","std::vector<int>",&superCluster_sim_fraction_min3_MatchedIndex);
         tree->Branch("superCluster_sim_rechit_diff_MatchedIndex","std::vector<int>",&superCluster_sim_rechit_diff_MatchedIndex);
         tree->Branch("superCluster_sim_rechit_fraction_MatchedIndex","std::vector<int>",&superCluster_sim_rechit_fraction_MatchedIndex);   
         tree->Branch("superCluster_global_sim_rechit_fraction_MatchedIndex","std::vector<int>",&superCluster_global_sim_rechit_fraction_MatchedIndex); 
      } 
      if(saveCaloParticles_ && saveScores_){
         tree->Branch("superCluster_dR_genScore","std::vector<std::vector<float> >",&superCluster_dR_genScore);
         tree->Branch("superCluster_dR_simScore","std::vector<std::vector<float> >",&superCluster_dR_simScore);
         tree->Branch("superCluster_n_shared_xtals","std::vector<std::vector<int> >",&superCluster_n_shared_xtals);
         tree->Branch("superCluster_sim_fraction","std::vector<std::vector<float> >",&superCluster_sim_fraction);
         tree->Branch("superCluster_sim_fraction_min1","std::vector<std::vector<float> >",&superCluster_sim_fraction_min1);
         tree->Branch("superCluster_sim_fraction_min3","std::vector<std::vector<float> >",&superCluster_sim_fraction_min3);
         tree->Branch("superCluster_sim_rechit_diff","std::vector<std::vector<float> >",&superCluster_sim_rechit_diff);
         tree->Branch("superCluster_sim_rechit_fraction","std::vector<std::vector<float> >",&superCluster_sim_rechit_fraction);   
         tree->Branch("superCluster_global_sim_rechit_fraction","std::vector<std::vector<float> >",&superCluster_global_sim_rechit_fraction); 
      }  
      if(saveShowerShapes_){ 
         tree->Branch("superCluster_r9","std::vector<float> ",&superCluster_r9);
         tree->Branch("superCluster_sigmaIetaIeta","std::vector<float> ",&superCluster_sigmaIetaIeta);
         tree->Branch("superCluster_sigmaIetaIphi","std::vector<float> ",&superCluster_sigmaIetaIphi);
         tree->Branch("superCluster_sigmaIphiIphi","std::vector<float> ",&superCluster_sigmaIphiIphi);
         tree->Branch("superCluster_full5x5_r9","std::vector<float> ",&superCluster_full5x5_r9);
         tree->Branch("superCluster_full5x5_sigmaIetaIeta","std::vector<float> ",&superCluster_full5x5_sigmaIetaIeta);
         tree->Branch("superCluster_full5x5_sigmaIetaIphi","std::vector<float> ",&superCluster_full5x5_sigmaIetaIphi);
         tree->Branch("superCluster_full5x5_sigmaIphiIphi","std::vector<float> ",&superCluster_full5x5_sigmaIphiIphi);
      }
   }
}

DeepClusteringDumper::~DeepClusteringDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void DeepClusteringDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
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

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();

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
  
   genParticle_id.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();
   genParticle_pfCluster_dR_genScore_MatchedIndex.clear();
   genParticle_superCluster_dR_genScore_MatchedIndex.clear();
   genParticle_pfCluster_dR_genScore_MatchedIndex.resize(nGenParticles);
   genParticle_superCluster_dR_genScore_MatchedIndex.resize(nGenParticles);

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
   caloParticle_pfCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_min1_MatchedIndex.clear();
   caloParticle_pfCluster_sim_fraction_min3_MatchedIndex.clear();
   caloParticle_pfCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_superCluster_dR_simScore_MatchedIndex.clear();
   caloParticle_superCluster_n_shared_xtals_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_min1_MatchedIndex.clear();
   caloParticle_superCluster_sim_fraction_min3_MatchedIndex.clear();
   caloParticle_superCluster_sim_rechit_diff_MatchedIndex.clear();
   caloParticle_superCluster_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex.clear();
   caloParticle_pfCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_min1_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_fraction_min3_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_dR_simScore_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_n_shared_xtals_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_min1_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_fraction_min3_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_rechit_diff_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex.resize(nCaloParticles);
   
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
   pfCluster_dR_genScore_MatchedIndex.clear();
   pfCluster_dR_simScore_MatchedIndex.clear();
   pfCluster_n_shared_xtals_MatchedIndex.clear();
   pfCluster_sim_fraction_MatchedIndex.clear();
   pfCluster_sim_fraction_min1_MatchedIndex.clear();
   pfCluster_sim_fraction_min3_MatchedIndex.clear();
   pfCluster_sim_rechit_diff_MatchedIndex.clear();
   pfCluster_sim_rechit_fraction_MatchedIndex.clear();
   pfCluster_global_sim_rechit_fraction_MatchedIndex.clear();  
   pfCluster_dR_genScore.clear();
   pfCluster_dR_simScore.clear();
   pfCluster_n_shared_xtals.clear();
   pfCluster_sim_fraction.clear();
   pfCluster_sim_fraction_min1.clear();
   pfCluster_sim_fraction_min3.clear();
   pfCluster_sim_rechit_diff.clear();
   pfCluster_sim_rechit_fraction.clear();
   pfCluster_global_sim_rechit_fraction.clear();   
   pfCluster_dR_genScore.resize(nPFClusters);
   pfCluster_dR_simScore.resize(nPFClusters);
   pfCluster_n_shared_xtals.resize(nPFClusters);
   pfCluster_sim_fraction.resize(nPFClusters);
   pfCluster_sim_fraction_min1.resize(nPFClusters);
   pfCluster_sim_fraction_min3.resize(nPFClusters);
   pfCluster_sim_rechit_diff.resize(nPFClusters);
   pfCluster_sim_rechit_fraction.resize(nPFClusters);
   pfCluster_global_sim_rechit_fraction.resize(nPFClusters); 
   pfCluster_superClustersIndex.resize(nPFClusters); 

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
   
   int nSuperClusters = (superClusterEB.product())->size() + (superClusterEE.product())->size();

   superCluster_energy.clear(); 
   superCluster_eta.clear(); 
   superCluster_phi.clear();  
   superCluster_ieta.clear(); 
   superCluster_iphi.clear();    
   superCluster_iz.clear(); 
   superCluster_pfClustersIndex.clear(); 
   superCluster_r9.clear(); 
   superCluster_sigmaIetaIeta.clear(); 
   superCluster_sigmaIetaIphi.clear(); 
   superCluster_sigmaIphiIphi.clear(); 
   superCluster_full5x5_r9.clear(); 
   superCluster_full5x5_sigmaIetaIeta.clear();
   superCluster_full5x5_sigmaIetaIphi.clear();
   superCluster_full5x5_sigmaIphiIphi.clear(); 
   superCluster_dR_genScore_MatchedIndex.clear();  
   superCluster_dR_simScore_MatchedIndex.clear();  
   superCluster_n_shared_xtals_MatchedIndex.clear();  
   superCluster_sim_fraction_MatchedIndex.clear();  
   superCluster_sim_fraction_min1_MatchedIndex.clear();  
   superCluster_sim_fraction_min3_MatchedIndex.clear();  
   superCluster_sim_rechit_diff_MatchedIndex.clear();  
   superCluster_sim_rechit_fraction_MatchedIndex.clear();  
   superCluster_global_sim_rechit_fraction_MatchedIndex.clear();   
   superCluster_dR_genScore.clear();
   superCluster_dR_simScore.clear();
   superCluster_n_shared_xtals.clear();
   superCluster_sim_fraction.clear();
   superCluster_sim_fraction_min1.clear();
   superCluster_sim_fraction_min3.clear();
   superCluster_sim_rechit_diff.clear();
   superCluster_sim_rechit_fraction.clear();
   superCluster_global_sim_rechit_fraction.clear();  
   superCluster_dR_genScore.resize(nSuperClusters);
   superCluster_dR_simScore.resize(nSuperClusters);
   superCluster_n_shared_xtals.resize(nSuperClusters);
   superCluster_sim_fraction.resize(nSuperClusters);
   superCluster_sim_fraction_min1.resize(nSuperClusters);
   superCluster_sim_fraction_min3.resize(nSuperClusters);
   superCluster_sim_rechit_diff.resize(nSuperClusters);
   superCluster_sim_rechit_fraction.resize(nSuperClusters);
   superCluster_global_sim_rechit_fraction.resize(nSuperClusters);  
   superCluster_pfClustersIndex.resize(nSuperClusters);
  
   hitsAndEnergies_CaloPart.clear();
   hitsAndEnergies_PFCluster.clear();
   hitsAndEnergies_SuperClusterEB.clear();
   hitsAndEnergies_SuperClusterEE.clear();

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

       caloParticle_simPt.push_back(reduceFloat(caloParts.at(iCalo).pt(),nBits_));
      
       hitsAndEnergies_CaloPart.push_back(*getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo))));
       GlobalPoint caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
       caloParticle_simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
       caloParticle_simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
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
     hits_SuperCluster.resize(nSuperClusters);
     for(const auto& iSuperCluster : *(superClusterEB.product())) 
         hitsAndEnergies_SuperClusterEB.push_back(*getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE));
     for(unsigned int iSC = 0; iSC < hitsAndEnergies_SuperClusterEB.size(); iSC++)
         for(unsigned int i = 0; i < hitsAndEnergies_SuperClusterEB.at(iSC).size(); i++)
             hits_SuperCluster[iSC].push_back(hitsAndEnergies_SuperClusterEB.at(iSC).at(i).first);   
     for(const auto& iSuperCluster : *(superClusterEE.product())) 
         hitsAndEnergies_SuperClusterEE.push_back(*getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE));
     for(unsigned int iSC = 0; iSC < hitsAndEnergies_SuperClusterEE.size(); iSC++)
         for(unsigned int i = 0; i < hitsAndEnergies_SuperClusterEE.at(iSC).size(); i++)
             hits_SuperCluster[iSC].push_back(hitsAndEnergies_SuperClusterEE.at(iSC).at(i).first);  
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
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_min1.clear();
          sim_fraction_min3.clear();
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();

          pfCluster_energy.push_back(reduceFloat(iPFCluster.energy(),nBits_));
          pfCluster_eta.push_back(reduceFloat(iPFCluster.eta(),nBits_));
          pfCluster_phi.push_back(reduceFloat(iPFCluster.phi(),nBits_));
          reco::CaloCluster caloBC(iPFCluster);

          math::XYZPoint caloPos = caloBC.position();
          if(std::abs(iPFCluster.eta()) < 1.479){  
             EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             pfCluster_ieta.push_back(eb_id.ieta());
             pfCluster_iphi.push_back(eb_id.iphi());
             pfCluster_iz.push_back(0); 
          }else{  
             int iz=-99;
             EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
             if(ee_id.zside()<0) iz=-1;
             if(ee_id.zside()>0) iz=1;     
             pfCluster_ieta.push_back(ee_id.ix());
             pfCluster_iphi.push_back(ee_id.iy());
             pfCluster_iz.push_back(iz); 
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
                    pfClusterHit_ieta[iPFCl].push_back(eb_id.ieta());
                    pfClusterHit_iphi[iPFCl].push_back(eb_id.iphi());
                    pfClusterHit_iz[iPFCl].push_back(0); 
                 }else if(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first.subdetId()==EcalEndcap){  
                    int iz=-99;
                    EEDetId ee_id(hitsAndEnergies_PFCluster.at(iPFCl).at(i).first);  
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
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_genScore.push_back(reduceFloat(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi()),nBits_)); 
                 else dR_genScore.push_back(reduceFloat(999.,nBits_));    
             }    
             if(saveScores_) pfCluster_dR_genScore[iPFCl] = dR_genScore;        
             if(std::equal(dR_genScore.begin() + 1, dR_genScore.end(), dR_genScore.begin())) pfCluster_dR_genScore_MatchedIndex.push_back(-1);
             else pfCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin());  
          } 
          if(saveCaloParticles_){ 
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<float> scores = getScores(&hitsAndEnergies_PFCluster.at(iPFCl),&hitsAndEnergies_CaloPart.at(iCalo));         
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_simScore.push_back(reduceFloat(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi()),nBits_)); 
                 else dR_simScore.push_back(reduceFloat(999.,nBits_)); 
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 sim_fraction_min1.push_back(scores[5]);  
                 sim_fraction_min3.push_back(scores[6]);                     
             } 
             if(saveScores_){
                pfCluster_dR_simScore[iPFCl] = dR_simScore;  
                pfCluster_n_shared_xtals[iPFCl] = n_shared_xtals;  
                pfCluster_sim_fraction[iPFCl] = sim_fraction;  
                pfCluster_sim_rechit_diff[iPFCl] = sim_rechit_diff; 
                pfCluster_sim_rechit_fraction[iPFCl] = sim_rechit_fraction;           
                pfCluster_global_sim_rechit_fraction[iPFCl] = global_sim_rechit_fraction;
                pfCluster_sim_fraction_min1[iPFCl] = sim_fraction_min1;  
                pfCluster_sim_fraction_min3[iPFCl] = sim_fraction_min3; 
             }
             if(std::equal(dR_simScore.begin() + 1, dR_simScore.end(), dR_simScore.begin())) pfCluster_dR_simScore_MatchedIndex.push_back(-1);
             else pfCluster_dR_simScore_MatchedIndex.push_back(std::min_element(dR_simScore.begin(),dR_simScore.end()) - dR_simScore.begin());  
             if(std::equal(n_shared_xtals.begin() + 1, n_shared_xtals.end(), n_shared_xtals.begin())) pfCluster_n_shared_xtals_MatchedIndex.push_back(-1);
             else pfCluster_n_shared_xtals_MatchedIndex.push_back(std::max_element(n_shared_xtals.begin(),n_shared_xtals.end()) - n_shared_xtals.begin());  
             if(std::equal(sim_fraction.begin() + 1, sim_fraction.end(), sim_fraction.begin())) pfCluster_sim_fraction_MatchedIndex.push_back(-1);
             else pfCluster_sim_fraction_MatchedIndex.push_back(std::max_element(sim_fraction.begin(),sim_fraction.end()) - sim_fraction.begin()); 
             if(std::equal(sim_fraction_min1.begin() + 1, sim_fraction_min1.end(), sim_fraction_min1.begin())) pfCluster_sim_fraction_min1_MatchedIndex.push_back(-1);
             else pfCluster_sim_fraction_min1_MatchedIndex.push_back(std::max_element(sim_fraction_min1.begin(),sim_fraction_min1.end()) - sim_fraction_min1.begin());
             if(std::equal(sim_fraction_min3.begin() + 1, sim_fraction_min3.end(), sim_fraction_min3.begin())) pfCluster_sim_fraction_min3_MatchedIndex.push_back(-1);
             else pfCluster_sim_fraction_min3_MatchedIndex.push_back(std::max_element(sim_fraction_min3.begin(),sim_fraction_min3.end()) - sim_fraction_min3.begin()); 
             if(std::equal(sim_rechit_diff.begin() + 1, sim_rechit_diff.end(), sim_rechit_diff.begin())) pfCluster_sim_rechit_diff_MatchedIndex.push_back(-1);
             else pfCluster_sim_rechit_diff_MatchedIndex.push_back(std::max_element(sim_rechit_diff.begin(),sim_rechit_diff.end()) - sim_rechit_diff.begin());  
             if(std::equal(sim_rechit_fraction.begin() + 1, sim_rechit_fraction.end(), sim_rechit_fraction.begin())) pfCluster_sim_rechit_fraction_MatchedIndex.push_back(-1);
             else pfCluster_sim_rechit_fraction_MatchedIndex.push_back(std::max_element(sim_rechit_fraction.begin(),sim_rechit_fraction.end()) - sim_rechit_fraction.begin()); 
             if(std::equal(global_sim_rechit_fraction.begin() + 1, global_sim_rechit_fraction.end(), global_sim_rechit_fraction.begin())) pfCluster_global_sim_rechit_fraction_MatchedIndex.push_back(-1);
             else pfCluster_global_sim_rechit_fraction_MatchedIndex.push_back(std::max_element(global_sim_rechit_fraction.begin(),global_sim_rechit_fraction.end()) - global_sim_rechit_fraction.begin());
          }    
    
          iPFCl++;        
      } 
   }

   //save inverse of matchings
   if(saveCaloParticles_ && savePFCluster_){ 
      for(unsigned int iPF=0; iPF<pfCluster_dR_simScore_MatchedIndex.size(); iPF++)
          if(pfCluster_dR_genScore_MatchedIndex.at(iPF)>=0) genParticle_pfCluster_dR_genScore_MatchedIndex[pfCluster_dR_genScore_MatchedIndex.at(iPF)].push_back(iPF);
      for(int iGen=0; iGen<nGenParticles; iGen++) 
          if(genParticle_pfCluster_dR_genScore_MatchedIndex[iGen].empty()) genParticle_pfCluster_dR_genScore_MatchedIndex[iGen].push_back(-1); 
   } 
   
   if(saveCaloParticles_ && savePFCluster_){ 
      for(unsigned int iPF=0; iPF<pfCluster_dR_simScore_MatchedIndex.size(); iPF++)
          if(pfCluster_dR_simScore_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_dR_simScore_MatchedIndex[pfCluster_dR_simScore_MatchedIndex.at(iPF)].push_back(iPF);
      for(unsigned int iPF=0; iPF<pfCluster_n_shared_xtals_MatchedIndex.size(); iPF++)
          if(pfCluster_n_shared_xtals_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_n_shared_xtals_MatchedIndex[pfCluster_n_shared_xtals_MatchedIndex.at(iPF)].push_back(iPF);
      for(unsigned int iPF=0; iPF<pfCluster_sim_fraction_MatchedIndex.size(); iPF++)
          if(pfCluster_sim_fraction_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_sim_fraction_MatchedIndex[pfCluster_sim_fraction_MatchedIndex.at(iPF)].push_back(iPF);
      for(unsigned int iPF=0; iPF<pfCluster_sim_fraction_min1_MatchedIndex.size(); iPF++)
          if(pfCluster_sim_fraction_min1_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_sim_fraction_min1_MatchedIndex[pfCluster_sim_fraction_min1_MatchedIndex.at(iPF)].push_back(iPF);
      for(unsigned int iPF=0; iPF<pfCluster_sim_fraction_min3_MatchedIndex.size(); iPF++)
          if(pfCluster_sim_fraction_min3_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_sim_fraction_min3_MatchedIndex[pfCluster_sim_fraction_min3_MatchedIndex.at(iPF)].push_back(iPF);
      for(unsigned int iPF=0; iPF<pfCluster_sim_rechit_diff_MatchedIndex.size(); iPF++)
          if(pfCluster_sim_rechit_diff_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_sim_rechit_diff_MatchedIndex[pfCluster_sim_rechit_diff_MatchedIndex.at(iPF)].push_back(iPF);
      for(unsigned int iPF=0; iPF<pfCluster_sim_rechit_fraction_MatchedIndex.size(); iPF++)
          if(pfCluster_sim_rechit_fraction_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex[pfCluster_sim_rechit_fraction_MatchedIndex.at(iPF)].push_back(iPF); 
      for(unsigned int iPF=0; iPF<pfCluster_global_sim_rechit_fraction_MatchedIndex.size(); iPF++)
         if(pfCluster_global_sim_rechit_fraction_MatchedIndex.at(iPF)>=0) caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex[pfCluster_global_sim_rechit_fraction_MatchedIndex.at(iPF)].push_back(iPF);   

      for(int iCalo=0; iCalo<nCaloParticles; iCalo++){ 
          if(caloParticle_pfCluster_dR_simScore_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_dR_simScore_MatchedIndex[iCalo].push_back(-1);  
          if(caloParticle_pfCluster_n_shared_xtals_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_n_shared_xtals_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_pfCluster_sim_fraction_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_sim_fraction_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_pfCluster_sim_fraction_min1_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_sim_fraction_min1_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_pfCluster_sim_fraction_min3_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_sim_fraction_min3_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_pfCluster_sim_rechit_diff_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_sim_rechit_diff_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_sim_rechit_fraction_MatchedIndex[iCalo].push_back(-1);
          if(caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex[iCalo].empty()) caloParticle_pfCluster_global_sim_rechit_fraction_MatchedIndex[iCalo].push_back(-1); 
      }    
   } 
  
   //Save SuperClusters 
   locCov.clear();
   full5x5_locCov.clear();
   if(saveSuperCluster_){
      int iSC=0;
      //std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
      for(const auto& iSuperCluster : *(superClusterEB.product())){  

          dR_genScore.clear();
          dR_simScore.clear();
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_min1.clear();
          sim_fraction_min3.clear(); 
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();

          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          math::XYZPoint caloPos = iSuperCluster.seed()->position();
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          superCluster_ieta.push_back(eb_id.ieta());
          superCluster_iphi.push_back(eb_id.iphi());
          superCluster_iz.push_back(0);   
 
          if(saveShowerShapes_){
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             locCov = EcalClusterTools::localCovariances(caloBC, &(*(recHitsEB.product())), &(*topology));
             full5x5_locCov = noZS::EcalClusterTools::localCovariances(caloBC, &(*(recHitsEB.product())), &(*topology));
             superCluster_r9.push_back(reduceFloat(EcalClusterTools::e3x3(caloBC, &(*(recHitsEB.product())), &(*topology))/iSuperCluster.energy(),nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(noZS::EcalClusterTools::e3x3(caloBC, &(*(recHitsEB.product())), &(*topology))/iSuperCluster.energy(),nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(sqrt(locCov[0]),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(sqrt(full5x5_locCov[0]),nBits_));
             superCluster_sigmaIetaIphi.push_back(reduceFloat(locCov[1],nBits_));
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(full5x5_locCov[1],nBits_));
             superCluster_sigmaIphiIphi.push_back(reduceFloat((!edm::isFinite(locCov[2]) ? 0. : sqrt(locCov[2])),nBits_));
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat((!edm::isFinite(full5x5_locCov[2]) ? 0. : sqrt(full5x5_locCov[2])),nBits_));
          } 
         
          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_genScore.push_back(reduceFloat(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_)); 
                 else dR_genScore.push_back(reduceFloat(999.,nBits_));  
             }  
             if(saveScores_) superCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::equal(dR_genScore.begin() + 1, dR_genScore.end(), dR_genScore.begin())) superCluster_dR_genScore_MatchedIndex.push_back(-1);
             else superCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin());  
          } 
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<float> scores = getScores(&hitsAndEnergies_SuperClusterEB.at(iSC),&hitsAndEnergies_CaloPart.at(iCalo));         
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_simScore.push_back(reduceFloat(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_)); 
                 else dR_simScore.push_back(999.);
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 sim_fraction_min1.push_back(scores[5]);  
                 sim_fraction_min3.push_back(scores[6]);  
             } 
             if(saveScores_){
                superCluster_dR_simScore[iSC] = dR_simScore;  
                superCluster_n_shared_xtals[iSC] = n_shared_xtals;  
                superCluster_sim_fraction[iSC] = sim_fraction;  
                superCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
                superCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
                superCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
                superCluster_sim_fraction_min1[iSC] = sim_fraction_min1;  
                superCluster_sim_fraction_min3[iSC] = sim_fraction_min3; 
             }
             if(std::equal(dR_simScore.begin() + 1, dR_simScore.end(), dR_simScore.begin())) superCluster_dR_simScore_MatchedIndex.push_back(-1);
             else superCluster_dR_simScore_MatchedIndex.push_back(std::min_element(dR_simScore.begin(),dR_simScore.end()) - dR_simScore.begin());  
             if(std::equal(n_shared_xtals.begin() + 1, n_shared_xtals.end(), n_shared_xtals.begin())) superCluster_n_shared_xtals_MatchedIndex.push_back(-1);
             else superCluster_n_shared_xtals_MatchedIndex.push_back(std::max_element(n_shared_xtals.begin(),n_shared_xtals.end()) - n_shared_xtals.begin());  
             if(std::equal(sim_fraction.begin() + 1, sim_fraction.end(), sim_fraction.begin())) superCluster_sim_fraction_MatchedIndex.push_back(-1);
             else superCluster_sim_fraction_MatchedIndex.push_back(std::max_element(sim_fraction.begin(),sim_fraction.end()) - sim_fraction.begin()); 
             if(std::equal(sim_fraction_min1.begin() + 1, sim_fraction_min1.end(), sim_fraction_min1.begin())) superCluster_sim_fraction_min1_MatchedIndex.push_back(-1);
             else superCluster_sim_fraction_min1_MatchedIndex.push_back(std::max_element(sim_fraction_min1.begin(),sim_fraction_min1.end()) - sim_fraction_min1.begin());
             if(std::equal(sim_fraction_min3.begin() + 1, sim_fraction_min3.end(), sim_fraction_min3.begin())) superCluster_sim_fraction_min3_MatchedIndex.push_back(-1);
             else superCluster_sim_fraction_min3_MatchedIndex.push_back(std::max_element(sim_fraction_min3.begin(),sim_fraction_min3.end()) - sim_fraction_min3.begin());  
             if(std::equal(sim_rechit_diff.begin() + 1, sim_rechit_diff.end(), sim_rechit_diff.begin())) superCluster_sim_rechit_diff_MatchedIndex.push_back(-1);
             else superCluster_sim_rechit_diff_MatchedIndex.push_back(std::max_element(sim_rechit_diff.begin(),sim_rechit_diff.end()) - sim_rechit_diff.begin());  
             if(std::equal(sim_rechit_fraction.begin() + 1, sim_rechit_fraction.end(), sim_rechit_fraction.begin())) superCluster_sim_rechit_fraction_MatchedIndex.push_back(-1);
             else superCluster_sim_rechit_fraction_MatchedIndex.push_back(std::max_element(sim_rechit_fraction.begin(),sim_rechit_fraction.end()) - sim_rechit_fraction.begin()); 
             if(std::equal(global_sim_rechit_fraction.begin() + 1, global_sim_rechit_fraction.end(), global_sim_rechit_fraction.begin())) superCluster_global_sim_rechit_fraction_MatchedIndex.push_back(-1);
             else superCluster_global_sim_rechit_fraction_MatchedIndex.push_back(std::max_element(global_sim_rechit_fraction.begin(),global_sim_rechit_fraction.end()) - global_sim_rechit_fraction.begin());   
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
          n_shared_xtals.clear();
          sim_fraction.clear();
          sim_fraction_min1.clear();
          sim_fraction_min3.clear();  
          sim_rechit_diff.clear();
          sim_rechit_fraction.clear();
          global_sim_rechit_fraction.clear();
          iSC_tmp++;
        
          superCluster_energy.push_back(reduceFloat(iSuperCluster.energy(),nBits_));
          superCluster_eta.push_back(reduceFloat(iSuperCluster.eta(),nBits_));
          superCluster_phi.push_back(reduceFloat(iSuperCluster.phi(),nBits_));
          math::XYZPoint caloPos = iSuperCluster.seed()->position(); 
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          superCluster_ieta.push_back(ee_id.ix());
          superCluster_iphi.push_back(ee_id.iy());
          superCluster_iz.push_back(ee_id.zside());   

          if(saveShowerShapes_){ 
             reco::CaloCluster caloBC(*iSuperCluster.seed());  
             locCov = EcalClusterTools::localCovariances(caloBC, &(*(recHitsEE.product())), &(*topology));
             full5x5_locCov = noZS::EcalClusterTools::localCovariances(caloBC, &(*(recHitsEE.product())), &(*topology));
             superCluster_r9.push_back(reduceFloat(EcalClusterTools::e3x3(caloBC, &(*(recHitsEE.product())), &(*topology))/iSuperCluster.energy(),nBits_));
             superCluster_full5x5_r9.push_back(reduceFloat(noZS::EcalClusterTools::e3x3(caloBC, &(*(recHitsEE.product())), &(*topology))/iSuperCluster.energy(),nBits_));
             superCluster_sigmaIetaIeta.push_back(reduceFloat(sqrt(locCov[0]),nBits_));
             superCluster_full5x5_sigmaIetaIeta.push_back(reduceFloat(sqrt(full5x5_locCov[0]),nBits_));
             superCluster_sigmaIetaIphi.push_back(reduceFloat(locCov[1],nBits_));
             superCluster_full5x5_sigmaIetaIphi.push_back(reduceFloat(full5x5_locCov[1],nBits_));
             superCluster_sigmaIphiIphi.push_back(reduceFloat((!edm::isFinite(locCov[2]) ? 0. : sqrt(locCov[2])),nBits_));
             superCluster_full5x5_sigmaIphiIphi.push_back(reduceFloat((!edm::isFinite(full5x5_locCov[2]) ? 0. : sqrt(full5x5_locCov[2])),nBits_));
          }

          //compute scores  
          if(saveGenParticles_){
             for(unsigned int iGen=0; iGen<genParts.size(); iGen++){
                 if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_genScore.push_back(reduceFloat(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_)); 
                 else dR_genScore.push_back(reduceFloat(999.,nBits_));  
             }  
             if(saveScores_) superCluster_dR_genScore[iSC] = dR_genScore; 
             if(std::equal(dR_genScore.begin() + 1, dR_genScore.end(), dR_genScore.begin())) superCluster_dR_genScore_MatchedIndex.push_back(-1);
             else superCluster_dR_genScore_MatchedIndex.push_back(std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin());  
          }
          if(saveCaloParticles_){
             for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 caloParticle_position = calculateAndSetPositionActual(&hitsAndEnergies_CaloPart.at(iCalo), 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
                 std::vector<float> scores = getScores(&hitsAndEnergies_SuperClusterEE.at(iSC_tmp),&hitsAndEnergies_CaloPart.at(iCalo));         
                 if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi())<0.1) dR_simScore.push_back(reduceFloat(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_)); 
                 else dR_simScore.push_back(999.);
                 n_shared_xtals.push_back(scores[0]);  
                 sim_fraction.push_back(scores[1]);  
                 sim_rechit_diff.push_back(scores[2]); 
                 sim_rechit_fraction.push_back(scores[3]);           
                 global_sim_rechit_fraction.push_back(scores[4]);
                 sim_fraction_min1.push_back(scores[5]);  
                 sim_fraction_min3.push_back(scores[6]);  
             }
             if(saveScores_){
                superCluster_dR_simScore[iSC] = dR_simScore;  
                superCluster_n_shared_xtals[iSC] = n_shared_xtals;  
                superCluster_sim_fraction[iSC] = sim_fraction;  
                superCluster_sim_rechit_diff[iSC] = sim_rechit_diff; 
                superCluster_sim_rechit_fraction[iSC] = sim_rechit_fraction;           
                superCluster_global_sim_rechit_fraction[iSC] = global_sim_rechit_fraction;
                superCluster_sim_fraction_min1[iSC] = sim_fraction_min1;  
                superCluster_sim_fraction_min3[iSC] = sim_fraction_min3; 
             }   
             if(std::equal(dR_simScore.begin() + 1, dR_simScore.end(), dR_simScore.begin())) superCluster_dR_simScore_MatchedIndex.push_back(-1);
             else superCluster_dR_simScore_MatchedIndex.push_back(std::min_element(dR_simScore.begin(),dR_simScore.end()) - dR_simScore.begin());  
             if(std::equal(n_shared_xtals.begin() + 1, n_shared_xtals.end(), n_shared_xtals.begin())) superCluster_n_shared_xtals_MatchedIndex.push_back(-1);
             else superCluster_n_shared_xtals_MatchedIndex.push_back(std::max_element(n_shared_xtals.begin(),n_shared_xtals.end()) - n_shared_xtals.begin());  
             if(std::equal(sim_fraction.begin() + 1, sim_fraction.end(), sim_fraction.begin())) superCluster_sim_fraction_MatchedIndex.push_back(-1);
             else superCluster_sim_fraction_MatchedIndex.push_back(std::max_element(sim_fraction.begin(),sim_fraction.end()) - sim_fraction.begin()); 
             if(std::equal(sim_fraction_min1.begin() + 1, sim_fraction_min1.end(), sim_fraction_min1.begin())) superCluster_sim_fraction_min1_MatchedIndex.push_back(-1);
             else superCluster_sim_fraction_min1_MatchedIndex.push_back(std::max_element(sim_fraction_min1.begin(),sim_fraction_min1.end()) - sim_fraction_min1.begin());
             if(std::equal(sim_fraction_min3.begin() + 1, sim_fraction_min3.end(), sim_fraction_min3.begin())) superCluster_sim_fraction_min3_MatchedIndex.push_back(-1);
             else superCluster_sim_fraction_min3_MatchedIndex.push_back(std::max_element(sim_fraction_min3.begin(),sim_fraction_min3.end()) - sim_fraction_min3.begin());  
             if(std::equal(sim_rechit_diff.begin() + 1, sim_rechit_diff.end(), sim_rechit_diff.begin())) superCluster_sim_rechit_diff_MatchedIndex.push_back(-1);
             else superCluster_sim_rechit_diff_MatchedIndex.push_back(std::max_element(sim_rechit_diff.begin(),sim_rechit_diff.end()) - sim_rechit_diff.begin());  
             if(std::equal(sim_rechit_fraction.begin() + 1, sim_rechit_fraction.end(), sim_rechit_fraction.begin())) superCluster_sim_rechit_fraction_MatchedIndex.push_back(-1);
             else superCluster_sim_rechit_fraction_MatchedIndex.push_back(std::max_element(sim_rechit_fraction.begin(),sim_rechit_fraction.end()) - sim_rechit_fraction.begin()); 
             if(std::equal(global_sim_rechit_fraction.begin() + 1, global_sim_rechit_fraction.end(), global_sim_rechit_fraction.begin())) superCluster_global_sim_rechit_fraction_MatchedIndex.push_back(-1);
             else superCluster_global_sim_rechit_fraction_MatchedIndex.push_back(std::max_element(global_sim_rechit_fraction.begin(),global_sim_rechit_fraction.end()) - global_sim_rechit_fraction.begin());   
          }

          iSC++;  
      }
   }

   //save inverse of matchings
   if(saveCaloParticles_ && saveSuperCluster_){ 
      for(unsigned int iSC=0; iSC<superCluster_dR_simScore_MatchedIndex.size(); iSC++)
          if(superCluster_dR_genScore_MatchedIndex.at(iSC)>=0) genParticle_superCluster_dR_genScore_MatchedIndex[superCluster_dR_genScore_MatchedIndex.at(iSC)].push_back(iSC);
      for(int iGen=0; iGen<nGenParticles; iGen++) 
          if(genParticle_superCluster_dR_genScore_MatchedIndex[iGen].empty()) genParticle_superCluster_dR_genScore_MatchedIndex[iGen].push_back(-1); 
   } 
   
   if(saveCaloParticles_ && saveSuperCluster_){ 
      for(unsigned int iSC=0; iSC<superCluster_dR_simScore_MatchedIndex.size(); iSC++)
          if(superCluster_dR_simScore_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_dR_simScore_MatchedIndex[superCluster_dR_simScore_MatchedIndex.at(iSC)].push_back(iSC);
      for(unsigned int iSC=0; iSC<superCluster_n_shared_xtals_MatchedIndex.size(); iSC++)
          if(superCluster_n_shared_xtals_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_n_shared_xtals_MatchedIndex[superCluster_n_shared_xtals_MatchedIndex.at(iSC)].push_back(iSC);
      for(unsigned int iSC=0; iSC<superCluster_sim_fraction_MatchedIndex.size(); iSC++)
          if(superCluster_sim_fraction_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_sim_fraction_MatchedIndex[superCluster_sim_fraction_MatchedIndex.at(iSC)].push_back(iSC);
      for(unsigned int iSC=0; iSC<superCluster_sim_fraction_min1_MatchedIndex.size(); iSC++)
          if(superCluster_sim_fraction_min1_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_sim_fraction_min1_MatchedIndex[superCluster_sim_fraction_min1_MatchedIndex.at(iSC)].push_back(iSC);
      for(unsigned int iSC=0; iSC<superCluster_sim_fraction_min3_MatchedIndex.size(); iSC++)
          if(superCluster_sim_fraction_min3_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_sim_fraction_min3_MatchedIndex[superCluster_sim_fraction_min3_MatchedIndex.at(iSC)].push_back(iSC);
      for(unsigned int iSC=0; iSC<superCluster_sim_rechit_diff_MatchedIndex.size(); iSC++)
          if(superCluster_sim_rechit_diff_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_sim_rechit_diff_MatchedIndex[superCluster_sim_rechit_diff_MatchedIndex.at(iSC)].push_back(iSC);
      for(unsigned int iSC=0; iSC<superCluster_sim_rechit_fraction_MatchedIndex.size(); iSC++)
          if(superCluster_sim_rechit_fraction_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_sim_rechit_fraction_MatchedIndex[superCluster_sim_rechit_fraction_MatchedIndex.at(iSC)].push_back(iSC); 
      for(unsigned int iSC=0; iSC<superCluster_global_sim_rechit_fraction_MatchedIndex.size(); iSC++)
         if(superCluster_global_sim_rechit_fraction_MatchedIndex.at(iSC)>=0) caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex[superCluster_global_sim_rechit_fraction_MatchedIndex.at(iSC)].push_back(iSC);   

      for(int iCalo=0; iCalo<nCaloParticles; iCalo++){ 
          if(caloParticle_superCluster_dR_simScore_MatchedIndex[iCalo].empty()) caloParticle_superCluster_dR_simScore_MatchedIndex[iCalo].push_back(-1);  
          if(caloParticle_superCluster_n_shared_xtals_MatchedIndex[iCalo].empty()) caloParticle_superCluster_n_shared_xtals_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_superCluster_sim_fraction_MatchedIndex[iCalo].empty()) caloParticle_superCluster_sim_fraction_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_superCluster_sim_fraction_min1_MatchedIndex[iCalo].empty()) caloParticle_superCluster_sim_fraction_min1_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_superCluster_sim_fraction_min3_MatchedIndex[iCalo].empty()) caloParticle_superCluster_sim_fraction_min3_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_superCluster_sim_rechit_diff_MatchedIndex[iCalo].empty()) caloParticle_superCluster_sim_rechit_diff_MatchedIndex[iCalo].push_back(-1); 
          if(caloParticle_superCluster_sim_rechit_fraction_MatchedIndex[iCalo].empty()) caloParticle_superCluster_sim_rechit_fraction_MatchedIndex[iCalo].push_back(-1);
          if(caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex[iCalo].empty()) caloParticle_superCluster_global_sim_rechit_fraction_MatchedIndex[iCalo].push_back(-1); 
      }    
   }

   //save clusters and superClusters mutual info
   if(savePFCluster_ && saveSuperCluster_){
      for(unsigned int iSC=0; iSC<hits_SuperCluster.size(); iSC++){
          for(unsigned int iPF=0; iPF<hits_PFCluster.size(); iPF++){
              bool isMatched = false; 
              for(unsigned int i=0; i<hits_SuperCluster.at(iSC).size(); i++){
                  std::vector<DetId>::iterator it = std::find(hits_PFCluster.at(iPF).begin(), hits_PFCluster.at(iPF).end(), hits_SuperCluster.at(iSC).at(i)); 
                  if (it != hits_PFCluster.at(iPF).end()){
                      isMatched = true;    
                      break;
                  }  
              }
              if(isMatched) superCluster_pfClustersIndex[iSC].push_back(iPF);
          }
          if(superCluster_pfClustersIndex[iSC].empty()) superCluster_pfClustersIndex[iSC].push_back(-1);    
      }
 
      for(unsigned int iSC=0; iSC<superCluster_pfClustersIndex.size(); iSC++)
          for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex.at(iSC).size(); iPF++)
              if(superCluster_pfClustersIndex[iSC].at(iPF)>=0) pfCluster_superClustersIndex[superCluster_pfClustersIndex[iSC].at(iPF)].push_back(iSC);
      for(unsigned int iPF=0; iPF<pfCluster_superClustersIndex.size(); iPF++)
          if(pfCluster_superClustersIndex[iPF].empty()) pfCluster_superClustersIndex[iPF].push_back(-1);
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

void DeepClusteringDumper::beginJob()
{

}

void DeepClusteringDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<std::pair<DetId, float> >* DeepClusteringDumper::getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CaloPart_tmp = new std::vector<std::pair<DetId, float> >;
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_tmp = new std::vector<std::pair<DetId, float> >;
    std::map<DetId, float> HitsAndEnergies_map;
    
    const auto& simClusters = iCaloParticle->simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int i = 0; i < hits_and_energies.size(); i++){  
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

std::vector<std::pair<DetId, float> >* DeepClusteringDumper::getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
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


std::vector<std::pair<DetId, float> >* DeepClusteringDumper::getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
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

std::vector<float> DeepClusteringDumper::getScores(const std::vector<std::pair<DetId, float> >*hits_and_energies_Cluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart)
{
    std::vector<float> scores;
    scores.resize(5);

    float nSharedXtals=-1.;
    float simFraction=-1.;
    float sim_rechit_diff=0.;
    float sim_rechit_fraction=0.;     
    float global_sim_rechit_fraction=-1.;       
   
    float rechits_tot_CaloPart = 0.;
    float rechits_tot_CaloPart_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart) {
        rechits_tot_CaloPart+=hit_CaloPart.second;
        rechits_tot_CaloPart_noEnergy+=1.;
    }

    float rechits_tot_Cluster = 0.;
    float rechits_tot_Cluster_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster) {
        rechits_tot_Cluster+=hit_Cluster.second;
        rechits_tot_Cluster_noEnergy+=1.;
    }
   
    float rechits_match_Cluster = 0.;
    float rechits_match_CaloPart = 0.;
    float rechits_match_CaloPart_noEnergy = 0.;
    for(const std::pair<DetId, float>& hit_Cluster : *hits_and_energies_Cluster){
        float reco_ratio=0.; 
        if(rechits_tot_Cluster!=0.) reco_ratio = hit_Cluster.second/rechits_tot_Cluster;
        float sim_ratio = 0.;        
        for(const std::pair<DetId, float>& hit_CaloPart : *hits_and_energies_CaloPart){  
            if(hit_CaloPart.first.rawId() == hit_Cluster.first.rawId()){

               rechits_match_Cluster += hit_Cluster.second;
               rechits_match_CaloPart += hit_CaloPart.second;    
               rechits_match_CaloPart_noEnergy += 1.0;

               sim_rechit_diff += fabs(hit_CaloPart.second-hit_Cluster.second);
               if(rechits_tot_CaloPart!=0.) sim_ratio = hit_CaloPart.second/rechits_tot_CaloPart; 
            }         
        } 
        sim_rechit_fraction += fabs(sim_ratio - reco_ratio);  
    }

    if(rechits_tot_CaloPart_noEnergy!=0.) nSharedXtals = rechits_match_CaloPart_noEnergy;

    if(rechits_tot_CaloPart!=0.) simFraction = rechits_match_CaloPart/rechits_tot_CaloPart;

    if(rechits_match_CaloPart_noEnergy!=0.) sim_rechit_diff = 1-(1./rechits_match_CaloPart_noEnergy)*sim_rechit_diff;
    else sim_rechit_diff=-1.; 

    if(sim_rechit_fraction!=0.) sim_rechit_fraction = 1-sim_rechit_fraction;
    else sim_rechit_fraction=-1.;
    
    if(rechits_tot_CaloPart!=0. && rechits_tot_Cluster!=0. && rechits_match_CaloPart/rechits_tot_CaloPart!=0. && rechits_match_Cluster/rechits_tot_Cluster!=0.) global_sim_rechit_fraction = 1-fabs(rechits_match_CaloPart/rechits_tot_CaloPart - rechits_match_Cluster/rechits_tot_Cluster);
    
    scores[0] = nSharedXtals;
    scores[1] = simFraction;
    scores[2] = sim_rechit_diff;
    scores[3] = sim_rechit_fraction;     
    scores[4] = global_sim_rechit_fraction;  
    if(simFraction>0.01) scores[5] = simFraction; 
    else scores[5] = -1.; 
    if(simFraction>0.03) scores[6] = simFraction; 
    else scores[6] = -1.; 

    return scores;
}

GlobalPoint DeepClusteringDumper::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES)
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

float DeepClusteringDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(DeepClusteringDumper);

