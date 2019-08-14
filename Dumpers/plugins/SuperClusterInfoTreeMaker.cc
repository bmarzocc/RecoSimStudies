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

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoSimStudies/Dumpers/plugins/SuperClusterInfoTreeMaker.h"

#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

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
SuperClusterTreeMaker::SuperClusterTreeMaker(const edm::ParameterSet& iConfig)
{

   vtxToken_                = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   caloPartToken_           = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   
   doCompression_           = iConfig.getParameter<bool>("doCompression");
   nBits_                   = iConfig.getParameter<int>("nBits");
   genID_                   = iConfig.getParameter<std::vector<int>>("genID");
   doSimMatch_              = iConfig.getParameter<bool>("doSimMatch"); 

   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("SuperClusterTree","Dump of all available SC info"); 
   tree->Branch("N_ECALClusters", &N_ECALClusters, "N_ECALClusters/I");
   tree->Branch("N_PSClusters", &N_PSClusters, "N_PSClusters/I");
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   tree->Branch("scRawEnergy", &scRawEnergy, "scRawEnergy/F");
   tree->Branch("scCalibratedEnergy", &scCalibratedEnergy, "scCalibratedEnergy/F");
   tree->Branch("scPreshowerEnergy", &scPreshowerEnergy, "scPreshowerEnergy/F");
   tree->Branch("scEta", &scEta, "scEta/F");
   tree->Branch("scPhi", &scPhi, "scPhi/F");
   tree->Branch("scR", &scR, "scR/F");
   tree->Branch("scPhiWidth", &scPhiWidth, "scPhiWidth/F");
   tree->Branch("scEtaWidth", &scEtaWidth, "scEtaWidth/F");
   tree->Branch("scSeedRawEnergy", &scSeedRawEnergy, "scSeedRawEnergy/F");
   tree->Branch("scSeedCalibratedEnergy", &scSeedCalibratedEnergy, "scSeedCalibratedEnergy/F");
   tree->Branch("scSeedEta", &scSeedEta, "scSeedEta/F");
   tree->Branch("scSeedPhi", &scSeedPhi, "scSeedPhi/F");
   // ecal cluster information
   tree->Branch("clusterRawEnergy", "std::vector<float>", &clusterRawEnergy);
   tree->Branch("clusterCalibEnergy", "std::vector<float>", &clusterCalibEnergy);
   tree->Branch("clusterEta", "std::vector<float>", &clusterEta);
   tree->Branch("clusterPhi", "std::vector<float>", &clusterPhi);
   tree->Branch("clusterDPhiToSeed", "std::vector<float>", &clusterDPhiToSeed);
   tree->Branch("clusterDEtaToSeed", "std::vector<float>", &clusterDEtaToSeed);
   tree->Branch("clusterDPhiToCentroid", "std::vector<float>", &clusterDPhiToCentroid);
   tree->Branch("clusterDEtaToCentroid", "std::vector<float>", &clusterDEtaToCentroid);
   tree->Branch("clusterHitFractionSharedWithSeed","std::vector<float>", &clusterHitFractionSharedWithSeed);
   tree->Branch("clusterLeakage","std::vector<float>", &clusterLeakage);
   tree->Branch("clusterInMustache", "std::vector<int>", &clusterInMustache);
   tree->Branch("clusterInDynDPhi", "std::vector<int>", &clusterInDynDPhi);
   // preshower information
   tree->Branch("psClusterRawEnergy", "std::vector<float>", &psClusterRawEnergy);
   tree->Branch("psClusterEta", "std::vector<float>", &psClusterEta);
   tree->Branch("psClusterPhi", "std::vector<float>", &psClusterPhi);

   if(doSimMatch_) {
     tree->Branch("genEta", &genEta, "genEta/F");
     tree->Branch("genPhi", &genPhi, "genPhi/F");
     tree->Branch("genEnergy", &genEnergy, "genEnergy/F");
     tree->Branch("genDRToCentroid", &genDRToCentroid, "genDRToCentroid/F");
     tree->Branch("genDRToSeed", &genDRToSeed, "genDRToSeed/F");
     tree->Branch("clusterDPhiToGen", "std::vector<float>", &clusterDPhiToGen);
     tree->Branch("clusterDEtaToGen", "std::vector<float>", &clusterDEtaToGen);
     tree->Branch("clusterLeakageWrtSim", "std::vector<float>", &clusterLeakageWrtSim);
   }
}

SuperClusterTreeMaker::~SuperClusterTreeMaker()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void SuperClusterTreeMaker::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<CaloParticle> > caloParticles;
   ev.getByToken(caloPartToken_,caloParticles);
   if (!caloParticles.isValid()) {
       std::cerr << "Analyze --> caloParticles not found" << std::endl; 
       return;
   }

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   ev.getByToken(ebSuperClusterToken_, superClusterEB);
   if (!superClusterEB.isValid()) {
       std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
       return;
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   ev.getByToken(eeSuperClusterToken_, superClusterEE);
   if (!superClusterEE.isValid()) {
       std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
       return;
   } 

   nVtx = -1;
   nVtx = vertices->size();
   
   scRawEnergy = -1.;
   scCalibratedEnergy = -1.;
   scPreshowerEnergy = -1.;
   scEta = -99.;
   scPhi = -999.;
   scR = -999.;
   scPhiWidth = -1.;
   scEtaWidth = -1.;
   scSeedRawEnergy = -1.;
   scSeedCalibratedEnergy = -1.;
   scSeedEta = -99.;
   scSeedPhi = -999.;
   genEnergy = -1.;
   genEta = -99.;
   genPhi = -999.;
   genDRToCentroid = -1.;
   genDRToSeed = -1.;
   N_ECALClusters = -1;
   N_PSClusters = -1;

   clusterRawEnergy.clear();;
   clusterCalibEnergy.clear();
   clusterEta.clear();
   clusterPhi.clear();
   clusterDPhiToSeed.clear();
   clusterDEtaToSeed.clear();
   clusterDPhiToCentroid.clear();
   clusterDEtaToCentroid.clear();
   clusterDPhiToGen.clear();
   clusterDEtaToGen.clear();
   clusterHitFractionSharedWithSeed.clear();
   clusterLeakage.clear();
   clusterLeakageWrtSim.clear();
   clusterInMustache.clear();
   clusterInDynDPhi.clear();
   
   simMatched_.clear();
   findBestSimMatches(caloParticles,superClusterEB, &genID_);

   std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
   for(const auto& iSuperCluster : *(superClusterEB.product())){   

       scRawEnergy = reduceFloat(iSuperCluster.rawEnergy(),nBits_);
       scCalibratedEnergy = reduceFloat(iSuperCluster.energy(),nBits_);
       scPreshowerEnergy = reduceFloat(iSuperCluster.preshowerEnergy(),nBits_);
       scEta = reduceFloat(iSuperCluster.eta(),nBits_);
       scPhi = reduceFloat(iSuperCluster.phi(),nBits_);
       scR = reduceFloat(iSuperCluster.position().R(),nBits_);
       scPhiWidth = reduceFloat(iSuperCluster.phiWidth(),nBits_);
       scEtaWidth = reduceFloat(iSuperCluster.etaWidth(),nBits_);
       if(doSimMatch_){ 
          genEnergy = reduceFloat(caloPartEnergy(&simMatched_[iSuperCluster]),nBits_);
          genEta = reduceFloat(simMatched_[iSuperCluster].eta(),nBits_);
          genPhi = reduceFloat(simMatched_[iSuperCluster].phi(),nBits_);
          genDRToCentroid = reduceFloat(deltaR(simMatched_[iSuperCluster].eta(),simMatched_[iSuperCluster].phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_);
       } 
       if(iSuperCluster.seed().isAvailable()){ 
          scSeedRawEnergy = reduceFloat(iSuperCluster.seed()->energy(),nBits_);
          scSeedCalibratedEnergy = reduceFloat(iSuperCluster.seed()->correctedEnergy(),nBits_);
          scSeedEta = reduceFloat(iSuperCluster.seed()->eta(),nBits_);
          scSeedPhi = reduceFloat(iSuperCluster.seed()->phi(),nBits_);
          if(doSimMatch_) genDRToSeed = reduceFloat( deltaR(simMatched_[iSuperCluster].eta(), simMatched_[iSuperCluster].phi(), iSuperCluster.seed()->eta(), iSuperCluster.seed()->phi()),nBits_);
       }
       if(iSuperCluster.clusters().isAvailable()){
          N_ECALClusters = iSuperCluster.clusters().size();
          clusterRawEnergy.resize(N_ECALClusters);
          clusterCalibEnergy.resize(N_ECALClusters);
          clusterEta.resize(N_ECALClusters);
          clusterPhi.resize(N_ECALClusters);
          clusterDPhiToSeed.resize(N_ECALClusters);
          clusterDEtaToSeed.resize(N_ECALClusters);
          clusterDPhiToCentroid.resize(N_ECALClusters);
          clusterDEtaToCentroid.resize(N_ECALClusters);
          clusterDPhiToGen.resize(N_ECALClusters);
          clusterDEtaToGen.resize(N_ECALClusters);
          clusterHitFractionSharedWithSeed.resize(N_ECALClusters);
          clusterLeakage.resize(N_ECALClusters);
          clusterLeakageWrtSim.resize(N_ECALClusters);
          clusterInMustache.resize(N_ECALClusters);
          clusterInDynDPhi.resize(N_ECALClusters);  
          int iClus=0;
          for(unsigned int iBC=0; iBC<(unsigned int)N_ECALClusters; iBC++){
              if(!iSuperCluster.clusters()[iBC].isAvailable()) { continue; } 
              if(iSuperCluster.clusters()[iBC] != iSuperCluster.seed()) { continue; } 
              clusterRawEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->energy(),nBits_);
              clusterCalibEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->correctedEnergy(),nBits_);
              clusterEta[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta(),nBits_);
              clusterPhi[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->phi(),nBits_);
              clusterDPhiToSeed[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.seed()->phi()),nBits_);  
              clusterDEtaToSeed[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.seed()->eta(),nBits_);    
              clusterDPhiToCentroid[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.phi()),nBits_); 
              clusterDEtaToCentroid[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.eta(),nBits_);
              clusterHitFractionSharedWithSeed[iClus] = reduceFloat(getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),&iSuperCluster.seed()->hitsAndFractions()),nBits_); 
              clusterLeakage[iClus] = reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),getHitsAndFractionsSC(iSuperCluster)),nBits_);
              clusterInMustache[iClus] = (int)reco::MustacheKernel::inMustache(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(),iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              clusterInDynDPhi[iClus] = (int)reco::MustacheKernel::inDynamicDPhiWindow(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(), iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              if(doSimMatch_) clusterDPhiToGen[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - simMatched_[iSuperCluster].phi()),nBits_);    
              if(doSimMatch_) clusterDEtaToGen[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - simMatched_[iSuperCluster].eta(),nBits_);    
              if(doSimMatch_) clusterLeakageWrtSim[iClus] = reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),getHitsAndFractionsCaloPart(simMatched_[iSuperCluster])),nBits_);  
              iClus++;
          }
       }

       tree->Fill();
   }

   scRawEnergy = -1.;
   scCalibratedEnergy = -1.;
   scPreshowerEnergy = -1.;
   scEta = -99.;
   scPhi = -999.;
   scR = -999.;
   scPhiWidth = -1.;
   scEtaWidth = -1.;
   scSeedRawEnergy = -1.;
   scSeedCalibratedEnergy = -1.;
   scSeedEta = -99.;
   scSeedPhi = -999.;
   genEnergy = -1.;
   genEta = -99.;
   genPhi = -999.;
   genDRToCentroid = -1.;
   genDRToSeed = -1.;
   N_ECALClusters = -1;
   N_PSClusters = -1;

   clusterRawEnergy.clear();;
   clusterCalibEnergy.clear();
   clusterEta.clear();
   clusterPhi.clear();
   clusterDPhiToSeed.clear();
   clusterDEtaToSeed.clear();
   clusterDPhiToCentroid.clear();
   clusterDEtaToCentroid.clear();
   clusterDPhiToGen.clear();
   clusterDEtaToGen.clear();
   clusterHitFractionSharedWithSeed.clear();
   clusterLeakage.clear();
   clusterLeakageWrtSim.clear();
   clusterInMustache.clear();
   clusterInDynDPhi.clear();
   psClusterRawEnergy.clear();
   psClusterEta.clear();
   psClusterPhi.clear();

   simMatched_.clear();
   findBestSimMatches(caloParticles,superClusterEE, &genID_);
 
   std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;
   for(const auto& iSuperCluster : *(superClusterEE.product())){   

       scRawEnergy = reduceFloat(iSuperCluster.rawEnergy(),nBits_);
       scCalibratedEnergy = reduceFloat(iSuperCluster.energy(),nBits_);
       scPreshowerEnergy = reduceFloat(iSuperCluster.preshowerEnergy(),nBits_);
       scEta = reduceFloat(iSuperCluster.eta(),nBits_);
       scPhi = reduceFloat(iSuperCluster.phi(),nBits_);
       scR = reduceFloat(iSuperCluster.position().R(),nBits_);
       scPhiWidth = reduceFloat(iSuperCluster.phiWidth(),nBits_);
       scEtaWidth = reduceFloat(iSuperCluster.etaWidth(),nBits_);
       if(doSimMatch_){ 
          genEnergy = reduceFloat(caloPartEnergy(&simMatched_[iSuperCluster]),nBits_);
          genEta = reduceFloat(simMatched_[iSuperCluster].eta(),nBits_);
          genPhi = reduceFloat(simMatched_[iSuperCluster].phi(),nBits_);
          genDRToCentroid = reduceFloat(deltaR(simMatched_[iSuperCluster].eta(),simMatched_[iSuperCluster].phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_);
       } 
       if(iSuperCluster.seed().isAvailable()){ 
          scSeedRawEnergy = reduceFloat(iSuperCluster.seed()->energy(),nBits_);
          scSeedCalibratedEnergy = reduceFloat(iSuperCluster.seed()->correctedEnergy(),nBits_);
          scSeedEta = reduceFloat(iSuperCluster.seed()->eta(),nBits_);
          scSeedPhi = reduceFloat(iSuperCluster.seed()->phi(),nBits_);
          if(doSimMatch_) genDRToSeed = reduceFloat( deltaR(simMatched_[iSuperCluster].eta(), simMatched_[iSuperCluster].phi(), iSuperCluster.seed()->eta(), iSuperCluster.seed()->phi()),nBits_);
       }
       if(iSuperCluster.clusters().isAvailable()){
          N_ECALClusters = iSuperCluster.clusters().size();
          clusterRawEnergy.resize(N_ECALClusters);
          clusterCalibEnergy.resize(N_ECALClusters);
          clusterEta.resize(N_ECALClusters);
          clusterPhi.resize(N_ECALClusters);
          clusterDPhiToSeed.resize(N_ECALClusters);
          clusterDEtaToSeed.resize(N_ECALClusters);
          clusterDPhiToCentroid.resize(N_ECALClusters);
          clusterDEtaToCentroid.resize(N_ECALClusters);
          clusterDPhiToGen.resize(N_ECALClusters);
          clusterDEtaToGen.resize(N_ECALClusters);
          clusterHitFractionSharedWithSeed.resize(N_ECALClusters);
          clusterLeakage.resize(N_ECALClusters);
          clusterLeakageWrtSim.resize(N_ECALClusters);
          clusterInMustache.resize(N_ECALClusters);
          clusterInDynDPhi.resize(N_ECALClusters);  
          int iClus=0;
          for(unsigned int iBC=0; iBC<(unsigned int)N_ECALClusters; iBC++){
              if(!iSuperCluster.clusters()[iBC].isAvailable()) { continue; } 
              if(iSuperCluster.clusters()[iBC] != iSuperCluster.seed()) { continue; } 
              clusterRawEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->energy(),nBits_);
              clusterCalibEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->correctedEnergy(),nBits_);
              clusterEta[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta(),nBits_);
              clusterPhi[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->phi(),nBits_);
              clusterDPhiToSeed[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.seed()->phi()),nBits_);  
              clusterDEtaToSeed[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.seed()->eta(),nBits_);    
              clusterDPhiToCentroid[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.phi()),nBits_); 
              clusterDEtaToCentroid[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.eta(),nBits_);
              clusterHitFractionSharedWithSeed[iClus] = reduceFloat(getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),&iSuperCluster.seed()->hitsAndFractions()),nBits_); 
              clusterLeakage[iClus] = reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),getHitsAndFractionsSC(iSuperCluster)),nBits_);
              clusterInMustache[iClus] = (int)reco::MustacheKernel::inMustache(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(),iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              clusterInDynDPhi[iClus] = (int)reco::MustacheKernel::inDynamicDPhiWindow(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(), iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              if(doSimMatch_) clusterDPhiToGen[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - simMatched_[iSuperCluster].phi()),nBits_);    
              if(doSimMatch_) clusterDEtaToGen[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - simMatched_[iSuperCluster].eta(),nBits_);    
              if(doSimMatch_) clusterLeakageWrtSim[iClus] = reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),getHitsAndFractionsCaloPart(simMatched_[iSuperCluster])),nBits_);  
              iClus++;
          }
       }

       if(iSuperCluster.preshowerClusters().isAvailable()){
          N_PSClusters = iSuperCluster.preshowerClusters().size();
          psClusterRawEnergy.resize(N_PSClusters);
          psClusterEta.resize(N_PSClusters);
          psClusterPhi.resize(N_PSClusters); 
          int iClus=0;
          for(unsigned int iPC=0; iPC<(unsigned int)N_PSClusters; iPC++){
              if(!iSuperCluster.preshowerClusters()[iPC].isAvailable()) { continue; } 
              psClusterRawEnergy[iClus] = reduceFloat(iSuperCluster.preshowerClusters()[iPC]->energy(),nBits_);
              psClusterEta[iClus] = reduceFloat(iSuperCluster.preshowerClusters()[iPC]->eta(),nBits_);;
              psClusterPhi[iClus] = reduceFloat(iSuperCluster.preshowerClusters()[iPC]->phi(),nBits_);;   
              iClus++;
          }
       }

       tree->Fill();
   }
}

void SuperClusterTreeMaker::beginJob()
{

}

void SuperClusterTreeMaker::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SuperClusterTreeMaker::findBestSimMatches(edm::Handle<std::vector<CaloParticle> > caloParticles, edm::Handle<std::vector<reco::SuperCluster> > superClusters, std::vector<int>* genID_)
{
    std::vector<CaloParticle> caloParts;
    for(const auto& iCalo : *(caloParticles.product()))
    {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_->size(); id++)
           if(iCalo.pdgId()==genID_->at(id) || genID_->at(id)==0) isGoodParticle=true;
      
       if(!isGoodParticle) continue; 
       caloParts.push_back(iCalo); 
    }

    for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++)
    {
        double dE_min=999.;
        reco::SuperCluster bestmatch;
        for(const auto& iSuperCluster : *(superClusters.product())){   
            double dE = fabs(iSuperCluster.energy() - caloPartEnergy(&caloParts[iCalo]));
            if(dE < dE_min) {
               dE_min = dE;
               bestmatch = iSuperCluster;
            }
        }
        simMatched_[bestmatch]=caloParts[iCalo]; 
    }   
}

float SuperClusterTreeMaker::caloPartEnergy(CaloParticle* caloPart)
{
    float energy = 0.;
    const auto& simClusters = caloPart->simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_energies = simCluster->hits_and_energies();
        for(unsigned int iHit= 0; iHit<hits_and_energies.size(); iHit++){
            energy += hits_and_energies[iHit].second;
        }
    }

    return energy;
}

float SuperClusterTreeMaker::getSharedRecHitFraction(const std::vector<std::pair<DetId, float> >*hits_and_fractions_BC, const std::vector<std::pair<DetId, float> > *hits_and_fractions_Seed)
{
    float fraction = -1.;
    
    float rechits_tot_BC = hits_and_fractions_BC->size();
    float rechits_match_BC = 0.0;

    for(const std::pair<DetId, float>& hit_Seed : *hits_and_fractions_Seed) {
        for(const std::pair<DetId, float>& hit_BC : *hits_and_fractions_BC) {
            if(hit_Seed.first.rawId() == hit_BC.first.rawId()) {
               rechits_match_BC += 1.0;
            }
        }
    }

    fraction = rechits_match_BC/rechits_tot_BC;
    return fraction;
}

std::vector<std::pair<DetId, float> >* SuperClusterTreeMaker::getHitsAndFractionsSC(reco::SuperCluster iSuperCluster)
{
    std::vector<std::pair<DetId, float> >* HitsAndFractions_SC = new std::vector<std::pair<DetId, float> >;
    
    for(reco::CaloCluster_iterator iBC = iSuperCluster.clustersBegin(); iBC != iSuperCluster.clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < seedrechits.size(); i++)     
            HitsAndFractions_SC->push_back(seedrechits[i]);                      
    } 

    return HitsAndFractions_SC;
}

std::vector<std::pair<DetId, float> >* SuperClusterTreeMaker::getHitsAndFractionsCaloPart(CaloParticle iCaloParticle)
{
    std::vector<std::pair<DetId, float> >* HitsAndFractions_CP = new std::vector<std::pair<DetId, float> >;
    
    const auto& simClusters = iCaloParticle.simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
        auto simCluster = simClusters[iSC];  
        auto hits_and_fractions = simCluster->hits_and_fractions();
        for(unsigned int i = 0; i < hits_and_fractions.size(); i++)  
            HitsAndFractions_CP->push_back(hits_and_fractions[i]);   
    }

    return HitsAndFractions_CP;
}

float SuperClusterTreeMaker::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SuperClusterTreeMaker);

