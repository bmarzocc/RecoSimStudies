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
   genToken_                = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_           = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
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
     tree->Branch("genEta", "std::vector<float>", &genEta);
     tree->Branch("genPhi", "std::vector<float>", &genPhi);
     tree->Branch("genEnergy", "std::vector<float>", &genEnergy);
     tree->Branch("simEta", "std::vector<float>", &simEta);
     tree->Branch("simPhi", "std::vector<float>", &simPhi);
     tree->Branch("simEnergy", "std::vector<float>", &simEnergy);
     tree->Branch("simDRToCentroid", "std::vector<float>", &simDRToCentroid);
     tree->Branch("simDRToSeed", "std::vector<float>", &simDRToSeed);
     tree->Branch("simFractionBCtoBC", "std::vector<float>", &simFractionBCtoBC);
     tree->Branch("simFractionBCtoCP", "std::vector<float>", &simFractionBCtoCP);
     tree->Branch("simFractionCPtoBC", "std::vector<float>", &simFractionCPtoBC);
     tree->Branch("simFractionCPtoCP", "std::vector<float>", &simFractionCPtoCP);   
     tree->Branch("clusterDPhiToSim", "std::vector<std::vector<float> >", &clusterDPhiToSim);
     tree->Branch("clusterDEtaToSim", "std::vector<std::vector<float> >", &clusterDEtaToSim);
     tree->Branch("clusterLeakageWrtSim", "std::vector<std::vector<float> >", &clusterLeakageWrtSim);
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

   edm::ESHandle<CaloGeometry> caloGeometry;
   iSetup.get<CaloGeometryRecord>().get(caloGeometry);
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
   if (!recHitsEB.isValid()) {
       std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if (!recHitsEE.isValid()) {
       std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
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
   
   std::vector<CaloParticle> caloParts;
   for(const auto& iCalo : *(caloParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++)
           if(iCalo.pdgId()==genID_.at(id) || genID_.at(id)==0) isGoodParticle=true;
      
       if(!isGoodParticle) continue; 
       caloParts.push_back(iCalo); 
   } 
   std::vector<GenParticle> genParts;
   for(const auto& iGen : *(genParticles.product()))
   {
       bool isGoodParticle = false; 
       for(unsigned int id=0; id<genID_.size(); id++)
           if((iGen.pdgId()==genID_.at(id) || genID_.at(id)==0) && iGen.status()==1) isGoodParticle=true;
      
       if(!isGoodParticle) continue; 
       genParts.push_back(iGen); 
   } 

   std::vector<std::pair<DetId, float> >* hitsAndEnergies_SC;
   std::vector<std::pair<DetId, float> >* hitsAndEnergies_CP;
   GlobalPoint caloParticle_position;

   std::cout << "SuperClustersEB size: " << (superClusterEB.product())->size() << std::endl;
   for(const auto& iSuperCluster : *(superClusterEB.product())){   

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
       N_ECALClusters = -1;
       genEnergy.clear();
       genEta.clear();
       genPhi.clear();
       simEnergy.clear();
       simEta.clear();
       simPhi.clear();
       simDRToCentroid.clear();
       simDRToSeed.clear();
       simFractionBCtoBC.clear();
       simFractionBCtoCP.clear();
       simFractionCPtoBC.clear();
       simFractionCPtoCP.clear();
       clusterRawEnergy.clear();;
       clusterCalibEnergy.clear();
       clusterEta.clear();
       clusterPhi.clear();
       clusterDPhiToSeed.clear();
       clusterDEtaToSeed.clear();
       clusterDPhiToCentroid.clear();
       clusterDEtaToCentroid.clear();
       clusterDPhiToSim.clear();
       clusterDEtaToSim.clear();
       clusterHitFractionSharedWithSeed.clear();
       clusterLeakage.clear();
       clusterLeakageWrtSim.clear();
       clusterInMustache.clear();
       clusterInDynDPhi.clear();
       psClusterRawEnergy.clear();
       psClusterEta.clear();
       psClusterPhi.clear();
       
       hitsAndEnergies_SC = getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE);

       scRawEnergy = reduceFloat(iSuperCluster.rawEnergy(),nBits_);
       scCalibratedEnergy = reduceFloat(iSuperCluster.energy(),nBits_);
       scPreshowerEnergy = reduceFloat(iSuperCluster.preshowerEnergy(),nBits_);
       scEta = reduceFloat(iSuperCluster.eta(),nBits_);
       scPhi = reduceFloat(iSuperCluster.phi(),nBits_);
       scR = reduceFloat(iSuperCluster.position().R(),nBits_);
       scPhiWidth = reduceFloat(iSuperCluster.phiWidth(),nBits_);
       scEtaWidth = reduceFloat(iSuperCluster.etaWidth(),nBits_);
       GlobalPoint caloParticle_position;  
       if(doSimMatch_){
          for(unsigned int iGen=0; iGen<genParts.size(); iGen++)
          {
              genEnergy.push_back(reduceFloat(genParts.at(iGen).energy(),nBits_));
              genEta.push_back(reduceFloat(genParts.at(iGen).eta(),nBits_));
              genPhi.push_back(reduceFloat(genParts.at(iGen).phi(),nBits_));
          } 
          for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
              float simEnergy_tmp=0.;  
              hitsAndEnergies_CP = getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)));
              
              for(const std::pair<DetId, float>& hit_CP : *hitsAndEnergies_CP) 
                  simEnergy_tmp+=hit_CP.second;
        
              std::vector<float> fractionEnergy = getSharedRecHitFraction(hitsAndEnergies_SC,hitsAndEnergies_CP,true);
              
              caloParticle_position = calculateAndSetPositionActual(hitsAndEnergies_CP, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
              simEnergy.push_back(reduceFloat(simEnergy_tmp,nBits_));
              simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
              simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
              simDRToCentroid.push_back(reduceFloat(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_)); 
              simFractionBCtoBC.push_back(reduceFloat(fractionEnergy[0],nBits_)); 
              simFractionBCtoCP.push_back(reduceFloat(fractionEnergy[1],nBits_)); 
              simFractionCPtoBC.push_back(reduceFloat(fractionEnergy[2],nBits_)); 
              simFractionCPtoCP.push_back(reduceFloat(fractionEnergy[3],nBits_));  
              if(iSuperCluster.seed().isAvailable() && doSimMatch_) simDRToSeed.push_back(reduceFloat( deltaR(caloParticle_position.eta(),caloParticle_position.phi(), iSuperCluster.seed()->eta(), iSuperCluster.seed()->phi()),nBits_)); 
          }  
       } 
       if(iSuperCluster.seed().isAvailable()){ 
          scSeedRawEnergy = reduceFloat(iSuperCluster.seed()->energy(),nBits_);
          scSeedCalibratedEnergy = reduceFloat(iSuperCluster.seed()->correctedEnergy(),nBits_);
          scSeedEta = reduceFloat(iSuperCluster.seed()->eta(),nBits_);
          scSeedPhi = reduceFloat(iSuperCluster.seed()->phi(),nBits_);
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
          clusterDPhiToSim.resize(N_ECALClusters);
          clusterDEtaToSim.resize(N_ECALClusters);
          clusterHitFractionSharedWithSeed.resize(N_ECALClusters);
          clusterLeakage.resize(N_ECALClusters);
          clusterLeakageWrtSim.resize(N_ECALClusters);
          clusterInMustache.resize(N_ECALClusters);
          clusterInDynDPhi.resize(N_ECALClusters);  
          int iClus=0;
          for(unsigned int iBC=0; iBC<(unsigned int)N_ECALClusters; iBC++){
              if(!iSuperCluster.clusters()[iBC].isAvailable()) { continue; } 
              if(iSuperCluster.clusters()[iBC] == iSuperCluster.seed()) { continue; } 
              clusterRawEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->energy(),nBits_);
              clusterCalibEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->correctedEnergy(),nBits_);
              clusterEta[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta(),nBits_);
              clusterPhi[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->phi(),nBits_);
              clusterDPhiToSeed[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.seed()->phi()),nBits_);  
              clusterDEtaToSeed[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.seed()->eta(),nBits_);    
              clusterDPhiToCentroid[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.phi()),nBits_); 
              clusterDEtaToCentroid[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.eta(),nBits_);
              clusterHitFractionSharedWithSeed[iClus] = reduceFloat(getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),&iSuperCluster.seed()->hitsAndFractions(),true)[1],nBits_);  
              clusterLeakage[iClus] = reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE),true)[1],nBits_);
              clusterInMustache[iClus] = (int)reco::MustacheKernel::inMustache(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(),iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              clusterInDynDPhi[iClus] = (int)reco::MustacheKernel::inDynamicDPhiWindow(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(), iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              if(doSimMatch_){ 
                 for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 
                     hitsAndEnergies_CP = getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)));
                     caloParticle_position = calculateAndSetPositionActual(hitsAndEnergies_CP, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);

                     clusterDPhiToSim[iClus].push_back(reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - caloParticle_position.phi()),nBits_));    
                     clusterDEtaToSim[iClus].push_back(reduceFloat(iSuperCluster.clusters()[iBC]->eta() - caloParticle_position.eta(),nBits_));    
                     clusterLeakageWrtSim[iClus].push_back(reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),hitsAndEnergies_CP,false)[1],nBits_));    
                 }
              }
  
              iClus++;
          }
       }

       tree->Fill();
   }

   std::cout << "SuperClustersEE size: " << (superClusterEE.product())->size() << std::endl;
   for(const auto& iSuperCluster : *(superClusterEE.product())){   

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
       N_ECALClusters = -1;
       N_PSClusters = -1;  
       genEnergy.clear();
       genEta.clear();
       genPhi.clear();
       simEnergy.clear();
       simEta.clear();
       simPhi.clear();
       simDRToCentroid.clear();
       simDRToSeed.clear();
       simFractionBCtoBC.clear();
       simFractionBCtoCP.clear();
       simFractionCPtoBC.clear();
       simFractionCPtoCP.clear();
       clusterRawEnergy.clear();;
       clusterCalibEnergy.clear();
       clusterEta.clear();
       clusterPhi.clear();
       clusterDPhiToSeed.clear();
       clusterDEtaToSeed.clear();
       clusterDPhiToCentroid.clear();
       clusterDEtaToCentroid.clear();
       clusterDPhiToSim.clear();
       clusterDEtaToSim.clear();
       clusterHitFractionSharedWithSeed.clear();
       clusterLeakage.clear();
       clusterLeakageWrtSim.clear();
       clusterInMustache.clear();
       clusterInDynDPhi.clear();
       psClusterRawEnergy.clear();
       psClusterEta.clear();
       psClusterPhi.clear();
       
       hitsAndEnergies_SC = getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE);

       scRawEnergy = reduceFloat(iSuperCluster.rawEnergy(),nBits_);
       scCalibratedEnergy = reduceFloat(iSuperCluster.energy(),nBits_);
       scPreshowerEnergy = reduceFloat(iSuperCluster.preshowerEnergy(),nBits_);
       scEta = reduceFloat(iSuperCluster.eta(),nBits_);
       scPhi = reduceFloat(iSuperCluster.phi(),nBits_);
       scR = reduceFloat(iSuperCluster.position().R(),nBits_);
       scPhiWidth = reduceFloat(iSuperCluster.phiWidth(),nBits_);
       scEtaWidth = reduceFloat(iSuperCluster.etaWidth(),nBits_);
       GlobalPoint caloParticle_position;  
       if(doSimMatch_){ 
          for(const auto& iGen : *(genParticles.product()))
          {
              bool isGoodParticle = false; 
              for(unsigned int id=0; id<genID_.size(); id++)
                  if((iGen.pdgId()==genID_.at(id) || genID_.at(id)==0) && iGen.status()==1) isGoodParticle=true;
      
              if(!isGoodParticle) continue; 
              genEnergy.push_back(reduceFloat(iGen.energy(),nBits_));
              genEta.push_back(reduceFloat(iGen.eta(),nBits_));
              genPhi.push_back(reduceFloat(iGen.phi(),nBits_));
          }    
          for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
              float simEnergy_tmp=0.;  
              hitsAndEnergies_CP = getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)));
              
              for(const std::pair<DetId, float>& hit_CP : *hitsAndEnergies_CP) 
                  simEnergy_tmp+=hit_CP.second;
        
              std::vector<float> fractionEnergy = getSharedRecHitFraction(hitsAndEnergies_SC,hitsAndEnergies_CP,true);
              
              caloParticle_position = calculateAndSetPositionActual(hitsAndEnergies_CP, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
              simEnergy.push_back(reduceFloat(simEnergy_tmp,nBits_));
              simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
              simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
              simDRToCentroid.push_back(reduceFloat(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iSuperCluster.eta(),iSuperCluster.phi()),nBits_)); 
              simFractionBCtoBC.push_back(reduceFloat(fractionEnergy[0],nBits_)); 
              simFractionBCtoCP.push_back(reduceFloat(fractionEnergy[1],nBits_)); 
              simFractionCPtoBC.push_back(reduceFloat(fractionEnergy[2],nBits_)); 
              simFractionCPtoCP.push_back(reduceFloat(fractionEnergy[3],nBits_));  
              if(iSuperCluster.seed().isAvailable() && doSimMatch_) simDRToSeed.push_back(reduceFloat( deltaR(caloParticle_position.eta(),caloParticle_position.phi(), iSuperCluster.seed()->eta(), iSuperCluster.seed()->phi()),nBits_)); 
          }  
       } 
       if(iSuperCluster.seed().isAvailable()){ 
          scSeedRawEnergy = reduceFloat(iSuperCluster.seed()->energy(),nBits_);
          scSeedCalibratedEnergy = reduceFloat(iSuperCluster.seed()->correctedEnergy(),nBits_);
          scSeedEta = reduceFloat(iSuperCluster.seed()->eta(),nBits_);
          scSeedPhi = reduceFloat(iSuperCluster.seed()->phi(),nBits_);
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
          clusterDPhiToSim.resize(N_ECALClusters);
          clusterDEtaToSim.resize(N_ECALClusters);
          clusterHitFractionSharedWithSeed.resize(N_ECALClusters);
          clusterLeakage.resize(N_ECALClusters);
          clusterLeakageWrtSim.resize(N_ECALClusters);
          clusterInMustache.resize(N_ECALClusters);
          clusterInDynDPhi.resize(N_ECALClusters);  
          int iClus=0;
          for(unsigned int iBC=0; iBC<(unsigned int)N_ECALClusters; iBC++){
              if(!iSuperCluster.clusters()[iBC].isAvailable()) { continue; } 
              if(iSuperCluster.clusters()[iBC] == iSuperCluster.seed()) { continue; } 
              clusterRawEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->energy(),nBits_);
              clusterCalibEnergy[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->correctedEnergy(),nBits_);
              clusterEta[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta(),nBits_);
              clusterPhi[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->phi(),nBits_);
              clusterDPhiToSeed[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.seed()->phi()),nBits_);  
              clusterDEtaToSeed[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.seed()->eta(),nBits_);    
              clusterDPhiToCentroid[iClus] = reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - iSuperCluster.phi()),nBits_); 
              clusterDEtaToCentroid[iClus] = reduceFloat(iSuperCluster.clusters()[iBC]->eta() - iSuperCluster.eta(),nBits_);
              clusterHitFractionSharedWithSeed[iClus] = reduceFloat(getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),&iSuperCluster.seed()->hitsAndFractions(),true)[1],nBits_);  
              clusterLeakage[iClus] = reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),getHitsAndEnergiesSC(&iSuperCluster,recHitsEB,recHitsEE),true)[1],nBits_);
              clusterInMustache[iClus] = (int)reco::MustacheKernel::inMustache(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(),iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              clusterInDynDPhi[iClus] = (int)reco::MustacheKernel::inDynamicDPhiWindow(iSuperCluster.seed()->eta(),iSuperCluster.seed()->phi(), iSuperCluster.clusters()[iBC]->energy(),iSuperCluster.clusters()[iBC]->eta(),iSuperCluster.clusters()[iBC]->phi()); 
              if(doSimMatch_){ 
                 for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
                 
                     hitsAndEnergies_CP = getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)));
                     caloParticle_position = calculateAndSetPositionActual(hitsAndEnergies_CP, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);

                     clusterDPhiToSim[iClus].push_back(reduceFloat(TVector2::Phi_mpi_pi(iSuperCluster.clusters()[iBC]->phi() - caloParticle_position.phi()),nBits_));    
                     clusterDEtaToSim[iClus].push_back(reduceFloat(iSuperCluster.clusters()[iBC]->eta() - caloParticle_position.eta(),nBits_));    
                     clusterLeakageWrtSim[iClus].push_back(reduceFloat(1.-getSharedRecHitFraction(&iSuperCluster.clusters()[iBC]->hitsAndFractions(),hitsAndEnergies_CP,false)[1],nBits_));    
                 }
              }
  
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
std::vector<float> SuperClusterTreeMaker::getSharedRecHitFraction(const std::vector<std::pair<DetId, float> >*hits_and_energies_BC, const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, bool useEnergy)
{
    std::vector<float> fraction;
    fraction.resize(4);
   
    float rechits_tot_CP = 0.;
    for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) {
        if(useEnergy)  rechits_tot_CP+=hit_CP.second;
        if(!useEnergy) rechits_tot_CP+=1.;
    }

    float rechits_tot_BC = 0.;
    for(const std::pair<DetId, float>& hit_BC : *hits_and_energies_BC) {
        if(useEnergy)  rechits_tot_BC+=hit_BC.second;
        if(!useEnergy) rechits_tot_BC+=1.;
    }
   
    float rechits_match_BC = 0.;
    float rechits_match_CP = 0.;
    for(const std::pair<DetId, float>& hit_CP : *hits_and_energies_CP) 
        for(const std::pair<DetId, float>& hit_BC : *hits_and_energies_BC){      
            if(hit_CP.first.rawId() == hit_BC.first.rawId()){
               if(useEnergy){  
                  rechits_match_BC += hit_BC.second;
                  rechits_match_CP += hit_CP.second;    
               }
               if(!useEnergy){
                  rechits_match_BC += 1.0;
                  rechits_match_CP += 1.0;
               }
            } 
        } 
    
    if(rechits_tot_BC!=0.) fraction[0] = rechits_match_BC/rechits_tot_BC;
    else fraction[0]=-1.;
    if(rechits_tot_CP!=0.) fraction[1] = rechits_match_BC/rechits_tot_CP;
    else fraction[1]=-1.; 
    if(rechits_tot_BC!=0.) fraction[2] = rechits_match_CP/rechits_tot_BC;
    else fraction[2]=-1.;
    if(rechits_tot_CP!=0.) fraction[3] = rechits_match_CP/rechits_tot_CP;
    else fraction[3]=-1.;

    return fraction;
}

std::vector<std::pair<DetId, float> >* SuperClusterTreeMaker::getHitsAndEnergiesSC(const reco::SuperCluster* iSuperCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_SC = new std::vector<std::pair<DetId, float> >;
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
                    HitsAndEnergies_map[seedrechits.at(i).first]=HitsAndEnergies_map[seedrechits.at(i).first]+seedrechits.at(i).second * (*(recHitsEB.product())->find(seedrechits.at(i).first)).energy();
                } 
            }
        }                      
    } 

    for(auto const& hit : HitsAndEnergies_map) 
         HitsAndEnergies_SC->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_SC;
}

std::vector<std::pair<DetId, float> >* SuperClusterTreeMaker::getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CP = new std::vector<std::pair<DetId, float> >;
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
         HitsAndEnergies_CP->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_CP;
}

GlobalPoint SuperClusterTreeMaker::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES)
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

float SuperClusterTreeMaker::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SuperClusterTreeMaker);

