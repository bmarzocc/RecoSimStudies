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
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

#include "DataFormats/Math/interface/deltaR.h"

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

class SuperClusterTreeMaker : public edm::EDAnalyzer
{
      public:
         explicit SuperClusterTreeMaker(const edm::ParameterSet&);
	 ~SuperClusterTreeMaker();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
         virtual void endJob() ;
        
      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      float caloPartEnergy(CaloParticle* caloPart);
      void findBestSimMatches(edm::Handle<std::vector<CaloParticle> > caloParticles, edm::Handle<std::vector<reco::SuperCluster> > superCluster, std::vector<int>* genID_);
      float getSharedRecHitFraction(const std::vector<std::pair<DetId, float> >*hits_and_fractions_BC, const std::vector<std::pair<DetId, float> > *hits_and_fractions_Seed);
      std::vector<std::pair<DetId, float> >* getHitsAndFractionsSC(reco::SuperCluster iSuperCluster);  
      std::vector<std::pair<DetId, float> >* getHitsAndFractionsCaloPart(CaloParticle iCaloParticle);
       
      // ----------collection tokens-------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeSuperClusterToken_; 

      edm::Service<TFileService> iFile;

      // ----------config inputs-------------------
      bool doCompression_;
      int nBits_;
      std::vector<int> genID_;
      bool doSimMatch_;
      std::map<reco::SuperCluster, CaloParticle> simMatched_;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      int nVtx;
      float scRawEnergy;
      float scCalibratedEnergy;
      float scPreshowerEnergy;
      float scEta;
      float scPhi;
      float scR;
      float scPhiWidth;
      float scEtaWidth;
      float scSeedRawEnergy;
      float scSeedCalibratedEnergy;
      float scSeedEta;
      float scSeedPhi;
      float genEnergy;
      float genEta;
      float genPhi;
      float genDRToCentroid;
      float genDRToSeed;
      int N_ECALClusters;
      std::vector<float> clusterRawEnergy;
      std::vector<float> clusterCalibEnergy;
      std::vector<float> clusterEta;
      std::vector<float> clusterPhi;
      std::vector<float> clusterDPhiToSeed;
      std::vector<float> clusterDEtaToSeed;
      std::vector<float> clusterDPhiToCentroid;
      std::vector<float> clusterDEtaToCentroid;
      std::vector<float> clusterDPhiToGen;
      std::vector<float> clusterDEtaToGen;
      std::vector<float> clusterHitFractionSharedWithSeed;
      std::vector<float> clusterLeakage;
      std::vector<float> clusterLeakageWrtSim;
      std::vector<int> clusterInMustache;
      std::vector<int> clusterInDynDPhi;
      int N_PSClusters;
      std::vector<float> psClusterRawEnergy;
      std::vector<float> psClusterEta;
      std::vector<float> psClusterPhi;
};

#endif

