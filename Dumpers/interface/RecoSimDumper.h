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
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
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
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"

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
      void setDefaultValues();
      std::map<DetId,float> superClusterXtalInfo(reco::SuperCluster iSC);
      float reduceFloat(float val, int bits);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_; 
      edm::EDGetTokenT<std::vector<PCaloHit>  > PCaloHitEBToken_;
      edm::EDGetTokenT<std::vector<PCaloHit>  > PCaloHitEEToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFRecHit>  > pfRecHitToken_; 
      edm::EDGetTokenT<std::vector<reco::CaloCluster> > pfClusterToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeSuperClusterToken_; 

      edm::Service<TFileService> iFile;

      // ----------config inputs-------------------
      bool doCompression_;
      int nBits_;
      bool saveCalohits_;
      bool saveSimhits_;
      bool saveRechits_;
      bool savePFRechits_; 
      bool savePFCluster_;
      bool saveSuperCluster_;
      bool useEnergyRegression_;
      int motherID_;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      int genParticle_id;
      float genParticle_energy;
      float genParticle_pt;
      float genParticle_eta;
      float genParticle_phi;
      float caloParticle_energy;
      float caloParticle_pt;
      float caloParticle_eta;
      float caloParticle_phi;
      std::vector<float> caloHit_energy;
      std::vector<float> caloHit_time;
      std::vector<float> simHit_energy;
      std::vector<float> simHit_eta;
      std::vector<float> simHit_phi;
      std::vector<int> simHit_ieta;
      std::vector<int> simHit_iphi;
      std::vector<int> simHit_iz;
      std::vector<float> recHit_energy;
      std::vector<float> recHit_eta;
      std::vector<float> recHit_phi;
      std::vector<int> recHit_ieta;
      std::vector<int> recHit_iphi;
      std::vector<int> recHit_iz;
      std::vector<bool> pfRecHit_isMatched;
      std::vector<float> pfRecHit_energy;
      std::vector<float> pfRecHit_eta;
      std::vector<float> pfRecHit_phi;
      std::vector<int> pfRecHit_ieta;
      std::vector<int> pfRecHit_iphi;
      std::vector<int> pfRecHit_iz;
      std::vector<float> pfClusterHit_energy;
      std::vector<float> pfClusterHit_eta;
      std::vector<float> pfClusterHit_phi;
      std::vector<int> pfClusterHit_ieta;
      std::vector<int> pfClusterHit_iphi;
      std::vector<int> pfClusterHit_iz;
      std::vector<float> pfCluster_energy;
      std::vector<float> pfCluster_eta;
      std::vector<float> pfCluster_phi;
      std::vector<float> superClusterHit_energy;
      std::vector<float> superClusterHit_eta;
      std::vector<float> superClusterHit_phi;  
      std::vector<int> superClusterHit_ieta;
      std::vector<int> superClusterHit_iphi;    
      std::vector<int> superClusterHit_iz;    
      std::vector<float> superCluster_energy;
      std::vector<float> superCluster_eta;
      std::vector<float> superCluster_phi;
      std::map<int,int> map_simHit_pfCluster; 
      std::map<int,int> map_recHit_pfCluster;
      std::map<int,int> map_pfRecHit_pfCluster;
      std::map<int,int> map_pfClusterHit_pfCluster;   
      std::map<int,int> map_simHit_superCluster; 
      std::map<int,int> map_recHit_superCluster;
      std::map<int,int> map_pfRecHit_superCluster;
      std::map<int,int> map_superClusterHit_superCluster;    
};

#endif

