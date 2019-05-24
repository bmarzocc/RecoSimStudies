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
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
      std::vector<reco::SuperCluster> matchedSuperClusters(std::vector<uint32_t>* simHit_detIds, const std::vector<reco::SuperCluster> superClusters);
      std::vector<DetId> superClusterXtalInfo(reco::SuperCluster iSC);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > ebSuperClusterToken_;
      edm::EDGetTokenT<std::vector<reco::SuperCluster> > eeSuperClusterToken_; 

      edm::Service<TFileService> iFile;

      // ----------config inputs-------------------
      bool useRechits_;
      bool useSuperCluster_;
      int motherID_;

      std::vector<uint32_t> simHit_detIds_;
      std::vector<float> recHit_energy_;
      std::vector<float> recHit_time_;
      std::vector<float> recHit_eta_;
      std::vector<float> recHit_phi_;
      std::vector<int> recHit_ieta_;
      std::vector<int> recHit_iphi_;
      std::vector<int> recHit_iz_;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      int genPart_id;
      float genPart_energy;
      float genPart_pt;
      float genPart_eta;
      float genPart_phi;
      float caloPart_energy;
      float caloPart_pt;
      float caloPart_eta;
      float caloPart_phi;
      std::vector<float> simHit_energy;
      std::vector<float> simHit_eta;
      std::vector<float> simHit_phi;
      std::vector<int> simHit_ieta;
      std::vector<int> simHit_iphi;
      std::vector<int> simHit_iz;
      std::vector<float> recHit_energy;
      std::vector<float> recHit_time;
      std::vector<float> superCluster_energy;
      std::vector<float> superCluster_eta;
      std::vector<float> superCluster_phi;
      std::vector<int> superCluster_ieta;
      std::vector<int> superCluster_iphi;
      std::vector<int> superCluster_iz;
      std::vector<std::vector<float> > superCluster_recHit_energy;
      std::vector<std::vector<float> > superCluster_recHit_time;
      std::vector<std::vector<float> > superCluster_recHit_eta;
      std::vector<std::vector<float> > superCluster_recHit_phi;
      std::vector<std::vector<int> > superCluster_recHit_ieta;
      std::vector<std::vector<int> > superCluster_recHit_iphi;
      std::vector<std::vector<int> > superCluster_recHit_iz;
};

#endif

