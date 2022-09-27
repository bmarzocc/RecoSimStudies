#ifndef RecoSimStudies_Dumpers_RechitDumper_H
#define RecoSimStudies_Dumpers_RechitDumper_H

// system include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

#include "DataFormats/Math/interface/libminifloat.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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

class RechitDumper : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
      public:
         explicit RechitDumper(const edm::ParameterSet&);
	 ~RechitDumper();
  
  
      private:
	 void beginJob() override;
	 void analyze(const edm::Event&, const edm::EventSetup&) override;
         void endJob() override;
        
      // ----------additional functions-------------------
      float reduceFloat(float val, int bits);
      void setTree(TTree* tree);
      void clearVectors(); 
      
      // ----------collection tokens-------------------
      edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
      edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
      edm::ESGetToken<EcalChannelStatus, EcalChannelStatusRcd> channelStatusToken_;

      const CaloGeometry* geometry;  
      const CaloTopology* topology;
      const EcalChannelStatus* chStatus;
      const CaloSubdetectorGeometry* ebGeometry;
      const CaloSubdetectorGeometry* eeGeometry;
      const CaloSubdetectorGeometry* hbGeometry;
      const CaloSubdetectorGeometry* heGeometry;
      const CaloSubdetectorGeometry* esGeometry;
      
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> ecalEBRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> ecalEERechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> ecalESRecHitToken_; 
      edm::EDGetTokenT<HBHERecHitCollection> hcalHBHERecHitToken_;   
    
      edm::Service<TFileService> iFile;
      
      // ----------config inputs-------------------
      bool isMC_;
      bool doCompression_;
      int nBits_;
      bool saveEB_;       
      bool saveEE_;      
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      long int eventId;
      int lumiId;
      int runId; 
      double truePU;
      double obsPU;
      int nVtx;
      float rho; 
      std::vector<uint32_t> ecalRecHit_rawId;
      std::vector<int> ecalRecHit_chStatus;
      std::vector<float> ecalRecHit_energy;
      std::vector<float> ecalRecHit_eta; 
      std::vector<float> ecalRecHit_phi;
      std::vector<int> ecalRecHit_ieta; 
      std::vector<int> ecalRecHit_iphi;
      std::vector<int> ecalRecHit_iz; 
      std::vector<std::vector<uint32_t> > matchedHcalRecHit_rawId;
      std::vector<std::vector<float> > matchedHcalRecHit_energy;
      std::vector<std::vector<float> > matchedHcalRecHit_eta; 
      std::vector<std::vector<float> > matchedHcalRecHit_phi;
      std::vector<std::vector<int> > matchedHcalRecHit_ieta; 
      std::vector<std::vector<int> > matchedHcalRecHit_iphi;
      std::vector<std::vector<int> > matchedHcalRecHit_iz; 
      std::vector<std::vector<int> > matchedHcalRecHit_depth; 
      std::vector<std::vector<uint32_t> > matchedESRecHit_rawId;
      std::vector<std::vector<float> > matchedESRecHit_energy;
      std::vector<std::vector<float> > matchedESRecHit_eta; 
      std::vector<std::vector<float> > matchedESRecHit_phi;
      std::vector<std::vector<int> > matchedESRecHit_ix; 
      std::vector<std::vector<int> > matchedESRecHit_iy;
      std::vector<std::vector<int> > matchedESRecHit_iz;
      std::vector<std::vector<int> > matchedESRecHit_strip; 
      std::vector<std::vector<int> > matchedESRecHit_plane; 
        
  
};

#endif
