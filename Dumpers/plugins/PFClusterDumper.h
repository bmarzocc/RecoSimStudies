#ifndef RecoSimStudies_Dumpers_PFClusterDumper_H
#define RecoSimStudies_Dumpers_PFClusterDumper_H

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

class PFClusterDumper : public edm::EDAnalyzer
{
      public:
         explicit PFClusterDumper(const edm::ParameterSet&);
	 ~PFClusterDumper();
  
  
      private:
	 virtual void beginJob() ;
	 virtual void analyze(const edm::Event&, const edm::EventSetup&);
         virtual void endJob() ;
        
      // ----------additional functions-------------------
      std::vector<float> getScores(const std::vector<std::pair<DetId, float> >*hits_and_energies_Cluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle);
      std::vector<std::pair<DetId, float> >* getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE);
      GlobalPoint calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES);
      float reduceFloat(float val, int bits);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > caloPartToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<reco::PFCluster> > pfClusterToken_; 
      
      edm::Service<TFileService> iFile;
      const CaloSubdetectorGeometry* _ebGeom;
      const CaloSubdetectorGeometry* _eeGeom;
      const CaloSubdetectorGeometry* _esGeom;
      bool _esPlus;
      bool _esMinus;

      // ----------config inputs-------------------
      bool doCompression_;
      int nBits_;
      bool simMatch_energyFraction; 
      bool saveHitsPosition_;
      bool useES_;
      std::vector<int> genID_;
      bool saveScores_;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      float pfCluster_energy;
      float pfCluster_eta;
      float pfCluster_phi;
      int pfCluster_ieta;
      int pfCluster_iphi;
      int pfCluster_iz;
      std::vector<float> genEnergy;
      std::vector<float> genEta;
      std::vector<float> genPhi;
      std::vector<float> dR_genScore;
      int dR_genScore_MatchedIndex;
      std::vector<float> simEnergy;
      std::vector<float> simEta;
      std::vector<float> simPhi;
      std::vector<float> simIeta;
      std::vector<float> simIphi;
      std::vector<float> simIz;
      std::vector<float> dR_simScore;
      int dR_simScore_MatchedIndex;
      std::vector<float> n_shared_xtals;
      int n_shared_xtals_MatchedIndex;
      std::vector<float> sim_fraction;
      int sim_fraction_MatchedIndex;
      std::vector<float> sim_fraction_min1;
      int sim_fraction_min1_MatchedIndex;
      std::vector<float> sim_fraction_min3;
      int sim_fraction_min3_MatchedIndex;
      std::vector<float> sim_rechit_diff;
      int sim_rechit_diff_MatchedIndex;
      std::vector<float> sim_rechit_fraction;
      int sim_rechit_fraction_MatchedIndex;
      std::vector<float> global_sim_rechit_fraction;
      int global_sim_rechit_fraction_MatchedIndex;

      std::vector<std::vector<std::pair<DetId, float>>> hitsAndEnergies_CaloPart; 
};

#endif

