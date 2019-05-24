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
#include "RecoSimStudies/Dumpers/interface/RecoSimDumper.h"

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
   genToken_                = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   caloPartToken_           = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   ebRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   ebSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   
   useRechits_              = iConfig.getParameter<bool>("useRechits");
   useSuperCluster_         = iConfig.getParameter<bool>("useSuperCluster");
   motherID_                = iConfig.getParameter<int>("motherID");

   //output file, historgrams and trees
   tree                     = iFile->make<TTree>("caloTree","caloTree"); 

   tree->Branch("genParticle_id", &genParticle_id,"genParticle_id/i");
   tree->Branch("genParticle_energy", &genParticle_energy,"genParticle_energy/F");
   tree->Branch("genParticle_pt", &genParticle_pt,"genParticle_pt/F");
   tree->Branch("genParticle_eta", &genParticle_eta,"genParticle_eta/F");
   tree->Branch("genParticle_phi", &genParticle_phi,"genParticle_phi/F");
   tree->Branch("caloParticle_energy", &caloParticle_energy,"caloParticle_energy/F");
   tree->Branch("caloParticle_pt", &caloParticle_pt,"caloParticle_pt/F");
   tree->Branch("caloParticle_eta", &caloParticle_eta,"caloParticle_eta/F");
   tree->Branch("caloParticle_phi", &caloParticle_phi,"caloParticle_phi/F");
   tree->Branch("caloParticle_simHit_energy","std::vector<float>",&caloParticle_simHit_energy);
   tree->Branch("caloParticle_simHit_eta","std::vector<float>",&caloParticle_simHit_eta);
   tree->Branch("caloParticle_simHit_phi","std::vector<float>",&caloParticle_simHit_phi);
   tree->Branch("caloParticle_simHit_ieta","std::vector<int>",&caloParticle_simHit_ieta);
   tree->Branch("caloParticle_simHit_iphi","std::vector<int>",&caloParticle_simHit_iphi);
   tree->Branch("caloParticle_simHit_iz","std::vector<int>",&caloParticle_simHit_iz);
   if(useRechits_){ 
      tree->Branch("caloParticle_recHit_energy","std::vector<float>",&caloParticle_recHit_energy);
      tree->Branch("caloParticle_recHit_time","std::vector<float>",&caloParticle_recHit_time);  
   }
   if(useSuperCluster_){
      tree->Branch("superCluster_energy","std::vector<float>",&superCluster_energy);
      tree->Branch("superCluster_eta","std::vector<float>",&superCluster_eta);
      tree->Branch("superCluster_phi","std::vector<float>",&superCluster_phi);
      tree->Branch("superCluster_ieta","std::vector<int>",&superCluster_ieta);
      tree->Branch("superCluster_iphi","std::vector<int>",&superCluster_iphi);
      tree->Branch("superCluster_iz","std::vector<int>",&superCluster_iz);
      tree->Branch("superCluster_recHit_energy","std::vector<std::vector<float> >",&superCluster_recHit_energy);
      tree->Branch("superCluster_recHit_time","std::vector<std::vector<float> >",&superCluster_recHit_time);
      tree->Branch("superCluster_recHit_eta","std::vector<std::vector<float> >",&superCluster_recHit_eta);
      tree->Branch("superCluster_recHit_phi","std::vector<std::vector<float> >",&superCluster_recHit_phi);
      tree->Branch("superCluster_recHit_ieta","std::vector<std::vector<int> >",&superCluster_recHit_ieta);
      tree->Branch("superCluster_recHit_iphi","std::vector<std::vector<int> >",&superCluster_recHit_iphi);
      tree->Branch("superCluster_recHit_iz","std::vector<std::vector<int> >",&superCluster_recHit_iz);
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
   if(useRechits_) {
      if (!recHitsEB.isValid()) {
          std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
          return;
      }
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if(useRechits_) {
      if (!recHitsEE.isValid()) {
          std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEB;
   ev.getByToken(ebSuperClusterToken_, superClusterEB);
   if(useSuperCluster_) {
      if (!superClusterEB.isValid()) {
          std::cerr << "Analyze --> superClusterEB not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::SuperCluster> > superClusterEE;
   ev.getByToken(eeSuperClusterToken_, superClusterEE);
   if(useSuperCluster_) {
      if (!superClusterEE.isValid()) {
          std::cerr << "Analyze --> superClusterEE not found" << std::endl; 
          return;
      }
   } 

   for(const auto& iCalo : *(caloParticles.product()))
   {
       if(iCalo.pdgId()!=motherID_ && motherID_!=0) continue; 

       //SetBranch values to default
       setDefaultValues(); 
       caloParticle_simHit_energy.clear();
       caloParticle_simHit_eta.clear();
       caloParticle_simHit_phi.clear();
       caloParticle_simHit_ieta.clear();
       caloParticle_simHit_iphi.clear();
       caloParticle_simHit_iz.clear();
       caloParticle_recHit_energy.clear();
       caloParticle_recHit_time.clear();
       superCluster_energy.clear();
       superCluster_eta.clear();
       superCluster_phi.clear();
       superCluster_ieta.clear();
       superCluster_iphi.clear();
       superCluster_iz.clear();
       superCluster_recHit_energy.clear();
       superCluster_recHit_time.clear();
       superCluster_recHit_eta.clear();
       superCluster_recHit_phi.clear();
       superCluster_recHit_ieta.clear();
       superCluster_recHit_iphi.clear();
       superCluster_recHit_iz.clear();
       
       simHit_detIds_.clear();
       GlobalPoint cell;

       /*for(reco::GenParticleRefVector::iterator iGen = iCalo.genParticle_begin(); iGen != iCalo.genParticle_end(); ++iGen ){
           cout <<"iGen energy: " << (*iGen)->energy() << endl;
       }*/  
       
       const auto& genParticles = iCalo.genParticles();
       genParticle_id = iCalo.pdgId();
       if(genParticles.empty()){
          cout << "WARNING: no associated genParticle found" << endl;
          genParticle_energy = iCalo.energy();
          genParticle_pt = iCalo.pt();
          genParticle_eta = iCalo.eta();
          genParticle_phi = iCalo.phi();
       }else{
          genParticle_energy = (*genParticles.begin())->energy();
          genParticle_pt = (*genParticles.begin())->pt();
          genParticle_eta = (*genParticles.begin())->eta();
          genParticle_phi = (*genParticles.begin())->phi();
       }
 
       caloParticle_energy = iCalo.energy();
       caloParticle_pt = iCalo.pt();
       caloParticle_eta = iCalo.eta();
       caloParticle_phi = iCalo.phi();

       //Get hits from simClusters, and associated rechits 
       const auto& simClusters = iCalo.simClusters();
       for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
            auto simCluster = simClusters[iSC];  
            auto hits_and_fractions = simCluster->hits_and_fractions();
            for(unsigned int iHit= 0; iHit<hits_and_fractions.size(); iHit++){
                DetId id(hits_and_fractions[iHit].first);
                cell = geometry->getPosition(id);
                simHit_detIds_.push_back(id.rawId());
                caloParticle_simHit_energy.push_back(hits_and_fractions[iHit].second*simCluster->energy());
                caloParticle_simHit_eta.push_back(cell.eta());
                caloParticle_simHit_phi.push_back(cell.phi());
                if(id.subdetId()==EcalBarrel){
                   EBDetId eb_id(id);
                   caloParticle_simHit_ieta.push_back(eb_id.ieta());
                   caloParticle_simHit_iphi.push_back(eb_id.iphi());
                   caloParticle_simHit_iz.push_back(0);
                   if(useRechits_){
                      caloParticle_recHit_energy.push_back((*(recHitsEB.product())->find(id)).energy());
                      caloParticle_recHit_time.push_back((*(recHitsEB.product())->find(id)).time());
                   }
                }else if(id.subdetId()==EcalEndcap){
                   EEDetId ee_id(id);
                   int iz=0;
                   if(ee_id.zside()<0) iz=-1;
                   if(ee_id.zside()>0) iz=1;
                   caloParticle_simHit_ieta.push_back(ee_id.ix());
                   caloParticle_simHit_iphi.push_back(ee_id.iy());
                   caloParticle_simHit_iz.push_back(iz);
                   if(useRechits_){
                      caloParticle_recHit_energy.push_back((*(recHitsEE.product())->find(id)).energy());
                      caloParticle_recHit_time.push_back((*(recHitsEE.product())->find(id)).time());
                   } 
                }
            }    
       }
       
       //Match superclusters with caloParticles and get associated rechits 
       if(useSuperCluster_){
          if((*superClusterEB.product()).size() == 0 && (*superClusterEE.product()).size() == 0) cout << "WARNING: no superclusters" << endl;
          else{
             std::vector<reco::SuperCluster> skimmedSuperClustersEB = matchedSuperClusters(&simHit_detIds_,*superClusterEB.product());
             std::vector<reco::SuperCluster> skimmedSuperClustersEE = matchedSuperClusters(&simHit_detIds_,*superClusterEE.product());
             if(skimmedSuperClustersEB.size()==0 && skimmedSuperClustersEE.size()==0) cout << "WARNING: no matched superclusters" << endl;
             else{
                for(unsigned int iSC=0; iSC<skimmedSuperClustersEB.size(); iSC++){

                    const DetId& seed_id = skimmedSuperClustersEB[iSC].seed()->hitsAndFractions().at(0).first;
                    EBDetId eb_id(seed_id);
                    superCluster_energy.push_back(skimmedSuperClustersEB[iSC].energy());
                    superCluster_eta.push_back(skimmedSuperClustersEB[iSC].eta());
                    superCluster_phi.push_back(skimmedSuperClustersEB[iSC].phi());
                    superCluster_ieta.push_back(eb_id.ieta());
                    superCluster_iphi.push_back(eb_id.iphi());
                    superCluster_iz.push_back(0);
                    
                    recHit_energy_.clear();
                    recHit_time_.clear();
                    recHit_eta_.clear();
                    recHit_phi_.clear();
                    recHit_ieta_.clear();
                    recHit_iphi_.clear();
                    recHit_iz_.clear();

                    std::vector<DetId> recHits_id = superClusterXtalInfo(skimmedSuperClustersEB[iSC]);
                    for(unsigned int iId=0; iId<recHits_id.size(); iId++){
                        EBDetId recHit_EBid(recHits_id[iId]);
                        cell = geometry->getPosition(recHits_id[iId]);
                        recHit_energy_.push_back((*(recHitsEB.product())->find(recHits_id[iId])).energy());
                        recHit_time_.push_back((*(recHitsEB.product())->find(recHits_id[iId])).time());
                        recHit_eta_.push_back(cell.eta());
                        recHit_phi_.push_back(cell.phi());
                        recHit_ieta_.push_back(recHit_EBid.ieta());
                        recHit_iphi_.push_back(recHit_EBid.iphi());
                        recHit_iz_.push_back(0);
                    }

                    superCluster_recHit_energy.push_back(recHit_energy_);
                    superCluster_recHit_time.push_back(recHit_time_);
                    superCluster_recHit_eta.push_back(recHit_eta_);
                    superCluster_recHit_phi.push_back(recHit_phi_);
                    superCluster_recHit_ieta.push_back(recHit_ieta_);
                    superCluster_recHit_iphi.push_back(recHit_iphi_);
                    superCluster_recHit_iz.push_back(recHit_iz_);  
                } 
                for(unsigned int iSC=0; iSC<skimmedSuperClustersEE.size(); iSC++){

                    const DetId& seed_id = skimmedSuperClustersEE[iSC].seed()->hitsAndFractions().at(0).first;
                    EEDetId ee_id(seed_id);
                    int iz=0;
                    if(ee_id.zside()<0) iz=-1;
                    if(ee_id.zside()>0) iz=1;
                    superCluster_energy.push_back(skimmedSuperClustersEE[iSC].energy());
                    superCluster_eta.push_back(skimmedSuperClustersEE[iSC].eta());
                    superCluster_phi.push_back(skimmedSuperClustersEE[iSC].phi());
                    superCluster_ieta.push_back(ee_id.ix());
                    superCluster_iphi.push_back(ee_id.iy());
                    superCluster_iz.push_back(iz);

                    recHit_energy_.clear();
                    recHit_time_.clear();
                    recHit_eta_.clear();
                    recHit_phi_.clear();
                    recHit_ieta_.clear();
                    recHit_iphi_.clear();
                    recHit_iz_.clear();

                    std::vector<DetId> recHits_id = superClusterXtalInfo(skimmedSuperClustersEE[iSC]);
                    for(unsigned int iId=0; iId<recHits_id.size(); iId++){
                        EEDetId recHit_EEid(recHits_id[iId]);
                        if(recHit_EEid.zside()<0) iz=-1;
                        if(recHit_EEid.zside()>0) iz=1;
                        cell = geometry->getPosition(recHits_id[iId]);
                        recHit_energy_.push_back((*(recHitsEE.product())->find(recHits_id[iId])).energy());
                        recHit_time_.push_back((*(recHitsEE.product())->find(recHits_id[iId])).time());
                        recHit_eta_.push_back(cell.eta());
                        recHit_phi_.push_back(cell.phi());
                        recHit_ieta_.push_back(recHit_EEid.ix());
                        recHit_iphi_.push_back(recHit_EEid.iy());
                        recHit_iz_.push_back(iz);
                    }

                    superCluster_recHit_energy.push_back(recHit_energy_);
                    superCluster_recHit_time.push_back(recHit_time_);
                    superCluster_recHit_eta.push_back(recHit_eta_);
                    superCluster_recHit_phi.push_back(recHit_phi_);
                    superCluster_recHit_ieta.push_back(recHit_ieta_);
                    superCluster_recHit_iphi.push_back(recHit_iphi_);
                    superCluster_recHit_iz.push_back(recHit_iz_); 
                }  
             }
          }
       }
       tree->Fill();
   }
}

void RecoSimDumper::beginJob()
{

}

void RecoSimDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void RecoSimDumper::setDefaultValues() 
{
    genParticle_id=0;
    genParticle_energy=-999.;
    genParticle_pt=-999.;
    genParticle_eta=-999.;
    genParticle_phi=-999.; 

    caloParticle_energy=-999.;
    caloParticle_pt=-999.;
    caloParticle_eta=-999.;
    caloParticle_phi=-999.; 
}

std::vector<reco::SuperCluster> RecoSimDumper::matchedSuperClusters(std::vector<uint32_t>* caloParticle_simHit_detIds, const std::vector<reco::SuperCluster> superClusters)
{
    std::vector<reco::SuperCluster> skimmedSuperClusters;
    for(unsigned int iSC=0; iSC<superClusters.size(); iSC++){
        reco::CaloCluster caloBC = *superClusters[iSC].seed();
        const DetId& seed_id = superClusters[iSC].seed()->hitsAndFractions().at(0).first;
        if(std::find(caloParticle_simHit_detIds->begin(),caloParticle_simHit_detIds->end(),seed_id.rawId())!= caloParticle_simHit_detIds->end()) skimmedSuperClusters.push_back(superClusters[iSC]); 
    } 

    return skimmedSuperClusters;
}

std::vector<DetId> RecoSimDumper::superClusterXtalInfo(reco::SuperCluster iSC)
{
    std::vector<DetId> hitsAndFractions;
    for(reco::CaloCluster_iterator iBC = iSC.clustersBegin(); iBC != iSC.clustersEnd(); ++iBC){
        const std::vector< std::pair<DetId, float> > &seedrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < seedrechits.size(); i++ ){      
            hitsAndFractions.push_back(seedrechits[i].first);
        }
    }
    return hitsAndFractions;
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumper);



