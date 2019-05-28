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
   pfRecHitToken_           = consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag>("pfRechitCollection")); 
   pfClusterToken_          = consumes<std::vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   ebSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("ebSuperClusterCollection"));
   eeSuperClusterToken_     = consumes<std::vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("eeSuperClusterCollection"));
   
   useRechits_              = iConfig.getParameter<bool>("useRechits");
   usePFRechits_            = iConfig.getParameter<bool>("usePFRechits"); 
   usePFCluster_            = iConfig.getParameter<bool>("usePFCluster");
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
   tree->Branch("simHit_energy","std::vector<float>",&simHit_energy);
   tree->Branch("simHit_eta","std::vector<float>",&simHit_eta);
   tree->Branch("simHit_phi","std::vector<float>",&simHit_phi);
   tree->Branch("simHit_ieta","std::vector<int>",&simHit_ieta);
   tree->Branch("simHit_iphi","std::vector<int>",&simHit_iphi);
   tree->Branch("simHit_iz","std::vector<int>",&simHit_iz);
   if(useRechits_) tree->Branch("recHit_energy","std::vector<float>",&recHit_energy);
   if(usePFRechits_) tree->Branch("pfRecHit_energy","std::vector<float>",&pfRecHit_energy);
   if(usePFCluster_){
      tree->Branch("pfClusterHit_energy","std::vector<float>",&pfClusterHit_energy);
      tree->Branch("pfCluster_energy","std::vector<float>",&pfCluster_energy);
      tree->Branch("pfCluster_eta","std::vector<float>",&pfCluster_eta);
      tree->Branch("pfCluster_phi","std::vector<float>",&pfCluster_phi);
      tree->Branch("map_simHit_pfCLuster","std::map<int,int>",&map_simHit_pfCLuster); 
   } 
   if(useSuperCluster_){
      tree->Branch("superClusterHit_energy","std::vector<float>",&superClusterHit_energy);
      tree->Branch("superCluster_energy","std::vector<float>",&superCluster_energy);
      tree->Branch("superCluster_eta","std::vector<float>",&superCluster_eta);
      tree->Branch("superCluster_phi","std::vector<float>",&superCluster_phi);
      tree->Branch("map_simHit_superCLuster","std::map<int,int>",&map_simHit_superCLuster); 
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

   edm::Handle<std::vector<reco::PFRecHit> > pfRecHits;
   ev.getByToken(pfRecHitToken_, pfRecHits);
   if(usePFRechits_) {
      if (!pfRecHits.isValid()) {
          std::cerr << "Analyze --> pfRecHits not found" << std::endl; 
          return;
      }
   } 

   edm::Handle<std::vector<reco::CaloCluster> > pfClusters;
   ev.getByToken(pfClusterToken_, pfClusters);
   if(usePFCluster_) {
      if (!pfClusters.isValid()) {
          std::cerr << "Analyze --> pfClusters not found" << std::endl; 
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
       simHit_energy.clear();
       simHit_eta.clear();
       simHit_phi.clear();
       simHit_ieta.clear();
       simHit_iphi.clear();
       simHit_iz.clear();
       recHit_energy.clear();
       pfRecHit_energy.clear();
       pfClusterHit_energy.clear();
       superClusterHit_energy.clear();
       pfCluster_energy.clear();
       pfCluster_eta.clear();
       pfCluster_phi.clear();
       superCluster_energy.clear(); 
       superCluster_eta.clear(); 
       superCluster_phi.clear(); 
      
       GlobalPoint cell;

       const auto& genParticles_caloPart = iCalo.genParticles();
       genParticle_id = iCalo.pdgId();
       if(genParticles_caloPart.empty()){
          cout << "WARNING: no associated genParticle found, making standard dR matching" << endl;
          float dR=999.;
          int igen_tmp=-1; 
          int igen=0; 
          for(const auto& iGen : *(genParticles.product()))
          {
              float dR_tmp = deltaR(iCalo.eta(),iCalo.phi(),iGen.eta(),iGen.phi());  
              if(dR_tmp<dR && iGen.status()==1){
                 dR=dR_tmp;
                 igen_tmp=igen;
              }  
              igen++;
          } 
          const auto& genParticles_tmp = *(genParticles.product());
          auto genParticle = genParticles_tmp[igen_tmp]; 
          genParticle_energy = genParticle.energy();
          genParticle_pt = genParticle.pt();
          genParticle_eta = genParticle.eta();
          genParticle_phi = genParticle.phi();
       }else{
          genParticle_energy = (*genParticles_caloPart.begin())->energy();
          genParticle_pt = (*genParticles_caloPart.begin())->pt();
          genParticle_eta = (*genParticles_caloPart.begin())->eta();
          genParticle_phi = (*genParticles_caloPart.begin())->phi();
       }
 
       caloParticle_energy = iCalo.energy();
       caloParticle_pt = iCalo.pt();
       caloParticle_eta = iCalo.eta();
       caloParticle_phi = iCalo.phi();

       //Get hits from simClusters, and associated recHits, pfRechits, PFClusterhit and superClusterhit 
       int simHit_index=-1; 
       const auto& simClusters = iCalo.simClusters();
       for( unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
            auto simCluster = simClusters[iSC];  
            auto hits_and_fractions = simCluster->hits_and_fractions();
            for(unsigned int iHit= 0; iHit<hits_and_fractions.size(); iHit++){
                DetId id(hits_and_fractions[iHit].first);
                if(id.subdetId()!=EcalBarrel && id.subdetId()!=EcalEndcap) continue;

                simHit_index++;
                int pfCluster_index=-999;
                int superCluster_index=-999;

                cell = geometry->getPosition(id);
                simHit_energy.push_back(hits_and_fractions[iHit].second*simCluster->energy());
                simHit_eta.push_back(cell.eta());
                simHit_phi.push_back(cell.phi());

                if(id.subdetId()==EcalBarrel){
                   EBDetId eb_id(id);
                   simHit_ieta.push_back(eb_id.ieta());
                   simHit_iphi.push_back(eb_id.iphi());
                   simHit_iz.push_back(0);

                   //Save associated recHit energy
                   float recHit_energy_ = -999.;
                   if(useRechits_){
                      if((*(recHitsEB.product())->find(id)).energy()>0.) recHit_energy_ = (*(recHitsEB.product())->find(id)).energy();
                      recHit_energy.push_back(recHit_energy_);
                   } 

                   //Save associated SuperClusterHit energy
                   float superClusterHit_energy_ = -999.;
                   if(useSuperCluster_){
                      for(const auto& iSuperCluster : *(superClusterEB.product())){
                          std::map<DetId,float> scInfos = superClusterXtalInfo(iSuperCluster);  
                          for(std::map<DetId,float>::iterator iter = scInfos.begin(); iter != scInfos.end(); ++iter)
                          {
                              if(iter->first.rawId() == id.rawId()){      
                                 superClusterHit_energy_ = iSuperCluster.energy()*iter->second;
                                 break;
                              }  
                          }
                      }
                      superClusterHit_energy.push_back(superClusterHit_energy_);
                   }

                }else if(id.subdetId()==EcalEndcap){
                   EEDetId ee_id(id);
                   int iz=0;
                   if(ee_id.zside()<0) iz=-1;
                   if(ee_id.zside()>0) iz=1;
                   simHit_ieta.push_back(ee_id.ix());
                   simHit_iphi.push_back(ee_id.iy());
                   simHit_iz.push_back(iz);

                   //Save associated recHit energy
                   float recHit_energy_ = -999.;
                   if(useRechits_){
                      if((*(recHitsEE.product())->find(id)).energy()>0.) recHit_energy_ = (*(recHitsEE.product())->find(id)).energy();
                      recHit_energy.push_back(recHit_energy_);
                   }
                   
                   //Save associated SuperClusterHit energy
                   float superClusterHit_energy_ = -999.;
                   int superCluster_index_tmp=0;
                   if(useSuperCluster_){
                      for(const auto& iSuperCluster : *(superClusterEE.product())){
                          std::map<DetId,float> scInfos = superClusterXtalInfo(iSuperCluster);  
                          for(std::map<DetId,float>::iterator iter = scInfos.begin(); iter != scInfos.end(); ++iter)
                          {
                              if(iter->first.rawId() == id.rawId()){      
                                 superClusterHit_energy_ = iSuperCluster.energy()*iter->second;
                                 superCluster_index = superCluster_index_tmp;
                                 break;
                              }  
                          }
                          superCluster_index_tmp++;
                      }
                      superClusterHit_energy.push_back(superClusterHit_energy_);
                   }
                }
                
                //Save associated pfRechit energy
                float pfRecHit_energy_ = -999.;
                if(usePFRechits_){
                   for(const auto& iPFRechit : *(pfRecHits.product())){
                       if(iPFRechit.detId() == id.rawId()){
                          pfRecHit_energy_ = iPFRechit.energy();
                          break;
                       }   
                   }  
                   pfRecHit_energy.push_back(pfRecHit_energy_);
                } 
                
                //Save associated PFClusterHit energy
                float pfClusterHit_energy_ = -999.;
                int pfCluster_index_tmp=0;
                if(usePFCluster_){
                   for(const auto& iPFCluster : *(pfClusters.product())){
                       reco::CaloCluster caloBC(iPFCluster);
                       const std::vector< std::pair<DetId,float> > &hitsAndFractions = caloBC.hitsAndFractions();
                       for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
                           if(hitsAndFractions[i].first.rawId() == id.rawId()){      
                              pfClusterHit_energy_ = iPFCluster.energy()*hitsAndFractions[i].second;
                              pfCluster_index = pfCluster_index_tmp;
                              break;
                           }  
                       }
                       pfCluster_index_tmp++;  
                   }  
                   pfClusterHit_energy.push_back(pfClusterHit_energy_);
                }
                //map simHit to PFCluster and to SuperCluster 
                map_simHit_pfCLuster.insert(pair<int,int>(simHit_index,pfCluster_index));
                map_simHit_superCLuster.insert(pair<int,int>(simHit_index,superCluster_index));
            }    
       }

       //Save PFClusters
       if(usePFCluster_){
          for(const auto& iPFCluster : *(pfClusters.product())){      
               pfCluster_energy.push_back(iPFCluster.energy());
               pfCluster_eta.push_back(iPFCluster.eta());
               pfCluster_phi.push_back(iPFCluster.phi());
          } 
       }
       
       //Save SuperClusters  
       if(useSuperCluster_){
          for(const auto& iSuperCluster : *(superClusterEE.product())){    
               superCluster_energy.push_back(iSuperCluster.energy());
               superCluster_eta.push_back(iSuperCluster.eta());
               superCluster_phi.push_back(iSuperCluster.phi());
          } 
       }

       //fill tree for each caloParticle
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

std::map<DetId,float> RecoSimDumper::superClusterXtalInfo(reco::SuperCluster iSC)
{
    std::map<DetId,float> hitsAndFractions;
    for(reco::CaloCluster_iterator iBC = iSC.clustersBegin(); iBC != iSC.clustersEnd(); ++iBC){
        const std::vector<std::pair<DetId,float> > &seedrechits = ( *iBC )->hitsAndFractions();
        for(unsigned int i = 0; i < seedrechits.size(); i++){      
                 hitsAndFractions[seedrechits[i].first] += seedrechits[i].second; 
        }
    }
    return hitsAndFractions;
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RecoSimDumper);



