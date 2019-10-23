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
#include "RecoSimStudies/Dumpers/plugins/PFClusterDumper.h"

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
PFClusterDumper::PFClusterDumper(const edm::ParameterSet& iConfig)
{

   caloPartToken_           = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
   genToken_                = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   ebRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_           = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   pfClusterToken_          = consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag>("pfClusterCollection")); 
   
   doCompression_           = iConfig.getParameter<bool>("doCompression");
   nBits_                   = iConfig.getParameter<int>("nBits");
   saveHitsPosition_        = iConfig.getParameter<bool>("saveHitsPosition");
   useES_                   = iConfig.getParameter<bool>("useES");
   genID_                   = iConfig.getParameter<std::vector<int>>("genID");
   saveScores_              = iConfig.getParameter<bool>("saveScores"); 

   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   tree->Branch("genEta", "std::vector<float>", &genEta);
   tree->Branch("genPhi", "std::vector<float>", &genPhi);
   tree->Branch("genEnergy", "std::vector<float>", &genEnergy);
   if(saveScores_) tree->Branch("dR_genScore", "std::vector<float>", &dR_genScore);   
   tree->Branch("simEta", "std::vector<float>", &simEta);
   tree->Branch("simPhi", "std::vector<float>", &simPhi);
   tree->Branch("simIeta", "std::vector<float>", &simIeta);
   tree->Branch("simIphi", "std::vector<float>", &simIphi);
   tree->Branch("simIz", "std::vector<float>", &simIz);
   tree->Branch("simEnergy", "std::vector<float>", &simEnergy);
   if(saveScores_){
      tree->Branch("dR_simScore", "std::vector<float>", &dR_simScore);
      tree->Branch("n_shared_xtals", "std::vector<float>", &n_shared_xtals);
      tree->Branch("sim_fraction", "std::vector<float>", &sim_fraction);
      tree->Branch("sim_fraction_min1", "std::vector<float>", &sim_fraction_min1);
      tree->Branch("sim_fraction_min3", "std::vector<float>", &sim_fraction_min3);
      tree->Branch("sim_rechit_diff", "std::vector<float>", &sim_rechit_diff);
      tree->Branch("sim_rechit_fraction", "std::vector<float>", &sim_rechit_fraction);   
      tree->Branch("global_sim_rechit_fraction", "std::vector<float>", &global_sim_rechit_fraction); 
   }  
   tree->Branch("dR_genScore_MatchedIndex",&dR_genScore_MatchedIndex,"dR_genScore_MatchedIndex/I");
   tree->Branch("dR_simScore_MatchedIndex",&dR_simScore_MatchedIndex,"dR_simScore_MatchedIndex/I");
   tree->Branch("n_shared_xtals_MatchedIndex",&n_shared_xtals_MatchedIndex,"n_shared_xtals_MatchedIndex/I");
   tree->Branch("sim_fraction_MatchedIndex",&sim_fraction_MatchedIndex,"sim_fraction_MatchedIndex/I");
   tree->Branch("sim_fraction_min1_MatchedIndex",&sim_fraction_min1_MatchedIndex,"sim_fraction_min1_MatchedIndex/I");
   tree->Branch("sim_fraction_min3_MatchedIndex",&sim_fraction_min3_MatchedIndex,"sim_fraction_min3_MatchedIndex/I");
   tree->Branch("sim_rechit_diff_MatchedIndex",&sim_rechit_diff_MatchedIndex,"sim_rechit_diff_MatchedIndex/I");
   tree->Branch("sim_rechit_fraction_MatchedIndex",&sim_rechit_fraction_MatchedIndex,"sim_rechit_fraction_MatchedIndex/I");   
   tree->Branch("global_sim_rechit_fraction_MatchedIndex",&global_sim_rechit_fraction_MatchedIndex,"global_sim_rechit_fraction_MatchedIndex/I"); 
   tree->Branch("pfCluster_energy",&pfCluster_energy,"pfCluster_energy/F");
   tree->Branch("pfCluster_eta",&pfCluster_eta,"pfCluster_eta/F");
   tree->Branch("pfCluster_phi",&pfCluster_phi,"pfCluster_phi/F"); 
   tree->Branch("pfCluster_ieta",&pfCluster_ieta,"pfCluster_ieta/I");
   tree->Branch("pfCluster_iphi",&pfCluster_iphi,"pfCluster_iphi/I");
   tree->Branch("pfCluster_iz",&pfCluster_iz,"pfCluster_iz/I");      
  
}

PFClusterDumper::~PFClusterDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void PFClusterDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{

   //calo geometry
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

   edm::Handle<std::vector<reco::PFCluster> > pfClusters;
   ev.getByToken(pfClusterToken_, pfClusters);
   if (!pfClusters.isValid()) {
       std::cerr << "Analyze --> pfClusters not found" << std::endl; 
       return;
   } 
  
   GlobalPoint cell;

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

   std::vector<std::pair<DetId, float> >* hitsAndEnergies_CP;
   std::vector<std::pair<DetId, float> >* hitsAndEnergies_BC;
   GlobalPoint caloParticle_position;

   //Save PFClusters
   std::cout << "PFClusters size: " << (pfClusters.product())->size() << std::endl;
   for(const auto& iPFCluster : *(pfClusters.product())){   

       pfCluster_energy=-1.;
       pfCluster_eta=-99.;
       pfCluster_phi=-99.;
       pfCluster_ieta=-99.;
       pfCluster_iphi=-99.;
       pfCluster_iz=-99.;
       genEnergy.clear();
       genEta.clear();
       genPhi.clear();
       dR_genScore.clear();
       simEnergy.clear();
       simEta.clear();
       simPhi.clear();
       simIeta.clear();
       simIphi.clear();
       simIz.clear();
       dR_simScore.clear();
       n_shared_xtals.clear();
       sim_fraction.clear();
       sim_fraction_min1.clear();
       sim_fraction_min3.clear(); 
       sim_rechit_diff.clear();
       sim_rechit_fraction.clear();
       global_sim_rechit_fraction.clear();
       dR_genScore_MatchedIndex=-1;
       dR_simScore_MatchedIndex=-1;
       n_shared_xtals_MatchedIndex=-1;
       sim_fraction_MatchedIndex=-1;
       sim_fraction_min1_MatchedIndex=-1;
       sim_fraction_min3_MatchedIndex=-1;   
       sim_rechit_diff_MatchedIndex=-1;
       sim_rechit_fraction_MatchedIndex=-1;
       global_sim_rechit_fraction_MatchedIndex=-1;
       
       pfCluster_energy=reduceFloat(iPFCluster.energy(),nBits_);
       pfCluster_eta=reduceFloat(iPFCluster.eta(),nBits_);
       pfCluster_phi=reduceFloat(iPFCluster.phi(),nBits_);

       reco::CaloCluster caloBC(iPFCluster);
       hitsAndEnergies_BC = getHitsAndEnergiesBC(&caloBC,recHitsEB,recHitsEE);

       math::XYZPoint caloPos = caloBC.position();
       if(std::abs(iPFCluster.eta()) < 1.479){  
          EBDetId eb_id(_ebGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));  
          pfCluster_ieta=eb_id.ieta();
          pfCluster_iphi=eb_id.iphi();
          pfCluster_iz=0; 
       }else{  
          EEDetId ee_id(_eeGeom->getClosestCell(GlobalPoint(caloPos.x(),caloPos.y(),caloPos.z())));   
          pfCluster_ieta=ee_id.ix();
          pfCluster_iphi=ee_id.iy();
          pfCluster_iz=ee_id.zside(); 
       }          

       for(unsigned int iGen=0; iGen<genParts.size(); iGen++)
       {
           genEnergy.push_back(reduceFloat(genParts.at(iGen).energy(),nBits_));
           genEta.push_back(reduceFloat(genParts.at(iGen).eta(),nBits_));
           genPhi.push_back(reduceFloat(genParts.at(iGen).phi(),nBits_));
           if(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_genScore.push_back(reduceFloat(deltaR(genParts.at(iGen).eta(),genParts.at(iGen).phi(),iPFCluster.eta(),iPFCluster.phi()),nBits_));
           else dR_genScore.push_back(reduceFloat(999.,nBits_));
       }
       if(std::equal(dR_genScore.begin() + 1, dR_genScore.end(), dR_genScore.begin())) dR_genScore_MatchedIndex = -1;
       else dR_genScore_MatchedIndex = std::min_element(dR_genScore.begin(),dR_genScore.end()) - dR_genScore.begin();  

       for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
           float simEnergy_tmp=0.;  
           hitsAndEnergies_CP = getHitsAndEnergiesCaloPart(&(caloParts.at(iCalo)));
              
           for(const std::pair<DetId, float>& hit_CP : *hitsAndEnergies_CP) 
               simEnergy_tmp+=hit_CP.second;

           caloParticle_position = calculateAndSetPositionActual(hitsAndEnergies_CP, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,false);
           simEnergy.push_back(reduceFloat(simEnergy_tmp,nBits_));
           simEta.push_back(reduceFloat(caloParticle_position.eta(),nBits_));
           simPhi.push_back(reduceFloat(caloParticle_position.phi(),nBits_));
           if(std::abs(iPFCluster.eta()) < 1.479){  
              EBDetId eb_id(_ebGeom->getClosestCell(caloParticle_position));  
              simIeta.push_back(eb_id.ieta());
              simIphi.push_back(eb_id.iphi());
              simIz.push_back(0); 
           }else{  
             EEDetId ee_id(_eeGeom->getClosestCell(caloParticle_position));   
             simIeta.push_back(ee_id.ix());
             simIphi.push_back(ee_id.iy());
             simIz.push_back(ee_id.zside());  
           }   

           std::vector<float> scores = getScores(hitsAndEnergies_BC,hitsAndEnergies_CP);         
           if(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi())<0.1) dR_simScore.push_back(reduceFloat(deltaR(caloParticle_position.eta(),caloParticle_position.phi(),iPFCluster.eta(),iPFCluster.phi()),nBits_)); 
           else dR_simScore.push_back(reduceFloat(999.,nBits_)); 
           n_shared_xtals.push_back(scores[0]);  
           sim_fraction.push_back(reduceFloat(scores[1],nBits_));  
           sim_rechit_diff.push_back(reduceFloat(scores[2],nBits_)); 
           sim_rechit_fraction.push_back(reduceFloat(scores[3],nBits_));           
           global_sim_rechit_fraction.push_back(reduceFloat(scores[4],nBits_));
           sim_fraction_min1.push_back(reduceFloat(scores[5],nBits_));  
           sim_fraction_min3.push_back(reduceFloat(scores[6],nBits_));  
       } 
       if(std::equal(dR_simScore.begin() + 1, dR_simScore.end(), dR_simScore.begin())) dR_simScore_MatchedIndex = -1;
       else dR_simScore_MatchedIndex = std::min_element(dR_simScore.begin(),dR_simScore.end()) - dR_simScore.begin();  
       if(std::equal(n_shared_xtals.begin() + 1, n_shared_xtals.end(), n_shared_xtals.begin())) n_shared_xtals_MatchedIndex = -1;
       else n_shared_xtals_MatchedIndex = std::max_element(n_shared_xtals.begin(),n_shared_xtals.end()) - n_shared_xtals.begin();  
       if(std::equal(sim_fraction.begin() + 1, sim_fraction.end(), sim_fraction.begin())) sim_fraction_MatchedIndex = -1;
       else sim_fraction_MatchedIndex = std::max_element(sim_fraction.begin(),sim_fraction.end()) - sim_fraction.begin(); 
       if(std::equal(sim_fraction_min1.begin() + 1, sim_fraction_min1.end(), sim_fraction_min1.begin())) sim_fraction_min1_MatchedIndex = -1;
       else sim_fraction_min1_MatchedIndex = std::max_element(sim_fraction_min1.begin(),sim_fraction_min1.end()) - sim_fraction_min1.begin(); 
       if(std::equal(sim_fraction_min3.begin() + 1, sim_fraction_min3.end(), sim_fraction_min3.begin())) sim_fraction_min3_MatchedIndex = -1;
       else sim_fraction_min3_MatchedIndex = std::max_element(sim_fraction_min3.begin(),sim_fraction_min3.end()) - sim_fraction_min3.begin();    
       if(std::equal(sim_rechit_diff.begin() + 1, sim_rechit_diff.end(), sim_rechit_diff.begin())) sim_rechit_diff_MatchedIndex = -1;
       else sim_rechit_diff_MatchedIndex = std::max_element(sim_rechit_diff.begin(),sim_rechit_diff.end()) - sim_rechit_diff.begin();  
       if(std::equal(sim_rechit_fraction.begin() + 1, sim_rechit_fraction.end(), sim_rechit_fraction.begin())) sim_rechit_fraction_MatchedIndex = -1;
       else sim_rechit_fraction_MatchedIndex = std::max_element(sim_rechit_fraction.begin(),sim_rechit_fraction.end()) - sim_rechit_fraction.begin(); 
       if(std::equal(global_sim_rechit_fraction.begin() + 1, global_sim_rechit_fraction.end(), global_sim_rechit_fraction.begin())) global_sim_rechit_fraction_MatchedIndex = -1;
       else global_sim_rechit_fraction_MatchedIndex = std::max_element(global_sim_rechit_fraction.begin(),global_sim_rechit_fraction.end()) - global_sim_rechit_fraction.begin();    

       tree->Fill();       
   }

}

void PFClusterDumper::beginJob()
{

}

void PFClusterDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
std::vector<float> PFClusterDumper::getScores(const std::vector<std::pair<DetId, float> >*hits_and_energies_Cluster, const std::vector<std::pair<DetId, float> > *hits_and_energies_CaloPart)
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

std::vector<std::pair<DetId, float> >* PFClusterDumper::getHitsAndEnergiesBC(reco::CaloCluster* iPFCluster, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_BC = new std::vector<std::pair<DetId, float> >;
    
    const std::vector<std::pair<DetId,float> > &hitsAndFractions = iPFCluster->hitsAndFractions();
    for(unsigned int i = 0; i < hitsAndFractions.size(); i++){
        if(hitsAndFractions.at(i).first.subdetId()==EcalBarrel){
           HitsAndEnergies_BC->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*(recHitsEB.product())->find(hitsAndFractions[i].first)).energy()));
        }else if(hitsAndFractions.at(i).first.subdetId()==EcalEndcap){
           HitsAndEnergies_BC->push_back(make_pair(hitsAndFractions.at(i).first,hitsAndFractions.at(i).second*(*(recHitsEE.product())->find(hitsAndFractions[i].first)).energy()));
        }
    }

    return HitsAndEnergies_BC;
}


std::vector<std::pair<DetId, float> >* PFClusterDumper::getHitsAndEnergiesCaloPart(CaloParticle* iCaloParticle)
{
    std::vector<std::pair<DetId, float> >* HitsAndEnergies_CaloPart = new std::vector<std::pair<DetId, float> >;
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
         HitsAndEnergies_CaloPart->push_back(make_pair(hit.first,hit.second));

    return HitsAndEnergies_CaloPart;
}


GlobalPoint PFClusterDumper::calculateAndSetPositionActual(const std::vector<std::pair<DetId, float> > *hits_and_energies_CP, double _param_T0_EB, double _param_T0_EE, double _param_T0_ES, double _param_W0, double _param_X0, double _minAllowedNorm, bool useES)
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

float PFClusterDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFClusterDumper);
