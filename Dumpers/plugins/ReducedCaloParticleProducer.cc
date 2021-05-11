#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace edm;
using namespace std;

class ReducedCaloParticleProducer : public EDProducer
{

 public:
    ReducedCaloParticleProducer( const ParameterSet & );
 private:
    void produce( Event &, const EventSetup & ) override;
    edm::EDGetTokenT<View<CaloParticle> > caloPartToken_;
    bool saveSignal_;
    bool savePU_;
    bool saveOOTPU_;
    std::map<DetId, float> HitsAndEnergiesMap_;
    std::map<int,vector<int>> SimClusterMap_;
    int iSC_;
};

ReducedCaloParticleProducer::ReducedCaloParticleProducer( const ParameterSet &iConfig ) : 
   caloPartToken_( consumes<View<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleTag")) )
{
   saveSignal_ = iConfig.getParameter<bool>("saveSignal");
   savePU_ = iConfig.getParameter<bool>("savePU");
   saveOOTPU_ = iConfig.getParameter<bool>("saveOOTPU"); 
   produces<vector<CaloParticle>>();
   produces<vector<SimCluster>>();
}

void ReducedCaloParticleProducer::produce( Event &evt, const EventSetup & )
{
   edm::Handle<View<CaloParticle> > caloParticles;
   evt.getByToken(caloPartToken_,caloParticles);
   if (!caloParticles.isValid()) {
       std::cerr << "ReducedCaloParticleProducer --> caloParticles not found" << std::endl; 
       return;
   }

   const std::vector<edm::Ptr<CaloParticle> > &caloPointers = caloParticles->ptrs();
   
   unique_ptr<vector<CaloParticle> > outCaloParticle( new vector<CaloParticle> );   
   unique_ptr<vector<SimCluster> > outSimCluster( new vector<SimCluster> );  

   HitsAndEnergiesMap_.clear();
   SimClusterMap_.clear();
   iSC_=0;

   for( unsigned int i = 0 ; i < caloPointers.size() ; i++ ) {
        edm::Ptr<CaloParticle> calo = caloPointers[i];

        if( saveSignal_) {

            if( calo->g4Tracks()[0].eventId().event()==0 && calo->g4Tracks()[0].eventId().bunchCrossing()==0 ) { 

                const auto& simClusters = calo->simClusters();
                for( unsigned int j = 0; j < simClusters.size() ; j++ ) {
                     outSimCluster->push_back( *simClusters[j] );
                     SimClusterMap_[i].push_back(iSC_);
                     iSC_++; 
                }
            } 

        } else {

            if(calo->g4Tracks()[0].eventId().bunchCrossing()!=0 && savePU_) continue; 
            if(calo->g4Tracks()[0].eventId().bunchCrossing()==0 && saveOOTPU_) continue;  

            const auto& simClusters = calo->simClusters();
            for( unsigned int iSC = 0; iSC < simClusters.size() ; iSC++ ) {
                 auto simCluster = simClusters[iSC];  
                 auto hits_and_energies = simCluster->hits_and_energies();
                 for( unsigned int i = 0; i < hits_and_energies.size(); i++ ) {
                      if ( HitsAndEnergiesMap_.find(DetId(hits_and_energies[i].first)) == HitsAndEnergiesMap_.end() ) {
                           HitsAndEnergiesMap_[DetId(hits_and_energies[i].first)]=hits_and_energies[i].second;      
                      }else{
                           HitsAndEnergiesMap_[DetId(hits_and_energies[i].first)]=HitsAndEnergiesMap_[DetId(hits_and_energies[i].first)]+hits_and_energies[i].second; 
                      }
                 }
            }
        }
   }

   if(  saveSignal_ ) {

        auto scHandle = evt.put( std::move( outSimCluster ) ); 
        for( auto const &map : SimClusterMap_ ) {
             CaloParticle reducedCalo = *caloPointers[map.first];
             reducedCalo.clearSimClusters();
             for( unsigned int i=0; i<map.second.size(); i++ ) {
                  edm::Ref<vector<SimCluster>> ref(scHandle, map.second.at(i));
                  reducedCalo.addSimCluster(ref);
             } 
             outCaloParticle->push_back( reducedCalo );
        }
        evt.put( std::move( outCaloParticle ) );

   } else {

      CaloParticle reducedCalo;
      SimCluster reducedCluster;
      for( auto const &hit_and_energy : HitsAndEnergiesMap_ ) {
           reducedCluster.addRecHitAndFraction(hit_and_energy.first, hit_and_energy.second);
           //reducedCluster.addHitEnergy(hit_and_energy.second);
      }

      outSimCluster->push_back(reducedCluster);    
      auto scHandle = evt.put( std::move( outSimCluster ) ); 
      edm::Ref<vector<SimCluster>> ref(scHandle, 0);
      reducedCalo.addSimCluster(ref);
      outCaloParticle->push_back( reducedCalo );
      evt.put( std::move( outCaloParticle ) );
   }
}

DEFINE_FWK_MODULE( ReducedCaloParticleProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
