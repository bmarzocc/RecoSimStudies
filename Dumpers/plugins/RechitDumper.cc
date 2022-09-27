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

#include "DataFormats/Math/interface/libminifloat.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoSimStudies/Dumpers/plugins/RechitDumper.h"

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
RechitDumper::RechitDumper(const edm::ParameterSet& iConfig):
   caloTopologyToken_(esConsumes()),
   caloGeometryToken_(esConsumes()),
   channelStatusToken_(esConsumes())
{
   usesResource(TFileService::kSharedResource);

   pileupSummaryToken_            = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   ecalEBRechitToken_             = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalEBRechitCollection"));
   ecalEERechitToken_             = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalEERechitCollection"));
   ecalESRecHitToken_             = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ecalESRechitCollection"));
   hcalHBHERecHitToken_           = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("hcalEBHERechitCollection"));
   
   nBits_                         = iConfig.getParameter<int>("nBits");
   doCompression_                 = iConfig.getParameter<bool>("doCompression");
   isMC_                          = iConfig.getParameter<bool>("isMC");
   saveEB_                        = iConfig.getParameter<bool>("saveEB");
   saveEE_                        = iConfig.getParameter<bool>("saveEE");
  
   if(nBits_>23 && doCompression_){
      cout << "WARNING: float compression bits > 23 ---> Using 23 (i.e. no compression) instead!" << endl;
      nBits_=23;
   }

   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   setTree(tree);
}

RechitDumper::~RechitDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void RechitDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
   //calo geometry
   edm::ESHandle<CaloGeometry> pCaloGeometry = iSetup.getHandle(caloGeometryToken_);
   geometry = pCaloGeometry.product();

   ebGeometry = geometry->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
   eeGeometry = geometry->getSubdetectorGeometry(DetId::Ecal, EcalEndcap); 
   hbGeometry = geometry->getSubdetectorGeometry(DetId::Hcal, HcalBarrel); 
   heGeometry = geometry->getSubdetectorGeometry(DetId::Hcal, HcalEndcap); 
   esGeometry = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower); 

   //calo topology
   edm::ESHandle<CaloTopology> pCaloTopology = iSetup.getHandle(caloTopologyToken_); 
   topology = pCaloTopology.product();

   // channel status
   edm::ESHandle<EcalChannelStatus> pChannelStatus = iSetup.getHandle(channelStatusToken_); 
   chStatus = pChannelStatus.product();
   
   //MC-only info and collections
   truePU=-1.;
   obsPU=-1.;
   edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   
   if(isMC_){    
    ev.getByToken(pileupSummaryToken_, PupInfo);
    if (PupInfo.isValid()) 
    {
       for(auto &pu : *PupInfo){
           if(pu.getBunchCrossing() == 0 ){
              truePU = pu.getTrueNumInteractions();
              obsPU = pu.getPU_NumInteractions();
              break;
           } 
       } 
    }else{
       std::cerr << "Analyze --> PupInfo not found" << std::endl;
    }
   }

   //Other collections
   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }
    
   edm::Handle<EcalRecHitCollection> ecalRecHitsEB;
   ev.getByToken(ecalEBRechitToken_, ecalRecHitsEB);
   if (!ecalRecHitsEB.isValid()) {
       std::cerr << "Analyze --> ecalRecHitsEB not found" << std::endl; 
       return;
   }
   const EcalRecHitCollection* ecalEBHits = ecalRecHitsEB.product(); 
  
   edm::Handle<EcalRecHitCollection> ecalRecHitsEE;
   ev.getByToken(ecalEERechitToken_, ecalRecHitsEE);
   if (!ecalRecHitsEE.isValid()) {
       std::cerr << "Analyze --> ecalRecHitsEE not found" << std::endl; 
       return;
   } 
   const EcalRecHitCollection* ecalEEHits = ecalRecHitsEE.product();  
   
   edm::Handle<EcalRecHitCollection> ecalRecHitsES;
   ev.getByToken(ecalESRecHitToken_, ecalRecHitsES);
   if (!ecalRecHitsES.isValid()) {
       std::cerr << "Analyze --> ecalRecHitsES not found" << std::endl; 
       return;
   } 
   const EcalRecHitCollection* ecalESHits = ecalRecHitsES.product(); 
   
   edm::Handle<HBHERecHitCollection> hcalRecHitsHBHE;
   ev.getByToken(hcalHBHERecHitToken_, hcalRecHitsHBHE);
   if (!hcalRecHitsHBHE.isValid()) {
       std::cerr << "Analyze --> hcalRecHitsHBHE not found" << std::endl; 
       return;
   } 
   const HBHERecHitCollection* hcalHBHEHits = hcalRecHitsHBHE.product(); 
   
   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   nVtx = vertices->size();
   rho = *(rhos.product());

   std::vector<uint32_t> hcalRecHit_rawId;
   std::vector<float> hcalRecHit_energy;
   std::vector<float> hcalRecHit_eta; 
   std::vector<float> hcalRecHit_phi;
   std::vector<int> hcalRecHit_ieta; 
   std::vector<int> hcalRecHit_iphi;
   std::vector<int> hcalRecHit_iz;
   std::vector<int> hcalRecHit_depth;
   std::vector<uint32_t> esRecHit_rawId;
   std::vector<float> esRecHit_energy;
   std::vector<float> esRecHit_eta; 
   std::vector<float> esRecHit_phi;
   std::vector<int> esRecHit_ix; 
   std::vector<int> esRecHit_iy;
   std::vector<int> esRecHit_iz;
   std::vector<int> esRecHit_strip;
   std::vector<int> esRecHit_plane;

   GlobalPoint notFound(0, 0, 0);
   GlobalPoint ecalCell;
   GlobalPoint hcalCell;
   GlobalPoint esCell; 

   //Save ecalRechits 
   if(saveEB_){

      for(const auto& iRechit : *ecalEBHits){
          
          DetId rechit_id(iRechit.detid());

          uint32_t rawId = rechit_id.rawId();
          ecalRecHit_rawId.push_back(rawId);

          int status = (*chStatus->getMap().find(rechit_id.rawId())).getStatusCode();
          ecalRecHit_chStatus.push_back(status);

          ecalCell = ebGeometry->getGeometry(rechit_id)->getPosition(); 
          ecalRecHit_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          ecalRecHit_eta.push_back(reduceFloat(ecalCell.eta(),nBits_));  
          ecalRecHit_phi.push_back(reduceFloat(ecalCell.phi(),nBits_)); 

          EBDetId eb_id(rechit_id);  
          ecalRecHit_ieta.push_back(eb_id.ieta());  
          ecalRecHit_iphi.push_back(eb_id.iphi());  
          ecalRecHit_iz.push_back(0);     
 
          hcalRecHit_rawId.clear();
          hcalRecHit_energy.clear();
          hcalRecHit_eta.clear(); 
          hcalRecHit_phi.clear();
          hcalRecHit_ieta.clear(); 
          hcalRecHit_iphi.clear();
          hcalRecHit_iz.clear();
          hcalRecHit_depth.clear();

          for (auto hcalRechit: *hcalHBHEHits){
 
              HcalDetId hcal_id(hcalRechit.id());

              hcalCell = hbGeometry->getGeometry(hcal_id)->getPosition();
              double dEta = hbGeometry->deltaEta(hcal_id);
              double dPhi = hbGeometry->deltaPhi(hcal_id);
 
              if(hcalCell==notFound){
                 hcalCell = heGeometry->getGeometry(hcal_id)->getPosition();
                 dEta = heGeometry->deltaEta(hcal_id);
                 dPhi = heGeometry->deltaPhi(hcal_id);
              }
   
              if(fabs(ecalCell.eta()-hcalCell.eta()) > 1.1 * (ebGeometry->deltaEta(rechit_id) + dEta)) continue;
              if(fabs(deltaPhi(double(ecalCell.phi()),double(hcalCell.phi()))) > 1.1 * (ebGeometry->deltaPhi(rechit_id) + dPhi)) continue;   
              
              int iz=-99;
              if(hcal_id.zside()==0) iz=0;
              if(hcal_id.zside()<0) iz=-1;
              if(hcal_id.zside()>0) iz=1; 

              hcalRecHit_rawId.push_back(hcal_id.rawId());
              hcalRecHit_energy.push_back(hcalRechit.energy());
              hcalRecHit_eta.push_back(hcalCell.eta()); 
              hcalRecHit_phi.push_back(hcalCell.phi());
              hcalRecHit_ieta.push_back(hcal_id.ieta()); 
              hcalRecHit_iphi.push_back(hcal_id.iphi());
              hcalRecHit_iz.push_back(iz); 
              hcalRecHit_depth.push_back(hcal_id.depth()); 
          }
           
          if(hcalRecHit_rawId.size()==0) std::cout << "WARNING: no matching HEHB Rechit for EB ecalHit: (" << ecalCell.eta() << "," << ecalCell.phi() << ") - ( " << eb_id.ieta() << "," << eb_id.iphi() << ")" << std::endl;
 
          matchedHcalRecHit_rawId.push_back(hcalRecHit_rawId);
          matchedHcalRecHit_energy.push_back(hcalRecHit_energy);
          matchedHcalRecHit_eta.push_back(hcalRecHit_eta);
          matchedHcalRecHit_phi.push_back(hcalRecHit_phi);
          matchedHcalRecHit_ieta.push_back(hcalRecHit_ieta); 
          matchedHcalRecHit_iphi.push_back(hcalRecHit_iphi);
          matchedHcalRecHit_iz.push_back(hcalRecHit_iz); 
          matchedHcalRecHit_depth.push_back(hcalRecHit_depth); 
      }
   }

   if(saveEE_){
  
      for(const auto& iRechit : *ecalEEHits){

          DetId rechit_id(iRechit.detid());

          uint32_t rawId = rechit_id.rawId();
          ecalRecHit_rawId.push_back(rawId);

          int status = (*chStatus->getMap().find(rechit_id.rawId())).getStatusCode();
          ecalRecHit_chStatus.push_back(status);

          ecalCell = eeGeometry->getGeometry(rechit_id)->getPosition();
          ecalRecHit_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          ecalRecHit_eta.push_back(reduceFloat(ecalCell.eta(),nBits_));  
          ecalRecHit_phi.push_back(reduceFloat(ecalCell.phi(),nBits_)); 

          int iz=-99;
          EEDetId ee_id(rechit_id);  
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1; 
          ecalRecHit_ieta.push_back(ee_id.ix());  
          ecalRecHit_iphi.push_back(ee_id.iy());  
          ecalRecHit_iz.push_back(iz);  

          hcalRecHit_rawId.clear();
          hcalRecHit_energy.clear();
          hcalRecHit_eta.clear(); 
          hcalRecHit_phi.clear();
          hcalRecHit_ieta.clear(); 
          hcalRecHit_iphi.clear();
          hcalRecHit_iz.clear();
          hcalRecHit_depth.clear();

          for (auto hcalRechit: *hcalHBHEHits){
 
              HcalDetId hcal_id(hcalRechit.id());

              hcalCell = heGeometry->getGeometry(hcal_id)->getPosition();
              double dEta = heGeometry->deltaEta(hcal_id);
              double dPhi = heGeometry->deltaPhi(hcal_id);
 
              if(hcalCell==notFound){
                 hcalCell = hbGeometry->getGeometry(hcal_id)->getPosition();
                 dEta = hbGeometry->deltaEta(hcal_id);
                 dPhi = hbGeometry->deltaPhi(hcal_id);
              }
              
              if(hcal_id.zside()!=ee_id.zside()) continue;
              if(fabs(ecalCell.eta()-hcalCell.eta()) > 1.1 * (eeGeometry->deltaEta(rechit_id) + dEta)) continue;
              if(fabs(deltaPhi(double(ecalCell.phi()),double(hcalCell.phi()))) > 1.1 * (eeGeometry->deltaPhi(rechit_id) + dPhi)) continue;       

              int iz=-99;
              if(hcal_id.zside()<0) iz=-1;
              if(hcal_id.zside()>0) iz=1; 

              hcalRecHit_rawId.push_back(hcal_id.rawId());
              hcalRecHit_energy.push_back(hcalRechit.energy());
              hcalRecHit_eta.push_back(hcalCell.eta()); 
              hcalRecHit_phi.push_back(hcalCell.phi());
              hcalRecHit_ieta.push_back(hcal_id.ieta()); 
              hcalRecHit_iphi.push_back(hcal_id.iphi());
              hcalRecHit_iz.push_back(iz); 
              hcalRecHit_depth.push_back(hcal_id.depth());
          }
 
          if(hcalRecHit_rawId.size()==0) std::cout << "WARNING: no matching HEHB Rechit for EE ecalHit: (" << ecalCell.eta() << "," << ecalCell.phi() << ") - ( " << ee_id.ix() << "," << ee_id.iy() << "," << ee_id.zside() << ")" << std::endl;

          matchedHcalRecHit_rawId.push_back(hcalRecHit_rawId);
          matchedHcalRecHit_energy.push_back(hcalRecHit_energy);
          matchedHcalRecHit_eta.push_back(hcalRecHit_eta);
          matchedHcalRecHit_phi.push_back(hcalRecHit_phi);
          matchedHcalRecHit_ieta.push_back(hcalRecHit_ieta); 
          matchedHcalRecHit_iphi.push_back(hcalRecHit_iphi);
          matchedHcalRecHit_iz.push_back(hcalRecHit_iz); 
          matchedHcalRecHit_depth.push_back(hcalRecHit_depth);  
  
          esRecHit_rawId.clear();
          esRecHit_energy.clear();
          esRecHit_eta.clear(); 
          esRecHit_phi.clear();
          esRecHit_ix.clear(); 
          esRecHit_iy.clear();
          esRecHit_iz.clear();
          esRecHit_strip.clear();
          esRecHit_plane.clear();

          for (auto esRechit: *ecalESHits){

              ESDetId es_id(esRechit.id()); 
              esCell = esGeometry->getGeometry(es_id)->getPosition();

              if(es_id.zside()!=ee_id.zside()) continue;
              if(fabs(ecalCell.eta()-esCell.eta()) > 1.1 * (eeGeometry->deltaEta(rechit_id) + esGeometry->deltaEta(es_id))) continue;
              if(fabs(deltaPhi(double(ecalCell.phi()),double(esCell.phi()))) > 1.1 * (eeGeometry->deltaPhi(rechit_id) + esGeometry->deltaPhi(es_id))) continue;     
             
              esRecHit_rawId.push_back(es_id.rawId());
              esRecHit_energy.push_back(esRechit.energy());
              esRecHit_eta.push_back(esCell.eta()); 
              esRecHit_phi.push_back(esCell.phi());
              esRecHit_ix.push_back(es_id.six()); 
              esRecHit_iy.push_back(es_id.siy());
              esRecHit_iz.push_back(es_id.zside());
              esRecHit_strip.push_back(es_id.strip()); 
              esRecHit_plane.push_back(es_id.plane());  
          }
 
          if(hcalRecHit_rawId.size()==0) std::cout << "WARNING: no matching ES Rechit for EE ecalHit: (" << ecalCell.eta() << "," << ecalCell.phi() << ") - ( " << ee_id.ix() << "," << ee_id.iy() << "," << ee_id.zside() << ")" << std::endl;

          matchedESRecHit_rawId.push_back(esRecHit_rawId);
          matchedESRecHit_energy.push_back(esRecHit_energy);
          matchedESRecHit_eta.push_back(esRecHit_eta);
          matchedESRecHit_phi.push_back(esRecHit_phi);
          matchedESRecHit_ix.push_back(esRecHit_ix); 
          matchedESRecHit_iy.push_back(esRecHit_iy);
          matchedESRecHit_iz.push_back(esRecHit_iz); 
          matchedESRecHit_strip.push_back(esRecHit_strip); 
          matchedESRecHit_plane.push_back(esRecHit_plane);  
      }   
   }  

   //fill tree for each event
   tree->Fill();
}

void RechitDumper::beginJob()
{

}

void RechitDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void RechitDumper::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("rho", &rho, "rho/F"); 
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   if(isMC_){ 
    tree->Branch("truePU", &truePU, "truePU/D");
    tree->Branch("obsPU", &obsPU, "obsPU/D");
   }   
   tree->Branch("ecalRecHit_rawId","std::vector<uint32_t>",&ecalRecHit_rawId);
   tree->Branch("ecalRecHit_chStatus","std::vector<int>",&ecalRecHit_chStatus);
   tree->Branch("ecalRecHit_energy","std::vector<float>",&ecalRecHit_energy);
   tree->Branch("ecalRecHit_eta","std::vector<float>",&ecalRecHit_eta); 
   tree->Branch("ecalRecHit_phi","std::vector<float>",&ecalRecHit_phi);
   if(!saveEB_){
      tree->Branch("ecalRecHit_ix","std::vector<int>",&ecalRecHit_ieta); 
      tree->Branch("ecalRecHit_iy","std::vector<int>",&ecalRecHit_iphi);
   }else{  
      tree->Branch("ecalRecHit_ieta","std::vector<int>",&ecalRecHit_ieta); 
      tree->Branch("ecalRecHit_iphi","std::vector<int>",&ecalRecHit_iphi);
   }
   tree->Branch("ecalRecHit_iz","std::vector<int>",&ecalRecHit_iz); 
   tree->Branch("matchedHcalRecHit_rawId","std::vector<std::vector<uint32_t> >",&matchedHcalRecHit_rawId);
   tree->Branch("matchedHcalRecHit_energy","std::vector<std::vector<float> >",&matchedHcalRecHit_energy);
   tree->Branch("matchedHcalRecHit_eta","std::vector<std::vector<float> >",&matchedHcalRecHit_eta); 
   tree->Branch("matchedHcalRecHit_phi","std::vector<std::vector<float> >",&matchedHcalRecHit_phi);
   tree->Branch("matchedHcalRecHit_ieta","std::vector<std::vector<int> >",&matchedHcalRecHit_ieta); 
   tree->Branch("matchedHcalRecHit_iphi","std::vector<std::vector<int> >",&matchedHcalRecHit_iphi);
   tree->Branch("matchedHcalRecHit_iz","std::vector<std::vector<int> >",&matchedHcalRecHit_iz); 
   tree->Branch("matchedHcalRecHit_depth","std::vector<std::vector<int> >",&matchedHcalRecHit_depth);  
   tree->Branch("matchedESRecHit_rawId","std::vector<std::vector<uint32_t> >",&matchedESRecHit_rawId);
   tree->Branch("matchedESRecHit_energy","std::vector<std::vector<float> >",&matchedESRecHit_energy);
   tree->Branch("matchedESRecHit_eta","std::vector<std::vector<float> >",&matchedESRecHit_eta); 
   tree->Branch("matchedESRecHit_phi","std::vector<std::vector<float> >",&matchedESRecHit_phi);
   tree->Branch("matchedESRecHit_ix","std::vector<std::vector<int> >",&matchedESRecHit_ix); 
   tree->Branch("matchedESRecHit_iy","std::vector<std::vector<int> >",&matchedESRecHit_iy);
   tree->Branch("matchedESRecHit_iz","std::vector<std::vector<int> >",&matchedESRecHit_iz);      
   tree->Branch("matchedESRecHit_strip","std::vector<std::vector<int> >",&matchedESRecHit_strip);   
   tree->Branch("matchedESRecHit_plane","std::vector<std::vector<int> >",&matchedESRecHit_plane);   
}

void RechitDumper::clearVectors()
{
   ecalRecHit_rawId.clear();
   ecalRecHit_chStatus.clear();
   ecalRecHit_energy.clear();
   ecalRecHit_eta.clear(); 
   ecalRecHit_phi.clear();
   ecalRecHit_ieta.clear(); 
   ecalRecHit_iphi.clear();
   ecalRecHit_iz.clear(); 
   matchedHcalRecHit_rawId.clear();
   matchedHcalRecHit_energy.clear();
   matchedHcalRecHit_eta.clear(); 
   matchedHcalRecHit_phi.clear();
   matchedHcalRecHit_ieta.clear(); 
   matchedHcalRecHit_iphi.clear();
   matchedHcalRecHit_iz.clear(); 
   matchedHcalRecHit_depth.clear(); 
   matchedESRecHit_rawId.clear();
   matchedESRecHit_energy.clear();
   matchedESRecHit_eta.clear(); 
   matchedESRecHit_phi.clear();
   matchedESRecHit_ix.clear(); 
   matchedESRecHit_iy.clear();
   matchedESRecHit_iz.clear(); 
   matchedESRecHit_strip.clear(); 
   matchedESRecHit_plane.clear(); 
}

float RechitDumper::reduceFloat(float val, int bits)
{
    if(!doCompression_) return val;
    else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RechitDumper);

