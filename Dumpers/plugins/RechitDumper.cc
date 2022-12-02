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

   matchingScale_                 = iConfig.getParameter<double>("matchingScale"); 
   neighborXtalsMatrix_           = iConfig.getParameter<std::vector<int> >("neighborXtalsMatrix");
   deadXtalsEB_                   = iConfig.getParameter<std::vector<uint32_t> >("deadXtalsEB");
   deadXtalsEE_                   = iConfig.getParameter<std::vector<uint32_t> >("deadXtalsEE");
  
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
   CaloCellGeometry::CornersVec ecalCornerVec;
   CaloCellGeometry::CornersVec hcalCornerVec;
   CaloCellGeometry::CornersVec esCornerVec;
   std::vector<Point> ecalPolygon;
   std::vector<Point> hcalPolygon;
   std::vector<Point> esPolygon; 
   Point ecalBar;
   Point hcalBar;
   Point esBar;
   
   clearVectors();
   saveXtalsEB_.clear();
   saveXtalsEE_.clear();
  
   if(saveEB_){

      //Find the xtals close to the deadXtals
      for(auto& xtal : deadXtalsEB_){
          DetId deadId(xtal);
          std::vector<DetId> neighborXtals = topology->getSubdetectorTopology(DetId::Ecal, EcalBarrel)->getWindow(deadId,neighborXtalsMatrix_.at(0), neighborXtalsMatrix_.at(1));  
          for(auto&  id: neighborXtals)
              saveXtalsEB_.push_back(id.rawId());      
      }

      //Remove the duplicates
      saveXtalsEB_.insert(saveXtalsEB_.end(),deadXtalsEB_.begin(),deadXtalsEB_.end());
      std::sort(saveXtalsEB_.begin(),saveXtalsEB_.end());
      saveXtalsEB_.erase(std::unique(saveXtalsEB_.begin(),saveXtalsEB_.end()),saveXtalsEB_.end());

      //std::cout << "saveXtalsEB_.size() = " << saveXtalsEB_.size() << std::endl; 

      for(const auto& iRechit : *ecalEBHits){
          
          DetId rechit_id(iRechit.detid());
          uint32_t rawId = rechit_id.rawId();

          if(std::find(saveXtalsEB_.begin(),saveXtalsEB_.end(), rawId) == saveXtalsEB_.end() && saveXtalsEB_.size()!=0) continue;

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

          ecalCornerVec = ebGeometry->getGeometry(rechit_id)->getCorners(); 

          //ECAL crystal rear face
          ecalPolygon = makePolygon(ecalCornerVec,std::vector<int>{4,5,6,7});
          ecalBar = computeBarycenter(ecalPolygon);   
          double dRmaxEcal = computeMaxDR(ecalPolygon, ecalBar);

          for (auto hcalRechit: *hcalHBHEHits){
 
              HcalDetId hcal_id(hcalRechit.id());

              hcalCell = hbGeometry->getGeometry(hcal_id)->getPosition();
              hcalCornerVec = hbGeometry->getGeometry(hcal_id)->getCorners(); 
              if(hcalCell==notFound){
                 hcalCell = heGeometry->getGeometry(hcal_id)->getPosition();
                 hcalCornerVec = heGeometry->getGeometry(hcal_id)->getCorners(); 
              }

              //HCAL crystal front face 
              hcalPolygon = makePolygon(hcalCornerVec,std::vector<int>{0,1,2,3});
              hcalBar = computeBarycenter(hcalPolygon);   
              double dRmaxHcal = computeMaxDR(hcalPolygon, hcalBar);

              if(deltaR(ecalBar.eta,ecalBar.phi,hcalBar.eta,hcalBar.phi)>matchingScale_*(dRmaxEcal+dRmaxHcal)) continue;  
              
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
  
      //Find the xtals close to the deadXtals
      for(auto& xtal : deadXtalsEE_){
          DetId deadId(xtal);
          std::vector<DetId> neighborXtals = topology->getSubdetectorTopology(DetId::Ecal, EcalEndcap)->getWindow(deadId,neighborXtalsMatrix_.at(0), neighborXtalsMatrix_.at(1));   
          for(auto&  id: neighborXtals)
              saveXtalsEE_.push_back(id.rawId());      
      }

      //Remove the duplicates
      saveXtalsEE_.insert(saveXtalsEE_.end(),deadXtalsEE_.begin(),deadXtalsEE_.end());
      std::sort(saveXtalsEE_.begin(),saveXtalsEE_.end());
      saveXtalsEE_.erase(std::unique(saveXtalsEE_.begin(),saveXtalsEE_.end()),saveXtalsEE_.end());

      //std::cout << "saveXtalsEE_.size() = " << saveXtalsEE_.size() << std::endl; 
 
      for(const auto& iRechit : *ecalEEHits){

          DetId rechit_id(iRechit.detid());
          uint32_t rawId = rechit_id.rawId();

          if(std::find(saveXtalsEE_.begin(), saveXtalsEE_.end(), rawId) == saveXtalsEE_.end() && saveXtalsEE_.size()!=0) continue;

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

          ecalCornerVec = eeGeometry->getGeometry(rechit_id)->getCorners(); 

          //ECAL crystal rear face
          ecalPolygon = makePolygon(ecalCornerVec,std::vector<int>{4,5,6,7});
          ecalBar = computeBarycenter(ecalPolygon); 
          double dRmaxEcal = computeMaxDR(ecalPolygon, ecalBar); 

          for (auto hcalRechit: *hcalHBHEHits){
 
              HcalDetId hcal_id(hcalRechit.id());
              if(hcal_id.zside()!=ee_id.zside()) continue;
 
              hcalCell = heGeometry->getGeometry(hcal_id)->getPosition();
              hcalCornerVec = heGeometry->getGeometry(hcal_id)->getCorners(); 
              if(hcalCell==notFound){
                 hcalCell = hbGeometry->getGeometry(hcal_id)->getPosition();
                 hcalCornerVec = hbGeometry->getGeometry(hcal_id)->getCorners(); 
              }

              //HCAL crystal front face 
              hcalPolygon = makePolygon(hcalCornerVec,std::vector<int>{0,1,2,3});
              hcalBar = computeBarycenter(hcalPolygon);   
              double dRmaxHcal = computeMaxDR(hcalPolygon, hcalBar); 
   
              if(deltaR(ecalBar.eta,ecalBar.phi,hcalBar.eta,hcalBar.phi)>matchingScale_*(dRmaxEcal+dRmaxHcal)) continue;  
              
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
 
          if(hcalRecHit_rawId.size()==0) std::cout << "WARNING: no matching HEHB Rechit for EE ecalHit: (" << (double)ecalCell.eta() << "," << (double)ecalCell.phi() << ") - ( " << (int)ee_id.ix() << "," << (int)ee_id.iy() << "," << (int)ee_id.zside() << ")" << std::endl;

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

          //ECAL crystal front face
          ecalPolygon = makePolygon(ecalCornerVec,std::vector<int>{0,1,2,3});
          ecalBar = computeBarycenter(ecalPolygon); 
          dRmaxEcal = computeMaxDR(ecalPolygon, ecalBar);

          for (auto esRechit: *ecalESHits){

              ESDetId es_id(esRechit.id()); 
              if(es_id.zside()!=ee_id.zside()) continue;

              esCell = esGeometry->getGeometry(es_id)->getPosition();
              esCornerVec = esGeometry->getGeometry(es_id)->getCorners(); 

              //ES strips rear face
              esPolygon = makePolygon(esCornerVec,std::vector<int>{4,5,6,7});
              esBar = computeBarycenter(esPolygon);     
              double dRmaxES = computeMaxDR(esPolygon, esBar); 
   
              if(deltaR(ecalBar.eta,ecalBar.phi,esBar.eta,esBar.phi)>matchingScale_*(dRmaxEcal+dRmaxES)) continue; 

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
 
          if(hcalRecHit_rawId.size()==0) std::cout << "WARNING: no matching ES Rechit for EE ecalHit: (" << (double)ecalCell.eta() << "," << (double)ecalCell.phi() << ") - ( " << (int)ee_id.ix() << "," << (int)ee_id.iy() << "," << (int)ee_id.zside() << ")" << std::endl;

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

std::vector<Point> RechitDumper::makePolygon(CaloCellGeometry::CornersVec& corners, std::vector<int> indices)
{
    std::vector<Point> polygon;
    polygon.resize(indices.size());
    for(unsigned int i=0; i<indices.size(); i++){
        int index = indices[i];
        polygon[i] = Point{corners[index].x(),corners[index].y(),corners[index].z(),corners[index].eta(),corners[index].phi()};     
    }
    return polygon;
}

Point RechitDumper::computeBarycenter(std::vector<Point>& polygon)
{
    double barX=0.;
    double barY=0.; 
    double barZ=0.; 
    double barEta=0.;
    double barPhi=0.; 
    for(unsigned int i=0; i<polygon.size(); i++){
       barX += polygon[i].x;
       barY += polygon[i].y;
       barZ += polygon[i].z;
       barEta += polygon[i].eta;
       barPhi += polygon[i].phi;
   }
   return Point{barX/polygon.size(),barY/polygon.size(),barZ/polygon.size(),barEta/polygon.size(),barPhi/polygon.size()}; 
}

double RechitDumper::computeMaxDR(std::vector<Point>& polygon, Point& barycenter)
{
   std::vector<double> dRs;
   for(unsigned int i=0; i<polygon.size(); i++)
       dRs.push_back(deltaR(polygon[i].eta,polygon[i].phi,barycenter.eta,barycenter.phi));

   return *max_element(dRs.begin(),dRs.end()); 
}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(RechitDumper);

