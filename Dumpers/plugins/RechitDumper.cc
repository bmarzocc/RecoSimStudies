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
 
   //calo topology
   edm::ESHandle<CaloTopology> pCaloTopology = iSetup.getHandle(caloTopologyToken_); 
   topology = pCaloTopology.product();

   
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

   GlobalPoint cell;

   //Save ecalRechits 
   if(saveEB_){
      for(const auto& iRechit : *ecalEBHits){
          
          DetId rechit_id(iRechit.detid());

          uint32_t rawId = rechit_id.rawId();
          ecalRecHit_rawId.push_back(rawId);

          int status = (*chStatus->getMap().find(rechit_id.rawId())).getStatusCode();
          ecalRecHit_chStatus.push_back(status);

          cell = geometry->getPosition(rechit_id); 
          ecalRecHit_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          ecalRecHit_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          ecalRecHit_phi.push_back(reduceFloat(cell.phi(),nBits_)); 

          EBDetId eb_id(rechit_id);  
          ecalRecHit_ieta.push_back(eb_id.ieta());  
          ecalRecHit_iphi.push_back(eb_id.iphi());  
          ecalRecHit_iz.push_back(0);     
      
      }
   }

   if(saveEE_){
      for(const auto& iRechit : *ecalEEHits){

          DetId rechit_id(iRechit.detid());

          uint32_t rawId = rechit_id.rawId();
          ecalRecHit_rawId.push_back(rawId);

          int status = (*chStatus->getMap().find(rechit_id.rawId())).getStatusCode();
          ecalRecHit_chStatus.push_back(status);

          cell = geometry->getPosition(rechit_id); 
          ecalRecHit_energy.push_back(reduceFloat(iRechit.energy(),nBits_));    
          ecalRecHit_eta.push_back(reduceFloat(cell.eta(),nBits_));  
          ecalRecHit_phi.push_back(reduceFloat(cell.phi(),nBits_)); 

          int iz=-99;
          EEDetId ee_id(rechit_id);  
          if(ee_id.zside()<0) iz=-1;
          if(ee_id.zside()>0) iz=1; 
          ecalRecHit_ieta.push_back(ee_id.ix());  
          ecalRecHit_iphi.push_back(ee_id.iy());  
          ecalRecHit_iz.push_back(iz);    
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
   tree->Branch("ecalRecHit_ieta","std::vector<int>",&ecalRecHit_ieta); 
   tree->Branch("ecalRecHit_iphi","std::vector<int>",&ecalRecHit_iphi);
   tree->Branch("ecalRecHit_iz","std::vector<int>",&ecalRecHit_iz); 
   tree->Branch("ecalRecHit_relatedHcalRecHit_energy","std::vector<float>",&ecalRecHit_relatedHcalRecHit_energy);
   tree->Branch("ecalRecHit_relatedHcalRecHit_eta","std::vector<float>",&ecalRecHit_relatedHcalRecHit_eta); 
   tree->Branch("ecalRecHit_relatedHcalRecHit_phi","std::vector<float>",&ecalRecHit_relatedHcalRecHit_phi);
   tree->Branch("ecalRecHit_relatedHcalRecHit_ieta","std::vector<int>",&ecalRecHit_relatedHcalRecHit_ieta); 
   tree->Branch("ecalRecHit_relatedHcalRecHit_iphi","std::vector<int>",&ecalRecHit_relatedHcalRecHit_iphi);
   tree->Branch("ecalRecHit_relatedHcalRecHit_iz","std::vector<int>",&ecalRecHit_relatedHcalRecHit_iz);      
   
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
   ecalRecHit_relatedHcalRecHit_energy.clear();
   ecalRecHit_relatedHcalRecHit_eta.clear(); 
   ecalRecHit_relatedHcalRecHit_phi.clear();
   ecalRecHit_relatedHcalRecHit_ieta.clear(); 
   ecalRecHit_relatedHcalRecHit_iphi.clear();
   ecalRecHit_relatedHcalRecHit_iz.clear(); 
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

