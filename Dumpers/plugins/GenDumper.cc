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

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoSimStudies/Dumpers/plugins/GenDumper.h"

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
GenDumper::GenDumper(const edm::ParameterSet& iConfig)
{
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   
   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   setTree(tree);
}

GenDumper::~GenDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void GenDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   ev.getByToken(genToken_,genParticles);
   if (!genParticles.isValid()) {
       std::cerr << "Analyze --> genParticles not found" << std::endl; 
       return;
   }

   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   
   genParticle_pdgId.clear();
   genParticle_status.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();

   genMother_pdgId.clear();
   genMother_status.clear();
   genMother_energy.clear();
   genMother_pt.clear();
   genMother_eta.clear();
   genMother_phi.clear();
   
   int genPart_index = 0;
   int genMother_index = 0;

   for(const auto& iGen : *(genParticles.product()))
   {
       
       genParticle_pdgId.push_back(iGen.pdgId()); 
       genParticle_status.push_back(iGen.status()); 
       genParticle_energy.push_back(iGen.energy()); 
       genParticle_pt.push_back(iGen.pt());
       genParticle_eta.push_back(iGen.eta());
       genParticle_phi.push_back(iGen.phi());
       genPart_index++;

       if(iGen.numberOfMothers()!=0) continue;

       genMother_pdgId.push_back(iGen.pdgId()); 
       genMother_status.push_back(iGen.status()); 
       genMother_energy.push_back(iGen.energy()); 
       genMother_pt.push_back(iGen.pt());
       genMother_eta.push_back(iGen.eta());
       genMother_phi.push_back(iGen.phi());
       genMother_index++;
       
   } 

   genParticle_size = genPart_index; 
   genMother_size = genMother_index; 
   
   //fill tree for each event
   tree->Fill();
}

void GenDumper::beginJob()
{

}

void GenDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void GenDumper::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("genParticle_size", &genParticle_size, "genParticle_size/I");  
   tree->Branch("genParticle_pdgId","std::vector<int>",&genParticle_pdgId);
   tree->Branch("genParticle_status","std::vector<int>",&genParticle_status); 
   tree->Branch("genParticle_energy","std::vector<float>",&genParticle_energy);
   tree->Branch("genParticle_pt","std::vector<float>",&genParticle_pt);
   tree->Branch("genParticle_eta","std::vector<float>",&genParticle_eta);
   tree->Branch("genParticle_phi","std::vector<float>",&genParticle_phi);
   tree->Branch("genMother_size", &genMother_size, "genMother_size/I");  
   tree->Branch("genMother_pdgId","std::vector<int>",&genMother_pdgId);
   tree->Branch("genMother_status","std::vector<int>",&genMother_status); 
   tree->Branch("genMother_energy","std::vector<float>",&genMother_energy);
   tree->Branch("genMother_pt","std::vector<float>",&genMother_pt);
   tree->Branch("genMother_eta","std::vector<float>",&genMother_eta);
   tree->Branch("genMother_phi","std::vector<float>",&genMother_phi);
}   


///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenDumper);

