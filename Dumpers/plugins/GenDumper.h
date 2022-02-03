#ifndef RecoSimStudies_Dumpers_GenDumper_H
#define RecoSimStudies_Dumpers_GenDumper_H

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

class GenDumper : public edm::one::EDAnalyzer<edm::one::SharedResources>
{
      public:
         explicit GenDumper(const edm::ParameterSet&);
	 ~GenDumper();
  
      private:
	 void beginJob() override;
	 void analyze(const edm::Event&, const edm::EventSetup&) override;
         void endJob() override;
        
      // ----------additional functions-------------------
      void setTree(TTree* tree);
     
      // ----------collection tokens-------------------
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::Service<TFileService> iFile;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      int runId;
      int lumiId;
      int eventId;

      int genParticle_size; 
      std::vector<int> genParticle_pdgId;
      std::vector<int> genParticle_status; 
      std::vector<float> genParticle_energy;
      std::vector<float> genParticle_pt;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;  
 
      int genMother_size; 
      std::vector<int> genMother_pdgId;
      std::vector<int> genMother_status; 
      std::vector<float> genMother_energy;
      std::vector<float> genMother_pt;
      std::vector<float> genMother_eta;
      std::vector<float> genMother_phi;

};

#endif
