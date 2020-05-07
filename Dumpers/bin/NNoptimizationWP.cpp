#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/NNoptimizationWP.h"
#include "RecoSimStudies/Dumpers/interface/CruijffPdf.h"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include <algorithm> 
#include <iostream>
#include <utility>

#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooRealProxy.h"

using namespace std;

int main(int argc, char** argv)
{
    const edm::ParameterSet &process         = edm::readPSetsFrom( argv[1] )->getParameter<edm::ParameterSet>( "process" );
    const edm::ParameterSet &filesOpt        = process.getParameter<edm::ParameterSet>( "ioFilesOpt" );
    
    //Read config
    string inputFileEB_    = filesOpt.getParameter<string>( "inputFile_EB" );
    string inputFileEE_    = filesOpt.getParameter<string>( "inputFile_EE" );
    string inputFileMusEB_ = filesOpt.getParameter<string>( "inputFileMustache_EB" );
    string inputFileMusEE_ = filesOpt.getParameter<string>( "inputFileMustache_EE" );  

    TFile* inFileEB_       = TFile::Open(inputFileEB_.c_str());
    TFile* inFileEE_       = TFile::Open(inputFileEE_.c_str());
    TFile* inFileMusEB_    = TFile::Open(inputFileMusEB_.c_str());
    TFile* inFileMusEE_    = TFile::Open(inputFileMusEE_.c_str());

    TTree* inputTreeEB_    = (TTree*)inFileEB_->Get("resolution_scan");
    TTree* inputTreeEE_    = (TTree*)inFileEE_->Get("resolution_scan");
    TTree* inputTreeMusEB_ = (TTree*)inFileMusEB_->Get("resolution_scan"); 
    TTree* inputTreeMusEE_ = (TTree*)inFileMusEE_->Get("resolution_scan");  

    int maxEvents_         = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    string fitFunction_    = filesOpt.getParameter<string>( "fitFunction" );
    string etBinning_      = filesOpt.getParameter<string>( "etBinning" ); 
    string etaBinning_     = filesOpt.getParameter<string>( "etaBinning" ); 
    string dnnBinning_     = filesOpt.getParameter<string>( "dnnBinning" ); 

    std::vector<std::string> etCuts      = split(etBinning_,' ');
    std::vector<std::string> etCutsName  = etCuts;
    std::vector<std::string> etaCuts     = split(etaBinning_,' ');
    std::vector<std::string> etaCutsName = etaCuts;
    std::vector<std::string> dnnCuts     = split(dnnBinning_,' ');
    std::vector<std::string> dnnCutsName = dnnCuts;
    
    setTreeBranches(inputTreeEB_);
    setTreeBranches(inputTreeEE_);
    setTreeBranches(inputTreeMusEB_);
    setTreeBranches(inputTreeMusEE_);
    setHistograms(dnnCutsName, etCutsName, etaCutsName);

    std::cout << "Fill Mustache histos in EB..." << std::endl;
    for(int entry = 0; entry < inputTreeMusEB_->GetEntries(); entry++){

        if(entry>maxEvents_ && maxEvents_>0) continue;
        if(entry%100000==0) std::cout << "--- Reading tree = " << entry << std::endl;
        inputTreeMusEB_->GetEntry(entry);
        
        int etBin = -1;
        for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
            float minEt = std::stof(etCuts.at(iBin));
            float maxEt = std::stof(etCuts.at(iBin+1)); 
            if(minEt<90. && et_seed>minEt && et_seed<=maxEt) etBin = iBin;
            else if(minEt>=90. && et_seed>minEt) etBin = iBin;
        }

        int etaBin = -1;
        for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
            float minEta = std::stof(etaCuts.at(iBin));
            float maxEta = std::stof(etaCuts.at(iBin+1)); 
            if(seed_eta>=minEta && seed_eta<=maxEta) etaBin = iBin;
        }

        if(etBin>=0) EoEtrue_Mustache_seedEt_EB[etBin]->Fill(EoEtrue);
        if(etaBin>=0) EoEtrue_Mustache_seedEta[etaBin]->Fill(EoEtrue); 
        if(etaBin>=0 && etBin>=0) EoEtrue_Mustache_seedEta_seedEt[etBin][etaBin]->Fill(EoEtrue); 
    }
    
    std::cout << "Fill Mustache histos in EE..." << std::endl;
    for(int entry = 0; entry < inputTreeMusEE_->GetEntries(); entry++){

        if(entry>maxEvents_ && maxEvents_>0) continue;
        if(entry%100000==0) std::cout << "--- Reading tree = " << entry << std::endl;
        inputTreeMusEE_->GetEntry(entry);
        
        int etBin = -1;
        for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
            float minEt = std::stof(etCuts.at(iBin));
            float maxEt = std::stof(etCuts.at(iBin+1)); 
            if(minEt<90. && et_seed>minEt && et_seed<=maxEt) etBin = iBin;
            else if(minEt>=90. && et_seed>minEt) etBin = iBin;
        }

        int etaBin = -1;
        for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
            float minEta = std::stof(etaCuts.at(iBin));
            float maxEta = std::stof(etaCuts.at(iBin+1)); 
            if(seed_eta>=minEta && seed_eta<=maxEta) etaBin = iBin;
        }

        if(etBin>=0) EoEtrue_Mustache_seedEt_EE[etBin]->Fill(EoEtrue);
        if(etaBin>=0) EoEtrue_Mustache_seedEta[etaBin]->Fill(EoEtrue); 
        if(etaBin>=0 && etBin>=0) EoEtrue_Mustache_seedEta_seedEt[etBin][etaBin]->Fill(EoEtrue); 
    } 
    
    std::cout << "Fill DeepSC histos in EB..." << std::endl;
    for(int entry = 0; entry < inputTreeEB_->GetEntries(); entry++){

        if(entry>maxEvents_ && maxEvents_>0) continue;
        if(entry%1000000==0) std::cout << "--- Reading tree = " << entry << std::endl;
        inputTreeEB_->GetEntry(entry);

        int etBin = -1;
        for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
            float minEt = std::stof(etCuts.at(iBin));
            float maxEt = std::stof(etCuts.at(iBin+1)); 
            if(minEt<90. && et_seed>minEt && et_seed<=maxEt) etBin = iBin;
            else if(minEt>=90. && et_seed>minEt) etBin = iBin;
        }

        int etaBin = -1;
        for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
            float minEta = std::stof(etaCuts.at(iBin));
            float maxEta = std::stof(etaCuts.at(iBin+1)); 
            if(seed_eta>=minEta && seed_eta<=maxEta) etaBin = iBin;
        }

        int dnnBin = -1;
        for(unsigned int iBin=0; iBin<dnnCuts.size(); iBin++){
            float minDNN = std::stod(dnnCuts.at(iBin))-0.001;
            float maxDNN = std::stod(dnnCuts.at(iBin))+0.001;
            if(dnn_thre>=minDNN && dnn_thre<=maxDNN) dnnBin = iBin;
        }
        
        if(etBin>=0 && dnnBin>=0) EoEtrue_vs_DnnThreshold_seedEt_EB[etBin][dnnBin]->Fill(EoEtrue);
        if(etaBin>=0 && dnnBin>=0) EoEtrue_vs_DnnThreshold_seedEta[etaBin][dnnBin]->Fill(EoEtrue);
        if(etaBin>=0 && etBin>=0  && dnnBin>=0) EoEtrue_vs_DnnThreshold_seedEta_seedEt[dnnBin][etBin][etaBin]->Fill(EoEtrue);
    }
    
    std::cout << "Fill DeepSC histos in EE..." << std::endl;
    for(int entry = 0; entry < inputTreeEE_->GetEntries(); entry++){

        if(entry>maxEvents_ && maxEvents_>0) continue;
        if(entry%1000000==0) std::cout << "--- Reading tree = " << entry << std::endl;
        inputTreeEE_->GetEntry(entry);
        
        int etBin = -1;
        for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
            float minEt = std::stof(etCuts.at(iBin));
            float maxEt = std::stof(etCuts.at(iBin+1)); 
            if(minEt<90. && et_seed>minEt && et_seed<=maxEt) etBin = iBin;
            else if(minEt>=90. && et_seed>minEt) etBin = iBin;
        }

        int etaBin = -1;
        for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
            float minEta = std::stof(etaCuts.at(iBin));
            float maxEta = std::stof(etaCuts.at(iBin+1)); 
            if(seed_eta>=minEta && seed_eta<=maxEta) etaBin = iBin;
        }

        int dnnBin = -1;
        for(unsigned int iBin=0; iBin<dnnCuts.size(); iBin++){
            float minDNN = std::stod(dnnCuts.at(iBin))-0.001;
            float maxDNN = std::stod(dnnCuts.at(iBin))+0.001;
            if(dnn_thre>=minDNN && dnn_thre<=maxDNN) dnnBin = iBin;
        }
        
        if(etBin>=0 && dnnBin>=0) EoEtrue_vs_DnnThreshold_seedEt_EE[etBin][dnnBin]->Fill(EoEtrue);
        if(etaBin>=0 && dnnBin>=0) EoEtrue_vs_DnnThreshold_seedEta[etaBin][dnnBin]->Fill(EoEtrue);
        if(etaBin>=0 && etBin>=0  && dnnBin>=0) EoEtrue_vs_DnnThreshold_seedEta_seedEt[dnnBin][etBin][etaBin]->Fill(EoEtrue);
    } 

    drawPlots(fitFunction_, etCuts, etaCuts, dnnCuts);
   
}
