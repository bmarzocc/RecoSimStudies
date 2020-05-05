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

    string fitFunction_    = filesOpt.getParameter<string>( "fitFunction" );
    string dnnBinning_     = filesOpt.getParameter<string>( "dnnBinning" ); 

    std::vector<std::string> dnnCuts = split(dnnBinning_,' ');
    std::vector<std::string> dnnCutsName = dnnCuts;

    setTreeBranches(inputTreeEB_);
    setTreeBranches(inputTreeEE_);
    setTreeBranches(inputTreeMusEB_);
    setTreeBranches(inputTreeMusEE_);
    setHistograms(dnnCutsName);

    std::cout << "Fill Mustache histos..." << std::endl;
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_0_10_EB","et_seed>0 && et_seed<=10");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_10_20_EB","et_seed>10 && et_seed<=20");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_20_30_EB","et_seed>20 && et_seed<=30");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_30_40_EB","et_seed>30 && et_seed<=40");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_40_50_EB","et_seed>40 && et_seed<=50");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_50_60_EB","et_seed>50 && et_seed<=60");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_60_70_EB","et_seed>60 && et_seed<=70");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_70_80_EB","et_seed>70 && et_seed<=80");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_80_90_EB","et_seed>80 && et_seed<=90");
    inputTreeMusEB_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_90_100_EB","et_seed>90");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_0_10_EE","et_seed>0 && et_seed<=10");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_10_20_EE","et_seed>10 && et_seed<=20");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_20_30_EE","et_seed>20 && et_seed<=30");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_30_40_EE","et_seed>30 && et_seed<=40");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_40_50_EE","et_seed>40 && et_seed<=50");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_50_60_EE","et_seed>50 && et_seed<=60");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_60_70_EE","et_seed>60 && et_seed<=70");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_70_80_EE","et_seed>70 && et_seed<=80");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_80_90_EE","et_seed>80 && et_seed<=90");
    inputTreeMusEE_->Draw("EoEtrue >> EoEtrue_Mustache_seedEt_90_100_EE","et_seed>90");
    
    std::cout << "Fill DeepSC histos..." << std::endl;
    for(unsigned int iBin=0; iBin<dnnCuts.size(); iBin++)
    {
        dnnCutsName.at(iBin).replace(1,1,string("_")); 
        string cut_min = std::to_string(std::stod(dnnCuts.at(iBin))-0.001);
        string cut_max = std::to_string(std::stod(dnnCuts.at(iBin))+0.001);
        
        std::cout << " DNN score = " << dnnCuts.at(iBin) << std::endl;
   
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_0_10_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>0 && et_seed<=10 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_10_20_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>10 && et_seed<=20 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_20_30_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>20 && et_seed<=30 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_30_40_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>30 && et_seed<=40 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_40_50_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>40 && et_seed<=50 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_50_60_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>50 && et_seed<=60 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_60_70_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>60 && et_seed<=70 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_70_80_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>70 && et_seed<=80 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());  
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_80_90_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>80 && et_seed<=90 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEB_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_90_100_DNN_"+dnnCutsName.at(iBin)+"_EB").c_str(),string("et_seed>90 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str()); 
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_0_10_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>0 && et_seed<=10 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_10_20_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>10 && et_seed<=20 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_20_30_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>20 && et_seed<=30 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_30_40_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>30 && et_seed<=40 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_40_50_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>40 && et_seed<=50 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_50_60_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>50 && et_seed<=60 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_60_70_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>60 && et_seed<=70 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_70_80_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>70 && et_seed<=80 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());  
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_80_90_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>80 && et_seed<=90 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());
        inputTreeEE_->Draw(string("EoEtrue >> EoEtrue_vs_DnnThreshold_seedEt_90_100_DNN_"+dnnCutsName.at(iBin)+"_EE").c_str(),string("et_seed>90 && dnn_thre>="+cut_min+" && dnn_thre<="+cut_max).c_str());   
    }  

    drawPlots(fitFunction_, dnnCuts);
   
}
