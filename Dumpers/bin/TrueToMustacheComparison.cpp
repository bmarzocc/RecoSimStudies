#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/TrueToMustacheComparison.h"
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
    string inputFile_    = filesOpt.getParameter<string>( "inputFile" );
    TFile* inFile_       = TFile::Open(inputFile_.c_str());
    TTree* inputTree_    = (TTree*)inFile_->Get("recosimdumper/caloTree");
    
    int maxEvents_         = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    string fitFunction_    = filesOpt.getParameter<string>( "fitFunction" );
    string matching_       = filesOpt.getParameter<string>( "matching" ); 
    string etBinning_      = filesOpt.getParameter<string>( "etBinning" ); 
    string etaBinning_     = filesOpt.getParameter<string>( "etaBinning" ); 
    string scoreBinning_   = filesOpt.getParameter<string>( "scoreBinning" ); 

    std::vector<std::string> etCuts        = split(etBinning_,' ');
    std::vector<std::string> etCutsName    = etCuts;
    std::vector<std::string> etaCuts       = split(etaBinning_,' ');
    std::vector<std::string> etaCutsName   = etaCuts;
    std::vector<std::string> scoreCuts     = split(scoreBinning_,' ');
    std::vector<std::string> scoreCutsName = scoreCuts;
    
    setTreeBranches(inputTree_, matching_);
    setHistograms(scoreCutsName, etCutsName, etaCutsName);

    std::cout << "Fill histos..." << std::endl;
    for(int entry = 0; entry < inputTree_->GetEntries(); entry++){

        if(entry>maxEvents_ && maxEvents_>0) continue;
        if(entry%10000==0) std::cout << "--- Reading tree = " << entry << std::endl;
        inputTree_->GetEntry(entry);

        for(unsigned int iSeed=0; iSeed<superCluster_seedIndex->size(); iSeed++)
        {
            float seed_eta = pfCluster_eta->at(iSeed);  
            float et_seed = pfCluster_energy->at(iSeed)/TMath::CosH(pfCluster_eta->at(iSeed)); 
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
       
            vector<vector<double>> pfCluster_simScore;
            bool useMax = true;

            if(matching_ == "dR_genScore"){ 
               pfCluster_simScore = *pfCluster_dR_genScore;
               useMax = false;  
            }else if(matching_ == "dR_simScore"){
               pfCluster_simScore = *pfCluster_dR_simScore;
               useMax = false;
            }else if(matching_ == "sim_fraction_old"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_old;
               useMax = false;
            }else if(matching_ == "pfCluster_n_shared_xtals"){
               pfCluster_simScore = *pfCluster_n_shared_xtals;
               useMax = true;
            }else if(matching_ == "sim_fraction"){ 
               pfCluster_simScore = *pfCluster_sim_fraction;
               useMax = true;
            }else if(matching_ == "sim_fraction_1MeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_1MeVCut;
               useMax = true;
            }else if(matching_ == "sim_fraction_5MeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_5MeVCut;
               useMax = true;
            }else if(matching_ == "sim_fraction_10MeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_10MeVCut;
               useMax = true;
            }else if(matching_ == "sim_fraction_50MeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_50MeVCut; 
               useMax = true;
            }else if(matching_ == "sim_fraction_100MeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_100MeVCut;
               useMax = true;
            }else if(matching_ == "sim_fraction_500MeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_500MeVCut;
               useMax = true;
            }else if(matching_ == "sim_fraction_1GeVCut"){ 
               pfCluster_simScore = *pfCluster_sim_fraction_1GeVCut;
               useMax = true;
            }else if(matching_ == "sim_rechit_diff"){
               pfCluster_simScore = *pfCluster_sim_rechit_diff;
               useMax = false;
            }else if(matching_ == "sim_rechit_fraction"){
               pfCluster_simScore = *pfCluster_sim_rechit_fraction;
               useMax = false;
            }else if(matching_ == "global_sim_rechit_fraction"){ 
               pfCluster_simScore = *pfCluster_global_sim_rechit_fraction; 
               useMax = false; 
            }
 
            for(unsigned int scoreBin=0; scoreBin<scoreCuts.size(); scoreBin++){  

                float minScore = std::stod(scoreCuts.at(scoreBin)); 
      
                int seed_matchedIndex = getMatchedIndex(&pfCluster_simScore, minScore, useMax, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), superCluster_seedIndex->at(iSeed)); 
                if(seed_matchedIndex<0) continue;

                //for(unsigned int iCalo=0; iCalo<pfCluster_simScore.at(superCluster_seedIndex->at(iSeed)).size() ; iCalo++)
                //    std::cout << superCluster_seedIndex->at(iSeed) << " - " << iCalo << " - " << pfCluster_simScore.at(superCluster_seedIndex->at(iSeed)).at(iCalo) << " - " << seed_matchedIndex << std::endl;
 
                float energy_tot = 0.; 
                for(unsigned int iPF=0; iPF<pfCluster_energy->size(); iPF++){
                    int pfCluster_matchedIndex = getMatchedIndex(&pfCluster_simScore, minScore, useMax, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPF);
 
                    if(pfCluster_matchedIndex==seed_matchedIndex){ 
                       energy_tot += pfCluster_rawEnergy->at(iPF); 
                       float dphi = TVector2::Phi_mpi_pi(pfCluster_phi->at(superCluster_seedIndex->at(iSeed))-pfCluster_phi->at(iPF)); 
                       float deta = pfCluster_eta->at(superCluster_seedIndex->at(iSeed))-pfCluster_eta->at(iPF);
                       if(etaBin>=0 && etBin>=0 && dphi!=0. && deta!=0.) dEtadPhi_vs_scoreThreshold_seedEta_seedEt[scoreBin][etBin][etaBin]->Fill(dphi,deta);
                    } 
                }

                float EoEtrue_must = superCluster_rawEnergy->at(iSeed)/caloParticle_simEnergy->at(seed_matchedIndex); 
                if(etBin>=0 && superCluster_iz->at(iSeed)==0) EoEtrue_vs_Mustache_seedEt_EB[etBin][scoreBin]->Fill(EoEtrue_must);
                if(etBin>=0 && superCluster_iz->at(iSeed)!=0) EoEtrue_vs_Mustache_seedEt_EE[etBin][scoreBin]->Fill(EoEtrue_must);
                if(etaBin>=0) EoEtrue_vs_Mustache_seedEta[etaBin][scoreBin]->Fill(EoEtrue_must);
                if(etaBin>=0 && etBin>=0) EoEtrue_vs_Mustache_seedEta_seedEt[scoreBin][etBin][etaBin]->Fill(EoEtrue_must);

                float EoEtrue = energy_tot/caloParticle_simEnergy->at(seed_matchedIndex); 
                if(etBin>=0 && superCluster_iz->at(iSeed)==0) EoEtrue_vs_scoreThreshold_seedEt_EB[etBin][scoreBin]->Fill(EoEtrue);
                if(etBin>=0 && superCluster_iz->at(iSeed)!=0) EoEtrue_vs_scoreThreshold_seedEt_EE[etBin][scoreBin]->Fill(EoEtrue);
                if(etaBin>=0) EoEtrue_vs_scoreThreshold_seedEta[etaBin][scoreBin]->Fill(EoEtrue);
                if(etaBin>=0 && etBin>=0) EoEtrue_vs_scoreThreshold_seedEta_seedEt[scoreBin][etBin][etaBin]->Fill(EoEtrue);


                //std::cout << etaBin << " - " << etBin << " - " << scoreBin << std::endl; 
            }
        } 
    }
    
    drawPlots(fitFunction_, etCuts, etaCuts, scoreCuts, matching_);
   
}
