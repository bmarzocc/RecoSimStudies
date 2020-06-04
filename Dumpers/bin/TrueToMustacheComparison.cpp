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
    
    string outDir_           = filesOpt.getParameter<string>( "outDir" );
    int maxEvents_           = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    string fitFunction_      = filesOpt.getParameter<string>( "fitFunction" );        
    string etBinning_        = filesOpt.getParameter<string>( "etBinning" ); 
    string etaBinning_       = filesOpt.getParameter<string>( "etaBinning" );

    std::vector<std::string> etCuts        = split(etBinning_,' ');
    std::vector<std::string> etCutsName    = etCuts;
    std::vector<std::string> etaCuts       = split(etaBinning_,' ');
    std::vector<std::string> etaCutsName   = etaCuts;
    
    setTreeBranches(inputTree_);
    setHistograms(etCutsName, etaCutsName);

    std::vector<int> pfCluster_caloIndex;

    std::cout << "Fill histos..." << std::endl;
    for(int entry = 0; entry < inputTree_->GetEntries(); entry++){

        if(entry>maxEvents_ && maxEvents_>0) continue;
        if(entry%10000==0) std::cout << "--- Reading tree = " << entry << std::endl;
        inputTree_->GetEntry(entry);

        for(unsigned int iSC=0; iSC<superCluster_seedIndex->size(); iSC++)
        {

            int seed_index = superCluster_seedIndex->at(iSC);
            float seed_eta = pfCluster_eta->at(seed_index);  
            float seed_et = pfCluster_energy->at(seed_index)/TMath::CosH(pfCluster_eta->at(seed_index)); 
            int etBin = -1;
            for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
                float minEt = std::stof(etCuts.at(iBin));
                float maxEt = std::stof(etCuts.at(iBin+1)); 
                if(minEt<90. && seed_et>minEt && seed_et<=maxEt) etBin = iBin;
                else if(minEt>=90. && seed_et>minEt) etBin = iBin;
            }
            int etaBin = -1;
            for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
                float minEta = std::stof(etaCuts.at(iBin));
                float maxEta = std::stof(etaCuts.at(iBin+1)); 
                if(fabs(seed_eta)>minEta && fabs(seed_eta)<=maxEta) etaBin = iBin;
            }
       
            pfCluster_caloIndex.clear();
            pfCluster_caloIndex.resize(pfCluster_rawEnergy->size());
            for(unsigned int iPF=0; iPF<pfCluster_rawEnergy->size(); iPF++){  
                int cluster_caloIndex = getMatchedIndex(pfCluster_sim_fraction, 0.002, true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPF);  
                pfCluster_caloIndex[iPF]= cluster_caloIndex;  
            }  
 
            int superCluster_caloIndex = superCluster_sim_fraction_MatchedIndex->at(iSC);
            int seed_caloIndex = getMatchedIndex(pfCluster_sim_fraction, 0., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), seed_index);
           
            if(superCluster_caloIndex<0) continue;  
            if(seed_caloIndex<0) continue; 
            if(superCluster_caloIndex!=seed_caloIndex) continue;   

            for(unsigned int iPF=0; iPF<pfCluster_rawEnergy->size(); iPF++){

                int cluster_caloIndex = pfCluster_caloIndex[iPF];
                if(cluster_caloIndex!=seed_caloIndex) continue;

                float dphi = TVector2::Phi_mpi_pi(pfCluster_phi->at(seed_index)-pfCluster_phi->at(iPF)); 
                float deta = pfCluster_eta->at(seed_index)-pfCluster_eta->at(iPF); 
                if(pfCluster_eta->at(seed_index)>0.) deta = -deta; 
                 
                bool pass_dphi_must = reco::MustacheKernel::inDynamicDPhiWindow(pfCluster_eta->at(seed_index), pfCluster_phi->at(seed_index), pfCluster_rawEnergy->at(iPF), pfCluster_eta->at(iPF), pfCluster_phi->at(iPF)); 
                bool pass_deta_must = reco::MustacheKernel::inMustache(pfCluster_eta->at(seed_index), pfCluster_phi->at(seed_index), pfCluster_rawEnergy->at(iPF), pfCluster_eta->at(iPF), pfCluster_phi->at(iPF)); 

                double simScore = pfCluster_sim_fraction->at(iPF).at(seed_caloIndex); 
                double simScore_withFraction = pfCluster_sim_fraction_withFraction->at(iPF).at(seed_caloIndex); 
                double RecoToCalo = pfCluster_rawEnergy->at(iPF)/(simScore*caloParticle_simEnergy->at(seed_caloIndex));
                double recoToSeed = pfCluster_rawEnergy->at(iPF)/pfCluster_rawEnergy->at(seed_index); 
                if(simScore<=0. || simScore_withFraction<=0.) continue;
   
                if(etaBin>=0 && etBin>=0)
                { 
                   if(dphi!=0. && deta!=0.) dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta);
                   SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,simScore); 
                   SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,simScore_withFraction); 
                   RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,RecoToCalo); 
                   RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,recoToSeed);
               
                   if(pass_dphi_must && pass_deta_must){ 
                      if(dphi!=0. && deta!=0.) dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin]->Fill(dphi,deta);
                      SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin]->Fill(dphi,deta,simScore); 
                      SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin]->Fill(dphi,deta,simScore_withFraction);  
                      RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin]->Fill(dphi,deta,RecoToCalo); 
                      RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin]->Fill(dphi,deta,recoToSeed);  
                   } 
                   if(!pass_dphi_must || !pass_deta_must){ 
                      if(dphi!=0. && deta!=0.) dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin]->Fill(dphi,deta);
                      SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin]->Fill(dphi,deta,simScore);
                      SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin]->Fill(dphi,deta,simScore_withFraction);
                      RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin]->Fill(dphi,deta,RecoToCalo); 
                      RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin]->Fill(dphi,deta,recoToSeed);  
                   } 
                   
                   for(unsigned int hit=0; hit<pfClusterHit_fraction->at(iPF).size(); hit++)
                   {
                       if(pfClusterHit_phi->at(iPF).at(hit)<=0.) continue;
                       //dphi = TVector2::Phi_mpi_pi(pfCluster_phi->at(seed_index)-pfClusterHit_phi->at(iPF).at(hit)); 
                       //deta = pfCluster_eta->at(seed_index)-pfClusterHit_eta->at(iPF).at(hit); 
                       //if(pfCluster_eta->at(seed_index)>0.) deta = -deta; 
             
                       HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,pfClusterHit_fraction->at(iPF).at(hit));
                       if(pass_dphi_must && pass_deta_must) HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin]->Fill(dphi,deta,pfClusterHit_fraction->at(iPF).at(hit));
                       if(!pass_dphi_must || !pass_deta_must) HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin]->Fill(dphi,deta,pfClusterHit_fraction->at(iPF).at(hit));
                   } 
                } 
            }    
  
         
            for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex->at(iSC).size(); iPF++)
            {  
                int pfCluster_index = superCluster_pfClustersIndex->at(iSC).at(iPF);          
                float dphi = TVector2::Phi_mpi_pi(pfCluster_phi->at(seed_index)-pfCluster_phi->at(pfCluster_index)); 
                float deta = pfCluster_eta->at(seed_index)-pfCluster_eta->at(pfCluster_index); 
                if(pfCluster_eta->at(seed_index)>0.) deta = -deta;
       
                int cluster_caloIndex = pfCluster_caloIndex[pfCluster_index];        
                if(etaBin>=0 && etBin>=0 && dphi!=0. && deta!=0.){  
                   dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta); 
                   if(cluster_caloIndex==seed_caloIndex) dEtadPhi_Mustache_vs_seedEt_seedEta_inCalo[etBin][etaBin]->Fill(dphi,deta); 
                   if(cluster_caloIndex!=seed_caloIndex) dEtadPhi_Mustache_vs_seedEt_seedEta_notCalo[etBin][etaBin]->Fill(dphi,deta); 
                }
            }
        } 
    }
    
    drawPlots(etCuts, etaCuts, outDir_);
   
}
