#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/SimScoreOptimization.h"
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
    string inputFiles_         = filesOpt.getParameter<string>( "inputFiles" );
    vector<string> inputFiles  = split(inputFiles_,',');
    
    int maxEvents_           = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    double simEnergyCut      = filesOpt.getParameter<double>( "simEnergyCut" );

    string outDir_           = filesOpt.getParameter<string>( "outDir" );
    string matching_         = filesOpt.getParameter<string>( "matching" );        
    string etBinning_        = filesOpt.getParameter<string>( "etBinning" ); 
    string etaBinning_       = filesOpt.getParameter<string>( "etaBinning" ); 
    string scoreBinning_     = filesOpt.getParameter<string>( "scoreBinning" ); 
 
    std::vector<std::string> etCuts            = split(etBinning_,' ');
    std::vector<std::string> etCutsName        = etCuts;
    std::vector<std::string> etaCuts           = split(etaBinning_,' ');
    std::vector<std::string> etaCutsName       = etaCuts;
    std::vector<std::string> scoreCuts         = split(scoreBinning_,' ');
    std::vector<std::string> scoreCutsName     = scoreCuts;
    
    std::cout << "Set histograms..." << std::endl;
    setHistograms(scoreCutsName, etCutsName, etaCutsName);

    std::vector<std::vector<int>> pfCluster_caloIndex;
    std::map<int,float> energy_simScore; 
    
    std::vector<double> scoreVec;
    for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
       scoreVec.push_back(std::stod(scoreCuts.at(iBin)));
 
    std::cout << "Fill histos inMustache..." << std::endl;
    for(unsigned int iFile=0; iFile<inputFiles.size(); ++iFile)
    {
       TFile* inFile = TFile::Open(inputFiles[iFile].c_str());
       TTree* inTree = (TTree*)inFile->Get("recosimdumper/caloTree");
       if(inFile)
       {
          cout << "\n--- Reading file: " << inputFiles[iFile].c_str() << endl;

          //setBranches
          setTreeBranches(inTree, matching_);
  
          for(int entry = 0; entry < inTree->GetEntries(); entry++){

              if(entry>maxEvents_ && maxEvents_>0) continue;
              if(entry%10000==0) std::cout << "--- Reading tree = " << entry << std::endl;
              inTree->GetEntry(entry);

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
              }else if(matching_ == "sim_fraction_withFraction"){ 
                 pfCluster_simScore = *pfCluster_sim_fraction_withFraction;
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

              pfCluster_caloIndex.clear();
              pfCluster_caloIndex.resize(pfCluster_rawEnergy->size());
              for(unsigned int iPF=0; iPF<pfCluster_rawEnergy->size(); iPF++){  
                  pfCluster_caloIndex[iPF].resize(scoreCuts.size()); 
                  for(unsigned int scoreBin=0; scoreBin<scoreCuts.size(); scoreBin++) 
                  {
                      float minScore = std::stod(scoreCuts.at(scoreBin)); 
                      int cluster_caloIndex = getMatchedIndex(&pfCluster_simScore, minScore, useMax, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), iPF);  
                      pfCluster_caloIndex[iPF][scoreBin] = cluster_caloIndex;  
                  }
              }  
        
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
               
                  if(etBin<0 || etaBin<0) continue; 

                 int superCluster_caloIndex = superCluster_sim_fraction_MatchedIndex->at(iSC);
                 int seed_caloIndex = getMatchedIndex(&pfCluster_simScore, 0., true, std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), std::vector<std::vector<std::vector<double>>>({}), std::vector<double>({}), seed_index);
             
                 if(superCluster_caloIndex<0) continue;  
                 if(seed_caloIndex<0) continue; 
                 if(superCluster_caloIndex!=seed_caloIndex) continue;   

                 for(unsigned int scoreBin = 0; scoreBin<pfCluster_caloIndex[seed_index].size(); scoreBin++)
                 {
                     MustOEtrue_vs_scoreThreshold_seedEt_seedEta[etBin][etaBin][scoreBin]->Fill(superCluster_rawEnergy->at(iSC)/caloParticle_simEnergy->at(seed_caloIndex));
                     energy_simScore[scoreBin] = pfCluster_rawEnergy->at(seed_index);   
 
                     for(unsigned int iPF=0; iPF<pfCluster_rawEnergy->size(); iPF++)  
                     {
                         if((int)iPF==seed_index) continue;

                         int cluster_caloIndex = pfCluster_caloIndex[iPF][scoreBin];  
                         if(cluster_caloIndex!=seed_caloIndex) continue;

                         float pfCluster_simEnergy = pfCluster_simScore.at(iPF).at(seed_caloIndex)*caloParticle_simEnergy->at(seed_caloIndex); 
                         if(pfCluster_simEnergy<simEnergyCut) continue;  
                         energy_simScore[scoreBin] += pfCluster_rawEnergy->at(iPF); 
                     }

                     float EoEtrue = energy_simScore[scoreBin]/caloParticle_simEnergy->at(seed_caloIndex);
                     if(EoEtrue!=0.) CaloOEtrue_vs_scoreThreshold_seedEt_seedEta[etBin][etaBin][scoreBin]->Fill(EoEtrue);
                 }
              }       
          }  
       }
    }
     
    scoreOptimization(etCuts, etaCuts, scoreCuts, matching_, outDir_);
    
    TFile* optimization_scores = new TFile("simScore_Minima_ElectronsOnly.root","RECREATE"); 
    optimization_scores->cd();
    h2_Minimum_simScore->Write();
    h2_Minimum_Ratio_simScore->Write();
    h2_GainToMustache_simScore->Write();
    h2_GainToMustache_Ratio_simScore->Write();
    optimization_scores->Close();
}
