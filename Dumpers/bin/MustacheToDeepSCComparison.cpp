#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/MustacheToDeepSCComparison.h"
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
    string inputFiles_        = filesOpt.getParameter<string>( "inputFiles" );
    vector<string> inputFiles = split(inputFiles_,',');
    
    string outDir_           = filesOpt.getParameter<string>( "outDir" );
    int maxEvents_           = filesOpt.getUntrackedParameter<int>( "maxEvents" );       
    string etBinning_        = filesOpt.getParameter<string>( "etBinning" ); 
    string etaBinning_       = filesOpt.getParameter<string>( "etaBinning" );

    std::vector<std::string> etCuts        = split(etBinning_,' ');
    std::vector<std::string> etCutsName    = etCuts;
    std::vector<std::string> etaCuts       = split(etaBinning_,' ');
    std::vector<std::string> etaCutsName   = etaCuts;
    
    setHistograms(etCutsName, etaCutsName);

    std::vector<int> pfCluster_caloIndex;

    std::cout << "Fill histos..." << std::endl;
    for(unsigned int iFile=0; iFile<inputFiles.size(); ++iFile)
    {
       TFile* inFile = TFile::Open(inputFiles[iFile].c_str());
       TTree* inTree = (TTree*)inFile->Get("recosimdumper/caloTree");
       if(inFile)
       {
          cout << "\n--- Reading file: " << inputFiles[iFile].c_str() << endl;

          //setBranches
          setTreeBranches(inTree);

          // Loop over all entries of the TTree
          for(int entry = 0; entry < inTree->GetEntries(); entry++)
          {
          
              if(entry>maxEvents_ && maxEvents_>0) continue;
              if(entry%10000==0) std::cout << "--- Reading tree = " << entry << std::endl;
              inTree->GetEntry(entry);

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
       
                  for(unsigned int iPF=0; iPF<superCluster_pfClustersIndex->at(iSC).size(); iPF++)
                  {  
                      int pfCluster_index = superCluster_pfClustersIndex->at(iSC).at(iPF);          
                      float dphi = deltaPhi(pfCluster_phi->at(seed_index),pfCluster_phi->at(pfCluster_index)); 
                      float deta = deltaEta(pfCluster_eta->at(seed_index),pfCluster_eta->at(pfCluster_index)); 
               
                      if(etaBin>=0 && etBin>=0){  
                         Energy_dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,pfCluster_rawEnergy->at(pfCluster_index));
                         ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,pfCluster_rawEnergy->at(pfCluster_index)/pfCluster_rawEnergy->at(seed_index));
                         if(dphi!=0. && deta!=0.) dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta); 
                      } 
                  }
              }

              for(unsigned int iSC=0; iSC<deepSuperCluster_seedIndex->size(); iSC++)
              {

                  int seed_index = deepSuperCluster_seedIndex->at(iSC);
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
       
                  for(unsigned int iPF=0; iPF<deepSuperCluster_pfClustersIndex->at(iSC).size(); iPF++)
                  {  
                      int pfCluster_index = deepSuperCluster_pfClustersIndex->at(iSC).at(iPF);          
                      float dphi = deltaPhi(pfCluster_phi->at(seed_index),pfCluster_phi->at(pfCluster_index)); 
                      float deta = deltaEta(pfCluster_eta->at(seed_index),pfCluster_eta->at(pfCluster_index)); 
               
                      if(etaBin>=0 && etBin>=0){  
                         Energy_dEtadPhi_DeepSC_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,pfCluster_rawEnergy->at(pfCluster_index));
                         ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta,pfCluster_rawEnergy->at(pfCluster_index)/pfCluster_rawEnergy->at(seed_index));
                         if(dphi!=0. && deta!=0.) dEtadPhi_DeepSC_vs_seedEt_seedEta[etBin][etaBin]->Fill(dphi,deta); 
                      }  
                  }
              } 
           }   
        }   
    }
    
    drawPlots(etCuts, etaCuts, outDir_);
   
}
