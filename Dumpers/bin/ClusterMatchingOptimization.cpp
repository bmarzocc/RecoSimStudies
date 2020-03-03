#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/ClusterMatchingOptimization.h"
#include "RecoSimStudies/Dumpers/interface/CruijffPdf.h"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
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
    string inputFiles_ = filesOpt.getParameter<string>( "inputFiles" );
    vector<string> inputFiles = split(inputFiles_,',');
    string outputDir_ = filesOpt.getParameter<string>( "outputDir" );
    if( outputDir_ == "" ) outputDir_ = "output/"; 
    int maxEvents_ = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    fitFunction_ = filesOpt.getParameter<string>( "fitFunction" );
    
    bool gotoMain = false;
    for(unsigned int iFile=0; iFile<inputFiles.size() && !gotoMain; ++iFile)
    {
       TFile* inFile = TFile::Open(inputFiles[iFile].c_str());
       TTree* inTree = (TTree*)inFile->Get("recosimdumper/caloTree");
       if(inFile)
       {
          cout << "\n--- Reading file: " << inputFiles[iFile].c_str() << endl;

          //setBranches
          setTreeBranches(inTree);

          // Loop over all entries of the TTree
          for(int entry = 0; entry < inTree->GetEntries(); entry++){
     
             if(entry>maxEvents_ && maxEvents_>0) continue;
             if(entry%10000==0) std::cout << "--- Reading tree = " << entry << std::endl;
             inTree->GetEntry(entry);
              
             fillHisto(pfCluster_sim_fraction_old, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_01, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_01, false);  
             fillHisto(pfCluster_sim_fraction_old, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_02, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_02, false);  
             fillHisto(pfCluster_sim_fraction_old, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_03, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_03, false);  
             fillHisto(pfCluster_sim_fraction_old, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_04, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_04, false);  
     	     fillHisto(pfCluster_sim_fraction_old, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_05, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_05, false);     
             fillHisto(pfCluster_sim_fraction_old, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_06, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_06, false);     
             fillHisto(pfCluster_sim_fraction_old, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_07, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_07, false);     
             fillHisto(pfCluster_sim_fraction_old, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_08, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_08, false);     
             fillHisto(pfCluster_sim_fraction_old, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_09, true);      
             fillHisto(pfCluster_sim_fraction_old, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_09, false);  
             fillHisto(pfCluster_sim_fraction_old, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_old_10, true);      
             fillHisto(pfCluster_sim_fraction_old, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_old_10, false);            
                            
             fillHisto(pfCluster_sim_fraction, 0.01, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_001, true); 
             fillHisto(pfCluster_sim_fraction, 0.01, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_001, false); 
             fillHisto(pfCluster_sim_fraction, 0.02, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_002, true); 
             fillHisto(pfCluster_sim_fraction, 0.02, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_002, true); 
       	     fillHisto(pfCluster_sim_fraction, 0.03, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_003, true); 
             fillHisto(pfCluster_sim_fraction, 0.03, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_003, true);   
             fillHisto(pfCluster_sim_fraction, 0.04, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_004, true); 
             fillHisto(pfCluster_sim_fraction, 0.04, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_004, true); 
             fillHisto(pfCluster_sim_fraction, 0.05, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_005, true); 
             fillHisto(pfCluster_sim_fraction, 0.05, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_005, true); 
       	     fillHisto(pfCluster_sim_fraction, 0.06, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_006, true); 
             fillHisto(pfCluster_sim_fraction, 0.06, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_006, true);    
             fillHisto(pfCluster_sim_fraction, 0.07, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_007, true); 
             fillHisto(pfCluster_sim_fraction, 0.07, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_007, true);   
             fillHisto(pfCluster_sim_fraction, 0.08, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_008, true); 
             fillHisto(pfCluster_sim_fraction, 0.08, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_008, true); 
             fillHisto(pfCluster_sim_fraction, 0.09, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_009, true); 
             fillHisto(pfCluster_sim_fraction, 0.09, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_009, true); 
       	     fillHisto(pfCluster_sim_fraction, 0.1, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_01, true); 
             fillHisto(pfCluster_sim_fraction, 0.1, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_01, true);         	     
             fillHisto(pfCluster_sim_fraction, 0.2, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_02, true); 
             fillHisto(pfCluster_sim_fraction, 0.2, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_02, true);  
             fillHisto(pfCluster_sim_fraction, 0.3, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_03, true); 
             fillHisto(pfCluster_sim_fraction, 0.3, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_03, true);  
             fillHisto(pfCluster_sim_fraction, 0.4, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_04, true); 
             fillHisto(pfCluster_sim_fraction, 0.4, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_04, true);     
             fillHisto(pfCluster_sim_fraction, 0.5, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_05, true); 
             fillHisto(pfCluster_sim_fraction, 0.5, true, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_05, true);   
             
             fillHisto(pfCluster_sim_rechit_diff, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_01, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_01, false);  
             fillHisto(pfCluster_sim_rechit_diff, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_02, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_02, false);  
             fillHisto(pfCluster_sim_rechit_diff, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_03, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_03, false);  
             fillHisto(pfCluster_sim_rechit_diff, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_04, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_04, false);  
     	     fillHisto(pfCluster_sim_rechit_diff, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_05, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_05, false);     
             fillHisto(pfCluster_sim_rechit_diff, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_06, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_06, false);     
             fillHisto(pfCluster_sim_rechit_diff, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_07, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_07, false);     
             fillHisto(pfCluster_sim_rechit_diff, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_08, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_08, false);     
             fillHisto(pfCluster_sim_rechit_diff, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_09, true);      
             fillHisto(pfCluster_sim_rechit_diff, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_09, false);  
             fillHisto(pfCluster_sim_rechit_diff, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_diff_10, true);      
             fillHisto(pfCluster_sim_rechit_diff, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_diff_10, false);  

             fillHisto(pfCluster_sim_rechit_fraction, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_01, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_01, false);  
             fillHisto(pfCluster_sim_rechit_fraction, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_02, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_02, false);  
             fillHisto(pfCluster_sim_rechit_fraction, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_03, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_03, false);  
             fillHisto(pfCluster_sim_rechit_fraction, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_04, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_04, false);  
     	     fillHisto(pfCluster_sim_rechit_fraction, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_05, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_05, false);     
             fillHisto(pfCluster_sim_rechit_fraction, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_06, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_06, false);     
             fillHisto(pfCluster_sim_rechit_fraction, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_07, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_07, false);     
             fillHisto(pfCluster_sim_rechit_fraction, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_08, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_08, false);     
             fillHisto(pfCluster_sim_rechit_fraction, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_09, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_09, false);  
             fillHisto(pfCluster_sim_rechit_fraction, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_fraction_10, true);      
             fillHisto(pfCluster_sim_rechit_fraction, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_fraction_10, false);  

             fillHisto(pfCluster_global_sim_rechit_fraction, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_01, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.1, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_01, false);  
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_02, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.2, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_02, false);  
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_03, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.3, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_03, false);  
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_04, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.4, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_04, false);  
     	     fillHisto(pfCluster_global_sim_rechit_fraction, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_05, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.5, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_05, false);     
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_06, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.6, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_06, false);     
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_07, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.7, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_07, false);     
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_08, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.8, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_08, false);     
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_09, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 0.9, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_09, false);  
             fillHisto(pfCluster_global_sim_rechit_fraction, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_global_sim_rechit_fraction_10, true);      
             fillHisto(pfCluster_global_sim_rechit_fraction, 1.0, false, vector<vector<vector<double>>>({}), vector<double>({}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_global_sim_rechit_fraction_10, false);

             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.1}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_combined_fraction01, true);
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.1}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_combined_fraction01, false);
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.2}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_combined_fraction02, true); 
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.2}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_combined_fraction02, false); 
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.3}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_combined_fraction03, true);
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.3}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_combined_fraction03, false); 
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.4}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_combined_fraction04, false);   
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.4}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_combined_fraction04, true);     
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.5}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_rechit_combined_fraction05, true);  
             fillHisto(pfCluster_sim_fraction, 0., true, vector<vector<vector<double>>>({*pfCluster_sim_rechit_combined_fraction}), vector<double>({0.5}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_rechit_combined_fraction05, false);                 
   
             fillHisto(pfCluster_sim_fraction_1MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_1MeVCut, true); 
             fillHisto(pfCluster_sim_fraction_1MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_1MeVCut, false);
             fillHisto(pfCluster_sim_fraction_5MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_5MeVCut, true); 
             fillHisto(pfCluster_sim_fraction_5MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_5MeVCut, false);
             fillHisto(pfCluster_sim_fraction_10MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_10MeVCut, true); 
             fillHisto(pfCluster_sim_fraction_10MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_10MeVCut, false);    
             fillHisto(pfCluster_sim_fraction_50MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_50MeVCut, true); 
             fillHisto(pfCluster_sim_fraction_50MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_50MeVCut, false); 
             fillHisto(pfCluster_sim_fraction_100MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_100MeVCut, true); 
             fillHisto(pfCluster_sim_fraction_100MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_100MeVCut, false); 
             fillHisto(pfCluster_sim_fraction_500MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_500MeVCut, true); 
             fillHisto(pfCluster_sim_fraction_500MeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_500MeVCut, false); 
             fillHisto(pfCluster_sim_fraction_1GeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_sim_fraction_1GeVCut, true); 
             fillHisto(pfCluster_sim_fraction_1GeVCut, 0.01, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction}), vector<double>({0.01}), vector<vector<vector<double>>>({}), vector<double>({}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_sim_fraction_1GeVCut, false); 

             fillHisto(pfCluster_sim_fraction, 0.04, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction_100MeVCut}), vector<double>({0.01}), vector<vector<vector<double>>>({*pfCluster_sim_fraction_old,*pfCluster_global_sim_rechit_fraction}), vector<double>({0.8,0.5}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_simScore_final_combination, true); 
             //fillHisto(pfCluster_sim_fraction, 0.04, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction_100MeVCut}), vector<double>({0.01}), vector<vector<vector<double>>>({*pfCluster_sim_fraction_old,*pfCluster_global_sim_rechit_fraction}), vector<double>({0.8,1.0}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EB_simScore_final_combination, true); 
             //fillHisto(pfCluster_sim_fraction, 0.04, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction_100MeVCut}), vector<double>({0.01}), vector<vector<vector<double>>>({*pfCluster_sim_fraction_old,*pfCluster_global_sim_rechit_fraction}), vector<double>({0.8,0.5}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_simScore_final_combination, false); 
             fillHisto(pfCluster_sim_fraction, 0.04, true, vector<vector<vector<double>>>({*pfCluster_sim_fraction_100MeVCut}), vector<double>({0.01}), vector<vector<vector<double>>>({*pfCluster_sim_fraction_old,*pfCluster_global_sim_rechit_fraction}), vector<double>({0.1,0.5}), pfCluster_energy, caloParticle_simEnergy, caloParticle_simEta, h_EoEtrue_EE_simScore_final_combination, false); 
          }
       }
    }

    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_01,fitHisto(h_EoEtrue_EB_sim_fraction_old_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_01"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_01,fitHisto(h_EoEtrue_EE_sim_fraction_old_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_01"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_02,fitHisto(h_EoEtrue_EB_sim_fraction_old_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_02"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_02,fitHisto(h_EoEtrue_EE_sim_fraction_old_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_02"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_03,fitHisto(h_EoEtrue_EB_sim_fraction_old_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_03"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_03,fitHisto(h_EoEtrue_EE_sim_fraction_old_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_03")); 
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_04,fitHisto(h_EoEtrue_EB_sim_fraction_old_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_04"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_04,fitHisto(h_EoEtrue_EE_sim_fraction_old_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_04")); 
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_05,fitHisto(h_EoEtrue_EB_sim_fraction_old_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_05"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_05,fitHisto(h_EoEtrue_EE_sim_fraction_old_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_05")); 
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_06,fitHisto(h_EoEtrue_EB_sim_fraction_old_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_06"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_06,fitHisto(h_EoEtrue_EE_sim_fraction_old_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_06"));  
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_07,fitHisto(h_EoEtrue_EB_sim_fraction_old_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_07"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_07,fitHisto(h_EoEtrue_EE_sim_fraction_old_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_07"));  
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_08,fitHisto(h_EoEtrue_EB_sim_fraction_old_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_08"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_08,fitHisto(h_EoEtrue_EE_sim_fraction_old_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_08"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_09,fitHisto(h_EoEtrue_EB_sim_fraction_old_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_09"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_09,fitHisto(h_EoEtrue_EE_sim_fraction_old_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_09"));    
    drawHistFunc(h_EoEtrue_EB_sim_fraction_old_10,fitHisto(h_EoEtrue_EB_sim_fraction_old_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_old_10"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_old_10,fitHisto(h_EoEtrue_EE_sim_fraction_old_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_old_10")); 
    drawHistFunc(h_EoEtrue_EB_sim_fraction_001,fitHisto(h_EoEtrue_EB_sim_fraction_001), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_001"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_001,fitHisto(h_EoEtrue_EE_sim_fraction_001), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_001"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_002,fitHisto(h_EoEtrue_EB_sim_fraction_002), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_002"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_002,fitHisto(h_EoEtrue_EE_sim_fraction_002), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_002"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_003,fitHisto(h_EoEtrue_EB_sim_fraction_003), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_003"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_003,fitHisto(h_EoEtrue_EE_sim_fraction_003), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_003"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_004,fitHisto(h_EoEtrue_EB_sim_fraction_004), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_004"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_004,fitHisto(h_EoEtrue_EE_sim_fraction_004), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_004"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_005,fitHisto(h_EoEtrue_EB_sim_fraction_005), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_005"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_005,fitHisto(h_EoEtrue_EE_sim_fraction_005), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_005"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_006,fitHisto(h_EoEtrue_EB_sim_fraction_006), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_006"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_006,fitHisto(h_EoEtrue_EE_sim_fraction_006), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_006"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_007,fitHisto(h_EoEtrue_EB_sim_fraction_007), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_007"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_007,fitHisto(h_EoEtrue_EE_sim_fraction_007), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_007"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_008,fitHisto(h_EoEtrue_EB_sim_fraction_008), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_008"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_008,fitHisto(h_EoEtrue_EE_sim_fraction_008), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_008"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_009,fitHisto(h_EoEtrue_EB_sim_fraction_009), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_009"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_009,fitHisto(h_EoEtrue_EE_sim_fraction_009), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_009"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_01,fitHisto(h_EoEtrue_EB_sim_fraction_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_01"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_01,fitHisto(h_EoEtrue_EE_sim_fraction_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_01"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_02,fitHisto(h_EoEtrue_EB_sim_fraction_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_02"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_02,fitHisto(h_EoEtrue_EE_sim_fraction_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_02"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_03,fitHisto(h_EoEtrue_EB_sim_fraction_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_03"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_03,fitHisto(h_EoEtrue_EE_sim_fraction_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_03"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_04,fitHisto(h_EoEtrue_EB_sim_fraction_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_04"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_04,fitHisto(h_EoEtrue_EE_sim_fraction_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_04"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_05,fitHisto(h_EoEtrue_EB_sim_fraction_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_05"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_05,fitHisto(h_EoEtrue_EE_sim_fraction_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_05")); 
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_01,fitHisto(h_EoEtrue_EB_sim_rechit_diff_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_01"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_01,fitHisto(h_EoEtrue_EE_sim_rechit_diff_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_01"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_02,fitHisto(h_EoEtrue_EB_sim_rechit_diff_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_02"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_02,fitHisto(h_EoEtrue_EE_sim_rechit_diff_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_02"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_03,fitHisto(h_EoEtrue_EB_sim_rechit_diff_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_03"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_03,fitHisto(h_EoEtrue_EE_sim_rechit_diff_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_03")); 
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_04,fitHisto(h_EoEtrue_EB_sim_rechit_diff_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_04"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_04,fitHisto(h_EoEtrue_EE_sim_rechit_diff_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_04")); 
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_05,fitHisto(h_EoEtrue_EB_sim_rechit_diff_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_05"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_05,fitHisto(h_EoEtrue_EE_sim_rechit_diff_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_05")); 
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_06,fitHisto(h_EoEtrue_EB_sim_rechit_diff_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_06"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_06,fitHisto(h_EoEtrue_EE_sim_rechit_diff_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_06"));  
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_07,fitHisto(h_EoEtrue_EB_sim_rechit_diff_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_07"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_07,fitHisto(h_EoEtrue_EE_sim_rechit_diff_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_07"));  
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_08,fitHisto(h_EoEtrue_EB_sim_rechit_diff_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_08"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_08,fitHisto(h_EoEtrue_EE_sim_rechit_diff_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_08"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_09,fitHisto(h_EoEtrue_EB_sim_rechit_diff_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_09"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_09,fitHisto(h_EoEtrue_EE_sim_rechit_diff_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_09"));    
    drawHistFunc(h_EoEtrue_EB_sim_rechit_diff_10,fitHisto(h_EoEtrue_EB_sim_rechit_diff_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_diff_10"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_diff_10,fitHisto(h_EoEtrue_EE_sim_rechit_diff_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_diff_10"));   
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_01,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_01"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_01,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_01"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_02,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_02"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_02,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_02"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_03,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_03"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_03,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_03"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_04,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_04"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_04,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_04"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_05,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_05"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_05,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_05"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_06,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_06"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_06,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_06"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_07,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_07"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_07,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_07"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_08,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_08"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_08,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_08"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_09,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_09"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_09,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_09"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_fraction_10,fitHisto(h_EoEtrue_EB_sim_rechit_fraction_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_fraction_10"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_fraction_10,fitHisto(h_EoEtrue_EE_sim_rechit_fraction_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_fraction_10"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_01,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_01"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_01,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_01"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_02,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_02"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_02,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_02"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_03,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_03"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_03,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_03"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_04,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_04"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_04,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_04"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_05,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_05"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_05,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_05"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_06,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_06"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_06,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_06), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_06"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_07,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_07"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_07,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_07), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_07"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_08,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_08"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_08,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_08), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_08"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_09,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_09"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_09,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_09), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_09"));
    drawHistFunc(h_EoEtrue_EB_global_sim_rechit_fraction_10,fitHisto(h_EoEtrue_EB_global_sim_rechit_fraction_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_global_sim_rechit_fraction_10"));
    drawHistFunc(h_EoEtrue_EE_global_sim_rechit_fraction_10,fitHisto(h_EoEtrue_EE_global_sim_rechit_fraction_10), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_global_sim_rechit_fraction_10"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_combined_fraction01,fitHisto(h_EoEtrue_EB_sim_rechit_combined_fraction01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_combined_fraction01"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_combined_fraction01,fitHisto(h_EoEtrue_EE_sim_rechit_combined_fraction01), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_combined_fraction01"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_combined_fraction02,fitHisto(h_EoEtrue_EB_sim_rechit_combined_fraction02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_combined_fraction02"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_combined_fraction02,fitHisto(h_EoEtrue_EE_sim_rechit_combined_fraction02), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_combined_fraction02"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_combined_fraction03,fitHisto(h_EoEtrue_EB_sim_rechit_combined_fraction03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_combined_fraction03"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_combined_fraction03,fitHisto(h_EoEtrue_EE_sim_rechit_combined_fraction03), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_combined_fraction03"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_combined_fraction04,fitHisto(h_EoEtrue_EB_sim_rechit_combined_fraction04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_combined_fraction04"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_combined_fraction04,fitHisto(h_EoEtrue_EE_sim_rechit_combined_fraction04), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_combined_fraction04"));
    drawHistFunc(h_EoEtrue_EB_sim_rechit_combined_fraction05,fitHisto(h_EoEtrue_EB_sim_rechit_combined_fraction05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_rechit_combined_fraction05"));
    drawHistFunc(h_EoEtrue_EE_sim_rechit_combined_fraction05,fitHisto(h_EoEtrue_EE_sim_rechit_combined_fraction05), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_rechit_combined_fraction05"));
    //drawHistFunc(h_EoEtrue_EB_hgcal_clusterToCalo,fitHisto(h_EoEtrue_EB_hgcal_clusterToCalo), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_hgcal_clusterToCalo"));
    //drawHistFunc(h_EoEtrue_EE_hgcal_clusterToCalo,fitHisto(h_EoEtrue_EE_hgcal_clusterToCalo), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_hgcal_clusterToCalo"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_1MeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_1MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_1MeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_1MeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_1MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_1MeVCut"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_5MeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_5MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_5MeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_5MeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_5MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_5MeVCut"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_10MeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_10MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_10MeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_10MeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_10MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_10MeVCut"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_50MeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_50MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_50MeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_50MeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_50MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_50MeVCut"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_100MeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_100MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_100MeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_100MeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_100MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_100MeVCut"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_500MeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_500MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_500MeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_500MeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_500MeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_500MeVCut"));
    drawHistFunc(h_EoEtrue_EB_sim_fraction_1GeVCut,fitHisto(h_EoEtrue_EB_sim_fraction_1GeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_sim_fraction_1GeVCut"));
    drawHistFunc(h_EoEtrue_EE_sim_fraction_1GeVCut,fitHisto(h_EoEtrue_EE_sim_fraction_1GeVCut), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_sim_fraction_1GeVCut"));
    drawHistFunc(h_EoEtrue_EB_simScore_final_combination,fitHisto(h_EoEtrue_EB_simScore_final_combination), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EB_simScore_final_combination"));
    drawHistFunc(h_EoEtrue_EE_simScore_final_combination,fitHisto(h_EoEtrue_EE_simScore_final_combination), std::string("#sum pfCluster_E_{Reco}/E_{SIM}"), std::string("h_EoEtrue_EE_simScore_final_combination"));
}
