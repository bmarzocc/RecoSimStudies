#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "PhysicsTools/Utilities/macros/setTDRStyle.C"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "RecoSimStudies/Dumpers/interface/SCAlgoValidation.h"
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
    const edm::ParameterSet &histOpt        = process.getParameter<edm::ParameterSet>( "histOpt" );  

    //Read config
    string inputFiles_ = filesOpt.getParameter<string>( "inputFiles" );
    vector<string> inputFiles = split(inputFiles_,',');

    string outputDir_ = filesOpt.getParameter<string>( "outputDir" );
    if( outputDir_ == "" ) outputDir_ = "output/"; 

    int maxEvents_ = filesOpt.getUntrackedParameter<int>( "maxEvents" );
    string fitFunction_ = filesOpt.getParameter<string>( "fitFunction" );
    string superClusterRef_ = filesOpt.getParameter<string>( "superClusterRef" );
    string superClusterVal_ = filesOpt.getParameter<string>( "superClusterVal" );

    std::vector<std::pair<std::string,std::vector<double>>> binOpts = getBinOpts(histOpt); 
    setHistograms(binOpts);
    binOpts = getBinOpts(histOpt); //reset binOpts 

    //useful vars 
    vector<float> simEnergy;
    std::map<int,float> SuperCluster_RecoEnergy_simScore;
    std::map<int,float> DeepSuperCluster_RecoEnergy_simScore;
    std::map<int,double> Calo_SuperCluster_simScore;
    std::map<int,int> Calo_SuperCluster_index;
    std::map<int,double> Calo_DeepSuperCluster_simScore;
    std::map<int,int> Calo_DeepSuperCluster_index;
    std::vector<int> superCluster_seeds;
    std::vector<int> deepSuperCluster_seeds;
    std::vector<int> common_seeds; 

    bool gotoMain = false;
    for(unsigned int iFile=0; iFile<inputFiles.size() && !gotoMain; ++iFile)
    {
       TFile* inFile = TFile::Open(inputFiles[iFile].c_str());
       TTree* inTree = (TTree*)inFile->Get("recosimdumper/caloTree");
       if(inFile)
       {
          cout << "\n--- Reading file: " << inputFiles[iFile].c_str() << endl;

          //setBranches
          setTreeBranches(inTree,superClusterRef_,superClusterVal_);

          // Loop over all entries of the TTree
          for(int entry = 0; entry < inTree->GetEntries(); entry++){
     
             if(entry>maxEvents_ && maxEvents_>0) continue;
             if(entry%1000==0) std::cout << "--- Reading tree = " << entry << std::endl;
             inTree->GetEntry(entry);
             
             //get SC with highest score
             Calo_SuperCluster_simScore.clear();
             for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){
                 Calo_SuperCluster_simScore[iCalo]=0.;
                 Calo_SuperCluster_index[iCalo]=-1;
             }
             Calo_DeepSuperCluster_simScore.clear();
             for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){
                 Calo_DeepSuperCluster_simScore[iCalo]=0.;
                 Calo_DeepSuperCluster_index[iCalo]=-1;
             }
             for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){
                 for(unsigned int iSC=0; iSC<superCluster_simScore->size(); iSC++){
                     if(superCluster_simScore->at(iSC).at(iCalo) > Calo_SuperCluster_simScore[iCalo]){
                        Calo_SuperCluster_simScore[iCalo]=superCluster_simScore->at(iSC).at(iCalo);
                        Calo_SuperCluster_index[iCalo]=iSC;
                     }
                 }
                 for(unsigned int iSC=0; iSC<deepSuperCluster_simScore->size(); iSC++){
                     if(deepSuperCluster_simScore->at(iSC).at(iCalo) > Calo_DeepSuperCluster_simScore[iCalo]){
                        Calo_DeepSuperCluster_simScore[iCalo]=deepSuperCluster_simScore->at(iSC).at(iCalo);
                        Calo_DeepSuperCluster_index[iCalo]=iSC;
                     }
                 }
                 if(Calo_SuperCluster_index[iCalo]>=(int)superCluster_energy->size() && Calo_SuperCluster_index[iCalo]>-1) std::cout << "WARNING SC---> iCalo = " << iCalo << " - " << Calo_SuperCluster_index[iCalo] << " - " << superCluster_energy->size() << std::endl;  
                 if(Calo_DeepSuperCluster_index[iCalo]>=(int)deepSuperCluster_energy->size() && Calo_DeepSuperCluster_index[iCalo]>-1) std::cout << "WARNING DeepSC iCalo = " << iCalo << " - " << Calo_DeepSuperCluster_index[iCalo] << " - " << deepSuperCluster_energy->size() << std::endl;  
             }

             //match SC and DeepSC seeds         
             superCluster_seeds.clear();
             for(unsigned int iPF=0; iPF<superCluster_energy->size(); iPF++)
                 superCluster_seeds.push_back(superCluster_seedIndex->at(iPF));      
             deepSuperCluster_seeds.clear();
             for(unsigned int iPF=0; iPF<deepSuperCluster_energy->size(); iPF++)
                 deepSuperCluster_seeds.push_back(deepSuperCluster_seedIndex->at(iPF));
             common_seeds.clear(); 
             for(unsigned int iPF=0; iPF<superCluster_seeds.size(); iPF++)
                 for(unsigned int iDeepPF=0; iDeepPF<deepSuperCluster_seeds.size(); iDeepPF++)
                     if(deepSuperCluster_seeds.at(iDeepPF)==superCluster_seeds.at(iPF)) common_seeds.push_back(superCluster_seeds.at(iPF));
             
             //fill total and seeMatched histograms
             for(unsigned int iPF=0; iPF<superCluster_simScore_MatchedIndex->size(); iPF++){
                
                if(superCluster_simScore_MatchedIndex->at(iPF)>=0)   
                   SuperCluster_RecoEnergy_simScore[superCluster_simScore_MatchedIndex->at(iPF)]+=superCluster_energy->at(iPF); 
   
                h_Eta_old->Fill(superCluster_eta->at(iPF));
                h_nPFClusters_old->Fill(superCluster_nPFClusters->at(iPF));  
                if(fabs(superCluster_eta->at(iPF))<1.4442){ 
                   h_EtaWidth_EB_old->Fill(superCluster_etaWidth->at(iPF)); 
                   h_Phi_EB_old->Fill(superCluster_phi->at(iPF)); 
                   h_PhiWidth_EB_old->Fill(superCluster_phiWidth->at(iPF));   
             	   h_Energy_EB_old->Fill(superCluster_energy->at(iPF));  
                   h_nPFClusters_EB_old->Fill(superCluster_nPFClusters->at(iPF));   
                   h_R9_EB_old->Fill(superCluster_r9->at(iPF));    
                   h_full5x5_R9_EB_old->Fill(superCluster_full5x5_r9->at(iPF));   
                   h_sigmaIetaIeta_EB_old->Fill(superCluster_sigmaIetaIeta->at(iPF));   
                   h_full5x5_sigmaIetaIeta_EB_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));      
                   h_sigmaIetaIphi_EB_old->Fill(superCluster_sigmaIetaIphi->at(iPF));   
                   h_full5x5_sigmaIetaIphi_EB_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EB_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EB_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));             
                }else if(fabs(superCluster_eta->at(iPF))>1.566){
                   h_EtaWidth_EE_old->Fill(superCluster_etaWidth->at(iPF));
                   h_Phi_EE_old->Fill(superCluster_phi->at(iPF));
                   h_PhiWidth_EE_old->Fill(superCluster_phiWidth->at(iPF)); 
                   h_Energy_EE_old->Fill(superCluster_energy->at(iPF)); 
                   h_nPFClusters_EE_old->Fill(superCluster_nPFClusters->at(iPF));  
                   h_R9_EE_old->Fill(superCluster_r9->at(iPF));             
                   h_full5x5_R9_EE_old->Fill(superCluster_full5x5_r9->at(iPF));
                   h_sigmaIetaIeta_EE_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                   h_full5x5_sigmaIetaIeta_EE_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                   h_sigmaIetaIphi_EE_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                   h_full5x5_sigmaIetaIphi_EE_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EE_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EE_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                } 
         
                bool sameSeed=false;  
                for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
                    if(superCluster_seedIndex->at(iPF) == common_seeds.at(iSeed)){
                       sameSeed=true;
                       break;
                    }    
          
                if(sameSeed==false) continue;
          
                h_Eta_seedMatched_old->Fill(superCluster_eta->at(iPF));
                h_nPFClusters_seedMatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                if(fabs(superCluster_eta->at(iPF))<1.4442){ 
                   h_EtaWidth_EB_seedMatched_old->Fill(superCluster_etaWidth->at(iPF));  
                   h_Phi_EB_seedMatched_old->Fill(superCluster_phi->at(iPF)); 
                   h_PhiWidth_EB_seedMatched_old->Fill(superCluster_phiWidth->at(iPF));  
                   h_Energy_EB_seedMatched_old->Fill(superCluster_energy->at(iPF));  
                   h_nPFClusters_EB_seedMatched_old->Fill(superCluster_nPFClusters->at(iPF));   
                   h_R9_EB_seedMatched_old->Fill(superCluster_r9->at(iPF));    
                   h_full5x5_R9_EB_seedMatched_old->Fill(superCluster_full5x5_r9->at(iPF));   
                   h_sigmaIetaIeta_EB_seedMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                   h_full5x5_sigmaIetaIeta_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                   h_sigmaIetaIphi_EB_seedMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                   h_full5x5_sigmaIetaIphi_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EB_seedMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));              
                }else if(fabs(superCluster_eta->at(iPF))>1.566){
                   h_EtaWidth_EE_seedMatched_old->Fill(superCluster_etaWidth->at(iPF));
                   h_Phi_EE_seedMatched_old->Fill(superCluster_phi->at(iPF));
                   h_PhiWidth_EE_seedMatched_old->Fill(superCluster_phiWidth->at(iPF)); 
                   h_Energy_EE_seedMatched_old->Fill(superCluster_energy->at(iPF)); 
                   h_nPFClusters_EE_seedMatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                   h_R9_EE_seedMatched_old->Fill(superCluster_r9->at(iPF));             
                   h_full5x5_R9_EE_seedMatched_old->Fill(superCluster_full5x5_r9->at(iPF));
                   h_sigmaIetaIeta_EE_seedMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                   h_full5x5_sigmaIetaIeta_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                   h_sigmaIetaIphi_EE_seedMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                   h_full5x5_sigmaIetaIphi_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                   h_sigmaIphiIphi_EE_seedMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                   h_full5x5_sigmaIphiIphi_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                }  
             } 

             for(unsigned int iPF=0; iPF<deepSuperCluster_simScore_MatchedIndex->size(); iPF++){
                 if(deepSuperCluster_simScore_MatchedIndex->at(iPF)>=0)   
                 DeepSuperCluster_RecoEnergy_simScore[deepSuperCluster_simScore_MatchedIndex->at(iPF)]+=deepSuperCluster_energy->at(iPF); 
          
                 h_Eta_new->Fill(deepSuperCluster_eta->at(iPF));
                 h_nPFClusters_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                 if(fabs(deepSuperCluster_eta->at(iPF))<1.4442){ 
                    h_EtaWidth_EB_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                    h_Phi_EB_new->Fill(deepSuperCluster_phi->at(iPF)); 
                    h_PhiWidth_EB_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                    h_Energy_EB_new->Fill(deepSuperCluster_energy->at(iPF));  
                    h_nPFClusters_EB_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                    h_R9_EB_new->Fill(deepSuperCluster_r9->at(iPF));    
                    h_full5x5_R9_EB_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                    h_sigmaIetaIeta_EB_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EB_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EB_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EB_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EB_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EB_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));              
                 }else if(fabs(deepSuperCluster_eta->at(iPF))>1.566){
                    h_EtaWidth_EE_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                    h_Phi_EE_new->Fill(deepSuperCluster_phi->at(iPF));
                    h_PhiWidth_EE_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                    h_Energy_EE_new->Fill(deepSuperCluster_energy->at(iPF)); 
                    h_nPFClusters_EE_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                    h_R9_EE_new->Fill(deepSuperCluster_r9->at(iPF));             
                    h_full5x5_R9_EE_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                    h_sigmaIetaIeta_EE_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EE_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EE_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EE_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EE_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EE_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                 }  
        
                 bool sameSeed=false;  
                 for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
                     if(deepSuperCluster_seedIndex->at(iPF) == common_seeds.at(iSeed)){
                        sameSeed=true;
                        break;
                     }    
          
                 if(sameSeed==false) continue;

                 h_Eta_seedMatched_new->Fill(deepSuperCluster_eta->at(iPF));
                 h_nPFClusters_seedMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                 if(fabs(deepSuperCluster_eta->at(iPF))<1.4442){ 
                    h_EtaWidth_EB_seedMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                    h_Phi_EB_seedMatched_new->Fill(deepSuperCluster_phi->at(iPF)); 
                    h_PhiWidth_EB_seedMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                    h_Energy_EB_seedMatched_new->Fill(deepSuperCluster_energy->at(iPF));  
                    h_nPFClusters_EB_seedMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                    h_R9_EB_seedMatched_new->Fill(deepSuperCluster_r9->at(iPF));    
                    h_full5x5_R9_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                    h_sigmaIetaIeta_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));              
                 }else if(fabs(deepSuperCluster_eta->at(iPF))>1.566){
                    h_EtaWidth_EE_seedMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                    h_Phi_EE_seedMatched_new->Fill(deepSuperCluster_phi->at(iPF));
                    h_PhiWidth_EE_seedMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                    h_Energy_EE_seedMatched_new->Fill(deepSuperCluster_energy->at(iPF)); 
                    h_nPFClusters_EE_seedMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                    h_R9_EE_seedMatched_new->Fill(deepSuperCluster_r9->at(iPF));             
                    h_full5x5_R9_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                    h_sigmaIetaIeta_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                    h_full5x5_sigmaIetaIeta_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                    h_sigmaIetaIphi_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                    h_full5x5_sigmaIetaIphi_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                    h_sigmaIphiIphi_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                    h_full5x5_sigmaIphiIphi_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                 }  
              }

              //fill caloMatched histograms
              for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){

                  h_Eta_Calo_Denum->Fill(caloParticle_simEta->at(iCalo));
                  if(fabs(caloParticle_simEta->at(iCalo))<1.479) h_Et_Calo_EB_Denum->Fill(caloParticle_simEt->at(iCalo));
                  else h_Et_Calo_EE_Denum->Fill(caloParticle_simEt->at(iCalo)); 
  
                  int iPF=Calo_SuperCluster_index.at(iCalo); 
                  if(iPF==-1) continue; 

                  h_Eta_Calo_old->Fill(caloParticle_simEta->at(iCalo)); 
                  if(fabs(caloParticle_simEta->at(iCalo))<1.479) h_Et_Calo_EB_old->Fill(caloParticle_simEt->at(iCalo)); 
                  else h_Et_Calo_EE_old->Fill(caloParticle_simEt->at(iCalo)); 

                  h_Eta_caloMatched_old->Fill(superCluster_eta->at(iPF));
                  h_nPFClusters_caloMatched_old->Fill(superCluster_nPFClusters->at(iPF));
                  prof_EoEtrue_vs_Eta_Calo_old->Fill(caloParticle_simEta->at(iCalo),superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
   
                  int nBins_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[0]; 
                  float min_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[1]; 
                  float max_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[2]; 
                  float delta = fabs(min_Eta-max_Eta)/nBins_Eta;
                  int iBin = int((caloParticle_simEta->at(iCalo)-min_Eta)/delta); 
                  if(caloParticle_simEta->at(iCalo)>=min_Eta && caloParticle_simEta->at(iCalo)<max_Eta)
                     EoEtrue_vs_Eta_Calo_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
        
                  int nBins_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[0]; 
                  float min_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[1]; 
                  float max_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);             
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_Et_Calo_EB_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  nBins_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[0]; 
                  min_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[1]; 
                  max_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_Et_Calo_EE_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));      
        
                  int nBins_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[0]; 
                  float min_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[1]; 
                  float max_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_Energy_Calo_EB_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  nBins_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[0]; 
                  min_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[1]; 
                  max_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_Energy_Calo_EE_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                  int nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[0]; 
                  float min_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[1]; 
                  float max_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_nVtx_EB_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[0]; 
                  min_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[1]; 
                  max_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_nVtx_EE_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    

                  int nBins_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[0]; 
                  float min_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[1]; 
                  float max_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_Rho_EB_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));      
                  nBins_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[0]; 
                  min_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[1]; 
                  max_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_Rho_EE_old[iBin]->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    

                  if(fabs(caloParticle_simEta->at(iCalo))<1.479){
                     h_EoEtrue_EB_old->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EtaWidth_EB_caloMatched_old->Fill(superCluster_etaWidth->at(iPF));  
                     h_Phi_EB_caloMatched_old->Fill(superCluster_phi->at(iPF)); 
                     h_PhiWidth_EB_caloMatched_old->Fill(superCluster_phiWidth->at(iPF));  
                     h_Energy_EB_caloMatched_old->Fill(superCluster_energy->at(iPF));  
                     h_nPFClusters_EB_caloMatched_old->Fill(superCluster_nPFClusters->at(iPF));   
                     h_R9_EB_caloMatched_old->Fill(superCluster_r9->at(iPF));    
                     h_full5x5_R9_EB_caloMatched_old->Fill(superCluster_full5x5_r9->at(iPF));   
                     h_sigmaIetaIeta_EB_caloMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EB_caloMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EB_caloMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF)); 
                     prof_EoEtrue_vs_Et_Calo_EB_old->Fill(caloParticle_simEt->at(iCalo),superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EB_old->Fill(caloParticle_simEnergy->at(iCalo),superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EB_old->Fill(nVtx,superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EB_old->Fill(rho,superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  }else if(fabs(caloParticle_simEta->at(iCalo))>=1.479){
                     h_EoEtrue_EE_old->Fill(superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EtaWidth_EE_caloMatched_old->Fill(superCluster_etaWidth->at(iPF));
                     h_Phi_EE_caloMatched_old->Fill(superCluster_phi->at(iPF));
                     h_PhiWidth_EE_caloMatched_old->Fill(superCluster_phiWidth->at(iPF)); 
                     h_Energy_EE_caloMatched_old->Fill(superCluster_energy->at(iPF)); 
                     h_nPFClusters_EE_caloMatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                     h_R9_EE_caloMatched_old->Fill(superCluster_r9->at(iPF));             
                     h_full5x5_R9_EE_caloMatched_old->Fill(superCluster_full5x5_r9->at(iPF));
                     h_sigmaIetaIeta_EE_caloMatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EE_caloMatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EE_caloMatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));          
                     prof_EoEtrue_vs_Et_Calo_EE_old->Fill(caloParticle_simEt->at(iCalo),superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EE_old->Fill(caloParticle_simEnergy->at(iCalo),superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EE_old->Fill(nVtx,superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EE_old->Fill(rho,superCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  }
               }

               for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++){

                  int iPF=Calo_DeepSuperCluster_index.at(iCalo); 
                  if(iPF==-1) continue; 

                  h_Eta_Calo_new->Fill(caloParticle_simEta->at(iCalo)); 
                  if(fabs(caloParticle_simEta->at(iCalo))<1.479) h_Et_Calo_EB_new->Fill(caloParticle_simEt->at(iCalo)); 
                  else h_Et_Calo_EE_new->Fill(caloParticle_simEt->at(iCalo)); 

                  h_Eta_caloMatched_new->Fill(deepSuperCluster_eta->at(iPF));
                  h_nPFClusters_caloMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));
                  prof_EoEtrue_vs_Eta_Calo_new->Fill(caloParticle_simEta->at(iCalo),deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));

                  int nBins_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[0]; 
                  float min_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[1]; 
                  float max_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[2]; 
                  float delta = fabs(min_Eta-max_Eta)/nBins_Eta;
                  int iBin = int((caloParticle_simEta->at(iCalo)-min_Eta)/delta); 
                  if(caloParticle_simEta->at(iCalo)>=min_Eta && caloParticle_simEta->at(iCalo)<max_Eta)
                     EoEtrue_vs_Eta_Calo_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
         
                  int nBins_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[0]; 
                  float min_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[1]; 
                  float max_Et = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_Et_Calo_EB_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  nBins_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[0]; 
                  min_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[1]; 
                  max_Et = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[2]; 
                  delta = fabs(min_Et-max_Et)/nBins_Et;
                  iBin = int((caloParticle_simEt->at(iCalo)-min_Et)/delta);
                  if(caloParticle_simEt->at(iCalo)>=min_Et && caloParticle_simEt->at(iCalo)<max_Et)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_Et_Calo_EE_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));      
                 
                  int nBins_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[0]; 
                  float min_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[1]; 
                  float max_Energy = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_Energy_Calo_EB_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  nBins_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[0]; 
                  min_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[1]; 
                  max_Energy = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Energy-max_Energy)/nBins_Energy;
                  iBin = int((caloParticle_simEnergy->at(iCalo)-min_Energy)/delta);
                  if(caloParticle_simEnergy->at(iCalo)>=min_Energy && caloParticle_simEnergy->at(iCalo)<max_Energy)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_Energy_Calo_EE_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  
                  int nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[0]; 
                  float min_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[1]; 
                  float max_nVtx = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_nVtx_EB_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  nBins_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[0]; 
                  min_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[1]; 
                  max_nVtx = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_nVtx-max_nVtx)/nBins_nVtx;
                  iBin = int((nVtx-min_nVtx)/delta);
                  if(nVtx>=min_nVtx && nVtx<max_nVtx)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_nVtx_EE_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    

                  int nBins_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[0]; 
                  float min_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[1]; 
                  float max_Rho = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho)
                     if(fabs(caloParticle_simEta->at(iCalo))<1.479) EoEtrue_vs_Rho_EB_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    
                  nBins_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[0]; 
                  min_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[1]; 
                  max_Rho = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[2];      
                  delta = fabs(min_Rho-max_Rho)/nBins_Rho;
                  iBin = int((rho-min_Rho)/delta);
                  if(rho>=min_Rho && rho<max_Rho)
                     if(fabs(caloParticle_simEta->at(iCalo))>1.479) EoEtrue_vs_Rho_EE_new[iBin]->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));    

                  if(fabs(caloParticle_simEta->at(iCalo))<1.479){
                     h_EoEtrue_EB_new->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EtaWidth_EB_caloMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                     h_Phi_EB_caloMatched_new->Fill(deepSuperCluster_phi->at(iPF)); 
                     h_PhiWidth_EB_caloMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                     h_Energy_EB_caloMatched_new->Fill(deepSuperCluster_energy->at(iPF));  
                     h_nPFClusters_EB_caloMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                     h_R9_EB_caloMatched_new->Fill(deepSuperCluster_r9->at(iPF));    
                     h_full5x5_R9_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                     h_sigmaIetaIeta_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF)); 
                     prof_EoEtrue_vs_Et_Calo_EB_new->Fill(caloParticle_simEt->at(iCalo),deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EB_new->Fill(caloParticle_simEnergy->at(iCalo),deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EB_new->Fill(nVtx,deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EB_new->Fill(rho,deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  }else if(fabs(caloParticle_simEta->at(iCalo))>=1.479){
                     h_EoEtrue_EE_new->Fill(deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     h_EtaWidth_EE_caloMatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                     h_Phi_EE_caloMatched_new->Fill(deepSuperCluster_phi->at(iPF));
                     h_PhiWidth_EE_caloMatched_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                     h_Energy_EE_caloMatched_new->Fill(deepSuperCluster_energy->at(iPF)); 
                     h_nPFClusters_EE_caloMatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                     h_R9_EE_caloMatched_new->Fill(deepSuperCluster_r9->at(iPF));             
                     h_full5x5_R9_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                     h_sigmaIetaIeta_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                     h_full5x5_sigmaIetaIeta_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                     h_sigmaIetaIphi_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                     h_full5x5_sigmaIetaIphi_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                     h_sigmaIphiIphi_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                     h_full5x5_sigmaIphiIphi_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));   
                     prof_EoEtrue_vs_Et_Calo_EE_new->Fill(caloParticle_simEt->at(iCalo),deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo)); 
                     prof_EoEtrue_vs_Energy_Calo_EE_new->Fill(caloParticle_simEnergy->at(iCalo),deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));  
                     prof_EoEtrue_vs_nVtx_EE_new->Fill(nVtx,deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                     prof_EoEtrue_vs_Rho_EE_new->Fill(rho,deepSuperCluster_energy->at(iPF)/caloParticle_simEnergy->at(iCalo));
                  }
               }

               for(unsigned int iPF=0; iPF<superCluster_energy->size(); iPF++){

                   std::map<int,int>::iterator it = Calo_SuperCluster_index.find(iPF); 
                   if(it != Calo_SuperCluster_index.end()) continue;   
 
                   h_Eta_caloUnmatched_old->Fill(superCluster_eta->at(iPF));
                   h_nPFClusters_caloUnmatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                   if(fabs(superCluster_eta->at(iPF))<1.4442){ 
                      h_EtaWidth_EB_caloUnmatched_old->Fill(superCluster_etaWidth->at(iPF));  
                      h_Phi_EB_caloUnmatched_old->Fill(superCluster_phi->at(iPF)); 
                      h_PhiWidth_EB_caloUnmatched_old->Fill(superCluster_phiWidth->at(iPF));  
                      h_Energy_EB_caloUnmatched_old->Fill(superCluster_energy->at(iPF));  
                      h_nPFClusters_EB_caloUnmatched_old->Fill(superCluster_nPFClusters->at(iPF));   
                      h_R9_EB_caloUnmatched_old->Fill(superCluster_r9->at(iPF));    
                      h_full5x5_R9_EB_caloUnmatched_old->Fill(superCluster_full5x5_r9->at(iPF));   
                      h_sigmaIetaIeta_EB_caloUnmatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EB_caloUnmatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EB_caloUnmatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));              
                   }else if(fabs(superCluster_eta->at(iPF))>1.566){
                      h_EtaWidth_EE_caloUnmatched_old->Fill(superCluster_etaWidth->at(iPF));
             	      h_Phi_EE_caloUnmatched_old->Fill(superCluster_phi->at(iPF));
                      h_PhiWidth_EE_caloUnmatched_old->Fill(superCluster_phiWidth->at(iPF)); 
                      h_Energy_EE_caloUnmatched_old->Fill(superCluster_energy->at(iPF)); 
                      h_nPFClusters_EE_caloUnmatched_old->Fill(superCluster_nPFClusters->at(iPF));  
                      h_R9_EE_caloUnmatched_old->Fill(superCluster_r9->at(iPF));             
                      h_full5x5_R9_EE_caloUnmatched_old->Fill(superCluster_full5x5_r9->at(iPF));
                      h_sigmaIetaIeta_EE_caloUnmatched_old->Fill(superCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EE_caloUnmatched_old->Fill(superCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EE_caloUnmatched_old->Fill(superCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIphiIphi->at(iPF));                             
                   } 
               }

               for(unsigned int iPF=0; iPF<deepSuperCluster_energy->size(); iPF++){

                   std::map<int,int>::iterator it = Calo_DeepSuperCluster_index.find(iPF); 
                   if(it != Calo_DeepSuperCluster_index.end()) continue;  

                   h_Eta_caloUnmatched_new->Fill(deepSuperCluster_eta->at(iPF));
                   h_nPFClusters_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  

                   if(fabs(deepSuperCluster_eta->at(iPF))<1.4442){ 
                      h_EtaWidth_EB_caloUnmatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));  
                      h_Phi_EB_caloUnmatched_new->Fill(deepSuperCluster_phi->at(iPF)); 
                      h_PhiWidth_EB_caloUnmatched_new->Fill(deepSuperCluster_phiWidth->at(iPF));  
                      h_Energy_EB_caloUnmatched_new->Fill(deepSuperCluster_energy->at(iPF));  
                      h_nPFClusters_EB_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));   
                      h_R9_EB_caloUnmatched_new->Fill(deepSuperCluster_r9->at(iPF));    
                      h_full5x5_R9_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));   
                      h_sigmaIetaIeta_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));              
                   }else if(fabs(deepSuperCluster_eta->at(iPF))>1.566){
                      h_EtaWidth_EE_caloUnmatched_new->Fill(deepSuperCluster_etaWidth->at(iPF));
                      h_Phi_EE_caloUnmatched_new->Fill(deepSuperCluster_phi->at(iPF));
                      h_PhiWidth_EE_caloUnmatched_new->Fill(deepSuperCluster_phiWidth->at(iPF)); 
                      h_Energy_EE_caloUnmatched_new->Fill(deepSuperCluster_energy->at(iPF)); 
                      h_nPFClusters_EE_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters->at(iPF));  
                      h_R9_EE_caloUnmatched_new->Fill(deepSuperCluster_r9->at(iPF));             
                      h_full5x5_R9_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_r9->at(iPF));
                      h_sigmaIetaIeta_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIeta->at(iPF));    
                      h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta->at(iPF));       
                      h_sigmaIetaIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIphi->at(iPF));    
                      h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi->at(iPF));    
                      h_sigmaIphiIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIphiIphi->at(iPF));    
                      h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi->at(iPF));  
                   }           
                }

                simEnergy.clear();
                SuperCluster_RecoEnergy_simScore.clear(); 
                DeepSuperCluster_RecoEnergy_simScore.clear();       
          }
       }
    }
   
    setEfficiencies();
    drawPlots(fitFunction_,superClusterRef_,superClusterVal_);
    
}
