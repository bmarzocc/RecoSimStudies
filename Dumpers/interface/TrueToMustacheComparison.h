#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include <algorithm> 
#include <iostream>
#include <utility>
#include <string>
#include <vector>

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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "RecoSimStudies/Dumpers/interface/CruijffPdf.h"
#include "RecoSimStudies/Dumpers/interface/DoubleCBPdf.h"

#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"

#include <boost/algorithm/string/replace.hpp>

using namespace std;
using namespace edm;

//DEFINE HISTOGRAMS 
std::vector<std::vector<TH2F*>> dEtadPhi_Calo_vs_seedEt_seedEta;
std::vector<std::vector<TH2F*>> dEtadPhi_Calo_vs_seedEt_seedEta_inMustache;
std::vector<std::vector<TH2F*>> dEtadPhi_Calo_vs_seedEt_seedEta_notMustache;
std::vector<std::vector<TH2F*>> dEtadPhi_Mustache_vs_seedEt_seedEta;
std::vector<std::vector<TH2F*>> dEtadPhi_Mustache_vs_seedEt_seedEta_inCalo;
std::vector<std::vector<TH2F*>> dEtadPhi_Mustache_vs_seedEt_seedEta_notCalo;

std::vector<std::vector<TProfile2D*>> HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache;
std::vector<std::vector<TProfile2D*>> HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache;
std::vector<std::vector<TProfile2D*>> SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache;
std::vector<std::vector<TProfile2D*>> SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache;
std::vector<std::vector<TProfile2D*>> SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache;
std::vector<std::vector<TProfile2D*>> SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache;
std::vector<std::vector<TProfile2D*>> RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache;
std::vector<std::vector<TProfile2D*>> RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache;
std::vector<std::vector<TProfile2D*>> RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache;
std::vector<std::vector<TProfile2D*>> RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache;

//DEFINE BRANCHES
vector<float> *caloParticle_simEnergy;
vector<vector<double> > *pfCluster_dR_genScore;
vector<vector<double> > *pfCluster_dR_simScore;
vector<vector<double> > *pfCluster_sim_fraction_old;
vector<vector<double> > *pfCluster_n_shared_xtals;
vector<vector<double> > *pfCluster_sim_fraction;
vector<vector<double> > *pfCluster_sim_fraction_withFraction;
vector<vector<double> > *pfCluster_sim_fraction_1MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_5MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_10MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_50MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_100MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_500MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_1GeVCut;
vector<vector<double> > *pfCluster_sim_rechit_diff;
vector<vector<double> > *pfCluster_sim_rechit_fraction;
vector<vector<double> > *pfCluster_global_sim_rechit_fraction;
vector<vector<float> > *pfClusterHit_fraction;
vector<vector<float> > *pfClusterHit_eta;
vector<vector<float> > *pfClusterHit_phi;
vector<float> *pfCluster_rawEnergy;
vector<float> *pfCluster_energy;
vector<float> *pfCluster_eta;
vector<float> *pfCluster_phi;
vector<vector<int> > *superCluster_pfClustersIndex;
vector<int> *superCluster_sim_fraction_MatchedIndex;
vector<int> *superCluster_seedIndex;
vector<float> *superCluster_rawEnergy;
vector<int> *superCluster_iz;

TBranch *b_caloParticle_simEnergy;
TBranch *b_pfCluster_dR_genScore;
TBranch *b_pfCluster_dR_simScore;
TBranch *b_pfCluster_sim_fraction_old;
TBranch *b_pfCluster_n_shared_xtals;
TBranch *b_pfCluster_sim_fraction;
TBranch *b_pfCluster_sim_fraction_withFraction;
TBranch *b_pfCluster_sim_fraction_1MeVCut;
TBranch *b_pfCluster_sim_fraction_5MeVCut;
TBranch *b_pfCluster_sim_fraction_10MeVCut;
TBranch *b_pfCluster_sim_fraction_50MeVCut;
TBranch *b_pfCluster_sim_fraction_100MeVCut;
TBranch *b_pfCluster_sim_fraction_500MeVCut;
TBranch *b_pfCluster_sim_fraction_1GeVCut;
TBranch *b_pfCluster_sim_rechit_diff;
TBranch *b_pfCluster_sim_rechit_fraction;
TBranch *b_pfCluster_global_sim_rechit_fraction;
TBranch *b_pfClusterHit_fraction;
TBranch *b_pfClusterHit_eta;
TBranch *b_pfClusterHit_phi;
TBranch *b_pfCluster_rawEnergy;
TBranch *b_pfCluster_energy;
TBranch *b_pfCluster_eta;
TBranch *b_pfCluster_phi;
TBranch *b_superCluster_pfClustersIndex;
TBranch *b_superCluster_sim_fraction_MatchedIndex;
TBranch *b_superCluster_seedIndex;
TBranch *b_superCluster_rawEnergy;
TBranch *b_superCluster_iz;

//setTreeBranches
void setTreeBranches(TTree* tree)
{
   tree->SetBranchStatus("*",0);   
   tree->SetBranchStatus("caloParticle_simEnergy",1);       
   tree->SetBranchStatus("pfCluster_sim_fraction",1);
   tree->SetBranchStatus("pfCluster_sim_fraction_withFraction",1); 
   tree->SetBranchStatus("pfClusterHit_fraction",1);  
   tree->SetBranchStatus("pfClusterHit_eta",1); 
   tree->SetBranchStatus("pfClusterHit_phi",1); 
   tree->SetBranchStatus("pfCluster_rawEnergy",1);  
   tree->SetBranchStatus("pfCluster_energy",1);  
   tree->SetBranchStatus("pfCluster_eta",1); 
   tree->SetBranchStatus("pfCluster_phi",1); 
   tree->SetBranchStatus("superCluster_pfClustersIndex",1); 
   tree->SetBranchStatus("superCluster_sim_fraction_MatchedIndex",1); 
   tree->SetBranchStatus("superCluster_seedIndex",1);
   tree->SetBranchStatus("superCluster_rawEnergy",1); 
   tree->SetBranchStatus("superCluster_iz",1); 
  
   tree->SetBranchAddress("caloParticle_simEnergy", &caloParticle_simEnergy, &b_caloParticle_simEnergy); 
   tree->SetBranchAddress("pfCluster_sim_fraction", &pfCluster_sim_fraction, &b_pfCluster_sim_fraction);
   tree->SetBranchAddress("pfCluster_sim_fraction_withFraction", &pfCluster_sim_fraction_withFraction, &b_pfCluster_sim_fraction_withFraction);
   tree->SetBranchAddress("pfClusterHit_fraction", &pfClusterHit_fraction, &b_pfClusterHit_fraction); 
   tree->SetBranchAddress("pfClusterHit_eta", &pfClusterHit_eta, &b_pfClusterHit_eta); 
   tree->SetBranchAddress("pfClusterHit_phi", &pfClusterHit_phi, &b_pfClusterHit_phi);  
   tree->SetBranchAddress("pfCluster_rawEnergy", &pfCluster_rawEnergy, &b_pfCluster_rawEnergy);
   tree->SetBranchAddress("pfCluster_energy", &pfCluster_energy, &b_pfCluster_energy); 
   tree->SetBranchAddress("pfCluster_eta", &pfCluster_eta, &b_pfCluster_eta);
   tree->SetBranchAddress("pfCluster_phi", &pfCluster_phi, &b_pfCluster_phi); 
   tree->SetBranchAddress("superCluster_pfClustersIndex", &superCluster_pfClustersIndex, &b_superCluster_pfClustersIndex);
   tree->SetBranchAddress("superCluster_sim_fraction_MatchedIndex", &superCluster_sim_fraction_MatchedIndex, &b_superCluster_sim_fraction_MatchedIndex);
   tree->SetBranchAddress("superCluster_seedIndex", &superCluster_seedIndex, &b_superCluster_seedIndex);
   tree->SetBranchAddress("superCluster_rawEnergy", &superCluster_rawEnergy, &b_superCluster_rawEnergy);
   tree->SetBranchAddress("superCluster_iz", &superCluster_iz, &b_superCluster_iz); 
}

//split
std::vector<std::string> split(const string &text, char sep)
{
   vector<string> tokens;
   size_t start = 0, end = 0;
   while ((end = text.find(sep, start)) != string::npos) {
     tokens.push_back(text.substr(start, end - start));
     start = end + 1;
   }
   tokens.push_back(text.substr(start));
   return tokens;
}

//set histograms
void setHistograms(std::vector<std::string> etCuts, std::vector<std::string> etaCuts)
{
   
   dEtadPhi_Calo_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Calo_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Calo_vs_seedEt_seedEta[iBin][jBin] = new TH2F(std::string("dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   } 

   dEtadPhi_Calo_vs_seedEt_seedEta_inMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin][jBin] = new TH2F(std::string("dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), std::string("dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   }

   dEtadPhi_Calo_vs_seedEt_seedEta_notMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin][jBin] = new TH2F(std::string("dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), std::string("dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   }  

   dEtadPhi_Mustache_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Mustache_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_seedEt_seedEta[iBin][jBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   } 

   dEtadPhi_Mustache_vs_seedEt_seedEta_inCalo.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Mustache_vs_seedEt_seedEta_inCalo[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_seedEt_seedEta_inCalo[iBin][jBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inCalo").c_str(), std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inCalo").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   }

   dEtadPhi_Mustache_vs_seedEt_seedEta_notCalo.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Mustache_vs_seedEt_seedEta_notCalo[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_seedEt_seedEta_notCalo[iBin][jBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notCalo").c_str(), std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notCalo").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   }

   HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("HitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("HitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   } 

   HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin][jBin] = new TProfile2D(std::string("HitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), std::string("HitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   }

   HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin][jBin] = new TProfile2D(std::string("HitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), std::string("HitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   }

   SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("SimFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("SimFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   } 

   SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin][jBin] = new TProfile2D(std::string("SimFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), std::string("SimFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   }

   SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin][jBin] = new TProfile2D(std::string("SimFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), std::string("SimFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   }

   SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   } 

   SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin][jBin] = new TProfile2D(std::string("SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), std::string("SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   }

   SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin][jBin] = new TProfile2D(std::string("SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), std::string("SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
        
       }
   }

   RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("RecoToCalo_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("RecoToCalo_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
        
       }
   } 

   RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin][jBin] = new TProfile2D(std::string("RecoToCalo_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), std::string("RecoToCalo_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
        
       }
   }

   RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin][jBin] = new TProfile2D(std::string("RecoToCalo_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), std::string("RecoToCalo_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
        
       }
   }

   RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("RecoToSeed_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("RecoToSeed_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
        
       }
   } 

   RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[iBin][jBin] = new TProfile2D(std::string("RecoToSeed_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), std::string("RecoToSeed_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
        
       }
   }

   RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[iBin][jBin] = new TProfile2D(std::string("RecoToSeed_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), std::string("RecoToSeed_dEtadPhi_Calo_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
        
       }
   } 
}

int getMatchedIndex(std::vector<std::vector<double>>* score, double selection, bool useMax, std::vector<std::vector<std::vector<double>>> scoreSelMax, std::vector<double> selectionMax, std::vector<std::vector<std::vector<double>>> scoreSelMin, std::vector<double> selectionMin, int iCl)
{
   int matchedIndex = -1; 
   if(!useMax){ 
      std::replace(score->at(iCl).begin(),score->at(iCl).end(), -999., 999.);
      if(std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::min_element(score->at(iCl).begin(),score->at(iCl).end()) - score->at(iCl).begin();
   }else{
      std::replace(score->at(iCl).begin(),score->at(iCl).end(), 999., -999.); 
      if(std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iCl).begin(),score->at(iCl).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::max_element(score->at(iCl).begin(),score->at(iCl).end()) - score->at(iCl).begin();
   }

   if(matchedIndex==-1) return -1;
   
   //std::cout << "getMatchedIndex - Max - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMax.size() << " - " << selectionMax.size() << std::endl;
   //std::cout << "getMatchedIndex - Min - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMin.size() << " - " << selectionMin.size() << std::endl;
    
   bool passSelection = true;
   for(unsigned int iSelMax=0; iSelMax < scoreSelMax.size(); iSelMax++)
       if(scoreSelMax.at(iSelMax).at(iCl).at(matchedIndex) < selectionMax.at(iSelMax)) passSelection = false;
      
   for(unsigned int iSelMin=0; iSelMin < scoreSelMin.size(); iSelMin++)
       if(scoreSelMin.at(iSelMin).at(iCl).at(matchedIndex) > selectionMin.at(iSelMin)) passSelection = false;
 
   if(useMax && score->at(iCl).at(matchedIndex) > selection && passSelection) return matchedIndex;
   if(!useMax && score->at(iCl).at(matchedIndex) < selection && passSelection) return matchedIndex; 
   
   return -1; 
}


void drawH2(TH2F* h2, std::string xtitle, std::string ytitle, string ztitle, std::string Name)
{
   gStyle->SetOptStat(0000); 

   h2->GetXaxis()->SetTitle(xtitle.c_str()); 
   h2->GetYaxis()->SetTitle(ytitle.c_str()); 
   h2->GetZaxis()->SetTitle(ztitle.c_str());

   TCanvas* c = new TCanvas();
   h2->Draw("COLZ");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");

   gStyle->SetOptStat(1110); 
}

void drawProfile2(TProfile2D* h2, std::string xtitle, std::string ytitle, string ztitle, std::string Name, float z_min=0., float z_max=1.)
{
   gStyle->SetOptStat(0000); 

   h2->GetXaxis()->SetTitle(xtitle.c_str()); 
   h2->GetYaxis()->SetTitle(ytitle.c_str()); 
   h2->GetZaxis()->SetTitle(ztitle.c_str());
   if(z_min!=-1. && z_max!=-1.) h2->GetZaxis()->SetRangeUser(z_min,z_max);

   TCanvas* c = new TCanvas();
   c->SetLogz(); 
   h2->Draw("PROFCOLZ");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");

   gStyle->SetOptStat(1110); 
}

//draw plots
void drawPlots(std::vector<std::string> etCuts, std::vector<std::string> etaCuts, std::string outDir_)
{
   for(unsigned int etBin=0; etBin<etCuts.size()-1; etBin++){
       for(unsigned int etaBin=0; etaBin<etaCuts.size()-1; etaBin++){
           std::string output_dir = outDir_+"/seedEt_"+etCuts.at(etBin)+"_"+etCuts.at(etBin+1)+"_seedEta_"+etaCuts.at(etaBin)+"_"+etaCuts.at(etaBin+1)+"/";
           system(string("mkdir -p "+output_dir).c_str()); 
           system(string("cp index.php "+output_dir).c_str()); 
           drawH2(dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Calo"));
           drawH2(dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Calo_inMustache"));
           drawH2(dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Calo_notMustache"));
           drawH2(dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache")); 
           drawH2(dEtadPhi_Mustache_vs_seedEt_seedEta_inCalo[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache_inCalo")); 
           drawH2(dEtadPhi_Mustache_vs_seedEt_seedEta_notCalo[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache_notCalo")); 
           drawProfile2(HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"HitFraction_dEtadPhi_Calo"));
           drawProfile2(HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"HitFraction_dEtadPhi_Calo_inMustache"));
           drawProfile2(HitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"HitFraction_dEtadPhi_Calo_notMustache"));
           drawProfile2(SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"SimFraction_dEtadPhi_Calo"));
           drawProfile2(SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"SimFraction_dEtadPhi_Calo_inMustache"));
           drawProfile2(SimFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"SimFraction_dEtadPhi_Calo_notMustache"));
           drawProfile2(SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"SimFraction_withHitFraction_dEtadPhi_Calo"));
           drawProfile2(SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"SimFraction_withHitFraction_dEtadPhi_Calo_inMustache"));
           drawProfile2(SimFraction_withHitFraction_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"SimFraction_withHitFraction_dEtadPhi_Calo_notMustache"));
           drawProfile2(RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoToCalo_dEtadPhi_Calo"), -1., -1.);
           drawProfile2(RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoToCalo_dEtadPhi_Calo_inMustache"), -1., 1.); 
           drawProfile2(RecoToCalo_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoToCalo_dEtadPhi_Calo_notMustache"), -1., 1.);
           drawProfile2(RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoToSeed_dEtadPhi_Calo"), 0.0001, 1.);
           drawProfile2(RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_inMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoToSeed_dEtadPhi_Calo_inMustache"), 0.0001, 1.); 
           drawProfile2(RecoToSeed_dEtadPhi_Calo_vs_seedEt_seedEta_notMustache[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoToSeed_dEtadPhi_Calo_notMustache"), 0.0001, 1.);
       }  
   } 
}

