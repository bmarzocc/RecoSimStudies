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
std::vector<std::vector<TH1F*>> MustOEtrue_vs_scoreThreshold_seedEt_EB;
std::vector<std::vector<TH1F*>> MustOEtrue_vs_scoreThreshold_seedEt_EE;
std::vector<std::vector<TH1F*>> MustOEtrue_vs_scoreThreshold_seedEta;
std::vector<std::vector<std::vector<TH1F*>>> MustOEtrue_vs_scoreThreshold_seedEta_seedEt;

std::vector<std::vector<TH1F*>> EtrueOEtrue_vs_scoreThreshold_seedEt_EB;
std::vector<std::vector<TH1F*>> EtrueOEtrue_vs_scoreThreshold_seedEt_EE;
std::vector<std::vector<TH1F*>> EtrueOEtrue_vs_scoreThreshold_seedEta;
std::vector<std::vector<std::vector<TH1F*>>> EtrueOEtrue_vs_scoreThreshold_seedEta_seedEt;

std::vector<std::vector<TH1F*>> ErecoOEtrue_vs_scoreThreshold_seedEt_EB;
std::vector<std::vector<TH1F*>> ErecoOEtrue_vs_scoreThreshold_seedEt_EE;
std::vector<std::vector<TH1F*>> ErecoOEtrue_vs_scoreThreshold_seedEta;
std::vector<std::vector<std::vector<TH1F*>>> ErecoOEtrue_vs_scoreThreshold_seedEta_seedEt;

std::vector<std::vector<std::vector<TH2F*>>> dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt;
std::vector<std::vector<std::vector<TH2F*>>> dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache;
std::vector<std::vector<std::vector<TH2F*>>> dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache;
std::vector<std::vector<std::vector<TH2F*>>> dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt;
std::vector<std::vector<std::vector<TH2F*>>> dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_inCalo;
std::vector<std::vector<std::vector<TH2F*>>> dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_notCalo;

std::vector<std::vector<std::vector<TProfile2D*>>> TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt;
std::vector<std::vector<std::vector<TProfile2D*>>> TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache;
std::vector<std::vector<std::vector<TProfile2D*>>> TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache;
std::vector<std::vector<std::vector<TProfile2D*>>> RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt;
std::vector<std::vector<std::vector<TProfile2D*>>> RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache;
std::vector<std::vector<std::vector<TProfile2D*>>> RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache;

//DEFINE BRANCHES
vector<float> *caloParticle_simEnergy;
vector<vector<double> > *pfCluster_dR_genScore;
vector<vector<double> > *pfCluster_dR_simScore;
vector<vector<double> > *pfCluster_sim_fraction_old;
vector<vector<double> > *pfCluster_n_shared_xtals;
vector<vector<double> > *pfCluster_sim_fraction;
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
vector<float> *pfCluster_rawEnergy;
vector<float> *pfCluster_energy;
vector<float> *pfCluster_eta;
vector<float> *pfCluster_phi;
vector<vector<int> > *superCluster_pfClustersIndex;
vector<int> *superCluster_seedIndex;
vector<float> *superCluster_rawEnergy;
vector<int> *superCluster_iz;

TBranch *b_caloParticle_simEnergy;
TBranch *b_pfCluster_dR_genScore;
TBranch *b_pfCluster_dR_simScore;
TBranch *b_pfCluster_sim_fraction_old;
TBranch *b_pfCluster_n_shared_xtals;
TBranch *b_pfCluster_sim_fraction;
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
TBranch *b_pfCluster_rawEnergy;
TBranch *b_pfCluster_energy;
TBranch *b_pfCluster_eta;
TBranch *b_pfCluster_phi;
TBranch *b_superCluster_pfClustersIndex;
TBranch *b_superCluster_seedIndex;
TBranch *b_superCluster_rawEnergy;
TBranch *b_superCluster_iz;

//setTreeBranches
void setTreeBranches(TTree* tree, string matching_)
{
   tree->SetBranchStatus("*",0);   
   tree->SetBranchStatus("caloParticle_simEnergy",1);       
   if(matching_ == "dR_genScore") tree->SetBranchStatus("pfCluster_dR_genScore",1); 
   else if(matching_ == "dR_simScore") tree->SetBranchStatus("pfCluster_dR_simScore",1); 
   else if(matching_ == "sim_fraction_old") tree->SetBranchStatus("pfCluster_sim_fraction_old",1); 
   else if(matching_ == "n_shared_xtals") tree->SetBranchStatus("pfCluster_n_shared_xtals",1); 
   else if(matching_ == "sim_fraction") tree->SetBranchStatus("pfCluster_sim_fraction",1); 
   else if(matching_ == "sim_fraction_1MeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_1MeVCut",1); 
   else if(matching_ == "sim_fraction_5MeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_5MeVCut",1); 
   else if(matching_ == "sim_fraction_10MeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_10MeVCut",1); 
   else if(matching_ == "sim_fraction_50MeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_50MeVCut",1);  
   else if(matching_ == "sim_fraction_100MeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_100MeVCut",1); 
   else if(matching_ == "sim_fraction_500MeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_500MeVCut",1); 
   else if(matching_ == "sim_fraction_1GeVCut") tree->SetBranchStatus("pfCluster_sim_fraction_1GeVCut",1); 
   else if(matching_ == "sim_rechit_diff") tree->SetBranchStatus("pfCluster_sim_rechit_diff",1);  
   else if(matching_ == "sim_rechit_fraction") tree->SetBranchStatus("pfCluster_sim_rechit_fraction",1); 
   else if(matching_ == "global_sim_rechit_fraction") tree->SetBranchStatus("pfCluster_global_sim_rechit_fraction",1); 
   tree->SetBranchStatus("pfCluster_rawEnergy",1); 
   tree->SetBranchStatus("pfCluster_energy",1);  
   tree->SetBranchStatus("pfCluster_eta",1); 
   tree->SetBranchStatus("pfCluster_phi",1); 
   tree->SetBranchStatus("superCluster_pfClustersIndex",1); 
   tree->SetBranchStatus("superCluster_seedIndex",1);
   tree->SetBranchStatus("superCluster_rawEnergy",1); 
   tree->SetBranchStatus("superCluster_iz",1); 
  
   tree->SetBranchAddress("caloParticle_simEnergy", &caloParticle_simEnergy, &b_caloParticle_simEnergy); 
   if(matching_ == "dR_genScore") tree->SetBranchAddress("pfCluster_dR_genScore", &pfCluster_dR_genScore, &b_pfCluster_dR_genScore);
   else if(matching_ == "dR_simScore") tree->SetBranchAddress("pfCluster_dR_simScore", &pfCluster_dR_simScore, &b_pfCluster_dR_simScore);
   else if(matching_ == "sim_fraction_old") tree->SetBranchAddress("pfCluster_sim_fraction_old", &pfCluster_sim_fraction_old, &b_pfCluster_sim_fraction_old);
   else if(matching_ == "n_shared_xtals") tree->SetBranchAddress("pfCluster_n_shared_xtals", &pfCluster_n_shared_xtals, &b_pfCluster_n_shared_xtals);
   else if(matching_ == "sim_fraction") tree->SetBranchAddress("pfCluster_sim_fraction", &pfCluster_sim_fraction, &b_pfCluster_sim_fraction);
   else if(matching_ == "sim_fraction_1MeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_1MeVCut", &pfCluster_sim_fraction_1MeVCut, &b_pfCluster_sim_fraction_1MeVCut);
   else if(matching_ == "sim_fraction_5MeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_5MeVCut", &pfCluster_sim_fraction_5MeVCut, &b_pfCluster_sim_fraction_5MeVCut);
   else if(matching_ == "sim_fraction_10MeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_10MeVCut", &pfCluster_sim_fraction_10MeVCut, &b_pfCluster_sim_fraction_10MeVCut);
   else if(matching_ == "sim_fraction_50MeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_50MeVCut", &pfCluster_sim_fraction_50MeVCut, &b_pfCluster_sim_fraction_50MeVCut);
   else if(matching_ == "sim_fraction_100MeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_100MeVCut", &pfCluster_sim_fraction_100MeVCut, &b_pfCluster_sim_fraction_100MeVCut);
   else if(matching_ == "sim_fraction_500MeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_500MeVCut", &pfCluster_sim_fraction_500MeVCut, &b_pfCluster_sim_fraction_500MeVCut);
   else if(matching_ == "sim_fraction_1GeVCut") tree->SetBranchAddress("pfCluster_sim_fraction_1GeVCut", &pfCluster_sim_fraction_1GeVCut, &b_pfCluster_sim_fraction_1GeVCut);
   else if(matching_ == "sim_rechit_diff") tree->SetBranchAddress("pfCluster_sim_rechit_diff", &pfCluster_sim_rechit_diff, &b_pfCluster_sim_rechit_diff);
   else if(matching_ == "sim_rechit_fraction") tree->SetBranchAddress("pfCluster_sim_rechit_fraction", &pfCluster_sim_rechit_fraction, &b_pfCluster_sim_rechit_fraction);
   else if(matching_ == "global_sim_rechit_fraction") tree->SetBranchAddress("pfCluster_global_sim_rechit_fraction", &pfCluster_global_sim_rechit_fraction, &b_pfCluster_global_sim_rechit_fraction);
   tree->SetBranchAddress("pfCluster_rawEnergy", &pfCluster_rawEnergy, &b_pfCluster_rawEnergy);
   tree->SetBranchAddress("pfCluster_energy", &pfCluster_energy, &b_pfCluster_energy); 
   tree->SetBranchAddress("pfCluster_eta", &pfCluster_eta, &b_pfCluster_eta);
   tree->SetBranchAddress("pfCluster_phi", &pfCluster_phi, &b_pfCluster_phi); 
   tree->SetBranchAddress("superCluster_pfClustersIndex", &superCluster_pfClustersIndex, &b_superCluster_pfClustersIndex);
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
void setHistograms(std::vector<std::string> scoreCuts, std::vector<std::string> etCuts, std::vector<std::string> etaCuts)
{
   
   MustOEtrue_vs_scoreThreshold_seedEt_EB.resize(etCuts.size()-1); 
   MustOEtrue_vs_scoreThreshold_seedEt_EE.resize(etCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       MustOEtrue_vs_scoreThreshold_seedEt_EB[iBin].resize(scoreCuts.size());
       MustOEtrue_vs_scoreThreshold_seedEt_EE[iBin].resize(scoreCuts.size());
       for(unsigned int jBin=0; jBin<scoreCuts.size(); jBin++)
       {
           boost::replace_all(scoreCuts.at(jBin), ".", "_");
           MustOEtrue_vs_scoreThreshold_seedEt_EB[iBin][jBin] = new TH1F(std::string("MustOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EB").c_str(), std::string("MustOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EB").c_str(), 10000, 0., 10.);
           MustOEtrue_vs_scoreThreshold_seedEt_EE[iBin][jBin] = new TH1F(std::string("MustOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EE").c_str(), std::string("MustOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EE").c_str(), 10000, 0., 10.);    
       }     
   }

   MustOEtrue_vs_scoreThreshold_seedEta.resize(etaCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++)
   {    
       MustOEtrue_vs_scoreThreshold_seedEta[iBin].resize(scoreCuts.size());
       for(unsigned int jBin=0; jBin<scoreCuts.size(); jBin++)
       { 
           MustOEtrue_vs_scoreThreshold_seedEta[iBin][jBin] = new TH1F(std::string("MustOEtrue_vs_scoreThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)).c_str(), std::string("MustOEtrue_vs_scoreThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)).c_str(), 10000, 0., 10.);
       }     
   }

   MustOEtrue_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       MustOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           MustOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               MustOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TH1F(std::string("MustOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("MustOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 10000, 0., 10.);
           }
       }
   }

   ErecoOEtrue_vs_scoreThreshold_seedEt_EB.resize(etCuts.size()-1); 
   ErecoOEtrue_vs_scoreThreshold_seedEt_EE.resize(etCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       ErecoOEtrue_vs_scoreThreshold_seedEt_EB[iBin].resize(scoreCuts.size());
       ErecoOEtrue_vs_scoreThreshold_seedEt_EE[iBin].resize(scoreCuts.size());
       for(unsigned int jBin=0; jBin<scoreCuts.size(); jBin++)
       {
           boost::replace_all(scoreCuts.at(jBin), ".", "_");
           ErecoOEtrue_vs_scoreThreshold_seedEt_EB[iBin][jBin] = new TH1F(std::string("ErecoOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EB").c_str(), std::string("ErecoOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EB").c_str(), 10000, 0., 10.);
           ErecoOEtrue_vs_scoreThreshold_seedEt_EE[iBin][jBin] = new TH1F(std::string("ErecoOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EE").c_str(), std::string("ErecoOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EE").c_str(), 10000, 0., 10.);    
       }     
   }

   ErecoOEtrue_vs_scoreThreshold_seedEta.resize(etaCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++)
   {    
       ErecoOEtrue_vs_scoreThreshold_seedEta[iBin].resize(scoreCuts.size());
       for(unsigned int jBin=0; jBin<scoreCuts.size(); jBin++)
       { 
           ErecoOEtrue_vs_scoreThreshold_seedEta[iBin][jBin] = new TH1F(std::string("ErecoOEtrue_vs_scoreThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)).c_str(), std::string("ErecoOEtrue_vs_scoreThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)).c_str(), 10000, 0., 10.);
       }     
   }

   ErecoOEtrue_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       ErecoOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           ErecoOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               ErecoOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TH1F(std::string("ErecoOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("ErecoOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 10000, 0., 10.);
           }
       }
   }

   EtrueOEtrue_vs_scoreThreshold_seedEt_EB.resize(etCuts.size()-1); 
   EtrueOEtrue_vs_scoreThreshold_seedEt_EE.resize(etCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       EtrueOEtrue_vs_scoreThreshold_seedEt_EB[iBin].resize(scoreCuts.size());
       EtrueOEtrue_vs_scoreThreshold_seedEt_EE[iBin].resize(scoreCuts.size());
       for(unsigned int jBin=0; jBin<scoreCuts.size(); jBin++)
       {
           boost::replace_all(scoreCuts.at(jBin), ".", "_");
           EtrueOEtrue_vs_scoreThreshold_seedEt_EB[iBin][jBin] = new TH1F(std::string("EtrueOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EB").c_str(), std::string("EtrueOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EB").c_str(), 10000, 0., 10.);
           EtrueOEtrue_vs_scoreThreshold_seedEt_EE[iBin][jBin] = new TH1F(std::string("EtrueOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EE").c_str(), std::string("EtrueOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)+"_EE").c_str(), 10000, 0., 10.);    
       }     
   }

   EtrueOEtrue_vs_scoreThreshold_seedEta.resize(etaCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++)
   {    
       EtrueOEtrue_vs_scoreThreshold_seedEta[iBin].resize(scoreCuts.size());
       for(unsigned int jBin=0; jBin<scoreCuts.size(); jBin++)
       { 
           EtrueOEtrue_vs_scoreThreshold_seedEta[iBin][jBin] = new TH1F(std::string("EtrueOEtrue_vs_scoreThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)).c_str(), std::string("EtrueOEtrue_vs_scoreThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_simScore_"+scoreCuts.at(jBin)).c_str(), 10000, 0., 10.);
       }     
   }

   EtrueOEtrue_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       EtrueOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           EtrueOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               EtrueOEtrue_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TH1F(std::string("EtrueOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("EtrueOEtrue_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 10000, 0., 10.);
           }
       }
   }

   dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TH2F(std::string("dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
           }
       }
   } 

   dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin][jBin][kBin] = new TH2F(std::string("dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inMustache").c_str(), std::string("dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
           }
       }
   }

   dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin][jBin][kBin] = new TH2F(std::string("dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notMustache").c_str(), std::string("dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
           }
       }
   }

   dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("dEtadPhi_Mustache_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
           }
       }
   }

   dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_inCalo.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_inCalo[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_inCalo[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_inCalo[iBin][jBin][kBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inCalo").c_str(), std::string("dEtadPhi_Mustache_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inCalo").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
           }
       }
   }

   dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_notCalo.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_notCalo[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_notCalo[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_notCalo[iBin][jBin][kBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notCalo").c_str(), std::string("dEtadPhi_Mustache_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notCalo").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
           }
       }
   } 

   TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TProfile2D(std::string("TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
           }
       }
   }

   TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin][jBin][kBin] = new TProfile2D(std::string("TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inMustache").c_str(), std::string("TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
           }
       }
   }

   TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin][jBin][kBin] = new TProfile2D(std::string("TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notMustache").c_str(), std::string("TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.);
           }
       }
   }

   RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TProfile2D(std::string("RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), std::string("RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
           }
       }
   }

   RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[iBin][jBin][kBin] = new TProfile2D(std::string("RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inMustache").c_str(), std::string("RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_inMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
           }
       }
   }

   RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache.resize(scoreCuts.size());
   for(unsigned int iBin=0; iBin<scoreCuts.size(); iBin++)
   {
       RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[iBin][jBin][kBin] = new TProfile2D(std::string("RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notMustache").c_str(), std::string("RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_simScore_"+scoreCuts.at(iBin)+"_notMustache").c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 10000.);
           }
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

double my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  double xcur = x[0];
  double alphaL = par[0];
  double nL = par[1];
  double mu = par[2];
  double sigma = par[3];
  double N = par[4];
  double alphaR = par[5];
  double nR = par[6];
  double t = (xcur-mu)/sigma;
  double absAlphaL = fabs((double)alphaL);
  double invAbsAlphaL = 1./absAlphaL;
  double absAlphaR = fabs((double)alphaR);
  double invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    double AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    double BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    double AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    double BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

double mycruijff(double* x, double* par) 
{
  double m = x[0];
  double m0 = par[0];
  double sigmaL = par[1];
  double sigmaR = par[2];
  double alphaL = par[3];
  double alphaR = par[4];
  double N = par[5];
  double dx =  (m-m0) ;
  double sigma = dx<0 ? sigmaL: sigmaR ;
  double alpha = dx<0 ? alphaL: alphaR ;
  double f = 2*sigma*sigma + alpha*dx*dx ; 
  
  return N*exp(-dx*dx/f) ;
}

TF1* makeCruijffFit(TH1* hist, float xmin, float xmax, float meanSet=1., float sigmaLSet=0.1, float sigmaRSet=0.1, float alphaLSet=0.1, float alphaRSet=0.1)
{
  //hist->Scale(1./hist->GetEntries()); 
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;
  float nrEntries = hist->Integral();
  RooRealVar  nsig("N_{S}", "#signal events", nrEntries, nrEntries*0.1, nrEntries*10.);
  RooRealVar mean( "#DeltaE", "mean_{cb}", meanSet ,meanSet-3.*hist->GetMean(),meanSet+3.*hist->GetMean(),"");
  RooRealVar sigmaL("#sigma_{L}","#sigma_{L}", sigmaLSet, 0., 0.5);
  RooRealVar sigmaR("#sigma_{R}","#sigma_{R}", sigmaRSet, 0., 0.5);
  RooRealVar alphaL( "alpha_{L}", "alpha_{L}", alphaLSet, 0., 20.);
  RooRealVar alphaR( "alpha_{R}", "alpha_{R}", alphaRSet, 0., 20.);
  CruijffPdf cruijff("cruijff","cruijff",res,mean,sigmaL,sigmaR,alphaL,alphaR);

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
   
  RooAddPdf model("model", "model", RooArgList(cruijff), RooArgList(nsig));
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1),RooFit::PrintEvalErrors(-1),RooFit::SumW2Error(true));
  auto fitRes = model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1));
  if( (sigmaL.getValV()<=0.001 || sigmaR.getValV()<=0.001) || fitRes->edm()>10){
    std::cout <<" fit status "<<fitRes->status()<<" sigma "<<sigmaL.getValV()<<" "<<sigmaR.getValV()<<" nrsig "<<nsig.getValV()<<" edm "<<fitRes->edm()<<std::endl;
    double maxSigma = std::max(sigmaL.getValV(),sigmaR.getValV());
    sigmaL.setVal(maxSigma);
    sigmaR.setVal(maxSigma);
    fitRes = model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1));
    std::cout <<"trying again "<<fitRes->status()<<" sigma "<<sigmaL.getValV()<<" "<<sigmaR.getValV()<<" nrsig "<<nsig.getValV()<<" edm "<<fitRes->edm()<<std::endl;
  }
  
  double errors[6] ={(mean.getAsymErrorHi()+mean.getAsymErrorLo())/2., (sigmaL.getAsymErrorHi()+sigmaL.getAsymErrorLo())/2., (sigmaR.getAsymErrorHi()+sigmaR.getAsymErrorLo())/2., (alphaL.getAsymErrorHi()+alphaL.getAsymErrorLo())/2., (alphaR.getAsymErrorHi()+alphaR.getAsymErrorLo())/2., 0.};
  TF1* Cruijff = new TF1("Cruijff",&mycruijff,xmin,xmax,6);
  Cruijff->SetParameters(mean.getVal(), sigmaL.getVal(), sigmaR.getVal(), alphaL.getVal(), alphaR.getVal(), hist->GetMaximum());
  Cruijff->SetParErrors(errors);

  return Cruijff;
}

TF1* makeDoubleCBFit(TH1* hist,float xmin,float xmax, float meanSet=1., float sigmaSet=0.1, float alpha1Set=0.1, float n1Set=2., float alpha2Set=0.1, float n2Set=2.)
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  float nrEntries = hist->Integral();
  RooRealVar nsig("N_{S}", "#signal events", nrEntries,nrEntries*0.5,nrEntries*2.);
  RooRealVar mean( "#DeltaE", "mean_{cb}", meanSet ,meanSet-3.*hist->GetMean(),meanSet+3.*hist->GetMean(),"");
  RooRealVar cbSigma("#sigma_{CB}","CB Width", sigmaSet, 0., 0.5);
  RooRealVar alpha1( "alpha_{1}", "alpha_{1}", alpha1Set, 0., 20.);
  RooRealVar alpha2( "alpha_{2}", "alpha_{2}", alpha2Set,0.,20.);
  RooRealVar n1( "n_{1}", "n_{1}", n1Set,0.,5000.);
  RooRealVar n2( "n_{2}", "n_{2}", n2Set,0.,5000.);

  DoubleCBPdf doubleCB("doubleCB","doubleCB",res,mean,cbSigma,alpha1,n1,alpha2,n2);
  RooAddPdf model("model", "model", RooArgList(doubleCB), RooArgList(nsig));

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
 
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1)); 
  
  double errors[7] ={(mean.getAsymErrorHi()+mean.getAsymErrorLo())/2., (cbSigma.getAsymErrorHi()+cbSigma.getAsymErrorLo())/2., (alpha1.getAsymErrorHi()+alpha1.getAsymErrorLo())/2., (n1.getAsymErrorHi()+n1.getAsymErrorLo())/2., (alpha2.getAsymErrorHi()+alpha2.getAsymErrorLo())/2., (n2.getAsymErrorHi()+n2.getAsymErrorLo())/2., 0.};
   
  hist->GetXaxis()->SetRangeUser(0.2,1.8); 
  double maximum = hist->GetMaximum();
  hist->GetXaxis()->SetRangeUser(0.,2.); 
  TF1* DoubleCB = new TF1("DoubleCB",&my2sideCrystalBall,xmin,xmax,7);
  DoubleCB->SetParameters(mean.getVal(), cbSigma.getVal(), alpha1.getVal(), n1.getVal(), alpha2.getVal(), n2.getVal(), maximum);
  DoubleCB->SetParErrors(errors);

  return DoubleCB;
}

TF1* fitHisto(TH1* hist, std::string fitFunction_="cruijff")
{
   hist->GetXaxis()->SetRangeUser(0.2,1.8); 
   int binmax = hist->GetMaximumBin(); 
   hist->GetXaxis()->SetRangeUser(0.,2.); 
   double xMAX = hist->GetXaxis()->GetBinCenter(binmax);  

   TF1* fit1;
   TF1* fit2;
   TF1* fit3; 
   //std::cout << "fitHisto: " << hist->GetName() << std::endl;
   if(fitFunction_=="doubleCB"){ 
      fit1 = makeDoubleCBFit(hist,xMAX-1.5*hist->GetRMS(),xMAX+1.5*hist->GetRMS());
      fit2 = makeDoubleCBFit(hist, fit1->GetParameter(0)-5.*fit1->GetParameter(1), fit1->GetParameter(0)+5.*fit1->GetParameter(1), fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2), fit1->GetParameter(3), fit1->GetParameter(4), fit1->GetParameter(5));
      fit3 = makeDoubleCBFit(hist, fit2->GetParameter(0)-5.*fit2->GetParameter(1), fit2->GetParameter(0)+5.*fit2->GetParameter(1), fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(2), fit2->GetParameter(3), fit2->GetParameter(4), fit2->GetParameter(5));
      return fit3;
   }else{
      float sigma = 0;
      fit1 = makeCruijffFit(hist,xMAX-1.5*hist->GetRMS(),xMAX+1.5*hist->GetRMS(), 1., hist->GetRMS()/2., hist->GetRMS()/2., 0.1, 0.1);
      sigma = (fit1->GetParameter(1)+fit1->GetParameter(2))/2.;
      fit2 = makeCruijffFit(hist, fit1->GetParameter(0)-5.*sigma, fit1->GetParameter(0)+5.*sigma, fit1->GetParameter(0), fit1->GetParameter(1), fit1->GetParameter(2), fit1->GetParameter(3), fit1->GetParameter(4)); 
      sigma = (fit2->GetParameter(1)+fit2->GetParameter(2))/2.;
      fit3 = makeCruijffFit(hist, fit2->GetParameter(0)-3.*sigma, fit2->GetParameter(0)+3.*sigma, fit2->GetParameter(0), fit2->GetParameter(1), fit2->GetParameter(2), fit2->GetParameter(3), fit2->GetParameter(4));  
      return fit3;
   }
}

void drawHistFunc(TH1F* hist, TF1* func, std::string x_label, std::string Name)
{
   gStyle->SetOptStat(0000); 
   hist->SetMaximum(hist->GetMaximum()*1.05);
   hist->SetLineColor(kRed+1);
   hist->SetMarkerColor(kRed+1);
   hist->SetLineWidth(2);
   hist->GetXaxis()->SetTitle(x_label.c_str());
   hist->GetXaxis()->SetRangeUser(hist->GetMean()-1.,hist->GetMean()+1.);   

   func->SetLineColor(kRed+1);

   TLegend* legend = new TLegend(0.57, 0.77, 0.72, 0.89);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(hist,std::string("#mu = "+to_string(func->GetParameter(0))+" +/- "+to_string(func->GetParError(0))).c_str(),"");
   legend -> AddEntry(hist,std::string("#sigma = "+to_string((func->GetParameter(2)+func->GetParameter(1))/2.)+" +/- "+to_string(0.5*sqrt(func->GetParError(2)*func->GetParError(2) +func->GetParError(1)*func->GetParError(1)))).c_str(),"");
     
   TCanvas* c = new TCanvas();
   hist->Draw("E");    
   func->Draw("L,same");
   legend -> Draw("same");  
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf"); 
   gStyle->SetOptStat(1110); 
   
}

void drawHist(TH1F* hist, std::string x_label, std::string Name)
{
   gStyle->SetOptStat(0000); 
   hist->SetMaximum(hist->GetMaximum()*1.05);
   hist->SetLineColor(kRed+1);
   hist->SetMarkerColor(kRed+1);
   hist->SetLineWidth(2);
   hist->GetXaxis()->SetTitle(x_label.c_str());
   hist->GetXaxis()->SetRangeUser(hist->GetMean()-1.,hist->GetMean()+1.);   

   TCanvas* c = new TCanvas();
   hist->Draw("E");    
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf"); 
   gStyle->SetOptStat(1110); 
   
}

std::pair<float,float> computeRange(TGraphErrors* graph)
{
   std::vector<double> y_points;
   double x,y;
   for(int i = 0; i<graph->GetN(); i++){ 
       graph->GetPoint(i,x,y);   
       if(y>=0.) y_points.push_back(y); 
   } 
   std::sort(y_points.begin(),y_points.end());
   if(y_points.size() == 0) return make_pair(0.,2.); 
   else return make_pair(y_points.at(0),y_points.at(y_points.size()-1)); 
}

double computeMean(TH1 * hist, int imin, int imax)
{
   if(imin<1) imin = 1;
   if(imax>hist->GetNbinsX()) imax = hist->GetNbinsX();
   
   double val = 0.;
   double total = 0.;
   for(int ibin=imin; ibin<imax+1; ibin++){
       val+=hist->GetXaxis()->GetBinCenter(ibin) *hist->GetBinContent(ibin);
       total+=hist->GetBinContent(ibin); 
   } 

   if(total==0) return -1.;
   else return val/total;   
}

std::pair<double,double> computeEffectiveSigma(TH1 * hist)
{
    TAxis *xaxis = hist->GetXaxis();
    int nb = xaxis->GetNbins();
    if(nb < 10) {
       cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
       return std::make_pair(-1.,-1.);
    }

    double bwid = xaxis->GetBinWidth(1);
    if(bwid == 0) {
       cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
       return std::make_pair(-1.,-1.);
    }

    //double xmax = xaxis->GetXmax();
    double xmin = xaxis->GetXmin();
    double ave = hist->GetMean();
    double rms = hist->GetRMS();
    //int aveBin = hist->FindBin(hist->GetMean());

    double total=0.;
    for(int i=0; i<nb+2; i++) {
        total+=hist->GetBinContent(i);
    }
    
    int ierr=0;
    int ismin=999;

    double rlim=0.68269*total;
    int nrms=rms/(bwid);    // Set scan size to +/- rms
    if(nrms > nb/10) nrms=nb/10; // Could be tuned...

    double widmin=9999999.;
    int jbin = nb;
    int kbin = 1;  
     

    for(int iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre

        int ibm=(ave-xmin)/bwid+1+iscan;
        double x=(ibm-0.5)*bwid+xmin;
        double xj=x;
        double xk=x;
        int jbm=ibm;
        int kbm=ibm;
        double bin=hist->GetBinContent(ibm);
        total=bin;

        for(int j=1;j<nb;j++){
            if(jbm < nb) {
               jbm++;
               xj+=bwid;
               bin=hist->GetBinContent(jbm);
               total+=bin;
               jbin = jbm; 
               if(total > rlim) break;
            }else ierr=1;

            if(kbm > 0) {
               kbm--;
               xk-=bwid;
               bin=hist->GetBinContent(kbm);
               total+=bin;
               kbin = kbm;  
               if(total > rlim) break;
            }else ierr=1;
        }

        double dxf=(total-rlim)*bwid/bin;
        double wid=(xj-xk+bwid-dxf)*0.5;

        if(wid < widmin){
           widmin=wid;
           ismin=iscan;
        }
    }

    if(ismin == nrms || ismin == -nrms) ierr=3;
    if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;

    //std::cout << "EFFECTIVE SIGMA: " << kbin << " - " << aveBin << " - " << jbin << " - " << nb << " - " << (float)hist->Integral(kbin,jbin)/(float)hist->Integral() << " - " << (float)hist->Integral(kbin,aveBin)/(float)hist->Integral() << " - " << (float)hist->Integral(aveBin,jbin)/(float)hist->Integral() << " - " << widmin << std::endl;
    return std::make_pair(computeMean(hist,kbin,jbin),widmin);
}

std::vector<double> evaluateHisto(TH1F* hist, std::string fitFunction_="cruijff", bool doEffective = false)
{
   double mean = 0.;
   double mean_error = 0.;
   double sigma = 0.;
   double sigma_error = 0.; 

   TF1* doubleCB;
   if(hist->Integral()>20){
      if(!doEffective){
         doubleCB = fitHisto(hist, fitFunction_);
         mean = doubleCB->GetParameter(0);
         mean_error = doubleCB->GetParError(0);   
         sigma = (doubleCB->GetParameter(1)+doubleCB->GetParameter(2))/2.;
         sigma_error = sqrt(doubleCB->GetParError(1)*doubleCB->GetParError(1)+doubleCB->GetParError(2)*doubleCB->GetParError(2))/2.;
         if(mean_error>0.1) mean = -1.; 
         if(sigma_error>0.03) sigma = -1.; 
      }else{
         std::pair<double,double> effective = computeEffectiveSigma(hist);
         mean = effective.first;
         mean_error = 0.;
         sigma =  effective.second;
         sigma_error = 0.; 
         if(effective.first>10.) mean = -1.; 
         if(effective.second>10.) sigma = -1.; 
      }
   }else{
     std::cout << "WARNING: too few events in " << hist->GetName() << std::endl;
     mean = -1.;
     sigma = -1.;
  } 
  
  return std::vector<double>({mean,mean_error,sigma,sigma_error});       
}

std::pair<TGraphErrors*,TGraphErrors*> makeFitProfile(std::vector<TH1F*>* vecHist, std::vector<std::string> scoreCuts, std::string xTitle, std::string fitFunction_="cruijff", bool doEffective = false)
{
   TF1* doubleCB;
   int nBins = vecHist->size();
   double x[nBins], x_error[nBins], y_mean[nBins], y_meanError[nBins], y_sigma[nBins], y_sigmaError[nBins];
   
   for(unsigned int iBin=0; iBin<vecHist->size(); iBin++)
   {
       x[iBin] = std::stod(scoreCuts.at(iBin)); 
       x_error[iBin] = 0.;
       //std::cout << "makeFit: " << iBin << " - " << vecHist->at(iBin)->Integral() << std::endl;
       if(vecHist->at(iBin)->Integral()>20){
          if(!doEffective){
             doubleCB = fitHisto(vecHist->at(iBin), fitFunction_);
             y_mean[iBin] = doubleCB->GetParameter(0);
             y_meanError[iBin] = doubleCB->GetParError(0);   
             if(y_meanError[iBin]>0.1){
                y_mean[iBin] = -1.;
                y_meanError[iBin] = 0.;     
             }  
             y_sigma[iBin] = (doubleCB->GetParameter(1)+doubleCB->GetParameter(2))/2.;
             y_sigmaError[iBin] = sqrt(doubleCB->GetParError(1)*doubleCB->GetParError(1)+doubleCB->GetParError(2)*doubleCB->GetParError(2))/2.;
             if(y_sigmaError[iBin]>0.03){
                y_sigma[iBin] = -1.;
                y_sigmaError[iBin] = 0.;     
             }    
             if(y_meanError[iBin]<0.1 && y_sigmaError[iBin]<0.03) drawHistFunc(vecHist->at(iBin),doubleCB, std::string(""), std::string(vecHist->at(iBin)->GetName())); 
          }else{
            //drawHist(vecHist->at(iBin), std::string("EoEtrue"), std::string(string("h1_")+string(vecHist->at(iBin)->GetName())));
            std::pair<double,double> effective = computeEffectiveSigma(vecHist->at(iBin));
            y_mean[iBin] = effective.first;
            y_meanError[iBin] = 0.;
            y_sigma[iBin] =  effective.second;
            y_sigmaError[iBin] = 0.; 
            if(effective.first>10.) y_mean[iBin] = -1.; 
            if(effective.second>10.) y_sigma[iBin] = -1.; 
          }   
       }else{
          y_mean[iBin] = -1.;
          y_meanError[iBin] = 0.;
          y_sigma[iBin] = -1.;
          y_sigmaError[iBin] = 0.;
       } 
   }
         
   TGraphErrors* gr_Mean = new TGraphErrors(nBins,x,y_mean,x_error,y_meanError);
   TGraphErrors* gr_Sigma = new TGraphErrors(nBins,x,y_sigma,x_error,y_sigmaError);;

   //delete doubleCB;
   return std::make_pair(gr_Mean,gr_Sigma);
   
}

std::pair<TH2F*,TH2F*> makeFitHist2D(int bin, std::vector<std::vector<TH1F*>> vecHist, std::vector<float> etCuts, std::vector<float> etaCuts, std::string fitFunction_="cruijff", bool doEffective = false)
{
   TF1* doubleCB;
   TH2F* h2_mean = new TH2F(std::string("h2_mean_"+std::to_string(bin)).c_str(), std::string("h2_mean_"+std::to_string(bin)).c_str(), etCuts.size()-1, &etCuts[0], etaCuts.size()-1, &etaCuts[0]);  
   TH2F* h2_res = new TH2F(std::string("h2_res_"+std::to_string(bin)).c_str(), std::string("h2_res_"+std::to_string(bin)).c_str(), etCuts.size()-1, &etCuts[0], etaCuts.size()-1, &etaCuts[0]); 
   
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++){
           if(vecHist[iBin][jBin]->Integral()<=20) continue;
           if(!doEffective){
              doubleCB = fitHisto(vecHist[iBin][jBin], fitFunction_);
              float mean_error = doubleCB->GetParError(0);
              float sigma_error = sqrt(doubleCB->GetParError(1)*doubleCB->GetParError(1)+doubleCB->GetParError(2)*doubleCB->GetParError(2))/2.;
              if(mean_error<=0.1 && sigma_error<=0.03){
                 h2_mean->SetBinContent(iBin+1,jBin+1,doubleCB->GetParameter(0)); 
                 h2_res->SetBinContent(iBin+1,jBin+1,(doubleCB->GetParameter(1)+doubleCB->GetParameter(2))/2.); 
              } 
           }else{
              drawHist(vecHist[iBin][jBin], std::string("EoEtrue"), std::string(string("h1_")+string(vecHist[iBin][jBin]->GetName())));
              std::pair<double,double> effective = computeEffectiveSigma(vecHist[iBin][jBin]);
              if(effective.first<=10. && effective.second<=10.){
                 h2_mean->SetBinContent(iBin+1,jBin+1,effective.first); 
                 h2_res->SetBinContent(iBin+1,jBin+1,effective.second); 
              } 
           }
       }
   } 
 
   return std::make_pair(h2_mean,h2_res); 
}

TGraphErrors* scaledGraph(TGraphErrors* gr_SuperCluster_ErecoOEtrue, TGraphErrors* gr_SuperCluster_MustOEtrue)
{
   TGraphErrors* gr = new TGraphErrors();
   double x_true,y_true;
   double x_must,y_must;
   for(int i=0; i<gr_SuperCluster_ErecoOEtrue->GetN(); i++)
   {
       gr_SuperCluster_ErecoOEtrue->GetPoint(i,x_true,y_true);
       gr_SuperCluster_MustOEtrue->GetPoint(i,x_must,y_must);
       gr->SetPoint(i,x_true,y_true/y_must);
       gr->SetPointError(i,0.,0.);
   }  
   return gr;
}

void drawGraph(TGraphErrors* gr_SuperCluster_ErecoOEtrue, TGraphErrors* gr_SuperCluster_MustOEtrue, TGraphErrors* gr_SuperCluster_EtrueOEtrue, std::string xtitle, std::string ytitle, std::string Name, bool log = false, std::string outDir = "")
{ 
   gStyle->SetOptStat(0000);  

   gr_SuperCluster_ErecoOEtrue->SetTitle(Name.c_str());
   gr_SuperCluster_ErecoOEtrue->SetLineColor(kBlue+1);
   gr_SuperCluster_ErecoOEtrue->SetMarkerStyle(20);
   gr_SuperCluster_ErecoOEtrue->SetMarkerSize(0.5);
   gr_SuperCluster_ErecoOEtrue->SetMarkerColor(kBlue+1);
   gr_SuperCluster_ErecoOEtrue->SetLineWidth(2);
   gr_SuperCluster_ErecoOEtrue->GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_SuperCluster_ErecoOEtrue->GetYaxis()->SetTitle(ytitle.c_str());

   gr_SuperCluster_MustOEtrue->SetLineColor(kRed+1);
   gr_SuperCluster_MustOEtrue->SetMarkerStyle(20);
   gr_SuperCluster_MustOEtrue->SetMarkerSize(0.5);
   gr_SuperCluster_MustOEtrue->SetMarkerColor(kRed+1);
   gr_SuperCluster_MustOEtrue->SetLineWidth(2);   
  
   gr_SuperCluster_EtrueOEtrue->SetLineColor(kGreen+1);
   gr_SuperCluster_EtrueOEtrue->SetMarkerStyle(20);
   gr_SuperCluster_EtrueOEtrue->SetMarkerSize(0.5);
   gr_SuperCluster_EtrueOEtrue->SetMarkerColor(kGreen+1);
   gr_SuperCluster_EtrueOEtrue->SetLineWidth(2); 
  
   TGraphErrors* gr_SuperCluster_ratio = scaledGraph(gr_SuperCluster_ErecoOEtrue, gr_SuperCluster_MustOEtrue);
   gr_SuperCluster_ratio->SetLineColor(kGreen+1);
   gr_SuperCluster_ratio->SetMarkerStyle(20);
   gr_SuperCluster_ratio->SetMarkerSize(0.5);
   gr_SuperCluster_ratio->SetMarkerColor(kGreen+1);
   gr_SuperCluster_ratio->SetLineWidth(2); 
   gr_SuperCluster_ratio->GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_SuperCluster_ratio->GetYaxis()->SetTitle("Reco/Mustache");

   std::vector<float> ranges;
   ranges.push_back(0.99*computeRange(gr_SuperCluster_ErecoOEtrue).first); 
   ranges.push_back(1.01*computeRange(gr_SuperCluster_ErecoOEtrue).second);  
   ranges.push_back(0.99*computeRange(gr_SuperCluster_MustOEtrue).first); 
   ranges.push_back(1.01*computeRange(gr_SuperCluster_MustOEtrue).second);  
   //ranges.push_back(0.99*computeRange(gr_SuperCluster_EtrueOEtrue).first); 
   //ranges.push_back(1.01*computeRange(gr_SuperCluster_EtrueOEtrue).second); 
   std::sort(ranges.begin(), ranges.end());

   float max = ranges.at(ranges.size()-1);
   float min = ranges.at(0);
   gr_SuperCluster_ErecoOEtrue->GetYaxis()->SetRangeUser(min,max);    

   TLegend* legend = new TLegend(0.799, 0.73, 0.999, 0.95);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.03);
   legend -> AddEntry(gr_SuperCluster_ErecoOEtrue,"Reco/True","L");
   legend -> AddEntry(gr_SuperCluster_MustOEtrue,"Mustache/True","L");
   //legend -> AddEntry(gr_SuperCluster_EtrueOEtrue,"True/True","L");

   TCanvas* c = new TCanvas();
   if(log) c->SetLogx();
   gr_SuperCluster_ErecoOEtrue->Draw("APL");
   gr_SuperCluster_MustOEtrue->Draw("PL, same");
   //gr_SuperCluster_EtrueOEtrue->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(outDir+Name+".png").c_str(),"png");
   c->SaveAs(std::string(outDir+Name+".pdf").c_str(),"pdf");

   c->cd();
   if(log) c->SetLogx();
   gr_SuperCluster_ratio->Draw("APL");
   c->SaveAs(std::string(outDir+Name+"_ratio.png").c_str(),"png");
   c->SaveAs(std::string(outDir+Name+"_ratio.pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
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
void drawPlots(std::string fitFunction_, std::vector<std::string> etCuts, std::vector<std::string> etaCuts, std::vector<std::string> scoreCuts, std::string matching_, std::string outDir_)
{
   
   //68% Quantile
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_MustOEtrue_vs_simScore_seedEt_eff_EB;
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_MustOEtrue_vs_simScore_seedEt_eff_EE; 
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_MustOEtrue_vs_simScore_seedEta_eff; 
   gr_MustOEtrue_vs_simScore_seedEt_eff_EB.resize(etCuts.size()-1);
   gr_MustOEtrue_vs_simScore_seedEt_eff_EE.resize(etCuts.size()-1); 
   gr_MustOEtrue_vs_simScore_seedEta_eff.resize(etaCuts.size()-1); 
 
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_ErecoOEtrue_vs_simScore_seedEt_eff_EB;
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_ErecoOEtrue_vs_simScore_seedEt_eff_EE; 
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_ErecoOEtrue_vs_simScore_seedEta_eff; 
   gr_ErecoOEtrue_vs_simScore_seedEt_eff_EB.resize(etCuts.size()-1);
   gr_ErecoOEtrue_vs_simScore_seedEt_eff_EE.resize(etCuts.size()-1); 
   gr_ErecoOEtrue_vs_simScore_seedEta_eff.resize(etaCuts.size()-1); 

   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_EtrueOEtrue_vs_simScore_seedEt_eff_EB;
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_EtrueOEtrue_vs_simScore_seedEt_eff_EE; 
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_EtrueOEtrue_vs_simScore_seedEta_eff; 
   gr_EtrueOEtrue_vs_simScore_seedEt_eff_EB.resize(etCuts.size()-1);
   gr_EtrueOEtrue_vs_simScore_seedEt_eff_EE.resize(etCuts.size()-1); 
   gr_EtrueOEtrue_vs_simScore_seedEta_eff.resize(etaCuts.size()-1);  

   system(string("mkdir -p "+outDir_).c_str()); 
   system(string("cp index.php "+outDir_).c_str()); 

   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
       gr_ErecoOEtrue_vs_simScore_seedEt_eff_EB[iBin] = makeFitProfile(&ErecoOEtrue_vs_scoreThreshold_seedEt_EB[iBin],scoreCuts,matching_,fitFunction_,true); 
       gr_EtrueOEtrue_vs_simScore_seedEt_eff_EB[iBin] = makeFitProfile(&EtrueOEtrue_vs_scoreThreshold_seedEt_EB[iBin],scoreCuts,matching_,fitFunction_,true);    
       gr_MustOEtrue_vs_simScore_seedEt_eff_EB[iBin] = makeFitProfile(&MustOEtrue_vs_scoreThreshold_seedEt_EB[iBin],scoreCuts,matching_,fitFunction_,true); 
       drawGraph(gr_ErecoOEtrue_vs_simScore_seedEt_eff_EB[iBin].first, gr_MustOEtrue_vs_simScore_seedEt_eff_EB[iBin].first, gr_EtrueOEtrue_vs_simScore_seedEt_eff_EB[iBin].first, matching_, std::string("#mu"), std::string("EoEtrue_vs_simScore_Mean_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EB"), true, outDir_); 
       drawGraph(gr_ErecoOEtrue_vs_simScore_seedEt_eff_EB[iBin].second, gr_MustOEtrue_vs_simScore_seedEt_eff_EB[iBin].second, gr_EtrueOEtrue_vs_simScore_seedEt_eff_EB[iBin].second, matching_, std::string("#sigma"), std::string("EoEtrue_vs_simScore_Resolution_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EB"), true, outDir_);

       gr_ErecoOEtrue_vs_simScore_seedEt_eff_EE[iBin] = makeFitProfile(&ErecoOEtrue_vs_scoreThreshold_seedEt_EE[iBin],scoreCuts,matching_,fitFunction_,true); 
       gr_EtrueOEtrue_vs_simScore_seedEt_eff_EE[iBin] = makeFitProfile(&EtrueOEtrue_vs_scoreThreshold_seedEt_EE[iBin],scoreCuts,matching_,fitFunction_,true); 
       gr_MustOEtrue_vs_simScore_seedEt_eff_EE[iBin] = makeFitProfile(&MustOEtrue_vs_scoreThreshold_seedEt_EE[iBin],scoreCuts,matching_,fitFunction_,true); 
       drawGraph(gr_ErecoOEtrue_vs_simScore_seedEt_eff_EE[iBin].first, gr_MustOEtrue_vs_simScore_seedEt_eff_EE[iBin].first, gr_EtrueOEtrue_vs_simScore_seedEt_eff_EE[iBin].first, matching_, std::string("#mu"), std::string("EoEtrue_vs_simScore_Mean_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EE"), true, outDir_); 
       drawGraph(gr_ErecoOEtrue_vs_simScore_seedEt_eff_EE[iBin].second, gr_MustOEtrue_vs_simScore_seedEt_eff_EE[iBin].second, gr_EtrueOEtrue_vs_simScore_seedEt_eff_EE[iBin].second, matching_, std::string("#sigma"), std::string("EoEtrue_vs_simScore_Resolution_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EE"), true, outDir_); 
   }

   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
       gr_ErecoOEtrue_vs_simScore_seedEta_eff[iBin] = makeFitProfile(&ErecoOEtrue_vs_scoreThreshold_seedEta[iBin],scoreCuts,matching_,fitFunction_,true); 
       gr_EtrueOEtrue_vs_simScore_seedEta_eff[iBin] = makeFitProfile(&EtrueOEtrue_vs_scoreThreshold_seedEta[iBin],scoreCuts,matching_,fitFunction_,true); 
       gr_MustOEtrue_vs_simScore_seedEta_eff[iBin] = makeFitProfile(&MustOEtrue_vs_scoreThreshold_seedEta[iBin],scoreCuts,matching_,fitFunction_,true); 
       drawGraph(gr_ErecoOEtrue_vs_simScore_seedEta_eff[iBin].first, gr_MustOEtrue_vs_simScore_seedEta_eff[iBin].first, gr_EtrueOEtrue_vs_simScore_seedEta_eff[iBin].first, matching_, std::string("#mu"), std::string("EoEtrue_vs_simScore_Mean_seedEt_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_Effective"), true, outDir_); 
       drawGraph(gr_ErecoOEtrue_vs_simScore_seedEta_eff[iBin].second, gr_MustOEtrue_vs_simScore_seedEta_eff[iBin].second, gr_EtrueOEtrue_vs_simScore_seedEta_eff[iBin].second, matching_, std::string("#sigma"), std::string("EoEtrue_vs_simScore_Resolution_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_Effective"), true, outDir_);
   }

   for(unsigned int scoreBin=0; scoreBin<scoreCuts.size(); scoreBin++){
       if(scoreCuts.at(scoreBin)!="0.002") continue;
       for(unsigned int etBin=0; etBin<etCuts.size()-1; etBin++){
           for(unsigned int etaBin=0; etaBin<etaCuts.size()-1; etaBin++){
               std::string output_dir = outDir_+"/seedEt_"+etCuts.at(etBin)+"_"+etCuts.at(etBin+1)+"_seedEta_"+etaCuts.at(etaBin)+"_"+etaCuts.at(etaBin+1)+"/";
               system(string("mkdir -p "+output_dir).c_str()); 
               system(string("cp index.php "+output_dir).c_str()); 
               drawH2(dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Calo"));
               drawH2(dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Calo_inMustache"));
               drawH2(dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Calo_notMustache"));
               drawH2(dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache")); 
               drawH2(dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_inCalo[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache_inCalo")); 
               drawH2(dEtadPhi_Mustache_vs_scoreThreshold_seedEta_seedEt_notCalo[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache_notCalo")); 
               drawProfile2(TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"TrueRatio_dEtadPhi_Calo"));
               drawProfile2(TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"TrueRatio_dEtadPhi_Calo_inMustache"));
               drawProfile2(TrueRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"TrueRatio_dEtadPhi_Calo_notMustache"));
               drawProfile2(RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoRatio_dEtadPhi_Calo"), -1., -1.);
               drawProfile2(RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_inMustache[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoRatio_dEtadPhi_Calo_inMustache"), -1., 1.); 
               drawProfile2(RecoRatio_dEtadPhi_Calo_vs_scoreThreshold_seedEta_seedEt_notMustache[scoreBin][etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"RecoRatio_dEtadPhi_Calo_notMustache"), -1., 1.);
           }
       }  
   } 
}

