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
std::vector<std::vector<TH2F*>> dEtadPhi_Mustache_vs_seedEt_seedEta;
std::vector<std::vector<TH2F*>> dEtadPhi_DeepSC_vs_seedEt_seedEta;

std::vector<std::vector<TProfile2D*>> Energy_dEtadPhi_Mustache_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> Energy_dEtadPhi_DeepSC_vs_seedEt_seedEta;

std::vector<std::vector<TProfile2D*>> ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_seedEta;
std::vector<std::vector<TProfile2D*>> ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_seedEta;

//DEFINE BRANCHES
vector<float> *pfCluster_rawEnergy;
vector<float> *pfCluster_energy;
vector<float> *pfCluster_eta;
vector<float> *pfCluster_phi;
vector<vector<int> > *superCluster_pfClustersIndex;
vector<int> *superCluster_seedIndex;
vector<float> *superCluster_rawEnergy;
vector<int> *superCluster_iz;
vector<vector<int> > *deepSuperCluster_pfClustersIndex;
vector<int> *deepSuperCluster_seedIndex;
vector<float> *deepSuperCluster_rawEnergy;
vector<int> *deepSuperCluster_iz;

TBranch *b_pfCluster_rawEnergy;
TBranch *b_pfCluster_energy;
TBranch *b_pfCluster_eta;
TBranch *b_pfCluster_phi;
TBranch *b_superCluster_pfClustersIndex;
TBranch *b_superCluster_seedIndex;
TBranch *b_superCluster_rawEnergy;
TBranch *b_superCluster_iz;
TBranch *b_deepSuperCluster_pfClustersIndex;
TBranch *b_deepSuperCluster_seedIndex;
TBranch *b_deepSuperCluster_rawEnergy;
TBranch *b_deepSuperCluster_iz;

//setTreeBranches
void setTreeBranches(TTree* tree)
{
   tree->SetBranchStatus("*",0);   
   tree->SetBranchStatus("pfCluster_rawEnergy",1);  
   tree->SetBranchStatus("pfCluster_energy",1);  
   tree->SetBranchStatus("pfCluster_eta",1); 
   tree->SetBranchStatus("pfCluster_phi",1); 
   tree->SetBranchStatus("superCluster_pfClustersIndex",1); 
   tree->SetBranchStatus("superCluster_seedIndex",1);
   tree->SetBranchStatus("superCluster_rawEnergy",1); 
   tree->SetBranchStatus("superCluster_iz",1); 
   tree->SetBranchStatus("deepSuperCluster_pfClustersIndex",1); 
   tree->SetBranchStatus("deepSuperCluster_seedIndex",1);
   tree->SetBranchStatus("deepSuperCluster_rawEnergy",1); 
   tree->SetBranchStatus("deepSuperCluster_iz",1);  
  
   tree->SetBranchAddress("pfCluster_rawEnergy", &pfCluster_rawEnergy, &b_pfCluster_rawEnergy);
   tree->SetBranchAddress("pfCluster_energy", &pfCluster_energy, &b_pfCluster_energy); 
   tree->SetBranchAddress("pfCluster_eta", &pfCluster_eta, &b_pfCluster_eta);
   tree->SetBranchAddress("pfCluster_phi", &pfCluster_phi, &b_pfCluster_phi); 
   tree->SetBranchAddress("superCluster_pfClustersIndex", &superCluster_pfClustersIndex, &b_superCluster_pfClustersIndex);
   tree->SetBranchAddress("superCluster_seedIndex", &superCluster_seedIndex, &b_superCluster_seedIndex);
   tree->SetBranchAddress("superCluster_rawEnergy", &superCluster_rawEnergy, &b_superCluster_rawEnergy);
   tree->SetBranchAddress("superCluster_iz", &superCluster_iz, &b_superCluster_iz); 
   tree->SetBranchAddress("deepSuperCluster_pfClustersIndex", &deepSuperCluster_pfClustersIndex, &b_deepSuperCluster_pfClustersIndex);
   tree->SetBranchAddress("deepSuperCluster_seedIndex", &deepSuperCluster_seedIndex, &b_deepSuperCluster_seedIndex);
   tree->SetBranchAddress("deepSuperCluster_rawEnergy", &deepSuperCluster_rawEnergy, &b_deepSuperCluster_rawEnergy);
   tree->SetBranchAddress("deepSuperCluster_iz", &deepSuperCluster_iz, &b_deepSuperCluster_iz); 
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

double deltaPhi(double seed_phi, double cluster_phi)
{
     double dphi = seed_phi - cluster_phi;
     if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
     if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
     return dphi;
} 

double deltaEta(double seed_eta, double cluster_eta)
{
     double deta = 0.;
     if(seed_eta > 0.) deta = cluster_eta - seed_eta;
     if(seed_eta < 0.) deta = seed_eta - cluster_eta;
     return deta;
}

//set histograms
void setHistograms(std::vector<std::string> etCuts, std::vector<std::string> etaCuts)
{
   dEtadPhi_Mustache_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_Mustache_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_Mustache_vs_seedEt_seedEta[iBin][jBin] = new TH2F(std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   }

   dEtadPhi_DeepSC_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       dEtadPhi_DeepSC_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           dEtadPhi_DeepSC_vs_seedEt_seedEta[iBin][jBin] = new TH2F(std::string("dEtadPhi_DeepSC_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("dEtadPhi_DeepSC_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2);
        
       }
   } 

   Energy_dEtadPhi_Mustache_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       Energy_dEtadPhi_Mustache_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           Energy_dEtadPhi_Mustache_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("Energy_dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("Energy_dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 99999999.);
        
       }
   } 

   Energy_dEtadPhi_DeepSC_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       Energy_dEtadPhi_DeepSC_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           Energy_dEtadPhi_DeepSC_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("Energy_dEtadPhi_DeepSC_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("Energy_dEtadPhi_DeepSC_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 99999999.);
        
       }
   }

   ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.1);
        
       }
   } 

   ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_seedEta.resize(etCuts.size()-1); 
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_seedEta[iBin].resize(etaCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_seedEta[iBin][jBin] = new TProfile2D(std::string("ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), std::string("ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(jBin+1)).c_str(), 100,-0.6,0.6, 100, -0.2, 0.2, 0., 1.1);
        
       }
   } 
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
           drawH2(dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_Mustache")); 
           drawH2(dEtadPhi_DeepSC_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"dEtadPhi_DeepSC"));
           drawProfile2(Energy_dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"Energy_dEtadPhi_Mustache"), -1., -1.);   
           drawProfile2(Energy_dEtadPhi_DeepSC_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"Energy_dEtadPhi_DeepSC"), -1., -1.);   
           drawProfile2(ClusterToSeed_dEtadPhi_Mustache_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"ClusterToSeed_dEtadPhi_Mustache"), 0., 1.);   
           drawProfile2(ClusterToSeed_dEtadPhi_DeepSC_vs_seedEt_seedEta[etBin][etaBin], std::string("#Delta#Phi"), std::string("#Delta#eta"), std::string(""), std::string(output_dir+"ClusterToSeed_dEtadPhi_DeepSC"), 0., 1.);      
       }  
   } 
}

