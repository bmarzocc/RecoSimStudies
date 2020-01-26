#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include <algorithm> 
#include <iostream>
#include "TStyle.h"

double my2sideCrystalBall(double* x, double* par);
void drawHisto(TH1F* h_old, TH1F* h_new, std::string x_label, std::string drawType, std::string Name, bool log);
void drawEfficiency(TEfficiency* eff_SuperCluster, TEfficiency* eff_DeepSuperCluster, std::string xtitle, std::string Name);
void drawProfile(TProfile* prof_SuperCluster, TProfile* prof_DeepSuperCluster, std::string xtitle, std::string ytitle, std::string Name);

void SCAlgo_ValidationPlots(){

   gStyle->SetOptStat(1110); 
   //gStyle->SetOptStat(0000); 

   TFile* inFile = TFile::Open("/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourGammasGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_and_DeepSC_v2_showervars_Dumper.root");
   if (inFile == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open files!\n");
      return;
   }

   TH1F* h_Energy_EB_old = new TH1F("h_Energy_EB_old","SC_Energy_EB",150,0.,300.);
   TH1F* h_Energy_EE_old = new TH1F("h_Energy_EE_old","SC_Energy_EE",150,0.,300.);
   TH1F* h_Energy_EB_new = new TH1F("h_Energy_EB_new","SC_Energy_EB",150,0.,300.);
   TH1F* h_Energy_EE_new = new TH1F("h_Energy_EE_new","SC_Energy_EE",150,0.,300.);

   TH1F* h_EoEtrue_EB_old = new TH1F("h_EoEtrue_EB_old","SC_EoEtrue_EB",100,0.7,1.3);
   TH1F* h_EoEtrue_EE_old = new TH1F("h_EoEtrue_EE_old","SC_EoEtrue_EE",100,0.5,1.5);
   TH1F* h_EoEtrue_EB_new = new TH1F("h_EoEtrue_EB_new","SC_EoEtrue_EB",100,0.7,1.3);
   TH1F* h_EoEtrue_EE_new = new TH1F("h_EoEtrue_EE_new","SC_EoEtrue_EE",100,0.5,1.5);

   TH1F* h_Eta_Calo_Denum = new TH1F("h_Eta_Calo_Denum","h_Eta_Calo_Denum",100,-3.1,3.1);
   TH1F* h_Eta_Calo_SuperCluster_Num = new TH1F("h_Eta_Calo_SuperCluster_Num","h_Eta_Calo_SuperCluster_Num",100,-3.1,3.1);
   TH1F* h_Eta_Calo_DeepSuperCluster_Num = new TH1F("h_Eta_Calo_DeepSuperCluster_Num","h_Eta_Calo_DeepSuperCluster_Num",100,-3.1,3.1);
   TH1F* h_Et_Calo_Denum = new TH1F("h_Et_Calo_Denum","h_Et_Calo_Denum",101,0.,101.);
   TH1F* h_Et_Calo_SuperCluster_Num = new TH1F("h_Et_Calo_SuperCluster_Num","h_Et_Calo_SuperCluster_Num",101,0.,101.);
   TH1F* h_Et_Calo_DeepSuperCluster_Num = new TH1F("h_Et_Calo_DeepSuperCluster_Num","h_Et_Calo_DeepSuperCluster_Num",101,0.,101.);

   TProfile* prof_EoEtrue_vs_Eta_Calo_old  = new TProfile("prof_EoEtrue_vs_Eta_Calo_old","prof_EoEtrue_vs_Eta_Calo",100,-3.1,3.1,0.5,1.5);
   TProfile* prof_EoEtrue_vs_Et_Calo_old  = new TProfile("prof_EoEtrue_vs_Et_Calo_old","prof_EoEtrue_vs_Et_Calo",101,0.,101.,0.5,1.5); 
   TProfile* prof_EoEtrue_vs_Eta_Calo_new  = new TProfile("prof_EoEtrue_vs_Eta_Calo_new","prof_EoEtrue_vs_Eta_Calo",100,-3.1,3.1,0.5,1.5);
   TProfile* prof_EoEtrue_vs_Et_Calo_new  = new TProfile("prof_EoEtrue_vs_Et_Calo_new","prof_EoEtrue_vs_Et_Calo",101,0.,101.,0.5,1.5); 
   TProfile* prof_EoEtrue_vs_nVtx_old  = new TProfile("prof_EoEtrue_vs_nVtx_old","prof_EoEtrue_vs_nVtx",100,0.,100.,0.5,1.5);
   TProfile* prof_EoEtrue_vs_Rho_old  = new TProfile("prof_EoEtrue_vs_Rho_old","prof_EoEtrue_vs_Rho",100,0.,100.,0.5,1.5); 
   TProfile* prof_EoEtrue_vs_nVtx_new  = new TProfile("prof_EoEtrue_vs_nVtx_new","prof_EoEtrue_vs_nVtx",100,0.,100.,0.5,1.5);
   TProfile* prof_EoEtrue_vs_Rho_new  = new TProfile("prof_EoEtrue_vs_Rho_new","prof_EoEtrue_vs_Rho",100,0.,100.,0.5,1.5); 

   TH1F* h_Eta_old = new TH1F("h_Eta_old","SC_#eta",100,-3.1,3.1);
   TH1F* h_Eta_new = new TH1F("h_Eta_new","SC_#eta",100,-3.1,3.1);
   TH1F* h_Phi_EB_old = new TH1F("h_Phi_EB_old","SC_#phi_EB",100,-3.3,3.3);
   TH1F* h_Phi_EE_old = new TH1F("h_Phi_EE_old","SC_#phi_EE",100,-3.3,3.3);
   TH1F* h_Phi_EB_new = new TH1F("h_Phi_EB_new","SC_#phi_EB",100,-3.3,3.3);
   TH1F* h_Phi_EE_new = new TH1F("h_Phi_EE_new","SC_#phi_EE",100,-3.3,3.3);

   TH1F* h_EtaWidth_EB_old = new TH1F("h_EtaWidth_EB_old","SC_#etaWidth_EB",100,0.,0.1);
   TH1F* h_EtaWidth_EE_old = new TH1F("h_EtaWidth_EE_old","SC_#etaWidth_EE",100,0.,0.1);
   TH1F* h_EtaWidth_EB_new = new TH1F("h_EtaWidth_EB_new","SC_#etaWidth_EB",100,0.,0.1);
   TH1F* h_EtaWidth_EE_new = new TH1F("h_EtaWidth_EE_new","SC_#etaWidth_EE",100,0.,0.1); 

   TH1F* h_PhiWidth_EB_old = new TH1F("h_PhiWidth_EB_old","SC_#phiWidth_EB",100,0.,0.5);
   TH1F* h_PhiWidth_EE_old = new TH1F("h_PhiWidth_EE_old","SC_#phiWidth_EE",100,0.,0.5);
   TH1F* h_PhiWidth_EB_new = new TH1F("h_PhiWidth_EB_new","SC_#phiWidth_EB",100,0.,0.5);
   TH1F* h_PhiWidth_EE_new = new TH1F("h_PhiWidth_EE_new","SC_#phiWidth_EE",100,0.,0.5); 

   TH1F* h_nPFClusters_old = new TH1F("h_nPFClusters_old","SC_nPFClusters",15,0.,15.);
   TH1F* h_nPFClusters_new = new TH1F("h_nPFClusters_new","SC_nPFClusters",15,0.,15.);
   TH1F* h_nPFClusters_EB_old = new TH1F("h_nPFClusters_EB_old","SC_nPFClusters_EB",15,0.,15.);
   TH1F* h_nPFClusters_EE_old = new TH1F("h_nPFClusters_EE_old","SC_nPFClusters_EE",15,0.,15.);
   TH1F* h_nPFClusters_EB_new = new TH1F("h_nPFClusters_EB_new","SC_nPFClusters_EB",15,0.,15.);
   TH1F* h_nPFClusters_EE_new = new TH1F("h_nPFClusters_EE_new","SC_nPFClusters_EE",15,0.,15.);

   TH1F* h_R9_EB_old = new TH1F("h_R9_EB_old","SC_R9_EB",200,0.,1.2);
   TH1F* h_R9_EE_old = new TH1F("h_R9_EE_old","SC_R9_EE",200,0.,1.2);
   TH1F* h_R9_EB_new = new TH1F("h_R9_EB_new","SC_R9_EB",200,0.,1.2);
   TH1F* h_R9_EE_new = new TH1F("h_R9_EE_new","SC_R9_EE",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_old = new TH1F("h_full5x5_R9_EB_old","SC_R9_EB",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_old = new TH1F("h_full5x5_R9_EE_old","SC_R9_EE",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_new = new TH1F("h_full5x5_R9_EB_new","SC_R9_EB",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_new = new TH1F("h_full5x5_R9_EE_new","SC_R9_EE",200,0.,1.2);

   TH1F* h_sigmaIetaIeta_EB_old = new TH1F("h_sigmaIetaIeta_EB_old","SC_#sigmai#etai#eta_EB",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_old = new TH1F("h_sigmaIetaIeta_EE_old","SC_#sigmai#etai#eta_EE",200,0.,0.04);
   TH1F* h_sigmaIetaIeta_EB_new = new TH1F("h_sigmaIetaIeta_EB_new","SC_#sigmai#etai#eta_EB",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_new = new TH1F("h_sigmaIetaIeta_EE_new","SC_#sigmai#etai#eta_EE",200,0.,0.04);
   TH1F* h_full5x5_sigmaIetaIeta_EB_old = new TH1F("h_full5x5_sigmaIetaIeta_EB_old","SC_full5x5_#sigmai#etai#eta_EB",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_old = new TH1F("h_full5x5_sigmaIetaIeta_EE_old","SC_full5x5_#sigmai#etai#eta_EE",400,0.,0.08);
   TH1F* h_full5x5_sigmaIetaIeta_EB_new = new TH1F("h_full5x5_sigmaIetaIeta_EB_new","SC_full5x5_#sigmai#etai#eta_EB",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_new = new TH1F("h_full5x5_sigmaIetaIeta_EE_new","SC_full5x5_#sigmai#etai#eta_EE",400,0.,0.08);

   TH1F* h_sigmaIetaIphi_EB_old = new TH1F("h_sigmaIetaIphi_EB_old","SC_#sigmai#etai#phi_EB",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_old = new TH1F("h_sigmaIetaIphi_EE_old","SC_#sigmai#etai#phi_EE",400,0.,0.005);
   TH1F* h_sigmaIetaIphi_EB_new = new TH1F("h_sigmaIetaIphi_EB_new","SC_#sigmai#etai#phi_EB",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_new = new TH1F("h_sigmaIetaIphi_EE_new","SC_#sigmai#etai#phi_EE",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_old = new TH1F("h_full5x5_sigmaIetaIphi_EB_old","SC_full5x5_#sigmai#etai#phi_EB",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_old = new TH1F("h_full5x5_sigmaIetaIphi_EE_old","SC_full5x5_#sigmai#etai#phi_EE",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_new = new TH1F("h_full5x5_sigmaIetaIphi_EB_new","SC_full5x5_#sigmai#etai#phi_EB",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_new = new TH1F("h_full5x5_sigmaIetaIphi_EE_new","SC_full5x5_#sigmai#etai#phi_EE",400,0.,0.005);

   TH1F* h_sigmaIphiIphi_EB_old = new TH1F("h_sigmaIphiIphi_EB_old","SC_#sigmai#phii#phi_EB",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_old = new TH1F("h_sigmaIphiIphi_EE_old","SC_#sigmai#phii#phi_EE",300,0.,0.08);
   TH1F* h_sigmaIphiIphi_EB_new = new TH1F("h_sigmaIphiIphi_EB_new","SC_#sigmai#phii#phi_EB",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_new = new TH1F("h_sigmaIphiIphi_EE_new","SC_#sigmai#phii#phi_EE",300,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_old = new TH1F("h_full5x5_sigmaIphiIphi_EB_old","SC_full5x5_#sigmai#phii#phi_EB",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_old = new TH1F("h_full5x5_sigmaIphiIphi_EE_old","SC_full5x5_#sigmai#phii#phi_EE",400,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_new = new TH1F("h_full5x5_sigmaIphiIphi_EB_new","SC_full5x5_#sigmai#phii#phi_EB",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_new = new TH1F("h_full5x5_sigmaIphiIphi_EE_new","SC_full5x5_#sigmai#phii#phi_EE",400,0.,0.08);

   //seed matched plots
   TH1F* h_Energy_EB_seedMatched_old = new TH1F("h_Energy_EB_seedMatched_old","SC_Energy_EB_seedMatched",150,0.,300.);
   TH1F* h_Energy_EE_seedMatched_old = new TH1F("h_Energy_EE_seedMatched_old","SC_Energy_EE_seedMatched",150,0.,300.);
   TH1F* h_Energy_EB_seedMatched_new = new TH1F("h_Energy_EB_seedMatched_new","SC_Energy_EB_seedMatched",150,0.,300.);
   TH1F* h_Energy_EE_seedMatched_new = new TH1F("h_Energy_EE_seedMatched_new","SC_Energy_EE_seedMatched",150,0.,300.);

   TH1F* h_EoEtrue_EB_seedMatched_old = new TH1F("h_EoEtrue_EB_seedMatched_old","SC_EoEtrue_EB_seedMatched",100,0.8,1.2);
   TH1F* h_EoEtrue_EE_seedMatched_old = new TH1F("h_EoEtrue_EE_seedMatched_old","SC_EoEtrue_EE_seedMatched",100,0.8,1.5);
   TH1F* h_EoEtrue_EB_seedMatched_new = new TH1F("h_EoEtrue_EB_seedMatched_new","SC_EoEtrue_EB_seedMatched",100,0.8,1.2);
   TH1F* h_EoEtrue_EE_seedMatched_new = new TH1F("h_EoEtrue_EE_seedMatched_new","SC_EoEtrue_EE_seedMatched",100,0.8,1.5);

   TH1F* h_Eta_seedMatched_old = new TH1F("h_Eta_seedMatched_old","SC_#eta_seedMatched",100,-3.1,3.1);
   TH1F* h_Eta_seedMatched_new = new TH1F("h_Eta_seedMatched_new","SC_#eta_seedMatched",100,-3.1,3.1);
   TH1F* h_Phi_EB_seedMatched_old = new TH1F("h_Phi_EB_seedMatched_old","SC_#phi_EB_seedMatched",100,-3.3,3.3);
   TH1F* h_Phi_EE_seedMatched_old = new TH1F("h_Phi_EE_seedMatched_old","SC_#phi_EE_seedMatched",100,-3.3,3.3);
   TH1F* h_Phi_EB_seedMatched_new = new TH1F("h_Phi_EB_seedMatched_new","SC_#phi_EB_seedMatched",100,-3.3,3.3);
   TH1F* h_Phi_EE_seedMatched_new = new TH1F("h_Phi_EE_seedMatched_new","SC_#phi_EE_seedMatched",100,-3.3,3.3);

   TH1F* h_EtaWidth_EB_seedMatched_old = new TH1F("h_EtaWidth_EB_seedMatched_old","SC_#etaWidth_EB_seedMatched",100,0.,0.1);
   TH1F* h_EtaWidth_EE_seedMatched_old = new TH1F("h_EtaWidth_EE_seedMatched_old","SC_#etaWidth_EE_seedMatched",100,0.,0.1);
   TH1F* h_EtaWidth_EB_seedMatched_new = new TH1F("h_EtaWidth_EB_seedMatched_new","SC_#etaWidth_EB_seedMatched",100,0.,0.1);
   TH1F* h_EtaWidth_EE_seedMatched_new = new TH1F("h_EtaWidth_EE_seedMatched_new","SC_#etaWidth_EE_seedMatched",100,0.,0.1); 

   TH1F* h_PhiWidth_EB_seedMatched_old = new TH1F("h_PhiWidth_EB_seedMatched_old","SC_#phiWidth_EB_seedMatched",100,0.,0.5);
   TH1F* h_PhiWidth_EE_seedMatched_old = new TH1F("h_PhiWidth_EE_seedMatched_old","SC_#phiWidth_EE_seedMatched",100,0.,0.5);
   TH1F* h_PhiWidth_EB_seedMatched_new = new TH1F("h_PhiWidth_EB_seedMatched_new","SC_#phiWidth_EB_seedMatched",100,0.,0.5);
   TH1F* h_PhiWidth_EE_seedMatched_new = new TH1F("h_PhiWidth_EE_seedMatched_new","SC_#phiWidth_EE_seedMatched",100,0.,0.5); 

   TH1F* h_nPFClusters_seedMatched_old = new TH1F("h_nPFClusters_seedMatched_old","SC_nPFClusters_seedMatched",15,0.,15.);
   TH1F* h_nPFClusters_seedMatched_new = new TH1F("h_nPFClusters_seedMatched_new","SC_nPFClusters_seedMatched",15,0.,15.);
   TH1F* h_nPFClusters_EB_seedMatched_old = new TH1F("h_nPFClusters_EB_seedMatched_old","SC_nPFClusters_EB_seedMatched",15,0.,15.);
   TH1F* h_nPFClusters_EE_seedMatched_old = new TH1F("h_nPFClusters_EE_seedMatched_old","SC_nPFClusters_EE_seedMatched",15,0.,15.);
   TH1F* h_nPFClusters_EB_seedMatched_new = new TH1F("h_nPFClusters_EB_seedMatched_new","SC_nPFClusters_EB_seedMatched",15,0.,15.);
   TH1F* h_nPFClusters_EE_seedMatched_new = new TH1F("h_nPFClusters_EE_seedMatched_new","SC_nPFClusters_EE_seedMatched",15,0.,15.);

   TH1F* h_R9_EB_seedMatched_old = new TH1F("h_R9_EB_seedMatched_old","SC_R9_EB_seedMatched",200,0.,1.2);
   TH1F* h_R9_EE_seedMatched_old = new TH1F("h_R9_EE_seedMatched_old","SC_R9_EE_seedMatched",200,0.,1.2);
   TH1F* h_R9_EB_seedMatched_new = new TH1F("h_R9_EB_seedMatched_new","SC_R9_EB_seedMatched",200,0.,1.2);
   TH1F* h_R9_EE_seedMatched_new = new TH1F("h_R9_EE_seedMatched_new","SC_R9_EE_seedMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_seedMatched_old = new TH1F("h_full5x5_R9_EB_seedMatched_old","SC_R9_EB_seedMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_seedMatched_old = new TH1F("h_full5x5_R9_EE_seedMatched_old","SC_R9_EE_seedMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_seedMatched_new = new TH1F("h_full5x5_R9_EB_seedMatched_new","SC_R9_EB_seedMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_seedMatched_new = new TH1F("h_full5x5_R9_EE_seedMatched_new","SC_R9_EE_seedMatched",200,0.,1.2);

   TH1F* h_sigmaIetaIeta_EB_seedMatched_old = new TH1F("h_sigmaIetaIeta_EB_seedMatched_old","SC_#sigmai#etai#eta_EB_seedMatched",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_seedMatched_old = new TH1F("h_sigmaIetaIeta_EE_seedMatched_old","SC_#sigmai#etai#eta_EE_seedMatched",200,0.,0.04);
   TH1F* h_sigmaIetaIeta_EB_seedMatched_new = new TH1F("h_sigmaIetaIeta_EB_seedMatched_new","SC_#sigmai#etai#eta_EB_seedMatched",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_seedMatched_new = new TH1F("h_sigmaIetaIeta_EE_seedMatched_new","SC_#sigmai#etai#eta_EE_seedMatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIetaIeta_EB_seedMatched_old = new TH1F("h_full5x5_sigmaIetaIeta_EB_seedMatched_old","SC_full5x5_#sigmai#etai#eta_EB_seedMatched",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_seedMatched_old = new TH1F("h_full5x5_sigmaIetaIeta_EE_seedMatched_old","SC_full5x5_#sigmai#etai#eta_EE_seedMatched",400,0.,0.08);
   TH1F* h_full5x5_sigmaIetaIeta_EB_seedMatched_new = new TH1F("h_full5x5_sigmaIetaIeta_EB_seedMatched_new","SC_full5x5_#sigmai#etai#eta_EB_seedMatched",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_seedMatched_new = new TH1F("h_full5x5_sigmaIetaIeta_EE_seedMatched_new","SC_full5x5_#sigmai#etai#eta_EE_seedMatched",400,0.,0.08);

   TH1F* h_sigmaIetaIphi_EB_seedMatched_old = new TH1F("h_sigmaIetaIphi_EB_seedMatched_old","SC_#sigmai#etai#phi_EB_seedMatched",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_seedMatched_old = new TH1F("h_sigmaIetaIphi_EE_seedMatched_old","SC_#sigmai#etai#phi_EE_seedMatched",400,0.,0.005);
   TH1F* h_sigmaIetaIphi_EB_seedMatched_new = new TH1F("h_sigmaIetaIphi_EB_seedMatched_new","SC_#sigmai#etai#phi_EB_seedMatched",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_seedMatched_new = new TH1F("h_sigmaIetaIphi_EE_seedMatched_new","SC_#sigmai#etai#phi_EE_seedMatched",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_seedMatched_old = new TH1F("h_full5x5_sigmaIetaIphi_EB_seedMatched_old","SC_full5x5_#sigmai#etai#phi_EB_seedMatched",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_seedMatched_old = new TH1F("h_full5x5_sigmaIetaIphi_EE_seedMatched_old","SC_full5x5_#sigmai#etai#phi_EE_seedMatched",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_seedMatched_new = new TH1F("h_full5x5_sigmaIetaIphi_EB_seedMatched_new","SC_full5x5_#sigmai#etai#phi_EB_seedMatched",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_seedMatched_new = new TH1F("h_full5x5_sigmaIetaIphi_EE_seedMatched_new","SC_full5x5_#sigmai#etai#phi_EE_seedMatched",400,0.,0.005);

   TH1F* h_sigmaIphiIphi_EB_seedMatched_old = new TH1F("h_sigmaIphiIphi_EB_seedMatched_old","SC_#sigmai#phii#phi_EB_seedMatched",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_seedMatched_old = new TH1F("h_sigmaIphiIphi_EE_seedMatched_old","SC_#sigmai#phii#phi_EE_seedMatched",300,0.,0.08);
   TH1F* h_sigmaIphiIphi_EB_seedMatched_new = new TH1F("h_sigmaIphiIphi_EB_seedMatched_new","SC_#sigmai#phii#phi_EB_seedMatched",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_seedMatched_new = new TH1F("h_sigmaIphiIphi_EE_seedMatched_new","SC_#sigmai#phii#phi_EE_seedMatched",300,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_seedMatched_old = new TH1F("h_full5x5_sigmaIphiIphi_EB_seedMatched_old","SC_full5x5_#sigmai#phii#phi_EB_seedMatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_seedMatched_old = new TH1F("h_full5x5_sigmaIphiIphi_EE_seedMatched_old","SC_full5x5_#sigmai#phii#phi_EE_seedMatched",400,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_seedMatched_new = new TH1F("h_full5x5_sigmaIphiIphi_EB_seedMatched_new","SC_full5x5_#sigmai#phii#phi_EB_seedMatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_seedMatched_new = new TH1F("h_full5x5_sigmaIphiIphi_EE_seedMatched_new","SC_full5x5_#sigmai#phii#phi_EE_seedMatched",400,0.,0.08);

   //caloMatched
   TH1F* h_Energy_EB_caloMatched_old = new TH1F("h_Energy_EB_caloMatched_old","SC_Energy_EB_caloMatched",150,0.,300.);
   TH1F* h_Energy_EE_caloMatched_old = new TH1F("h_Energy_EE_caloMatched_old","SC_Energy_EE_caloMatched",150,0.,300.);
   TH1F* h_Energy_EB_caloMatched_new = new TH1F("h_Energy_EB_caloMatched_new","SC_Energy_EB_caloMatched",150,0.,300.);
   TH1F* h_Energy_EE_caloMatched_new = new TH1F("h_Energy_EE_caloMatched_new","SC_Energy_EE_caloMatched",150,0.,300.);

   TH1F* h_EoEtrue_EB_caloMatched_old = new TH1F("h_EoEtrue_EB_caloMatched_old","SC_EoEtrue_EB_caloMatched",100,0.8,1.2);
   TH1F* h_EoEtrue_EE_caloMatched_old = new TH1F("h_EoEtrue_EE_caloMatched_old","SC_EoEtrue_EE_caloMatched",100,0.8,1.5);
   TH1F* h_EoEtrue_EB_caloMatched_new = new TH1F("h_EoEtrue_EB_caloMatched_new","SC_EoEtrue_EB_caloMatched",100,0.8,1.2);
   TH1F* h_EoEtrue_EE_caloMatched_new = new TH1F("h_EoEtrue_EE_caloMatched_new","SC_EoEtrue_EE_caloMatched",100,0.8,1.5);

   TH1F* h_Eta_caloMatched_old = new TH1F("h_Eta_caloMatched_old","SC_#eta_caloMatched",100,-3.1,3.1);
   TH1F* h_Eta_caloMatched_new = new TH1F("h_Eta_caloMatched_new","SC_#eta_caloMatched",100,-3.1,3.1);
   TH1F* h_Phi_EB_caloMatched_old = new TH1F("h_Phi_EB_caloMatched_old","SC_#phi_EB_caloMatched",100,-3.3,3.3);
   TH1F* h_Phi_EE_caloMatched_old = new TH1F("h_Phi_EE_caloMatched_old","SC_#phi_EE_caloMatched",100,-3.3,3.3);
   TH1F* h_Phi_EB_caloMatched_new = new TH1F("h_Phi_EB_caloMatched_new","SC_#phi_EB_caloMatched",100,-3.3,3.3);
   TH1F* h_Phi_EE_caloMatched_new = new TH1F("h_Phi_EE_caloMatched_new","SC_#phi_EE_caloMatched",100,-3.3,3.3);

   TH1F* h_EtaWidth_EB_caloMatched_old = new TH1F("h_EtaWidth_EB_caloMatched_old","SC_#etaWidth_EB_caloMatched",100,0.,0.1);
   TH1F* h_EtaWidth_EE_caloMatched_old = new TH1F("h_EtaWidth_EE_caloMatched_old","SC_#etaWidth_EE_caloMatched",100,0.,0.1);
   TH1F* h_EtaWidth_EB_caloMatched_new = new TH1F("h_EtaWidth_EB_caloMatched_new","SC_#etaWidth_EB_caloMatched",100,0.,0.1);
   TH1F* h_EtaWidth_EE_caloMatched_new = new TH1F("h_EtaWidth_EE_caloMatched_new","SC_#etaWidth_EE_caloMatched",100,0.,0.1); 

   TH1F* h_PhiWidth_EB_caloMatched_old = new TH1F("h_PhiWidth_EB_caloMatched_old","SC_#phiWidth_EB_caloMatched",100,0.,0.5);
   TH1F* h_PhiWidth_EE_caloMatched_old = new TH1F("h_PhiWidth_EE_caloMatched_old","SC_#phiWidth_EE_caloMatched",100,0.,0.5);
   TH1F* h_PhiWidth_EB_caloMatched_new = new TH1F("h_PhiWidth_EB_caloMatched_new","SC_#phiWidth_EB_caloMatched",100,0.,0.5);
   TH1F* h_PhiWidth_EE_caloMatched_new = new TH1F("h_PhiWidth_EE_caloMatched_new","SC_#phiWidth_EE_caloMatched",100,0.,0.5); 

   TH1F* h_nPFClusters_caloMatched_old = new TH1F("h_nPFClusters_caloMatched_old","SC_nPFClusters_caloMatched",15,0.,15.);
   TH1F* h_nPFClusters_caloMatched_new = new TH1F("h_nPFClusters_caloMatched_new","SC_nPFClusters_caloMatched",15,0.,15.);
   TH1F* h_nPFClusters_EB_caloMatched_old = new TH1F("h_nPFClusters_EB_caloMatched_old","SC_nPFClusters_EB_caloMatched",15,0.,15.);
   TH1F* h_nPFClusters_EE_caloMatched_old = new TH1F("h_nPFClusters_EE_caloMatched_old","SC_nPFClusters_EE_caloMatched",15,0.,15.);
   TH1F* h_nPFClusters_EB_caloMatched_new = new TH1F("h_nPFClusters_EB_caloMatched_new","SC_nPFClusters_EB_caloMatched",15,0.,15.);
   TH1F* h_nPFClusters_EE_caloMatched_new = new TH1F("h_nPFClusters_EE_caloMatched_new","SC_nPFClusters_EE_caloMatched",15,0.,15.);

   TH1F* h_R9_EB_caloMatched_old = new TH1F("h_R9_EB_caloMatched_old","SC_R9_EB_caloMatched",200,0.,1.2);
   TH1F* h_R9_EE_caloMatched_old = new TH1F("h_R9_EE_caloMatched_old","SC_R9_EE_caloMatched",200,0.,1.2);
   TH1F* h_R9_EB_caloMatched_new = new TH1F("h_R9_EB_caloMatched_new","SC_R9_EB_caloMatched",200,0.,1.2);
   TH1F* h_R9_EE_caloMatched_new = new TH1F("h_R9_EE_caloMatched_new","SC_R9_EE_caloMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_caloMatched_old = new TH1F("h_full5x5_R9_EB_caloMatched_old","SC_R9_EB_caloMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_caloMatched_old = new TH1F("h_full5x5_R9_EE_caloMatched_old","SC_R9_EE_caloMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_caloMatched_new = new TH1F("h_full5x5_R9_EB_caloMatched_new","SC_R9_EB_caloMatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_caloMatched_new = new TH1F("h_full5x5_R9_EE_caloMatched_new","SC_R9_EE_caloMatched",200,0.,1.2);

   TH1F* h_sigmaIetaIeta_EB_caloMatched_old = new TH1F("h_sigmaIetaIeta_EB_caloMatched_old","SC_#sigmai#etai#eta_EB_caloMatched",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_caloMatched_old = new TH1F("h_sigmaIetaIeta_EE_caloMatched_old","SC_#sigmai#etai#eta_EE_caloMatched",200,0.,0.04);
   TH1F* h_sigmaIetaIeta_EB_caloMatched_new = new TH1F("h_sigmaIetaIeta_EB_caloMatched_new","SC_#sigmai#etai#eta_EB_caloMatched",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_caloMatched_new = new TH1F("h_sigmaIetaIeta_EE_caloMatched_new","SC_#sigmai#etai#eta_EE_caloMatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloMatched_old = new TH1F("h_full5x5_sigmaIetaIeta_EB_caloMatched_old","SC_full5x5_#sigmai#etai#eta_EB_caloMatched",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloMatched_old = new TH1F("h_full5x5_sigmaIetaIeta_EE_caloMatched_old","SC_full5x5_#sigmai#etai#eta_EE_caloMatched",400,0.,0.08);
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloMatched_new = new TH1F("h_full5x5_sigmaIetaIeta_EB_caloMatched_new","SC_full5x5_#sigmai#etai#eta_EB_caloMatched",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloMatched_new = new TH1F("h_full5x5_sigmaIetaIeta_EE_caloMatched_new","SC_full5x5_#sigmai#etai#eta_EE_caloMatched",400,0.,0.08);

   TH1F* h_sigmaIetaIphi_EB_caloMatched_old = new TH1F("h_sigmaIetaIphi_EB_caloMatched_old","SC_#sigmai#etai#phi_EB_caloMatched",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_caloMatched_old = new TH1F("h_sigmaIetaIphi_EE_caloMatched_old","SC_#sigmai#etai#phi_EE_caloMatched",400,0.,0.005);
   TH1F* h_sigmaIetaIphi_EB_caloMatched_new = new TH1F("h_sigmaIetaIphi_EB_caloMatched_new","SC_#sigmai#etai#phi_EB_caloMatched",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_caloMatched_new = new TH1F("h_sigmaIetaIphi_EE_caloMatched_new","SC_#sigmai#etai#phi_EE_caloMatched",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloMatched_old = new TH1F("h_full5x5_sigmaIetaIphi_EB_caloMatched_old","SC_full5x5_#sigmai#etai#phi_EB_caloMatched",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloMatched_old = new TH1F("h_full5x5_sigmaIetaIphi_EE_caloMatched_old","SC_full5x5_#sigmai#etai#phi_EE_caloMatched",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloMatched_new = new TH1F("h_full5x5_sigmaIetaIphi_EB_caloMatched_new","SC_full5x5_#sigmai#etai#phi_EB_caloMatched",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloMatched_new = new TH1F("h_full5x5_sigmaIetaIphi_EE_caloMatched_new","SC_full5x5_#sigmai#etai#phi_EE_caloMatched",400,0.,0.005);

   TH1F* h_sigmaIphiIphi_EB_caloMatched_old = new TH1F("h_sigmaIphiIphi_EB_caloMatched_old","SC_#sigmai#phii#phi_EB_caloMatched",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_caloMatched_old = new TH1F("h_sigmaIphiIphi_EE_caloMatched_old","SC_#sigmai#phii#phi_EE_caloMatched",300,0.,0.08);
   TH1F* h_sigmaIphiIphi_EB_caloMatched_new = new TH1F("h_sigmaIphiIphi_EB_caloMatched_new","SC_#sigmai#phii#phi_EB_caloMatched",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_caloMatched_new = new TH1F("h_sigmaIphiIphi_EE_caloMatched_new","SC_#sigmai#phii#phi_EE_caloMatched",300,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloMatched_old = new TH1F("h_full5x5_sigmaIphiIphi_EB_caloMatched_old","SC_full5x5_#sigmai#phii#phi_EB_caloMatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloMatched_old = new TH1F("h_full5x5_sigmaIphiIphi_EE_caloMatched_old","SC_full5x5_#sigmai#phii#phi_EE_caloMatched",400,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloMatched_new = new TH1F("h_full5x5_sigmaIphiIphi_EB_caloMatched_new","SC_full5x5_#sigmai#phii#phi_EB_caloMatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloMatched_new = new TH1F("h_full5x5_sigmaIphiIphi_EE_caloMatched_new","SC_full5x5_#sigmai#phii#phi_EE_caloMatched",400,0.,0.08);

   //caloUnmatched
   TH1F* h_Energy_EB_caloUnmatched_old = new TH1F("h_Energy_EB_caloUnmatched_old","SC_Energy_EB_caloUnmatched",150,0.,300.);
   TH1F* h_Energy_EE_caloUnmatched_old = new TH1F("h_Energy_EE_caloUnmatched_old","SC_Energy_EE_caloUnmatched",150,0.,300.);
   TH1F* h_Energy_EB_caloUnmatched_new = new TH1F("h_Energy_EB_caloUnmatched_new","SC_Energy_EB_caloUnmatched",150,0.,300.);
   TH1F* h_Energy_EE_caloUnmatched_new = new TH1F("h_Energy_EE_caloUnmatched_new","SC_Energy_EE_caloUnmatched",150,0.,300.);

   TH1F* h_EoEtrue_EB_caloUnmatched_old = new TH1F("h_EoEtrue_EB_caloUnmatched_old","SC_EoEtrue_EB_caloUnmatched",100,0.8,1.2);
   TH1F* h_EoEtrue_EE_caloUnmatched_old = new TH1F("h_EoEtrue_EE_caloUnmatched_old","SC_EoEtrue_EE_caloUnmatched",100,0.8,1.5);
   TH1F* h_EoEtrue_EB_caloUnmatched_new = new TH1F("h_EoEtrue_EB_caloUnmatched_new","SC_EoEtrue_EB_caloUnmatched",100,0.8,1.2);
   TH1F* h_EoEtrue_EE_caloUnmatched_new = new TH1F("h_EoEtrue_EE_caloUnmatched_new","SC_EoEtrue_EE_caloUnmatched",100,0.8,1.5);

   TH1F* h_Eta_caloUnmatched_old = new TH1F("h_Eta_caloUnmatched_old","SC_#eta_caloUnmatched",100,-3.1,3.1);
   TH1F* h_Eta_caloUnmatched_new = new TH1F("h_Eta_caloUnmatched_new","SC_#eta_caloUnmatched",100,-3.1,3.1);
   TH1F* h_Phi_EB_caloUnmatched_old = new TH1F("h_Phi_EB_caloUnmatched_old","SC_#phi_EB_caloUnmatched",100,-3.3,3.3);
   TH1F* h_Phi_EE_caloUnmatched_old = new TH1F("h_Phi_EE_caloUnmatched_old","SC_#phi_EE_caloUnmatched",100,-3.3,3.3);
   TH1F* h_Phi_EB_caloUnmatched_new = new TH1F("h_Phi_EB_caloUnmatched_new","SC_#phi_EB_caloUnmatched",100,-3.3,3.3);
   TH1F* h_Phi_EE_caloUnmatched_new = new TH1F("h_Phi_EE_caloUnmatched_new","SC_#phi_EE_caloUnmatched",100,-3.3,3.3);

   TH1F* h_EtaWidth_EB_caloUnmatched_old = new TH1F("h_EtaWidth_EB_caloUnmatched_old","SC_#etaWidth_EB_caloUnmatched",100,0.,0.1);
   TH1F* h_EtaWidth_EE_caloUnmatched_old = new TH1F("h_EtaWidth_EE_caloUnmatched_old","SC_#etaWidth_EE_caloUnmatched",100,0.,0.1);
   TH1F* h_EtaWidth_EB_caloUnmatched_new = new TH1F("h_EtaWidth_EB_caloUnmatched_new","SC_#etaWidth_EB_caloUnmatched",100,0.,0.1);
   TH1F* h_EtaWidth_EE_caloUnmatched_new = new TH1F("h_EtaWidth_EE_caloUnmatched_new","SC_#etaWidth_EE_caloUnmatched",100,0.,0.1); 

   TH1F* h_PhiWidth_EB_caloUnmatched_old = new TH1F("h_PhiWidth_EB_caloUnmatched_old","SC_#phiWidth_EB_caloUnmatched",100,0.,0.5);
   TH1F* h_PhiWidth_EE_caloUnmatched_old = new TH1F("h_PhiWidth_EE_caloUnmatched_old","SC_#phiWidth_EE_caloUnmatched",100,0.,0.5);
   TH1F* h_PhiWidth_EB_caloUnmatched_new = new TH1F("h_PhiWidth_EB_caloUnmatched_new","SC_#phiWidth_EB_caloUnmatched",100,0.,0.5);
   TH1F* h_PhiWidth_EE_caloUnmatched_new = new TH1F("h_PhiWidth_EE_caloUnmatched_new","SC_#phiWidth_EE_caloUnmatched",100,0.,0.5); 

   TH1F* h_nPFClusters_caloUnmatched_old = new TH1F("h_nPFClusters_caloUnmatched_old","SC_nPFClusters_caloUnmatched",15,0.,15.);
   TH1F* h_nPFClusters_caloUnmatched_new = new TH1F("h_nPFClusters_caloUnmatched_new","SC_nPFClusters_caloUnmatched",15,0.,15.);
   TH1F* h_nPFClusters_EB_caloUnmatched_old = new TH1F("h_nPFClusters_EB_caloUnmatched_old","SC_nPFClusters_EB_caloUnmatched",15,0.,15.);
   TH1F* h_nPFClusters_EE_caloUnmatched_old = new TH1F("h_nPFClusters_EE_caloUnmatched_old","SC_nPFClusters_EE_caloUnmatched",15,0.,15.);
   TH1F* h_nPFClusters_EB_caloUnmatched_new = new TH1F("h_nPFClusters_EB_caloUnmatched_new","SC_nPFClusters_EB_caloUnmatched",15,0.,15.);
   TH1F* h_nPFClusters_EE_caloUnmatched_new = new TH1F("h_nPFClusters_EE_caloUnmatched_new","SC_nPFClusters_EE_caloUnmatched",15,0.,15.);

   TH1F* h_R9_EB_caloUnmatched_old = new TH1F("h_R9_EB_caloUnmatched_old","SC_R9_EB_caloUnmatched",200,0.,1.2);
   TH1F* h_R9_EE_caloUnmatched_old = new TH1F("h_R9_EE_caloUnmatched_old","SC_R9_EE_caloUnmatched",200,0.,1.2);
   TH1F* h_R9_EB_caloUnmatched_new = new TH1F("h_R9_EB_caloUnmatched_new","SC_R9_EB_caloUnmatched",200,0.,1.2);
   TH1F* h_R9_EE_caloUnmatched_new = new TH1F("h_R9_EE_caloUnmatched_new","SC_R9_EE_caloUnmatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_caloUnmatched_old = new TH1F("h_full5x5_R9_EB_caloUnmatched_old","SC_R9_EB_caloUnmatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_caloUnmatched_old = new TH1F("h_full5x5_R9_EE_caloUnmatched_old","SC_R9_EE_caloUnmatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EB_caloUnmatched_new = new TH1F("h_full5x5_R9_EB_caloUnmatched_new","SC_R9_EB_caloUnmatched",200,0.,1.2);
   TH1F* h_full5x5_R9_EE_caloUnmatched_new = new TH1F("h_full5x5_R9_EE_caloUnmatched_new","SC_R9_EE_caloUnmatched",200,0.,1.2);

   TH1F* h_sigmaIetaIeta_EB_caloUnmatched_old = new TH1F("h_sigmaIetaIeta_EB_caloUnmatched_old","SC_#sigmai#etai#eta_EB_caloUnmatched",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_caloUnmatched_old = new TH1F("h_sigmaIetaIeta_EE_caloUnmatched_old","SC_#sigmai#etai#eta_EE_caloUnmatched",200,0.,0.04);
   TH1F* h_sigmaIetaIeta_EB_caloUnmatched_new = new TH1F("h_sigmaIetaIeta_EB_caloUnmatched_new","SC_#sigmai#etai#eta_EB_caloUnmatched",100,0.,0.02);
   TH1F* h_sigmaIetaIeta_EE_caloUnmatched_new = new TH1F("h_sigmaIetaIeta_EE_caloUnmatched_new","SC_#sigmai#etai#eta_EE_caloUnmatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old = new TH1F("h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old","SC_full5x5_#sigmai#etai#eta_EB_caloUnmatched",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old = new TH1F("h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old","SC_full5x5_#sigmai#etai#eta_EE_caloUnmatched",400,0.,0.08);
   TH1F* h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new = new TH1F("h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new","SC_full5x5_#sigmai#etai#eta_EB_caloUnmatched",100,0.,0.02);
   TH1F* h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new = new TH1F("h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new","SC_full5x5_#sigmai#etai#eta_EE_caloUnmatched",400,0.,0.08);

   TH1F* h_sigmaIetaIphi_EB_caloUnmatched_old = new TH1F("h_sigmaIetaIphi_EB_caloUnmatched_old","SC_#sigmai#etai#phi_EB_caloUnmatched",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_caloUnmatched_old = new TH1F("h_sigmaIetaIphi_EE_caloUnmatched_old","SC_#sigmai#etai#phi_EE_caloUnmatched",400,0.,0.005);
   TH1F* h_sigmaIetaIphi_EB_caloUnmatched_new = new TH1F("h_sigmaIetaIphi_EB_caloUnmatched_new","SC_#sigmai#etai#phi_EB_caloUnmatched",200,0.,0.001);
   TH1F* h_sigmaIetaIphi_EE_caloUnmatched_new = new TH1F("h_sigmaIetaIphi_EE_caloUnmatched_new","SC_#sigmai#etai#phi_EE_caloUnmatched",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old = new TH1F("h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old","SC_full5x5_#sigmai#etai#phi_EB_caloUnmatched",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old = new TH1F("h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old","SC_full5x5_#sigmai#etai#phi_EE_caloUnmatched",400,0.,0.005);
   TH1F* h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new = new TH1F("h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new","SC_full5x5_#sigmai#etai#phi_EB_caloUnmatched",200,0.,0.001);
   TH1F* h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new = new TH1F("h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new","SC_full5x5_#sigmai#etai#phi_EE_caloUnmatched",400,0.,0.005);

   TH1F* h_sigmaIphiIphi_EB_caloUnmatched_old = new TH1F("h_sigmaIphiIphi_EB_caloUnmatched_old","SC_#sigmai#phii#phi_EB_caloUnmatched",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_caloUnmatched_old = new TH1F("h_sigmaIphiIphi_EE_caloUnmatched_old","SC_#sigmai#phii#phi_EE_caloUnmatched",300,0.,0.08);
   TH1F* h_sigmaIphiIphi_EB_caloUnmatched_new = new TH1F("h_sigmaIphiIphi_EB_caloUnmatched_new","SC_#sigmai#phii#phi_EB_caloUnmatched",200,0.,0.04);
   TH1F* h_sigmaIphiIphi_EE_caloUnmatched_new = new TH1F("h_sigmaIphiIphi_EE_caloUnmatched_new","SC_#sigmai#phii#phi_EE_caloUnmatched",300,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old = new TH1F("h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old","SC_full5x5_#sigmai#phii#phi_EB_caloUnmatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old = new TH1F("h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old","SC_full5x5_#sigmai#phii#phi_EE_caloUnmatched",400,0.,0.08);
   TH1F* h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new = new TH1F("h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new","SC_full5x5_#sigmai#phii#phi_EB_caloUnmatched",200,0.,0.04);
   TH1F* h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new = new TH1F("h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new","SC_full5x5_#sigmai#phii#phi_EE_caloUnmatched",400,0.,0.08);
   
   int nVtx_tmp;
   float rho_tmp;
   vector<vector<int>> caloParticle_superCluster_simScore_MatchedIndex_tmp;
   vector<vector<int>> caloParticle_deepSuperCluster_simScore_MatchedIndex_tmp;
   vector<float> caloParticle_simEnergy_tmp;
   vector<float> caloParticle_simEta_tmp;
   vector<float> caloParticle_simEt_tmp; 
   vector<float> superCluster_energy_tmp;
   vector<float> superCluster_eta_tmp;
   vector<float> superCluster_phi_tmp;
   vector<float> superCluster_etaWidth_tmp;
   vector<float> superCluster_phiWidth_tmp;
   vector<int> superCluster_nPFClusters_tmp;
   vector<int> superCluster_seedIndex_tmp;
   vector<float> superCluster_R9_tmp;
   vector<float> superCluster_sigmaIetaIeta_tmp;
   vector<float> superCluster_sigmaIetaIphi_tmp;
   vector<float> superCluster_sigmaIphiIphi_tmp;
   vector<float> superCluster_full5x5_R9_tmp;
   vector<float> superCluster_full5x5_sigmaIetaIeta_tmp;
   vector<float> superCluster_full5x5_sigmaIetaIphi_tmp;
   vector<float> superCluster_full5x5_sigmaIphiIphi_tmp;
   vector<int> superCluster_simScore_MatchedIndex_tmp;
   vector<vector<double>> superCluster_simScore_tmp;
   vector<float> deepSuperCluster_energy_tmp;
   vector<float> deepSuperCluster_eta_tmp;
   vector<float> deepSuperCluster_phi_tmp;
   vector<float> deepSuperCluster_etaWidth_tmp;
   vector<float> deepSuperCluster_phiWidth_tmp;
   vector<int> deepSuperCluster_nPFClusters_tmp;
   vector<int> deepSuperCluster_seedIndex_tmp;
   vector<float> deepSuperCluster_R9_tmp;
   vector<float> deepSuperCluster_sigmaIetaIeta_tmp;
   vector<float> deepSuperCluster_sigmaIetaIphi_tmp;
   vector<float> deepSuperCluster_sigmaIphiIphi_tmp;
   vector<float> deepSuperCluster_full5x5_R9_tmp;
   vector<float> deepSuperCluster_full5x5_sigmaIetaIeta_tmp;
   vector<float> deepSuperCluster_full5x5_sigmaIetaIphi_tmp;
   vector<float> deepSuperCluster_full5x5_sigmaIphiIphi_tmp;
   vector<int> deepSuperCluster_simScore_MatchedIndex_tmp;
   vector<vector<double>> deepSuperCluster_simScore_tmp;
   vector<float> pfCluster_eta_tmp;
   vector<float> pfCluster_phi_tmp;
   
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
   
   // Create tyhe tree reader and its data containers
   TTreeReader myReader("recosimdumper/caloTree", inFile); 
   TTreeReaderValue<int> nVtx(myReader, "nVtx"); 
   TTreeReaderValue<float> rho(myReader, "rho");    
   TTreeReaderValue<vector<vector<int>>> caloParticle_superCluster_simScore_MatchedIndex(myReader, "caloParticle_superCluster_simScore_MatchedIndex");
   TTreeReaderValue<vector<vector<int>>> caloParticle_deepSuperCluster_simScore_MatchedIndex(myReader, "caloParticle_deepSuperCluster_simScore_MatchedIndex");
   TTreeReaderValue<vector<float>> caloParticle_simEnergy(myReader, "caloParticle_simEnergy");
   TTreeReaderValue<vector<float>> caloParticle_simEta(myReader, "caloParticle_simEta");
   TTreeReaderValue<vector<float>> caloParticle_simEt(myReader, "caloParticle_simPt");     
   TTreeReaderValue<vector<float>> pfCluster_eta(myReader, "pfCluster_eta");
   TTreeReaderValue<vector<float>> pfCluster_phi(myReader, "pfCluster_phi");
   TTreeReaderValue<vector<float>> superCluster_energy_old(myReader, "superCluster_energy");
   TTreeReaderValue<vector<float>> superCluster_eta_old(myReader, "superCluster_eta");
   TTreeReaderValue<vector<float>> superCluster_etaWidth_old(myReader, "superCluster_etaWidth");
   TTreeReaderValue<vector<float>> superCluster_phi_old(myReader, "superCluster_phi");
   TTreeReaderValue<vector<float>> superCluster_phiWidth_old(myReader, "superCluster_phiWidth");
   TTreeReaderValue<vector<int>> superCluster_nPFClusters_old(myReader, "superCluster_nPFClusters");
   TTreeReaderValue<vector<int>> superCluster_seedIndex_old(myReader, "superCluster_seedIndex");
   TTreeReaderValue<vector<float>> superCluster_R9_old(myReader, "superCluster_r9");
   TTreeReaderValue<vector<float>> superCluster_sigmaIetaIeta_old(myReader, "superCluster_sigmaIetaIeta");
   TTreeReaderValue<vector<float>> superCluster_sigmaIetaIphi_old(myReader, "superCluster_sigmaIetaIphi");
   TTreeReaderValue<vector<float>> superCluster_sigmaIphiIphi_old(myReader, "superCluster_sigmaIphiIphi");
   TTreeReaderValue<vector<float>> superCluster_full5x5_R9_old(myReader, "superCluster_full5x5_r9");
   TTreeReaderValue<vector<float>> superCluster_full5x5_sigmaIetaIeta_old(myReader, "superCluster_full5x5_sigmaIetaIeta");
   TTreeReaderValue<vector<float>> superCluster_full5x5_sigmaIetaIphi_old(myReader, "superCluster_full5x5_sigmaIetaIphi");
   TTreeReaderValue<vector<float>> superCluster_full5x5_sigmaIphiIphi_old(myReader, "superCluster_full5x5_sigmaIphiIphi");
   TTreeReaderValue<vector<int>> superCluster_simScore_MatchedIndex_old(myReader, "superCluster_simScore_MatchedIndex");
   TTreeReaderValue<vector<vector<double>>> superCluster_simScore_old(myReader, "superCluster_simScore");
   TTreeReaderValue<vector<float>> superCluster_energy_new(myReader, "deepSuperCluster_energy");
   TTreeReaderValue<vector<float>> superCluster_eta_new(myReader, "deepSuperCluster_eta");
   TTreeReaderValue<vector<float>> superCluster_etaWidth_new(myReader, "deepSuperCluster_etaWidth");
   TTreeReaderValue<vector<float>> superCluster_phi_new(myReader, "deepSuperCluster_phi");
   TTreeReaderValue<vector<float>> superCluster_phiWidth_new(myReader, "deepSuperCluster_phiWidth");
   TTreeReaderValue<vector<int>> superCluster_nPFClusters_new(myReader, "deepSuperCluster_nPFClusters");
   TTreeReaderValue<vector<int>> superCluster_seedIndex_new(myReader, "deepSuperCluster_seedIndex");
   TTreeReaderValue<vector<float>> superCluster_R9_new(myReader, "deepSuperCluster_r9");
   TTreeReaderValue<vector<float>> superCluster_sigmaIetaIeta_new(myReader, "deepSuperCluster_sigmaIetaIeta");
   TTreeReaderValue<vector<float>> superCluster_sigmaIetaIphi_new(myReader, "deepSuperCluster_sigmaIetaIphi");
   TTreeReaderValue<vector<float>> superCluster_sigmaIphiIphi_new(myReader, "deepSuperCluster_sigmaIphiIphi");
   TTreeReaderValue<vector<float>> superCluster_full5x5_R9_new(myReader, "deepSuperCluster_full5x5_r9");
   TTreeReaderValue<vector<float>> superCluster_full5x5_sigmaIetaIeta_new(myReader, "deepSuperCluster_full5x5_sigmaIetaIeta");
   TTreeReaderValue<vector<float>> superCluster_full5x5_sigmaIetaIphi_new(myReader, "deepSuperCluster_full5x5_sigmaIetaIphi");
   TTreeReaderValue<vector<float>> superCluster_full5x5_sigmaIphiIphi_new(myReader, "deepSuperCluster_full5x5_sigmaIphiIphi");
   TTreeReaderValue<vector<int>> superCluster_simScore_MatchedIndex_new(myReader, "deepSuperCluster_simScore_MatchedIndex");
   TTreeReaderValue<vector<vector<double>>> superCluster_simScore_new(myReader, "deepSuperCluster_simScore");
   
   // Loop over all entries of the TTree or TChain.
   int  entry = 0;
   while (myReader.Next()) {
       
      if(entry%100==0) std::cout<<"--- Reading tree = "<<entry<<std::endl;
      if(entry>20000) continue;

      nVtx_tmp = *nVtx;
      rho_tmp = *rho;
      caloParticle_superCluster_simScore_MatchedIndex_tmp = *caloParticle_superCluster_simScore_MatchedIndex;
      caloParticle_deepSuperCluster_simScore_MatchedIndex_tmp = *caloParticle_deepSuperCluster_simScore_MatchedIndex;
      caloParticle_simEnergy_tmp = *caloParticle_simEnergy;
      caloParticle_simEta_tmp = *caloParticle_simEta;
      caloParticle_simEt_tmp = *caloParticle_simEt;
      pfCluster_eta_tmp = *pfCluster_eta;
      pfCluster_phi_tmp = *pfCluster_phi;
      superCluster_energy_tmp = *superCluster_energy_old; 
      superCluster_eta_tmp = *superCluster_eta_old; 
      superCluster_etaWidth_tmp = *superCluster_etaWidth_old; 
      superCluster_phi_tmp = *superCluster_phi_old;   
      superCluster_phiWidth_tmp = *superCluster_phiWidth_old;    
      superCluster_nPFClusters_tmp = *superCluster_nPFClusters_old;   
      superCluster_seedIndex_tmp = *superCluster_seedIndex_old;     
      superCluster_R9_tmp = *superCluster_R9_old;   
      superCluster_full5x5_R9_tmp = *superCluster_full5x5_R9_old;  
      superCluster_sigmaIetaIeta_tmp = *superCluster_sigmaIetaIeta_old;   
      superCluster_full5x5_sigmaIetaIeta_tmp = *superCluster_full5x5_sigmaIetaIeta_old;   
      superCluster_sigmaIetaIphi_tmp = *superCluster_sigmaIetaIphi_old;   
      superCluster_full5x5_sigmaIetaIphi_tmp = *superCluster_full5x5_sigmaIetaIphi_old;       
      superCluster_sigmaIphiIphi_tmp = *superCluster_sigmaIphiIphi_old;   
      superCluster_full5x5_sigmaIphiIphi_tmp = *superCluster_full5x5_sigmaIphiIphi_old;          
      superCluster_simScore_MatchedIndex_tmp = *superCluster_simScore_MatchedIndex_old; 
      superCluster_simScore_tmp = *superCluster_simScore_old; 
      deepSuperCluster_energy_tmp = *superCluster_energy_new; 
      deepSuperCluster_eta_tmp = *superCluster_eta_new; 
      deepSuperCluster_etaWidth_tmp = *superCluster_etaWidth_new; 
      deepSuperCluster_phi_tmp = *superCluster_phi_new;   
      deepSuperCluster_phiWidth_tmp = *superCluster_phiWidth_new;    
      deepSuperCluster_nPFClusters_tmp = *superCluster_nPFClusters_new; 
      deepSuperCluster_seedIndex_tmp = *superCluster_seedIndex_new;       
      deepSuperCluster_R9_tmp = *superCluster_R9_new;   
      deepSuperCluster_full5x5_R9_tmp = *superCluster_full5x5_R9_new;  
      deepSuperCluster_sigmaIetaIeta_tmp = *superCluster_sigmaIetaIeta_new;   
      deepSuperCluster_full5x5_sigmaIetaIeta_tmp = *superCluster_full5x5_sigmaIetaIeta_new;   
      deepSuperCluster_sigmaIetaIphi_tmp = *superCluster_sigmaIetaIphi_new;   
      deepSuperCluster_full5x5_sigmaIetaIphi_tmp = *superCluster_full5x5_sigmaIetaIphi_new;       
      deepSuperCluster_sigmaIphiIphi_tmp = *superCluster_sigmaIphiIphi_new;   
      deepSuperCluster_full5x5_sigmaIphiIphi_tmp = *superCluster_full5x5_sigmaIphiIphi_new;          
      deepSuperCluster_simScore_MatchedIndex_tmp = *superCluster_simScore_MatchedIndex_new; 
      deepSuperCluster_simScore_tmp = *superCluster_simScore_new; 

      Calo_SuperCluster_simScore.clear();
      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){
          Calo_SuperCluster_simScore[iCalo]=0.;
          Calo_SuperCluster_index[iCalo]=-1;
      }
      Calo_DeepSuperCluster_simScore.clear();
      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){
          Calo_DeepSuperCluster_simScore[iCalo]=0.;
          Calo_DeepSuperCluster_index[iCalo]=-1;
      }

      simEnergy.resize(caloParticle_simEnergy_tmp.size());  
    
      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){
          for(unsigned int iSC=0; iSC<superCluster_simScore_tmp.size(); iSC++){
              if(superCluster_simScore_tmp.at(iSC).at(iCalo) > Calo_SuperCluster_simScore[iCalo]){
                 Calo_SuperCluster_simScore[iCalo]=superCluster_simScore_tmp.at(iSC).at(iCalo);
                 Calo_SuperCluster_index[iCalo]=iSC;
              }
          }
          for(unsigned int iSC=0; iSC<deepSuperCluster_simScore_tmp.size(); iSC++){
              if(deepSuperCluster_simScore_tmp.at(iSC).at(iCalo) > Calo_DeepSuperCluster_simScore[iCalo]){
                 Calo_DeepSuperCluster_simScore[iCalo]=deepSuperCluster_simScore_tmp.at(iSC).at(iCalo);
                 Calo_DeepSuperCluster_index[iCalo]=iSC;
              }
          }

          if(Calo_SuperCluster_index[iCalo]>=superCluster_energy_tmp.size() && Calo_SuperCluster_index[iCalo]>-1)std::cout << "WARNING SC---> iCalo = " << iCalo << " - " << Calo_SuperCluster_index[iCalo] << " - " << superCluster_energy_tmp.size() << std::endl;  
          if(Calo_DeepSuperCluster_index[iCalo]>=deepSuperCluster_energy_tmp.size() && Calo_DeepSuperCluster_index[iCalo]>-1)std::cout << "WARNING DeepSC iCalo = " << iCalo << " - " << Calo_DeepSuperCluster_index[iCalo] << " - " << deepSuperCluster_energy_tmp.size() << std::endl;  
      }

      superCluster_seeds.clear();
      for(unsigned int iPF=0; iPF<superCluster_energy_tmp.size(); iPF++)
          superCluster_seeds.push_back(superCluster_seedIndex_tmp.at(iPF));
      
      deepSuperCluster_seeds.clear();
      for(unsigned int iPF=0; iPF<deepSuperCluster_energy_tmp.size(); iPF++)
          deepSuperCluster_seeds.push_back(deepSuperCluster_seedIndex_tmp.at(iPF));

      common_seeds.clear(); 
      for(unsigned int iPF=0; iPF<superCluster_seeds.size(); iPF++)
          for(unsigned int iDeepPF=0; iDeepPF<deepSuperCluster_seeds.size(); iDeepPF++)
          {
              if(deepSuperCluster_seeds.at(iDeepPF)==superCluster_seeds.at(iPF)) common_seeds.push_back(superCluster_seeds.at(iPF));
          }    
 
      for(unsigned int iPF=0; iPF<superCluster_simScore_MatchedIndex_tmp.size(); iPF++){
          if(superCluster_simScore_MatchedIndex_tmp.at(iPF)>=0)   
             SuperCluster_RecoEnergy_simScore[superCluster_simScore_MatchedIndex_tmp.at(iPF)]+=superCluster_energy_tmp.at(iPF); 
          
          h_Eta_old->Fill(superCluster_eta_tmp.at(iPF));
          h_nPFClusters_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
          if(fabs(superCluster_eta_tmp.at(iPF))<1.442){ 
             h_EtaWidth_EB_old->Fill(superCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_old->Fill(superCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_old->Fill(superCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_old->Fill(superCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_old->Fill(superCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_old->Fill(superCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(superCluster_eta_tmp.at(iPF))>1.5){
             h_EtaWidth_EE_old->Fill(superCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_old->Fill(superCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_old->Fill(superCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_old->Fill(superCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_old->Fill(superCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));                             
          } 
         
          bool sameSeed=false;  
          for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
              if(superCluster_seedIndex_tmp.at(iPF) == common_seeds.at(iSeed)){
                 sameSeed=true;
                 break;
              }    
          
          if(sameSeed==false) continue;
           
          h_Eta_seedMatched_old->Fill(superCluster_eta_tmp.at(iPF));
          h_nPFClusters_seedMatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
          if(fabs(superCluster_eta_tmp.at(iPF))<1.442){ 
             h_EtaWidth_EB_seedMatched_old->Fill(superCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_seedMatched_old->Fill(superCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_seedMatched_old->Fill(superCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_seedMatched_old->Fill(superCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_seedMatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_seedMatched_old->Fill(superCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_seedMatched_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_seedMatched_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_seedMatched_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_seedMatched_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_seedMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(superCluster_eta_tmp.at(iPF))>1.5){
             h_EtaWidth_EE_seedMatched_old->Fill(superCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_seedMatched_old->Fill(superCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_seedMatched_old->Fill(superCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_seedMatched_old->Fill(superCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_seedMatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_seedMatched_old->Fill(superCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_seedMatched_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_seedMatched_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_seedMatched_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_seedMatched_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_seedMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));                             
          }  
            
      }

      for(unsigned int iPF=0; iPF<deepSuperCluster_simScore_MatchedIndex_tmp.size(); iPF++){
          if(deepSuperCluster_simScore_MatchedIndex_tmp.at(iPF)>=0)   
             DeepSuperCluster_RecoEnergy_simScore[deepSuperCluster_simScore_MatchedIndex_tmp.at(iPF)]+=deepSuperCluster_energy_tmp.at(iPF); 
          
          h_Eta_new->Fill(deepSuperCluster_eta_tmp.at(iPF));
          h_nPFClusters_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
          if(fabs(deepSuperCluster_eta_tmp.at(iPF))<1.442){ 
             h_EtaWidth_EB_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_new->Fill(deepSuperCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_new->Fill(deepSuperCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_new->Fill(deepSuperCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(deepSuperCluster_eta_tmp.at(iPF))>1.5){
             h_EtaWidth_EE_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_new->Fill(deepSuperCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_new->Fill(deepSuperCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_new->Fill(deepSuperCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));                             
          }  
         
          bool sameSeed=false;  
          for(unsigned int iSeed=0; iSeed<common_seeds.size(); iSeed++)
              if(deepSuperCluster_seedIndex_tmp.at(iPF) == common_seeds.at(iSeed)){
                 sameSeed=true;
                 break;
              }    
          
          if(sameSeed==false) continue;

          h_Eta_seedMatched_new->Fill(deepSuperCluster_eta_tmp.at(iPF));
          h_nPFClusters_seedMatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
          if(fabs(deepSuperCluster_eta_tmp.at(iPF))<1.442){ 
             h_EtaWidth_EB_seedMatched_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_seedMatched_new->Fill(deepSuperCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_seedMatched_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_seedMatched_new->Fill(deepSuperCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_seedMatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_seedMatched_new->Fill(deepSuperCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_seedMatched_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(deepSuperCluster_eta_tmp.at(iPF))>1.5){
             h_EtaWidth_EE_seedMatched_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_seedMatched_new->Fill(deepSuperCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_seedMatched_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_seedMatched_new->Fill(deepSuperCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_seedMatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_seedMatched_new->Fill(deepSuperCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_seedMatched_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_seedMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));                             
          }  
      }

      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){

          h_Eta_Calo_Denum->Fill(caloParticle_simEta_tmp.at(iCalo));
          h_Et_Calo_Denum->Fill(caloParticle_simEt_tmp.at(iCalo));
  
          int iPF=Calo_SuperCluster_index.at(iCalo); 
          if(iPF==-1) continue; 
          h_Eta_Calo_SuperCluster_Num->Fill(caloParticle_simEta_tmp.at(iCalo)); 
          h_Et_Calo_SuperCluster_Num->Fill(caloParticle_simEt_tmp.at(iCalo));  

          h_Eta_caloMatched_old->Fill(superCluster_eta_tmp.at(iPF));
          h_nPFClusters_caloMatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));
          prof_EoEtrue_vs_Eta_Calo_old->Fill(caloParticle_simEta_tmp.at(iCalo),superCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
          prof_EoEtrue_vs_Et_Calo_old->Fill(caloParticle_simEt_tmp.at(iCalo),superCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo)); 
          prof_EoEtrue_vs_nVtx_old->Fill(nVtx_tmp,superCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
          prof_EoEtrue_vs_Rho_old->Fill(rho_tmp,superCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));

          if(fabs(caloParticle_simEta_tmp.at(iCalo))<1.442){
             h_EoEtrue_EB_old->Fill(superCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
             h_EtaWidth_EB_caloMatched_old->Fill(superCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_caloMatched_old->Fill(superCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_caloMatched_old->Fill(superCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_caloMatched_old->Fill(superCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_caloMatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_caloMatched_old->Fill(superCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_caloMatched_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_caloMatched_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_caloMatched_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_caloMatched_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_caloMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF)); 
          }else if(fabs(caloParticle_simEta_tmp.at(iCalo))>1.5){
             h_EoEtrue_EE_old->Fill(superCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
             h_EtaWidth_EE_caloMatched_old->Fill(superCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_caloMatched_old->Fill(superCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_caloMatched_old->Fill(superCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_caloMatched_old->Fill(superCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_caloMatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_caloMatched_old->Fill(superCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_caloMatched_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_caloMatched_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_caloMatched_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_caloMatched_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_caloMatched_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));   
          }
      }
       
      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){

          int iPF=Calo_DeepSuperCluster_index.at(iCalo); 
          if(iPF==-1) continue; 

          h_Eta_Calo_DeepSuperCluster_Num->Fill(caloParticle_simEta_tmp.at(iCalo));
          h_Et_Calo_DeepSuperCluster_Num->Fill(caloParticle_simEt_tmp.at(iCalo));

          h_Eta_caloMatched_new->Fill(deepSuperCluster_eta_tmp.at(iPF));
          h_nPFClusters_caloMatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF)); 
          prof_EoEtrue_vs_Eta_Calo_new->Fill(caloParticle_simEta_tmp.at(iCalo),deepSuperCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
          prof_EoEtrue_vs_Et_Calo_new->Fill(caloParticle_simEt_tmp.at(iCalo),deepSuperCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo)); 
          prof_EoEtrue_vs_nVtx_new->Fill(nVtx_tmp,deepSuperCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
          prof_EoEtrue_vs_Rho_new->Fill(rho_tmp,deepSuperCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));   
 
          if(fabs(caloParticle_simEta_tmp.at(iCalo))<1.442){
             h_EoEtrue_EB_new->Fill(deepSuperCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
             h_EtaWidth_EB_caloMatched_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_caloMatched_new->Fill(deepSuperCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_caloMatched_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_caloMatched_new->Fill(deepSuperCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_caloMatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_caloMatched_new->Fill(deepSuperCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_caloMatched_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(caloParticle_simEta_tmp.at(iCalo))>1.5){
             h_EoEtrue_EE_new->Fill(deepSuperCluster_energy_tmp.at(iPF)/caloParticle_simEnergy_tmp.at(iCalo));
             h_EtaWidth_EE_caloMatched_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_caloMatched_new->Fill(deepSuperCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_caloMatched_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_caloMatched_new->Fill(deepSuperCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_caloMatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_caloMatched_new->Fill(deepSuperCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_caloMatched_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_caloMatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));  
          } 

      }

      for(unsigned int iPF=0; iPF<superCluster_energy_tmp.size(); iPF++){

          std::map<int,int>::iterator it = Calo_SuperCluster_index.find(iPF); 
          if(it != Calo_SuperCluster_index.end()) continue;   
 
          h_Eta_caloUnmatched_old->Fill(superCluster_eta_tmp.at(iPF));
          h_nPFClusters_caloUnmatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
          if(fabs(superCluster_eta_tmp.at(iPF))<1.442){ 
             h_EtaWidth_EB_caloUnmatched_old->Fill(superCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_caloUnmatched_old->Fill(superCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_caloUnmatched_old->Fill(superCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_caloUnmatched_old->Fill(superCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_caloUnmatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_caloUnmatched_old->Fill(superCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_caloUnmatched_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_caloUnmatched_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_caloUnmatched_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_caloUnmatched_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(superCluster_eta_tmp.at(iPF))>1.5){
             h_EtaWidth_EE_caloUnmatched_old->Fill(superCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_caloUnmatched_old->Fill(superCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_caloUnmatched_old->Fill(superCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_caloUnmatched_old->Fill(superCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_caloUnmatched_old->Fill(superCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_caloUnmatched_old->Fill(superCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_caloUnmatched_old->Fill(superCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_caloUnmatched_old->Fill(superCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_caloUnmatched_old->Fill(superCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_caloUnmatched_old->Fill(superCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old->Fill(superCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));                             
          } 
      }

      for(unsigned int iPF=0; iPF<deepSuperCluster_energy_tmp.size(); iPF++){
          std::map<int,int>::iterator it = Calo_DeepSuperCluster_index.find(iPF); 
          if(it != Calo_DeepSuperCluster_index.end()) continue;  
          h_Eta_caloUnmatched_new->Fill(deepSuperCluster_eta_tmp.at(iPF));
          h_nPFClusters_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
          if(fabs(deepSuperCluster_eta_tmp.at(iPF))<1.442){ 
             h_EtaWidth_EB_caloUnmatched_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));  
             h_Phi_EB_caloUnmatched_new->Fill(deepSuperCluster_phi_tmp.at(iPF)); 
             h_PhiWidth_EB_caloUnmatched_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF));  
             h_Energy_EB_caloUnmatched_new->Fill(deepSuperCluster_energy_tmp.at(iPF));  
             h_nPFClusters_EB_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));   
             h_R9_EB_caloUnmatched_new->Fill(deepSuperCluster_R9_tmp.at(iPF));    
             h_full5x5_R9_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));   
             h_sigmaIetaIeta_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));              
          }else if(fabs(deepSuperCluster_eta_tmp.at(iPF))>1.5){
             h_EtaWidth_EE_caloUnmatched_new->Fill(deepSuperCluster_etaWidth_tmp.at(iPF));
             h_Phi_EE_caloUnmatched_new->Fill(deepSuperCluster_phi_tmp.at(iPF));
             h_PhiWidth_EE_caloUnmatched_new->Fill(deepSuperCluster_phiWidth_tmp.at(iPF)); 
             h_Energy_EE_caloUnmatched_new->Fill(deepSuperCluster_energy_tmp.at(iPF)); 
             h_nPFClusters_EE_caloUnmatched_new->Fill(deepSuperCluster_nPFClusters_tmp.at(iPF));  
             h_R9_EE_caloUnmatched_new->Fill(deepSuperCluster_R9_tmp.at(iPF));             
             h_full5x5_R9_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_R9_tmp.at(iPF));
             h_sigmaIetaIeta_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIeta_tmp.at(iPF));    
             h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIeta_tmp.at(iPF));       
             h_sigmaIetaIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIetaIphi_tmp.at(iPF));    
             h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIetaIphi_tmp.at(iPF));    
             h_sigmaIphiIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_sigmaIphiIphi_tmp.at(iPF));    
             h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new->Fill(deepSuperCluster_full5x5_sigmaIphiIphi_tmp.at(iPF));  
          }           
      } 

      entry++;
      simEnergy.clear();
      SuperCluster_RecoEnergy_simScore.clear(); 
      DeepSuperCluster_RecoEnergy_simScore.clear();   
   }

   TEfficiency* eff_SuperCluster_vs_EtaCalo = new TEfficiency(*h_Eta_Calo_SuperCluster_Num,*h_Eta_Calo_Denum);
   TEfficiency* eff_DeepSuperCluster_vs_EtaCalo = new TEfficiency(*h_Eta_Calo_DeepSuperCluster_Num,*h_Eta_Calo_Denum);
   drawEfficiency(eff_SuperCluster_vs_EtaCalo, eff_DeepSuperCluster_vs_EtaCalo, std::string("caloParticle_#eta"), std::string("Efficiency_vs_CaloEta")); 
  
   TEfficiency* eff_SuperCluster_vs_EtCalo = new TEfficiency(*h_Et_Calo_SuperCluster_Num,*h_Et_Calo_Denum);
   TEfficiency* eff_DeepSuperCluster_vs_EtCalo = new TEfficiency(*h_Et_Calo_DeepSuperCluster_Num,*h_Et_Calo_Denum);
   drawEfficiency(eff_SuperCluster_vs_EtCalo, eff_DeepSuperCluster_vs_EtCalo, std::string("caloParticle_Et (GeV)"), std::string("Efficiency_vs_CaloEt")); 

   drawProfile(prof_EoEtrue_vs_Eta_Calo_old, prof_EoEtrue_vs_Eta_Calo_new, std::string("caloParticle_#eta"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEta"));
   drawProfile(prof_EoEtrue_vs_Et_Calo_old, prof_EoEtrue_vs_Et_Calo_new, std::string("caloParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEt"));
   drawProfile(prof_EoEtrue_vs_nVtx_old, prof_EoEtrue_vs_nVtx_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_nVtx"));
   drawProfile(prof_EoEtrue_vs_Rho_old, prof_EoEtrue_vs_Rho_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_Rho"));
   
   /*double gaussEdgeL = 0.97;
   double gaussEdgeR = 1.03;
   double gaussNorm = h_EoEtrue_EB_old->GetEntries();
   double gaussMean = h_EoEtrue_EB_old->GetMean();
   double gaussSigma = h_EoEtrue_EB_old->GetRMS();

   TF1* cb1 = new TF1("cb1",&my2sideCrystalBall,gaussEdgeL,gaussEdgeR,7);
   cb1->SetParNames("alphaL","nL","mu","sigma","N","alphaR","nR");  
   cb1->SetParLimits(cb1->GetParNumber("nL"),0.1,15);
   cb1->SetParLimits(cb1->GetParNumber("nR"),0.1,15);
   cb1->SetParLimits(cb1->GetParNumber("alphaL"),-10,-0.01); 
   cb1->SetParLimits(cb1->GetParNumber("alphaR"),0.01,10);
   cb1->SetParameters((gaussEdgeL-gaussMean)/gaussSigma,5,gaussMean,gaussSigma,gaussNorm,(gaussEdgeR-gaussMean)/gaussSigma,5);
   cb1->SetLineColor(kRed);

   TFitResultPtr frp1 = h_EoEtrue_EB_old->Fit(cb1,"E L I S Q B R","HE",gaussEdgeL,gaussEdgeR);*/

   drawHisto(h_EoEtrue_EB_old, h_EoEtrue_EB_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("hist"), std::string("h_SC_EoEtrue_EB"), 0); 
   drawHisto(h_EoEtrue_EB_old, h_EoEtrue_EB_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("hist"), std::string("h_SC_EoEtrue_EB"), 1); 
   drawHisto(h_EoEtrue_EE_old, h_EoEtrue_EE_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("hist"), std::string("h_SC_EoEtrue_EE"), 0);
   drawHisto(h_EoEtrue_EE_old, h_EoEtrue_EE_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("hist"), std::string("h_SC_EoEtrue_EE"), 1);
   drawHisto(h_Energy_EB_old, h_Energy_EB_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB"), 0); 
   drawHisto(h_Energy_EB_old, h_Energy_EB_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB"), 1); 
   drawHisto(h_Energy_EE_old, h_Energy_EE_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE"), 0); 
   drawHisto(h_Energy_EE_old, h_Energy_EE_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE"), 1); 
   drawHisto(h_Eta_old, h_Eta_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta"), 0); 
   drawHisto(h_Eta_old, h_Eta_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta"), 1); 
   drawHisto(h_Phi_EB_old, h_Phi_EB_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB"), 0);
   drawHisto(h_Phi_EB_old, h_Phi_EB_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB"), 1); 
   drawHisto(h_Phi_EE_old, h_Phi_EE_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE"), 0);  
   drawHisto(h_Phi_EE_old, h_Phi_EE_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE"), 1);  
   drawHisto(h_EtaWidth_EB_old, h_EtaWidth_EB_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB"), 0);
   drawHisto(h_EtaWidth_EB_old, h_EtaWidth_EB_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB"), 1); 
   drawHisto(h_EtaWidth_EE_old, h_EtaWidth_EE_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE"), 0);  
   drawHisto(h_EtaWidth_EE_old, h_EtaWidth_EE_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE"), 1);  
   drawHisto(h_PhiWidth_EB_old, h_PhiWidth_EB_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB"), 0);
   drawHisto(h_PhiWidth_EB_old, h_PhiWidth_EB_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB"), 1); 
   drawHisto(h_PhiWidth_EE_old, h_PhiWidth_EE_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE"), 0);  
   drawHisto(h_PhiWidth_EE_old, h_PhiWidth_EE_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE"), 1);  
   drawHisto(h_nPFClusters_old, h_nPFClusters_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters"), 0); 
   drawHisto(h_nPFClusters_old, h_nPFClusters_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters"), 1); 
   drawHisto(h_nPFClusters_EB_old, h_nPFClusters_EB_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB"), 0); 
   drawHisto(h_nPFClusters_EB_old, h_nPFClusters_EB_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB"), 1); 
   drawHisto(h_nPFClusters_EE_old, h_nPFClusters_EE_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE"), 0); 
   drawHisto(h_nPFClusters_EE_old, h_nPFClusters_EE_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE"), 1); 
   drawHisto(h_R9_EB_old, h_R9_EB_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB"), 0); 
   drawHisto(h_R9_EB_old, h_R9_EB_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB"), 1); 
   drawHisto(h_R9_EE_old, h_R9_EE_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE"), 0); 
   drawHisto(h_R9_EE_old, h_R9_EE_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE"), 1); 
   drawHisto(h_full5x5_R9_EB_old, h_full5x5_R9_EB_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB"), 0); 
   drawHisto(h_full5x5_R9_EB_old, h_full5x5_R9_EB_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB"), 1); 
   drawHisto(h_full5x5_R9_EE_old, h_full5x5_R9_EE_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE"), 0); 
   drawHisto(h_full5x5_R9_EE_old, h_full5x5_R9_EE_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE"), 1); 
   drawHisto(h_sigmaIetaIeta_EB_old, h_sigmaIetaIeta_EB_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB"), 0); 
   drawHisto(h_sigmaIetaIeta_EB_old, h_sigmaIetaIeta_EB_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB"), 1); 
   drawHisto(h_sigmaIetaIeta_EE_old, h_sigmaIetaIeta_EE_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE"), 0); 
   drawHisto(h_sigmaIetaIeta_EE_old, h_sigmaIetaIeta_EE_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_old, h_full5x5_sigmaIetaIeta_EB_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_old, h_full5x5_sigmaIetaIeta_EB_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_old, h_full5x5_sigmaIetaIeta_EE_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_old, h_full5x5_sigmaIetaIeta_EE_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE"), 1); 
   drawHisto(h_sigmaIetaIphi_EB_old, h_sigmaIetaIphi_EB_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB"), 0); 
   drawHisto(h_sigmaIetaIphi_EB_old, h_sigmaIetaIphi_EB_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB"), 1); 
   drawHisto(h_sigmaIetaIphi_EE_old, h_sigmaIetaIphi_EE_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE"), 0); 
   drawHisto(h_sigmaIetaIphi_EE_old, h_sigmaIetaIphi_EE_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_old, h_full5x5_sigmaIetaIphi_EB_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_old, h_full5x5_sigmaIetaIphi_EB_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_old, h_full5x5_sigmaIetaIphi_EE_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_old, h_full5x5_sigmaIetaIphi_EE_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE"), 1); 
   drawHisto(h_sigmaIphiIphi_EB_old, h_sigmaIphiIphi_EB_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB"), 0); 
   drawHisto(h_sigmaIphiIphi_EB_old, h_sigmaIphiIphi_EB_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB"), 1); 
   drawHisto(h_sigmaIphiIphi_EE_old, h_sigmaIphiIphi_EE_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE"), 0); 
   drawHisto(h_sigmaIphiIphi_EE_old, h_sigmaIphiIphi_EE_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_old, h_full5x5_sigmaIphiIphi_EB_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_old, h_full5x5_sigmaIphiIphi_EB_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_old, h_full5x5_sigmaIphiIphi_EE_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_old, h_full5x5_sigmaIphiIphi_EE_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE"), 1); 

   //make seedMatched plots
   drawHisto(h_Energy_EB_seedMatched_old, h_Energy_EB_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_seedMatched"), 0); 
   drawHisto(h_Energy_EB_seedMatched_old, h_Energy_EB_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_seedMatched"), 1); 
   drawHisto(h_Energy_EE_seedMatched_old, h_Energy_EE_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_seedMatched"), 0); 
   drawHisto(h_Energy_EE_seedMatched_old, h_Energy_EE_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_seedMatched"), 1); 
   drawHisto(h_Eta_seedMatched_old, h_Eta_seedMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_seedMatched"), 0); 
   drawHisto(h_Eta_seedMatched_old, h_Eta_seedMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_seedMatched"), 1); 
   drawHisto(h_Phi_EB_seedMatched_old, h_Phi_EB_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_seedMatched"), 0);
   drawHisto(h_Phi_EB_seedMatched_old, h_Phi_EB_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_seedMatched"), 1); 
   drawHisto(h_Phi_EE_seedMatched_old, h_Phi_EE_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_seedMatched"), 0);  
   drawHisto(h_Phi_EE_seedMatched_old, h_Phi_EE_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_seedMatched"), 1);  
   drawHisto(h_EtaWidth_EB_seedMatched_old, h_EtaWidth_EB_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_seedMatched"), 0);
   drawHisto(h_EtaWidth_EB_seedMatched_old, h_EtaWidth_EB_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_seedMatched"), 1); 
   drawHisto(h_EtaWidth_EE_seedMatched_old, h_EtaWidth_EE_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_seedMatched"), 0);  
   drawHisto(h_EtaWidth_EE_seedMatched_old, h_EtaWidth_EE_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_seedMatched"), 1);  
   drawHisto(h_PhiWidth_EB_seedMatched_old, h_PhiWidth_EB_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_seedMatched"), 0);
   drawHisto(h_PhiWidth_EB_seedMatched_old, h_PhiWidth_EB_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_seedMatched"), 1); 
   drawHisto(h_PhiWidth_EE_seedMatched_old, h_PhiWidth_EE_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_seedMatched"), 0);  
   drawHisto(h_PhiWidth_EE_seedMatched_old, h_PhiWidth_EE_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_seedMatched"), 1);  
   drawHisto(h_nPFClusters_seedMatched_old, h_nPFClusters_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_seedMatched"), 0); 
   drawHisto(h_nPFClusters_seedMatched_old, h_nPFClusters_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_seedMatched"), 1); 
   drawHisto(h_nPFClusters_EB_seedMatched_old, h_nPFClusters_EB_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_seedMatched"), 0); 
   drawHisto(h_nPFClusters_EB_seedMatched_old, h_nPFClusters_EB_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_seedMatched"), 1); 
   drawHisto(h_nPFClusters_EE_seedMatched_old, h_nPFClusters_EE_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_seedMatched"), 0); 
   drawHisto(h_nPFClusters_EE_seedMatched_old, h_nPFClusters_EE_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_seedMatched"), 1); 
   drawHisto(h_R9_EB_seedMatched_old, h_R9_EB_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_seedMatched"), 0); 
   drawHisto(h_R9_EB_seedMatched_old, h_R9_EB_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_seedMatched"), 1); 
   drawHisto(h_R9_EE_seedMatched_old, h_R9_EE_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_seedMatched"), 0); 
   drawHisto(h_R9_EE_seedMatched_old, h_R9_EE_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_seedMatched"), 1); 
   drawHisto(h_full5x5_R9_EB_seedMatched_old, h_full5x5_R9_EB_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_seedMatched"), 0); 
   drawHisto(h_full5x5_R9_EB_seedMatched_old, h_full5x5_R9_EB_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_seedMatched"), 1); 
   drawHisto(h_full5x5_R9_EE_seedMatched_old, h_full5x5_R9_EE_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_seedMatched"), 0); 
   drawHisto(h_full5x5_R9_EE_seedMatched_old, h_full5x5_R9_EE_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_seedMatched"), 1); 
   drawHisto(h_sigmaIetaIeta_EB_seedMatched_old, h_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_seedMatched"), 0); 
   drawHisto(h_sigmaIetaIeta_EB_seedMatched_old, h_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_seedMatched"), 1); 
   drawHisto(h_sigmaIetaIeta_EE_seedMatched_old, h_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_seedMatched"), 0); 
   drawHisto(h_sigmaIetaIeta_EE_seedMatched_old, h_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_seedMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_seedMatched_old, h_full5x5_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_seedMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_seedMatched_old, h_full5x5_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_seedMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_seedMatched_old, h_full5x5_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_seedMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_seedMatched_old, h_full5x5_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_seedMatched"), 1); 
   drawHisto(h_sigmaIetaIphi_EB_seedMatched_old, h_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_seedMatched"), 0); 
   drawHisto(h_sigmaIetaIphi_EB_seedMatched_old, h_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_seedMatched"), 1); 
   drawHisto(h_sigmaIetaIphi_EE_seedMatched_old, h_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_seedMatched"), 0); 
   drawHisto(h_sigmaIetaIphi_EE_seedMatched_old, h_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_seedMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_seedMatched_old, h_full5x5_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_seedMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_seedMatched_old, h_full5x5_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_seedMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_seedMatched_old, h_full5x5_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_seedMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_seedMatched_old, h_full5x5_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_seedMatched"), 1); 
   drawHisto(h_sigmaIphiIphi_EB_seedMatched_old, h_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_seedMatched"), 0); 
   drawHisto(h_sigmaIphiIphi_EB_seedMatched_old, h_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_seedMatched"), 1); 
   drawHisto(h_sigmaIphiIphi_EE_seedMatched_old, h_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_seedMatched"), 0); 
   drawHisto(h_sigmaIphiIphi_EE_seedMatched_old, h_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_seedMatched"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_seedMatched_old, h_full5x5_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_seedMatched"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_seedMatched_old, h_full5x5_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_seedMatched"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_seedMatched_old, h_full5x5_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_seedMatched"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_seedMatched_old, h_full5x5_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_seedMatched"), 1);

   //make caloMatched plots
   drawHisto(h_Energy_EB_caloMatched_old, h_Energy_EB_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloMatched"), 0); 
   drawHisto(h_Energy_EB_caloMatched_old, h_Energy_EB_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloMatched"), 1); 
   drawHisto(h_Energy_EE_caloMatched_old, h_Energy_EE_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloMatched"), 0); 
   drawHisto(h_Energy_EE_caloMatched_old, h_Energy_EE_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloMatched"), 1); 
   drawHisto(h_Eta_caloMatched_old, h_Eta_caloMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloMatched"), 0); 
   drawHisto(h_Eta_caloMatched_old, h_Eta_caloMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloMatched"), 1); 
   drawHisto(h_Phi_EB_caloMatched_old, h_Phi_EB_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloMatched"), 0);
   drawHisto(h_Phi_EB_caloMatched_old, h_Phi_EB_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloMatched"), 1); 
   drawHisto(h_Phi_EE_caloMatched_old, h_Phi_EE_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloMatched"), 0);  
   drawHisto(h_Phi_EE_caloMatched_old, h_Phi_EE_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloMatched"), 1);  
   drawHisto(h_EtaWidth_EB_caloMatched_old, h_EtaWidth_EB_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloMatched"), 0);
   drawHisto(h_EtaWidth_EB_caloMatched_old, h_EtaWidth_EB_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloMatched"), 1); 
   drawHisto(h_EtaWidth_EE_caloMatched_old, h_EtaWidth_EE_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloMatched"), 0);  
   drawHisto(h_EtaWidth_EE_caloMatched_old, h_EtaWidth_EE_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloMatched"), 1);  
   drawHisto(h_PhiWidth_EB_caloMatched_old, h_PhiWidth_EB_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloMatched"), 0);
   drawHisto(h_PhiWidth_EB_caloMatched_old, h_PhiWidth_EB_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloMatched"), 1); 
   drawHisto(h_PhiWidth_EE_caloMatched_old, h_PhiWidth_EE_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloMatched"), 0);  
   drawHisto(h_PhiWidth_EE_caloMatched_old, h_PhiWidth_EE_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloMatched"), 1);  
   drawHisto(h_nPFClusters_caloMatched_old, h_nPFClusters_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloMatched"), 0); 
   drawHisto(h_nPFClusters_caloMatched_old, h_nPFClusters_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloMatched"), 1); 
   drawHisto(h_nPFClusters_EB_caloMatched_old, h_nPFClusters_EB_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloMatched"), 0); 
   drawHisto(h_nPFClusters_EB_caloMatched_old, h_nPFClusters_EB_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloMatched"), 1); 
   drawHisto(h_nPFClusters_EE_caloMatched_old, h_nPFClusters_EE_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloMatched"), 0); 
   drawHisto(h_nPFClusters_EE_caloMatched_old, h_nPFClusters_EE_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloMatched"), 1); 
   drawHisto(h_R9_EB_caloMatched_old, h_R9_EB_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloMatched"), 0); 
   drawHisto(h_R9_EB_caloMatched_old, h_R9_EB_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloMatched"), 1); 
   drawHisto(h_R9_EE_caloMatched_old, h_R9_EE_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloMatched"), 0); 
   drawHisto(h_R9_EE_caloMatched_old, h_R9_EE_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloMatched"), 1); 
   drawHisto(h_full5x5_R9_EB_caloMatched_old, h_full5x5_R9_EB_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloMatched"), 0); 
   drawHisto(h_full5x5_R9_EB_caloMatched_old, h_full5x5_R9_EB_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloMatched"), 1); 
   drawHisto(h_full5x5_R9_EE_caloMatched_old, h_full5x5_R9_EE_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloMatched"), 0); 
   drawHisto(h_full5x5_R9_EE_caloMatched_old, h_full5x5_R9_EE_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloMatched"), 1); 
   drawHisto(h_sigmaIetaIeta_EB_caloMatched_old, h_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloMatched"), 0); 
   drawHisto(h_sigmaIetaIeta_EB_caloMatched_old, h_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloMatched"), 1); 
   drawHisto(h_sigmaIetaIeta_EE_caloMatched_old, h_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloMatched"), 0); 
   drawHisto(h_sigmaIetaIeta_EE_caloMatched_old, h_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloMatched_old, h_full5x5_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloMatched_old, h_full5x5_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloMatched_old, h_full5x5_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloMatched_old, h_full5x5_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloMatched"), 1); 
   drawHisto(h_sigmaIetaIphi_EB_caloMatched_old, h_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloMatched"), 0); 
   drawHisto(h_sigmaIetaIphi_EB_caloMatched_old, h_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloMatched"), 1); 
   drawHisto(h_sigmaIetaIphi_EE_caloMatched_old, h_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloMatched"), 0); 
   drawHisto(h_sigmaIetaIphi_EE_caloMatched_old, h_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloMatched_old, h_full5x5_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloMatched_old, h_full5x5_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloMatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloMatched_old, h_full5x5_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloMatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloMatched_old, h_full5x5_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloMatched"), 1); 
   drawHisto(h_sigmaIphiIphi_EB_caloMatched_old, h_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloMatched"), 0); 
   drawHisto(h_sigmaIphiIphi_EB_caloMatched_old, h_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloMatched"), 1); 
   drawHisto(h_sigmaIphiIphi_EE_caloMatched_old, h_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloMatched"), 0); 
   drawHisto(h_sigmaIphiIphi_EE_caloMatched_old, h_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloMatched"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloMatched_old, h_full5x5_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloMatched"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloMatched_old, h_full5x5_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloMatched"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloMatched_old, h_full5x5_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloMatched"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloMatched_old, h_full5x5_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloMatched"), 1);

   //make caloUnmatched plots
   drawHisto(h_Energy_EB_caloUnmatched_old, h_Energy_EB_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloUnmatched"), 0); 
   drawHisto(h_Energy_EB_caloUnmatched_old, h_Energy_EB_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloUnmatched"), 1); 
   drawHisto(h_Energy_EE_caloUnmatched_old, h_Energy_EE_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloUnmatched"), 0); 
   drawHisto(h_Energy_EE_caloUnmatched_old, h_Energy_EE_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloUnmatched"), 1); 
   drawHisto(h_Eta_caloUnmatched_old, h_Eta_caloUnmatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloUnmatched"), 0); 
   drawHisto(h_Eta_caloUnmatched_old, h_Eta_caloUnmatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloUnmatched"), 1); 
   drawHisto(h_Phi_EB_caloUnmatched_old, h_Phi_EB_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloUnmatched"), 0);
   drawHisto(h_Phi_EB_caloUnmatched_old, h_Phi_EB_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloUnmatched"), 1); 
   drawHisto(h_Phi_EE_caloUnmatched_old, h_Phi_EE_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloUnmatched"), 0);  
   drawHisto(h_Phi_EE_caloUnmatched_old, h_Phi_EE_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloUnmatched"), 1);  
   drawHisto(h_EtaWidth_EB_caloUnmatched_old, h_EtaWidth_EB_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloUnmatched"), 0);
   drawHisto(h_EtaWidth_EB_caloUnmatched_old, h_EtaWidth_EB_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloUnmatched"), 1); 
   drawHisto(h_EtaWidth_EE_caloUnmatched_old, h_EtaWidth_EE_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloUnmatched"), 0);  
   drawHisto(h_EtaWidth_EE_caloUnmatched_old, h_EtaWidth_EE_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloUnmatched"), 1);  
   drawHisto(h_PhiWidth_EB_caloUnmatched_old, h_PhiWidth_EB_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloUnmatched"), 0);
   drawHisto(h_PhiWidth_EB_caloUnmatched_old, h_PhiWidth_EB_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloUnmatched"), 1); 
   drawHisto(h_PhiWidth_EE_caloUnmatched_old, h_PhiWidth_EE_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloUnmatched"), 0);  
   drawHisto(h_PhiWidth_EE_caloUnmatched_old, h_PhiWidth_EE_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloUnmatched"), 1);  
   drawHisto(h_nPFClusters_caloUnmatched_old, h_nPFClusters_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloUnmatched"), 0); 
   drawHisto(h_nPFClusters_caloUnmatched_old, h_nPFClusters_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloUnmatched"), 1); 
   drawHisto(h_nPFClusters_EB_caloUnmatched_old, h_nPFClusters_EB_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloUnmatched"), 0); 
   drawHisto(h_nPFClusters_EB_caloUnmatched_old, h_nPFClusters_EB_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloUnmatched"), 1); 
   drawHisto(h_nPFClusters_EE_caloUnmatched_old, h_nPFClusters_EE_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloUnmatched"), 0); 
   drawHisto(h_nPFClusters_EE_caloUnmatched_old, h_nPFClusters_EE_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloUnmatched"), 1); 
   drawHisto(h_R9_EB_caloUnmatched_old, h_R9_EB_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloUnmatched"), 0); 
   drawHisto(h_R9_EB_caloUnmatched_old, h_R9_EB_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloUnmatched"), 1); 
   drawHisto(h_R9_EE_caloUnmatched_old, h_R9_EE_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloUnmatched"), 0); 
   drawHisto(h_R9_EE_caloUnmatched_old, h_R9_EE_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloUnmatched"), 1); 
   drawHisto(h_full5x5_R9_EB_caloUnmatched_old, h_full5x5_R9_EB_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloUnmatched"), 0); 
   drawHisto(h_full5x5_R9_EB_caloUnmatched_old, h_full5x5_R9_EB_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloUnmatched"), 1); 
   drawHisto(h_full5x5_R9_EE_caloUnmatched_old, h_full5x5_R9_EE_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloUnmatched"), 0); 
   drawHisto(h_full5x5_R9_EE_caloUnmatched_old, h_full5x5_R9_EE_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloUnmatched"), 1); 
   drawHisto(h_sigmaIetaIeta_EB_caloUnmatched_old, h_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloUnmatched"), 0); 
   drawHisto(h_sigmaIetaIeta_EB_caloUnmatched_old, h_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloUnmatched"), 1); 
   drawHisto(h_sigmaIetaIeta_EE_caloUnmatched_old, h_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloUnmatched"), 0); 
   drawHisto(h_sigmaIetaIeta_EE_caloUnmatched_old, h_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloUnmatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloUnmatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloUnmatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloUnmatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloUnmatched"), 1); 
   drawHisto(h_sigmaIetaIphi_EB_caloUnmatched_old, h_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloUnmatched"), 0); 
   drawHisto(h_sigmaIetaIphi_EB_caloUnmatched_old, h_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloUnmatched"), 1); 
   drawHisto(h_sigmaIetaIphi_EE_caloUnmatched_old, h_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloUnmatched"), 0); 
   drawHisto(h_sigmaIetaIphi_EE_caloUnmatched_old, h_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloUnmatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloUnmatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloUnmatched"), 1); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloUnmatched"), 0); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloUnmatched"), 1); 
   drawHisto(h_sigmaIphiIphi_EB_caloUnmatched_old, h_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloUnmatched"), 0); 
   drawHisto(h_sigmaIphiIphi_EB_caloUnmatched_old, h_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloUnmatched"), 1); 
   drawHisto(h_sigmaIphiIphi_EE_caloUnmatched_old, h_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloUnmatched"), 0); 
   drawHisto(h_sigmaIphiIphi_EE_caloUnmatched_old, h_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloUnmatched"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloUnmatched"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloUnmatched"), 1); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloUnmatched"), 0); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloUnmatched"), 1); 
}

double my2sideCrystalBall(double* x, double* par) {

  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 

  Double_t xcur = x[0];
  Double_t alphaL = par[0];
  Double_t nL = par[1];
  Double_t mu = par[2];
  Double_t sigma = par[3];
  Double_t N = par[4];
  Double_t alphaR = par[5];
  Double_t nR = par[6];
  Double_t t = (xcur-mu)/sigma;
  Double_t absAlphaL = fabs((Double_t)alphaL);
  Double_t invAbsAlphaL = 1./absAlphaL;
  Double_t absAlphaR = fabs((Double_t)alphaR);
  Double_t invAbsAlphaR = 1./absAlphaR;

  
  if ( t<-absAlphaL ) {
    //cout<<"checkpoint dscb left"<<endl;
    Double_t AL = TMath::Power(nL*invAbsAlphaL,nL)*exp(-0.5*absAlphaL*absAlphaL);
    Double_t BL = nL*invAbsAlphaL - absAlphaL;
    return N*AL*TMath::Power(BL-t,-nL);
  } else if ( t <= absAlphaR )  {
    //cout<<"checkpoint dscb gaussian"<<endl;
    return N*exp(-0.5*t*t);
  } else {
    //cout<<"checkpoint dscb right"<<endl;
    Double_t AR = TMath::Power(nR*invAbsAlphaR,nR)*exp(-0.5*absAlphaR*absAlphaR);
    Double_t BR = nR*invAbsAlphaR - absAlphaR;
    return N*AR*TMath::Power(BR+t,-nR);
  }

}

void drawHisto(TH1F* h_old, TH1F* h_new, std::string x_label, std::string drawType, std::string Name, bool log)
{

   h_old->Scale(1./h_old->GetEntries());
   h_new->Scale(1./h_new->GetEntries());
   //h_old->SetTitle("");
   //h_new->SetTitle("");

   h_old->SetLineColor(kRed+1);
   h_old->SetLineWidth(2);
   h_old->GetXaxis()->SetTitle(x_label.c_str());

   h_new->SetLineColor(kBlue+1);
   h_new->SetLineWidth(2);

   float maximum = -1.;
   if( h_old -> GetMaximum() > h_new -> GetMaximum()) maximum=h_old -> GetMaximum();
   else maximum=h_new -> GetMaximum();  
   h_old -> SetMaximum( 1.1*maximum );
    
   TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(h_old,"Mustache","L");
   legend -> AddEntry(h_new,"DeepSC","L");

   TPaveStats* st_old = new TPaveStats();
   TPaveStats* st_new = new TPaveStats();
   TPaveStats* st_ratio = new TPaveStats();
   
   TCanvas* c = new TCanvas();
 
   TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
   TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
   
   cUp->SetBottomMargin(0.01); 
   cDown->SetTopMargin(0.01); 
   cDown->SetBottomMargin(0.2); 
    
   cUp->Draw();
   if(log) cUp->SetLogy();
   cDown->Draw();
     
   cUp->cd();
   h_old->Draw(drawType.c_str());
   gPad -> Update();
   st_old= (TPaveStats*)(h_old->GetListOfFunctions()->FindObject("stats"));
   st_old->SetX1NDC(0.82); //new x start position
   st_old->SetX2NDC(0.99); //new x end position
   st_old->SetY1NDC(0.82); //new y start position
   st_old->SetY2NDC(0.94); //new y end position
   st_old->SetTextColor(kRed+1);
   st_old->Draw("sames");
   h_new->Draw(std::string(drawType+",sames").c_str());
   gPad -> Update();
   st_new= (TPaveStats*)(h_new->GetListOfFunctions()->FindObject("stats"));
   st_new->SetX1NDC(0.82); //new x start position
   st_new->SetX2NDC(0.99); //new x end position
   st_new->SetY1NDC(0.68); //new y start position
   st_new->SetY2NDC(0.80); //new y end position
   st_new->SetTextColor(kBlue+1);
   st_new->Draw("sames");
   legend -> Draw("same");

   cDown->cd();
    
   TH1F* histo_ratio=(TH1F*)h_new->Clone("histo_ratio");
   histo_ratio->Sumw2();
   histo_ratio->Divide(h_old);
    
   histo_ratio -> GetXaxis() -> SetTitle(x_label.c_str());
   histo_ratio -> GetYaxis() -> SetTitle("DeespSC/Mustache");
   histo_ratio -> SetMaximum(2);
   histo_ratio -> SetMinimum(0);
   histo_ratio -> SetMarkerColor(kBlack);
   histo_ratio -> SetMarkerSize(0.5);
   histo_ratio ->SetTitle("");
   histo_ratio -> GetXaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetYaxis() -> SetLabelSize(0.07);
   histo_ratio -> GetXaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleSize(0.07);
   histo_ratio -> GetYaxis() -> SetTitleOffset(0.7);
   histo_ratio -> Draw("e");
   TF1* f_const = new TF1("f_1", "[0]",histo_ratio->GetBinCenter(1)-histo_ratio->GetBinWidth(1)/2, histo_ratio->GetBinCenter(histo_ratio->GetNbinsX())+histo_ratio->GetBinWidth(histo_ratio->GetNbinsX())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kRed);
   f_const -> SetLineWidth(2);
   f_const -> Draw("same");
    
   gPad -> Update();
   st_ratio= (TPaveStats*)(histo_ratio->GetListOfFunctions()->FindObject("stats"));
   st_ratio->SetX1NDC(0.); //new x start position
   st_ratio->SetX2NDC(0.); //new x end position
   st_ratio->SetY1NDC(0.); //new y start position
   st_ratio->SetY2NDC(0.); //new y end position
   st_ratio->SetTextColor(kBlack);
   st_ratio->Draw("sames");

   if(!log){
      c->SaveAs(std::string(Name+".png").c_str(),"png");
      c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	
   }else{
      c->SaveAs(std::string(Name+"_logY.png").c_str(),"png");
      c->SaveAs(std::string(Name+"_logY.pdf").c_str(),"pdf");
   }
}

void drawEfficiency(TEfficiency* eff_SuperCluster, TEfficiency* eff_DeepSuperCluster, std::string xtitle, std::string Name)
{
   eff_SuperCluster->SetLineColor(kRed+1);
   eff_SuperCluster->SetLineWidth(2);
   eff_SuperCluster->SetTitle(std::string(Name+"; "+xtitle+" ; Efficiency").c_str()); 
   
   eff_DeepSuperCluster->SetLineColor(kBlue+1);
   eff_DeepSuperCluster->SetLineWidth(2);

   TLegend* legend = new TLegend(0.365, 0.12, 0.635, 0.34);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(eff_SuperCluster,"Mustache","L");
   legend -> AddEntry(eff_DeepSuperCluster,"DeepSC","L");

   TCanvas* c = new TCanvas();
   eff_SuperCluster->Draw("APL");
   eff_DeepSuperCluster->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	
}

void drawProfile(TProfile* prof_SuperCluster, TProfile* prof_DeepSuperCluster, std::string xtitle, std::string ytitle, std::string Name)
{ 
   gStyle->SetOptStat(0000);  
   prof_SuperCluster->SetLineColor(kRed+1);
   prof_SuperCluster->SetLineWidth(2);
   prof_SuperCluster->GetXaxis()->SetTitle(xtitle.c_str()); 
   prof_SuperCluster->GetYaxis()->SetTitle(ytitle.c_str());  
   
   prof_DeepSuperCluster->SetLineColor(kBlue+1);
   prof_DeepSuperCluster->SetLineWidth(2);

   TLegend* legend = new TLegend(0.365, 0.12, 0.635, 0.34);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(prof_SuperCluster,"Mustache","L");
   legend -> AddEntry(prof_DeepSuperCluster,"DeepSC","L");

   TCanvas* c = new TCanvas();
   prof_SuperCluster->Draw("PL");
   prof_DeepSuperCluster->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
}

