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
#include "TStyle.h"
#include "TTreeReader.h"
#include <algorithm> 
#include <iostream>
#include <utility>
#include <string>

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

using namespace std;
using namespace edm;

//DEFINE HISTOGRAMS

//Total
TH1F* h_Energy_EB_old;
TH1F* h_Energy_EE_old;
TH1F* h_Energy_EB_new;
TH1F* h_Energy_EE_new;
TH1F* h_Et_EB_old;
TH1F* h_Et_EE_old;
TH1F* h_Et_EB_new;
TH1F* h_Et_EE_new;
TH1F* h_Eta_old;
TH1F* h_Eta_new;
TH1F* h_Phi_EB_old;
TH1F* h_Phi_EE_old;
TH1F* h_Phi_EB_new;
TH1F* h_Phi_EE_new;
TH1F* h_EtaWidth_EB_old;
TH1F* h_EtaWidth_EE_old;
TH1F* h_EtaWidth_EB_new;
TH1F* h_EtaWidth_EE_new;
TH1F* h_PhiWidth_EB_old;
TH1F* h_PhiWidth_EE_old;
TH1F* h_PhiWidth_EB_new;
TH1F* h_PhiWidth_EE_new;
TH1F* h_nPFClusters_old;
TH1F* h_nPFClusters_new;
TH1F* h_nPFClusters_EB_old;
TH1F* h_nPFClusters_EE_old;
TH1F* h_nPFClusters_EB_new;
TH1F* h_nPFClusters_EE_new;
TH1F* h_R9_EB_old;
TH1F* h_R9_EE_old;
TH1F* h_R9_EB_new;
TH1F* h_R9_EE_new;
TH1F* h_full5x5_R9_EB_old;
TH1F* h_full5x5_R9_EE_old;
TH1F* h_full5x5_R9_EB_new;
TH1F* h_full5x5_R9_EE_new;
TH1F* h_sigmaIetaIeta_EB_old;
TH1F* h_sigmaIetaIeta_EE_old;
TH1F* h_sigmaIetaIeta_EB_new;
TH1F* h_sigmaIetaIeta_EE_new;
TH1F* h_full5x5_sigmaIetaIeta_EB_old;
TH1F* h_full5x5_sigmaIetaIeta_EE_old;
TH1F* h_full5x5_sigmaIetaIeta_EB_new;
TH1F* h_full5x5_sigmaIetaIeta_EE_new;
TH1F* h_sigmaIetaIphi_EB_old;
TH1F* h_sigmaIetaIphi_EE_old;
TH1F* h_sigmaIetaIphi_EB_new;
TH1F* h_sigmaIetaIphi_EE_new;
TH1F* h_full5x5_sigmaIetaIphi_EB_old;
TH1F* h_full5x5_sigmaIetaIphi_EE_old;
TH1F* h_full5x5_sigmaIetaIphi_EB_new;
TH1F* h_full5x5_sigmaIetaIphi_EE_new;
TH1F* h_sigmaIphiIphi_EB_old;
TH1F* h_sigmaIphiIphi_EE_old;
TH1F* h_sigmaIphiIphi_EB_new;
TH1F* h_sigmaIphiIphi_EE_new;
TH1F* h_full5x5_sigmaIphiIphi_EB_old;
TH1F* h_full5x5_sigmaIphiIphi_EE_old;
TH1F* h_full5x5_sigmaIphiIphi_EB_new;
TH1F* h_full5x5_sigmaIphiIphi_EE_new;

//seed matched plots
TH1F* h_Energy_EB_seedMatched_old;
TH1F* h_Energy_EE_seedMatched_old;
TH1F* h_Energy_EB_seedMatched_new;
TH1F* h_Energy_EE_seedMatched_new;
TH1F* h_Et_EB_seedMatched_old;
TH1F* h_Et_EE_seedMatched_old;
TH1F* h_Et_EB_seedMatched_new;
TH1F* h_Et_EE_seedMatched_new;
TH1F* h_Eta_seedMatched_old;
TH1F* h_Eta_seedMatched_new;
TH1F* h_Phi_EB_seedMatched_old;
TH1F* h_Phi_EE_seedMatched_old;
TH1F* h_Phi_EB_seedMatched_new;
TH1F* h_Phi_EE_seedMatched_new;
TH1F* h_EtaWidth_EB_seedMatched_old;
TH1F* h_EtaWidth_EE_seedMatched_old;
TH1F* h_EtaWidth_EB_seedMatched_new;
TH1F* h_EtaWidth_EE_seedMatched_new;
TH1F* h_PhiWidth_EB_seedMatched_old;
TH1F* h_PhiWidth_EE_seedMatched_old;
TH1F* h_PhiWidth_EB_seedMatched_new;
TH1F* h_PhiWidth_EE_seedMatched_new;
TH1F* h_nPFClusters_seedMatched_old;
TH1F* h_nPFClusters_seedMatched_new;
TH1F* h_nPFClusters_EB_seedMatched_old;
TH1F* h_nPFClusters_EE_seedMatched_old;
TH1F* h_nPFClusters_EB_seedMatched_new;
TH1F* h_nPFClusters_EE_seedMatched_new;
TH1F* h_R9_EB_seedMatched_old;
TH1F* h_R9_EE_seedMatched_old;
TH1F* h_R9_EB_seedMatched_new;
TH1F* h_R9_EE_seedMatched_new;
TH1F* h_full5x5_R9_EB_seedMatched_old;
TH1F* h_full5x5_R9_EE_seedMatched_old;
TH1F* h_full5x5_R9_EB_seedMatched_new;
TH1F* h_full5x5_R9_EE_seedMatched_new;
TH1F* h_sigmaIetaIeta_EB_seedMatched_old;
TH1F* h_sigmaIetaIeta_EE_seedMatched_old;
TH1F* h_sigmaIetaIeta_EB_seedMatched_new;
TH1F* h_sigmaIetaIeta_EE_seedMatched_new;
TH1F* h_full5x5_sigmaIetaIeta_EB_seedMatched_old;
TH1F* h_full5x5_sigmaIetaIeta_EE_seedMatched_old;
TH1F* h_full5x5_sigmaIetaIeta_EB_seedMatched_new;
TH1F* h_full5x5_sigmaIetaIeta_EE_seedMatched_new;
TH1F* h_sigmaIetaIphi_EB_seedMatched_old;
TH1F* h_sigmaIetaIphi_EE_seedMatched_old;
TH1F* h_sigmaIetaIphi_EB_seedMatched_new;
TH1F* h_sigmaIetaIphi_EE_seedMatched_new;
TH1F* h_full5x5_sigmaIetaIphi_EB_seedMatched_old;
TH1F* h_full5x5_sigmaIetaIphi_EE_seedMatched_old;
TH1F* h_full5x5_sigmaIetaIphi_EB_seedMatched_new;
TH1F* h_full5x5_sigmaIetaIphi_EE_seedMatched_new;
TH1F* h_sigmaIphiIphi_EB_seedMatched_old;
TH1F* h_sigmaIphiIphi_EE_seedMatched_old;
TH1F* h_sigmaIphiIphi_EB_seedMatched_new;
TH1F* h_sigmaIphiIphi_EE_seedMatched_new;
TH1F* h_full5x5_sigmaIphiIphi_EB_seedMatched_old;
TH1F* h_full5x5_sigmaIphiIphi_EE_seedMatched_old;
TH1F* h_full5x5_sigmaIphiIphi_EB_seedMatched_new;
TH1F* h_full5x5_sigmaIphiIphi_EE_seedMatched_new;

//caloMatched
TH1F* h_Energy_EB_caloMatched_old;
TH1F* h_Energy_EE_caloMatched_old;
TH1F* h_Energy_EB_caloMatched_new;
TH1F* h_Energy_EE_caloMatched_new;
TH1F* h_EoEtrue_EB_old;
TH1F* h_EoEtrue_EE_old;
TH1F* h_EoEtrue_EB_new;
TH1F* h_EoEtrue_EE_new;
TH1F* h_EoEgen_EB_old;
TH1F* h_EoEgen_EE_old;
TH1F* h_EoEgen_EB_new;
TH1F* h_EoEgen_EE_new;
TH1F* h_Et_EB_caloMatched_old;
TH1F* h_Et_EE_caloMatched_old;
TH1F* h_Et_EB_caloMatched_new;
TH1F* h_Et_EE_caloMatched_new;
TH1F* h_Eta_caloMatched_old;
TH1F* h_Eta_caloMatched_new;
TH1F* h_Phi_EB_caloMatched_old;
TH1F* h_Phi_EE_caloMatched_old;
TH1F* h_Phi_EB_caloMatched_new;
TH1F* h_Phi_EE_caloMatched_new;
TH1F* h_EtaWidth_EB_caloMatched_old;
TH1F* h_EtaWidth_EE_caloMatched_old;
TH1F* h_EtaWidth_EB_caloMatched_new;
TH1F* h_EtaWidth_EE_caloMatched_new;
TH1F* h_PhiWidth_EB_caloMatched_old;
TH1F* h_PhiWidth_EE_caloMatched_old;
TH1F* h_PhiWidth_EB_caloMatched_new;
TH1F* h_PhiWidth_EE_caloMatched_new;
TH1F* h_nPFClusters_caloMatched_old;
TH1F* h_nPFClusters_caloMatched_new;
TH1F* h_nPFClusters_EB_caloMatched_old;
TH1F* h_nPFClusters_EE_caloMatched_old;
TH1F* h_nPFClusters_EB_caloMatched_new;
TH1F* h_nPFClusters_EE_caloMatched_new;
TH1F* h_R9_EB_caloMatched_old;
TH1F* h_R9_EE_caloMatched_old;
TH1F* h_R9_EB_caloMatched_new;
TH1F* h_R9_EE_caloMatched_new;
TH1F* h_full5x5_R9_EB_caloMatched_old;
TH1F* h_full5x5_R9_EE_caloMatched_old;
TH1F* h_full5x5_R9_EB_caloMatched_new;
TH1F* h_full5x5_R9_EE_caloMatched_new;
TH1F* h_sigmaIetaIeta_EB_caloMatched_old;
TH1F* h_sigmaIetaIeta_EE_caloMatched_old;
TH1F* h_sigmaIetaIeta_EB_caloMatched_new;
TH1F* h_sigmaIetaIeta_EE_caloMatched_new;
TH1F* h_full5x5_sigmaIetaIeta_EB_caloMatched_old;
TH1F* h_full5x5_sigmaIetaIeta_EE_caloMatched_old;
TH1F* h_full5x5_sigmaIetaIeta_EB_caloMatched_new;
TH1F* h_full5x5_sigmaIetaIeta_EE_caloMatched_new;
TH1F* h_sigmaIetaIphi_EB_caloMatched_old;
TH1F* h_sigmaIetaIphi_EE_caloMatched_old;
TH1F* h_sigmaIetaIphi_EB_caloMatched_new;
TH1F* h_sigmaIetaIphi_EE_caloMatched_new;
TH1F* h_full5x5_sigmaIetaIphi_EB_caloMatched_old;
TH1F* h_full5x5_sigmaIetaIphi_EE_caloMatched_old;
TH1F* h_full5x5_sigmaIetaIphi_EB_caloMatched_new;
TH1F* h_full5x5_sigmaIetaIphi_EE_caloMatched_new;
TH1F* h_sigmaIphiIphi_EB_caloMatched_old;
TH1F* h_sigmaIphiIphi_EE_caloMatched_old;
TH1F* h_sigmaIphiIphi_EB_caloMatched_new;
TH1F* h_sigmaIphiIphi_EE_caloMatched_new;
TH1F* h_full5x5_sigmaIphiIphi_EB_caloMatched_old;
TH1F* h_full5x5_sigmaIphiIphi_EE_caloMatched_old;
TH1F* h_full5x5_sigmaIphiIphi_EB_caloMatched_new;
TH1F* h_full5x5_sigmaIphiIphi_EE_caloMatched_new;

//caloMatched && seedMatched
TH1F* h_EoEtrue_EB_seedMatched_old;
TH1F* h_EoEtrue_EE_seedMatched_old;
TH1F* h_EoEtrue_EB_seedMatched_new;
TH1F* h_EoEtrue_EE_seedMatched_new;
TH1F* h_EoEgen_EB_seedMatched_old;
TH1F* h_EoEgen_EE_seedMatched_old;
TH1F* h_EoEgen_EB_seedMatched_new;
TH1F* h_EoEgen_EE_seedMatched_new;

//caloUnmatched
TH1F* h_Energy_EB_caloUnmatched_old;
TH1F* h_Energy_EE_caloUnmatched_old;
TH1F* h_Energy_EB_caloUnmatched_new;
TH1F* h_Energy_EE_caloUnmatched_new;
TH1F* h_Et_EB_caloUnmatched_old;
TH1F* h_Et_EE_caloUnmatched_old;
TH1F* h_Et_EB_caloUnmatched_new;
TH1F* h_Et_EE_caloUnmatched_new;
TH1F* h_Eta_caloUnmatched_old;
TH1F* h_Eta_caloUnmatched_new;
TH1F* h_Phi_EB_caloUnmatched_old;
TH1F* h_Phi_EE_caloUnmatched_old;
TH1F* h_Phi_EB_caloUnmatched_new;
TH1F* h_Phi_EE_caloUnmatched_new;
TH1F* h_EtaWidth_EB_caloUnmatched_old;
TH1F* h_EtaWidth_EE_caloUnmatched_old;
TH1F* h_EtaWidth_EB_caloUnmatched_new;
TH1F* h_EtaWidth_EE_caloUnmatched_new;
TH1F* h_PhiWidth_EB_caloUnmatched_old;
TH1F* h_PhiWidth_EE_caloUnmatched_old;
TH1F* h_PhiWidth_EB_caloUnmatched_new;
TH1F* h_PhiWidth_EE_caloUnmatched_new;
TH1F* h_nPFClusters_caloUnmatched_old;
TH1F* h_nPFClusters_caloUnmatched_new;
TH1F* h_nPFClusters_EB_caloUnmatched_old;
TH1F* h_nPFClusters_EE_caloUnmatched_old;
TH1F* h_nPFClusters_EB_caloUnmatched_new;
TH1F* h_nPFClusters_EE_caloUnmatched_new;
TH1F* h_R9_EB_caloUnmatched_old;
TH1F* h_R9_EE_caloUnmatched_old;
TH1F* h_R9_EB_caloUnmatched_new;
TH1F* h_R9_EE_caloUnmatched_new;
TH1F* h_full5x5_R9_EB_caloUnmatched_old;
TH1F* h_full5x5_R9_EE_caloUnmatched_old;
TH1F* h_full5x5_R9_EB_caloUnmatched_new;
TH1F* h_full5x5_R9_EE_caloUnmatched_new;
TH1F* h_sigmaIetaIeta_EB_caloUnmatched_old;
TH1F* h_sigmaIetaIeta_EE_caloUnmatched_old;
TH1F* h_sigmaIetaIeta_EB_caloUnmatched_new;
TH1F* h_sigmaIetaIeta_EE_caloUnmatched_new;
TH1F* h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old;
TH1F* h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old;
TH1F* h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new;
TH1F* h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new;
TH1F* h_sigmaIetaIphi_EB_caloUnmatched_old;
TH1F* h_sigmaIetaIphi_EE_caloUnmatched_old;
TH1F* h_sigmaIetaIphi_EB_caloUnmatched_new;
TH1F* h_sigmaIetaIphi_EE_caloUnmatched_new;
TH1F* h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old;
TH1F* h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old;
TH1F* h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new;
TH1F* h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new;
TH1F* h_sigmaIphiIphi_EB_caloUnmatched_old;
TH1F* h_sigmaIphiIphi_EE_caloUnmatched_old;
TH1F* h_sigmaIphiIphi_EB_caloUnmatched_new;
TH1F* h_sigmaIphiIphi_EE_caloUnmatched_new;
TH1F* h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old;
TH1F* h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old;
TH1F* h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new;
TH1F* h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new;

//efficiencies
TH1F* h_Eta_Calo_Denum;
TH1F* h_Eta_Calo_old;
TH1F* h_Eta_Calo_new;
TH1F* h_Eta_Gen_Denum;
TH1F* h_Eta_Gen_old;
TH1F* h_Eta_Gen_new;
TH1F* h_Eta_Seed_Denum;
TH1F* h_Eta_Seed_old;
TH1F* h_Eta_Seed_new;
TH1F* h_Et_Calo_EB_Denum;
TH1F* h_Et_Calo_EB_old;
TH1F* h_Et_Calo_EB_new;
TH1F* h_Et_Calo_EE_Denum;
TH1F* h_Et_Calo_EE_old;
TH1F* h_Et_Calo_EE_new;
TH1F* h_Et_Gen_EB_Denum;
TH1F* h_Et_Gen_EB_old;
TH1F* h_Et_Gen_EB_new;
TH1F* h_Et_Gen_EE_Denum;
TH1F* h_Et_Gen_EE_old;
TH1F* h_Et_Gen_EE_new;
TH1F* h_Et_Seed_EB_Denum;
TH1F* h_Et_Seed_EB_old;
TH1F* h_Et_Seed_EB_new;
TH1F* h_Et_Seed_EE_Denum;
TH1F* h_Et_Seed_EE_old;
TH1F* h_Et_Seed_EE_new;
TEfficiency* eff_SuperCluster_vs_EtaCalo;
TEfficiency* eff_DeepSuperCluster_vs_EtaCalo;    
TEfficiency* eff_SuperCluster_vs_EtaSeed;
TEfficiency* eff_DeepSuperCluster_vs_EtaSeed;    
TEfficiency* eff_SuperCluster_vs_EtCalo_EB;
TEfficiency* eff_DeepSuperCluster_vs_EtCalo_EB;
TEfficiency* eff_SuperCluster_vs_EtSeed_EB;
TEfficiency* eff_DeepSuperCluster_vs_EtSeed_EB;
TEfficiency* eff_SuperCluster_vs_EtCalo_EE;
TEfficiency* eff_DeepSuperCluster_vs_EtCalo_EE;
TEfficiency* eff_SuperCluster_vs_EtSeed_EE;
TEfficiency* eff_DeepSuperCluster_vs_EtSeed_EE;

//DEFINE PROFILES
TProfile* prof_EoEtrue_vs_Eta_Calo_old;
TProfile* prof_EoEtrue_vs_Eta_Calo_new;
TProfile* prof_EoEtrue_vs_Eta_Seed_old;
TProfile* prof_EoEtrue_vs_Eta_Seed_new;
TProfile* prof_EoEtrue_vs_Et_Calo_EB_old;
TProfile* prof_EoEtrue_vs_Et_Calo_EB_new;
TProfile* prof_EoEtrue_vs_Et_Calo_EE_old;
TProfile* prof_EoEtrue_vs_Et_Calo_EE_new;
TProfile* prof_EoEtrue_vs_Et_Seed_EB_old;
TProfile* prof_EoEtrue_vs_Et_Seed_EB_new;
TProfile* prof_EoEtrue_vs_Et_Seed_EE_old;
TProfile* prof_EoEtrue_vs_Et_Seed_EE_new;
TProfile* prof_EoEtrue_vs_Energy_Calo_EB_old;
TProfile* prof_EoEtrue_vs_Energy_Calo_EB_new;
TProfile* prof_EoEtrue_vs_Energy_Calo_EE_old;
TProfile* prof_EoEtrue_vs_Energy_Calo_EE_new;
TProfile* prof_EoEtrue_vs_nVtx_EB_old;
TProfile* prof_EoEtrue_vs_nVtx_EB_new;
TProfile* prof_EoEtrue_vs_nVtx_EE_old;
TProfile* prof_EoEtrue_vs_nVtx_EE_new;
TProfile* prof_EoEtrue_vs_Rho_EB_old;
TProfile* prof_EoEtrue_vs_Rho_EB_new;
TProfile* prof_EoEtrue_vs_Rho_EE_old;
TProfile* prof_EoEtrue_vs_Rho_EE_new;
TProfile* prof_EoEgen_vs_Eta_Gen_old;
TProfile* prof_EoEgen_vs_Eta_Gen_new;
TProfile* prof_EoEgen_vs_Eta_Seed_old;
TProfile* prof_EoEgen_vs_Eta_Seed_new;
TProfile* prof_EoEgen_vs_Et_Gen_EB_old;
TProfile* prof_EoEgen_vs_Et_Gen_EB_new;
TProfile* prof_EoEgen_vs_Et_Gen_EE_old;
TProfile* prof_EoEgen_vs_Et_Gen_EE_new;
TProfile* prof_EoEgen_vs_Et_Seed_EB_old;
TProfile* prof_EoEgen_vs_Et_Seed_EB_new;
TProfile* prof_EoEgen_vs_Et_Seed_EE_old;
TProfile* prof_EoEgen_vs_Et_Seed_EE_new;
TProfile* prof_EoEgen_vs_Energy_Gen_EB_old;
TProfile* prof_EoEgen_vs_Energy_Gen_EB_new;
TProfile* prof_EoEgen_vs_Energy_Gen_EE_old;
TProfile* prof_EoEgen_vs_Energy_Gen_EE_new;
TProfile* prof_EoEgen_vs_nVtx_EB_old;
TProfile* prof_EoEgen_vs_nVtx_EB_new;
TProfile* prof_EoEgen_vs_nVtx_EE_old;
TProfile* prof_EoEgen_vs_nVtx_EE_new;
TProfile* prof_EoEgen_vs_Rho_EB_old;
TProfile* prof_EoEgen_vs_Rho_EB_new;
TProfile* prof_EoEgen_vs_Rho_EE_old;
TProfile* prof_EoEgen_vs_Rho_EE_new;

//DEFINE HISTOGRAM VECTORS 
std::vector<TH1F*> EoEtrue_vs_Eta_Calo_old;
std::vector<TH1F*> EoEtrue_vs_Eta_Calo_new;
std::vector<TH1F*> EoEtrue_vs_Eta_Seed_old;
std::vector<TH1F*> EoEtrue_vs_Eta_Seed_new;
std::vector<TH1F*> EoEtrue_vs_Et_Calo_EB_old;
std::vector<TH1F*> EoEtrue_vs_Et_Calo_EB_new;
std::vector<TH1F*> EoEtrue_vs_Et_Seed_EB_old;
std::vector<TH1F*> EoEtrue_vs_Et_Seed_EB_new;
std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EB_old;
std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EB_new;
std::vector<TH1F*> EoEtrue_vs_nVtx_EB_old;
std::vector<TH1F*> EoEtrue_vs_nVtx_EB_new;
std::vector<TH1F*> EoEtrue_vs_Rho_EB_old;
std::vector<TH1F*> EoEtrue_vs_Rho_EB_new;
std::vector<TH1F*> EoEtrue_vs_Et_Calo_EE_old;
std::vector<TH1F*> EoEtrue_vs_Et_Calo_EE_new;
std::vector<TH1F*> EoEtrue_vs_Et_Seed_EE_old;
std::vector<TH1F*> EoEtrue_vs_Et_Seed_EE_new;
std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EE_old;
std::vector<TH1F*> EoEtrue_vs_Energy_Calo_EE_new;
std::vector<TH1F*> EoEtrue_vs_nVtx_EE_old;
std::vector<TH1F*> EoEtrue_vs_nVtx_EE_new;
std::vector<TH1F*> EoEtrue_vs_Rho_EE_old;
std::vector<TH1F*> EoEtrue_vs_Rho_EE_new;
std::vector<TH1F*> EoEgen_vs_Eta_Gen_old;
std::vector<TH1F*> EoEgen_vs_Eta_Gen_new;
std::vector<TH1F*> EoEgen_vs_Eta_Seed_old;
std::vector<TH1F*> EoEgen_vs_Eta_Seed_new;
std::vector<TH1F*> EoEgen_vs_Et_Gen_EB_old;
std::vector<TH1F*> EoEgen_vs_Et_Gen_EB_new;
std::vector<TH1F*> EoEgen_vs_Et_Seed_EB_old;
std::vector<TH1F*> EoEgen_vs_Et_Seed_EB_new;
std::vector<TH1F*> EoEgen_vs_Energy_Gen_EB_old;
std::vector<TH1F*> EoEgen_vs_Energy_Gen_EB_new;
std::vector<TH1F*> EoEgen_vs_nVtx_EB_old;
std::vector<TH1F*> EoEgen_vs_nVtx_EB_new;
std::vector<TH1F*> EoEgen_vs_Rho_EB_old;
std::vector<TH1F*> EoEgen_vs_Rho_EB_new;
std::vector<TH1F*> EoEgen_vs_Et_Gen_EE_old;
std::vector<TH1F*> EoEgen_vs_Et_Gen_EE_new;
std::vector<TH1F*> EoEgen_vs_Et_Seed_EE_old;
std::vector<TH1F*> EoEgen_vs_Et_Seed_EE_new;
std::vector<TH1F*> EoEgen_vs_Energy_Gen_EE_old;
std::vector<TH1F*> EoEgen_vs_Energy_Gen_EE_new;
std::vector<TH1F*> EoEgen_vs_nVtx_EE_old;
std::vector<TH1F*> EoEgen_vs_nVtx_EE_new;
std::vector<TH1F*> EoEgen_vs_Rho_EE_old;
std::vector<TH1F*> EoEgen_vs_Rho_EE_new;


//DEFINE BRANCHES
int nVtx;
float rho;
vector<float>   *genParticle_energy;
vector<float>   *genParticle_eta;
vector<float>   *genParticle_phi;
vector<vector<int> > *genParticle_superCluster_dR_genScore_MatchedIndex;
vector<vector<int> > *genParticle_deepSuperCluster_dR_genScore_MatchedIndex;
vector<vector<int>> *caloParticle_superCluster_simScore_MatchedIndex;
vector<vector<int>> *caloParticle_deepSuperCluster_simScore_MatchedIndex;
vector<float> *caloParticle_simEnergy;
vector<float> *caloParticle_simEta;
vector<float> *caloParticle_simEt;
vector<float> *caloParticle_simPhi;
vector<float> *caloParticle_genEnergy;
vector<float> *caloParticle_genEta;
vector<float> *caloParticle_genEt;
vector<float> *caloParticle_genPhi;
vector<float> *superCluster_energy;
vector<float> *superCluster_eta;
vector<float> *superCluster_phi;
vector<float> *superCluster_etaWidth;
vector<float> *superCluster_phiWidth;
vector<int> *superCluster_nPFClusters;
vector<int> *superCluster_seedIndex;
vector<float> *superCluster_swissCross;
vector<float> *superCluster_r9;
vector<float> *superCluster_sigmaIetaIeta;
vector<float> *superCluster_sigmaIetaIphi;
vector<float> *superCluster_sigmaIphiIphi;
vector<float> *superCluster_full5x5_r9;
vector<float> *superCluster_full5x5_sigmaIetaIeta;
vector<float> *superCluster_full5x5_sigmaIetaIphi;
vector<float> *superCluster_full5x5_sigmaIphiIphi;
vector<int> *superCluster_simScore_MatchedIndex;
vector<int> *superCluster_dR_genScore_MatchedIndex;
vector<vector<double>> *superCluster_simScore;
vector<vector<double>> *superCluster_dR_genScore;
vector<float> *deepSuperCluster_energy;
vector<float> *deepSuperCluster_eta;
vector<float> *deepSuperCluster_phi;
vector<float> *deepSuperCluster_etaWidth;
vector<float> *deepSuperCluster_phiWidth;
vector<int> *deepSuperCluster_nPFClusters;
vector<int> *deepSuperCluster_seedIndex;
vector<float> *deepSuperCluster_swissCross;
vector<float> *deepSuperCluster_r9;
vector<float> *deepSuperCluster_sigmaIetaIeta;
vector<float> *deepSuperCluster_sigmaIetaIphi;
vector<float> *deepSuperCluster_sigmaIphiIphi;
vector<float> *deepSuperCluster_full5x5_r9;
vector<float> *deepSuperCluster_full5x5_sigmaIetaIeta;
vector<float> *deepSuperCluster_full5x5_sigmaIetaIphi;
vector<float> *deepSuperCluster_full5x5_sigmaIphiIphi;
vector<int> *deepSuperCluster_simScore_MatchedIndex;
vector<int> *deepSuperCluster_dR_genScore_MatchedIndex;
vector<vector<double>> *deepSuperCluster_dR_genScore;
vector<vector<double>> *deepSuperCluster_simScore;
vector<float> *pfCluster_eta;
vector<float> *pfCluster_energy;
vector<float> *pfCluster_phi;

TBranch *b_nVtx;
TBranch *b_rho;
TBranch *b_genParticle_energy;
TBranch *b_genParticle_eta;
TBranch *b_genParticle_phi;
TBranch *b_genParticle_superCluster_dR_genScore_MatchedIndex;
TBranch *b_genParticle_deepSuperCluster_dR_genScore_MatchedIndex;
TBranch *b_caloParticle_superCluster_simScore_MatchedIndex;
TBranch *b_caloParticle_deepSuperCluster_simScore_MatchedIndex;
TBranch *b_caloParticle_simEnergy;
TBranch *b_caloParticle_simEta;
TBranch *b_caloParticle_simEt;
TBranch *b_caloParticle_simPhi;
TBranch *b_caloParticle_genEnergy;
TBranch *b_caloParticle_genEta;
TBranch *b_caloParticle_genEt;
TBranch *b_caloParticle_genPhi;
TBranch *b_superCluster_energy;
TBranch *b_superCluster_eta;
TBranch *b_superCluster_phi;
TBranch *b_superCluster_etaWidth;
TBranch *b_superCluster_phiWidth;
TBranch *b_superCluster_nPFClusters;
TBranch *b_superCluster_seedIndex;
TBranch *b_superCluster_swissCross;
TBranch *b_superCluster_r9;
TBranch *b_superCluster_sigmaIetaIeta;
TBranch *b_superCluster_sigmaIetaIphi;
TBranch *b_superCluster_sigmaIphiIphi;
TBranch *b_superCluster_full5x5_r9;
TBranch *b_superCluster_full5x5_sigmaIetaIeta;
TBranch *b_superCluster_full5x5_sigmaIetaIphi;
TBranch *b_superCluster_full5x5_sigmaIphiIphi;
TBranch *b_superCluster_simScore_MatchedIndex;
TBranch *b_superCluster_dR_genScore_MatchedIndex;
TBranch *b_superCluster_simScore;
TBranch *b_superCluster_dR_genScore;
TBranch *b_deepSuperCluster_energy;
TBranch *b_deepSuperCluster_eta;
TBranch *b_deepSuperCluster_phi;
TBranch *b_deepSuperCluster_etaWidth;
TBranch *b_deepSuperCluster_phiWidth;
TBranch *b_deepSuperCluster_nPFClusters;
TBranch *b_deepSuperCluster_seedIndex;
TBranch *b_deepSuperCluster_swissCross;
TBranch *b_deepSuperCluster_r9;
TBranch *b_deepSuperCluster_sigmaIetaIeta;
TBranch *b_deepSuperCluster_sigmaIetaIphi;
TBranch *b_deepSuperCluster_sigmaIphiIphi;
TBranch *b_deepSuperCluster_full5x5_r9;
TBranch *b_deepSuperCluster_full5x5_sigmaIetaIeta;
TBranch *b_deepSuperCluster_full5x5_sigmaIetaIphi;
TBranch *b_deepSuperCluster_full5x5_sigmaIphiIphi;
TBranch *b_deepSuperCluster_simScore_MatchedIndex;
TBranch *b_deepSuperCluster_dR_genScore_MatchedIndex;
TBranch *b_deepSuperCluster_dR_genScore;
TBranch *b_deepSuperCluster_simScore;
TBranch *b_pfCluster_eta;
TBranch *b_pfCluster_energy;
TBranch *b_pfCluster_phi;

//setTreeBranches
void setTreeBranches(TTree* tree, std::string superClusterRef, std::string superClusterVal)
{
   tree->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   tree->SetBranchAddress("rho", &rho, &b_rho);
   tree->SetBranchAddress("genParticle_energy", &genParticle_energy, &b_genParticle_energy);
   tree->SetBranchAddress("genParticle_eta", &genParticle_eta, &b_genParticle_eta);
   tree->SetBranchAddress("genParticle_phi", &genParticle_phi, &b_genParticle_phi);
   tree->SetBranchAddress(std::string("genParticle_"+superClusterRef+"_dR_genScore_MatchedIndex").c_str(), &genParticle_superCluster_dR_genScore_MatchedIndex, &b_genParticle_superCluster_dR_genScore_MatchedIndex);
   tree->SetBranchAddress(std::string("genParticle_"+superClusterVal+"_dR_genScore_MatchedIndex").c_str(), &genParticle_deepSuperCluster_dR_genScore_MatchedIndex, &b_genParticle_deepSuperCluster_dR_genScore_MatchedIndex);
   tree->SetBranchAddress("caloParticle_simEnergy", &caloParticle_simEnergy, &b_caloParticle_simEnergy);
   tree->SetBranchAddress("caloParticle_simPt", &caloParticle_simEt, &b_caloParticle_simEt);
   tree->SetBranchAddress("caloParticle_simEta", &caloParticle_simEta, &b_caloParticle_simEta);
   tree->SetBranchAddress("caloParticle_simPhi", &caloParticle_simPhi, &b_caloParticle_simPhi);
   tree->SetBranchAddress("caloParticle_genEnergy", &caloParticle_genEnergy, &b_caloParticle_genEnergy);
   tree->SetBranchAddress("caloParticle_genPt", &caloParticle_genEt, &b_caloParticle_genEt);
   tree->SetBranchAddress("caloParticle_genEta", &caloParticle_genEta, &b_caloParticle_genEta);
   tree->SetBranchAddress("caloParticle_genPhi", &caloParticle_genPhi, &b_caloParticle_genPhi);
   tree->SetBranchAddress(std::string("caloParticle_"+superClusterRef+"_simScore_MatchedIndex").c_str(), &caloParticle_superCluster_simScore_MatchedIndex, &b_caloParticle_superCluster_simScore_MatchedIndex);
   tree->SetBranchAddress(std::string("caloParticle_"+superClusterVal+"_simScore_MatchedIndex").c_str(), &caloParticle_deepSuperCluster_simScore_MatchedIndex, &b_caloParticle_deepSuperCluster_simScore_MatchedIndex);
   tree->SetBranchAddress(std::string(superClusterRef+"_energy").c_str(), &superCluster_energy, &b_superCluster_energy);
   tree->SetBranchAddress(std::string(superClusterRef+"_eta").c_str(), &superCluster_eta, &b_superCluster_eta);
   tree->SetBranchAddress(std::string(superClusterRef+"_phi").c_str(), &superCluster_phi, &b_superCluster_phi);
   tree->SetBranchAddress(std::string(superClusterRef+"_etaWidth").c_str(), &superCluster_etaWidth, &b_superCluster_etaWidth);
   tree->SetBranchAddress(std::string(superClusterRef+"_phiWidth").c_str(), &superCluster_phiWidth, &b_superCluster_phiWidth);
   tree->SetBranchAddress(std::string(superClusterRef+"_nPFClusters").c_str(), &superCluster_nPFClusters, &b_superCluster_nPFClusters);
   tree->SetBranchAddress(std::string(superClusterRef+"_seedIndex").c_str(), &superCluster_seedIndex, &b_superCluster_seedIndex);
   tree->SetBranchAddress(std::string(superClusterVal+"_energy").c_str(), &deepSuperCluster_energy, &b_deepSuperCluster_energy);
   tree->SetBranchAddress(std::string(superClusterVal+"_eta").c_str(), &deepSuperCluster_eta, &b_deepSuperCluster_eta);
   tree->SetBranchAddress(std::string(superClusterVal+"_phi").c_str(), &deepSuperCluster_phi, &b_deepSuperCluster_phi);
   tree->SetBranchAddress(std::string(superClusterVal+"_etaWidth").c_str(), &deepSuperCluster_etaWidth, &b_deepSuperCluster_etaWidth);
   tree->SetBranchAddress(std::string(superClusterVal+"_phiWidth").c_str(), &deepSuperCluster_phiWidth, &b_deepSuperCluster_phiWidth);
   tree->SetBranchAddress(std::string(superClusterVal+"_nPFClusters").c_str(), &deepSuperCluster_nPFClusters, &b_deepSuperCluster_nPFClusters);
   tree->SetBranchAddress(std::string(superClusterVal+"_seedIndex").c_str(), &deepSuperCluster_seedIndex, &b_deepSuperCluster_seedIndex);
   tree->SetBranchAddress(std::string(superClusterRef+"_dR_genScore_MatchedIndex").c_str(), &superCluster_dR_genScore_MatchedIndex, &b_superCluster_dR_genScore_MatchedIndex);
   tree->SetBranchAddress(std::string(superClusterRef+"_simScore_MatchedIndex").c_str(), &superCluster_simScore_MatchedIndex, &b_superCluster_simScore_MatchedIndex);
   tree->SetBranchAddress(std::string(superClusterVal+"_dR_genScore_MatchedIndex").c_str(), &deepSuperCluster_dR_genScore_MatchedIndex, &b_deepSuperCluster_dR_genScore_MatchedIndex);
   tree->SetBranchAddress(std::string(superClusterVal+"_simScore_MatchedIndex").c_str(), &deepSuperCluster_simScore_MatchedIndex, &b_deepSuperCluster_simScore_MatchedIndex);
   tree->SetBranchAddress(std::string(superClusterRef+"_dR_genScore").c_str(), &superCluster_dR_genScore, &b_superCluster_dR_genScore);
   tree->SetBranchAddress(std::string(superClusterRef+"_simScore").c_str(), &superCluster_simScore, &b_superCluster_simScore);
   tree->SetBranchAddress(std::string(superClusterVal+"_dR_genScore").c_str(), &deepSuperCluster_dR_genScore, &b_deepSuperCluster_dR_genScore);
   tree->SetBranchAddress(std::string(superClusterVal+"_simScore").c_str(), &deepSuperCluster_simScore, &b_deepSuperCluster_simScore);
   tree->SetBranchAddress(std::string(superClusterRef+"_swissCross").c_str(), &superCluster_swissCross, &b_superCluster_swissCross);
   tree->SetBranchAddress(std::string(superClusterRef+"_r9").c_str(), &superCluster_r9, &b_superCluster_r9);
   tree->SetBranchAddress(std::string(superClusterRef+"_sigmaIetaIeta").c_str(), &superCluster_sigmaIetaIeta, &b_superCluster_sigmaIetaIeta);
   tree->SetBranchAddress(std::string(superClusterRef+"_sigmaIetaIphi").c_str(), &superCluster_sigmaIetaIphi, &b_superCluster_sigmaIetaIphi);
   tree->SetBranchAddress(std::string(superClusterRef+"_sigmaIphiIphi").c_str(), &superCluster_sigmaIphiIphi, &b_superCluster_sigmaIphiIphi);
   tree->SetBranchAddress(std::string(superClusterRef+"_full5x5_r9").c_str(), &superCluster_full5x5_r9, &b_superCluster_full5x5_r9);
   tree->SetBranchAddress(std::string(superClusterRef+"_full5x5_sigmaIetaIeta").c_str(), &superCluster_full5x5_sigmaIetaIeta, &b_superCluster_full5x5_sigmaIetaIeta);
   tree->SetBranchAddress(std::string(superClusterRef+"_full5x5_sigmaIetaIphi").c_str(), &superCluster_full5x5_sigmaIetaIphi, &b_superCluster_full5x5_sigmaIetaIphi);
   tree->SetBranchAddress(std::string(superClusterRef+"_full5x5_sigmaIphiIphi").c_str(), &superCluster_full5x5_sigmaIphiIphi, &b_superCluster_full5x5_sigmaIphiIphi);
   tree->SetBranchAddress(std::string(superClusterVal+"_swissCross").c_str(), &deepSuperCluster_swissCross, &b_deepSuperCluster_swissCross);
   tree->SetBranchAddress(std::string(superClusterVal+"_r9").c_str(), &deepSuperCluster_r9, &b_deepSuperCluster_r9);
   tree->SetBranchAddress(std::string(superClusterVal+"_sigmaIetaIeta").c_str(), &deepSuperCluster_sigmaIetaIeta, &b_deepSuperCluster_sigmaIetaIeta);
   tree->SetBranchAddress(std::string(superClusterVal+"_sigmaIetaIphi").c_str(), &deepSuperCluster_sigmaIetaIphi, &b_deepSuperCluster_sigmaIetaIphi);
   tree->SetBranchAddress(std::string(superClusterVal+"_sigmaIphiIphi").c_str(), &deepSuperCluster_sigmaIphiIphi, &b_deepSuperCluster_sigmaIphiIphi);
   tree->SetBranchAddress(std::string(superClusterVal+"_full5x5_r9").c_str(), &deepSuperCluster_full5x5_r9, &b_deepSuperCluster_full5x5_r9);
   tree->SetBranchAddress(std::string(superClusterVal+"_full5x5_sigmaIetaIeta").c_str(), &deepSuperCluster_full5x5_sigmaIetaIeta, &b_deepSuperCluster_full5x5_sigmaIetaIeta);
   tree->SetBranchAddress(std::string(superClusterVal+"_full5x5_sigmaIetaIphi").c_str(), &deepSuperCluster_full5x5_sigmaIetaIphi, &b_deepSuperCluster_full5x5_sigmaIetaIphi);
   tree->SetBranchAddress(std::string(superClusterVal+"_full5x5_sigmaIphiIphi").c_str(), &deepSuperCluster_full5x5_sigmaIphiIphi, &b_deepSuperCluster_full5x5_sigmaIphiIphi);
   tree->SetBranchAddress("pfCluster_phi", &pfCluster_phi, &b_pfCluster_phi);
   tree->SetBranchAddress("pfCluster_eta", &pfCluster_eta, &b_pfCluster_eta);
   tree->SetBranchAddress("pfCluster_energy", &pfCluster_energy, &b_pfCluster_energy);  
}

//split
vector<string> split(const string &text, char sep)
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

//erase string 
void removeSubstrs(std::string& s, const std::string& p) 
{
   std::string::size_type n = p.length();
   for (std::string::size_type i = s.find(p);
        i != std::string::npos;
        i = s.find(p))
      s.erase(i, n);
}

//getBinOpts
std::vector<std::pair<std::string,std::vector<double>>> getBinOpts(const edm::ParameterSet &histOpt)
{
    std::vector<std::pair<std::string,std::vector<double>>> bins;
    bins.resize(38); 
    bins[0] = std::make_pair(std::string("nVtxBins_Barrel"),histOpt.getParameter<std::vector<double>>( "nVtxBins_Barrel" ));
    bins[1] = std::make_pair(std::string("nVtxBins_Endcap"),histOpt.getParameter<std::vector<double>>( "nVtxBins_Endcap" ));
    bins[2] = std::make_pair(std::string("RhoBins_Barrel"),histOpt.getParameter<std::vector<double>>( "RhoBins_Barrel" ));
    bins[3] = std::make_pair(std::string("RhoBins_Endcap"),histOpt.getParameter<std::vector<double>>( "RhoBins_Endcap" )); 
    bins[4] = std::make_pair(std::string("nPFClustersBins"),histOpt.getParameter<std::vector<double>>( "nPFClustersBins" )); 
    bins[5] = std::make_pair(std::string("nPFClustersBins_Barrel"),histOpt.getParameter<std::vector<double>>( "nPFClustersBins_Barrel" ));
    bins[6] = std::make_pair(std::string("nPFClustersBins_Endcap"),histOpt.getParameter<std::vector<double>>( "nPFClustersBins_Endcap" )); 
    bins[7] = std::make_pair(std::string("EnergyBins_Barrel"),histOpt.getParameter<std::vector<double>>( "EnergyBins_Barrel" ));
    bins[8] = std::make_pair(std::string("EnergyBins_Endcap"),histOpt.getParameter<std::vector<double>>( "EnergyBins_Endcap" ));
    bins[9] = std::make_pair(std::string("EoEtrueBins_Barrel"),histOpt.getParameter<std::vector<double>>( "EoEtrueBins_Barrel" ));
    bins[10] = std::make_pair(std::string("EoEtrueBins_Endcap"),histOpt.getParameter<std::vector<double>>( "EoEtrueBins_Endcap" ));
    bins[11] = std::make_pair(std::string("EoEgenBins_Barrel"),histOpt.getParameter<std::vector<double>>( "EoEgenBins_Barrel" ));
    bins[12] = std::make_pair(std::string("EoEgenBins_Endcap"),histOpt.getParameter<std::vector<double>>( "EoEgenBins_Endcap" )); 
    bins[13] = std::make_pair(std::string("EtBins_Barrel"),histOpt.getParameter<std::vector<double>>( "EtBins_Barrel" ));
    bins[14] = std::make_pair(std::string("EtBins_Endcap"),histOpt.getParameter<std::vector<double>>( "EtBins_Endcap" )); 
    bins[15] = std::make_pair(std::string("EtaBins"),histOpt.getParameter<std::vector<double>>( "EtaBins" ));
    bins[16] = std::make_pair(std::string("PhiBins_Barrel"),histOpt.getParameter<std::vector<double>>( "PhiBins_Barrel" ));
    bins[17] = std::make_pair(std::string("PhiBins_Endcap"),histOpt.getParameter<std::vector<double>>( "PhiBins_Endcap" ));
    bins[18] = std::make_pair(std::string("EtaWidthBins_Barrel"),histOpt.getParameter<std::vector<double>>( "EtaWidthBins_Barrel" ));
    bins[19] = std::make_pair(std::string("EtaWidthBins_Endcap"),histOpt.getParameter<std::vector<double>>( "EtaWidthBins_Endcap" ));
    bins[20] = std::make_pair(std::string("PhiWidthBins_Barrel"),histOpt.getParameter<std::vector<double>>( "PhiWidthBins_Barrel" ));
    bins[21] = std::make_pair(std::string("PhiWidthBins_Endcap"),histOpt.getParameter<std::vector<double>>( "PhiWidthBins_Endcap" ));
    bins[22] = std::make_pair(std::string("R9Bins_Barrel"),histOpt.getParameter<std::vector<double>>( "R9Bins_Barrel" ));
    bins[23] = std::make_pair(std::string("R9Bins_Endcap"),histOpt.getParameter<std::vector<double>>( "R9Bins_Endcap" ));  
    bins[24] = std::make_pair(std::string("full5x5_R9Bins_Barrel"),histOpt.getParameter<std::vector<double>>( "full5x5_R9Bins_Barrel" ));
    bins[25] = std::make_pair(std::string("full5x5_R9Bins_Endcap"),histOpt.getParameter<std::vector<double>>( "full5x5_R9Bins_Endcap" )); 
    bins[26] = std::make_pair(std::string("SigmaIetaIetaBins_Barrel"),histOpt.getParameter<std::vector<double>>( "SigmaIetaIetaBins_Barrel" ));
    bins[27] = std::make_pair(std::string("SigmaIetaIetaBins_Endcap"),histOpt.getParameter<std::vector<double>>( "SigmaIetaIetaBins_Endcap" ));  
    bins[28] = std::make_pair(std::string("full5x5_SigmaIetaIetaBins_Barrel"),histOpt.getParameter<std::vector<double>>( "full5x5_SigmaIetaIetaBins_Barrel" ));
    bins[29] = std::make_pair(std::string("full5x5_SigmaIetaIetaBins_Endcap"),histOpt.getParameter<std::vector<double>>( "full5x5_SigmaIetaIetaBins_Endcap" ));  
    bins[30] = std::make_pair(std::string("SigmaIetaIphiBins_Barrel"),histOpt.getParameter<std::vector<double>>( "SigmaIetaIphiBins_Barrel" ));
    bins[31] = std::make_pair(std::string("SigmaIetaIphiBins_Endcap"),histOpt.getParameter<std::vector<double>>( "SigmaIetaIphiBins_Endcap" ));  
    bins[32] = std::make_pair(std::string("full5x5_SigmaIetaIphiBins_Barrel"),histOpt.getParameter<std::vector<double>>( "full5x5_SigmaIetaIphiBins_Barrel" ));
    bins[33] = std::make_pair(std::string("full5x5_SigmaIetaIphiBins_Endcap"),histOpt.getParameter<std::vector<double>>( "full5x5_SigmaIetaIphiBins_Endcap" )); 
    bins[34] = std::make_pair(std::string("SigmaIphiIphiBins_Barrel"),histOpt.getParameter<std::vector<double>>( "SigmaIphiIphiBins_Barrel" ));
    bins[35] = std::make_pair(std::string("SigmaIphiIphiBins_Endcap"),histOpt.getParameter<std::vector<double>>( "SigmaIphiIphiBins_Endcap" ));  
    bins[36] = std::make_pair(std::string("full5x5_SigmaIphiIphiBins_Barrel"),histOpt.getParameter<std::vector<double>>( "full5x5_SigmaIphiIphiBins_Barrel" ));
    bins[37] = std::make_pair(std::string("full5x5_SigmaIphiIphiBins_Endcap"),histOpt.getParameter<std::vector<double>>( "full5x5_SigmaIphiIphiBins_Endcap" )); 
    return bins; 
}

//findOption
int findOption(std::string var, std::vector<std::pair<std::string,std::vector<double>>> binOpts)
{
   for(unsigned int iBin=0; iBin<binOpts.size(); iBin++){
       if(binOpts[iBin].first.find(var.c_str())!=std::string::npos)
          return iBin;
   } 
   std::cout << "----> findOption: " << var << " ---> NO BIN!" << std::endl;
   return -1; 
}

//set efficiencies
void setEfficiencies()
{
   eff_SuperCluster_vs_EtaCalo = new TEfficiency(*h_Eta_Calo_old,*h_Eta_Calo_Denum);
   eff_DeepSuperCluster_vs_EtaCalo = new TEfficiency(*h_Eta_Calo_new,*h_Eta_Calo_Denum);   
   eff_SuperCluster_vs_EtCalo_EB = new TEfficiency(*h_Et_Calo_EB_old,*h_Et_Calo_EB_Denum);
   eff_DeepSuperCluster_vs_EtCalo_EB = new TEfficiency(*h_Et_Calo_EB_new,*h_Et_Calo_EB_Denum);
   eff_SuperCluster_vs_EtCalo_EE = new TEfficiency(*h_Et_Calo_EE_old,*h_Et_Calo_EE_Denum);
   eff_DeepSuperCluster_vs_EtCalo_EE = new TEfficiency(*h_Et_Calo_EE_new,*h_Et_Calo_EE_Denum);
}

//set histograms
void setHistograms(std::vector<std::pair<std::string,std::vector<double>>> binOpts)
{
   std::vector<std::pair<std::string,std::vector<double>>> binOpts_tmp = binOpts; 
   for(unsigned int iBin=0; iBin<binOpts.size(); iBin++)
   {
       if(binOpts[iBin].first.find("EtaBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")==std::string::npos && binOpts[iBin].first.find("Endcap")==std::string::npos){ 
             removeSubstrs(binOpts[iBin].first,std::string("Bins")); 
             h_Eta_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_Calo_Denum = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_Denum").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_Denum").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_Calo_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_SC").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_SC").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_Calo_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_DeepSC").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_DeepSC").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             prof_EoEtrue_vs_Eta_Calo_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Eta_Calo_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Eta_Gen_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Eta_Gen_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);    
             prof_EoEtrue_vs_Eta_Seed_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Eta_Seed_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Eta_Seed_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Eta_Seed_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             h_Eta_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Eta_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("nPFClustersBins")!=std::string::npos){    
          if(binOpts[iBin].first.find("Barrel")==std::string::npos && binOpts[iBin].first.find("Endcap")==std::string::npos){ 
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));  
             h_nPFClusters_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_nPFClusters_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_nPFClusters_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_nPFClusters_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }
       
       if(binOpts[iBin].first.find("nVtxBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             prof_EoEtrue_vs_nVtx_EB_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_nVtx_EB_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_nVtx_EB_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_nVtx_EB_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             prof_EoEtrue_vs_nVtx_EE_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_nVtx_EE_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_nVtx_EE_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_nVtx_EE_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
          }
       }

       if(binOpts[iBin].first.find("RhoBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             prof_EoEtrue_vs_Rho_EB_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Rho_EB_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Rho_EB_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Rho_EB_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             prof_EoEtrue_vs_Rho_EE_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Rho_EE_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Rho_EE_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Rho_EE_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
          }
       }
       if(binOpts[iBin].first.find("EnergyBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_Energy_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             prof_EoEtrue_vs_Energy_Calo_EB_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Energy_Calo_EB_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Energy_Gen_EB_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Energy_Gen_EB_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             h_Energy_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_Energy_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             prof_EoEtrue_vs_Energy_Calo_EE_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Energy_Calo_EE_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Energy_Gen_EE_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Energy_Gen_EE_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.); 
             h_Energy_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Energy_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("EoEtrueBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_EoEtrue_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEtrue_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEtrue_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEtrue_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_EoEtrue_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEtrue_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEtrue_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEtrue_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("EoEgenBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_EoEgen_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEgen_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEgen_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEgen_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_EoEgen_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEgen_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEgen_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EoEgen_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("EtBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_Et_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_Calo_EB_Denum = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_EB_Denum").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_EB_Denum").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_Calo_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_EB_SC").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_EB_SC").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_Calo_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_EB_DeepSC").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_EB_DeepSC").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             prof_EoEtrue_vs_Et_Calo_EB_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Et_Calo_EB_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Gen_EB_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Gen_EB_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.); 
             prof_EoEtrue_vs_Et_Seed_EB_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EB_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Et_Seed_EB_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EB_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Seed_EB_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EB_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Seed_EB_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EB_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EB").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);    
             h_Et_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_Et_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_Calo_EE_Denum = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_EE_Denum").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_EE_Denum").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_Calo_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_EE_SC").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_EE_SC").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_Calo_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_Calo_EE_DeepSC").c_str(),std::string("h_"+binOpts[iBin].first+"_Calo_EE_DeepSC").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             prof_EoEtrue_vs_Et_Calo_EE_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Et_Calo_EE_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Calo_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Gen_EE_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Gen_EE_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Gen_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);     
             prof_EoEtrue_vs_Et_Seed_EE_old  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EE_SC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEtrue_vs_Et_Seed_EE_new  = new TProfile(std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EE_DeepSC").c_str(), std::string("prof_EoEtrue_vs_"+binOpts[iBin].first+"_Seed_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Seed_EE_old  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EE_SC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             prof_EoEgen_vs_Et_Seed_EE_new  = new TProfile(std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EE_DeepSC").c_str(), std::string("prof_EoEgen_vs_"+binOpts[iBin].first+"_Seed_EE").c_str(),binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2],0.,10.);
             h_Et_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Et_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       } 

       if(binOpts[iBin].first.find("PhiBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_Phi_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_Phi_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_Phi_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("EtaWidthBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_EtaWidth_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_EtaWidth_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_EtaWidth_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("PhiWidthBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_PhiWidth_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_PhiWidth_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_PhiWidth_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("full5x5_R9Bins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins")); 
             h_full5x5_R9_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_R9_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_R9_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("R9Bins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_R9_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_R9_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_R9_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("full5x5_SigmaIetaIetaBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_sigmaIetaIeta_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_sigmaIetaIeta_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("full5x5_SigmaIetaIphiBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_sigmaIetaIphi_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_sigmaIetaIphi_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }  

       if(binOpts[iBin].first.find("SigmaIetaIetaBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_sigmaIetaIeta_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_sigmaIetaIeta_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIeta_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("SigmaIetaIphiBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_sigmaIetaIphi_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_sigmaIetaIphi_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIetaIphi_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("full5x5_SigmaIphiIphiBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_sigmaIphiIphi_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_full5x5_sigmaIphiIphi_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }

       if(binOpts[iBin].first.find("SigmaIphiIphiBins")!=std::string::npos){      
          if(binOpts[iBin].first.find("Barrel")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Barrel"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_sigmaIphiIphi_EB_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EB_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EB_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
          if(binOpts[iBin].first.find("Endcap")!=std::string::npos){
             removeSubstrs(binOpts[iBin].first,std::string("_Endcap"));
             removeSubstrs(binOpts[iBin].first,std::string("Bins"));
             h_sigmaIphiIphi_EE_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_SC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first).c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_seedMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_seedMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_seedMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_seedMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_caloMatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_caloMatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloMatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloMatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_caloUnmatched_old = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_SC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
             h_sigmaIphiIphi_EE_caloUnmatched_new = new TH1F(std::string("h_"+binOpts[iBin].first+"_EE_caloUnmatched_DeepSC").c_str(), std::string("h_"+binOpts[iBin].first+"_caloUnmatched").c_str(), binOpts[iBin].second[0], binOpts[iBin].second[1], binOpts[iBin].second[2]);
          }
       }
   } 

   binOpts = binOpts_tmp;
 
   int nBins_Eta = binOpts[findOption(std::string("EtaBins"),binOpts)].second[0];
   EoEtrue_vs_Eta_Calo_old.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEtrue_vs_Eta_Calo_old[iBin] = new TH1F(std::string("EoEtrue_vs_Eta_Calo_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Eta_Calo_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Eta_Calo_new.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEtrue_vs_Eta_Calo_new[iBin] = new TH1F(std::string("EoEtrue_vs_Eta_Calo_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Eta_Calo_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEgen_vs_Eta_Gen_old.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEgen_vs_Eta_Gen_old[iBin] = new TH1F(std::string("EoEgen_vs_Eta_Gen_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Eta_Gen_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Eta_Gen_new.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEgen_vs_Eta_Gen_new[iBin] = new TH1F(std::string("EoEgen_vs_Eta_Gen_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Eta_Gen_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEtrue_vs_Eta_Seed_old.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEtrue_vs_Eta_Seed_old[iBin] = new TH1F(std::string("EoEtrue_vs_Eta_Seed_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Eta_Seed_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Eta_Seed_new.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEtrue_vs_Eta_Seed_new[iBin] = new TH1F(std::string("EoEtrue_vs_Eta_Seed_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Eta_Seed_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEgen_vs_Eta_Seed_old.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEgen_vs_Eta_Seed_old[iBin] = new TH1F(std::string("EoEgen_vs_Eta_Seed_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Eta_Seed_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Eta_Seed_new.resize(nBins_Eta);
   for(int iBin=0; iBin<nBins_Eta; iBin++)
       EoEgen_vs_Eta_Seed_new[iBin] = new TH1F(std::string("EoEgen_vs_Eta_Seed_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Eta_Seed_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 

   int nBins_Et_EB = binOpts[findOption(std::string("EtBins_Barrel"),binOpts)].second[0];
   int nBins_Et_EE = binOpts[findOption(std::string("EtBins_Endcap"),binOpts)].second[0];
   EoEtrue_vs_Et_Calo_EB_old.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEtrue_vs_Et_Calo_EB_old[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Calo_EB_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Calo_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Et_Calo_EB_new.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEtrue_vs_Et_Calo_EB_new[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Calo_EB_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Calo_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Et_Calo_EE_old.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEtrue_vs_Et_Calo_EE_old[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Calo_EE_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Calo_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Et_Calo_EE_new.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEtrue_vs_Et_Calo_EE_new[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Calo_EE_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Calo_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEgen_vs_Et_Gen_EB_old.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEgen_vs_Et_Gen_EB_old[iBin] = new TH1F(std::string("EoEgen_vs_Et_Gen_EB_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Gen_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Et_Gen_EB_new.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEgen_vs_Et_Gen_EB_new[iBin] = new TH1F(std::string("EoEgen_vs_Et_Gen_EB_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Gen_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Et_Gen_EE_old.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEgen_vs_Et_Gen_EE_old[iBin] = new TH1F(std::string("EoEgen_vs_Et_Gen_EE_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Gen_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Et_Gen_EE_new.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEgen_vs_Et_Gen_EE_new[iBin] = new TH1F(std::string("EoEgen_vs_Et_Gen_EE_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Gen_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEtrue_vs_Et_Seed_EB_old.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEtrue_vs_Et_Seed_EB_old[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Seed_EB_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Seed_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Et_Seed_EB_new.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEtrue_vs_Et_Seed_EB_new[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Seed_EB_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Seed_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Et_Seed_EE_old.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEtrue_vs_Et_Seed_EE_old[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Seed_EE_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Seed_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Et_Seed_EE_new.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEtrue_vs_Et_Seed_EE_new[iBin] = new TH1F(std::string("EoEtrue_vs_Et_Seed_EE_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Et_Seed_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEgen_vs_Et_Seed_EB_old.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEgen_vs_Et_Seed_EB_old[iBin] = new TH1F(std::string("EoEgen_vs_Et_Seed_EB_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Seed_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Et_Seed_EB_new.resize(nBins_Et_EB);
   for(int iBin=0; iBin<nBins_Et_EB; iBin++)
       EoEgen_vs_Et_Seed_EB_new[iBin] = new TH1F(std::string("EoEgen_vs_Et_Seed_EB_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Seed_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Et_Seed_EE_old.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEgen_vs_Et_Seed_EE_old[iBin] = new TH1F(std::string("EoEgen_vs_Et_Seed_EE_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Seed_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Et_Seed_EE_new.resize(nBins_Et_EE);
   for(int iBin=0; iBin<nBins_Et_EE; iBin++)
       EoEgen_vs_Et_Seed_EE_new[iBin] = new TH1F(std::string("EoEgen_vs_Et_Seed_EE_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Et_Seed_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);  

   int nBins_Energy_EB = binOpts[findOption(std::string("EnergyBins_Barrel"),binOpts)].second[0];
   int nBins_Energy_EE = binOpts[findOption(std::string("EnergyBins_Endcap"),binOpts)].second[0];
   EoEtrue_vs_Energy_Calo_EB_old.resize(nBins_Energy_EB);
   for(int iBin=0; iBin<nBins_Energy_EB; iBin++)
       EoEtrue_vs_Energy_Calo_EB_old[iBin] = new TH1F(std::string("EoEtrue_vs_Energy_Calo_EB_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Energy_Calo_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Energy_Calo_EB_new.resize(nBins_Energy_EB);
   for(int iBin=0; iBin<nBins_Energy_EB; iBin++)
       EoEtrue_vs_Energy_Calo_EB_new[iBin] = new TH1F(std::string("EoEtrue_vs_Energy_Calo_EB_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Energy_Calo_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Energy_Calo_EE_old.resize(nBins_Energy_EE);
   for(int iBin=0; iBin<nBins_Energy_EE; iBin++)
       EoEtrue_vs_Energy_Calo_EE_old[iBin] = new TH1F(std::string("EoEtrue_vs_Energy_Calo_EE_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Energy_Calo_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Energy_Calo_EE_new.resize(nBins_Energy_EE);
   for(int iBin=0; iBin<nBins_Energy_EE; iBin++)
       EoEtrue_vs_Energy_Calo_EE_new[iBin] = new TH1F(std::string("EoEtrue_vs_Energy_Calo_EE_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Energy_Calo_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);   
   EoEgen_vs_Energy_Gen_EB_old.resize(nBins_Energy_EB);
   for(int iBin=0; iBin<nBins_Energy_EB; iBin++)
       EoEgen_vs_Energy_Gen_EB_old[iBin] = new TH1F(std::string("EoEgen_vs_Energy_Calo_EB_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Energy_Calo_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Energy_Gen_EB_new.resize(nBins_Energy_EB);
   for(int iBin=0; iBin<nBins_Energy_EB; iBin++)
       EoEgen_vs_Energy_Gen_EB_new[iBin] = new TH1F(std::string("EoEgen_vs_Energy_Calo_EB_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Energy_Calo_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Energy_Gen_EE_old.resize(nBins_Energy_EE);
   for(int iBin=0; iBin<nBins_Energy_EE; iBin++)
       EoEgen_vs_Energy_Gen_EE_old[iBin] = new TH1F(std::string("EoEgen_vs_Energy_Calo_EE_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Energy_Calo_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Energy_Gen_EE_new.resize(nBins_Energy_EE);
   for(int iBin=0; iBin<nBins_Energy_EE; iBin++)
       EoEgen_vs_Energy_Gen_EE_new[iBin] = new TH1F(std::string("EoEgen_vs_Energy_Calo_EE_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Energy_Calo_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);   
   
   int nBins_nVtx_EB = binOpts[findOption(std::string("nVtxBins_Barrel"),binOpts)].second[0];
   int nBins_nVtx_EE = binOpts[findOption(std::string("nVtxBins_Endcap"),binOpts)].second[0];
   EoEtrue_vs_nVtx_EB_old.resize(nBins_nVtx_EB);
   for(int iBin=0; iBin<nBins_nVtx_EB; iBin++)
       EoEtrue_vs_nVtx_EB_old[iBin] = new TH1F(std::string("EoEtrue_vs_nVtx_EB_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_nVtx_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_nVtx_EB_new.resize(nBins_nVtx_EB);
   for(int iBin=0; iBin<nBins_nVtx_EB; iBin++)
       EoEtrue_vs_nVtx_EB_new[iBin] = new TH1F(std::string("EoEtrue_vs_nVtx_EB_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_nVtx_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_nVtx_EE_old.resize(nBins_nVtx_EE);
   for(int iBin=0; iBin<nBins_nVtx_EE; iBin++)
       EoEtrue_vs_nVtx_EE_old[iBin] = new TH1F(std::string("EoEtrue_vs_nVtx_EE_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_nVtx_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_nVtx_EE_new.resize(nBins_nVtx_EE);
   for(int iBin=0; iBin<nBins_nVtx_EE; iBin++)
       EoEtrue_vs_nVtx_EE_new[iBin] = new TH1F(std::string("EoEtrue_vs_nVtx_EE_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_nVtx_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEgen_vs_nVtx_EB_old.resize(nBins_nVtx_EB);
   for(int iBin=0; iBin<nBins_nVtx_EB; iBin++)
       EoEgen_vs_nVtx_EB_old[iBin] = new TH1F(std::string("EoEgen_vs_nVtx_EB_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_nVtx_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_nVtx_EB_new.resize(nBins_nVtx_EB);
   for(int iBin=0; iBin<nBins_nVtx_EB; iBin++)
       EoEgen_vs_nVtx_EB_new[iBin] = new TH1F(std::string("EoEgen_vs_nVtx_EB_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_nVtx_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_nVtx_EE_old.resize(nBins_nVtx_EE);
   for(int iBin=0; iBin<nBins_nVtx_EE; iBin++)
       EoEgen_vs_nVtx_EE_old[iBin] = new TH1F(std::string("EoEgen_vs_nVtx_EE_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_nVtx_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_nVtx_EE_new.resize(nBins_nVtx_EE);
   for(int iBin=0; iBin<nBins_nVtx_EE; iBin++)
       EoEgen_vs_nVtx_EE_new[iBin] = new TH1F(std::string("EoEgen_vs_nVtx_EE_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_nVtx_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
  
   int nBins_Rho_EB = binOpts[findOption(std::string("RhoBins_Barrel"),binOpts)].second[0]; 
   int nBins_Rho_EE = binOpts[findOption(std::string("RhoBins_Endcap"),binOpts)].second[0];
   EoEtrue_vs_Rho_EB_old.resize(nBins_Rho_EB);
   for(int iBin=0; iBin<nBins_Rho_EB; iBin++)
       EoEtrue_vs_Rho_EB_old[iBin] = new TH1F(std::string("EoEtrue_vs_Rho_EB_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Rho_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Rho_EB_new.resize(nBins_Rho_EB);
   for(int iBin=0; iBin<nBins_Rho_EB; iBin++)
       EoEtrue_vs_Rho_EB_new[iBin] = new TH1F(std::string("EoEtrue_vs_Rho_EB_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Rho_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Rho_EE_old.resize(nBins_Rho_EE);
   for(int iBin=0; iBin<nBins_Rho_EE; iBin++)
       EoEtrue_vs_Rho_EE_old[iBin] = new TH1F(std::string("EoEtrue_vs_Rho_EE_old_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Rho_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEtrue_vs_Rho_EE_new.resize(nBins_Rho_EE);
   for(int iBin=0; iBin<nBins_Rho_EE; iBin++)
       EoEtrue_vs_Rho_EE_new[iBin] = new TH1F(std::string("EoEtrue_vs_Rho_EE_new_"+to_string(iBin)).c_str(), std::string("EoEtrue_vs_Rho_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.); 
   EoEgen_vs_Rho_EB_old.resize(nBins_Rho_EB);
   for(int iBin=0; iBin<nBins_Rho_EB; iBin++)
       EoEgen_vs_Rho_EB_old[iBin] = new TH1F(std::string("EoEgen_vs_Rho_EB_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Rho_EB_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Rho_EB_new.resize(nBins_Rho_EB);
   for(int iBin=0; iBin<nBins_Rho_EB; iBin++)
       EoEgen_vs_Rho_EB_new[iBin] = new TH1F(std::string("EoEgen_vs_Rho_EB_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Rho_EB_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Rho_EE_old.resize(nBins_Rho_EE);
   for(int iBin=0; iBin<nBins_Rho_EE; iBin++)
       EoEgen_vs_Rho_EE_old[iBin] = new TH1F(std::string("EoEgen_vs_Rho_EE_old_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Rho_EE_old_"+to_string(iBin)).c_str(), 10000, 0., 10.);
   EoEgen_vs_Rho_EE_new.resize(nBins_Rho_EE);
   for(int iBin=0; iBin<nBins_Rho_EE; iBin++)
       EoEgen_vs_Rho_EE_new[iBin] = new TH1F(std::string("EoEgen_vs_Rho_EE_new_"+to_string(iBin)).c_str(), std::string("EoEgen_vs_Rho_EE_new_"+to_string(iBin)).c_str(), 10000, 0., 10.);

   
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

double mycruijff(double* x, double* par) 
{
  Double_t m = x[0];
  Double_t m0 = par[0];
  Double_t sigmaL = par[1];
  Double_t sigmaR = par[2];
  Double_t alphaL = par[3];
  Double_t alphaR = par[4];
  Double_t N = par[5];
  Double_t dx =  (m-m0) ;
  Double_t sigma = dx<0 ? sigmaL: sigmaR ;
  Double_t alpha = dx<0 ? alphaL: alphaR ;
  Double_t f = 2*sigma*sigma + alpha*dx*dx ; 
  
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
   std::cout << "fitHisto: " << hist->GetName() << std::endl;
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

    //std::cout << "EFFECTIVE SIGMA: " << kbin << " - " << jbin << " - " << nb << " - " << (float)hist->Integral(kbin,jbin)/(float)hist->Integral() << std::endl;
    return std::make_pair(computeMean(hist,kbin,jbin),widmin);
}

std::pair<TGraphErrors*,TGraphErrors*> makeFitProfile(std::vector<TH1F*>* vecHist,double min, double max, std::string xTitle, std::string fitFunction_="cruijff", bool doEffective = false)
{
   TF1* doubleCB;
   int nBins = vecHist->size();
   double x[nBins], x_error[nBins], y_mean[nBins], y_meanError[nBins], y_sigma[nBins], y_sigmaError[nBins];
   float delta = fabs(max-min)/nBins; 

   for(unsigned int iBin=0; iBin<vecHist->size(); iBin++)
   {
       x[iBin] = min + (iBin+0.5)*delta; 
       x_error[iBin] = 0.;
       
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

void drawHisto(TH1F* h_old, TH1F* h_new, std::string x_label, std::string drawType, std::string Name, bool log, bool fit=false, std::string fitFunc_="cruijff", std::string refLegend="Mustache", std::string valLegend="DeepSC")
{

   h_old->Scale(1./h_old->GetEntries());
   h_new->Scale(1./h_new->GetEntries());

   TFitResultPtr frp_old;
   TFitResultPtr frp_new;
   TF1* doubleCB_old = new TF1();
   TF1* doubleCB_new = new TF1();

   if(fit){
      TH1F* h_old_clone = (TH1F*)h_old->Clone();  
      TH1F* h_new_clone = (TH1F*)h_new->Clone();    

      doubleCB_old = fitHisto(h_old_clone, fitFunc_);
      doubleCB_old->SetLineColor(kRed+1);

      doubleCB_new = fitHisto(h_new_clone, fitFunc_);
      doubleCB_new->SetLineColor(kBlue+1);
   }

   //h_old->SetTitle("");
   //h_new->SetTitle("");

   std::vector<float> maxima;
   maxima.resize(4);
   maxima[0] = h_old->GetMaximum(); 
   maxima[1] = h_new->GetMaximum();
   if(fit) maxima[2] = doubleCB_old->GetMaximum();
   if(fit) maxima[3] = doubleCB_new->GetMaximum(); 
   std::sort(maxima.begin(),maxima.end());  
   
   h_old->SetMaximum(maxima.at(maxima.size()-1)*1.05);
   h_old->SetLineColor(kRed+1);
   h_old->SetMarkerColor(kRed+1);
   h_old->SetLineWidth(2);
   h_old->GetXaxis()->SetTitle(x_label.c_str());

   h_new->SetLineColor(kBlue+1);
   h_new->SetMarkerColor(kBlue+1);
   h_new->SetLineWidth(2);

   /*float maximum = -1.;
   if( h_old -> GetMaximum() > h_new -> GetMaximum()) maximum=h_old -> GetMaximum();
   else maximum=h_new -> GetMaximum();  
   h_old -> SetMaximum( 1.1*maximum );*/
    
   TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);

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
   if(!fit){
      TLegend* legend = new TLegend(0.60, 0.82, 0.75, 0.94);
      legend -> SetFillColor(kWhite);
      legend -> SetFillStyle(1000);
      legend -> SetLineWidth(0);
      legend -> SetLineColor(kWhite);
      legend -> SetTextFont(42);  
      legend -> SetTextSize(0.04); 
      legend -> AddEntry(h_old,refLegend.c_str(),"L");
      legend -> AddEntry(h_new,valLegend.c_str(),"L");
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
   }else{
      TLegend* legend = new TLegend(0.57, 0.77, 0.72, 0.89);
      legend -> SetFillColor(kWhite);
      legend -> SetFillStyle(1000);
      legend -> SetLineWidth(0);
      legend -> SetLineColor(kWhite);
      legend -> SetTextFont(42);  
      legend -> SetTextSize(0.05);
      legend -> AddEntry(h_old,std::string(refLegend+": "+to_string(doubleCB_old->GetParameter(0))+" +/- "+to_string((doubleCB_old->GetParameter(2)+doubleCB_old->GetParameter(1))/2.)).c_str(),"L");
      legend -> AddEntry(h_new,std::string(valLegend+": "+to_string(doubleCB_new->GetParameter(0))+" +/- "+to_string((doubleCB_new->GetParameter(2)+doubleCB_new->GetParameter(1))/2.)).c_str(),"L");
      h_old->Draw(std::string(drawType).c_str());    
      h_new->Draw(std::string(drawType+",sames").c_str());
      doubleCB_old->Draw("L,same");
      doubleCB_new->Draw("L,same");
      legend -> Draw("same");   
      gPad -> Update();
      st_old= (TPaveStats*)(h_old->GetListOfFunctions()->FindObject("stats"));
      st_old->SetX1NDC(0.); //new x start position
      st_old->SetX2NDC(0.); //new x end position
      st_old->SetY1NDC(0.); //new y start position
      st_old->SetY2NDC(0.); //new y end position
      st_old->SetTextColor(kBlack);
      st_old->Draw("sames");
      gPad -> Update();
      st_new= (TPaveStats*)(h_new->GetListOfFunctions()->FindObject("stats"));
      st_new->SetX1NDC(0.); //new x start position
      st_new->SetX2NDC(0.); //new x end position
      st_new->SetY1NDC(0.); //new y start position
      st_new->SetY2NDC(0.); //new y end position
      st_new->SetTextColor(kBlack);
      st_new->Draw("sames");
   }
   cDown->cd();
    
   TH1F* histo_ratio=(TH1F*)h_new->Clone("histo_ratio");
   if(histo_ratio->GetSumw2N()<=0) histo_ratio->Sumw2();
   histo_ratio->Divide(h_old);
    
   histo_ratio -> GetXaxis() -> SetTitle(x_label.c_str());
   histo_ratio -> GetYaxis() -> SetTitle(std::string(valLegend+"/"+refLegend).c_str());
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

   delete histo_ratio;
}

TGraphErrors makeRatioGraph(TGraphErrors* gr_SuperCluster, TGraphErrors* gr_DeepSuperCluster)
{
   TGraphErrors* gr_ratio = new TGraphErrors();
   double x,y_ref,y_val,yErr_ref,yErr_val;
   for(int i = 0; i<gr_SuperCluster->GetN(); i++){ 
       gr_SuperCluster->GetPoint(i,x,y_ref); 
       yErr_ref = gr_SuperCluster->GetErrorY(i); 
       gr_DeepSuperCluster->GetPoint(i,x,y_val);  
       yErr_val = gr_DeepSuperCluster->GetErrorY(i);     
       if(y_ref>0. && y_val>0.){
          gr_ratio->SetPoint(i,x,y_val/y_ref);
          gr_ratio->SetPointError(i,0.,y_val/y_ref*sqrt((yErr_ref/y_ref)*(yErr_ref/y_ref)+(yErr_val/y_val)*(yErr_val/y_val)));
       }else{
          gr_ratio->SetPoint(i,x,-1.);
          gr_ratio->SetPointError(i,0.,0.);
       }   
   } 
   return *gr_ratio; 
}

void drawGraph(TGraphErrors* gr_SuperCluster, TGraphErrors* gr_DeepSuperCluster, std::string xtitle, std::string ytitle, std::string Name, std::string refLegend="Mustache", std::string valLegend="DeepSC",float y_min=-1., float y_max=-1.)
{ 
   gStyle->SetOptStat(0000);  
   gr_SuperCluster->SetTitle(Name.c_str());
   gr_SuperCluster->SetLineColor(kRed+1);
   gr_SuperCluster->SetMarkerStyle(20);
   gr_SuperCluster->SetMarkerSize(0.5);
   gr_SuperCluster->SetMarkerColor(kRed+1);
   gr_SuperCluster->SetLineWidth(2);
   gr_SuperCluster->GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_SuperCluster->GetYaxis()->SetTitle(ytitle.c_str()); 
   
   gr_DeepSuperCluster->SetLineColor(kBlue+1);
   gr_DeepSuperCluster->SetMarkerStyle(20);
   gr_DeepSuperCluster->SetMarkerSize(0.5);
   gr_DeepSuperCluster->SetMarkerColor(kBlue+1);
   gr_DeepSuperCluster->SetLineWidth(2);

   float min = y_min;
   float max = y_max;
   if(y_min<0. || y_max<0.){
      min = 0.9*computeRange(gr_SuperCluster).first; 
      max = 1.1*computeRange(gr_SuperCluster).second;
   }
   gr_SuperCluster->GetYaxis()->SetRangeUser(min,max);  

   TLegend* legend = new TLegend(0.799, 0.77, 0.999, 0.95);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.03);
   legend -> AddEntry(gr_SuperCluster,refLegend.c_str(),"L");
   legend -> AddEntry(gr_DeepSuperCluster,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();

   TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.36,1.00,1.00);
   TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.36);
   
   cUp->SetBottomMargin(0.01); 
   cDown->SetTopMargin(0.01); 
   cDown->SetBottomMargin(0.2); 
    
   cUp->Draw();
   cDown->Draw();
     
   cUp->cd();
   gr_SuperCluster->Draw("AP");
   gr_DeepSuperCluster->Draw("P, same");
   legend -> Draw("same");

   cDown->cd();
    
   TGraphErrors gr_ratio = makeRatioGraph(gr_SuperCluster,gr_DeepSuperCluster); 
   gr_ratio.GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_ratio.GetYaxis() -> SetTitle(std::string(valLegend+"/"+refLegend).c_str());
   gr_ratio.GetYaxis() -> SetRangeUser(0.5,1.5);
   gr_ratio.SetMarkerColor(kBlack);
   gr_ratio.SetMarkerStyle(20);
   gr_ratio.SetMarkerSize(0.5);
   gr_ratio.SetTitle("");
   gr_ratio.GetXaxis() -> SetLabelSize(0.07);
   gr_ratio.GetYaxis() -> SetLabelSize(0.07);
   gr_ratio.GetXaxis() -> SetTitleSize(0.07);
   gr_ratio.GetYaxis() -> SetTitleSize(0.07);
   gr_ratio.GetYaxis() -> SetTitleOffset(0.7);
   gr_ratio.Draw("AP");
   TF1* f_const = new TF1("f_1", "[0]",gr_ratio.GetXaxis()->GetBinCenter(1)-gr_ratio.GetXaxis()->GetBinWidth(1)/2, gr_ratio.GetXaxis()->GetBinCenter(gr_ratio.GetXaxis()->GetNbins())+gr_ratio.GetXaxis()->GetBinWidth(gr_ratio.GetXaxis()->GetNbins())/2);
   f_const -> FixParameter(0,1);
   f_const -> SetLineColor(kRed);
   //f_const -> SetLineWidth(2);
   f_const -> Draw("same");

   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
}

void drawEfficiency(TEfficiency* eff_SuperCluster, TEfficiency* eff_DeepSuperCluster, std::string xtitle, std::string Name, std::string refLegend="Mustache", std::string valLegend="DeepSC")
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
   legend -> AddEntry(eff_SuperCluster,refLegend.c_str(),"L");
   legend -> AddEntry(eff_DeepSuperCluster,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();
   eff_SuperCluster->Draw("APL");
   eff_DeepSuperCluster->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	
}

void drawProfile(TProfile* prof_SuperCluster, TProfile* prof_DeepSuperCluster, std::string xtitle, std::string ytitle, std::string Name, std::string refLegend="Mustache", std::string valLegend="DeepSC")
{ 
   gStyle->SetOptStat(0000);  
   prof_SuperCluster->SetLineColor(kRed+1);
   prof_SuperCluster->SetLineWidth(2);
   prof_SuperCluster->GetXaxis()->SetTitle(xtitle.c_str()); 
   prof_SuperCluster->GetYaxis()->SetTitle(ytitle.c_str());   
   
   prof_DeepSuperCluster->SetLineColor(kBlue+1);
   prof_DeepSuperCluster->SetLineWidth(2);

   float min = prof_SuperCluster->GetMinimum();
   if(prof_DeepSuperCluster->GetMinimum()<min) min = prof_DeepSuperCluster->GetMinimum(); 
   float max = prof_SuperCluster->GetMaximum();
   if(prof_DeepSuperCluster->GetMaximum()>max) max = prof_DeepSuperCluster->GetMaximum(); 
   min = min*0.9;
   max = max*1.1;

   prof_SuperCluster->GetYaxis()->SetRangeUser(min,max);

   TLegend* legend = new TLegend(0.365, 0.72, 0.635, 0.90);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.03);
   legend -> AddEntry(prof_SuperCluster,refLegend.c_str(),"L");
   legend -> AddEntry(prof_DeepSuperCluster,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();
   prof_SuperCluster->Draw("PL");
   prof_DeepSuperCluster->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
}

//draw plots
void drawPlots(std::string fitFunction_, string superClusterRef_, string superClusterVal_)
{
   drawEfficiency(eff_SuperCluster_vs_EtaCalo, eff_DeepSuperCluster_vs_EtaCalo, std::string("caloParticle_#eta"), std::string("Efficiency_vs_CaloEta"), superClusterRef_, superClusterVal_); 
   drawEfficiency(eff_SuperCluster_vs_EtCalo_EB, eff_DeepSuperCluster_vs_EtCalo_EB, std::string("caloParticle_Et (GeV)"), std::string("Efficiency_vs_CaloEt_EB"), superClusterRef_, superClusterVal_); 
   drawEfficiency(eff_SuperCluster_vs_EtCalo_EE, eff_DeepSuperCluster_vs_EtCalo_EE, std::string("caloParticle_Et (GeV)"), std::string("Efficiency_vs_CaloEt_EE"), superClusterRef_, superClusterVal_); 
   
   drawProfile(prof_EoEtrue_vs_Eta_Calo_old, prof_EoEtrue_vs_Eta_Calo_new, std::string("caloParticle_#eta"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Eta_Seed_old, prof_EoEtrue_vs_Eta_Seed_new, std::string("seed_#eta"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_SeedEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Et_Calo_EB_old, prof_EoEtrue_vs_Et_Calo_EB_new, std::string("caloParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEt_EB"), superClusterRef_, superClusterVal_);  
   drawProfile(prof_EoEtrue_vs_Et_Seed_EB_old, prof_EoEtrue_vs_Et_Seed_EB_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_SeedEt_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Energy_Calo_EB_old, prof_EoEtrue_vs_Energy_Calo_EB_new, std::string("caloParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEnergy_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_nVtx_EB_old, prof_EoEtrue_vs_nVtx_EB_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_nVtx_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Rho_EB_old, prof_EoEtrue_vs_Rho_EB_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_Rho_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Et_Calo_EE_old, prof_EoEtrue_vs_Et_Calo_EE_new, std::string("caloParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Et_Seed_EE_old, prof_EoEtrue_vs_Et_Seed_EE_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_SeedEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Energy_Calo_EE_old, prof_EoEtrue_vs_Energy_Calo_EE_new, std::string("caloParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_CaloEnergy_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_nVtx_EE_old, prof_EoEtrue_vs_nVtx_EE_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_nVtx_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEtrue_vs_Rho_EE_old, prof_EoEtrue_vs_Rho_EE_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("Profile_EoEtrue_vs_Rho_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Eta_Gen_old, prof_EoEgen_vs_Eta_Gen_new, std::string("genParticle_#eta"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Eta_Seed_old, prof_EoEgen_vs_Eta_Seed_new, std::string("seed_#eta"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_SeedEta"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Gen_EB_old, prof_EoEgen_vs_Et_Gen_EB_new, std::string("genParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEt_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Seed_EB_old, prof_EoEgen_vs_Et_Seed_EB_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_SeedEt_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Energy_Gen_EB_old, prof_EoEgen_vs_Energy_Gen_EB_new, std::string("genParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEnergy_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_nVtx_EB_old, prof_EoEgen_vs_nVtx_EB_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_nVtx_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Rho_EB_old, prof_EoEgen_vs_Rho_EB_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_Rho_EB"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Gen_EE_old, prof_EoEgen_vs_Et_Gen_EE_new, std::string("genParticle_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Et_Seed_EE_old, prof_EoEgen_vs_Et_Seed_EE_new, std::string("seed_Et (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_SeedEt_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Energy_Gen_EE_old, prof_EoEgen_vs_Energy_Gen_EE_new, std::string("genParticle_Energy (GeV)"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_GenEnergy_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_nVtx_EE_old, prof_EoEgen_vs_nVtx_EE_new, std::string("nVtx"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_nVtx_EE"), superClusterRef_, superClusterVal_);
   drawProfile(prof_EoEgen_vs_Rho_EE_old, prof_EoEgen_vs_Rho_EE_new, std::string("#rho"), std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("Profile_EoEgen_vs_Rho_EE"), superClusterRef_, superClusterVal_);

   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_old= makeFitProfile(&EoEtrue_vs_Eta_Calo_old,-3.,3.,std::string("caloParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_new= makeFitProfile(&EoEtrue_vs_Eta_Calo_new,-3.,3.,std::string("caloParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_old= makeFitProfile(&EoEtrue_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_new= makeFitProfile(&EoEtrue_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_old= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_old,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_new= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_new,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_old= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_new= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_old= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_old,0.,250.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_new= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_new,0.,250.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_old= makeFitProfile(&EoEtrue_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_new= makeFitProfile(&EoEtrue_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_old= makeFitProfile(&EoEtrue_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_new= makeFitProfile(&EoEtrue_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_old= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_old,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_new= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_new,0.,100.,std::string("caloParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_old= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_new= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_old= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_old,0.,1000.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_new= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_new,0.,1000.,std::string("caloParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_old= makeFitProfile(&EoEtrue_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_new= makeFitProfile(&EoEtrue_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_old= makeFitProfile(&EoEtrue_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_new= makeFitProfile(&EoEtrue_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_old= makeFitProfile(&EoEgen_vs_Eta_Gen_old,-3.,3.,std::string("genParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_new= makeFitProfile(&EoEgen_vs_Eta_Gen_new,-3.,3.,std::string("genParticle_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_old= makeFitProfile(&EoEgen_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_new= makeFitProfile(&EoEgen_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_old= makeFitProfile(&EoEgen_vs_Et_Gen_EB_old,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_new= makeFitProfile(&EoEgen_vs_Et_Gen_EB_new,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_old= makeFitProfile(&EoEgen_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_new= makeFitProfile(&EoEgen_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_old= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_old,0.,250.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_new= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_new,0.,250.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_old= makeFitProfile(&EoEgen_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_new= makeFitProfile(&EoEgen_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_old= makeFitProfile(&EoEgen_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_new= makeFitProfile(&EoEgen_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_old= makeFitProfile(&EoEgen_vs_Et_Gen_EE_old,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_new= makeFitProfile(&EoEgen_vs_Et_Gen_EE_new,0.,100.,std::string("genParticle_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_old= makeFitProfile(&EoEgen_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_new= makeFitProfile(&EoEgen_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_old= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_old,0.,1000.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_new= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_new,0.,1000.,std::string("genParticle_Energy"), fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_old= makeFitProfile(&EoEgen_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_new= makeFitProfile(&EoEgen_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_old= makeFitProfile(&EoEgen_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_new= makeFitProfile(&EoEgen_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_);

   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_old_eff = makeFitProfile(&EoEtrue_vs_Eta_Calo_old,-3.,3.,std::string("caloParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Eta_new_eff = makeFitProfile(&EoEtrue_vs_Eta_Calo_new,-3.,3.,std::string("caloParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_old_eff = makeFitProfile(&EoEtrue_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEta_new_eff = makeFitProfile(&EoEtrue_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_old_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_old,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EB_new_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EB_new,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_old_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EB_new_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_,true); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_old_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_old,0.,250.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EB_new_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EB_new,0.,250.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_old_eff= makeFitProfile(&EoEtrue_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EB_new_eff= makeFitProfile(&EoEtrue_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_old_eff= makeFitProfile(&EoEtrue_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EB_new_eff= makeFitProfile(&EoEtrue_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_old_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_old,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Et_EE_new_eff= makeFitProfile(&EoEtrue_vs_Et_Calo_EE_new,0.,100.,std::string("caloParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_old_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_SeedEt_EE_new_eff= makeFitProfile(&EoEtrue_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_,true); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_old_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_old,0.,1000.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Energy_EE_new_eff= makeFitProfile(&EoEtrue_vs_Energy_Calo_EE_new,0.,1000.,std::string("caloParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_old_eff= makeFitProfile(&EoEtrue_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_nVtx_EE_new_eff= makeFitProfile(&EoEtrue_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_old_eff= makeFitProfile(&EoEtrue_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_Rho_EE_new_eff= makeFitProfile(&EoEtrue_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true); 
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_old_eff= makeFitProfile(&EoEgen_vs_Eta_Gen_old,-3.,3.,std::string("genParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Eta_new_eff= makeFitProfile(&EoEgen_vs_Eta_Gen_new,-3.,3.,std::string("genParticle_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_old_eff= makeFitProfile(&EoEgen_vs_Eta_Seed_old,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEta_new_eff= makeFitProfile(&EoEgen_vs_Eta_Seed_new,-3.,3.,std::string("seed_#eta"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_old_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EB_old,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EB_new_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EB_new,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_old_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EB_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EB_new_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EB_new,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_old_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_old,0.,250.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EB_new_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EB_new,0.,250.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_old_eff= makeFitProfile(&EoEgen_vs_nVtx_EB_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EB_new_eff= makeFitProfile(&EoEgen_vs_nVtx_EB_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_old_eff= makeFitProfile(&EoEgen_vs_Rho_EB_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EB_new_eff= makeFitProfile(&EoEgen_vs_Rho_EB_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_old_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EE_old,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Et_EE_new_eff= makeFitProfile(&EoEgen_vs_Et_Gen_EE_new,0.,100.,std::string("genParticle_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_old_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EE_old,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_SeedEt_EE_new_eff= makeFitProfile(&EoEgen_vs_Et_Seed_EE_new,0.,100.,std::string("seed_Et"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_old_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_old,0.,1000.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Energy_EE_new_eff= makeFitProfile(&EoEgen_vs_Energy_Gen_EE_new,0.,1000.,std::string("genParticle_Energy"), fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_old_eff= makeFitProfile(&EoEgen_vs_nVtx_EE_old,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_nVtx_EE_new_eff= makeFitProfile(&EoEgen_vs_nVtx_EE_new,0.,120.,std::string("nVtx"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_old_eff= makeFitProfile(&EoEgen_vs_Rho_EE_old,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEgen_vs_Rho_EE_new_eff= makeFitProfile(&EoEgen_vs_Rho_EE_new,0.,80.,std::string("#rho (GeV)"),fitFunction_,true);

   drawGraph(gr_EoEtrue_vs_Eta_old.first, gr_EoEtrue_vs_Eta_new.first, std::string("caloParticle_#eta"), std::string("#mu"), std::string("EoEtrue_vs_caloEta_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old.first, gr_EoEtrue_vs_SeedEta_new.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEtrue_vs_seedEta_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Eta_old.second, gr_EoEtrue_vs_Eta_new.second, std::string("caloParticle_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_caloEta_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old.second, gr_EoEtrue_vs_Eta_new.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_seedEta_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Et_EB_old.first, gr_EoEtrue_vs_Et_EB_new.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old.first, gr_EoEtrue_vs_SeedEt_EB_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EB_old.second, gr_EoEtrue_vs_Et_EB_new.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old.second, gr_EoEtrue_vs_SeedEt_EB_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Et_EE_old.first, gr_EoEtrue_vs_Et_EE_new.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old.first, gr_EoEtrue_vs_SeedEt_EE_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EE_old.second, gr_EoEtrue_vs_Et_EE_new.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old.second, gr_EoEtrue_vs_SeedEt_EE_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old.first, gr_EoEtrue_vs_Energy_EB_new.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old.second, gr_EoEtrue_vs_Energy_EB_new.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Energy_EE_old.first, gr_EoEtrue_vs_Energy_EE_new.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EE_old.second, gr_EoEtrue_vs_Energy_EE_new.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old.first, gr_EoEtrue_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old.first, gr_EoEtrue_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old.second, gr_EoEtrue_vs_nVtx_EB_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old.first, gr_EoEtrue_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old.first, gr_EoEtrue_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.95,1.03);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old.second, gr_EoEtrue_vs_nVtx_EE_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old.first, gr_EoEtrue_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old.first, gr_EoEtrue_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old.second, gr_EoEtrue_vs_Rho_EB_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old.first, gr_EoEtrue_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old.first, gr_EoEtrue_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.95,1.03);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old.second, gr_EoEtrue_vs_Rho_EE_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEgen_vs_Eta_old.first, gr_EoEgen_vs_Eta_new.first, std::string("genParticle_#eta"), std::string("#mu"), std::string("EoEgen_vs_genEta_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old.first, gr_EoEgen_vs_SeedEta_new.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEgen_vs_seedEta_Mean"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Eta_old.second, gr_EoEgen_vs_Eta_new.second, std::string("genParticle_#eta"), std::string("#sigma"), std::string("EoEgen_vs_genEta_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old.second, gr_EoEgen_vs_SeedEta_new.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEgen_vs_seedEta_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EB_old.first, gr_EoEgen_vs_Et_EB_new.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old.first, gr_EoEgen_vs_Et_EB_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EB_old.second, gr_EoEgen_vs_Et_EB_new.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old.second, gr_EoEgen_vs_SeedEt_EB_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Et_EE_old.first, gr_EoEgen_vs_Et_EE_new.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old.first, gr_EoEgen_vs_SeedEt_EE_new.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EE_old.second, gr_EoEgen_vs_Et_EE_new.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old.second, gr_EoEgen_vs_SeedEt_EE_new.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EE_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Energy_EB_old.first, gr_EoEgen_vs_Energy_EB_new.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EB_old.second, gr_EoEgen_vs_Energy_EB_new.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EB_Resolution"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Energy_EE_old.first, gr_EoEgen_vs_Energy_EE_new.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EE_old.second, gr_EoEgen_vs_Energy_EE_new.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEgen_vs_nVtx_EB_old.first, gr_EoEgen_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EB_old.first, gr_EoEgen_vs_nVtx_EB_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEgen_vs_nVtx_EB_old.second, gr_EoEgen_vs_nVtx_EB_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old.first, gr_EoEgen_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old.first, gr_EoEgen_vs_nVtx_EE_new.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.93,1.);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old.second, gr_EoEgen_vs_nVtx_EE_new.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EE_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old.first, gr_EoEgen_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EB_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old.first, gr_EoEgen_vs_Rho_EB_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EB_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.97,1.03);
   drawGraph(gr_EoEgen_vs_Rho_EB_old.second, gr_EoEgen_vs_Rho_EB_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EB_Resolution"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old.first, gr_EoEgen_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EE_Mean"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old.first, gr_EoEgen_vs_Rho_EE_new.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EE_Mean_Zoomed"), superClusterRef_, superClusterVal_,0.93,1.0);
   drawGraph(gr_EoEgen_vs_Rho_EE_old.second, gr_EoEgen_vs_Rho_EE_new.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EE_Resolution"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_Eta_old_eff.first, gr_EoEtrue_vs_Eta_new_eff.first, std::string("caloParticle_#eta"), std::string("#mu"), std::string("EoEtrue_vs_caloEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old_eff.first, gr_EoEtrue_vs_SeedEta_new_eff.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEtrue_vs_seedEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Eta_old_eff.second, gr_EoEtrue_vs_Eta_new_eff.second, std::string("caloParticle_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_caloEta_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEta_old_eff.second, gr_EoEtrue_vs_SeedEta_new_eff.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEtrue_vs_seedEta_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Et_EB_old_eff.first, gr_EoEtrue_vs_Et_EB_new_eff.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old_eff.first, gr_EoEtrue_vs_SeedEt_EB_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EB_old_eff.second, gr_EoEtrue_vs_Et_EB_new_eff.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EB_old_eff.second, gr_EoEtrue_vs_SeedEt_EB_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_Et_EE_old_eff.first, gr_EoEtrue_vs_Et_EE_new_eff.first, std::string("caloParticle_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old_eff.first, gr_EoEtrue_vs_SeedEt_EE_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_seedEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Et_EE_old_eff.second, gr_EoEtrue_vs_Et_EE_new_eff.second, std::string("caloParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_SeedEt_EE_old_eff.second, gr_EoEtrue_vs_SeedEt_EE_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_seedEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old_eff.first, gr_EoEtrue_vs_Energy_EB_new_eff.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EB_old_eff.second, gr_EoEtrue_vs_Energy_EB_new_eff.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EB_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEtrue_vs_Energy_EE_old_eff.first, gr_EoEtrue_vs_Energy_EE_new_eff.first, std::string("caloParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_caloEnergy_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Energy_EE_old_eff.second, gr_EoEtrue_vs_Energy_EE_new_eff.second, std::string("caloParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_caloEnergy_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old_eff.first, gr_EoEtrue_vs_nVtx_EB_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EB_old_eff.second, gr_EoEtrue_vs_nVtx_EB_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old_eff.first, gr_EoEtrue_vs_nVtx_EE_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEtrue_vs_nVtx_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_nVtx_EE_old_eff.second, gr_EoEtrue_vs_nVtx_EE_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEtrue_vs_nVtx_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old_eff.first, gr_EoEtrue_vs_Rho_EB_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EB_old_eff.second, gr_EoEtrue_vs_Rho_EB_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old_eff.first, gr_EoEtrue_vs_Rho_EE_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEtrue_vs_Rho_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEtrue_vs_Rho_EE_old_eff.second, gr_EoEtrue_vs_Rho_EE_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEtrue_vs_Rho_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);  

   drawGraph(gr_EoEgen_vs_Eta_old_eff.first, gr_EoEgen_vs_Eta_new_eff.first, std::string("genParticle_#eta"), std::string("#mu"), std::string("EoEgen_vs_genEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old_eff.first, gr_EoEgen_vs_SeedEta_new_eff.first, std::string("seed_#eta"), std::string("#mu"), std::string("EoEgen_vs_seedEta_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Eta_old_eff.second, gr_EoEgen_vs_Eta_new_eff.second, std::string("genParticle_#eta"), std::string("#sigma"), std::string("EoEgen_vs_genEta_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEta_old_eff.second, gr_EoEgen_vs_SeedEta_new_eff.second, std::string("seed_#eta"), std::string("#sigma"), std::string("EoEgen_vs_seedEta_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Et_EB_old_eff.first, gr_EoEgen_vs_Et_EB_new_eff.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old_eff.first, gr_EoEgen_vs_SeedEt_EB_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EB_old_eff.second, gr_EoEgen_vs_Et_EB_new_eff.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EB_old_eff.second, gr_EoEgen_vs_SeedEt_EB_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EB_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Et_EE_old_eff.first, gr_EoEgen_vs_Et_EE_new_eff.first, std::string("genParticle_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old_eff.first, gr_EoEgen_vs_SeedEt_EE_new_eff.first, std::string("seed_Et (GeV)"), std::string("#mu"), std::string("EoEgen_vs_seedEt_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Et_EE_old_eff.second, gr_EoEgen_vs_Et_EE_new_eff.second, std::string("genParticle_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_SeedEt_EE_old_eff.second, gr_EoEgen_vs_SeedEt_EE_new_eff.second, std::string("seed_Et (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_seedEt_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EB_old_eff.first, gr_EoEgen_vs_Energy_EB_new_eff.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EB_old_eff.second, gr_EoEgen_vs_Energy_EB_new_eff.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EB_Resolution_Effective"), superClusterRef_, superClusterVal_); 
   drawGraph(gr_EoEgen_vs_Energy_EE_old_eff.first, gr_EoEgen_vs_Energy_EE_new_eff.first, std::string("genParticle_Energy (GeV)"), std::string("#mu"), std::string("EoEgen_vs_genEnergy_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Energy_EE_old_eff.second, gr_EoEgen_vs_Energy_EE_new_eff.second, std::string("genParticle_Energy (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_genEnergy_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);  
   drawGraph(gr_EoEgen_vs_nVtx_EB_old_eff.first, gr_EoEgen_vs_nVtx_EB_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EB_old_eff.second, gr_EoEgen_vs_nVtx_EB_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old_eff.first, gr_EoEgen_vs_nVtx_EE_new_eff.first, std::string("nVtx"), std::string("#mu"), std::string("EoEgen_vs_nVtx_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_nVtx_EE_old_eff.second, gr_EoEgen_vs_nVtx_EE_new_eff.second, std::string("nVtx"), std::string("#sigma"), std::string("EoEgen_vs_nVtx_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old_eff.first, gr_EoEgen_vs_Rho_EB_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EB_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EB_old_eff.second, gr_EoEgen_vs_Rho_EB_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EB_Resolution_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old_eff.first, gr_EoEgen_vs_Rho_EE_new_eff.first, std::string("#rho (GeV)"), std::string("#mu"), std::string("EoEgen_vs_Rho_EE_Mean_Effective"), superClusterRef_, superClusterVal_);
   drawGraph(gr_EoEgen_vs_Rho_EE_old_eff.second, gr_EoEgen_vs_Rho_EE_new_eff.second, std::string("#rho (GeV)"), std::string("#sigma"), std::string("EoEgen_vs_Rho_EE_Resolution_Effective"), superClusterRef_, superClusterVal_);

   drawHisto(h_EoEtrue_EB_old, h_EoEtrue_EB_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EB_old, h_EoEtrue_EB_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EE_old, h_EoEtrue_EE_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEtrue_EE_old, h_EoEtrue_EE_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE"), 1, true, fitFunction_);
   drawHisto(h_EoEtrue_EB_seedMatched_old, h_EoEtrue_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EB_seedMatched_old, h_EoEtrue_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EB_seedMatched"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEtrue_EE_seedMatched_old, h_EoEtrue_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEtrue_EE_seedMatched_old, h_EoEtrue_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{SIM}"), std::string("e"), std::string("h_SC_EoEtrue_EE_seedMatched"), 1, true, fitFunction_);
   drawHisto(h_EoEgen_EB_old, h_EoEgen_EB_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EB_old, h_EoEgen_EB_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EE_old, h_EoEgen_EE_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEgen_EE_old, h_EoEgen_EE_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE"), 1, true, fitFunction_);
   drawHisto(h_EoEgen_EB_seedMatched_old, h_EoEgen_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EB_seedMatched_old, h_EoEgen_EB_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EB_seedMatched"), 1, true, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EoEgen_EE_seedMatched_old, h_EoEgen_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE_seedMatched"), 0, true, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EoEgen_EE_seedMatched_old, h_EoEgen_EE_seedMatched_new, std::string("SC_E_{Reco}/SC_E_{GEN}"), std::string("e"), std::string("h_SC_EoEgen_EE_seedMatched"), 1, true, fitFunction_);
   drawHisto(h_Energy_EB_old, h_Energy_EB_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_old, h_Energy_EB_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_old, h_Energy_EE_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_old, h_Energy_EE_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_old, h_Eta_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_old, h_Eta_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_old, h_Phi_EB_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_old, h_Phi_EB_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_old, h_Phi_EE_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_old, h_Phi_EE_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_old, h_EtaWidth_EB_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_old, h_EtaWidth_EB_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_old, h_EtaWidth_EE_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_old, h_EtaWidth_EE_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_old, h_PhiWidth_EB_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_old, h_PhiWidth_EB_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_old, h_PhiWidth_EE_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_old, h_PhiWidth_EE_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_old, h_nPFClusters_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_old, h_nPFClusters_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_old, h_nPFClusters_EB_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_old, h_nPFClusters_EB_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_old, h_nPFClusters_EE_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_old, h_nPFClusters_EE_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_old, h_R9_EB_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_old, h_R9_EB_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_old, h_R9_EE_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_old, h_R9_EE_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_old, h_full5x5_R9_EB_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_old, h_full5x5_R9_EB_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_old, h_full5x5_R9_EE_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_old, h_full5x5_R9_EE_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_old, h_sigmaIetaIeta_EB_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_old, h_sigmaIetaIeta_EB_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_old, h_sigmaIetaIeta_EE_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_old, h_sigmaIetaIeta_EE_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_old, h_full5x5_sigmaIetaIeta_EB_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_old, h_full5x5_sigmaIetaIeta_EB_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_old, h_full5x5_sigmaIetaIeta_EE_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_old, h_full5x5_sigmaIetaIeta_EE_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_old, h_sigmaIetaIphi_EB_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_old, h_sigmaIetaIphi_EB_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_old, h_sigmaIetaIphi_EE_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_old, h_sigmaIetaIphi_EE_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_old, h_full5x5_sigmaIetaIphi_EB_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_old, h_full5x5_sigmaIetaIphi_EB_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_old, h_full5x5_sigmaIetaIphi_EE_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_old, h_full5x5_sigmaIetaIphi_EE_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_old, h_sigmaIphiIphi_EB_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_old, h_sigmaIphiIphi_EB_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_old, h_sigmaIphiIphi_EE_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_old, h_sigmaIphiIphi_EE_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_old, h_full5x5_sigmaIphiIphi_EB_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_old, h_full5x5_sigmaIphiIphi_EB_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_old, h_full5x5_sigmaIphiIphi_EE_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_old, h_full5x5_sigmaIphiIphi_EE_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);

   //make seedMatched plots
   drawHisto(h_Energy_EB_seedMatched_old, h_Energy_EB_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_seedMatched_old, h_Energy_EB_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_seedMatched_old, h_Energy_EE_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_seedMatched_old, h_Energy_EE_seedMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_seedMatched_old, h_Eta_seedMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_seedMatched_old, h_Eta_seedMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_seedMatched_old, h_Phi_EB_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_seedMatched_old, h_Phi_EB_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_seedMatched_old, h_Phi_EE_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_seedMatched_old, h_Phi_EE_seedMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_seedMatched_old, h_EtaWidth_EB_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_seedMatched_old, h_EtaWidth_EB_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_seedMatched_old, h_EtaWidth_EE_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_seedMatched_old, h_EtaWidth_EE_seedMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_seedMatched_old, h_PhiWidth_EB_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_seedMatched_old, h_PhiWidth_EB_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_seedMatched_old, h_PhiWidth_EE_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_seedMatched_old, h_PhiWidth_EE_seedMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_seedMatched_old, h_nPFClusters_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_seedMatched_old, h_nPFClusters_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_seedMatched_old, h_nPFClusters_EB_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_seedMatched_old, h_nPFClusters_EB_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_seedMatched_old, h_nPFClusters_EE_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_seedMatched_old, h_nPFClusters_EE_seedMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_seedMatched_old, h_R9_EB_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_seedMatched_old, h_R9_EB_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_seedMatched_old, h_R9_EE_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_seedMatched_old, h_R9_EE_seedMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_seedMatched_old, h_full5x5_R9_EB_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_seedMatched_old, h_full5x5_R9_EB_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_seedMatched_old, h_full5x5_R9_EE_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_seedMatched_old, h_full5x5_R9_EE_seedMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_seedMatched_old, h_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_seedMatched_old, h_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_seedMatched_old, h_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_seedMatched_old, h_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_seedMatched_old, h_full5x5_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_seedMatched_old, h_full5x5_sigmaIetaIeta_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_seedMatched_old, h_full5x5_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_seedMatched_old, h_full5x5_sigmaIetaIeta_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_seedMatched_old, h_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_seedMatched_old, h_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_seedMatched_old, h_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_seedMatched_old, h_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_seedMatched_old, h_full5x5_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_seedMatched_old, h_full5x5_sigmaIetaIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_seedMatched_old, h_full5x5_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_seedMatched_old, h_full5x5_sigmaIetaIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_seedMatched_old, h_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_seedMatched_old, h_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_seedMatched_old, h_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_seedMatched_old, h_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_seedMatched_old, h_full5x5_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_seedMatched_old, h_full5x5_sigmaIphiIphi_EB_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_seedMatched_old, h_full5x5_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_seedMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_seedMatched_old, h_full5x5_sigmaIphiIphi_EE_seedMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_seedMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);

   //make caloMatched plots
   drawHisto(h_Energy_EB_caloMatched_old, h_Energy_EB_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_caloMatched_old, h_Energy_EB_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloMatched_old, h_Energy_EE_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloMatched_old, h_Energy_EE_caloMatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloMatched_old, h_Eta_caloMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloMatched_old, h_Eta_caloMatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_caloMatched_old, h_Phi_EB_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_caloMatched_old, h_Phi_EB_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_caloMatched_old, h_Phi_EE_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_caloMatched_old, h_Phi_EE_caloMatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_caloMatched_old, h_EtaWidth_EB_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_caloMatched_old, h_EtaWidth_EB_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_caloMatched_old, h_EtaWidth_EE_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_caloMatched_old, h_EtaWidth_EE_caloMatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_caloMatched_old, h_PhiWidth_EB_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_caloMatched_old, h_PhiWidth_EB_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_caloMatched_old, h_PhiWidth_EE_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_caloMatched_old, h_PhiWidth_EE_caloMatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_caloMatched_old, h_nPFClusters_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_caloMatched_old, h_nPFClusters_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloMatched_old, h_nPFClusters_EB_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloMatched_old, h_nPFClusters_EB_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloMatched_old, h_nPFClusters_EE_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloMatched_old, h_nPFClusters_EE_caloMatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloMatched_old, h_R9_EB_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloMatched_old, h_R9_EB_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloMatched_old, h_R9_EE_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloMatched_old, h_R9_EE_caloMatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloMatched_old, h_full5x5_R9_EB_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloMatched_old, h_full5x5_R9_EB_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloMatched_old, h_full5x5_R9_EE_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloMatched_old, h_full5x5_R9_EE_caloMatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloMatched_old, h_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloMatched_old, h_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloMatched_old, h_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloMatched_old, h_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloMatched_old, h_full5x5_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloMatched_old, h_full5x5_sigmaIetaIeta_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloMatched_old, h_full5x5_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloMatched_old, h_full5x5_sigmaIetaIeta_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloMatched_old, h_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloMatched_old, h_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloMatched_old, h_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloMatched_old, h_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloMatched_old, h_full5x5_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloMatched_old, h_full5x5_sigmaIetaIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloMatched_old, h_full5x5_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloMatched_old, h_full5x5_sigmaIetaIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloMatched_old, h_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloMatched_old, h_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloMatched_old, h_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloMatched_old, h_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloMatched_old, h_full5x5_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloMatched_old, h_full5x5_sigmaIphiIphi_EB_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloMatched_old, h_full5x5_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloMatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloMatched_old, h_full5x5_sigmaIphiIphi_EE_caloMatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloMatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);

   //make caloUnmatched plots
   drawHisto(h_Energy_EB_caloUnmatched_old, h_Energy_EB_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EB_caloUnmatched_old, h_Energy_EB_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloUnmatched_old, h_Energy_EE_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Energy_EE_caloUnmatched_old, h_Energy_EE_caloUnmatched_new, std::string("SC_Energy (GeV)"), std::string("hist"), std::string("h_SC_Energy_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloUnmatched_old, h_Eta_caloUnmatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Eta_caloUnmatched_old, h_Eta_caloUnmatched_new, std::string("SC_#eta"), std::string("hist"), std::string("h_SC_Eta_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EB_caloUnmatched_old, h_Phi_EB_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_Phi_EB_caloUnmatched_old, h_Phi_EB_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_Phi_EE_caloUnmatched_old, h_Phi_EE_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_Phi_EE_caloUnmatched_old, h_Phi_EE_caloUnmatched_new, std::string("SC_#phi"), std::string("hist"), std::string("h_SC_Phi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EB_caloUnmatched_old, h_EtaWidth_EB_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_EtaWidth_EB_caloUnmatched_old, h_EtaWidth_EB_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_EtaWidth_EE_caloUnmatched_old, h_EtaWidth_EE_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_EtaWidth_EE_caloUnmatched_old, h_EtaWidth_EE_caloUnmatched_new, std::string("SC_#etaWidth"), std::string("hist"), std::string("h_SC_EtaWidth_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EB_caloUnmatched_old, h_PhiWidth_EB_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);
   drawHisto(h_PhiWidth_EB_caloUnmatched_old, h_PhiWidth_EB_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_PhiWidth_EE_caloUnmatched_old, h_PhiWidth_EE_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_PhiWidth_EE_caloUnmatched_old, h_PhiWidth_EE_caloUnmatched_new, std::string("SC_#phiWidth"), std::string("hist"), std::string("h_SC_PhiWidth_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_);  
   drawHisto(h_nPFClusters_caloUnmatched_old, h_nPFClusters_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_caloUnmatched_old, h_nPFClusters_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloUnmatched_old, h_nPFClusters_EB_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EB_caloUnmatched_old, h_nPFClusters_EB_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloUnmatched_old, h_nPFClusters_EE_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_nPFClusters_EE_caloUnmatched_old, h_nPFClusters_EE_caloUnmatched_new, std::string("SC_nPFClusters"), std::string("hist"), std::string("h_SC_nPFClusters_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloUnmatched_old, h_R9_EB_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EB_caloUnmatched_old, h_R9_EB_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloUnmatched_old, h_R9_EE_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_R9_EE_caloUnmatched_old, h_R9_EE_caloUnmatched_new, std::string("SC_R9"), std::string("hist"), std::string("h_SC_R9_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloUnmatched_old, h_full5x5_R9_EB_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EB_caloUnmatched_old, h_full5x5_R9_EB_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloUnmatched_old, h_full5x5_R9_EE_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_R9_EE_caloUnmatched_old, h_full5x5_R9_EE_caloUnmatched_new, std::string("SC_full5x5_R9"), std::string("hist"), std::string("h_SC_full5x5_R9_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloUnmatched_old, h_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EB_caloUnmatched_old, h_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloUnmatched_old, h_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIeta_EE_caloUnmatched_old, h_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_sigmaIetaIeta_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EB_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIeta_EE_caloUnmatched_old, h_full5x5_sigmaIetaIeta_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#eta"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIeta_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloUnmatched_old, h_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EB_caloUnmatched_old, h_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloUnmatched_old, h_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIetaIphi_EE_caloUnmatched_old, h_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_sigmaIetaIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EB_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIetaIphi_EE_caloUnmatched_old, h_full5x5_sigmaIetaIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#etai#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIetaIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloUnmatched_old, h_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EB_caloUnmatched_old, h_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloUnmatched_old, h_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_sigmaIphiIphi_EE_caloUnmatched_old, h_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_sigmaIphiIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EB_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EB_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EB_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloUnmatched"), 0, false, fitFunction_, superClusterRef_, superClusterVal_); 
   drawHisto(h_full5x5_sigmaIphiIphi_EE_caloUnmatched_old, h_full5x5_sigmaIphiIphi_EE_caloUnmatched_new, std::string("SC_full5x5_#sigmai#phii#phi"), std::string("hist"), std::string("h_SC_full5x5_sigmaIphiIphi_EE_caloUnmatched"), 1, false, fitFunction_, superClusterRef_, superClusterVal_); 
 
}

