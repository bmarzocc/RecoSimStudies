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

string fitFunction_;

//DEFINE HISTOGRAMS
TH1F* h_EoEtrue_EB_dR_genScore = new TH1F("h_EoEtrue_EB_dR_genScore","pfCLuster_EoEtrue_EB_dR_genScore",200,0.,2.);
TH1F* h_EoEtrue_EE_dR_genScore = new TH1F("h_EoEtrue_EE_dR_genScore","pfCluster_EoEtrue_EE_dR_genScore",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_01 = new TH1F("h_EoEtrue_EB_sim_fraction_old_01","pfCLuster_EoEtrue_EB_sim_fraction_old_01",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_01 = new TH1F("h_EoEtrue_EE_sim_fraction_old_01","pfCluster_EoEtrue_EE_sim_fraction_old_01",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_02 = new TH1F("h_EoEtrue_EB_sim_fraction_old_02","pfCLuster_EoEtrue_EB_sim_fraction_old_02",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_02 = new TH1F("h_EoEtrue_EE_sim_fraction_old_02","pfCluster_EoEtrue_EE_sim_fraction_old_02",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_03 = new TH1F("h_EoEtrue_EB_sim_fraction_old_03","pfCLuster_EoEtrue_EB_sim_fraction_old_03",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_03 = new TH1F("h_EoEtrue_EE_sim_fraction_old_03","pfCluster_EoEtrue_EE_sim_fraction_old_03",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_04 = new TH1F("h_EoEtrue_EB_sim_fraction_old_04","pfCLuster_EoEtrue_EB_sim_fraction_old_04",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_04 = new TH1F("h_EoEtrue_EE_sim_fraction_old_04","pfCluster_EoEtrue_EE_sim_fraction_old_04",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_05 = new TH1F("h_EoEtrue_EB_sim_fraction_old_05","pfCLuster_EoEtrue_EB_sim_fraction_old_05",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_05 = new TH1F("h_EoEtrue_EE_sim_fraction_old_05","pfCluster_EoEtrue_EE_sim_fraction_old_05",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_06 = new TH1F("h_EoEtrue_EB_sim_fraction_old_06","pfCLuster_EoEtrue_EB_sim_fraction_old_06",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_06 = new TH1F("h_EoEtrue_EE_sim_fraction_old_06","pfCluster_EoEtrue_EE_sim_fraction_old_06",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_07 = new TH1F("h_EoEtrue_EB_sim_fraction_old_07","pfCLuster_EoEtrue_EB_sim_fraction_old_07",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_07 = new TH1F("h_EoEtrue_EE_sim_fraction_old_07","pfCluster_EoEtrue_EE_sim_fraction_old_07",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_08 = new TH1F("h_EoEtrue_EB_sim_fraction_old_08","pfCLuster_EoEtrue_EB_sim_fraction_old_08",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_08 = new TH1F("h_EoEtrue_EE_sim_fraction_old_08","pfCluster_EoEtrue_EE_sim_fraction_old_08",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_09 = new TH1F("h_EoEtrue_EB_sim_fraction_old_09","pfCLuster_EoEtrue_EB_sim_fraction_old_09",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_09 = new TH1F("h_EoEtrue_EE_sim_fraction_old_09","pfCluster_EoEtrue_EE_sim_fraction_old_09",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_old_10 = new TH1F("h_EoEtrue_EB_sim_fraction_old_10","pfCLuster_EoEtrue_EB_sim_fraction_old_10",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_old_10 = new TH1F("h_EoEtrue_EE_sim_fraction_old_10","pfCluster_EoEtrue_EE_sim_fraction_old_10",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_001 = new TH1F("h_EoEtrue_EB_sim_fraction_001","pfCLuster_EoEtrue_EB_sim_fraction_001",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_001 = new TH1F("h_EoEtrue_EE_sim_fraction_001","pfCluster_EoEtrue_EE_sim_fraction_001",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_002 = new TH1F("h_EoEtrue_EB_sim_fraction_002","pfCLuster_EoEtrue_EB_sim_fraction_002",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_002 = new TH1F("h_EoEtrue_EE_sim_fraction_002","pfCluster_EoEtrue_EE_sim_fraction_002",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_003 = new TH1F("h_EoEtrue_EB_sim_fraction_003","pfCLuster_EoEtrue_EB_sim_fraction_003",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_003 = new TH1F("h_EoEtrue_EE_sim_fraction_003","pfCluster_EoEtrue_EE_sim_fraction_003",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_004 = new TH1F("h_EoEtrue_EB_sim_fraction_004","pfCLuster_EoEtrue_EB_sim_fraction_004",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_004 = new TH1F("h_EoEtrue_EE_sim_fraction_004","pfCluster_EoEtrue_EE_sim_fraction_004",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_005 = new TH1F("h_EoEtrue_EB_sim_fraction_005","pfCLuster_EoEtrue_EB_sim_fraction_005",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_005 = new TH1F("h_EoEtrue_EE_sim_fraction_005","pfCluster_EoEtrue_EE_sim_fraction_005",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_006 = new TH1F("h_EoEtrue_EB_sim_fraction_006","pfCLuster_EoEtrue_EB_sim_fraction_006",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_006 = new TH1F("h_EoEtrue_EE_sim_fraction_006","pfCluster_EoEtrue_EE_sim_fraction_006",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_007 = new TH1F("h_EoEtrue_EB_sim_fraction_007","pfCLuster_EoEtrue_EB_sim_fraction_007",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_007 = new TH1F("h_EoEtrue_EE_sim_fraction_007","pfCluster_EoEtrue_EE_sim_fraction_007",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_008 = new TH1F("h_EoEtrue_EB_sim_fraction_008","pfCLuster_EoEtrue_EB_sim_fraction_008",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_008 = new TH1F("h_EoEtrue_EE_sim_fraction_008","pfCluster_EoEtrue_EE_sim_fraction_008",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_009 = new TH1F("h_EoEtrue_EB_sim_fraction_009","pfCLuster_EoEtrue_EB_sim_fraction_009",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_009 = new TH1F("h_EoEtrue_EE_sim_fraction_009","pfCluster_EoEtrue_EE_sim_fraction_009",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_01 = new TH1F("h_EoEtrue_EB_sim_fraction_01","pfCLuster_EoEtrue_EB_sim_fraction_01",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_01 = new TH1F("h_EoEtrue_EE_sim_fraction_01","pfCluster_EoEtrue_EE_sim_fraction_01",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_02 = new TH1F("h_EoEtrue_EB_sim_fraction_02","pfCLuster_EoEtrue_EB_sim_fraction_02",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_02 = new TH1F("h_EoEtrue_EE_sim_fraction_02","pfCluster_EoEtrue_EE_sim_fraction_02",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_03 = new TH1F("h_EoEtrue_EB_sim_fraction_03","pfCLuster_EoEtrue_EB_sim_fraction_03",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_03 = new TH1F("h_EoEtrue_EE_sim_fraction_03","pfCluster_EoEtrue_EE_sim_fraction_03",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_04 = new TH1F("h_EoEtrue_EB_sim_fraction_04","pfCLuster_EoEtrue_EB_sim_fraction_04",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_04 = new TH1F("h_EoEtrue_EE_sim_fraction_04","pfCluster_EoEtrue_EE_sim_fraction_04",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_05 = new TH1F("h_EoEtrue_EB_sim_fraction_05","pfCLuster_EoEtrue_EB_sim_fraction_05",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_05 = new TH1F("h_EoEtrue_EE_sim_fraction_05","pfCluster_EoEtrue_EE_sim_fraction_05",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_01 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_01","pfCLuster_EoEtrue_EB_sim_rechit_diff_01",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_01 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_01","pfCluster_EoEtrue_EE_sim_rechit_diff_01",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_02 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_02","pfCLuster_EoEtrue_EB_sim_rechit_diff_02",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_02 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_02","pfCluster_EoEtrue_EE_sim_rechit_diff_02",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_03 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_03","pfCLuster_EoEtrue_EB_sim_rechit_diff_03",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_03 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_03","pfCluster_EoEtrue_EE_sim_rechit_diff_03",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_04 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_04","pfCLuster_EoEtrue_EB_sim_rechit_diff_04",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_04 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_04","pfCluster_EoEtrue_EE_sim_rechit_diff_04",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_05 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_05","pfCLuster_EoEtrue_EB_sim_rechit_diff_05",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_05 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_05","pfCluster_EoEtrue_EE_sim_rechit_diff_05",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_06 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_06","pfCLuster_EoEtrue_EB_sim_rechit_diff_06",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_06 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_06","pfCluster_EoEtrue_EE_sim_rechit_diff_06",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_07 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_07","pfCLuster_EoEtrue_EB_sim_rechit_diff_07",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_07 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_07","pfCluster_EoEtrue_EE_sim_rechit_diff_07",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_08 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_08","pfCLuster_EoEtrue_EB_sim_rechit_diff_08",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_08 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_08","pfCluster_EoEtrue_EE_sim_rechit_diff_08",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_09 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_09","pfCLuster_EoEtrue_EB_sim_rechit_diff_09",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_09 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_09","pfCluster_EoEtrue_EE_sim_rechit_diff_09",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_diff_10 = new TH1F("h_EoEtrue_EB_sim_rechit_diff_10","pfCLuster_EoEtrue_EB_sim_rechit_diff_10",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_diff_10 = new TH1F("h_EoEtrue_EE_sim_rechit_diff_10","pfCluster_EoEtrue_EE_sim_rechit_diff_10",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_01 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_01","pfCLuster_EoEtrue_EB_sim_rechit_fraction_01",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_01 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_01","pfCluster_EoEtrue_EE_sim_rechit_fraction_01",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_02 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_02","pfCLuster_EoEtrue_EB_sim_rechit_fraction_02",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_02 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_02","pfCluster_EoEtrue_EE_sim_rechit_fraction_02",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_03 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_03","pfCLuster_EoEtrue_EB_sim_rechit_fraction_03",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_03 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_03","pfCluster_EoEtrue_EE_sim_rechit_fraction_03",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_04 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_04","pfCLuster_EoEtrue_EB_sim_rechit_fraction_04",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_04 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_04","pfCluster_EoEtrue_EE_sim_rechit_fraction_04",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_05 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_05","pfCLuster_EoEtrue_EB_sim_rechit_fraction_05",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_05 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_05","pfCluster_EoEtrue_EE_sim_rechit_fraction_05",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_06 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_06","pfCLuster_EoEtrue_EB_sim_rechit_fraction_06",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_06 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_06","pfCluster_EoEtrue_EE_sim_rechit_fraction_06",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_07 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_07","pfCLuster_EoEtrue_EB_sim_rechit_fraction_07",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_07 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_07","pfCluster_EoEtrue_EE_sim_rechit_fraction_07",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_08 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_08","pfCLuster_EoEtrue_EB_sim_rechit_fraction_08",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_08 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_08","pfCluster_EoEtrue_EE_sim_rechit_fraction_08",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_09 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_09","pfCLuster_EoEtrue_EB_sim_rechit_fraction_09",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_09 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_09","pfCluster_EoEtrue_EE_sim_rechit_fraction_09",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_fraction_10 = new TH1F("h_EoEtrue_EB_sim_rechit_fraction_10","pfCLuster_EoEtrue_EB_sim_rechit_fraction_10",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_fraction_10 = new TH1F("h_EoEtrue_EE_sim_rechit_fraction_10","pfCluster_EoEtrue_EE_sim_rechit_fraction_10",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_01 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_01","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_01",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_01 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_01","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_01",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_02 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_02","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_02",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_02 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_02","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_02",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_03 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_03","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_03",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_03 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_03","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_03",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_04 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_04","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_04",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_04 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_04","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_04",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_05 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_05","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_05",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_05 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_05","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_05",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_06 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_06","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_06",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_06 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_06","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_06",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_07 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_07","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_07",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_07 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_07","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_07",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_08 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_08","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_08",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_08 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_08","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_08",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_09 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_09","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_09",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_09 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_09","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_09",200,0.,2.);
TH1F* h_EoEtrue_EB_global_sim_rechit_fraction_10 = new TH1F("h_EoEtrue_EB_global_sim_rechit_fraction_10","pfCLuster_EoEtrue_EB_global_sim_rechit_fraction_10",200,0.,2.);
TH1F* h_EoEtrue_EE_global_sim_rechit_fraction_10 = new TH1F("h_EoEtrue_EE_global_sim_rechit_fraction_10","pfCluster_EoEtrue_EE_global_sim_rechit_fraction_10",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_combined_fraction01 = new TH1F("h_EoEtrue_EB_sim_rechit_combined_fraction01","pfCLuster_EoEtrue_EB_sim_rechit_combined_fraction01",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_combined_fraction01 = new TH1F("h_EoEtrue_EE_sim_rechit_combined_fraction01","pfCluster_EoEtrue_EE_sim_rechit_combined_fraction01",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_combined_fraction02 = new TH1F("h_EoEtrue_EB_sim_rechit_combined_fraction02","pfCLuster_EoEtrue_EB_sim_rechit_combined_fraction02",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_combined_fraction02 = new TH1F("h_EoEtrue_EE_sim_rechit_combined_fraction02","pfCluster_EoEtrue_EE_sim_rechit_combined_fraction02",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_combined_fraction03 = new TH1F("h_EoEtrue_EB_sim_rechit_combined_fraction03","pfCLuster_EoEtrue_EB_sim_rechit_combined_fraction03",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_combined_fraction03 = new TH1F("h_EoEtrue_EE_sim_rechit_combined_fraction03","pfCluster_EoEtrue_EE_sim_rechit_combined_fraction03",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_combined_fraction04 = new TH1F("h_EoEtrue_EB_sim_rechit_combined_fraction04","pfCLuster_EoEtrue_EB_sim_rechit_combined_fraction04",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_combined_fraction04 = new TH1F("h_EoEtrue_EE_sim_rechit_combined_fraction04","pfCluster_EoEtrue_EE_sim_rechit_combined_fraction04",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_rechit_combined_fraction05 = new TH1F("h_EoEtrue_EB_sim_rechit_combined_fraction05","pfCLuster_EoEtrue_EB_sim_rechit_combined_fraction05",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_rechit_combined_fraction05 = new TH1F("h_EoEtrue_EE_sim_rechit_combined_fraction05","pfCluster_EoEtrue_EE_sim_rechit_combined_fraction05",200,0.,2.);
TH1F* h_EoEtrue_EB_hgcal_clusterToCalo = new TH1F("h_EoEtrue_EB_hgcal_clusterToCalo","pfCLuster_EoEtrue_EB_hgcal_clusterToCalo",200,0.,2.);
TH1F* h_EoEtrue_EE_hgcal_clusterToCalo = new TH1F("h_EoEtrue_EE_hgcal_clusterToCalo","pfCluster_EoEtrue_EE_hgcal_clusterToCalo",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_1MeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_1MeVCut","pfCLuster_EoEtrue_EB_sim_fraction_1MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_1MeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_1MeVCut","pfCluster_EoEtrue_EE_sim_fraction_1MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_5MeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_5MeVCut","pfCLuster_EoEtrue_EB_sim_fraction_5MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_5MeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_5MeVCut","pfCluster_EoEtrue_EE_sim_fraction_5MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_10MeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_10MeVCut","pfCLuster_EoEtrue_EB_sim_fraction_10MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_10MeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_10MeVCut","pfCluster_EoEtrue_EE_sim_fraction_10MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_50MeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_50MeVCut","pfCLuster_EoEtrue_EB_sim_fraction_50MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_50MeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_50MeVCut","pfCluster_EoEtrue_EE_sim_fraction_50MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_100MeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_100MeVCut","pfCLuster_EoEtrue_EB_sim_fraction_100MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_100MeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_100MeVCut","pfCluster_EoEtrue_EE_sim_fraction_100MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_500MeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_500MeVCut","pfCLuster_EoEtrue_EB_sim_fraction_500MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_500MeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_500MeVCut","pfCluster_EoEtrue_EE_sim_fraction_500MeVCut",200,0.,2.);
TH1F* h_EoEtrue_EB_sim_fraction_1GeVCut = new TH1F("h_EoEtrue_EB_sim_fraction_1GeVCut","pfCLuster_EoEtrue_EB_sim_fraction_1GeVCut",200,0.,2.);
TH1F* h_EoEtrue_EE_sim_fraction_1GeVCut = new TH1F("h_EoEtrue_EE_sim_fraction_1GeVCut","pfCluster_EoEtrue_EE_sim_fraction_1GeVCut",200,0.,2.);

TH1F* h_EoEtrue_EB_simScore_final_combination = new TH1F("h_EoEtrue_EB_simScore_final_combination","pfCLuster_EoEtrue_EB_simScore_final_combination",200,0.,2.);
TH1F* h_EoEtrue_EE_simScore_final_combination = new TH1F("h_EoEtrue_EE_simScore_final_combination","pfCluster_EoEtrue_EE_simScore_final_combination",200,0.,2.);

//DEFINE BRANCHES
vector<float> *genParticle_energy;
vector<float> *genParticle_eta;
vector<float> *genParticle_phi;
vector<float> *caloParticle_simEnergy;
vector<float> *caloParticle_simEta;
vector<float> *caloParticle_simPhi;
vector<float> *pfCluster_energy;
vector<float> *pfCluster_eta;
vector<float> *pfCluster_phi;
vector<vector<double> > *pfCluster_dR_genScore;
vector<vector<double> > *pfCluster_sim_fraction_old;
vector<vector<int> > *pfCluster_n_shared_xtals;
vector<vector<double> > *pfCluster_sim_fraction;
vector<vector<double> > *pfCluster_sim_fraction_1MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_5MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_10MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_50MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_100MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_500MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_1GeVCut;
vector<vector<double> > *pfCluster_sim_fraction_min1;
vector<vector<double> > *pfCluster_sim_fraction_min3;
vector<vector<double> > *pfCluster_sim_fraction_min3_1MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_min3_5MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_min3_10MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_min3_50MeVCut;
vector<vector<double> > *pfCluster_sim_fraction_min3_100MeVCut;
vector<vector<double> > *pfCluster_sim_rechit_diff;
vector<vector<double> > *pfCluster_sim_rechit_fraction;
vector<vector<double> > *pfCluster_global_sim_rechit_fraction;
vector<vector<double> > *pfCluster_hgcal_caloToCluster;
vector<vector<double> > *pfCluster_hgcal_clusterToCalo;
vector<vector<double> > *pfCluster_sim_rechit_combined_fraction;
vector<vector<double> > *pfCluster_rechit_sim_combined_fraction;

TBranch *b_genParticle_energy;   //!
TBranch *b_genParticle_eta;   //!
TBranch *b_genParticle_phi;   //!
TBranch *b_caloParticle_simEnergy;   //!
TBranch *b_caloParticle_simEta;   //!
TBranch *b_caloParticle_simPhi;   //!
TBranch *b_pfCluster_energy;   //!
TBranch *b_pfCluster_eta;   //!
TBranch *b_pfCluster_phi;   //!
TBranch *b_pfCluster_dR_genScore;   //!
TBranch *b_pfCluster_sim_fraction_old;   //!
TBranch *b_pfCluster_n_shared_xtals;   //!
TBranch *b_pfCluster_sim_fraction;   //!
TBranch *b_pfCluster_sim_fraction_1MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_5MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_10MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_50MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_100MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_500MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_1GeVCut;   //!
TBranch *b_pfCluster_sim_fraction_min1;   //!
TBranch *b_pfCluster_sim_fraction_min3;   //!
TBranch *b_pfCluster_sim_fraction_min3_1MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_min3_5MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_min3_10MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_min3_50MeVCut;   //!
TBranch *b_pfCluster_sim_fraction_min3_100MeVCut;   //!
TBranch *b_pfCluster_sim_rechit_diff;   //!
TBranch *b_pfCluster_sim_rechit_fraction;   //!
TBranch *b_pfCluster_global_sim_rechit_fraction;   //!
TBranch *b_pfCluster_hgcal_caloToCluster;   //!
TBranch *b_pfCluster_hgcal_clusterToCalo;   //!
TBranch *b_pfCluster_sim_rechit_combined_fraction;   //!
TBranch *b_pfCluster_rechit_sim_combined_fraction;   //!

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

//setTreeBranches
void setTreeBranches(TTree* tree)
{
   tree->SetBranchAddress("genParticle_energy", &genParticle_energy, &b_genParticle_energy);
   tree->SetBranchAddress("genParticle_eta", &genParticle_eta, &b_genParticle_eta);
   tree->SetBranchAddress("genParticle_phi", &genParticle_phi, &b_genParticle_phi);
   tree->SetBranchAddress("caloParticle_simEnergy", &caloParticle_simEnergy, &b_caloParticle_simEnergy);
   tree->SetBranchAddress("caloParticle_simEta", &caloParticle_simEta, &b_caloParticle_simEta);
   tree->SetBranchAddress("caloParticle_simPhi", &caloParticle_simPhi, &b_caloParticle_simPhi);
   tree->SetBranchAddress("pfCluster_energy", &pfCluster_energy, &b_pfCluster_energy);
   tree->SetBranchAddress("pfCluster_eta", &pfCluster_eta, &b_pfCluster_eta);
   tree->SetBranchAddress("pfCluster_phi", &pfCluster_phi, &b_pfCluster_phi);
   tree->SetBranchAddress("pfCluster_dR_genScore", &pfCluster_dR_genScore, &b_pfCluster_dR_genScore);
   tree->SetBranchAddress("pfCluster_sim_fraction_old", &pfCluster_sim_fraction_old, &b_pfCluster_sim_fraction_old);
   tree->SetBranchAddress("pfCluster_n_shared_xtals", &pfCluster_n_shared_xtals, &b_pfCluster_n_shared_xtals);
   tree->SetBranchAddress("pfCluster_sim_fraction", &pfCluster_sim_fraction, &b_pfCluster_sim_fraction);
   tree->SetBranchAddress("pfCluster_sim_fraction_1MeVCut", &pfCluster_sim_fraction_1MeVCut, &b_pfCluster_sim_fraction_1MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_5MeVCut", &pfCluster_sim_fraction_5MeVCut, &b_pfCluster_sim_fraction_5MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_10MeVCut", &pfCluster_sim_fraction_10MeVCut, &b_pfCluster_sim_fraction_10MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_50MeVCut", &pfCluster_sim_fraction_50MeVCut, &b_pfCluster_sim_fraction_50MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_100MeVCut", &pfCluster_sim_fraction_100MeVCut, &b_pfCluster_sim_fraction_100MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_500MeVCut", &pfCluster_sim_fraction_500MeVCut, &b_pfCluster_sim_fraction_500MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_1GeVCut", &pfCluster_sim_fraction_1GeVCut, &b_pfCluster_sim_fraction_1GeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_min1", &pfCluster_sim_fraction_min1, &b_pfCluster_sim_fraction_min1);
   tree->SetBranchAddress("pfCluster_sim_fraction_min3", &pfCluster_sim_fraction_min3, &b_pfCluster_sim_fraction_min3);
   tree->SetBranchAddress("pfCluster_sim_fraction_min3_1MeVCut", &pfCluster_sim_fraction_min3_1MeVCut, &b_pfCluster_sim_fraction_min3_1MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_min3_5MeVCut", &pfCluster_sim_fraction_min3_5MeVCut, &b_pfCluster_sim_fraction_min3_5MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_min3_10MeVCut", &pfCluster_sim_fraction_min3_10MeVCut, &b_pfCluster_sim_fraction_min3_10MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_min3_50MeVCut", &pfCluster_sim_fraction_min3_50MeVCut, &b_pfCluster_sim_fraction_min3_50MeVCut);
   tree->SetBranchAddress("pfCluster_sim_fraction_min3_100MeVCut", &pfCluster_sim_fraction_min3_100MeVCut, &b_pfCluster_sim_fraction_min3_100MeVCut);
   tree->SetBranchAddress("pfCluster_sim_rechit_diff", &pfCluster_sim_rechit_diff, &b_pfCluster_sim_rechit_diff);
   tree->SetBranchAddress("pfCluster_sim_rechit_fraction", &pfCluster_sim_rechit_fraction, &b_pfCluster_sim_rechit_fraction);
   tree->SetBranchAddress("pfCluster_global_sim_rechit_fraction", &pfCluster_global_sim_rechit_fraction, &b_pfCluster_global_sim_rechit_fraction);
   tree->SetBranchAddress("pfCluster_hgcal_caloToCluster", &pfCluster_hgcal_caloToCluster, &b_pfCluster_hgcal_caloToCluster);
   tree->SetBranchAddress("pfCluster_hgcal_clusterToCalo", &pfCluster_hgcal_clusterToCalo, &b_pfCluster_hgcal_clusterToCalo);
   tree->SetBranchAddress("pfCluster_sim_rechit_combined_fraction", &pfCluster_sim_rechit_combined_fraction, &b_pfCluster_sim_rechit_combined_fraction);
   tree->SetBranchAddress("pfCluster_rechit_sim_combined_fraction", &pfCluster_rechit_sim_combined_fraction, &b_pfCluster_rechit_sim_combined_fraction);
}

int getMatchedIndex(vector<vector<double>>* score, double selection, bool useMax, vector<vector<std::vector<double>>>* scoreSelMax, vector<double>* selectionMax, vector<vector<std::vector<double>>>* scoreSelMin, std::vector<double>* selectionMin, int iPF)
{
   int matchedIndex = -1; 
   if(!useMax){ 
      std::replace(score->at(iPF).begin(),score->at(iPF).end(), -999., 999.);
      if(std::all_of(score->at(iPF).begin(),score->at(iPF).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iPF).begin(),score->at(iPF).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::min_element(score->at(iPF).begin(),score->at(iPF).end()) - score->at(iPF).begin();
   }else{
      std::replace(score->at(iPF).begin(),score->at(iPF).end(), 999., -999.); 
      if(std::all_of(score->at(iPF).begin(),score->at(iPF).end(),[](double i){return i==-999.;}) || std::all_of(score->at(iPF).begin(),score->at(iPF).end(),[](double i){return i==999.;})) matchedIndex=-1;
      else matchedIndex = std::max_element(score->at(iPF).begin(),score->at(iPF).end()) - score->at(iPF).begin();
   }

   if(matchedIndex==-1) return -1;
   
   //std::cout << "getMatchedIndex - Max - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMax->size() << " - " << selectionMax->size() << std::endl;
   //std::cout << "getMatchedIndex - Min - iPF = " << iPF << " - " << matchedIndex << " - " << scoreSelMin->size() << " - " << selectionMin->size() << std::endl;
    
   bool passSelection = true;
   for(unsigned int iSelMax=0; iSelMax < scoreSelMax->size(); iSelMax++)
       if(scoreSelMax->at(iSelMax).at(iPF).at(matchedIndex) < selectionMax->at(iSelMax)) passSelection = false;
   for(unsigned int iSelMin=0; iSelMin < scoreSelMin->size(); iSelMin++)
       if(scoreSelMin->at(iSelMin).at(iPF).at(matchedIndex) > selectionMin->at(iSelMin)) passSelection = false;

   //if(!passSelection) return -1;
   if(useMax && score->at(iPF).at(matchedIndex) > selection && passSelection) return matchedIndex;
   if(!useMax && score->at(iPF).at(matchedIndex) < selection && passSelection) return matchedIndex; 
   
   return -1; 
}

void fillHisto(vector<vector<double>>* pfCluster_simScore, double selection, bool useMax, vector<vector<std::vector<double>>> scoreSelMax, vector<double> selectionMax,  vector<vector<std::vector<double>>> scoreSelMin, vector<double> selectionMin, vector<float> *pfCluster_energy, vector<float>* caloParticle_simEnergy, vector<float>* caloParticle_simEta, TH1F* h_EoEtrue, bool isEB)
{
   vector<vector<std::vector<double>>>* scoreSelMax_tmp = &scoreSelMax;
   vector<vector<std::vector<double>>>* scoreSelMin_tmp = &scoreSelMin;
   vector<double>* selectionMax_tmp = &selectionMax;
   vector<double>* selectionMin_tmp = &selectionMin;

   //std::cout << "fillHisto: Max - " << scoreSelMax_tmp->size() << " - " << selectionMax_tmp->size() << std::endl; 
   //std::cout << "fillHisto: Min - " << scoreSelMin_tmp->size() << " - " << selectionMin_tmp->size() << std::endl; 

   std::map<int,float> recoEnergy;
   for(unsigned int iPF=0; iPF<pfCluster_simScore->size(); iPF++)
   { 
       int matchedIndex = getMatchedIndex(pfCluster_simScore, selection, useMax, scoreSelMax_tmp, selectionMax_tmp, scoreSelMin_tmp, selectionMin_tmp, iPF);
       if(matchedIndex>=0) recoEnergy[matchedIndex] += pfCluster_energy->at(iPF);       
   }
   for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy->size(); iCalo++)
   { 
       if(fabs(caloParticle_simEta->at(iCalo))<1.479 && isEB) h_EoEtrue->Fill(recoEnergy[iCalo]/caloParticle_simEnergy->at(iCalo));
       if(fabs(caloParticle_simEta->at(iCalo))>1.479 && !isEB) h_EoEtrue->Fill(recoEnergy[iCalo]/caloParticle_simEnergy->at(iCalo)); 
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

TF1* makeCruijffFit(TH1* hist,float xmin,float xmax)
{
  //hist->Scale(1./hist->GetEntries()); 
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;
  float nrEntries = hist->Integral();
  RooRealVar  nsig("N_{S}", "#signal events", nrEntries, nrEntries*0.01, nrEntries*100.);
  RooRealVar mean( "#DeltaE", "mean_{cb}", 1. ,hist->GetMean()-1.,hist->GetMean()+1.,"");
  RooRealVar sigmaL("#sigma_{L}","#sigma_{L}", hist->GetRMS()/2., 0.00001, 1.0);
  RooRealVar sigmaR("#sigma_{R}","#sigma_{R}", hist->GetRMS()/2., 0.00001, 1.0);
  RooRealVar alphaL( "alpha_{L}", "alpha_{L}", 0.1,0.,20.);
  RooRealVar alphaR( "alpha_{R}", "alpha_{R}", 0.1,0.,20.);
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

  hist->GetXaxis()->SetRangeUser(0.2,1.8); 
  double maximum = hist->GetMaximum();
  hist->GetXaxis()->SetRangeUser(0.,2.); 
  TF1* Cruijff = new TF1("Cruijff",&mycruijff,xmin,xmax,6);
  Cruijff->SetParameters(mean.getVal(), sigmaL.getVal(), sigmaR.getVal(), alphaL.getVal(), alphaR.getVal(), maximum);
  Cruijff->SetParErrors(errors);

  return Cruijff;
}

double mydoublecb(double* x, double* par) 
{
  //a priori we allow for different shape of right and left tail, thus two values of alpha and n 
  Double_t xcur = x[0];
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t alphaL = par[2];
  Double_t nL = par[3];
  Double_t alphaR = par[4];
  Double_t nR = par[5];
  Double_t N = par[6];
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

TF1* makeDoubleCBFit(TH1* hist,float xmin,float xmax)
{
  RooRealVar  res("res","E^{reco}/E^{gen}", xmin,xmax,"");
  res.setBins(10000,"cache") ;
  res.setMin("cache",xmin) ;
  res.setMax("cache",xmax) ;

  float nrEntries = hist->Integral();
  RooRealVar nsig("N_{S}", "#signal events", nrEntries,nrEntries*0.5,nrEntries*2.);
  RooRealVar mean( "#DeltaE", "mean_{cb}",hist->GetMean()-1.,hist->GetMean()+1.,""); 
  RooRealVar cbSigma("#sigma_{CB}","CB Width", hist->GetRMS()/2., 0.0001, 0.5,"");
  RooRealVar alpha1( "alpha_{1}", "alpha_{1}", 0.1,0.,20.);
  //RooRealVar n1( "n_{1}", "n_{1}", 3 ,0,40);
  RooRealVar n1( "n_{1}", "n_{1}", 2 ,1.01,5000.);
  RooRealVar alpha2( "alpha_{2}", "alpha_{2}", 0.1,0.,20.);
  //  RooRealVar n2( "n_{2}", "n_{2}", 0.81 ,0,40);
  RooRealVar n2( "n_{2}", "n_{2}", 2 ,1.01,5000.);

  DoubleCBPdf doubleCB("doubleCB","doubleCB",res,mean,cbSigma,alpha1,n1,alpha2,n2);
  RooAddPdf model("model", "model", RooArgList(doubleCB), RooArgList(nsig));

  RooDataHist data("res","E^{reco}/E^{gen}",res,hist);
 
  model.fitTo(data,RooFit::FitOptions("mh"),RooFit::Optimize(0),RooFit::Timer(1));
  model.fitTo(data,RooFit::Optimize(1),RooFit::Timer(0),RooFit::PrintEvalErrors(-1),RooFit::Save(1)); 
  
  double errors[7] ={(mean.getAsymErrorHi()+mean.getAsymErrorLo())/2., (cbSigma.getAsymErrorHi()+cbSigma.getAsymErrorLo())/2., (alpha1.getAsymErrorHi()+alpha1.getAsymErrorLo())/2., (n1.getAsymErrorHi()+n1.getAsymErrorLo())/2., (alpha2.getAsymErrorHi()+alpha2.getAsymErrorLo())/2., (n2.getAsymErrorHi()+n2.getAsymErrorLo())/2., 0.};
   
  hist->GetXaxis()->SetRangeUser(0.2,1.8); 
  double maximum = hist->GetMaximum();
  hist->GetXaxis()->SetRangeUser(0.,2.); 
  TF1* DoubleCB = new TF1("DoubleCB",&mydoublecb,xmin,xmax,7);
  DoubleCB->SetParameters(mean.getVal(), cbSigma.getVal(), alpha1.getVal(), n1.getVal(), alpha2.getVal(), n2.getVal(), maximum);
  DoubleCB->SetParErrors(errors);

  return DoubleCB;
}

TF1* fitHisto(TH1* hist)
{
   hist->GetXaxis()->SetRangeUser(0.2,1.8); 
   int binmax = hist->GetMaximumBin(); 
   hist->GetXaxis()->SetRangeUser(0.,2.); 
   double xMAX = hist->GetXaxis()->GetBinCenter(binmax);  

   std::cout << "fitHisto: " << hist->GetName() << std::endl;
   if(fitFunction_=="doubleCB") return makeDoubleCBFit(hist,xMAX-1.5*hist->GetRMS(),xMAX+1.5*hist->GetRMS());
   return makeCruijffFit(hist,xMAX-1.5*hist->GetRMS(),xMAX+1.5*hist->GetRMS());
}

void drawHistFunc(TH1F* hist, TF1* func, std::string x_label, std::string Name)
{
   gStyle->SetOptStat(0000); 
   hist->SetMaximum(hist->GetMaximum()*1.05);
   hist->SetLineColor(kBlack);
   hist->SetMarkerColor(kBlack);
   hist->SetLineWidth(2);
   hist->GetXaxis()->SetTitle(x_label.c_str());
   hist->GetXaxis()->SetRangeUser(hist->GetMean()-1.,hist->GetMean()+1.);   

   func->SetLineColor(kRed+1);

   double mean = hist->GetMean();
   double mean_error = hist->GetMeanError();
   double rms = hist->GetRMS();
   double rms_error = hist->GetRMSError();
   double mu = func->GetParameter(0);
   double mu_error = func->GetParError(0);
   double sigma = (func->GetParameter(1)+func->GetParameter(2))/2.;
   double sigma_error = 0.5*sqrt(func->GetParError(2)*func->GetParError(2) +func->GetParError(1)*func->GetParError(1));
   double score = sigma*fabs(mu-1)*10000;
   double score_error = score*sqrt((mu_error/mu)*(mu_error/mu) + (sigma_error/sigma)*(sigma_error/sigma));
   double res = sigma/mu;
   double res_error = res*sqrt((mu_error/mu)*(mu_error/mu) + (sigma_error/sigma)*(sigma_error/sigma));
   
   TLegend* legend = new TLegend(0.57, 0.60, 0.72, 0.89);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(hist,std::string("mean = "+to_string(mean)+" +/- "+to_string(mean_error)).c_str(),"");
   legend -> AddEntry(hist,std::string("rms = "+to_string(rms)+" +/- "+to_string(rms_error)).c_str(),"");
   legend -> AddEntry(hist,std::string("#sigma = "+to_string(sigma)+" +/- "+to_string(sigma_error)).c_str(),"");
   legend -> AddEntry(hist,std::string("#mu = "+to_string(mu)+" +/- "+to_string(mu_error)).c_str(),"");
   legend -> AddEntry(hist,std::string("#sigma = "+to_string(sigma)+" +/- "+to_string(sigma_error)).c_str(),"");
   legend -> AddEntry(hist,std::string("res = "+to_string(res)+" +/- "+to_string(res_error)).c_str(),"");
   legend -> AddEntry(hist,std::string("score = "+to_string(score)+" +/- "+to_string(score_error)).c_str(),"");
     
   TCanvas* c = new TCanvas();
   hist->Draw("E");    
   func->Draw("L,same");
   legend -> Draw("same");  
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf"); 
   gStyle->SetOptStat(1110); 
   
}

