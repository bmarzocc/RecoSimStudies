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
TH1F* EoEtrue_Mustache_seedEt_0_10_EB;
TH1F* EoEtrue_Mustache_seedEt_10_20_EB;
TH1F* EoEtrue_Mustache_seedEt_20_30_EB;
TH1F* EoEtrue_Mustache_seedEt_30_40_EB;
TH1F* EoEtrue_Mustache_seedEt_40_50_EB;
TH1F* EoEtrue_Mustache_seedEt_50_60_EB;
TH1F* EoEtrue_Mustache_seedEt_60_70_EB;
TH1F* EoEtrue_Mustache_seedEt_70_80_EB;
TH1F* EoEtrue_Mustache_seedEt_80_90_EB;
TH1F* EoEtrue_Mustache_seedEt_90_100_EB;
TH1F* EoEtrue_Mustache_seedEt_0_10_EE;
TH1F* EoEtrue_Mustache_seedEt_10_20_EE;
TH1F* EoEtrue_Mustache_seedEt_20_30_EE;
TH1F* EoEtrue_Mustache_seedEt_30_40_EE;
TH1F* EoEtrue_Mustache_seedEt_40_50_EE;
TH1F* EoEtrue_Mustache_seedEt_50_60_EE;
TH1F* EoEtrue_Mustache_seedEt_60_70_EE;
TH1F* EoEtrue_Mustache_seedEt_70_80_EE;
TH1F* EoEtrue_Mustache_seedEt_80_90_EE;
TH1F* EoEtrue_Mustache_seedEt_90_100_EE;

std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_0_10_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_10_20_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_20_30_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_30_40_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_40_50_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_50_60_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_60_70_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_70_80_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_80_90_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_90_100_EB;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_0_10_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_10_20_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_20_30_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_30_40_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_40_50_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_50_60_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_60_70_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_70_80_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_80_90_EE;
std::vector<TH1F*> EoEtrue_vs_DnnThreshold_seedEt_90_100_EE;


//DEFINE BRANCHES
double en_cluster;
double en_true;
double et_seed;
double seed_eta;
double et_true;
double EoEtrue;
double dnn_thre;

TBranch *b_en_cluster;   //!
TBranch *b_en_true;   //!
TBranch *b_et_seed;   //!
TBranch *b_seed_eta;   //!
TBranch *b_et_true;   //!
TBranch *b_EoEtrue;   //!
TBranch *b_dnn_thre;   //!

//setTreeBranches
void setTreeBranches(TTree* tree)
{
   tree->SetBranchAddress("en_cluster", &en_cluster, &b_en_cluster);
   tree->SetBranchAddress("en_true", &en_true, &b_en_true);
   tree->SetBranchAddress("et_seed", &et_seed, &b_et_seed);
   tree->SetBranchAddress("seed_eta", &seed_eta, &b_seed_eta);
   tree->SetBranchAddress("et_true", &et_true, &b_et_true);
   tree->SetBranchAddress("EoEtrue", &EoEtrue, &b_EoEtrue);
   tree->SetBranchAddress("dnn_thre", &dnn_thre, &b_dnn_thre); 
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
void setHistograms(std::vector<std::string> dnnCuts)
{
   EoEtrue_Mustache_seedEt_0_10_EB = new TH1F("EoEtrue_Mustache_seedEt_0_10_EB","EoEtrue_Mustache_seedEt_0_10_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_10_20_EB = new TH1F("EoEtrue_Mustache_seedEt_10_20_EB","EoEtrue_Mustache_seedEt_10_20_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_20_30_EB = new TH1F("EoEtrue_Mustache_seedEt_20_30_EB","EoEtrue_Mustache_seedEt_20_30_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_30_40_EB = new TH1F("EoEtrue_Mustache_seedEt_30_40_EB","EoEtrue_Mustache_seedEt_30_40_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_40_50_EB = new TH1F("EoEtrue_Mustache_seedEt_40_50_EB","EoEtrue_Mustache_seedEt_40_50_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_50_60_EB = new TH1F("EoEtrue_Mustache_seedEt_50_60_EB","EoEtrue_Mustache_seedEt_50_60_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_60_70_EB = new TH1F("EoEtrue_Mustache_seedEt_60_70_EB","EoEtrue_Mustache_seedEt_60_70_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_70_80_EB = new TH1F("EoEtrue_Mustache_seedEt_70_80_EB","EoEtrue_Mustache_seedEt_70_80_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_80_90_EB = new TH1F("EoEtrue_Mustache_seedEt_80_90_EB","EoEtrue_Mustache_seedEt_80_90_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_90_100_EB = new TH1F("EoEtrue_Mustache_seedEt_90_100_EB","EoEtrue_Mustache_seedEt_90_100_EB",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_0_10_EE = new TH1F("EoEtrue_Mustache_seedEt_0_10_EE","EoEtrue_Mustache_seedEt_0_10_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_10_20_EE = new TH1F("EoEtrue_Mustache_seedEt_10_20_EE","EoEtrue_Mustache_seedEt_10_20_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_20_30_EE = new TH1F("EoEtrue_Mustache_seedEt_20_30_EE","EoEtrue_Mustache_seedEt_20_30_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_30_40_EE = new TH1F("EoEtrue_Mustache_seedEt_30_40_EE","EoEtrue_Mustache_seedEt_30_40_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_40_50_EE = new TH1F("EoEtrue_Mustache_seedEt_40_50_EE","EoEtrue_Mustache_seedEt_40_50_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_50_60_EE = new TH1F("EoEtrue_Mustache_seedEt_50_60_EE","EoEtrue_Mustache_seedEt_50_60_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_60_70_EE = new TH1F("EoEtrue_Mustache_seedEt_60_70_EE","EoEtrue_Mustache_seedEt_60_70_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_70_80_EE = new TH1F("EoEtrue_Mustache_seedEt_70_80_EE","EoEtrue_Mustache_seedEt_70_80_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_80_90_EE = new TH1F("EoEtrue_Mustache_seedEt_80_90_EE","EoEtrue_Mustache_seedEt_80_90_EE",10000, 0., 10.);
   EoEtrue_Mustache_seedEt_90_100_EE = new TH1F("EoEtrue_Mustache_seedEt_90_100_EE","EoEtrue_Mustache_seedEt_90_100_EE",10000, 0., 10.);
   
   EoEtrue_vs_DnnThreshold_seedEt_0_10_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_10_20_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_20_30_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_30_40_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_40_50_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_50_60_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_60_70_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_70_80_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_80_90_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_90_100_EB.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_0_10_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_10_20_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_20_30_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_30_40_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_40_50_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_50_60_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_60_70_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_70_80_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_80_90_EE.resize(dnnCuts.size());
   EoEtrue_vs_DnnThreshold_seedEt_90_100_EE.resize(dnnCuts.size());

   for(unsigned int iBin=0; iBin<dnnCuts.size(); iBin++)
   {
      dnnCuts.at(iBin).replace(1,1,string("_"));
      EoEtrue_vs_DnnThreshold_seedEt_0_10_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_0_10_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_0_10_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.);    
      EoEtrue_vs_DnnThreshold_seedEt_10_20_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_10_20_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_10_20_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_20_30_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_20_30_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_20_30_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_30_40_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_30_40_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_30_40_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_40_50_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_40_50_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_40_50_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_50_60_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_50_60_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_50_60_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_60_70_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_60_70_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_60_70_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_70_80_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_70_80_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_70_80_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_80_90_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_80_90_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_80_90_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_90_100_EB[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_90_100_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_90_100_DNN_"+dnnCuts.at(iBin)+"_EB").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_0_10_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_0_10_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_0_10_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.);    
      EoEtrue_vs_DnnThreshold_seedEt_10_20_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_10_20_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_10_20_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_20_30_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_20_30_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_20_30_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_30_40_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_30_40_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_30_40_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_40_50_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_40_50_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_40_50_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_50_60_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_50_60_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_50_60_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_60_70_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_60_70_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_60_70_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_70_80_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_70_80_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_70_80_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_80_90_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_80_90_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_80_90_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.); 
      EoEtrue_vs_DnnThreshold_seedEt_90_100_EE[iBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_90_100_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_90_100_DNN_"+dnnCuts.at(iBin)+"_EE").c_str(), 10000, 0., 10.);    
   }
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

std::pair<TGraphErrors*,TGraphErrors*> makeFitProfile(std::vector<TH1F*>* vecHist, std::vector<std::string> dnnCuts, std::string xTitle, std::string fitFunction_="cruijff", bool doEffective = false)
{
   TF1* doubleCB;
   int nBins = vecHist->size();
   double x[nBins], x_error[nBins], y_mean[nBins], y_meanError[nBins], y_sigma[nBins], y_sigmaError[nBins];
   
   for(unsigned int iBin=0; iBin<vecHist->size(); iBin++)
   {
       x[iBin] = std::stod(dnnCuts.at(iBin)); 
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

TGraphErrors* scaledGraph(TGraphErrors* gr_SuperCluster, float Mustache_val)
{
   TGraphErrors* gr = new TGraphErrors();
   double x,y;
   for(int i=0; i<gr_SuperCluster->GetN(); i++)
   {
       double y_err = gr_SuperCluster->GetErrorY(i);
       gr_SuperCluster->GetPoint(i,x,y);
       gr->SetPoint(i,x,y/Mustache_val);
       gr->SetPointError(i,0,y_err/Mustache_val);
   }  
   return gr;
}

void drawGraph(TGraphErrors* gr_SuperCluster_tmp, float Mustache_val, std::string xtitle, std::string ytitle, std::string Name, std::string refLegend="DeepSC", std::string valLegend="Mustache",float y_min=-1., float y_max=-1.)
{ 
   gStyle->SetOptStat(0000);  

   TGraphErrors* gr_SuperCluster; 
   if(Mustache_val!=-1) gr_SuperCluster = scaledGraph(gr_SuperCluster_tmp,Mustache_val);
   else gr_SuperCluster = gr_SuperCluster_tmp; 
  
   gr_SuperCluster->SetTitle(Name.c_str());
   gr_SuperCluster->SetLineColor(kBlue+1);
   gr_SuperCluster->SetMarkerStyle(20);
   gr_SuperCluster->SetMarkerSize(0.5);
   gr_SuperCluster->SetMarkerColor(kBlue+1);
   gr_SuperCluster->SetLineWidth(2);
   gr_SuperCluster->GetXaxis()->SetTitle(xtitle.c_str()); 
   gr_SuperCluster->GetYaxis()->SetTitle(ytitle.c_str()); 

   TLine* Mustache_line = new TLine(0.5,1., 1.,1.);
   Mustache_line->SetLineWidth(2);
   Mustache_line->SetLineColor(kGreen+1); 
   
   float min = y_min;
   float max = y_max;
   if(y_min<0. || y_max<0.){
      min = 0.98*computeRange(gr_SuperCluster).first; 
      max = 1.02*computeRange(gr_SuperCluster).second;
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
   legend -> AddEntry(Mustache_line,valLegend.c_str(),"L");

   TCanvas* c = new TCanvas();
   gr_SuperCluster->Draw("AL");
   if(1.>=min && 1.<=max) Mustache_line->Draw("L, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   gStyle->SetOptStat(1110);
}

//draw plots
void drawPlots(std::string fitFunction_, std::vector<std::string> dnnCuts)
{
   
   //68% Quantile
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_0_10_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_0_10_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_10_20_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_10_20_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_20_30_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_20_30_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_30_40_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_30_40_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_40_50_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_40_50_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_50_60_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_50_60_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_60_70_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_60_70_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_70_80_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_70_80_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_80_90_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_80_90_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_90_100_eff_EB = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_90_100_EB,dnnCuts,std::string("DNN Score"),fitFunction_,true); 

   drawGraph(gr_EoEtrue_vs_DNN_seedEt_0_10_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_0_10_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_0_10_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_0_10_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_0_10_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_0_10_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_10_20_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_10_20_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_10_20_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_10_20_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_10_20_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_10_20_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_20_30_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_20_30_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_20_30_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_20_30_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_20_30_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_20_30_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_30_40_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_30_40_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_30_40_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_30_40_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_30_40_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_30_40_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_40_50_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_40_50_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_40_50_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_40_50_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_40_50_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_40_50_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_50_60_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_50_60_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_50_60_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_50_60_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_50_60_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_50_60_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_60_70_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_60_70_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_60_70_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_60_70_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_60_70_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_60_70_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_70_80_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_70_80_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_70_80_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_70_80_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_70_80_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_70_80_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_80_90_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_80_90_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_80_90_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_80_90_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_80_90_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_80_90_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_90_100_eff_EB.first,evaluateHisto(EoEtrue_Mustache_seedEt_90_100_EB,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_90_100_Effective_EB"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_90_100_eff_EB.second,evaluateHisto(EoEtrue_Mustache_seedEt_90_100_EB,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_90_100_Effective_EB"));

   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_0_10_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_0_10_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_10_20_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_10_20_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_20_30_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_20_30_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_30_40_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_30_40_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_40_50_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_40_50_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_50_60_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_50_60_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_60_70_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_60_70_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_70_80_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_70_80_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_80_90_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_80_90_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true);
   std::pair<TGraphErrors*,TGraphErrors*> gr_EoEtrue_vs_DNN_seedEt_90_100_eff_EE = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_90_100_EE,dnnCuts,std::string("DNN Score"),fitFunction_,true); 

   drawGraph(gr_EoEtrue_vs_DNN_seedEt_0_10_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_0_10_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_0_10_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_0_10_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_0_10_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_0_10_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_10_20_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_10_20_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_10_20_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_10_20_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_10_20_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_10_20_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_20_30_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_20_30_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_20_30_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_20_30_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_20_30_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_20_30_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_30_40_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_30_40_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_30_40_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_30_40_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_30_40_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_30_40_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_40_50_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_40_50_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_40_50_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_40_50_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_40_50_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_40_50_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_50_60_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_50_60_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_50_60_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_50_60_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_50_60_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_50_60_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_60_70_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_60_70_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_60_70_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_60_70_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_60_70_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_60_70_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_70_80_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_70_80_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_70_80_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_70_80_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_70_80_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_70_80_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_80_90_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_80_90_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_80_90_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_80_90_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_80_90_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_80_90_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_90_100_eff_EE.first,evaluateHisto(EoEtrue_Mustache_seedEt_90_100_EE,fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_90_100_Effective_EE"));
   drawGraph(gr_EoEtrue_vs_DNN_seedEt_90_100_eff_EE.second,evaluateHisto(EoEtrue_Mustache_seedEt_90_100_EE,fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_90_100_Effective_EE"));


}

