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

#include <boost/algorithm/string/replace.hpp>

using namespace std;
using namespace edm;

//DEFINE HISTOGRAMS 
std::vector<TH1F*> EoEtrue_Mustache_seedEt_EB;
std::vector<TH1F*> EoEtrue_Mustache_seedEt_EE;
std::vector<TH1F*> EoEtrue_Mustache_seedEta;
std::vector<std::vector<TH1F*>> EoEtrue_Mustache_seedEta_seedEt;

std::vector<std::vector<TH1F*>> EoEtrue_vs_DnnThreshold_seedEt_EB;
std::vector<std::vector<TH1F*>> EoEtrue_vs_DnnThreshold_seedEt_EE;
std::vector<std::vector<TH1F*>> EoEtrue_vs_DnnThreshold_seedEta;
std::vector<std::vector<std::vector<TH1F*>>> EoEtrue_vs_DnnThreshold_seedEta_seedEt;

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
void setHistograms(std::vector<std::string> dnnCuts, std::vector<std::string> etCuts, std::vector<std::string> etaCuts)
{
   EoEtrue_Mustache_seedEt_EB.resize(etCuts.size()-1);
   EoEtrue_Mustache_seedEt_EE.resize(etCuts.size()-1);
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
       boost::replace_all(etCuts.at(iBin), ".", "_");
       boost::replace_all(etCuts.at(iBin+1), ".", "_");   
       EoEtrue_Mustache_seedEt_EB[iBin] = new TH1F(string("EoEtrue_Mustache_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_EB").c_str(), string("EoEtrue_Mustache_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_EB").c_str(),10000, 0., 10.);
       EoEtrue_Mustache_seedEt_EE[iBin] = new TH1F(string("EoEtrue_Mustache_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_EE").c_str(), string("EoEtrue_Mustache_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_EE").c_str(),10000, 0., 10.); 
   }   

   EoEtrue_Mustache_seedEta.resize(etaCuts.size()-1);
   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
       boost::replace_all(etaCuts.at(iBin), ".", "_");
       boost::replace_all(etaCuts.at(iBin+1), ".", "_");   
       EoEtrue_Mustache_seedEta[iBin] = new TH1F(string("EoEtrue_Mustache_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)).c_str(), string("EoEtrue_Mustache_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)).c_str(),10000, 0., 10.);
   } 

   EoEtrue_Mustache_seedEta_seedEt.resize(etCuts.size()-1);
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       EoEtrue_Mustache_seedEta_seedEt[iBin].resize(etaCuts.size()-1);
       for(unsigned int jBin=0; jBin<etaCuts.size()-1; jBin++)
       {
           EoEtrue_Mustache_seedEta_seedEt[iBin][jBin] = new TH1F(std::string("EoEtrue_Mustache_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(iBin+1)).c_str(), std::string("EoEtrue_Mustache_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_seedEta_"+etaCuts.at(jBin)+"_"+etaCuts.at(iBin+1)).c_str(), 10000, 0., 10.);  
       }    
   }
 
   EoEtrue_vs_DnnThreshold_seedEt_EB.resize(etCuts.size()-1); 
   EoEtrue_vs_DnnThreshold_seedEt_EE.resize(etCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++)
   {
       EoEtrue_vs_DnnThreshold_seedEt_EB[iBin].resize(dnnCuts.size());
       EoEtrue_vs_DnnThreshold_seedEt_EE[iBin].resize(dnnCuts.size());
       for(unsigned int jBin=0; jBin<dnnCuts.size(); jBin++)
       {
           boost::replace_all(dnnCuts.at(jBin), ".", "_");
           EoEtrue_vs_DnnThreshold_seedEt_EB[iBin][jBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_DNN_"+dnnCuts.at(jBin)+"_EB").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_DNN_"+dnnCuts.at(jBin)+"_EB").c_str(), 10000, 0., 10.);
           EoEtrue_vs_DnnThreshold_seedEt_EE[iBin][jBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_DNN_"+dnnCuts.at(jBin)+"_EE").c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_DNN_"+dnnCuts.at(jBin)+"_EE").c_str(), 10000, 0., 10.);    
       }     
   }

   EoEtrue_vs_DnnThreshold_seedEta.resize(etaCuts.size()-1);      
   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++)
   {    
       EoEtrue_vs_DnnThreshold_seedEta[iBin].resize(dnnCuts.size());
       for(unsigned int jBin=0; jBin<dnnCuts.size(); jBin++)
       { 
           EoEtrue_vs_DnnThreshold_seedEta[iBin][jBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_DNN_"+dnnCuts.at(jBin)).c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_DNN_"+dnnCuts.at(jBin)).c_str(), 10000, 0., 10.);
       }     
   }

   EoEtrue_vs_DnnThreshold_seedEta_seedEt.resize(dnnCuts.size());
   for(unsigned int iBin=0; iBin<dnnCuts.size(); iBin++)
   {
       EoEtrue_vs_DnnThreshold_seedEta_seedEt[iBin].resize(etCuts.size()-1); 
       for(unsigned int jBin=0; jBin<etCuts.size()-1; jBin++)
       {
           EoEtrue_vs_DnnThreshold_seedEta_seedEt[iBin][jBin].resize(etaCuts.size()-1); 
           for(unsigned int kBin=0; kBin<etaCuts.size()-1; kBin++)
           {
               EoEtrue_vs_DnnThreshold_seedEta_seedEt[iBin][jBin][kBin] = new TH1F(std::string("EoEtrue_vs_DnnThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_DNN_"+dnnCuts.at(iBin)).c_str(), std::string("EoEtrue_vs_DnnThreshold_seedEt_"+etCuts.at(jBin)+"_"+etCuts.at(jBin+1)+"_seedEta_"+etaCuts.at(kBin)+"_"+etaCuts.at(kBin+1)+"_DNN_"+dnnCuts.at(iBin)).c_str(), 10000, 0., 10.);
           }
       }
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

   TLine* Mustache_line = new TLine(0.3,1., 1.,1.);
   Mustache_line->SetLineWidth(2);
   Mustache_line->SetLineColor(kGreen+1); 
   
   float min = y_min;
   float max = y_max;
   if(y_min<0. || y_max<0.){
      min = 0.98*computeRange(gr_SuperCluster).first; 
      max = 1.02*computeRange(gr_SuperCluster).second;
   }
   if(min>1.) min = 0.99;
   if(max<1.) max = 1.01;
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
   Mustache_line->Draw("L, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

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

   TH2F* h2_clone = (TH2F*)h2->Clone();
   h2_clone->Reset(); 
   for(int binx=1; binx<=h2->GetNbinsX(); binx++){ 
       for(int biny=1; biny<=h2->GetNbinsY(); biny++){
           if(h2->GetBinContent(binx,biny)<=1.){ 
              h2_clone->SetBinContent(binx,biny,fabs(1.-h2->GetBinContent(binx,biny))*100.);
           } 
       }
   }

   //std::sort(vals.begin(),vals.begin()); 
   //h2_clone->GetZaxis()->SetRangeUser(0.99*vals.at(0),1.);

   c->cd();
   h2_clone->Draw("COLZ");
   c->SaveAs(std::string(Name+"_Gain.png").c_str(),"png");
   c->SaveAs(std::string(Name+"_Gain.pdf").c_str(),"pdf");

   gStyle->SetOptStat(1110); 
}

//draw plots
void drawPlots(std::string fitFunction_, std::vector<std::string> etCuts, std::vector<std::string> etaCuts, std::vector<std::string> dnnCuts)
{
   
   //68% Quantile
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_EoEtrue_vs_DNN_seedEt_eff_EB;
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_EoEtrue_vs_DNN_seedEt_eff_EE; 
   std::vector<std::pair<TGraphErrors*,TGraphErrors*>> gr_EoEtrue_vs_DNN_seedEta_eff; 
   gr_EoEtrue_vs_DNN_seedEt_eff_EB.resize(etCuts.size()-1);
   gr_EoEtrue_vs_DNN_seedEt_eff_EE.resize(etCuts.size()-1); 
   gr_EoEtrue_vs_DNN_seedEta_eff.resize(etaCuts.size()-1);  

   for(unsigned int iBin=0; iBin<etCuts.size()-1; iBin++){
       gr_EoEtrue_vs_DNN_seedEt_eff_EB[iBin] = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_EB[iBin],dnnCuts,std::string("DNN Score"),fitFunction_,true); 
       drawGraph(gr_EoEtrue_vs_DNN_seedEt_eff_EB[iBin].first,evaluateHisto(EoEtrue_Mustache_seedEt_EB[iBin],fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EB")); 
       drawGraph(gr_EoEtrue_vs_DNN_seedEt_eff_EB[iBin].second,evaluateHisto(EoEtrue_Mustache_seedEt_EB[iBin],fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EB"));

       gr_EoEtrue_vs_DNN_seedEt_eff_EE[iBin] = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEt_EE[iBin],dnnCuts,std::string("DNN Score"),fitFunction_,true); 
       drawGraph(gr_EoEtrue_vs_DNN_seedEt_eff_EE[iBin].first,evaluateHisto(EoEtrue_Mustache_seedEt_EE[iBin],fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EE")); 
       drawGraph(gr_EoEtrue_vs_DNN_seedEt_eff_EE[iBin].second,evaluateHisto(EoEtrue_Mustache_seedEt_EE[iBin],fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEt_"+etCuts.at(iBin)+"_"+etCuts.at(iBin+1)+"_Effective_EE")); 
   }

   for(unsigned int iBin=0; iBin<etaCuts.size()-1; iBin++){
       gr_EoEtrue_vs_DNN_seedEta_eff[iBin] = makeFitProfile(&EoEtrue_vs_DnnThreshold_seedEta[iBin],dnnCuts,std::string("DNN Score"),fitFunction_,true); 
       drawGraph(gr_EoEtrue_vs_DNN_seedEta_eff[iBin].first,evaluateHisto(EoEtrue_Mustache_seedEta[iBin],fitFunction_,true)[0], std::string("DNN Score"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("EoEtrue_vs_DNN_Mean_seedEt_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_Effective")); 
       drawGraph(gr_EoEtrue_vs_DNN_seedEta_eff[iBin].second,evaluateHisto(EoEtrue_Mustache_seedEta[iBin],fitFunction_,true)[2], std::string("DNN Score"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("EoEtrue_vs_DNN_Resolution_seedEta_"+etaCuts.at(iBin)+"_"+etaCuts.at(iBin+1)+"_Effective"));
   }
 
   TH2F* Ratio_deepMust;
   std::vector<float> vecEtCuts;
   std::vector<float> vecEtaCuts;

   for(unsigned int iBin=0; iBin<etCuts.size(); iBin++){
       boost::replace_all(etCuts.at(iBin), "_", ".");
       vecEtCuts.push_back(std::stof(etCuts.at(iBin)));
   } 
   for(unsigned int iBin=0; iBin<etaCuts.size(); iBin++){
       boost::replace_all(etaCuts.at(iBin), "_", ".");
       vecEtaCuts.push_back(std::stof(etaCuts.at(iBin)));
   } 
  
   std::pair<TH2F*,TH2F*> h2_must = makeFitHist2D(-1,EoEtrue_Mustache_seedEta_seedEt,vecEtCuts,vecEtaCuts,fitFunction_,true);
   for(unsigned int iBin=0; iBin<dnnCuts.size(); iBin++){
       std::pair<TH2F*,TH2F*> h2_deepSC = makeFitHist2D(iBin,EoEtrue_vs_DnnThreshold_seedEta_seedEt[iBin],vecEtCuts,vecEtaCuts,fitFunction_,true);
       Ratio_deepMust = (TH2F*)(h2_deepSC.first)->Clone();
       Ratio_deepMust->Divide(h2_must.first);
       drawH2(Ratio_deepMust, std::string("Et (GeV)"), std::string("#eta"), std::string("#mu_DeepSC/#mu_Mustache"), std::string("h2_EoEtrue_Mean_DNN_"+dnnCuts.at(iBin)));
       Ratio_deepMust = (TH2F*)(h2_deepSC.second)->Clone();
       Ratio_deepMust->Divide(h2_must.second);
       drawH2(Ratio_deepMust, std::string("Et (GeV)"), std::string("#eta"), std::string("#sigma_DeepSC/#sigma_Mustache"), std::string("h2_EoEtrue_Resolution_DNN_"+dnnCuts.at(iBin))); 
   }
   
}

