#!/usr/bin/python
import numpy as n
from ROOT import *
import sys, getopt
from array import array
from optparse import OptionParser
import operator
import math  

def drawHisto(h1,name,var):

   gStyle.SetOptStat(0000)

   h1.SetMarkerStyle(20)
   h1.SetMarkerColor(kBlack) 
   h1.SetLineColor(kBlack) 
   h1.SetLineWidth(2) 
   h1.Scale(1./h1.Integral())
   h1.GetXaxis().SetTitle(var)
   
   maximum = 1.01*h1.GetMaximum()
   
   h1.GetYaxis().SetRangeUser(0.,maximum)
   c = TCanvas()
   h1.Draw("HIST")
   c.SaveAs(name+".png","png") 
   c.SaveAs(name+".pdf","pdf") 

   h1.GetYaxis().SetRangeUser(0.001,maximum*20.)
   c.SetLogy() 
   h1.Draw("HIST")
   c.SaveAs(name+"_log.png","png") 
   c.SaveAs(name+"_log.pdf","pdf") 

def drawHistos(h1,h2,name,var):

   gStyle.SetOptStat(0000)

   h1.SetMarkerStyle(20)
   h1.SetMarkerColor(kRed+1) 
   h1.SetLineColor(kRed+1) 
   h1.SetLineWidth(2) 
   h1.Scale(1./h1.Integral())
   h1.GetXaxis().SetTitle(var)
   
   h2.SetMarkerStyle(20)
   h2.SetMarkerColor(kBlack) 
   h2.SetLineColor(kBlack) 
   h2.SetLineWidth(2) 
   h2.Scale(1./h2.Integral()) 

   legend = TLegend(0.70, 0.82, 0.95, 0.94)
   legend.SetFillColor(kWhite)
   legend.SetFillStyle(1000)
   legend.SetLineWidth(0)
   legend.SetLineColor(kWhite)
   legend.SetTextFont(42) 
   legend.SetTextSize(0.04)
   legend.AddEntry(h1,"caloMatched","L")
   legend.AddEntry(h2,"caloUnmatched","L")

   maximum = 1.01*h1.GetMaximum()
   if h2.GetMaximum()>h1.GetMaximum(): maximum = 1.01*h2.GetMaximum()  

   h1.GetYaxis().SetRangeUser(0.,maximum)
   c = TCanvas()
   h1.Draw("HIST")
   h2.Draw("HIST,same")
   legend.Draw("same")
   c.SaveAs(name+".png","png") 
   c.SaveAs(name+".pdf","pdf") 

   h1.GetYaxis().SetRangeUser(0.001,maximum*20.)
   c.SetLogy() 
   h1.Draw("HIST")
   h2.Draw("HIST,same")
   legend.Draw("same")
   c.SaveAs(name+"_log.png","png") 
   c.SaveAs(name+"_log.pdf","pdf") 
 

if __name__ == '__main__':

  inFile = '/eos/cms/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Clustering/FourElectronsGammasGunPt1-100_pythia8_withPU_withTracker_106X_mcRun3_2021_realistic_v3_RAW_StdSeedingGathering_Mustache_Dumper_v2_Total.root'
 
  print "inputFile = ",inFile
  
  f_in = TFile(inFile)
  tree = f_in.Get('recosimdumper/caloTree')

  h_seed_nXtals_EB = TH1F("h_seed_nXtals_EB","h_seed_nXtals",25,0.5,25.5)
  tree.Draw("pfCluster_nXtals[superCluster_seedIndex]>>h_seed_nXtals_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_nXtals_EB,"h_seed_nXtals_EB","nXtals")

  h_seed_nXtals_EE = TH1F("h_seed_nXtals_EE","h_seed_nXtals",25,0.5,25.5)
  tree.Draw("pfCluster_nXtals[superCluster_seedIndex]>>h_seed_nXtals_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_nXtals_EE,"h_seed_nXtals_EE","nXtals")

  h_seed_R9_EB = TH1F("h_seed_R9_EB","h_seed_R9",80,0.6,1.4)
  tree.Draw("pfCluster_full5x5_r9[superCluster_seedIndex]>>h_seed_R9_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_R9_EB,"h_seed_R9_EB","R9")

  h_seed_R9_EE = TH1F("h_seed_R9_EE","h_seed_R9",80,0.6,1.4)
  tree.Draw("pfCluster_full5x5_r9[superCluster_seedIndex]>>h_seed_R9_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_R9_EE,"h_seed_R9_EE","R9")

  h_seed_etaWidth_EB = TH1F("h_seed_etaWidth_EB","h_seed_etaWidth",100,0.,0.016)
  tree.Draw("pfCluster_etaWidth[superCluster_seedIndex]>>h_seed_etaWidth_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_etaWidth_EB,"h_seed_etaWidth_EB","#eta-Width")

  h_seed_etaWidth_EE = TH1F("h_seed_etaWidth_EE","h_seed_etaWidth",100,0.,0.06)
  tree.Draw("pfCluster_etaWidth[superCluster_seedIndex]>>h_seed_etaWidth_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_etaWidth_EE,"h_seed_etaWidth_EE","#eta-Width")

  h_seed_phiWidth_EB = TH1F("h_seed_phiWidth_EB","h_seed_phiWidth",100,0.,0.018)
  tree.Draw("pfCluster_phiWidth[superCluster_seedIndex]>>h_seed_phiWidth_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_phiWidth_EB,"h_seed_phiWidth_EB","#phi-Width")

  h_seed_phiWidth_EE = TH1F("h_seed_phiWidth_EE","h_seed_phiWidth",100,0.,0.06)
  tree.Draw("pfCluster_phiWidth[superCluster_seedIndex]>>h_seed_phiWidth_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_phiWidth_EE,"h_seed_phiWidth_EE","#phi-Width")
   
  h_seed_sigmaIetaIeta_EB = TH1F("h_seed_sigmaIetaIeta_EB","h_seed_sigmaIetaIeta",100,0.,0.025)
  tree.Draw("pfCluster_full5x5_sigmaIetaIeta[superCluster_seedIndex]>>h_seed_sigmaIetaIeta_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_sigmaIetaIeta_EB,"h_seed_sigmaIetaIeta_EB","#sigma_{i#etai#eta}")

  h_seed_sigmaIetaIeta_EE = TH1F("h_seed_sigmaIetaIeta_EE","h_seed_sigmaIetaIeta",300,0.,0.08)
  tree.Draw("pfCluster_full5x5_sigmaIetaIeta[superCluster_seedIndex]>>h_seed_sigmaIetaIeta_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_sigmaIetaIeta_EE,"h_seed_sigmaIetaIeta_EE","#sigma_{i#etai#eta}")

  h_seed_sigmaIetaIphi_EB = TH1F("h_seed_sigmaIetaIphi_EB","h_seed_sigmaIetaIphi",100,-0.001,0.001)
  tree.Draw("pfCluster_full5x5_sigmaIetaIphi[superCluster_seedIndex]>>h_seed_sigmaIetaIphi_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_sigmaIetaIphi_EB,"h_seed_sigmaIetaIphi_EB","#sigma_{i#etai#phi}")

  h_seed_sigmaIetaIphi_EE = TH1F("h_seed_sigmaIetaIphi_EE","h_seed_sigmaIetaIphi",100,-0.005,0.005)
  tree.Draw("pfCluster_full5x5_sigmaIetaIphi[superCluster_seedIndex]>>h_seed_sigmaIetaIphi_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_sigmaIetaIphi_EE,"h_seed_sigmaIetaIphi_EE","#sigma_{i#etai#phi}")

  h_seed_sigmaIphiIphi_EB = TH1F("h_seed_sigmaIphiIphi_EB","h_seed_sigmaIphiIphi",200,0.,0.05)
  tree.Draw("pfCluster_full5x5_sigmaIphiIphi[superCluster_seedIndex]>>h_seed_sigmaIphiIphi_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_sigmaIphiIphi_EB,"h_seed_sigmaIphiIphi_EB","#sigma_{i#phii#phi}")

  h_seed_sigmaIphiIphi_EE = TH1F("h_seed_sigmaIphiIphi_EE","h_seed_sigmaIphiIphi",400,0.,0.08)
  tree.Draw("pfCluster_full5x5_sigmaIphiIphi[superCluster_seedIndex]>>h_seed_sigmaIphiIphi_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_sigmaIphiIphi_EE,"h_seed_sigmaIphiIphi_EE","#sigma_{i#phii#phi}")

  h_seed_swissCross_EB = TH1F("h_seed_swissCross_EB","h_seed_swissCross",200,-1.,1.1)
  tree.Draw("pfCluster_full5x5_swissCross[superCluster_seedIndex]>>h_seed_swissCross_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_swissCross_EB,"h_seed_swissCross_EB","1.-e4/e1")

  h_seed_swissCross_EE = TH1F("h_seed_swissCross_EE","h_seed_swissCross",200,-1.,1.1)
  tree.Draw("pfCluster_full5x5_swissCross[superCluster_seedIndex]>>h_seed_swissCross_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_swissCross_EE,"h_seed_swissCross_EE","1.-e4/e1")
  
  h_seed_energy_EB = TH1F("h_seed_energy_EB","h_seed_energy",600,0.,300.)
  tree.Draw("pfCluster_energy[superCluster_seedIndex]>>h_seed_energy_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_energy_EB,"h_seed_energy_EB","Energy (GeV)")

  h_seed_energy_EE = TH1F("h_seed_energy_EE","h_seed_energy",600,0.,300.)
  tree.Draw("pfCluster_energy[superCluster_seedIndex]>>h_seed_energy_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_energy_EE,"h_seed_energy_EE","Energy (GeV)")

  h_seed_pt_EB = TH1F("h_seed_pt_EB","h_seed_pt",200,0.,100.)
  tree.Draw("pfCluster_energy[superCluster_seedIndex]/TMath::CosH(pfCluster_eta[superCluster_seedIndex])>>h_seed_pt_EB","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]==0")
  drawHisto(h_seed_pt_EB,"h_seed_pt_EB","Pt (GeV)")
 
  h_seed_pt_EE = TH1F("h_seed_pt_EE","h_seed_pt",200,0.,100.)
  tree.Draw("pfCluster_energy[superCluster_seedIndex]/TMath::CosH(pfCluster_eta[superCluster_seedIndex])>>h_seed_pt_EE","superCluster_seedIndex>=0 && pfCluster_iz[superCluster_seedIndex]!=0")
  drawHisto(h_seed_pt_EE,"h_seed_pt_EE","Pt (GeV)") 
  
  window_selection = " && superCluster_seedIndex>=0 && fabs(pfCluster_eta[superCluster_seedIndex]-pfCluster_eta)>0. && fabs(pfCluster_eta[superCluster_seedIndex]-pfCluster_eta)<0.2 && fabs(TVector2::Phi_mpi_pi(pfCluster_phi[superCluster_seedIndex]-pfCluster_phi))>0. && fabs(TVector2::Phi_mpi_pi(pfCluster_phi[superCluster_seedIndex]-pfCluster_phi))<0.6" 
  #window_selection = " && 1.>0."

  h_pfCluster_nXtals_caloMatched_EB = TH1F("h_pfCluster_nXtals_caloMatched_EB","h_pfCluster_nXtals",20,0.5,20.5)
  h_pfCluster_nXtals_caloUnmatched_EB = TH1F("h_pfCluster_nXtals_caloUnmatched_EB","h_pfCluster_nXtals",20,0.5,20.5)
  tree.Draw("pfCluster_nXtals>>h_pfCluster_nXtals_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_nXtals>>h_pfCluster_nXtals_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_nXtals_caloMatched_EB,h_pfCluster_nXtals_caloUnmatched_EB,"h_pfCluster_nXtals_EB","nXtals")

  h_pfCluster_nXtals_caloMatched_EE = TH1F("h_pfCluster_nXtals_caloMatched_EE","h_pfCluster_nXtals",20,0.5,20.5)
  h_pfCluster_nXtals_caloUnmatched_EE = TH1F("h_pfCluster_nXtals_caloUnmatched_EE","h_pfCluster_nXtals",20,0.5,20.5)
  tree.Draw("pfCluster_nXtals>>h_pfCluster_nXtals_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_nXtals>>h_pfCluster_nXtals_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_nXtals_caloMatched_EE,h_pfCluster_nXtals_caloUnmatched_EE,"h_pfCluster_nXtals_EE","nXtals") 

  h_pfCluster_R9_caloMatched_EB = TH1F("h_pfCluster_R9_caloMatched_EB","h_pfCluster_R9",80,0.6,1.4)
  h_pfCluster_R9_caloUnmatched_EB = TH1F("h_pfCluster_R9_caloUnmatched_EB","h_pfCluster_R9",80,0.6,1.4)
  tree.Draw("pfCluster_full5x5_r9>>h_pfCluster_R9_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_full5x5_r9>>h_pfCluster_R9_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_R9_caloMatched_EB,h_pfCluster_R9_caloUnmatched_EB,"h_pfCluster_R9_EB","R9")

  h_pfCluster_R9_caloMatched_EE = TH1F("h_pfCluster_R9_caloMatched_EE","h_pfCluster_R9",80,0.6,1.4)
  h_pfCluster_R9_caloUnmatched_EE = TH1F("h_pfCluster_R9_caloUnmatched_EE","h_pfCluster_R9",80,0.6,1.4)
  tree.Draw("pfCluster_full5x5_r9>>h_pfCluster_R9_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_full5x5_r9>>h_pfCluster_R9_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_R9_caloMatched_EE,h_pfCluster_R9_caloUnmatched_EE,"h_pfCluster_R9_EE","R9")

  h_pfCluster_etaWidth_caloMatched_EB = TH1F("h_pfCluster_etaWidth_caloMatched_EB","h_pfCluster_etaWidth",100,0.,0.016)
  h_pfCluster_etaWidth_caloUnmatched_EB = TH1F("h_pfCluster_etaWidth_caloUnmatched_EB","h_pfCluster_etaWidth",100,0.,0.016)
  tree.Draw("pfCluster_etaWidth>>h_pfCluster_etaWidth_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_etaWidth>>h_pfCluster_etaWidth_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_etaWidth_caloMatched_EB,h_pfCluster_etaWidth_caloUnmatched_EB,"h_pfCluster_etaWidth_EB","#eta-Width")

  h_pfCluster_etaWidth_caloMatched_EE = TH1F("h_pfCluster_etaWidth_caloMatched_EE","h_pfCluster_etaWidth",100,0.,0.06)
  h_pfCluster_etaWidth_caloUnmatched_EE = TH1F("h_pfCluster_etaWidth_caloUnmatched_EE","h_pfCluster_etaWidth",100,0.,0.06)
  tree.Draw("pfCluster_etaWidth>>h_pfCluster_etaWidth_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_etaWidth>>h_pfCluster_etaWidth_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_etaWidth_caloMatched_EE,h_pfCluster_etaWidth_caloUnmatched_EE,"h_pfCluster_etaWidth_EE","#eta-Width")

  h_pfCluster_phiWidth_caloMatched_EB = TH1F("h_pfCluster_phiWidth_caloMatched_EB","h_pfCluster_phiWidth",100,0.,0.018)
  h_pfCluster_phiWidth_caloUnmatched_EB = TH1F("h_pfCluster_phiWidth_caloUnmatched_EB","h_pfCluster_phiWidth",100,0.,0.018)
  tree.Draw("pfCluster_phiWidth>>h_pfCluster_phiWidth_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_phiWidth>>h_pfCluster_phiWidth_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_phiWidth_caloMatched_EB,h_pfCluster_phiWidth_caloUnmatched_EB,"h_pfCluster_phiWidth_EB","#phi-Width")

  h_pfCluster_phiWidth_caloMatched_EE = TH1F("h_pfCluster_phiWidth_caloMatched_EE","h_pfCluster_phiWidth",100,0.,0.06)
  h_pfCluster_phiWidth_caloUnmatched_EE = TH1F("h_pfCluster_phiWidth_caloUnmatched_EE","h_pfCluster_phiWidth",100,0.,0.06)
  tree.Draw("pfCluster_phiWidth>>h_pfCluster_phiWidth_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_phiWidth>>h_pfCluster_phiWidth_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_phiWidth_caloMatched_EE,h_pfCluster_phiWidth_caloUnmatched_EE,"h_pfCluster_phiWidth_EE","#phi-Width")
  
  h_pfCluster_sigmaIetaIeta_caloMatched_EB = TH1F("h_pfCluster_sigmaIetaIeta_caloMatched_EB","h_pfCluster_sigmaIetaIeta",100,0.,0.025)
  h_pfCluster_sigmaIetaIeta_caloUnmatched_EB = TH1F("h_pfCluster_sigmaIetaIeta_caloUnmatched_EB","h_pfCluster_sigmaIetaIeta",100,0.,0.025)
  tree.Draw("pfCluster_full5x5_sigmaIetaIeta>>h_pfCluster_sigmaIetaIeta_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_full5x5_sigmaIetaIeta>>h_pfCluster_sigmaIetaIeta_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_sigmaIetaIeta_caloMatched_EB,h_pfCluster_sigmaIetaIeta_caloUnmatched_EB,"h_pfCluster_sigmaIetaIeta_EB","#sigma_{i#etai#eta}")

  h_pfCluster_sigmaIetaIeta_caloMatched_EE = TH1F("h_pfCluster_sigmaIetaIeta_caloMatched_EE","h_pfCluster_sigmaIetaIeta",300,0.,0.08)
  h_pfCluster_sigmaIetaIeta_caloUnmatched_EE = TH1F("h_pfCluster_sigmaIetaIeta_caloUnmatched_EE","h_pfCluster_sigmaIetaIeta",300,0.,0.08)
  tree.Draw("pfCluster_full5x5_sigmaIetaIeta>>h_pfCluster_sigmaIetaIeta_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_full5x5_sigmaIetaIeta>>h_pfCluster_sigmaIetaIeta_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_sigmaIetaIeta_caloMatched_EE,h_pfCluster_sigmaIetaIeta_caloUnmatched_EE,"h_pfCluster_sigmaIetaIeta_EE","#sigma_{i#etai#eta}")
  
  h_pfCluster_sigmaIetaIphi_caloMatched_EB = TH1F("h_pfCluster_sigmaIetaIphi_caloMatched_EB","h_pfCluster_sigmaIetaIphi",100,-0.001,0.001)
  h_pfCluster_sigmaIetaIphi_caloUnmatched_EB = TH1F("h_pfCluster_sigmaIetaIphi_caloUnmatched_EB","h_pfCluster_sigmaIetaIphi",100,-0.001,0.001)
  tree.Draw("pfCluster_full5x5_sigmaIetaIphi>>h_pfCluster_sigmaIetaIphi_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_full5x5_sigmaIetaIphi>>h_pfCluster_sigmaIetaIphi_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_sigmaIetaIphi_caloMatched_EB,h_pfCluster_sigmaIetaIphi_caloUnmatched_EB,"h_pfCluster_sigmaIetaIphi_EB","#sigma_{i#etai#phi}")
  
  h_pfCluster_sigmaIetaIphi_caloMatched_EE = TH1F("h_pfCluster_sigmaIetaIphi_caloMatched_EE","h_pfCluster_sigmaIetaIphi",100,-0.005,0.005)
  h_pfCluster_sigmaIetaIphi_caloUnmatched_EE = TH1F("h_pfCluster_sigmaIetaIphi_caloUnmatched_EE","h_pfCluster_sigmaIetaIphi",100,-0.005,0.005)
  tree.Draw("pfCluster_full5x5_sigmaIetaIphi>>h_pfCluster_sigmaIetaIphi_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_full5x5_sigmaIetaIphi>>h_pfCluster_sigmaIetaIphi_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_sigmaIetaIphi_caloMatched_EE,h_pfCluster_sigmaIetaIphi_caloUnmatched_EE,"h_pfCluster_sigmaIetaIphi_EE","#sigma_{i#etai#phi}")

  h_pfCluster_sigmaIphiIphi_caloMatched_EB = TH1F("h_pfCluster_sigmaIphiIphi_caloMatched_EB","h_pfCluster_sigmaIphiIphi",200,0.,0.05)
  h_pfCluster_sigmaIphiIphi_caloUnmatched_EB = TH1F("h_pfCluster_sigmaIphiIphi_caloUnmatched_EB","h_pfCluster_sigmaIphiIphi",200,0.,0.05)
  tree.Draw("pfCluster_full5x5_sigmaIphiIphi>>h_pfCluster_sigmaIphiIphi_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_full5x5_sigmaIphiIphi>>h_pfCluster_sigmaIphiIphi_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_sigmaIphiIphi_caloMatched_EB,h_pfCluster_sigmaIphiIphi_caloUnmatched_EB,"h_pfCluster_sigmaIphiIphi_EB","#sigma_{i#phii#phi}")
  
  h_pfCluster_sigmaIphiIphi_caloMatched_EE = TH1F("h_pfCluster_sigmaIphiIphi_caloMatched_EE","h_pfCluster_sigmaIphiIphi",400,0.,0.08)
  h_pfCluster_sigmaIphiIphi_caloUnmatched_EE = TH1F("h_pfCluster_sigmaIphiIphi_caloUnmatched_EE","h_pfCluster_sigmaIphiIphi",400,0.,0.08)
  tree.Draw("pfCluster_full5x5_sigmaIphiIphi>>h_pfCluster_sigmaIphiIphi_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_full5x5_sigmaIphiIphi>>h_pfCluster_sigmaIphiIphi_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_sigmaIphiIphi_caloMatched_EE,h_pfCluster_sigmaIphiIphi_caloUnmatched_EE,"h_pfCluster_sigmaIphiIphi_EE","#sigma_{i#phii#phi}")
  
  h_pfCluster_swissCross_caloMatched_EB = TH1F("h_pfCluster_swissCross_caloMatched_EB","h_pfCluster_swissCross",200,-1.,1.1)
  h_pfCluster_swissCross_caloUnmatched_EB = TH1F("h_pfCluster_swissCross_caloUnmatched_EB","h_pfCluster_swissCross",200,-1.,1.1)
  tree.Draw("pfCluster_full5x5_swissCross>>h_pfCluster_swissCross_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_full5x5_swissCross>>h_pfCluster_swissCross_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_swissCross_caloMatched_EB,h_pfCluster_swissCross_caloUnmatched_EB,"h_pfCluster_swissCross_EB","1.-e4/e1")

  h_pfCluster_swissCross_caloMatched_EE = TH1F("h_pfCluster_swissCross_caloMatched_EE","h_pfCluster_swissCross",200,-1.,1.1)
  h_pfCluster_swissCross_caloUnmatched_EE = TH1F("h_pfCluster_swissCross_caloUnmatched_EE","h_pfCluster_swissCross",200,-1.,1.1)
  tree.Draw("pfCluster_full5x5_swissCross>>h_pfCluster_swissCross_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_full5x5_swissCross>>h_pfCluster_swissCross_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_swissCross_caloMatched_EE,h_pfCluster_swissCross_caloUnmatched_EE,"h_pfCluster_swissCross_EE","1.-e4/e1")

  h_pfCluster_energy_caloMatched_EB = TH1F("h_pfCluster_energy_caloMatched_EB","h_pfCluster_energy",600,0.,300.)
  h_pfCluster_energy_caloUnmatched_EB = TH1F("h_pfCluster_energy_caloUnmatched_EB","h_pfCluster_energy",600,0.,300.)
  tree.Draw("pfCluster_energy>>h_pfCluster_energy_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_energy>>h_pfCluster_energy_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_energy_caloMatched_EB,h_pfCluster_energy_caloUnmatched_EB,"h_pfCluster_energy_EB","Energy (GeV)")

  h_pfCluster_energy_caloMatched_EE = TH1F("h_pfCluster_energy_caloMatched_EE","h_pfCluster_energy",600,0.,300.)
  h_pfCluster_energy_caloUnmatched_EE = TH1F("h_pfCluster_energy_caloUnmatched_EE","h_pfCluster_energy",600,0.,300.)
  tree.Draw("pfCluster_energy>>h_pfCluster_energy_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_energy>>h_pfCluster_energy_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_energy_caloMatched_EE,h_pfCluster_energy_caloUnmatched_EE,"h_pfCluster_energy_EE","Energy (GeV)")

  h_pfCluster_pt_caloMatched_EB = TH1F("h_pfCluster_pt_caloMatched_EB","h_pfCluster_pt",200,0.,100.)
  h_pfCluster_pt_caloUnmatched_EB = TH1F("h_pfCluster_pt_caloUnmatched_EB","h_pfCluster_pt",200,0.,100.)
  tree.Draw("pfCluster_energy/TMath::CosH(pfCluster_eta)>>h_pfCluster_pt_caloMatched_EB","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz==0"+window_selection)
  tree.Draw("pfCluster_energy/TMath::CosH(pfCluster_eta)>>h_pfCluster_pt_caloUnmatched_EB","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz==0"+window_selection)
  drawHistos(h_pfCluster_pt_caloMatched_EB,h_pfCluster_pt_caloUnmatched_EB,"h_pfCluster_pt_EB","Pt (GeV)")
 
  h_pfCluster_pt_caloMatched_EE = TH1F("h_pfCluster_pt_caloMatched_EE","h_pfCluster_pt",200,0.,100.)
  h_pfCluster_pt_caloUnmatched_EE = TH1F("h_pfCluster_pt_caloUnmatched_EE","h_pfCluster_pt",200,0.,100.)
  tree.Draw("pfCluster_energy/TMath::CosH(pfCluster_eta)>>h_pfCluster_pt_caloMatched_EE","pfCluster_simScore_MatchedIndex>=0 && pfCluster_iz!=0"+window_selection)
  tree.Draw("pfCluster_energy/TMath::CosH(pfCluster_eta)>>h_pfCluster_pt_caloUnmatched_EE","pfCluster_simScore_MatchedIndex<0 && pfCluster_iz!=0"+window_selection)
  drawHistos(h_pfCluster_pt_caloMatched_EE,h_pfCluster_pt_caloUnmatched_EE,"h_pfCluster_pt_EE","Pt (GeV)") 


