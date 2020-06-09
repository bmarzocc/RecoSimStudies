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
   #h1.Scale(1./h1.Integral())
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

  #inFile = '/eos/cms/store/user/bmarzocc/Clustering/DeepCluster_FourGammasGun_RECO.root'
  #inFile = '/eos/cms/store/user/bmarzocc/Clustering/Mustache_FourGammasGun_RECO.root'
  #inFile = '/eos/cms/store/user/bmarzocc/Clustering/DeepCluster_FourElectronsGun_RECO.root'
  inFile = '/eos/cms/store/user/bmarzocc/Clustering/Mustache_FourElectronsGun_RECO.root'
 
  print "inputFile = ",inFile
  
  f_in = TFile(inFile)
  tree = f_in.Get('Events')

  print "Filling SuperClusters..."

  h_superCluster_Eta_EB = TH1F("h_superCluster_Eta_EB","h_superCluster_Eta_EB",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALBarrel_RECO.obj.eta()>>h_superCluster_Eta_EB")
  drawHisto(h_superCluster_Eta_EB,"h_superCluster_Eta_EB","#eta")

  h_superCluster_Eta_EE = TH1F("h_superCluster_Eta_EE","h_superCluster_Eta_EE",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALEndcapWithPreshower_RECO.obj.eta()>>h_superCluster_Eta_EE") 
  drawHisto(h_superCluster_Eta_EE,"h_superCluster_Eta_EE","#eta")

  h_superCluster_Phi_EB = TH1F("h_superCluster_Phi_EB","h_superCluster_Phi_EB",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALBarrel_RECO.obj.phi()>>h_superCluster_Phi_EB")
  drawHisto(h_superCluster_Phi_EB,"h_superCluster_Phi_EB","#phi")

  h_superCluster_Phi_EE = TH1F("h_superCluster_Phi_EE","h_superCluster_Phi_EE",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALEndcapWithPreshower_RECO.obj.phi()>>h_superCluster_Phi_EE")
  drawHisto(h_superCluster_Phi_EE,"h_superCluster_Phi_EE","#phi")

  h_superCluster_Energy_EB = TH1F("h_superCluster_Energy_EB","h_superCluster_Energy_EB",60,0.,300.)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALBarrel_RECO.obj.energy()>>h_superCluster_Energy_EB")
  drawHisto(h_superCluster_Energy_EB,"h_superCluster_Energy_EB","Energy (GeV)")

  h_superCluster_Energy_EE = TH1F("h_superCluster_Energy_EE","h_superCluster_Energy_EE",200,0.,1000.)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALEndcapWithPreshower_RECO.obj.energy()>>h_superCluster_Energy_EE")
  drawHisto(h_superCluster_Energy_EE,"h_superCluster_Energy_EE","Energy (GeV)")

  h_superCluster_Pt_EB = TH1F("h_superCluster_Pt_EB","h_superCluster_Pt_EB",30,0.,150.)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALBarrel_RECO.obj.energy() / TMath::CosH( recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALBarrel_RECO.obj.eta()) >> h_superCluster_Pt_EB")
  drawHisto(h_superCluster_Pt_EB,"h_superCluster_Pt_EB","Pt (GeV)")

  h_superCluster_Pt_EE = TH1F("h_superCluster_Pt_EE","h_superCluster_Pt_EE",30,0.,150.)
  tree.Draw("recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALEndcapWithPreshower_RECO.obj.energy() / TMath::CosH( recoSuperClusters_particleFlowSuperClusterECAL_particleFlowSuperClusterECALEndcapWithPreshower_RECO.obj.eta()) >> h_superCluster_Pt_EE")
  drawHisto(h_superCluster_Pt_EE,"h_superCluster_Pt_EE","Pt (GeV)")

  print "Filling RefinedSCs..."
  
  h_refinedSC_Eta_EB = TH1F("h_refinedSC_Eta_EB","h_refinedSC_Eta_EB",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.eta()>>h_refinedSC_Eta_EB","fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())<1.479")
  drawHisto(h_refinedSC_Eta_EB,"h_refinedSC_Eta_EB","#eta")

  h_refinedSC_Eta_EE = TH1F("h_refinedSC_Eta_EE","h_refinedSC_Eta_EE",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.eta()>>h_refinedSC_Eta_EE","fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())>1.479")
  drawHisto(h_refinedSC_Eta_EE,"h_refinedSC_Eta_EE","#eta")

  h_refinedSC_Phi_EB = TH1F("h_refinedSC_Phi_EB","h_refinedSC_Phi_EB",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.phi()>>h_refinedSC_Phi_EB","fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())<1.479")
  drawHisto(h_refinedSC_Phi_EB,"h_refinedSC_Phi_EB","#phi")

  h_refinedSC_Phi_EE = TH1F("h_refinedSC_Phi_EE","h_refinedSC_Phi_EE",70,-3.5,3.5)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.phi()>>h_refinedSC_Phi_EE","fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())>1.479")
  drawHisto(h_refinedSC_Phi_EE,"h_refinedSC_Phi_EE","#phi")

  h_refinedSC_Energy_EB = TH1F("h_refinedSC_Energy_EB","h_refinedSC_Energy_EB",60,0.,300.)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.energy()>>h_refinedSC_Energy_EB","fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())<1.479")
  drawHisto(h_refinedSC_Energy_EB,"h_refinedSC_Energy_EB","Energy (GeV)")

  h_refinedSC_Energy_EE = TH1F("h_refinedSC_Energy_EE","h_refinedSC_Energy_EE",200,0.,1000.)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.energy()>>h_refinedSC_Energy_EE","fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())>1.479")
  drawHisto(h_refinedSC_Energy_EE,"h_refinedSC_Energy_EE","Energy (GeV)")

  h_refinedSC_Pt_EB = TH1F("h_refinedSC_Pt_EB","h_refinedSC_Pt_EB",30,0.,150.)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.energy() / TMath::CosH(recoSuperClusters_particleFlowEGamma__RECO.obj.eta()) >> h_refinedSC_Pt_EB", "fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())<1.479")
  drawHisto(h_refinedSC_Pt_EB,"h_refinedSC_Pt_EB","Pt (GeV)")

  h_refinedSC_Pt_EE = TH1F("h_refinedSC_Pt_EE","h_refinedSC_Pt_EE",30,0.,150.)
  tree.Draw("recoSuperClusters_particleFlowEGamma__RECO.obj.energy() / TMath::CosH(recoSuperClusters_particleFlowEGamma__RECO.obj.eta()) >> h_refinedSC_Pt_EE", "fabs(recoSuperClusters_particleFlowEGamma__RECO.obj.eta())>1.479")
  drawHisto(h_refinedSC_Pt_EE,"h_refinedSC_Pt_EE","Pt (GeV)")

  if "Gammas" in inFile: 
     print "Filling gedPhotons..."
  
     h_photon_Eta_EB = TH1F("h_photon_Eta_EB","h_photon_Eta_EB",70,-3.5,3.5)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.eta()>>h_photon_Eta_EB","fabs(recoPhotons_gedPhotons__RECO.obj.eta())<1.479")
     drawHisto(h_photon_Eta_EB,"h_photon_Eta_EB","#eta")

     h_photon_Eta_EE = TH1F("h_photon_Eta_EE","h_photon_Eta_EE",70,-3.5,3.5)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.eta()>>h_photon_Eta_EE","fabs(recoPhotons_gedPhotons__RECO.obj.eta())>1.479")
     drawHisto(h_photon_Eta_EE,"h_photon_Eta_EE","#eta")

     h_photon_Phi_EB = TH1F("h_photon_Phi_EB","h_photon_Phi_EB",70,-3.5,3.5)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.phi()>>h_photon_Phi_EB","fabs(recoPhotons_gedPhotons__RECO.obj.eta())<1.479")
     drawHisto(h_photon_Phi_EB,"h_photon_Phi_EB","#phi")

     h_photon_Phi_EE = TH1F("h_photon_Phi_EE","h_photon_Phi_EE",70,-3.5,3.5)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.phi()>>h_photon_Phi_EE","fabs(recoPhotons_gedPhotons__RECO.obj.eta())>1.479")
     drawHisto(h_photon_Phi_EE,"h_photon_Phi_EE","#phi")

     h_photon_Energy_EB = TH1F("h_photon_Energy_EB","h_photon_Energy_EB",60,0.,300.)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.energy()>>h_photon_Energy_EB","fabs(recoPhotons_gedPhotons__RECO.obj.eta())<1.479")
     drawHisto(h_photon_Energy_EB,"h_photon_Energy_EB","Energy (GeV)")

     h_photon_Energy_EE = TH1F("h_photon_Energy_EE","h_photon_Energy_EE",200,0.,1000.)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.energy()>>h_photon_Energy_EE","fabs(recoPhotons_gedPhotons__RECO.obj.eta())>1.479")
     drawHisto(h_photon_Energy_EE,"h_photon_Energy_EE","Energy (GeV)")

     h_photon_Pt_EB = TH1F("h_photon_Pt_EB","h_photon_Pt_EB",30,0.,150.)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.energy() / TMath::CosH(recoPhotons_gedPhotons__RECO.obj.eta()) >> h_photon_Pt_EB", "fabs(recoPhotons_gedPhotons__RECO.obj.eta())<1.479")
     drawHisto(h_photon_Pt_EB,"h_photon_Pt_EB","Pt (GeV)")

     h_photon_Pt_EE = TH1F("h_photon_Pt_EE","h_photon_Pt_EE",30,0.,150.)
     tree.Draw("recoPhotons_gedPhotons__RECO.obj.energy() / TMath::CosH(recoPhotons_gedPhotons__RECO.obj.eta()) >> h_photon_Pt_EE", "fabs(recoPhotons_gedPhotons__RECO.obj.eta())>1.479")
     drawHisto(h_photon_Pt_EE,"h_photon_Pt_EE","Pt (GeV)")

  if "Electrons" in inFile: 
     print "Filling gedGsfElectrons..."
  
     h_electron_Eta_EB = TH1F("h_electron_Eta_EB","h_electron_Eta_EB",70,-3.5,3.5)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.eta()>>h_electron_Eta_EB","fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())<1.479")
     drawHisto(h_electron_Eta_EB,"h_electron_Eta_EB","#eta")

     h_electron_Eta_EE = TH1F("h_electron_Eta_EE","h_electron_Eta_EE",70,-3.5,3.5)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.eta()>>h_electron_Eta_EE","fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())>1.479")
     drawHisto(h_electron_Eta_EE,"h_electron_Eta_EE","#eta")

     h_electron_Phi_EB = TH1F("h_electron_Phi_EB","h_electron_Phi_EB",70,-3.5,3.5)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.phi()>>h_electron_Phi_EB","fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())<1.479")
     drawHisto(h_electron_Phi_EB,"h_electron_Phi_EB","#phi")

     h_electron_Phi_EE = TH1F("h_electron_Phi_EE","h_electron_Phi_EE",70,-3.5,3.5)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.phi()>>h_electron_Phi_EE","fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())>1.479")
     drawHisto(h_electron_Phi_EE,"h_electron_Phi_EE","#phi")

     h_electron_Energy_EB = TH1F("h_electron_Energy_EB","h_electron_Energy_EB",60,0.,300.)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.energy()>>h_electron_Energy_EB","fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())<1.479")
     drawHisto(h_electron_Energy_EB,"h_electron_Energy_EB","Energy (GeV)")

     h_electron_Energy_EE = TH1F("h_electron_Energy_EE","h_electron_Energy_EE",200,0.,1000.)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.energy()>>h_electron_Energy_EE","fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())>1.479")
     drawHisto(h_electron_Energy_EE,"h_electron_Energy_EE","Energy (GeV)")

     h_electron_Pt_EB = TH1F("h_electron_Pt_EB","h_electron_Pt_EB",30,0.,150.)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.energy() / TMath::CosH(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta()) >> h_electron_Pt_EB", "fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())<1.479")
     drawHisto(h_electron_Pt_EB,"h_electron_Pt_EB","Pt (GeV)")

     h_electron_Pt_EE = TH1F("h_electron_Pt_EE","h_electron_Pt_EE",30,0.,150.)
     tree.Draw("recoGsfElectrons_gedGsfElectrons__RECO.obj.energy() / TMath::CosH(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta()) >> h_electron_Pt_EE", "fabs(recoGsfElectrons_gedGsfElectrons__RECO.obj.eta())>1.479")
     drawHisto(h_electron_Pt_EE,"h_electron_Pt_EE","Pt (GeV)")  

  
  
