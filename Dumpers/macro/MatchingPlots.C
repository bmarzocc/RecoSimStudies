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
#include <algorithm> 
#include <iostream>
#include "TStyle.h"

void drawHisto(TH1F* h_dR_simScore, TH1F* h_n_shared_xtals, TH1F* h_sim_fraction, TH1F* h_sim_fraction_min1, std::string Name);
void drawEfficiency(TEfficiency* eff_dR_simScore, TEfficiency* eff_n_shared_xtals, TEfficiency* eff_sim_fraction, TEfficiency* eff_sim_fraction_min1, std::string Name);

void MatchingPlots(){

   //gStyle->SetOptStat(1110); 
   gStyle->SetOptStat(0000); 

   TFile* inFile = TFile::Open("../output.root");
   if (inFile == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open ../output.root !\n");
      return;
   }

   TH1F* h_PFCluster_EoverEtrue_dR_simScore_EB = new TH1F("h_PFCluster_EoverEtrue_dR_simScore_EB","PFCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_n_shared_xtals_EB = new TH1F("h_PFCluster_EoverEtrue_n_shared_xtals_EB","PFCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_sim_fraction_EB = new TH1F("h_PFCluster_EoverEtrue_sim_fraction_EB","PFCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_sim_fraction_min1_EB = new TH1F("h_PFCluster_EoverEtrue_sim_fraction_min1_EB","PFCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_dR_simScore_EE = new TH1F("h_PFCluster_EoverEtrue_dR_simScore_EE","PFCluster_EoverEtrue_EE",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_n_shared_xtals_EE = new TH1F("h_PFCluster_EoverEtrue_n_shared_xtals_EE","PFCluster_EoverEtrue_EE",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_sim_fraction_EE = new TH1F("h_PFCluster_EoverEtrue_sim_fraction_EE","PFCluster_EoverEtrue_EE",25,0.,1.5);
   TH1F* h_PFCluster_EoverEtrue_sim_fraction_min1_EE = new TH1F("h_PFCluster_EoverEtrue_sim_fraction_min1_EE","PFCluster_EoverEtrue_EE",25,0.,1.5);

   TH1F* h_GenParticle_Denum_vs_Eta = new TH1F("h_GenParticle_Denum_vs_Eta","",30,-3.,3.);
   TH1F* h_CaloParticle_Denum_vs_Eta = new TH1F("h_CaloParticle_Denum_vs_Eta","",30,-3.,3.);
   TH1F* h_PFCluster_dR_simScore_Num_vs_Eta = new TH1F("h_PFCluster_dR_simScore_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_PFCluster_n_shared_xtals_Num_vs_Eta = new TH1F("h_PFCluster_n_shared_xtals_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_PFCluster_sim_fraction_Num_vs_Eta = new TH1F("h_PFCluster_sim_fraction_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_PFCluster_sim_fraction_min1_Num_vs_Eta = new TH1F("h_PFCluster_sim_fraction_min1_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_SuperCluster_dR_simScore_Num_vs_Eta = new TH1F("h_SuperCluster_dR_simScore_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_SuperCluster_n_shared_xtals_Num_vs_Eta = new TH1F("h_SuperCluster_n_shared_xtals_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_SuperCluster_sim_fraction_Num_vs_Eta = new TH1F("h_SuperCluster_sim_fraction_Num_vs_Eta","",30,-3.,3.);
   TH1F* h_SuperCluster_sim_fraction_min1_Num_vs_Eta = new TH1F("h_SuperCluster_sim_fraction_min1_Num_vs_Eta","",30,-3.,3.);
   
   TH1F* h_SuperCluster_EoverEtrue_dR_simScore_EB = new TH1F("h_SuperCluster_EoverEtrue_dR_simScore_EB","SuperCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_n_shared_xtals_EB = new TH1F("h_SuperCluster_EoverEtrue_n_shared_xtals_EB","SuperCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_sim_fraction_EB = new TH1F("h_SuperCluster_EoverEtrue_sim_fraction_EB","SuperCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_sim_fraction_min1_EB = new TH1F("h_SuperCluster_EoverEtrue_sim_fraction_min1_EB","SuperCluster_EoverEtrue_EB",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_dR_simScore_EE = new TH1F("h_SuperCluster_EoverEtrue_dR_simScore_EE","SuperCluster_EoverEtrue_EE",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_n_shared_xtals_EE = new TH1F("h_SuperCluster_EoverEtrue_n_shared_xtals_EE","SuperCluster_EoverEtrue_EE",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_sim_fraction_EE = new TH1F("h_SuperCluster_EoverEtrue_sim_fraction_EE","SuperCluster_EoverEtrue_EE",25,0.,1.5);
   TH1F* h_SuperCluster_EoverEtrue_sim_fraction_min1_EE = new TH1F("h_SuperCluster_EoverEtrue_sim_fraction_min1_EE","SuperCluster_EoverEtrue_EE",25,0.,1.5);

   std::map<int,float> simEnergy;
   std::map<int,float> PFCluster_RecoEnergy_dR_simScore;
   std::map<int,int> PFCluster_nSimMatches_dR_simScore;
   std::map<int,float> PFCluster_RecoEnergy_n_shared_xtals;
   std::map<int,int> PFCluster_nSimMatches_n_shared_xtals;
   std::map<int,float> PFCluster_RecoEnergy_sim_fraction;
   std::map<int,int> PFCluster_nSimMatches_sim_fraction;
   std::map<int,float> PFCluster_RecoEnergy_sim_fraction_min1;
   std::map<int,int> PFCluster_nSimMatches_sim_fraction_min1;
   std::map<int,float> SuperCluster_RecoEnergy_dR_simScore;
   std::map<int,int> SuperCluster_nSimMatches_dR_simScore;
   std::map<int,float> SuperCluster_RecoEnergy_n_shared_xtals;
   std::map<int,int> SuperCluster_nSimMatches_n_shared_xtals;
   std::map<int,float> SuperCluster_RecoEnergy_sim_fraction;
   std::map<int,int> SuperCluster_nSimMatches_sim_fraction;
   std::map<int,float> SuperCluster_RecoEnergy_sim_fraction_min1;
   std::map<int,int> SuperCluster_nSimMatches_sim_fraction_min1;

   vector<float> genParticle_eta_tmp;
   vector<float> caloParticle_simEnergy_tmp;
   vector<float> caloParticle_simEta_tmp;
   vector<float> pfCluster_energy_tmp;
   vector<int> pfCluster_dR_simScore_MatchedIndex_tmp;
   vector<int> pfCluster_n_shared_xtals_MatchedIndex_tmp;
   vector<int> pfCluster_sim_fraction_MatchedIndex_tmp;
   vector<int> pfCluster_sim_fraction_min1_MatchedIndex_tmp;
   vector<float> superCluster_energy_tmp;
   vector<int> superCluster_dR_simScore_MatchedIndex_tmp;
   vector<int> superCluster_n_shared_xtals_MatchedIndex_tmp;
   vector<int> superCluster_sim_fraction_MatchedIndex_tmp;
   vector<int> superCluster_sim_fraction_min1_MatchedIndex_tmp;
   
   // Create tyhe tree reader and its data containers
   TTreeReader myReader("deepclusteringdumper/caloTree", inFile); 
   TTreeReaderValue<vector<float>> genParticle_eta(myReader, "genParticle_eta");	
   TTreeReaderValue<vector<float>> caloParticle_simEnergy(myReader, "caloParticle_simEnergy");
   TTreeReaderValue<vector<float>> caloParticle_simEta(myReader, "caloParticle_simEta");
   TTreeReaderValue<vector<float>> pfCluster_energy(myReader, "pfCluster_energy");
   TTreeReaderValue<vector<int>> pfCluster_dR_simScore_MatchedIndex(myReader, "pfCluster_dR_simScore_MatchedIndex");
   TTreeReaderValue<vector<int>> pfCluster_n_shared_xtals_MatchedIndex(myReader, "pfCluster_n_shared_xtals_MatchedIndex");
   TTreeReaderValue<vector<int>> pfCluster_sim_fraction_MatchedIndex(myReader, "pfCluster_sim_fraction_MatchedIndex"); 
   TTreeReaderValue<vector<int>> pfCluster_sim_fraction_min1_MatchedIndex(myReader, "pfCluster_sim_fraction_min1_MatchedIndex");
   TTreeReaderValue<vector<float>> superCluster_energy(myReader, "superCluster_energy");
   TTreeReaderValue<vector<int>> superCluster_dR_simScore_MatchedIndex(myReader, "superCluster_dR_simScore_MatchedIndex");
   TTreeReaderValue<vector<int>> superCluster_n_shared_xtals_MatchedIndex(myReader, "superCluster_n_shared_xtals_MatchedIndex");
   TTreeReaderValue<vector<int>> superCluster_sim_fraction_MatchedIndex(myReader, "superCluster_sim_fraction_MatchedIndex"); 
   TTreeReaderValue<vector<int>> superCluster_sim_fraction_min1_MatchedIndex(myReader, "superCluster_sim_fraction_min1_MatchedIndex");
   
   // Loop over all entries of the TTree or TChain.
   int  entry = 0;
   while (myReader.Next()) {
       
      if(entry%10==0) std::cout<<"--- Reading tree = "<<entry<<std::endl;
      //if(entry>1) continue;

      genParticle_eta_tmp = *genParticle_eta;
      caloParticle_simEnergy_tmp = *caloParticle_simEnergy;
      caloParticle_simEta_tmp = *caloParticle_simEta;
      pfCluster_energy_tmp = *pfCluster_energy; 
      superCluster_energy_tmp = *superCluster_energy; 

      pfCluster_dR_simScore_MatchedIndex_tmp = *pfCluster_dR_simScore_MatchedIndex; 
      pfCluster_n_shared_xtals_MatchedIndex_tmp = *pfCluster_n_shared_xtals_MatchedIndex; 
      pfCluster_sim_fraction_MatchedIndex_tmp = *pfCluster_sim_fraction_MatchedIndex; 
      pfCluster_sim_fraction_min1_MatchedIndex_tmp = *pfCluster_sim_fraction_min1_MatchedIndex; 
      superCluster_dR_simScore_MatchedIndex_tmp = *superCluster_dR_simScore_MatchedIndex; 
      superCluster_n_shared_xtals_MatchedIndex_tmp = *superCluster_n_shared_xtals_MatchedIndex; 
      superCluster_sim_fraction_MatchedIndex_tmp = *superCluster_sim_fraction_MatchedIndex; 
      superCluster_sim_fraction_min1_MatchedIndex_tmp = *superCluster_sim_fraction_min1_MatchedIndex; 

      for(unsigned int iGen=0; iGen<genParticle_eta_tmp.size(); iGen++)
          h_GenParticle_Denum_vs_Eta->Fill(genParticle_eta_tmp.at(iGen));
      
      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){
          h_CaloParticle_Denum_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
          simEnergy[iCalo] = caloParticle_simEnergy_tmp.at(iCalo);  
          PFCluster_nSimMatches_dR_simScore[iCalo]=0;
          PFCluster_nSimMatches_n_shared_xtals[iCalo]=0;
          PFCluster_nSimMatches_sim_fraction[iCalo]=0;
          PFCluster_nSimMatches_sim_fraction_min1[iCalo]=0;
      }
       
      for(unsigned int iPF=0; iPF<pfCluster_dR_simScore_MatchedIndex_tmp.size(); iPF++)
          if(pfCluster_dR_simScore_MatchedIndex_tmp.at(iPF)>=0){   
             PFCluster_RecoEnergy_dR_simScore[pfCluster_dR_simScore_MatchedIndex_tmp.at(iPF)]+=pfCluster_energy_tmp.at(iPF);
             PFCluster_nSimMatches_dR_simScore[pfCluster_dR_simScore_MatchedIndex_tmp.at(iPF)]++; 
          }

      for(unsigned int iPF=0; iPF<pfCluster_n_shared_xtals_MatchedIndex_tmp.size(); iPF++)
          if(pfCluster_n_shared_xtals_MatchedIndex_tmp.at(iPF)>=0){   
             PFCluster_RecoEnergy_n_shared_xtals[pfCluster_n_shared_xtals_MatchedIndex_tmp.at(iPF)]+=pfCluster_energy_tmp.at(iPF);
             PFCluster_nSimMatches_n_shared_xtals[pfCluster_n_shared_xtals_MatchedIndex_tmp.at(iPF)]++; 
          }

      for(unsigned int iPF=0; iPF<pfCluster_sim_fraction_MatchedIndex_tmp.size(); iPF++)
          if(pfCluster_sim_fraction_MatchedIndex_tmp.at(iPF)>=0){   
             PFCluster_RecoEnergy_sim_fraction[pfCluster_sim_fraction_MatchedIndex_tmp.at(iPF)]+=pfCluster_energy_tmp.at(iPF);
             PFCluster_nSimMatches_sim_fraction[pfCluster_sim_fraction_MatchedIndex_tmp.at(iPF)]++;
          }  

      for(unsigned int iPF=0; iPF<pfCluster_sim_fraction_min1_MatchedIndex_tmp.size(); iPF++)
          if(pfCluster_sim_fraction_min1_MatchedIndex_tmp.at(iPF)>=0){   
             PFCluster_RecoEnergy_sim_fraction_min1[pfCluster_sim_fraction_min1_MatchedIndex_tmp.at(iPF)]+=pfCluster_energy_tmp.at(iPF); 
             PFCluster_nSimMatches_sim_fraction_min1[pfCluster_sim_fraction_min1_MatchedIndex_tmp.at(iPF)]++;
          }

      for(unsigned int iPF=0; iPF<superCluster_dR_simScore_MatchedIndex_tmp.size(); iPF++)
          if(superCluster_dR_simScore_MatchedIndex_tmp.at(iPF)>=0){   
             SuperCluster_RecoEnergy_dR_simScore[superCluster_dR_simScore_MatchedIndex_tmp.at(iPF)]+=superCluster_energy_tmp.at(iPF);
             SuperCluster_nSimMatches_dR_simScore[superCluster_dR_simScore_MatchedIndex_tmp.at(iPF)]++;
          }

      for(unsigned int iPF=0; iPF<superCluster_n_shared_xtals_MatchedIndex_tmp.size(); iPF++)
          if(superCluster_n_shared_xtals_MatchedIndex_tmp.at(iPF)>=0){   
             SuperCluster_RecoEnergy_n_shared_xtals[superCluster_n_shared_xtals_MatchedIndex_tmp.at(iPF)]+=superCluster_energy_tmp.at(iPF);
             SuperCluster_nSimMatches_n_shared_xtals[superCluster_n_shared_xtals_MatchedIndex_tmp.at(iPF)]++;
          }  

      for(unsigned int iPF=0; iPF<superCluster_sim_fraction_MatchedIndex_tmp.size(); iPF++)
          if(superCluster_sim_fraction_MatchedIndex_tmp.at(iPF)>=0){   
             SuperCluster_RecoEnergy_sim_fraction[superCluster_sim_fraction_MatchedIndex_tmp.at(iPF)]+=superCluster_energy_tmp.at(iPF); 
             SuperCluster_nSimMatches_sim_fraction[superCluster_sim_fraction_MatchedIndex_tmp.at(iPF)]++;
          }

      for(unsigned int iPF=0; iPF<superCluster_sim_fraction_min1_MatchedIndex_tmp.size(); iPF++)
          if(superCluster_sim_fraction_min1_MatchedIndex_tmp.at(iPF)>=0){   
             SuperCluster_RecoEnergy_sim_fraction_min1[superCluster_sim_fraction_min1_MatchedIndex_tmp.at(iPF)]+=superCluster_energy_tmp.at(iPF);
             SuperCluster_nSimMatches_sim_fraction_min1[superCluster_sim_fraction_min1_MatchedIndex_tmp.at(iPF)]++;   
          }  
                   
      for(unsigned int iCalo=0; iCalo<caloParticle_simEnergy_tmp.size(); iCalo++){
          if(fabs(caloParticle_simEta_tmp.at(iCalo))<1.442){
            h_PFCluster_EoverEtrue_dR_simScore_EB->Fill(PFCluster_RecoEnergy_dR_simScore[iCalo]/simEnergy[iCalo]);
            h_PFCluster_EoverEtrue_n_shared_xtals_EB->Fill(PFCluster_RecoEnergy_n_shared_xtals[iCalo]/simEnergy[iCalo]);
            h_PFCluster_EoverEtrue_sim_fraction_EB->Fill(PFCluster_RecoEnergy_sim_fraction[iCalo]/simEnergy[iCalo]);
            h_PFCluster_EoverEtrue_sim_fraction_min1_EB->Fill(PFCluster_RecoEnergy_sim_fraction_min1[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_dR_simScore_EB->Fill(SuperCluster_RecoEnergy_dR_simScore[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_n_shared_xtals_EB->Fill(SuperCluster_RecoEnergy_n_shared_xtals[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_sim_fraction_EB->Fill(SuperCluster_RecoEnergy_sim_fraction[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_sim_fraction_min1_EB->Fill(SuperCluster_RecoEnergy_sim_fraction_min1[iCalo]/simEnergy[iCalo]);
         }else if(fabs(caloParticle_simEta_tmp.at(iCalo))>1.5){
            h_PFCluster_EoverEtrue_dR_simScore_EE->Fill(PFCluster_RecoEnergy_dR_simScore[iCalo]/simEnergy[iCalo]);
            h_PFCluster_EoverEtrue_n_shared_xtals_EE->Fill(PFCluster_RecoEnergy_n_shared_xtals[iCalo]/simEnergy[iCalo]);
            h_PFCluster_EoverEtrue_sim_fraction_EE->Fill(PFCluster_RecoEnergy_sim_fraction[iCalo]/simEnergy[iCalo]);
            h_PFCluster_EoverEtrue_sim_fraction_min1_EE->Fill(PFCluster_RecoEnergy_sim_fraction_min1[iCalo]/simEnergy[iCalo]);  
            h_SuperCluster_EoverEtrue_dR_simScore_EE->Fill(SuperCluster_RecoEnergy_dR_simScore[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_n_shared_xtals_EE->Fill(SuperCluster_RecoEnergy_n_shared_xtals[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_sim_fraction_EE->Fill(SuperCluster_RecoEnergy_sim_fraction[iCalo]/simEnergy[iCalo]);
            h_SuperCluster_EoverEtrue_sim_fraction_min1_EE->Fill(SuperCluster_RecoEnergy_sim_fraction_min1[iCalo]/simEnergy[iCalo]); 
         }

         if(PFCluster_nSimMatches_dR_simScore[iCalo]>0) h_PFCluster_dR_simScore_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(PFCluster_nSimMatches_n_shared_xtals[iCalo]>0) h_PFCluster_n_shared_xtals_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(PFCluster_nSimMatches_sim_fraction[iCalo]>0) h_PFCluster_sim_fraction_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(PFCluster_nSimMatches_sim_fraction_min1[iCalo]>0) h_PFCluster_sim_fraction_min1_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(SuperCluster_nSimMatches_dR_simScore[iCalo]>0) h_SuperCluster_dR_simScore_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(SuperCluster_nSimMatches_n_shared_xtals[iCalo]>0) h_SuperCluster_n_shared_xtals_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(SuperCluster_nSimMatches_sim_fraction[iCalo]>0) h_SuperCluster_sim_fraction_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));
         if(SuperCluster_nSimMatches_sim_fraction_min1[iCalo]>0) h_SuperCluster_sim_fraction_min1_Num_vs_Eta->Fill(caloParticle_simEta_tmp.at(iCalo));     
      }

      entry++;
      simEnergy.clear();
      PFCluster_RecoEnergy_dR_simScore.clear();
      PFCluster_nSimMatches_dR_simScore.clear();  
      PFCluster_RecoEnergy_n_shared_xtals.clear();
      PFCluster_nSimMatches_n_shared_xtals.clear(); 
      PFCluster_RecoEnergy_sim_fraction.clear();  
      PFCluster_nSimMatches_sim_fraction.clear(); 
      PFCluster_RecoEnergy_sim_fraction_min1.clear();
      PFCluster_nSimMatches_sim_fraction_min1.clear(); 
      SuperCluster_RecoEnergy_dR_simScore.clear();
      SuperCluster_nSimMatches_dR_simScore.clear();  
      SuperCluster_RecoEnergy_n_shared_xtals.clear();
      SuperCluster_nSimMatches_n_shared_xtals.clear();  
      SuperCluster_RecoEnergy_sim_fraction.clear();
      SuperCluster_nSimMatches_sim_fraction.clear();   
      SuperCluster_RecoEnergy_sim_fraction_min1.clear();
      SuperCluster_nSimMatches_sim_fraction_min1.clear();  
   }

   drawHisto(h_PFCluster_EoverEtrue_dR_simScore_EB, h_PFCluster_EoverEtrue_n_shared_xtals_EB, h_PFCluster_EoverEtrue_sim_fraction_EB, h_PFCluster_EoverEtrue_sim_fraction_min1_EB, std::string("h_PFCluster_EoverEtrue_EB"));
   drawHisto(h_PFCluster_EoverEtrue_dR_simScore_EE, h_PFCluster_EoverEtrue_n_shared_xtals_EE, h_PFCluster_EoverEtrue_sim_fraction_EE, h_PFCluster_EoverEtrue_sim_fraction_min1_EE, std::string("h_PFCluster_EoverEtrue_EE"));
   drawHisto(h_SuperCluster_EoverEtrue_dR_simScore_EB, h_SuperCluster_EoverEtrue_n_shared_xtals_EB, h_SuperCluster_EoverEtrue_sim_fraction_EB, h_SuperCluster_EoverEtrue_sim_fraction_min1_EB, std::string("h_SuperCluster_EoverEtrue_EB"));
   drawHisto(h_SuperCluster_EoverEtrue_dR_simScore_EE, h_SuperCluster_EoverEtrue_n_shared_xtals_EE, h_SuperCluster_EoverEtrue_sim_fraction_EE, h_SuperCluster_EoverEtrue_sim_fraction_min1_EE, std::string("h_SuperCluster_EoverEtrue_EE"));

   TEfficiency* eff_PFCluster_dR_simScore_vs_Eta = new TEfficiency(*h_PFCluster_dR_simScore_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   TEfficiency* eff_PFCluster_n_shared_xtals_vs_Eta = new TEfficiency(*h_PFCluster_n_shared_xtals_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   TEfficiency* eff_PFCluster_sim_fraction_vs_Eta = new TEfficiency(*h_PFCluster_sim_fraction_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   TEfficiency* eff_PFCluster_sim_fraction_min1_vs_Eta = new TEfficiency(*h_PFCluster_sim_fraction_min1_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   drawEfficiency(eff_PFCluster_dR_simScore_vs_Eta, eff_PFCluster_n_shared_xtals_vs_Eta, eff_PFCluster_sim_fraction_vs_Eta, eff_PFCluster_sim_fraction_min1_vs_Eta, std::string("eff_PFCluster_vs_Eta"));

   TEfficiency* eff_SuperCluster_dR_simScore_vs_Eta = new TEfficiency(*h_SuperCluster_dR_simScore_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   TEfficiency* eff_SuperCluster_n_shared_xtals_vs_Eta = new TEfficiency(*h_SuperCluster_n_shared_xtals_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   TEfficiency* eff_SuperCluster_sim_fraction_vs_Eta = new TEfficiency(*h_SuperCluster_sim_fraction_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   TEfficiency* eff_SuperCluster_sim_fraction_min1_vs_Eta = new TEfficiency(*h_SuperCluster_sim_fraction_min1_Num_vs_Eta,*h_CaloParticle_Denum_vs_Eta);
   drawEfficiency(eff_SuperCluster_dR_simScore_vs_Eta, eff_SuperCluster_n_shared_xtals_vs_Eta, eff_SuperCluster_sim_fraction_vs_Eta, eff_SuperCluster_sim_fraction_min1_vs_Eta, std::string("eff_SuperCluster_vs_Eta")); 
   
}

void drawHisto(TH1F* h_dR_simScore, TH1F* h_n_shared_xtals, TH1F* h_sim_fraction, TH1F* h_sim_fraction_min1, std::string Name)
{

   //h_dR_simScore->Scale(1./h_dR_simScore->GetEntries());
   for(int iBin=1; iBin<=h_dR_simScore->GetNbinsX();iBin++)
       if(h_dR_simScore->GetBinContent(iBin)==0.) h_dR_simScore->SetBinContent(iBin,0.8); 
   h_dR_simScore->SetLineColor(kBlue+1);
   h_dR_simScore->SetLineWidth(2);
   h_dR_simScore->GetXaxis()->SetTitle("E_{Reco}/E_{SIM}");

   //h_n_shared_xtals->Scale(1./h_n_shared_xtals->GetEntries());
   for(int iBin=1; iBin<=h_n_shared_xtals->GetNbinsX();iBin++)
       if(h_n_shared_xtals->GetBinContent(iBin)==0.) h_n_shared_xtals->SetBinContent(iBin,0.8); 
   h_n_shared_xtals->SetLineColor(kRed+1);
   h_n_shared_xtals->SetLineWidth(2);

   //h_sim_fraction->Scale(1./h_sim_fraction->GetEntries());
   for(int iBin=1; iBin<=h_sim_fraction->GetNbinsX();iBin++)
       if(h_sim_fraction->GetBinContent(iBin)==0.) h_sim_fraction->SetBinContent(iBin,0.8);
   h_sim_fraction->SetLineColor(kGreen+1);
   h_sim_fraction->SetLineWidth(2);

   //h_sim_fraction_min1->Scale(1./h_sim_fraction_min1->GetEntries());
   for(int iBin=1; iBin<=h_sim_fraction_min1->GetNbinsX();iBin++)
       if(h_sim_fraction_min1->GetBinContent(iBin)==0.) h_sim_fraction_min1->SetBinContent(iBin,0.8); 
   h_sim_fraction_min1->SetLineColor(kViolet);
   h_sim_fraction_min1->SetLineWidth(2);

   TLegend* legend = new TLegend(0.72, 0.72, 0.99, 0.94);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(h_dR_simScore,"deltaR","L");
   legend -> AddEntry(h_n_shared_xtals,"nxtals","L");
   legend -> AddEntry(h_sim_fraction,"sim_fraction","L");
   legend -> AddEntry(h_sim_fraction_min1,"sim_fraction_min1","L");

   TCanvas* c = new TCanvas();
   h_dR_simScore->Draw("hist");
   h_n_shared_xtals->Draw("hist same");
   h_sim_fraction->Draw("hist same");
   h_sim_fraction_min1->Draw("hist same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	

   c->SetLogy();
   h_dR_simScore->Draw("hist");
   h_n_shared_xtals->Draw("hist same");
   h_sim_fraction->Draw("hist same");
   h_sim_fraction_min1->Draw("hist same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+"_logY.png").c_str(),"png");
   c->SaveAs(std::string(Name+"_logY.pdf").c_str(),"pdf");	
}

void drawEfficiency(TEfficiency* eff_dR_simScore, TEfficiency* eff_n_shared_xtals, TEfficiency* eff_sim_fraction, TEfficiency* eff_sim_fraction_min1, std::string Name)
{
   eff_dR_simScore->SetLineColor(kBlue+1);
   eff_dR_simScore->SetLineWidth(2);
   eff_dR_simScore->SetTitle(std::string(Name+"; #eta ; Efficiency").c_str()); 
   
   eff_n_shared_xtals->SetLineColor(kRed+1);
   eff_n_shared_xtals->SetLineWidth(2);

   eff_sim_fraction->SetLineColor(kGreen+1);
   eff_sim_fraction->SetLineWidth(2);

   eff_sim_fraction_min1->SetLineColor(kViolet+1);
   eff_sim_fraction_min1->SetLineWidth(2);

   TLegend* legend = new TLegend(0.365, 0.12, 0.635, 0.34);
   legend -> SetFillColor(kWhite);
   legend -> SetFillStyle(1000);
   legend -> SetLineWidth(0);
   legend -> SetLineColor(kWhite);
   legend -> SetTextFont(42);  
   legend -> SetTextSize(0.04);
   legend -> AddEntry(eff_dR_simScore,"deltaR","L");
   legend -> AddEntry(eff_n_shared_xtals,"nxtals","L");
   legend -> AddEntry(eff_sim_fraction,"sim_fraction","L");
   legend -> AddEntry(eff_sim_fraction_min1,"sim_fraction_min1","L");

   TCanvas* c = new TCanvas();
   eff_dR_simScore->Draw("APL");
   eff_n_shared_xtals->Draw("PL, same");
   eff_sim_fraction->Draw("PL, same");
   eff_sim_fraction_min1->Draw("PL, same");
   legend -> Draw("same");
   c->SaveAs(std::string(Name+".png").c_str(),"png");
   c->SaveAs(std::string(Name+".pdf").c_str(),"pdf");	
}
