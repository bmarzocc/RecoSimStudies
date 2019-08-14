/*
 *  MustacheMacro.c
 *  
 *
 *  Created by Rishi Patel on 5/31/13.
 *  Copyright 2013 Rutgers University. All rights reserved.
 *
 */
 #include<TFile.h>
 #include<TTree.h>
 #include<TGraph.h>
 #include<TH2F.h>
 #include<TCanvas.h>
 #include<TVector2.h>
 #include<TMath.h>
 #include <algorithm> 
 #include<iostream>
#include<TStyle.h>
//Tree Variables

Int_t N_ECALClusters;
Int_t nVtx;
Float_t genEnergy;
Float_t genDRToCentroid;
Float_t genDRToSeed;
Float_t scSeedRawEnergy;
Float_t scSeedEta;
Float_t scSeedPhi;
Float_t scRawEnergy;
Float_t scEta;
//clusterVar
Float_t clusterRawEnergy[20];
Float_t clusterEta[20];
Float_t clusterPhi[20];
Float_t clusterDEtaToSeed[20];
Float_t clusterDPhiToSeed[20];
Int_t clusterInMustache[20];
Int_t clusterLeakage[20];
TTree* EventTree;
Int_t EvMax;
TH2F* axis=new TH2F("axis", "", 100, -0.6, 0.6, 100,-0.15, 0.15);
TCanvas* c1;
TCanvas* c2;
TH2F* EtaLeakage=new TH2F("EtaLeakage", "", 200,0,2.5,200, 0.0,0.05);
TH2F* EtaLeakageMust=new TH2F("EtaLeakageMust", "", 200,0,2.5,200, 0.0,0.05);
void InitTree(TString Filename){
	TFile* f1=new TFile(Filename);
	f1->cd("bigBoxSCTree");
	EventTree=(TTree*)gDirectory->Get("SuperClusterTree");
	EventTree->SetBranchAddress("nVtx", &nVtx);
	EventTree->SetBranchAddress("genEnergy", &genEnergy);
	EventTree->SetBranchAddress("scRawEnergy", &scRawEnergy);
	EventTree->SetBranchAddress("genDRToCentroid", &genDRToCentroid);
	EventTree->SetBranchAddress("genDRToSeed", &genDRToSeed);
	EventTree->SetBranchAddress("scSeedRawEnergy", &scSeedRawEnergy);
	EventTree->SetBranchAddress("scEta", &scEta);
	EventTree->SetBranchAddress("scSeedEta", &scSeedEta);
	EventTree->SetBranchAddress("scSeedPhi", &scSeedPhi);
	
	EventTree->SetBranchAddress("N_ECALClusters", &N_ECALClusters);
	EventTree->SetBranchAddress("clusterRawEnergy", clusterRawEnergy);
	EventTree->SetBranchAddress("clusterEta", clusterEta);
	EventTree->SetBranchAddress("clusterPhi", clusterPhi);
	EventTree->SetBranchAddress("clusterDEtaToSeed", clusterDEtaToSeed);
	EventTree->SetBranchAddress("clusterDPhiToSeed", clusterDPhiToSeed);
	EventTree->SetBranchAddress("clusterInMustache", clusterInMustache);
	EventTree->SetBranchAddress("clusterLeakage", clusterLeakage);
	EvMax=EventTree->GetEntries();
}
bool inMustache(const float maxEta, const float maxPhi, 
				const float ClustE, const float ClusEta, 
				const float ClusPhi){
	//bool inMust=false;
	//float eta0 = maxEta;
	//float phi0 = maxPhi;      
	
	float p00 = -0.107537;
	float p01 = 0.590969;
	float p02 = -0.076494;
	float p10 = -0.0268843;
	float p11 = 0.147742;
	float p12 = -0.0191235;
	
	float w00 = -0.00571429;
	float w01 = -0.002;
	float w10 = 0.0135714;
	float w11 = 0.001;
	
	const float sineta0 = sin(maxEta);
	const float eta0xsineta0 = maxEta*sineta0;
	
	
	//2 parabolas (upper and lower) 
	//of the form: y = a*x*x + b
	
	//b comes from a fit to the width
	//and has a slight dependence on E on the upper edge    
	// this only works because of fine tuning :-D
	const float sqrt_log10_clustE = sqrt(log10(ClustE)+1.1);
	// we need to have this in two steps, so that we don't improperly shift
	// the lower bound!
	float b_upper = w10*eta0xsineta0 + w11 / sqrt_log10_clustE;      
	float b_lower = w00*eta0xsineta0 + w01 / sqrt_log10_clustE; 
	const float midpoint =  0.5*( b_upper + b_lower );
	b_upper -= midpoint;
	b_lower -= midpoint;      
	
	//the curvature comes from a parabolic 
	//fit for many slices in eta given a 
	//slice -0.1 < log10(Et) < 0.1
	const float curv_up=eta0xsineta0*(p00*eta0xsineta0+p01)+p02;
	const float curv_low=eta0xsineta0*(p10*eta0xsineta0+p11)+p12;
	
	//solving for the curviness given the width of this particular point
	const float a_upper=(1/(4*curv_up))-fabs(b_upper);
	const float a_lower = (1/(4*curv_low))-fabs(b_lower);
	
	const double dphi=TVector2::Phi_mpi_pi(ClusPhi-maxPhi);
	const double dphi2 = dphi*dphi;
	// minimum offset is half a crystal width in either direction
	// because science.
	const float upper_cut=(1./(4.*a_upper))*dphi2+b_upper; 
	const float lower_cut=(1./(4.*a_lower))*dphi2+b_lower;
	
	//if(deta < upper_cut && deta > lower_cut) inMust=true;
	
	const float deta=(1-2*(maxEta<0))*(ClusEta-maxEta); // sign flip deta
	return (deta < upper_cut && deta > lower_cut);
}

void MustachePlot(float minEt, float maxEt){
	c1 = new TCanvas("c1","Mustache Et Bins",0,0,1000,1000);
	c1->Divide(4,3);
	for(int i=0; i<12; ++i){
		float maxEta=(i+1)*0.2;
		float minEta=(i)*0.2;
		int incl=1;
		int excl=1;
		int l=1;
		TGraph* grMust=new TGraph();
		TGraph* grExc=new TGraph();
		TGraph* grL=new TGraph();
		for(int iev=0; iev<EvMax;++iev){
		EventTree->GetEvent(iev);
		if(N_ECALClusters==0)continue;
		if((scSeedEta)>maxEta || (scSeedEta)<minEta)continue;
		for(int ic=0; ic<N_ECALClusters; ++ic){
			if(clusterRawEnergy[ic]/cosh(clusterEta[ic])>maxEt || clusterRawEnergy[ic]/cosh(clusterEta[ic])<minEt )continue;
			bool inMust=inMustache(scSeedEta, scSeedPhi, clusterRawEnergy[ic], clusterEta[ic], clusterPhi[ic]);
			if(inMust){
				grMust->Expand(incl);
				grMust->SetPoint(incl,clusterDPhiToSeed[ic], clusterDEtaToSeed[ic]);
				//cout<<"incl"<<clusterRawEnergy[ic]<<endl;
				++incl;
			}
			else{
				grExc->Expand(excl);
				grExc->SetPoint(excl, clusterDPhiToSeed[ic], clusterDEtaToSeed[ic]);
				++excl;
			}
			if(clusterLeakage[ic]!=0 && !inMust ){
				grL->Expand(l);
				grL->SetPoint(l, clusterDPhiToSeed[ic], clusterDEtaToSeed[ic]);
				++l;
			}
		  }
		}
		
		c1->cd(i+1);
		axis->Draw("");
		grMust->GetXaxis()->SetLimits(-0.6, 0.6);
		grMust->GetYaxis()->SetLimits(-0.15, 0.15);
		grMust->SetMarkerColor(kBlue);
		grMust->SetMarkerSize(0.5);
		grMust->Draw("P");	
		grExc->GetXaxis()->SetLimits(-0.6, 0.6);
		grExc->GetYaxis()->SetLimits(-0.15, 0.15);
		grExc->SetMarkerColor(kRed);
		grExc->SetMarkerSize(0.5);
		grExc->Draw("PSame");
		grL->GetXaxis()->SetLimits(-0.6, 0.6);
		grL->GetYaxis()->SetLimits(-0.15, 0.15);
		grL->SetMarkerSize(0.5);
		grL->SetMarkerColor(kGreen);
		grL->SetMarkerStyle(24);
		grL->Draw("PSame");
		
	}
	c1->Print(TString::Format("TestMustacheShapesEt%0.4fto%0.4f.gif", minEt, maxEt));
	
}

bool pairgreater(const std::pair<float,int>& one, const std::pair<float,int>& two) { return one.first > two.first; }

void ContainmentAnalyzer(){
	
	
	int evcount=1;
	std::vector< std::pair<float,int> > Pair;
	for(int iev=0; iev<EvMax;++iev){
		EventTree->GetEvent(iev);
		if(N_ECALClusters==0)continue;
		//std::vector<int>SortInd;
		for(int ic=0; ic<N_ECALClusters; ++ic){
			Pair.push_back(make_pair(fabs(clusterDEtaToSeed[ic]),ic));
			//SortInd[ic]=ic;
		}
		
		std::sort(Pair.begin(),Pair.end()+1,pairgreater);
		float ClusSumE=scSeedRawEnergy;
		float ClustDEta=0.;
		float ClustDEtaMust=0;
		for(int ic=0; ic<N_ECALClusters; ++ic){
			int index=Pair[ic].second;
			ClusSumE=ClusSumE+clusterRawEnergy[index];
			if(ClusSumE/scRawEnergy>0.995){
				ClustDEta=clusterDEtaToSeed[index];
				if(clusterInMustache[index]==1)ClustDEtaMust=clusterDEtaToSeed[index];
			}
		}
		//std::cout<<scSeedEta<<", "<<ClustDEta<<std::endl;
		//grLeakage->SetPoint(evcount,scSeedEta,ClustDEta);
		Double_t seedEta=scSeedEta;
		if(ClustDEta>0)EtaLeakage->Fill(scSeedEta,ClustDEta);
		if(ClustDEtaMust>0)EtaLeakageMust->Fill(scSeedEta,ClustDEtaMust);
		++evcount;
		Pair.clear();
	}
	c2=new TCanvas("c2","Mustache Containment",0,0,1000,1000);
	c2->Divide(3,1);
	c2->cd(1);
	//grLeakage->Draw("P");
	gStyle->SetPalette(1);
	EtaLeakage->Draw("colz");
	c2->cd(2);
	EtaLeakageMust->Draw("colz");
	Double_t DEta[200];
	Double_t DEtaMust[200];
	Double_t EtaSeed[200];
	Double_t Displacement[200];
	for(int i=1; i<200; ++i){
		EtaSeed[i-1]=EtaLeakage->ProjectionX()->GetBinLowEdge(i);
		double IntegralTot=EtaLeakage->Integral(i,i,1,200);
		for(int j=1; j<200; ++j){
			double IntegralDeta=EtaLeakage->Integral(i,i,1, j+100);
			//cout<<"IntegralDeta "<<IntegralDeta << "Total "<<IntegralTot<<endl;
			if(IntegralDeta/IntegralTot>0.995){
				//cout<<"Fraction "<<IntegralDeta/IntegralTot<<endl;
				DEta[i-1]=EtaLeakage->ProjectionY()->GetBinLowEdge(j+100);
				//return;
				break;
			}
		}
		for(int j=1; j<200; ++j){
			double IntegralDeta=EtaLeakageMust->Integral(i,i,1, j+100);
			//cout<<"IntegralDeta "<<IntegralDeta << "Total "<<IntegralTot<<endl;
			if(IntegralDeta/IntegralTot>0.995){
				//cout<<"Fraction "<<IntegralDeta/IntegralTot<<endl;
				DEtaMust[i-1]=EtaLeakageMust->ProjectionY()->GetBinLowEdge(j+100);
				//return;
				break;
			}
		}
		Displacement[i-1]=DEta[i-1]-DEtaMust[i-1];
	}
	c2->cd(3);
	TGraph* grLeakage=new TGraph(200,EtaSeed,Displacement);
	grLeakage->SetMarkerStyle(kOpenCircle);
	grLeakage->SetMarkerSize(1.0);
	grLeakage->Draw("AP");
	//Now Loop over to 
}

void MustacheMacro(){
	InitTree("MustacheCompTreeHighStats.root");
	//doSomething
	ContainmentAnalyzer();
	return;
	//1 Mustache Plot
	for(int e=0; e<4; ++e){
		float minEt=e*0.5;
		float maxEt=(e+1)*0.5;
		MustachePlot(minEt, maxEt);
	}
	
	//2 Dynamic Dphi
	
	//3 routine to find maximal containment (dEta adjustment)
	
	
}
