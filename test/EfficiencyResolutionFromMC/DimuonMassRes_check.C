#include "TFile.h"
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TPad.h"

void DimuonMassRes_check(){

gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetLegendTextSize(0.03);

gROOT->LoadMacro("cruijff.C+");

gROOT->Reset();
gROOT->SetBatch();

// 	    Double_t MASS_BINS[] = {120, 200, 300, 400, 600, 800, 1000, 1300, 1600, 2000, 2500, 3100, 3800, 4500, 5500};
		Double_t MASS_BINS[] = {50, 120, 200, 400, 800, 1400, 2300, 3500, 4500, 6000};
	    Int_t  binnum = sizeof(MASS_BINS)/sizeof(Double_t)-1;
	    
	    float BB_sigma[14] = {0};
	    float BB_sigma_err[14] = {0};
	    float BE_sigma[14] = {0};
	    float BE_sigma_err[14] = {0};
	    float EE_sigma[14] = {0};
	    float EE_sigma_err[14] = {0};
	    
          std::vector<float> SIGMA;
          std::vector<float> SIGMA_ERR;
		
     TString samples[39] =  {"dyInclusive50", 
	 						"qcd80to120", "qcd120to170", "qcd170to300", "qcd300to470", "qcd470to600", "qcd600to800", "qcd800to1000", "qcd1000to1400", "qcd1400to1800", "qcd1800to2400", "qcd2400to3200", "qcd3200", 
	 						"Wantitop", "tW", 
	 						"ttbar_lep50to500", "ttbar_lep_500to800", "ttbar_lep_800to1200", "ttbar_lep_1200to1800", "ttbar_lep1800toInf", 
	 						"Wjets", 
	 						"WWinclusive", "WW200to600", "WW600to1200", "WW1200to2500", "WW2500", 
	 						"ZZ", "WZ", "ZZ_ext", "WZ_ext", 
	 						"dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000"};


	float events[39] = {19385554,  6986740, 6708572, 6958708, 4150588, 3959986, 3896412, 3992112, 2999069, 396409, 397660, 399226, 391735, 
						6933094, 6952830,
						79092400, 200000, 199800, 200000, 40829, 
						29705748,
						1999000, 200000, 200000, 200000, 38969, 
						990064, 1000000, 998034, 2995828,
						2760000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000
						};
						
	float sigma[39] = {6025.2, 2762530, 471100, 117276, 7823, 648.2, 186.9, 32.3992112, 9.4183, 0.84265, 0.114943, 0.00682981, 0.000165445,
						35.6, 35.6,
						87.31, 0.32611, 0.03265, 0.00305, 0.00017, 
						61526.7,
						12.178, 1.385, 0.0566, 0.0035, 0.00005,
						16.523, 47.13, 16.523, 47.13,
						1975, 19.32,  2.731, 0.241, 0.01678, 0.00139, 0.00008948, 0.0000041, 4.56E-7
						};
						
	float LUMINOSITY = 41903;


	float weight[39] = {0};
	
	TString NOME;
	float  fit_min, fit_max;
	float yPos;
	TString longstring;
	

	
	
	Long64_t ne;
	Long64_t Nne;
  float dil_mass;
  float cos_angle;
  float vertex_chi2;
  int dil_chosen;
  
  Int_t event;
  Int_t run;
  unsigned lumi;

  float lep_pt[2];
  float lep_phi[2];
  int lep_id[2];
  float lep_pt_err[2];
  float lep_eta[2];
  float lep_tk_pt[2];
  float lep_tk_eta[2];
  float lep_tk_phi[2];
  float lep_std_pt[2];
  float lep_std_eta[2];
  float lep_std_phi[2];
  float lep_picky_pt[2];
  float lep_picky_eta[2];
  float lep_picky_phi[2];
  float lep_dyt_pt[2];
  float lep_dyt_eta[2];
  float lep_dyt_phi[2];
  float lep_glb_pt[2];
  float lep_tpfms_pt[2];
  float lep_dB[2];
  float lep_sumPt[2];
  float lep_triggerMatchPt[2];
  short lep_glb_numberOfValidTrackerLayers[2]; 
  short lep_glb_numberOfValidPixelHits[2];
  short lep_glb_numberOfValidMuonHits[2];
  short lep_numberOfMatchedStations[2];
  bool lep_isGlobalMuon[2];
  bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool GoodVtx;
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    float vertex_m;
    float gen_dil_mass;
    float gen_lep_qOverPt[2];
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
	
	float MASS = -1;
	float MASS_SCALE = -1;
	float MASS_GEN = -1;
	int i = 0;
	int n_count;
	int p_count;
	
	float gen_eta_first = -999;
	float gen_eta_second = -999;
	
	float gen_pt_first = -999;
	float gen_pt_second = -999;
	
	float pt_first = -999;
	float pt_second = -999;
	float max_pt = -999;
	
	float gen_pt_plus = -999;
	float gen_pt_minus = -999;
		
	int pp;
	bool next;
	
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 
	
  Int_t prev_event = -88;
  
	float count_event[9] = {0};
	float count_event_BB[9] = {0};
	float count_event_BE[9] = {0};

			
	for(int i = 0; i < 39; i++){
// 		weight[i] = LUMINOSITY * sigma[i] / events[i];
		weight[i] = 1;
	}
	
	double mass;
	double rdil;
	TH2F* DileptonMassResVMass_2d = new TH2F("(dil. mass - gen dil. mass)/(gen dil. mass)", "(dil. mass - gen dil. mass)/(gen dil. mass)", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F* DileptonMassResVMass_2d_BB = new TH2F("(dil. mass - gen dil. mass)/(gen dil. mass): BB", "(dil. mass - gen dil. mass)/(gen dil. mass): BB", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F* DileptonMassResVMass_2d_BE = new TH2F("(dil. mass - gen dil. mass)/(gen dil. mass): BE", "(dil. mass - gen dil. mass)/(gen dil. mass): BE", binnum, MASS_BINS, 400, -1., 1.0);
	TH2F* DileptonMassResVMass_2d_EE = new TH2F("(dil. mass - gen dil. mass)/(gen dil. mass: EE)", "(dil. mass - gen dil. mass)/(gen dil. mass): EE", binnum, MASS_BINS, 400, -1., 1.0);

	TH1F *res = new TH1F("Resolution", "Resolution", binnum, MASS_BINS);
	res->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res->GetYaxis()->SetTitle("Entries"); 
	res->SetTitle("Mass residuals");
	TH1F *res_bb = new TH1F("Resolution BB", "Resolution BB", binnum, MASS_BINS);
	res_bb->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res_bb->GetYaxis()->SetTitle("Entries"); 
	res_bb->SetTitle("BB mass residuals");
	TH1F *res_be = new TH1F("Resolution BE + EE", "Resolution BE + EE", binnum, MASS_BINS);
	res_be->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res_be->GetYaxis()->SetTitle("Entries"); 
	res_be->SetTitle("BE + EE mass residuals");
	TH1F *res_ee = new TH1F("Resolution EE", "Resolution EE", binnum, MASS_BINS);
	res_ee->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV]");
// 	res_ee->GetYaxis()->SetTitle("Entries"); 
	res_ee->SetTitle("EE mass residuals");
		                    	
  for(int j=30; j < 39; j++){   

	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    

     TChain *treeMC = new TChain("SimpleNtupler/t");

     treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/ana_datamc_"+samples[j]+".root");

     treeMC->SetBranchAddress("event",&event);
     treeMC->SetBranchAddress("run",&run);
     treeMC->SetBranchAddress("lumi",&lumi);
     treeMC->SetBranchAddress("dil_mass",&dil_mass);
     treeMC->SetBranchAddress("cos_angle",&cos_angle);
     treeMC->SetBranchAddress("vertex_chi2",&vertex_chi2);
     treeMC->SetBranchAddress("dil_chosen",&dil_chosen);
     treeMC->SetBranchAddress("lep_pt",lep_pt);
     treeMC->SetBranchAddress("lep_id",lep_id);
     treeMC->SetBranchAddress("lep_eta",lep_eta);
     treeMC->SetBranchAddress("lep_phi",lep_phi);
     treeMC->SetBranchAddress("lep_dB",lep_dB);
     treeMC->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
     treeMC->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
     treeMC->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
     treeMC->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeMC->SetBranchAddress("lep_tk_eta",lep_tk_eta);
     treeMC->SetBranchAddress("lep_tk_phi",lep_tk_phi);
     treeMC->SetBranchAddress("lep_std_pt",lep_std_pt);
     treeMC->SetBranchAddress("lep_std_eta",lep_std_eta);
     treeMC->SetBranchAddress("lep_std_phi",lep_std_phi);
     treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
     treeMC->SetBranchAddress("lep_picky_eta",lep_picky_eta);
     treeMC->SetBranchAddress("lep_picky_phi",lep_picky_phi);
     treeMC->SetBranchAddress("lep_dyt_pt",lep_dyt_pt);
     treeMC->SetBranchAddress("lep_dyt_eta",lep_dyt_eta);
     treeMC->SetBranchAddress("lep_dyt_phi",lep_dyt_phi);
     treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeMC->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
     treeMC->SetBranchAddress("lep_pt_err",lep_pt_err);
     treeMC->SetBranchAddress("vertex_m",&vertex_m);
     treeMC->SetBranchAddress("GoodDataRan", &GoodDataRan);
     treeMC->SetBranchAddress("GoodVtx",&GoodVtx);
     treeMC->SetBranchAddress("lep_stationMask",&lep_stationMask);
     treeMC->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
     treeMC->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
     treeMC->SetBranchAddress("gen_lep_qOverPt", gen_lep_qOverPt);
     treeMC->SetBranchAddress("gen_lep_eta", gen_lep_eta);
     treeMC->SetBranchAddress("gen_lep_pt", gen_lep_pt);	

	
	ne = treeMC->GetEntries();
	std::cout<<"START"<<std::endl;
	for ( int p=0; p < ne ;p++){
// 	for ( int p=0; p<1000 ;p++){
	
		if(p % 100000 == 0) std::cout<<p<<" su "<<ne<<std::endl;		
			
		treeMC->GetEntry(p);
				
		if(
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			( lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
			lep_glb_numberOfValidMuonHits[0]>0 && lep_glb_numberOfValidMuonHits[1]>0 && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){
			
				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
				prev_event = event;

			count_event[j-30]++;

			
// 			lep_1.SetPtEtaPhiM(lep_tk_pt[0], lep_tk_eta[0], lep_tk_phi[0], 0.105);
// 			lep_2.SetPtEtaPhiM(lep_tk_pt[1], lep_tk_eta[1], lep_tk_phi[1], 0.105);
			lep_1.SetPtEtaPhiM(lep_std_pt[0], lep_std_eta[0], lep_std_phi[0], 0.105);
			lep_2.SetPtEtaPhiM(lep_std_pt[1], lep_std_eta[1], lep_std_phi[1], 0.105);
// 			lep_1.SetPtEtaPhiM(lep_picky_pt[0], lep_picky_eta[0], lep_picky_phi[0], 0.105);
// 			lep_2.SetPtEtaPhiM(lep_picky_pt[1], lep_picky_eta[1], lep_picky_phi[1], 0.105);
// 			lep_1.SetPtEtaPhiM(lep_dyt_pt[0], lep_dyt_eta[0], lep_dyt_phi[0], 0.105);
// 			lep_2.SetPtEtaPhiM(lep_dyt_pt[1], lep_dyt_eta[1], lep_dyt_phi[1], 0.105);

			ZPrime = lep_1 + lep_2;

			mass 		 = ZPrime.M(); // resolution using Tracker track		
// 			mass         = dil_mass; //vertex_m

			rdil    = mass    /gen_dil_mass     - 1;
			DileptonMassResVMass_2d   ->Fill(gen_dil_mass,     rdil, weight[j]);

		  if (fabs(gen_lep_eta[0]) < 0.9 && fabs(gen_lep_eta[1]) < 0.9){ 
		  	DileptonMassResVMass_2d_BB ->Fill(gen_dil_mass,     rdil, weight[j]);
		  }
		  else if (fabs(gen_lep_eta[0]) > 1.2 && fabs(gen_lep_eta[1]) > 1.2){
		  	DileptonMassResVMass_2d_EE ->Fill(gen_dil_mass,     rdil, weight[j]);
		  }
		  else{
		  	DileptonMassResVMass_2d_BE ->Fill(gen_dil_mass,     rdil, weight[j]);
		  }

			/*
			if (gen_lep_eta[0]<=1.2 && gen_lep_eta[1]<=1.2 && gen_lep_eta[0]>=-1.2 && gen_lep_eta[1]>=-1.2){ 
				DileptonMassResVMass_2d_BB ->Fill(gen_dil_mass,     rdil, weight[j]);
				count_event_BB[j-30]++;
			}

			else{                                                                    
				DileptonMassResVMass_2d_BE ->Fill(gen_dil_mass,     rdil, weight[j]);
				count_event_BE[j-30]++;
			}
			if((gen_lep_eta[0] > 1.2 || gen_lep_eta[0] < -1.2) && (gen_lep_eta[1] > 1.2 || gen_lep_eta[1] < -1.2))                                                                         
				DileptonMassResVMass_2d_EE ->Fill(gen_dil_mass,     rdil, weight[j]);	
			*/

			
		} //if selection
	} // for entries
 } // for samples

//  TCanvas* c1 = new TCanvas("res", "res", 600, 600);
//  DileptonMassResVMass_2d->Draw();
//  TCanvas* c2 = new TCanvas("res_BB", "res_BB", 600, 600);
//  DileptonMassResVMass_2d_BB->Draw();
//  TCanvas* c3 = new TCanvas("res_BE", "res_BE", 600, 600);
//  DileptonMassResVMass_2d_BE->Draw();
//  TCanvas* c4 = new TCanvas("res_EE", "res_EE", 600, 600);
//  DileptonMassResVMass_2d_EE->Draw();


 bool save = true;

 
 for(int i = 1; i <= DileptonMassResVMass_2d_BB->GetNbinsX(); i++){
 	
 	NOME = Form("DimuonMass Resolution BB: %.0f < m < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMassResVMass_2d_BB->ProjectionY(NOME,i,i);
 	proiezione->Rebin(2);
 	proiezione->Draw();
 	
 	fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BB: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    proiezione->SetTitle(title);
    proiezione->GetXaxis()->SetTitle("Reco / Gen - 1");
        
    res_bb->SetBinContent(i, funct->GetParameter(2));
    res_bb->SetBinError(i, funct->GetParError(2));
    
    BB_sigma[i-1] = funct->GetParameter(2);
    BB_sigma_err[i-1] = funct->GetParError(2);
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar()+1; i++){
		if(i == funct->GetNpar()){
	    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    		longstring = Form("#chi^{2}/ndof = %5.3f/%d = %5.3f", funct->GetChisquare(), funct->GetNDF(), funct->GetChisquare()/(Double_t)funct->GetNDF());
        	latexFit->DrawLatex(-0.85, yPos, longstring);
        	continue;
		}
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(-0.85, yPos, longstring);
    }
 	 if(save){	
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/BB_resolution.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/BB_resolution.pdf");
 	if(i==DileptonMassResVMass_2d_BB->GetNbinsX())
 		canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/BB_resolution.pdf]");
    canvas->Write();
	}
}

 for(int i = 1; i <= DileptonMassResVMass_2d_BE->GetNbinsX(); i++){
 	
 	NOME = Form("DimuonMass Resolution BE: %.0f < m < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMassResVMass_2d_BE->ProjectionY(NOME,i,i);
 	proiezione->Rebin(2);
 	proiezione->Draw();
 	
 	fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("BE: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    proiezione->SetTitle(title); 
    proiezione->GetXaxis()->SetTitle("Reco / Gen - 1");
    
    res_be->SetBinContent(i, funct->GetParameter(2));
    res_be->SetBinError(i, funct->GetParError(2));
    
    BE_sigma[i-1] = funct->GetParameter(2);
    BE_sigma_err[i-1] = funct->GetParError(2);
    
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
// 	TString longstring_2;
    for(int i = 0; i < funct->GetNpar()+1; i++){
		if(i == funct->GetNpar()){
	    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    		longstring = Form("#chi^{2}/ndof = %5.3f/%d = %5.3f", funct->GetChisquare(), funct->GetNDF(), funct->GetChisquare()/(Double_t)funct->GetNDF());
        	latexFit->DrawLatex(-0.85, yPos, longstring);
        	continue;
		}
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(-0.85, yPos, longstring);
    }
 if(save){
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/BE_resolution.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/BE_resolution.pdf");
 	if(i==DileptonMassResVMass_2d_BE->GetNbinsX())
 		canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/BE_resolution.pdf]");
    canvas->Write();
	}
}

 for(int i = 1; i <= DileptonMassResVMass_2d_EE->GetNbinsX(); i++){
 	
 	NOME = Form("DimuonMass Resolution EE: %.0f < m < %.0f", MASS_BINS[i-1], MASS_BINS[i]);
 	TCanvas *canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	TH1D* proiezione = DileptonMassResVMass_2d_EE->ProjectionY(NOME,i,i);
 	proiezione->Rebin(2);
 	proiezione->Draw();
 	
 	fit_min = proiezione->GetMean() - 2.0*proiezione->GetRMS();
 	fit_max = proiezione->GetMean() + 1.7*proiezione->GetRMS();
 	TF1 *gaus = new TF1("gaus","gaus",fit_min,fit_max);
 	gaus->SetParameters(0, proiezione->GetMean(), proiezione->GetRMS());
    proiezione->Fit("gaus","M0R+");
    
    TF1* funct = new  TF1("cruijff", "cruijff", fit_min,fit_max, 5);
    funct->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2), 0., 0.);// #15, 0.001);             
    funct->SetParNames("Constant","Mean","Sigma","AlphaL","AlphaR");
    funct->SetLineColor(kBlue);
    funct->SetLineWidth(2);
    proiezione->Fit("cruijff","MR+"); 
    TString title = Form("EE: Mass resolution for %.0f < m_{ll} <%.0f", MASS_BINS[i-1], MASS_BINS[i]);
    proiezione->SetTitle(title);
    proiezione->GetXaxis()->SetTitle("Reco / Gen - 1");
    
    res_ee->SetBinContent(i, funct->GetParameter(2));
    res_ee->SetBinError(i, funct->GetParError(2));
 
    EE_sigma[i-1] = funct->GetParameter(2);
    EE_sigma_err[i-1] = funct->GetParError(2);
   
    TLatex* latexFit = new TLatex();
    latexFit->SetTextFont(42);
    latexFit->SetTextSize(0.030);
    for(int i = 0; i < funct->GetNpar()+1; i++){
		if(i == funct->GetNpar()){
	    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    		longstring = Form("#chi^{2}/ndof = %5.3f/%d = %5.3f", funct->GetChisquare(), funct->GetNDF(), funct->GetChisquare()/(Double_t)funct->GetNDF());
        	latexFit->DrawLatex(-0.85, yPos, longstring);
        	continue;
		}
    	yPos = proiezione->GetMaximum() - 5*i*proiezione->GetMaximum()/(float)100;
    	longstring = Form("%s = %5.3g #pm %5.3g", funct->GetParName(i),funct->GetParameter(i),funct->GetParError(i));
        latexFit->DrawLatex(-0.85, yPos, longstring);
    }
    
 if(save){	 	
 	if(i==1)
 		canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/EE_resolution.pdf[");
	canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/EE_resolution.pdf");
 	if(i==DileptonMassResVMass_2d_EE->GetNbinsX())
 		canvas->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/EE_resolution.pdf]");
    canvas->Write();
	}
}

 
 TCanvas* c1 = new TCanvas("res", "res", 500, 500);
 c1->SetGrid();
 res_bb->SetTitle("");
 res_bb->SetLineColor(kRed+1);
 res_bb->SetMarkerStyle(20);
 res_bb->SetMarkerSize(0.5);
 res_bb->SetMarkerColor(kRed+1);
 res_bb->GetYaxis()->SetRangeUser(0, 0.5);
 res_bb->Draw();
 res_be->SetTitle("");
 res_be->SetLineColor(kGreen+1);
 res_be->SetMarkerStyle(20);
 res_be->SetMarkerSize(0.5);
 res_be->SetMarkerColor(kGreen+1);
 res_be->GetYaxis()->SetRangeUser(0, 0.5);
 res_be->Draw("same");
 res_ee->SetTitle("");
 res_ee->SetLineColor(kBlue+1);
 res_ee->SetMarkerStyle(20);
 res_ee->SetMarkerSize(0.5);
 res_ee->SetMarkerColor(kBlue+1);
 res_ee->GetYaxis()->SetRangeUser(0, 0.5);
 res_ee->Draw("same");
 TLegend *legend_3 = new TLegend(0.2,0.75,0.7,0.85);
 legend_3->AddEntry(res_bb, "BB category: |#eta| < 0.9", "lep");
 legend_3->AddEntry(res_ee,"EE category: |#eta| > 1.2", "lep");
 legend_3->AddEntry(res_be,"BE category: elsewhere", "lep");
 legend_3->Draw(); 
 if(save){
 c1->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/Resolution.png");
 c1->Print("./Check_Resolution_Code_2017_Std/Split_3categorie/Resolution.pdf");
	}


	std::cout<<"MC:"<<std::endl;
	std::cout<<"Double_t BB_NUM[] = {";
	for(int i = 0; i < binnum; i++){
		std::cout<<BB_sigma[i];
		if(i != binnum-1) std::cout<<", ";
		else std::cout<<"};"<<std::endl;
	}
	std::cout<<"Double_t BB_err_NUM[] = {";
	for(int i = 0; i < binnum; i++){
		std::cout<<BB_sigma_err[i];
		if(i != binnum-1) std::cout<<", ";
		else std::cout<<"};"<<std::endl;
	}
	std::cout<<"Double_t BE_NUM[] = {";
	for(int i = 0; i < binnum; i++){
		std::cout<<BE_sigma[i];
		if(i != binnum-1) std::cout<<", ";
		else std::cout<<"};"<<std::endl;
	}
	std::cout<<"Double_t BE_err_NUM[] = {";
	for(int i = 0; i < binnum; i++){
		std::cout<<BE_sigma_err[i];
		if(i != binnum-1) std::cout<<", ";
		else std::cout<<"};"<<std::endl;
	}
	std::cout<<"Double_t EE_NUM[] = {";
	for(int i = 0; i < binnum; i++){
		std::cout<<EE_sigma[i];
		if(i != binnum-1) std::cout<<", ";
		else std::cout<<"};"<<std::endl;
	}
	std::cout<<"Double_t EE_err_NUM[] = {";
	for(int i = 0; i < binnum; i++){
		std::cout<<EE_sigma_err[i];
		if(i != binnum-1) std::cout<<", ";
		else std::cout<<"};"<<std::endl;
	}

//  std::cout<<"{"<<std::endl;
//  for(int i = 0; i < binnum; i++)
// 	std::cout<<BB_sigma[i]<<", ";
//  std::cout<<"};\n{"<<std::endl;
//  for(int i = 0; i < binnum; i++)
// 	std::cout<<BB_sigma_err[i]<<", ";	
//  std::cout<<"};\n{"<<std::endl;
//  for(int i = 0; i < binnum; i++)
// 	std::cout<<BE_sigma[i]<<", ";
//  std::cout<<"};\n{"<<std::endl;
//  for(int i = 0; i < binnum; i++)
// 	std::cout<<BE_sigma_err[i]<<", ";	    
//  std::cout<<"};"<<std::endl;

// 
//  for(int i = 0; i < 9; i++) 
// 	std::cout<<count_event[i]<<", ";
// std::cout<<""<<std::endl;	     
//  for(int i = 0; i < 9; i++) 
// 	std::cout<<count_event_BB[i]<<", ";
// std::cout<<""<<std::endl;
//  for(int i = 0; i < 9; i++) 
// 	std::cout<<count_event_BE[i]<<", ";
// std::cout<<""<<std::endl;
 
} // main function


// entries
// 0, 23908, 48176, 64112, 75469, 81654, 84220, 84185, 83458, 
// 0, 10303, 17285, 22187, 29176, 38058, 46160, 50871, 52300, 
// 0, 13605, 30891, 41925, 46293, 43596, 38060, 33314, 31158,

