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


void EventSelectionBase(){

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);
gStyle->SetLegendTextSize(0.03);
	
    TString samples[15] =  {
                            "Wantitop", "tW", 
                            "ttbar",
                            "WW",
                            "ZZ",
                            "WZ",
                            "dy50to120", "dy120to200", "dy200to400", "dy400to800", "dy800to1400", "dy1400to2300", "dy2300to3500", "dy3500to4500", "dy4500to6000",
                            };
	
	float events[15] = {
						7780870, 7581624,
						33844772,
						7791498,
						1949768,
						3928630,
						2961000, 100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000,
					};
						
	float sigma[15] = {
						35.6, 35.6,
						831.76,
						118.7,
						16.523,
						47.13,
						2112.905, 20.553, 2.8861, 0.25126, 0.017075, 1.366E-3, 8.178E-5, 3.191E-6, 2.787E-7,
						};         
         
	float LUMINOSITY = 41903.837;

	float weight[15] = {0};
	float weight_BB[15] = {0};
	float weight_BE[15] = {0};
	
	float gM;
	float kFactor;
	float kFactor_BB;
	float kFactor_BE;
	float NNPDF;


	
	Long64_t ne;
	float dil_mass;
	float cos_angle;
	float vertex_chi2;
	int dil_chosen;
	
	Int_t event;
	Int_t run;
	unsigned lumi;
	Int_t prev_event = -88;

	float lep_pt[2];
	float lep_phi[2];
	int lep_id[2];
	float lep_pt_err[2];
	float lep_eta[2];
	float lep_tk_pt[2];
	float lep_glb_pt[2];
	float lep_picky_pt[2];
	float lep_tpfms_pt[2];
	float lep_dB[2];
	float lep_sumPt[2];
	float lep_triggerMatchPt[2];
	short lep_glb_numberOfValidTrackerLayers[2]; 
	short lep_glb_numberOfValidPixelHits[2];
	short lep_glb_numberOfValidMuonHits[2];
	short lep_TuneP_numberOfValidMuonHits[2];
	short lep_numberOfMatchedStations[2];
  	short lep_expectedNnumberOfMatchedStations[2];
	bool lep_isGlobalMuon[2];
	bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool GoodVtx;
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    float vertex_m;
    float gen_dil_mass;
	float genWeight;    
    float gen_lep_qOverPt[2];
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
	int c_run = -1;
	int c_lumi = -1;
	int c_event = -1;
	float c_mass = -1;

	int n_run = -1;
	int n_lumi = -1;
	int n_event = -1;
	float n_mass = -1;
	
	for(int i = 0; i < 15; i++){
		weight[i] = LUMINOSITY * sigma[i] / events[i];
// 		weight[i] = 1;
		weight_BB[i] = weight[i];
		weight_BE[i] = weight[i];
	}
	
                        	
  for(int j=0; j < 15; j++){   

	 printf("openning.. %s %i  --- %.5f\n",samples[j].Data(),j , weight[j]);    

     TChain *treeMC = new TChain("SimpleNtupler/t");

     treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/MC/ana_datamc_" + samples[j] + ".root");        
//      treeMC->Add("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_2017/DATA_WITH_OR/ana_datamc_data.root");        

     treeMC->SetBranchAddress("genWeight",&genWeight);
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
     treeMC->SetBranchAddress("lep_expectedNnumberOfMatchedStations",lep_expectedNnumberOfMatchedStations);     
     treeMC->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
     treeMC->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
     treeMC->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
     treeMC->SetBranchAddress("lep_sumPt",lep_sumPt);
     treeMC->SetBranchAddress("lep_tk_pt",lep_tk_pt);
     treeMC->SetBranchAddress("lep_glb_pt",lep_glb_pt);
     treeMC->SetBranchAddress("lep_picky_pt",lep_picky_pt);
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

	int p_chosen = -99;

	std::cout<<"START: "<<ne<<std::endl;

	for ( int p=0; p<ne ;p++){
// 	for ( int p=0; p<1000 ;p++){

		if(p % 100000 == 0) std::cout<<p<<std::endl;		
		
		treeMC->GetEntry(p);
				
		if(
			GoodVtx && 
			fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
			lep_pt[0]>53. && lep_pt[1]>53. && 
			lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
			lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
			fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || (!(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2)))) &&
			(lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || (!(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2)))) &&
// 			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 			
			(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
			lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
			lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
			lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
			lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
			cos_angle>-0.9998 && 
			lep_id[0]*lep_id[1]<0 && 
			(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
			vertex_chi2 < 20
			){

				int counting = 0;
				bool trovato = false;
				int successivo = 0;

				if(dil_mass > - 99){			

					c_event = event;
					c_mass = dil_mass;
					int p_1 = p + 1;
// 					std::cout<<"A = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
					treeMC->GetEntry(p_1);
					n_event = event;
					while(c_event == n_event){	
						if(c_mass > 70 && c_mass < 110){
							successivo = 1;
							break;
						}
// 						std::cout<<"B = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
						if(
							GoodVtx && 
							fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 && 
							lep_pt[0]>53. && lep_pt[1]>53. && 
							lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 && 
							lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 && 
							fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 && 
							(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || (!(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2)))) &&
			(				lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || (!(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2)))) &&
// 			(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  && (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) && 
							(lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) && 
							lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 && 
							lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 && 
							lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 && 
							lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 && 
							cos_angle>-0.9998 && 
							lep_id[0]*lep_id[1]<0 && 
							(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && 
							vertex_chi2 < 20
						){					
							n_mass = dil_mass;
// 							counting_OK++;
							if(dil_mass > 70 && dil_mass < 110){
								p = p_1;
// 								std::cout<<"C = "<<p<<"\t"<<p_1<<"\t"<<c_event<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
								trovato = true;
								break;
							}
														
// 							std::cout<<"couting = "<<counting<<std::endl;
						}

						counting++;					
						p_1 = p_1 + 1;
						treeMC->GetEntry(p_1);
						n_event = event;
						
						if(p_1 > ne)
							break;

					}
				}

				treeMC->GetEntry(p);	

// 				std::cout<<"Selected = "<<p<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;

// 				if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
// 				prev_event = event;

// 				std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<"\t"<<lep_pt[0] + lep_pt[1]<<std::endl;
				
//				weight[j] *= genWeight;
//				weight_BB[j] *= genWeight;
//				weight_BE[j] *= genWeight;
				
//				gM = gen_dil_mass - 400;
//				kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
//			   	kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
//	 		   	kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
//	 		   	NNPDF = 0.9803 - 0.0001068 * gM + 1.142e-07 * pow(gM,2) -2.013e-11 * pow(gM,3)-5.904e-15 * pow(gM,4)+1.634e-18 * pow(gM,5);

//				if(j >= 6 && j <=14){
//					weight[j] *= kFactor*NNPDF;
//					weight_BB[j] *= kFactor_BB*NNPDF;
//					weight_BE[j] *= kFactor_BE*NNPDF;
//				}
				
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 




				std::cout<<dil_mass<<std::endl;














				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 
				////////////  INSERISCI QUI IL CODICE CHE SERVE DOPO AVER SELEZIONATO L'EVENTO //////////////
				//////////// //////////// //////////// //////////// //////////// //////////// //////////// 

//				if(j >= 6 && j <=14){
//					weight[j] /= kFactor;
//					weight_BB[j] /= kFactor_BB;
//					weight_BE[j] /= kFactor_BE;
//					weight[j] /= NNPDF;
//					weight_BB[j] /= NNPDF;
//					weight_BE[j] /= NNPDF;
//				}
				
//				weight[j] /= genWeight;
//				weight_BB[j] /= genWeight;
//				weight_BE[j] /= genWeight;
				
				if(!trovato) p += counting;
				p += successivo;

	       }//end of condition on event
	
		} // for event
	}
		
}
