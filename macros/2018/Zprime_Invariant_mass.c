
//
//  Zprime_Invariant_mass.c
//
//
//  Created by Kalpanie Liyanage on 10/28/19.
//


//  To Compare dimuon invariant mass distributions of all MC and DATA

#include <stdio.h>
#include "TFile.h"
#include <iostream>
#include <iomanip>
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
#include "TGaxis.h"
#include "TText.h"
#include <TString.h>



void DrawPlot(TH1F* DY, TH1F* ttbar, TH1F* tW, TH1F* Wantitop, TH1F* WW, TH1F* WZ, TH1F* ZZ, TH1F* Wjets,  /*Th1F* qcd, */TH1F* Zjets, TH1F* Zprime, TH1F* DATA, THStack* DATA_MC, TString title, TString name, bool logx);

void PrintKine(UInt_t data_kinemat[3500][20]);

void Check_ratio(TH1F* DY, TH1F* ttbar, TH1F* tW, TH1F* Wantitop, TH1F* WW, TH1F* WZ, TH1F* ZZ, TH1F* Wjets, TH1F* Zjets,  TH1F* DATA);

void Zprime_Invariant_mass(){
    
    TFile *f = new TFile("out.root", "RECREATE");
    //file1->cd();
    
    const int mc = 19;
    const int data = 4;
    const double PI = 3.141592654;
    
    
    UInt_t data_kinemat[3500][20]={0};
    
    
    
    
    
    
    float dy[8]={0.0},ttbar[8]={0.0},tW[8]={0.0},Wantitop[8]={0.0},WW[8]={0.0},WZ[8]={0.0},ZZ[8]={0.0},Zjets[8]={0.0},Wjets[8]={0.0},Zprime[8]={0.0};
    
    /////////////////////////////////////Data//////////////////////////////////
    
    TString samp[data] =  {
        "Run2018MuonsOnly_SingleMuonRun2018A_v2-17Sep2018",
        "Run2018MuonsOnly_SingleMuonRun2018B_v1-17Sep2018",
        "Run2018MuonsOnly_SingleMuonRun2018C_v1-17Sep2018",
        "Run2018MuonsOnly_SingleMuonRun2018D_v2-22Jan2019",
    };
    
    
    TString DATA_samples[data] =  {
        "Run2018A",
        "Run2018B",
        "Run2018C",
        "Run2018D",
    };
    
    ///////////////////////////////////////MC//////////////////////////////////
    
    TString MC_samples[mc] =  {
        "dy50to120",    //0
        "dy120to200",   //1
        "dy200to400",   //2
        "dy400to800",   //3
        "dy800to1400",  //4
        "dy1400to2300", //5
        "dy2300to3500", //6
        "dy3500to4500", //7
        "dy4500to6000", //8
        "dy6000toInf",  //9
        "ttbar",        //10
        "tW",           //11
        "Wantitop",     //12
        "WW",           //13
        "WZ",           //14
        "ZZ",           //15
        "Wjets",        //16
        "Zjets_filtered", //17
        "Zprime",       //18
        
    };
    
    double events[mc] = {
        2982000,    //0
        100000,     //1
        100000,     //2
        100000,     //3
        100000,     //4
        100000,     //5
        100000,     //6
        100000,     //7
        100000,     //8
        100000,     //9
        64310000,   //10    ttbar
        9598000,    //11    tW
        7623000,    //12    Wantitop
        7850000,    //13    WW
        3885000,    //14    WZ
        1979000,    //15    ZZ
        71026861,   //16    Wjets
        31470710,  //17    Zjets_filtered
        100000,     //18    Zprime
    };
    
    
    // using cross sections in AN_2018_011, XSDB & Alexander (most correct)
    float sigma[mc] = {
        2112.904,     //0
        20.55,      //1
        2.886,      //2
        0.2517,     //3
        0.01707,    //4
        1.366E-3,   //5
        8.178E-5,   //6
        3.191E-6,   //7
        2.787E-7,   //8
        9.569E-9,   //9
        88.29,      //10    ttbar
        35.6,      //11    tW
        35.6,      //12    Wantitop
        118.7,       //13    WW
        50.2,       //14    WZ
        16.523,      //15    ZZ
        61526.7,        //16    Wjets
        6077.22,        //17    Zjets
        6.76E-5,    //18
    };
    
    
    
    // From AN_2018_011
    /* float sigma[mc] = {
     2112.904,     //0
     20.553,      //1
     2.886,      //2
     0.2517,     //3
     0.01707,    //4
     1.366E-3,   //5
     8.178E-5,   //6
     3.191E-6,   //7
     2.787E-7,   //8
     9.569E-9,   //9
     87.31,      //10
     35.6,      //11
     35.6,      //12
     118.7,       //13
     47.13,       //14
     16.523,      //15
     //6529.0,
     //6525.42,
     6.76E-5,    //16
     };
     */
    
    
    
    
    float dil_mass;
    float vertex_m;
    float dil_pt;
    float dil_eta;
    float dil_rap;
    float dil_phi;
    float cos_angle;
    float vertex_chi2;
    int dil_chosen;
    float met_pt;
    float met_phi;
    
    float mt = -1.0; // for transverse mass
    float leading_pt = -1.0; // leading pt of the dimuon
    
    UInt_t event;
    UInt_t run;
    unsigned lumi;
    Int_t prev_event = -88;
    
    float lep_pt[2];
    float lep_phi[2];
    int lep_id[2];
    float lep_pt_err[2];
    float lep_eta[2];
    float lep_tk_pt[2];
    float lep_tk_dz[2];
    
    float lep_tuneP_pt[2];
    float lep_dyt_pt[2];
    float lep_std_pt[2];
    float lep_cocktail_pt[2];
    float lep_pfIso[2];
    
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
    //float vertex_m;
    float gen_dil_mass;
    float genWeight;
    float gen_lep_qOverPt[2];
    float lep_qOverPt[2]; //for data
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
    
    
    
    //ratios between DATA and MC in the range 60 - 120 GeV (for normalizing the background events to z peak), considering also the luminosity sections applied for all backgrounds
    float Z_peak = 1.0062;
    float Z_peak_BB = 1.0124;
    float Z_peak_BE = 1.0017;
    
    int mass_bins=200;
    
    
    //k-factor applied for DY (ttbar) samples only to accamodate NNPDF corrections and PI bkgs
    float gM;
    float kFactor;
    float kFactor_BB;
    float kFactor_BE;
    
    
    
    
    const double LUMINOSITY = 61298.77523; //in pb (14217.529374833 + 6867.781918123 + 6612.242975873 +  32717.848260849)
    
    float weight[mc] = {0.0};
    float weight_BB[mc] = {0.0};
    float weight_BE[mc] = {0.0};
    
    for(int i = 0; i<mc; i++){ //Normalizing MC to the luminosity of DATA
        weight[i] = LUMINOSITY * (float) sigma[i] / (float) events[i];
        weight[i] *=  (float) Z_peak;
        weight_BB[i] = weight[i] *  (float) Z_peak_BB /  (float) Z_peak;
        weight_BE[i] = weight[i] *  (float) Z_peak_BE/  (float) Z_peak;
        std::cout<<MC_samples[i]<<" "<<weight[i]<<"    "<<weight_BB[i]<<"    "<<weight_BE[i]<<std::endl;
    }
    
    
    const int    NMBINS = 51;
    const double MMIN = 70., MMAX =4000.;
    double logMbins[NMBINS+1];
    
    for (int ibin = 0; ibin <= NMBINS; ibin++){
        logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
        //       if(ibin != 0 && logMbins[ibin] <= logMbins[ibin-1]) std::cout<<"PROBLEMA"<<std::endl;
        //ss  std::cout<<logMbins[ibin]<<",";
    }
    
    
    
    
    //////////////////Stack plots////////////////////////////
    
    THStack* mass_linear = new THStack("mass_linear", "mass_linear");
    THStack* mass_linear_BB = new THStack("mass_linear_BB", "mass_linear_BB");
    THStack* mass_linear_BE = new THStack("mass_linear_BE", "mass_linear_BE");
    
    THStack* mass_log = new THStack("mass_log", "mass_log");
    THStack* mass_log_BB = new THStack("mass_log_BB", "mass_log_BB");
    THStack* mass_log_BE = new THStack("mass_log_BE", "mass_log_BE");
    
    THStack* mass_cumulative = new THStack("mass_cumulative", "mass_cumulative");
    THStack* mass_cumulative_BB = new THStack("mass_cumulative_BB", "mass_cumulative_BB");
    THStack* mass_cumulative_BE = new THStack("mass_cumulative_BE", "mass_cumulative_BE");
    
    THStack* mass_log_cumulative = new THStack("mass_log_cumulative", "mass_log_cumulative");
    THStack* mass_log_cumulative_BB = new THStack("mass_log_cumulative_BB", "mass_log_cumulative_BB");
    THStack* mass_log_cumulative_BE = new THStack("mass_log_cumulative_BE", "mass_log_cumulative_BE");
    
    
    
    
    
    TH1F* DY_mass_linear = new TH1F("DY_mass_linear", "DY_mass_linear;DY_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* DY_mass_linear_BB = new TH1F("DY_mass_linear_BB", "DY_mass_linear_BB;DY_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* DY_mass_linear_BE = new TH1F("DY_mass_linear_BE", "DY_mass_linear_BE;DY_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* DY_mass_log = new TH1F("DY_mass_log", "DY_mass_log;DY_mass_log;#Events", NMBINS, logMbins);
    TH1F* DY_mass_log_BB = new TH1F("DY_mass_log_BB", "DY_mass_log_BB;DY_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* DY_mass_log_BE = new TH1F("DY_mass_log_BE", "DY_mass_log_BE;DY_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* DY_mass_cumulative = new TH1F("DY_mass_cumulative", "DY_mass_cumulative;DY_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* DY_mass_cumulative_BB = new TH1F("DY_mass_cumulative_BB", "DY_mass_cumulative_BB;DY_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* DY_mass_cumulative_BE = new TH1F("DY_mass_cumulative_BE", "DY_mass_cumulative_BE;DY_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* DY_mass_log_cumulative = new TH1F("DY_mass_log_cumulative", "DY_mass_log_cumulative;DY_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* DY_mass_log_cumulative_BB = new TH1F("DY_mass_log_cumulative_BB", "DY_mass_log_cumulative_BB;DY_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* DY_mass_log_cumulative_BE = new TH1F("DY_mass_log_cumulative_BE", "DY_mass_log_cumulative_BE;DY_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    
    TH1F* ttbar_mass_linear = new TH1F("ttbar_mass_linear", "ttbar_mass_linear;ttbar_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* ttbar_mass_linear_BB = new TH1F("ttbar_mass_linear_BB", "ttbar_mass_linear_BB;ttbar_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* ttbar_mass_linear_BE = new TH1F("ttbar_mass_linear_BE", "ttbar_mass_linear_BE;ttbar_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* ttbar_mass_log = new TH1F("ttbar_mass_log", "ttbar_mass_log;ttbar_mass_log;#Events", NMBINS, logMbins);
    TH1F* ttbar_mass_log_BB = new TH1F("ttbar_mass_log_BB", "ttbar_mass_log_BB;ttbar_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* ttbar_mass_log_BE = new TH1F("ttbar_mass_log_BE", "ttbar_mass_log_BE;ttbar_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* ttbar_mass_cumulative = new TH1F("ttbar_mass_cumulative", "ttbar_mass_cumulative;ttbar_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* ttbar_mass_cumulative_BB = new TH1F("ttbar_mass_cumulative_BB", "ttbar_mass_cumulative_BB;ttbar_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* ttbar_mass_cumulative_BE = new TH1F("ttbar_mass_cumulative_BE", "ttbar_mass_cumulative_BE;ttbar_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* ttbar_mass_log_cumulative = new TH1F("ttbar_mass_log_cumulative", "ttbar_mass_log_cumulative;ttbar_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* ttbar_mass_log_cumulative_BB = new TH1F("ttbar_mass_log_cumulative_BB", "ttbar_mass_log_cumulative_BB;ttbar_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* ttbar_mass_log_cumulative_BE = new TH1F("ttbar_mass_log_cumulative_BE", "ttbar_mass_log_cumulative_BE;ttbar_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    TH1F* tW_mass_linear = new TH1F("tW_mass_linear", "tW_mass_linear;tW_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* tW_mass_linear_BB = new TH1F("tW_mass_linear_BB", "tW_mass_linear_BB;tW_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* tW_mass_linear_BE = new TH1F("tW_mass_linear_BE", "tW_mass_linear_BE;tW_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* tW_mass_log = new TH1F("tW_mass_log", "tW_mass_log;tW_mass_log;#Events", NMBINS, logMbins);
    TH1F* tW_mass_log_BB = new TH1F("tW_mass_log_BB", "tW_mass_log_BB;tW_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* tW_mass_log_BE = new TH1F("tW_mass_log_BE", "tW_mass_log_BE;tW_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* tW_mass_cumulative = new TH1F("tW_mass_cumulative", "tW_mass_cumulative;tW_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* tW_mass_cumulative_BB = new TH1F("tW_mass_cumulative_BB", "tW_mass_cumulative_BB;tW_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* tW_mass_cumulative_BE = new TH1F("tW_mass_cumulative_BE", "tW_mass_cumulative_BE;tW_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* tW_mass_log_cumulative = new TH1F("tW_mass_log_cumulative", "tW_mass_log_cumulative;tW_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* tW_mass_log_cumulative_BB = new TH1F("tW_mass_log_cumulative_BB", "tW_mass_log_cumulative_BB;tW_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* tW_mass_log_cumulative_BE = new TH1F("tW_mass_log_cumulative_BE", "tW_mass_log_cumulative_BE;tW_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    TH1F* Wantitop_mass_linear = new TH1F("Wantitop_mass_linear", "Wantitop_mass_linear;Wantitop_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* Wantitop_mass_linear_BB = new TH1F("Wantitop_mass_linear_BB", "Wantitop_mass_linear_BB;Wantitop_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* Wantitop_mass_linear_BE = new TH1F("Wantitop_mass_linear_BE", "Wantitop_mass_linear_BE;Wantitop_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Wantitop_mass_log = new TH1F("Wantitop_mass_log", "Wantitop_mass_log;Wantitop_mass_log;#Events", NMBINS, logMbins);
    TH1F* Wantitop_mass_log_BB = new TH1F("Wantitop_mass_log_BB", "Wantitop_mass_log_BB;Wantitop_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* Wantitop_mass_log_BE = new TH1F("Wantitop_mass_log_BE", "Wantitop_mass_log_BE;Wantitop_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* Wantitop_mass_cumulative = new TH1F("Wantitop_mass_cumulative", "Wantitop_mass_cumulative;Wantitop_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* Wantitop_mass_cumulative_BB = new TH1F("Wantitop_mass_cumulative_BB", "Wantitop_mass_cumulative_BB;Wantitop_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* Wantitop_mass_cumulative_BE = new TH1F("Wantitop_mass_cumulative_BE", "Wantitop_mass_cumulative_BE;Wantitop_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Wantitop_mass_log_cumulative = new TH1F("Wantitop_mass_log_cumulative", "Wantitop_mass_log_cumulative;Wantitop_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* Wantitop_mass_log_cumulative_BB = new TH1F("Wantitop_mass_log_cumulative_BB", "Wantitop_mass_log_cumulative_BB;Wantitop_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* Wantitop_mass_log_cumulative_BE = new TH1F("Wantitop_mass_log_cumulative_BE", "Wantitop_mass_log_cumulative_BE;Wantitop_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    
    TH1F* WW_mass_linear = new TH1F("WW_mass_linear", "WW_mass_linear;WW_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* WW_mass_linear_BB = new TH1F("WW_mass_linear_BB", "WW_mass_linear_BB;WW_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* WW_mass_linear_BE = new TH1F("WW_mass_linear_BE", "WW_mass_linear_BE;WW_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* WW_mass_log = new TH1F("WW_mass_log", "WW_mass_log;WW_mass_log;#Events", NMBINS, logMbins);
    TH1F* WW_mass_log_BB = new TH1F("WW_mass_log_BB", "WW_mass_log_BB;WW_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* WW_mass_log_BE = new TH1F("WW_mass_log_BE", "WW_mass_log_BE;WW_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* WW_mass_cumulative = new TH1F("WW_mass_cumulative", "WW_mass_cumulative;WW_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* WW_mass_cumulative_BB = new TH1F("WW_mass_cumulative_BB", "WW_mass_cumulative_BB;WW_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* WW_mass_cumulative_BE = new TH1F("WW_mass_cumulative_BE", "WW_mass_cumulative_BE;WW_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* WW_mass_log_cumulative = new TH1F("WW_mass_log_cumulative", "WW_mass_log_cumulative;WW_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* WW_mass_log_cumulative_BB = new TH1F("WW_mass_log_cumulative_BB", "WW_mass_log_cumulative_BB;WW_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* WW_mass_log_cumulative_BE = new TH1F("WW_mass_log_cumulative_BE", "WW_mass_log_cumulative_BE;WW_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    
    TH1F* WZ_mass_linear = new TH1F("WZ_mass_linear", "WZ_mass_linear;WZ_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* WZ_mass_linear_BB = new TH1F("WZ_mass_linear_BB", "WZ_mass_linear_BB;WZ_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* WZ_mass_linear_BE = new TH1F("WZ_mass_linear_BE", "WZ_mass_linear_BE;WZ_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* WZ_mass_log = new TH1F("WZ_mass_log", "WZ_mass_log;WZ_mass_log;#Events", NMBINS, logMbins);
    TH1F* WZ_mass_log_BB = new TH1F("WZ_mass_log_BB", "WZ_mass_log_BB;WZ_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* WZ_mass_log_BE = new TH1F("WZ_mass_log_BE", "WZ_mass_log_BE;WZ_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* WZ_mass_cumulative = new TH1F("WZ_mass_cumulative", "WZ_mass_cumulative;WZ_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* WZ_mass_cumulative_BB = new TH1F("WZ_mass_cumulative_BB", "WZ_mass_cumulative_BB;WZ_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* WZ_mass_cumulative_BE = new TH1F("WZ_mass_cumulative_BE", "WZ_mass_cumulative_BE;WZ_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* WZ_mass_log_cumulative = new TH1F("WZ_mass_log_cumulative", "WZ_mass_log_cumulative;WZ_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* WZ_mass_log_cumulative_BB = new TH1F("WZ_mass_log_cumulative_BB", "WZ_mass_log_cumulative_BB;WZ_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* WZ_mass_log_cumulative_BE = new TH1F("WZ_mass_log_cumulative_BE", "WZ_mass_log_cumulative_BE;WZ_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    
    TH1F* ZZ_mass_linear = new TH1F("ZZ_mass_linear", "ZZ_mass_linear;ZZ_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* ZZ_mass_linear_BB = new TH1F("ZZ_mass_linear_BB", "ZZ_mass_linear_BB;ZZ_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* ZZ_mass_linear_BE = new TH1F("ZZ_mass_linear_BE", "ZZ_mass_linear_BE;ZZ_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* ZZ_mass_log = new TH1F("ZZ_mass_log", "ZZ_mass_log;ZZ_mass_log;#Events", NMBINS, logMbins);
    TH1F* ZZ_mass_log_BB = new TH1F("ZZ_mass_log_BB", "ZZ_mass_log_BB;ZZ_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* ZZ_mass_log_BE = new TH1F("ZZ_mass_log_BE", "ZZ_mass_log_BE;ZZ_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* ZZ_mass_cumulative = new TH1F("ZZ_mass_cumulative", "ZZ_mass_cumulative;ZZ_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* ZZ_mass_cumulative_BB = new TH1F("ZZ_mass_cumulative_BB", "ZZ_mass_cumulative_BB;ZZ_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* ZZ_mass_cumulative_BE = new TH1F("ZZ_mass_cumulative_BE", "ZZ_mass_cumulative_BE;ZZ_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* ZZ_mass_log_cumulative = new TH1F("ZZ_mass_log_cumulative", "ZZ_mass_log_cumulative;ZZ_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* ZZ_mass_log_cumulative_BB = new TH1F("ZZ_mass_log_cumulative_BB", "ZZ_mass_log_cumulative_BB;ZZ_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* ZZ_mass_log_cumulative_BE = new TH1F("ZZ_mass_log_cumulative_BE", "ZZ_mass_log_cumulative_BE;ZZ_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    TH1F* Wjets_mass_linear = new TH1F("Wjets_mass_linear", "Wjets_mass_linear;Wjets_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* Wjets_mass_linear_BB = new TH1F("Wjets_mass_linear_BB", "Wjets_mass_linear_BB;Wjets_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* Wjets_mass_linear_BE = new TH1F("Wjets_mass_linear_BE", "Wjets_mass_linear_BE;Wjets_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Wjets_mass_log = new TH1F("Wjets_mass_log", "Wjets_mass_log;Wjets_mass_log;#Events", NMBINS, logMbins);
    TH1F* Wjets_mass_log_BB = new TH1F("Wjets_mass_log_BB", "Wjets_mass_log_BB;Wjets_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* Wjets_mass_log_BE = new TH1F("Wjets_mass_log_BE", "Wjets_mass_log_BE;Wjets_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* Wjets_mass_cumulative = new TH1F("Wjets_mass_cumulative", "Wjets_mass_cumulative;Wjets_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* Wjets_mass_cumulative_BB = new TH1F("Wjets_mass_cumulative_BB", "Wjets_mass_cumulative_BB;Wjets_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* Wjets_mass_cumulative_BE = new TH1F("Wjets_mass_cumulative_BE", "Wjets_mass_cumulative_BE;Wjets_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Wjets_mass_log_cumulative = new TH1F("Wjets_mass_log_cumulative", "Wjets_mass_log_cumulative;Wjets_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* Wjets_mass_log_cumulative_BB = new TH1F("Wjets_mass_log_cumulative_BB", "Wjets_mass_log_cumulative_BB;Wjets_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* Wjets_mass_log_cumulative_BE = new TH1F("Wjets_mass_log_cumulative_BE", "Wjets_mass_log_cumulative_BE;Wjets_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    TH1F* Zjets_mass_linear = new TH1F("Zjets_mass_linear", "Zjets_mass_linear;Zjets_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* Zjets_mass_linear_BB = new TH1F("Zjets_mass_linear_BB", "Zjets_mass_linear_BB;Zjets_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* Zjets_mass_linear_BE = new TH1F("Zjets_mass_linear_BE", "Zjets_mass_linear_BE;Zjets_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Zjets_mass_log = new TH1F("Zjets_mass_log", "Zjets_mass_log;Zjets_mass_log;#Events", NMBINS, logMbins);
    TH1F* Zjets_mass_log_BB = new TH1F("Zjets_mass_log_BB", "Zjets_mass_log_BB;Zjets_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* Zjets_mass_log_BE = new TH1F("Zjets_mass_log_BE", "Zjets_mass_log_BE;Zjets_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* Zjets_mass_cumulative = new TH1F("Zjets_mass_cumulative", "Zjets_mass_cumulative;Zjets_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* Zjets_mass_cumulative_BB = new TH1F("Zjets_mass_cumulative_BB", "Zjets_mass_cumulative_BB;Zjets_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* Zjets_mass_cumulative_BE = new TH1F("Zjets_mass_cumulative_BE", "Zjets_mass_cumulative_BE;Zjets_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Zjets_mass_log_cumulative = new TH1F("Zjets_mass_log_cumulative", "Zjets_mass_log_cumulative;Zjets_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* Zjets_mass_log_cumulative_BB = new TH1F("Zjets_mass_log_cumulative_BB", "Zjets_mass_log_cumulative_BB;Zjets_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* Zjets_mass_log_cumulative_BE = new TH1F("Zjets_mass_log_cumulative_BE", "Zjets_mass_log_cumulative_BE;Zjets_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    TH1F* Zprime_mass_linear = new TH1F("Zprime_mass_linear", "Zprime_mass_linear;Zprime_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* Zprime_mass_linear_BB = new TH1F("Zprime_mass_linear_BB", "Zprime_mass_linear_BB;Zprime_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* Zprime_mass_linear_BE = new TH1F("Zprime_mass_linear_BE", "Zprime_mass_linear_BE;Zprime_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Zprime_mass_log = new TH1F("Zprime_mass_log", "Zprime_mass_log;Zprime_mass_log;#Events", NMBINS, logMbins);
    TH1F* Zprime_mass_log_BB = new TH1F("Zprime_mass_log_BB", "Zprime_mass_log_BB;Zprime_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* Zprime_mass_log_BE = new TH1F("Zprime_mass_log_BE", "Zprime_mass_log_BE;Zprime_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* Zprime_mass_cumulative = new TH1F("Zprime_mass_cumulative", "Zprime_mass_cumulative;Zprime_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* Zprime_mass_cumulative_BB = new TH1F("Zprime_mass_cumulative_BB", "Zprime_mass_cumulative_BB;Zprime_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* Zprime_mass_cumulative_BE = new TH1F("Zprime_mass_cumulative_BE", "Zprime_mass_cumulative_BE;Zprime_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* Zprime_mass_log_cumulative = new TH1F("Zprime_mass_log_cumulative", "Zprime_mass_log_cumulative;Zprime_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* Zprime_mass_log_cumulative_BB = new TH1F("Zprime_mass_log_cumulative_BB", "Zprime_mass_log_cumulative_BB;Zprime_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* Zprime_mass_log_cumulative_BE = new TH1F("Zprime_mass_log_cumulative_BE", "Zprime_mass_log_cumulative_BE;Zprime_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    
    
    TH1F* DATA_mass_linear = new TH1F("DATA_mass_linear", "DATA_mass_linear;DATA_mass_linear;#Events", mass_bins, 70, 4000);
    TH1F* DATA_mass_linear_BB = new TH1F("DATA_mass_linear_BB", "DATA_mass_linear_BB;DATA_mass_linear_BB;#Events", mass_bins, 70, 4000);
    TH1F* DATA_mass_linear_BE = new TH1F("DATA_mass_linear_BE", "DATA_mass_linear_BE;DATA_mass_linear_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* DATA_mass_log = new TH1F("DATA_mass_log", "DATA_mass_log;DATA_mass_log;#Events", NMBINS, logMbins);
    TH1F* DATA_mass_log_BB = new TH1F("DATA_mass_log_BB", "DATA_mass_log_BB;DATA_mass_log_BB;#Events", NMBINS, logMbins);
    TH1F* DATA_mass_log_BE = new TH1F("DATA_mass_log_BE", "DATA_mass_log_BE;DATA_mass_log_BE;#Events", NMBINS, logMbins);
    
    TH1F* DATA_mass_cumulative = new TH1F("DATA_mass_cumulative", "DATA_mass_cumulative;DATA_mass_cumulative;#Events", mass_bins, 70, 4000);
    TH1F* DATA_mass_cumulative_BB = new TH1F("DATA_mass_cumulative_BB", "DATA_mass_cumulative_BB;DATA_mass_cumulative_BB;#Events", mass_bins, 70, 4000);
    TH1F* DATA_mass_cumulative_BE = new TH1F("DATA_mass_cumulative_BE", "DATA_mass_cumulative_BE;DATA_mass_cumulative_BE;#Events", mass_bins, 70, 4000);
    
    TH1F* DATA_mass_log_cumulative = new TH1F("DATA_mass_log_cumulative", "DATA_mass_log_cumulative;DATA_mass_log_cumulative;#Events", NMBINS, logMbins);
    TH1F* DATA_mass_log_cumulative_BB = new TH1F("DATA_mass_log_cumulative_BB", "DATA_mass_log_cumulative_BB;DATA_mass_log_cumulative_BB;#Events", NMBINS, logMbins);
    TH1F* DATA_mass_log_cumulative_BE = new TH1F("DATA_mass_log_cumulative_BE", "DATA_mass_log_cumulative_BE;DATA_mass_log_cumulative_BE;#Events", NMBINS, logMbins);
    
    
    
    
    int  x=0, q=0, r=0, s=0, t=0, u=0, v=0;
    
    
    
    
    //strat of looping over MC samples
    for(int j=0; j<mc; j++){
        
        std::cout<<"opening.. "<<MC_samples[j]<<" --- "<<" No of evnets: "<<events[j]<<" Cross section: "<<sigma[j]<<std::endl;
        
        TChain *treeMC = new TChain("SimpleNtupler/t");
        
        treeMC->Add("/eos/user/k/kaliyana/2018_MC/FileBased/"+ MC_samples[j] +"/"+ MC_samples[j] +".root");
        
        Long64_t ne = treeMC->GetEntries();
        
        treeMC->SetBranchAddress("genWeight",&genWeight);
        treeMC->SetBranchAddress("event",&event);
        treeMC->SetBranchAddress("run",&run);
        treeMC->SetBranchAddress("lumi",&lumi);
        treeMC->SetBranchAddress("dil_mass",&dil_mass);
        treeMC->SetBranchAddress("vertex_m",&vertex_m);
        treeMC->SetBranchAddress("dil_pt",&dil_pt);
        treeMC->SetBranchAddress("dil_eta",&dil_eta);
        treeMC->SetBranchAddress("dil_phi",&dil_phi);
        treeMC->SetBranchAddress("dil_rap",&dil_rap);
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
        
        treeMC->SetBranchAddress("met_pt",&met_pt);
        treeMC->SetBranchAddress("met_phi",&met_phi);
        treeMC->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
        treeMC->SetBranchAddress("lep_tk_dz",lep_tk_dz);
        
        
        std::cout<<"START: "<<ne<<std::endl;
        
        //start looping over entries
        for ( int p=0; p<ne ;p++){
            // if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            treeMC->GetEntry(p);
            
            //start selection on entries
            if(
               GoodVtx &&
               fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
               lep_pt[0]>53. && lep_pt[1]>53. &&  //offline reconstruction pt threshold
               lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 &&
               lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 &&
               fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 &&
               (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
               (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
               (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) &&
               lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 &&
               lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 &&
               lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
               lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 &&
               cos_angle>-0.9998 &&
               lep_id[0]*lep_id[1]<0 &&
               (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && //trigger pt threshold
               vertex_chi2 < 20
               //fabs(lep_tk_dz[0]) < 1.0 && fabs(lep_tk_dz[1]) < 1.0
               ){
                //std::cout<<vertex_m<<std::endl;
                if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
                prev_event = event;
                
                weight[j] *= genWeight;
                weight_BB[j] *= genWeight;
                weight_BE[j] *= genWeight;
                
                
                
                
                if(j==10){ //for ttbar
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<"\t"<<weight[j]<<std::endl;
                    
                    
                    gM = gen_dil_mass;
                    //kFactor =0.994078695151 + gM*2.64819793287e-05 - gM*gM*3.73996461024e-08 - gM*gM*gM*1.11452866827e-11;
                    kFactor = 0.991403 + 3.05593e-05 * gM - 2.21967e-07 * pow(gM,2) + 6.63658e-11 * pow(gM,3);
                    kFactor_BB = 0.990973 + 6.17124e-06 * gM - 3.31244e-07 * pow(gM,2) + 1.2125e-10 * pow(gM,3);
                    kFactor_BE = 0.990038 + 5.44269e-05 * gM - 2.43311e-07 * pow(gM,2) + 5.9748e-11 * pow(gM,3);
                    
                    weight[j] *=  (float) kFactor;
                    weight_BB[j] *=  (float) kFactor_BB;
                    weight_BE[j] *=  (float) kFactor_BE;
                    
                    ttbar_mass_linear->Fill(vertex_m, weight[j]);
                    ttbar_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        ttbar_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        ttbar_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        ttbar_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        ttbar_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    
                    if(vertex_m<=60.0) ttbar[0] = (float) ttbar[0]+ (float) weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) ttbar[1]= (float) ttbar[1]+ (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) ttbar[2]= (float) ttbar[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) ttbar[3]= (float) ttbar[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) ttbar[4]= (float) ttbar[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) ttbar[5]= (float) ttbar[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) ttbar[6]= (float) ttbar[6] + (float) weight[j];
                    if(vertex_m>=1800.0) ttbar[7]= (float) ttbar[7] + (float) weight[j];
                    
                    
                    weight[j] /=  (float) kFactor;
                    weight_BB[j] /=  (float) kFactor_BB;
                    weight_BE[j] /=  (float) kFactor_BE;
                }
                
                
                
                else if(j==11){ //for tW
                    // std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<"\t"<<weight[j]<<std::endl;
                    
                    
                    tW_mass_linear->Fill(vertex_m, weight[j]);
                    tW_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        tW_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        tW_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        tW_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        tW_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) tW[0] = tW[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) tW[1]= (float) tW[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) tW[2]= (float) tW[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) tW[3]= (float) tW[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) tW[4]= (float) tW[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) tW[5]= (float) tW[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) tW[6]= (float) tW[6] + (float) weight[j];
                    if(vertex_m>=1800.0) tW[7]= (float) tW[7] + (float) weight[j];
                }
                
                
                else if(j==12){ //for Wantitop
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<"\t"<<weight[j]<<std::endl;
                    
                    
                    Wantitop_mass_linear->Fill(vertex_m, weight[j]);
                    Wantitop_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Wantitop_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        Wantitop_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        Wantitop_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        Wantitop_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) Wantitop[0] = Wantitop[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) Wantitop[1]= (float) Wantitop[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) Wantitop[2]= (float) Wantitop[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) Wantitop[3]= (float) Wantitop[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) Wantitop[4]= (float) Wantitop[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) Wantitop[5]= (float) Wantitop[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) Wantitop[6]= (float) Wantitop[6] + (float) weight[j];
                    if(vertex_m>=1800.0) Wantitop[7]= (float) Wantitop[7] + (float) weight[j];
                }
                
                
                
                else if(j==13){ //for WW
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<std::endl;
                    
                    
                    WW_mass_linear->Fill(vertex_m, weight[j]);
                    WW_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        WW_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        WW_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        WW_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        WW_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) WW[0] = WW[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) WW[1]= (float) WW[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) WW[2]= (float) WW[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) WW[3]= (float) WW[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) WW[4]= (float) WW[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) WW[5]= (float) WW[5] +  (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) WW[6]= (float) WW[6] + (float) weight[j];
                    if(vertex_m>=1800.0) WW[7]= (float) WW[7] + (float) weight[j];
                }
                
                else if(j==14){ //for WZ
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<std::endl;
                    
                    
                    WZ_mass_linear->Fill(vertex_m, weight[j]);
                    WZ_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        WZ_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        WZ_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        WZ_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        WZ_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) WZ[0] = WZ[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) WZ[1]= (float) WZ[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) WZ[2]= (float) WZ[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) WZ[3]= (float) WZ[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) WZ[4]= (float) WZ[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) WZ[5]= (float) WZ[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) WZ[6]= (float) WZ[6] + (float) weight[j];
                    if(vertex_m>=1800.0) WZ[7]= (float) WZ[7] + (float) weight[j];
                }
                
                
                else if(j==15){ //for ZZ
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<std::endl;
                    
                    
                    ZZ_mass_linear->Fill(vertex_m, weight[j]);
                    ZZ_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        ZZ_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        ZZ_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        ZZ_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        ZZ_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) ZZ[0] = ZZ[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) ZZ[1]= (float) ZZ[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) ZZ[2]= (float) ZZ[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) ZZ[3]= (float) ZZ[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) ZZ[4]= (float) ZZ[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) ZZ[5]= (float) ZZ[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) ZZ[6]= (float) ZZ[6] + (float) weight[j];
                    if(vertex_m>=1800.0) ZZ[7]= (float) ZZ[7] + (float) weight[j];
                }
                
                
                else if(j==16){ //for Wjets
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<"\t"<<weight[j]<<std::endl;
                    
                    Wjets_mass_linear->Fill(vertex_m, weight[j]);
                    Wjets_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Wjets_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        Wjets_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        Wjets_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        Wjets_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) Wjets[0] = Wjets[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) Wjets[1]= (float) Wjets[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) Wjets[2]= (float) Wjets[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) Wjets[3]= (float) Wjets[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) Wjets[4]= (float) Wjets[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) Wjets[5]= (float) Wjets[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) Wjets[6]= (float) Wjets[6] + (float) weight[j];
                    if(vertex_m>=1800.0) Wjets[7]= (float) Wjets[7] + (float) weight[j];
                }
                
                
                else if(j==17){ //for Zjets
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<"\t"<<weight[j]<<std::endl;
                    
                    Zjets_mass_linear->Fill(vertex_m, weight[j]);
                    Zjets_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Zjets_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        Zjets_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        Zjets_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        Zjets_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) Zjets[0] = Zjets[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) Zjets[1]= (float) Zjets[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) Zjets[2]= (float) Zjets[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) Zjets[3]= (float) Zjets[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) Zjets[4]= (float) Zjets[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) Zjets[5]= (float) Zjets[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) Zjets[6]= (float) Zjets[6] + (float) weight[j];
                    if(vertex_m>=1800.0) Zjets[7]= (float) Zjets[7] + (float) weight[j];
                }
                
                
                else if(j==18){ //for Zprime
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<"\t"<<weight[j]<<std::endl;
                    
                    Zprime_mass_linear->Fill(vertex_m, weight[j]);
                    Zprime_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Zprime_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        Zprime_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        Zprime_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        Zprime_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) ttbar[0] = ttbar[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) Zprime[1]= (float) Zprime[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) Zprime[2]= (float) Zprime[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) Zprime[3]= (float) Zprime[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) Zprime[4]= (float) Zprime[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) Zprime[5]= (float) Zprime[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) Zprime[6]= (float) Zprime[6] + (float) weight[j];
                    if(vertex_m>=1800.0) Zprime[7]= (float) Zprime[7] + (float) weight[j];
                }
                
                
                else{ //for DY
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<vertex_m<<std::endl;
                    
                    
                    gM = gen_dil_mass - 400.0;
                    kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
                    kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
                    kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
                    
                    weight[j] *=  (float) kFactor;
                    weight_BB[j] *=  (float) kFactor_BB;
                    weight_BE[j] *=  (float) kFactor_BE;
                    
                    DY_mass_linear->Fill(vertex_m, weight[j]);
                    DY_mass_log->Fill(vertex_m, weight[j]);
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        DY_mass_linear_BB->Fill(vertex_m,  weight_BB[j]);
                        DY_mass_log_BB->Fill(vertex_m,  weight_BB[j]);
                    }
                    else{
                        DY_mass_linear_BE->Fill(vertex_m,  weight_BE[j]);
                        DY_mass_log_BE->Fill(vertex_m,  weight_BE[j]);
                    }
                    
                    if(vertex_m<=60.0) dy[0] = dy[0]+weight[j];
                    if(vertex_m>=60.0 && vertex_m<120.0) dy[1]= (float) dy[1] + (float) weight[j];
                    if(vertex_m>=120.0 && vertex_m<400.0) dy[2]= (float) dy[2] + (float) weight[j];
                    if(vertex_m>=400.0 && vertex_m<600.0) dy[3]= (float) dy[3] + (float) weight[j];
                    if(vertex_m>=600.0 && vertex_m<900.0) dy[4]= (float) dy[4] + (float) weight[j];
                    if(vertex_m>=900.0 && vertex_m<1300.0) dy[5]= (float) dy[5] + (float) weight[j];
                    if(vertex_m>=1300.0 && vertex_m<1800.0) dy[6]= (float) dy[6] + (float) weight[j];
                    if(vertex_m>=1800.0) dy[7]= (float) dy[7] + (float) weight[j];
                    
                    weight[j] /=  (float) kFactor;
                    weight_BB[j] /=  (float) kFactor_BB;
                    weight_BE[j] /=  (float) kFactor_BE;
                    
                }
                
                
                
                weight[j] /= genWeight;
                weight_BB[j] /= genWeight;
                weight_BE[j] /= genWeight;
                
                /*
                 // For yeild counting
                 if(vertex_m>=60.0 && vertex_m<120.0) x++;
                 if(vertex_m>=120.0 && vertex_m<400.0) q++;
                 if(vertex_m>=400.0 && vertex_m<600.0) r++;
                 if(vertex_m>=600.0 && vertex_m<900.0) s++;
                 if(vertex_m>=900.0 && vertex_m<1300.0) t++;
                 if(vertex_m>=1300.0 && vertex_m<1800.0) u++;
                 if(vertex_m>=1800.0) v++;
                 */
                
                
            } //end selection on entries
            
        }//end looping over entries
        
    } //end looping over MC samples
    
    
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = DY_mass_linear->Integral(i, mass_bins);
        DY_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = DY_mass_linear_BB->Integral(i, mass_bins);
        DY_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DY_mass_linear_BE->Integral(i, mass_bins);
        DY_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = ttbar_mass_linear->Integral(i, mass_bins);
        ttbar_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = ttbar_mass_linear_BB->Integral(i, mass_bins);
        ttbar_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ttbar_mass_linear_BE->Integral(i, mass_bins);
        ttbar_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = tW_mass_linear->Integral(i, mass_bins);
        tW_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = tW_mass_linear_BB->Integral(i, mass_bins);
        tW_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = tW_mass_log_BE->Integral(i, mass_bins);
        tW_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Wantitop_mass_linear->Integral(i, mass_bins);
        Wantitop_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wantitop_mass_linear_BB->Integral(i, mass_bins);
        Wantitop_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wantitop_mass_log_BE->Integral(i, mass_bins);
        Wantitop_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = WW_mass_linear->Integral(i, mass_bins);
        WW_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = WW_mass_linear_BB->Integral(i, mass_bins);
        WW_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WW_mass_linear_BE->Integral(i, mass_bins);
        WW_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = WZ_mass_linear->Integral(i, mass_bins);
        WZ_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = WZ_mass_linear_BB->Integral(i, mass_bins);
        WZ_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WZ_mass_linear_BE->Integral(i, mass_bins);
        WZ_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = ZZ_mass_linear->Integral(i, mass_bins);
        ZZ_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = ZZ_mass_linear_BB->Integral(i, mass_bins);
        ZZ_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ZZ_mass_linear_BE->Integral(i, mass_bins);
        ZZ_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Wjets_mass_linear->Integral(i, mass_bins);
        Wjets_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wjets_mass_linear_BB->Integral(i, mass_bins);
        Wjets_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wjets_mass_linear_BE->Integral(i, mass_bins);
        Wjets_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Zjets_mass_linear->Integral(i, mass_bins);
        Zjets_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zjets_mass_linear_BB->Integral(i, mass_bins);
        Zjets_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zjets_mass_linear_BE->Integral(i, mass_bins);
        Zjets_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = Zprime_mass_linear->Integral(i, mass_bins);
        Zprime_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zprime_mass_linear_BB->Integral(i, mass_bins);
        Zprime_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zprime_mass_linear_BE->Integral(i, mass_bins);
        Zprime_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    ////////////////////////////////////////////////////////////////////////
    
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = DY_mass_log->Integral(i, NMBINS);
        DY_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = DY_mass_log_BB->Integral(i, NMBINS);
        DY_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DY_mass_log_BE->Integral(i, NMBINS);
        DY_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = ttbar_mass_log->Integral(i, NMBINS);
        ttbar_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = ttbar_mass_log_BB->Integral(i, NMBINS);
        ttbar_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ttbar_mass_log_BE->Integral(i, NMBINS);
        ttbar_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = tW_mass_log->Integral(i, NMBINS);
        tW_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = tW_mass_log_BB->Integral(i, NMBINS);
        tW_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = tW_mass_log_BE->Integral(i, NMBINS);
        tW_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Wantitop_mass_log->Integral(i, NMBINS);
        Wantitop_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wantitop_mass_log_BB->Integral(i, NMBINS);
        Wantitop_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wantitop_mass_log_BE->Integral(i, NMBINS);
        Wantitop_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = WW_mass_log->Integral(i, NMBINS);
        WW_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = WW_mass_log_BB->Integral(i, NMBINS);
        WW_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WW_mass_log_BE->Integral(i, NMBINS);
        WW_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = WZ_mass_log->Integral(i, NMBINS);
        WZ_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = WZ_mass_log_BB->Integral(i, NMBINS);
        WZ_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = WZ_mass_log_BE->Integral(i, NMBINS);
        WZ_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = ZZ_mass_log->Integral(i, NMBINS);
        ZZ_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = ZZ_mass_log_BB->Integral(i, NMBINS);
        ZZ_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = ZZ_mass_log_BE->Integral(i, NMBINS);
        ZZ_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Wjets_mass_log->Integral(i, NMBINS);
        Wjets_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wjets_mass_log_BB->Integral(i, NMBINS);
        Wjets_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Wjets_mass_log_BE->Integral(i, NMBINS);
        Wjets_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Zjets_mass_log->Integral(i, NMBINS);
        Zjets_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zjets_mass_log_BB->Integral(i, NMBINS);
        Zjets_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zjets_mass_log_BE->Integral(i, NMBINS);
        Zjets_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = Zprime_mass_log->Integral(i, NMBINS);
        Zprime_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zprime_mass_log_BB->Integral(i, NMBINS);
        Zprime_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = Zprime_mass_log_BE->Integral(i, NMBINS);
        Zprime_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    
    
    
    
    
    
    ///////////////////////////////////DATA//////////////////////////
    
    int di = -1, a=0, b=0, c=0, d=0, e=0, h=0, g=0, k=0, l=0, m=0, n=0, w=0, y=0, z=0;
    
    //strat of looping over Data samples
    for(int j=0; j<data; j++){
        
        std::cout<<"opening.. "<<DATA_samples[j]<<std::endl;
        
        TChain *treeDATA = new TChain("SimpleNtupler/t");
        
        treeDATA->Add("/eos/user/k/kaliyana/2018_data/"+ samp[j] +"/"+ DATA_samples[j] +".root");
        Long64_t ne = treeDATA->GetEntries();
        
        //treeDATA->SetBranchAddress("genWeight",&genWeight);
        treeDATA->SetBranchAddress("event",&event);
        treeDATA->SetBranchAddress("run",&run);
        treeDATA->SetBranchAddress("lumi",&lumi);
        treeDATA->SetBranchAddress("dil_mass",&dil_mass);
        treeDATA->SetBranchAddress("vertex_m",&vertex_m);
        treeDATA->SetBranchAddress("dil_pt",&dil_pt);
        treeDATA->SetBranchAddress("dil_eta",&dil_eta);
        treeDATA->SetBranchAddress("dil_rap",&dil_rap);
        treeDATA->SetBranchAddress("dil_phi",&dil_phi);
        treeDATA->SetBranchAddress("cos_angle",&cos_angle);
        treeDATA->SetBranchAddress("vertex_chi2",&vertex_chi2);
        treeDATA->SetBranchAddress("dil_chosen",&dil_chosen);
        treeDATA->SetBranchAddress("lep_pt",lep_pt);
        treeDATA->SetBranchAddress("lep_id",lep_id);
        treeDATA->SetBranchAddress("lep_eta",lep_eta);
        treeDATA->SetBranchAddress("lep_phi",lep_phi);
        treeDATA->SetBranchAddress("lep_dB",lep_dB);
        treeDATA->SetBranchAddress("lep_triggerMatchPt",lep_triggerMatchPt);
        treeDATA->SetBranchAddress("lep_isTrackerMuon",lep_isTrackerMuon);
        treeDATA->SetBranchAddress("lep_isGlobalMuon",lep_isGlobalMuon);
        treeDATA->SetBranchAddress("lep_numberOfMatchedStations",lep_numberOfMatchedStations);
        treeDATA->SetBranchAddress("lep_expectedNnumberOfMatchedStations",lep_expectedNnumberOfMatchedStations);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidMuonHits",lep_glb_numberOfValidMuonHits);
        treeDATA->SetBranchAddress("lep_TuneP_numberOfValidMuonHits",lep_TuneP_numberOfValidMuonHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidPixelHits",lep_glb_numberOfValidPixelHits);
        treeDATA->SetBranchAddress("lep_glb_numberOfValidTrackerLayers",lep_glb_numberOfValidTrackerLayers);
        treeDATA->SetBranchAddress("lep_sumPt",lep_sumPt);
        treeDATA->SetBranchAddress("lep_tk_pt",lep_tk_pt);
        treeDATA->SetBranchAddress("lep_glb_pt",lep_glb_pt);
        treeDATA->SetBranchAddress("lep_picky_pt",lep_picky_pt);
        treeDATA->SetBranchAddress("lep_tpfms_pt",lep_tpfms_pt);
        treeDATA->SetBranchAddress("lep_pt_err",lep_pt_err);
        treeDATA->SetBranchAddress("vertex_m",&vertex_m);
        treeDATA->SetBranchAddress("GoodDataRan", &GoodDataRan);
        treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);
        treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
        treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
        //treeDATA->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
        treeDATA->SetBranchAddress("lep_qOverPt", lep_qOverPt);
        //treeDATA->SetBranchAddress("gen_lep_eta", gen_lep_eta);
        //treeDATA->SetBranchAddress("gen_lep_pt", gen_lep_pt);
        //treeDATA->SetBranchAddress("gen_lep_pt",gen_lep_pt);
        
        treeDATA->SetBranchAddress("met_pt",&met_pt);
        treeDATA->SetBranchAddress("met_phi",&met_phi);
        treeDATA->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
        treeDATA->SetBranchAddress("lep_tk_dz",lep_tk_dz);
        
        
        
        
        
        std::cout<<"START: "<<ne<<std::endl;
        
        //start looping over entries
        for ( int p=0; p<ne ;p++){
            //if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            treeDATA->GetEntry(p);
            
            //start selection on entries
            if(
                GoodVtx &&
               fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
               lep_pt[0]>53. && lep_pt[1]>53. && //offline reconstruction pt threshold
               lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 &&
               lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 &&
               fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 &&
               (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && (lep_expectedNnumberOfMatchedStations[0] < 2 || !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16) || lep_numberOfMatchedRPCLayers[0] > 2))) &&
               (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && (lep_expectedNnumberOfMatchedStations[1] < 2 || !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16) || lep_numberOfMatchedRPCLayers[1] > 2))) &&
               (lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0) &&
               lep_glb_numberOfValidPixelHits[0]>=1 && lep_glb_numberOfValidPixelHits[1]>=1 &&
               lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5 &&
               lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10 &&
               lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3 &&
               cos_angle>-0.9998 &&
               lep_id[0]*lep_id[1]<0 &&
               (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) && //offline reconstruction pt threshold
               vertex_chi2 < 20
               //fabs(lep_tk_dz[0]) < 1.0  && fabs(lep_tk_dz[1]) < 1.0
               ){
                
                if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
                
                //std::cout.setf(ios::fixed);
                //std::cout<<run<<"    "<<lumi<<"    "<<setprecision(0)<<event<<std::endl;
                
                prev_event = event;
                
                
                
                
                DATA_mass_linear->Fill(vertex_m);
                DATA_mass_log->Fill(vertex_m);
                
                if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                    DATA_mass_linear_BB->Fill(vertex_m);
                    DATA_mass_log_BB->Fill(vertex_m);
                }
                else{
                    DATA_mass_linear_BE->Fill(vertex_m);
                    DATA_mass_log_BE->Fill(vertex_m);
                }
                
                
                //For the table
                
                if(vertex_m>600){
                    
                    di++;
                    data_kinemat[di][0] = (UInt_t)di+1;
                    data_kinemat[di][1] = (UInt_t)vertex_m;
                    data_kinemat[di][14] = (UInt_t)met_pt;
                    data_kinemat[di][15] = (UInt_t)met_phi;
                    data_kinemat[di][16] = (UInt_t)run;
                    data_kinemat[di][17] = (UInt_t)lumi;
                    data_kinemat[di][18] = (UInt_t)event;
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        data_kinemat[di][19] = 1.0 ;
                    }
                    else data_kinemat[di][19] = 2.0;
                    
                    if(lep_id[0]<0 && lep_id[1]>0){
                        
                        data_kinemat[di][2] = (UInt_t)lep_id[0];
                        data_kinemat[di][3] = (UInt_t)lep_id[1];
                        
                        data_kinemat[di][4] = (UInt_t)lep_tuneP_pt[0];
                        data_kinemat[di][5] = (UInt_t)lep_tuneP_pt[1];
                        
                        data_kinemat[di][6] = (UInt_t)lep_tk_pt[0];
                        data_kinemat[di][7] = (UInt_t)lep_tk_pt[1];
                        
                        data_kinemat[di][8] = (UInt_t)lep_tuneP_pt[0]/lep_tk_pt[0];
                        data_kinemat[di][9] = (UInt_t)lep_tuneP_pt[1]/lep_tk_pt[1];
                        
                        data_kinemat[di][10] = (UInt_t)lep_eta[0];
                        data_kinemat[di][11] = (UInt_t)lep_phi[0];
                        
                        data_kinemat[di][12] = (UInt_t)lep_eta[1];
                        data_kinemat[di][13] = (UInt_t)lep_phi[1];
                        
                        
                    }
                    
                    else{
                        
                        data_kinemat[di][2] = (UInt_t)lep_id[1];
                        data_kinemat[di][3] = (UInt_t)lep_id[0];
                        
                        data_kinemat[di][4] = (UInt_t)lep_tuneP_pt[1];
                        data_kinemat[di][5] = (UInt_t)lep_tuneP_pt[0];
                        
                        data_kinemat[di][6] = (UInt_t)lep_tk_pt[1];
                        data_kinemat[di][7] = (UInt_t)lep_tk_pt[0];
                        
                        data_kinemat[di][8] = (UInt_t)lep_tuneP_pt[1]/lep_tk_pt[1];
                        data_kinemat[di][9] = (UInt_t)lep_tuneP_pt[0]/lep_tk_pt[0];
                        
                        data_kinemat[di][10] = (UInt_t)lep_eta[1];
                        data_kinemat[di][11] = (UInt_t)lep_phi[1];
                        
                        data_kinemat[di][12] = (UInt_t)lep_eta[0];
                        data_kinemat[di][13] = (UInt_t)lep_phi[0];
                    }
                } //end for table
                
                
                // For yeild counting
                if(vertex_m<60.0 && vertex_m<=120.0) a++;
                if(vertex_m>120.0 && vertex_m<=400.0) b++;
                if(vertex_m>400.0 && vertex_m<=600.0) c++;
                if(vertex_m>600.0 && vertex_m<=900.0) d++;
                if(vertex_m>900.0 && vertex_m<=1300.0) e++;
                if(vertex_m>1300.0 && vertex_m<=1800.0) h++;
                if(vertex_m>1800.0) g++;
                /*
                 // For yeild counting
                 if(dil_mass>60.0 && vertex_m<=120.0) k++;
                 if(dil_mass>120.0 && vertex_m<=400.0) l++;
                 if(dil_mass>400.0 && vertex_m<=600.0) m++;
                 if(dil_mass>600.0 && vertex_m<=900.0) n++;
                 if(dil_mass>900.0 && vertex_m<=1300.0) w++;
                 if(dil_mass>1300.0 && vertex_m<=1800.0) y++;
                 if(dil_mass>1800.0) z++;
                 */
                
            } //end selection on entries
            
        }//end looping over entries
        
    } //end looping over Data samples
    
    
    for(int i = 1; i <  mass_bins+1; i++){
        float contenuto_mc = DATA_mass_linear->Integral(i, mass_bins);
        DATA_mass_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = DATA_mass_linear_BB->Integral(i, mass_bins);
        DATA_mass_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DATA_mass_linear_BE->Integral(i, mass_bins);
        DATA_mass_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    for(int i = 1; i <  NMBINS+1; i++){
        float contenuto_mc = DATA_mass_log->Integral(i, NMBINS);
        DATA_mass_log_cumulative->SetBinContent(i, contenuto_mc);
        contenuto_mc = DATA_mass_log_BB->Integral(i, NMBINS);
        DATA_mass_log_cumulative_BB->SetBinContent(i, contenuto_mc);
        contenuto_mc = DATA_mass_log_BE->Integral(i, NMBINS);
        DATA_mass_log_cumulative_BE->SetBinContent(i, contenuto_mc);
    }
    
    
    f->Write();
    
    
    
    std::cout<<"Mass yeild"<<std::endl;
    
    std::cout<<"Data"<<std::endl;
    std::cout<<"a="<<a<<"    "<<"b="<<b<<"    "<<"c="<<c<<"    "<<"d="<<d<<"    "<<"e="<<e<<"    "<<"h="<<h<<"    "<<"g="<<g<<std::endl;
    //    std::cout<<"k="<<k<<"    "<<"l="<<l<<"    "<<"m="<<m<<"    "<<"n="<<n<<"    "<<"w="<<w<<"    "<<"y="<<y<<"    "<<"z="<<z<<std::endl;
    //std::cout<<"x="<<x<<"    "<<"q="<<q<<"    "<<"r="<<r<<"    "<<"s="<<s<<"    "<<"t="<<t<<"    "<<"u="<<u<<"    "<<"v="<<v<<std::endl;
    
    std::cout<<"DY"<<std::endl;
    std::cout<<dy[0]<<"    "<<dy[1]<<"    "<<dy[2]<<"    "<<dy[3]<<"    "<<dy[4]<<"    "<<dy[5]<<"    "<<dy[6]<<"    "<<dy[7]<<std::endl;
    
    std::cout<<"ttbar"<<std::endl;
    std::cout<<ttbar[0]<<"    "<<ttbar[1]<<"    "<<ttbar[2]<<"    "<<ttbar[3]<<"    "<<ttbar[4]<<"    "<<ttbar[5]<<"    "<<ttbar[6]<<"    "<<ttbar[7]<<std::endl;
    
    std::cout<<"tW"<<std::endl;
    std::cout<<tW[0]<<"    "<<tW[1]<<"    "<<tW[2]<<"    "<<tW[3]<<"    "<<tW[4]<<"    "<<tW[5]<<"    "<<tW[6]<<"    "<<tW[7]<<std::endl;
    
    std::cout<<"Wantitop"<<std::endl;
    std::cout<<Wantitop[0]<<"    "<<Wantitop[1]<<"    "<<Wantitop[2]<<"    "<<Wantitop[3]<<"    "<<Wantitop[4]<<"    "<<Wantitop[5]<<"    "<<Wantitop[6]<<"    "<<Wantitop[7]<<std::endl;
    
    std::cout<<"WW"<<std::endl;
    std::cout<<WW[0]<<"    "<<WW[1]<<"    "<<WW[2]<<"    "<<WW[3]<<"    "<<WW[4]<<"    "<<WW[5]<<"    "<<WW[6]<<"    "<<WW[7]<<std::endl;
    
    std::cout<<"WZ"<<std::endl;
    std::cout<<WZ[0]<<"    "<<WZ[1]<<"    "<<WZ[2]<<"    "<<WZ[3]<<"    "<<WZ[4]<<"    "<<WZ[5]<<"    "<<WZ[6]<<"    "<<WZ[7]<<std::endl;
    
    
    std::cout<<"ZZ"<<std::endl;
    std::cout<<ZZ[0]<<"    "<<ZZ[1]<<"    "<<ZZ[2]<<"    "<<ZZ[3]<<"    "<<ZZ[4]<<"    "<<ZZ[5]<<"    "<<ZZ[6]<<"    "<<ZZ[7]<<std::endl;
    
    std::cout<<"Wjets"<<std::endl;
    std::cout<<Wjets[0]<<"    "<<Wjets[1]<<"    "<<Wjets[2]<<"    "<<Wjets[3]<<"    "<<Wjets[4]<<"    "<<Wjets[5]<<"    "<<Wjets[6]<<"    "<<Wjets[7]<<std::endl;
    
    std::cout<<"Zjets"<<std::endl;
    std::cout<<Zjets[0]<<"    "<<Zjets[1]<<"    "<<Zjets[2]<<"    "<<Zjets[3]<<"    "<<Zjets[4]<<"    "<<Zjets[5]<<"    "<<Zjets[6]<<"    "<<Zjets[7]<<std::endl;
    
    std::cout<<"Zprime"<<std::endl;
    std::cout<<Zprime[0]<<"    "<<Zprime[1]<<"    "<<Zprime[2]<<"    "<<Zprime[3]<<"    "<<Zprime[4]<<"    "<<Zprime[5]<<"    "<<Zprime[6]<<"    "<<Zprime[7]<<std::endl;
    
    
    /////////////////////////////////////////////
    
    DrawPlot(DY_mass_log, ttbar_mass_log, tW_mass_log, Wantitop_mass_log, WW_mass_log, WZ_mass_log, ZZ_mass_log, Wjets_mass_log, Zjets_mass_log, Zprime_mass_log,  DATA_mass_log, mass_log, "Dimuon invariant mass(BB+BE+EE)",  "dimu_mass_log_scale",1);
    
    DrawPlot(DY_mass_log_BB, ttbar_mass_log_BB, tW_mass_log_BB, Wantitop_mass_log_BB, WW_mass_log_BB, WZ_mass_log_BB, ZZ_mass_log_BB, Wjets_mass_log_BB, Zjets_mass_log_BB, Zprime_mass_log_BB, DATA_mass_log_BB, mass_log_BB, "Dimuon invariant mass(BB)", "dimu_mass_log_scale_BB", 1);
    
    DrawPlot(DY_mass_log_BE, ttbar_mass_log_BE, tW_mass_log_BE, Wantitop_mass_log_BE, WW_mass_log_BE, WZ_mass_log_BE, ZZ_mass_log_BE, Wjets_mass_log_BE, Zjets_mass_log_BE, Zprime_mass_log_BE, DATA_mass_log_BE, mass_log_BE, "Dimuon invariant mass(BE+EE)", "dimu_mass_log_scale_BE+EE", 1);
    
    

    
    DrawPlot(DY_mass_linear, ttbar_mass_linear, tW_mass_linear, Wantitop_mass_linear, WW_mass_linear, WZ_mass_linear, ZZ_mass_linear, Wjets_mass_linear, Zjets_mass_linear, Zprime_mass_linear, DATA_mass_linear, mass_linear, "Dimuon invariant mass(BB+BE+EE)", "dimu_mass_linear_scale", 0);
    
    DrawPlot(DY_mass_linear_BB, ttbar_mass_linear_BB, tW_mass_linear_BB, Wantitop_mass_linear_BB, WW_mass_linear_BB, WZ_mass_linear_BB, ZZ_mass_linear_BB, Wjets_mass_linear_BB, Zjets_mass_linear_BB, Zprime_mass_linear_BB, DATA_mass_linear_BB, mass_linear_BB, "Dimuon invariant mass(BB)", "dimu_mass_linear_scale_BB", 0);
    
    DrawPlot(DY_mass_linear_BE, ttbar_mass_linear_BE, tW_mass_linear_BE, Wantitop_mass_linear_BE, WW_mass_linear_BE, WZ_mass_linear_BE, ZZ_mass_linear_BE, Wjets_mass_linear_BE, Zjets_mass_linear_BE, Zprime_mass_linear_BE, DATA_mass_linear_BE, mass_linear_BE, "Dimuon invariant mass(BE+EE)", "dimu_mass_linear_scale_BE+EE", 0);
    
    
    
    
    DrawPlot(DY_mass_cumulative, ttbar_mass_cumulative, tW_mass_cumulative, Wantitop_mass_cumulative, WW_mass_cumulative, WZ_mass_cumulative, ZZ_mass_cumulative, Wjets_mass_cumulative, Zjets_mass_cumulative, Zprime_mass_cumulative, DATA_mass_cumulative, mass_cumulative, "Dimuon invariant mass(BB+BE+EE)",  "dimu_mass_cumulative", 0);
    
    DrawPlot(DY_mass_cumulative_BB, ttbar_mass_cumulative_BB, tW_mass_cumulative_BB, Wantitop_mass_cumulative_BB, WW_mass_cumulative_BB, WZ_mass_cumulative_BB, ZZ_mass_cumulative_BB, Wjets_mass_cumulative_BB, Zjets_mass_cumulative_BB, Zprime_mass_cumulative_BB, DATA_mass_cumulative_BB, mass_cumulative_BB,  "Dimuon invariant mass(BB)","dimu_mass_cumulative_BB", 0);
    
    DrawPlot(DY_mass_cumulative_BE, ttbar_mass_cumulative_BE, tW_mass_cumulative_BE, Wantitop_mass_cumulative_BE, WW_mass_cumulative_BE, WZ_mass_cumulative_BE, ZZ_mass_cumulative_BE, Wjets_mass_cumulative_BE, Zjets_mass_cumulative_BE,  Zprime_mass_cumulative_BE,  DATA_mass_cumulative_BE, mass_cumulative_BE,  "Dimuon invariant mass(BE+EE)",  "dimu_mass_cumulative_BE+EE", 0);
    
    
    
    
    
    DrawPlot(DY_mass_log_cumulative, ttbar_mass_log_cumulative, tW_mass_log_cumulative, Wantitop_mass_log_cumulative, WW_mass_log_cumulative, WZ_mass_log_cumulative, ZZ_mass_log_cumulative, Wjets_mass_log_cumulative, Zjets_mass_log_cumulative, Zprime_mass_log_cumulative, DATA_mass_log_cumulative, mass_log_cumulative,  "Dimuon invariant mass(BB+BE+EE)", "dimu_mass_log_scale_cumulative", 1);
    
    DrawPlot(DY_mass_log_cumulative_BB, ttbar_mass_log_cumulative_BB, tW_mass_log_cumulative_BB, Wantitop_mass_log_cumulative_BB, WW_mass_log_cumulative_BB, WZ_mass_log_cumulative_BB, ZZ_mass_log_cumulative_BB, Wjets_mass_log_cumulative_BB, Zjets_mass_log_cumulative_BB, Zprime_mass_log_cumulative_BB, DATA_mass_log_cumulative_BB, mass_log_cumulative_BB, "Dimuon invariant mass(BB)", "dimu_mass_log_scale_cumulative_BB", 1);
    
    DrawPlot(DY_mass_log_cumulative_BE, ttbar_mass_log_cumulative_BE, tW_mass_log_cumulative_BE, Wantitop_mass_log_cumulative_BE, WW_mass_log_cumulative_BE, WZ_mass_log_cumulative_BE, ZZ_mass_log_cumulative_BE, Wjets_mass_log_cumulative_BE, Zjets_mass_log_cumulative_BE, Zprime_mass_log_cumulative_BE,  DATA_mass_log_cumulative_BE, mass_log_cumulative_BE, "Dimuon invariant mass(BE+EE)", "dimu_mass_log_scale_cumulative_BE+EE", 1);
    
    
    //    Check_ratio(DY_mass_log_BB, ttbar_mass_log_BB, tW_mass_log_BB, Wantitop_mass_log_BB, WW_mass_log_BB, WZ_mass_log_BB, ZZ_mass_log_BB, Zjets_mass_log_BB,  DATA_mass_log_BB);
    
    
    //    Check_ratio(DY_mass_log_BE, ttbar_mass_log_BE, tW_mass_log_BE, Wantitop_mass_log_BE, WW_mass_log_BE, WZ_mass_log_BE, ZZ_mass_log_BE, Zjets_mass_log_BE, DATA_mass_log_BE);
    
    
    PrintKine(data_kinemat);
    
    f->Close();
    
    
    
}





void DrawPlot(TH1F* DY, TH1F* ttbar, TH1F* tW, TH1F* Wantitop, TH1F* WW, TH1F* WZ, TH1F* ZZ, TH1F* Wjets,  /*Th1F* qcd, */TH1F* Zjets, TH1F* Zprime, TH1F* DATA, THStack* DATA_MC, TString title, TString name, bool logx){
    
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    gStyle->SetLegendTextSize(0.02);
    //gStyle->SetStatX(.9);
    // gStyle->SetStatY(.9);
    // gStyle->SetTitleFontSize(0.05);
    //   gStyle->SetLabelSize(0.01,"XY");
    
    
    //float r = 0.25;
    //float epsilon = 0.02;
    
    TCanvas *c3 = new TCanvas("DATA_MC", "DATA_MC",  1000, 1000);
    c3->cd();
    
    TPad* pad1 = new TPad("pad1", "pad1",  0, 0.3, 1, 1.0);
    pad1->SetLogy();
    if(logx) pad1->SetLogx();
    pad1->SetBottomMargin(0);// Upper and lower plot are joined
    pad1->Draw();
    pad1->cd();
    pad1->SetTicks();
    pad1->SetGrid();
    
    Zjets->SetLineColor(kBlack);
    Zjets->SetFillColor(kViolet-9);
    DATA_MC->Add(Zjets, "HIST");
    //Zjets->Draw("sameHIST");
    c3->Update();
    
    WW->SetLineColor(kBlack);
    WW->SetFillColor(kAzure+1);
    DATA_MC->Add(WW, "HIST");
    //WW->Draw("sameHIST");
    c3->Update();
    
    Wantitop->SetLineColor(kBlack);
    Wantitop->SetFillColor(kOrange-9);
    DATA_MC->Add(Wantitop, "HIST");
    //Wantitop->Draw("sameHIST");
    c3->Update();
    
    tW->SetLineColor(kBlack);
    tW->SetFillColor(kGreen-5);
    DATA_MC->Add(tW, "HIST");
    //tW->Draw("sameHIST");
    c3->Update();
    
    ZZ->SetLineColor(kBlack);
    ZZ->SetFillColor(kYellow-9);
    DATA_MC->Add(ZZ, "HIST");
    //ZZ->Draw("sameHIST");
    c3->Update();
    
    WZ->SetLineColor(kBlack);
    WZ->SetFillColor(kGreen);
    DATA_MC->Add(WZ, "HIST");
    //WZ->Draw("sameHIST");
    c3->Update();
    
    ttbar->SetLineColor(kBlack);
    ttbar->SetFillColor(kOrange-7);
    DATA_MC->Add(ttbar, "HIST");
    //ttbar->Draw("sameHIST");
    c3->Update();
    
    Wjets->SetLineColor(kBlack);
    Wjets->SetFillColor(kRed-4);
    DATA_MC->Add(Wjets, "HIST");
    //Zjets->Draw("sameHIST");
    c3->Update();
    
    DY->SetLineColor(kBlack);
    DY->SetFillColor(kYellow-8);
    DATA_MC->Add(DY, "HIST");
    //DY->Draw("sameHIST");
    c3->Update();
    
    
    DATA_MC->Draw();
    c3->Update();
    
    
    Zprime->SetLineColor(kRed);
    Zprime->SetLineWidth(2);
    Zprime->Draw("HIST SAME");
    //DATA_MC->Add(Zprime, "HIST");
    //Zprime->Draw("sameHIST");
    c3->Update();
    
    DATA->SetMarkerStyle(20);
    DATA->SetMarkerColor(kBlack);
    DATA->SetMarkerSize(0.95);
    DATA->Draw("SAMEPE");
    //DATA_MC->Add(DATA, "PE");
    c3->Update();
    
    
    DATA_MC->GetXaxis()->SetTitle(" ");
    DATA_MC->GetYaxis()->SetTitle("Events");
    DATA_MC->GetYaxis()->SetTitleSize(18);
    DATA_MC->GetYaxis()->SetTitleFont(43);
    DATA_MC->GetYaxis()->SetTitleOffset(1.55);
    DATA_MC->GetYaxis()->SetLabelFont(43);
    DATA_MC->GetYaxis()->SetLabelSize(15);
    DATA_MC->SetMinimum(0.00001);
    DATA_MC->GetXaxis()->SetLabelSize(0);
    DATA_MC->SetTitle(title);
    //DATA->Draw("samePE");
    gPad->RedrawAxis();//to redraw the axes on top of all histograms
    c3->Update();
    
    TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText2->SetBorderSize(0);
    tText2->SetFillColor(0);
    tText2->SetFillStyle(0);
    TText *t2 = tText2->AddText("CMS 2018 Data 61.3 fb^{-1} (13TeV)");
    tText2->SetTextSize(0.035);
    tText2->Draw();
    c3->Update();
    
    
    TLegend *l1 = new TLegend(0.65,0.55,0.85,0.85);
    //l1->SetBorderSize(0);
    l1->AddEntry(DATA, "DATA", "lep");
    l1->AddEntry(DY, "#gamma^{*}/Z #rightarrow #mu^{+}#mu^{-}", "f");
    l1->AddEntry(Wjets, "W+jets", "f");
    l1->AddEntry(ttbar, "t#bar{t}", "f");
    l1->AddEntry(WZ, "WZ", "f");
    l1->AddEntry(ZZ, "ZZ", "f");
    l1->AddEntry(tW, "tW", "f");
    l1->AddEntry(Wantitop, "#bar{t}W", "f");
    l1->AddEntry(WW, "WW", "f");
    l1->AddEntry(Zjets, "#gamma^{*}/Z #rightarrow #tau^{+}#tau^{-}", "f");
    l1->AddEntry(Zprime, "Zprime", "l");
    l1->Draw();
    c3->Update();
    c3->cd();
    
    TPad* pad2 = new TPad("pad2", "pad2",  0, 0.05, 1, 0.3);
    //pad2->SetFrameFillStyle(4000);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGrid();
    if(logx) pad2->SetLogx();
    pad2->Draw();
    pad2->cd();
    pad2->SetTicks();
    pad2->Update();
    
    
    
    
    TH1F* ratio = (TH1F*)DATA->Clone();
    TH1F* MC = (TH1F*)DY->Clone();
    MC->Add(Wjets);
    MC->Add(ttbar);
    MC->Add(tW);
    MC->Add(Wantitop);
    MC->Add(WW);
    MC->Add(WZ);
    MC->Add(ZZ);
    MC->Add(Zjets);
    ratio->Add(MC,-1);
    ratio->Divide(MC);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerSize(0.95);
    ratio->GetYaxis()->SetTitle("(DATA-Bkg) / Bkg");
    ratio->GetXaxis()->SetTitle("m(#mu^{+}#mu^{-})[GeV]");
    ratio->GetYaxis()->SetTitleSize(18);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(15);
    ratio->GetYaxis()->SetTitleOffset(1.55); //Distance between the axis and the title
    
    //ratio->GetXaxis()->SetTitle(title);
    ratio->GetXaxis()->SetTitleSize(18);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(15);
    ratio->GetXaxis()->SetTitleOffset(4.);
    ratio->GetYaxis()->SetRangeUser(-1.0, 1.0);
    //ratio->SetMaximum(2.0);
    //ratio->SetMinimum(-1.0);
    
    ratio->SetTitle("");
    ratio->Draw("PE");
    pad2->Update();
    c3->Update();
    
    //For X axis settings
    Float_t Xmin=pad2->GetUxmin();
    Float_t Xmax=pad2->GetUxmax();
    Float_t Ymin=pad2->GetUymin();
    Float_t Ymax=pad2->GetUymax();
    
    
    
    
    TLine* line = new TLine(Xmin, 0, Xmax, 0);
    line->SetLineColor(kGreen);
    line->SetLineWidth(1);
    if(logx) line->DrawLine(pow(10,Xmin), 0, pow(10,Xmax), 0);
    else line->DrawLine(Xmin, 0, Xmax, 0);
    pad2->Update();
    
    
    
    c3->Update();
    
    c3->SaveAs(name+".pdf","pdf");
    c3->SaveAs(name+".png","png");
    
    
    
}



void PrintKine(UInt_t array[3500][20]){
    
    std::cout<<"Printing muon kinematics"<<std::endl;
    
    std::cout<<"No"<<";"<<"Vertex Mass[GeV]"<<";"<<"lep_id_muon(-)"<<";"<<"lep_id_muon(+)"<<";"<<"TuneP_pt[GeV]_muon(-)"<<";"<<"TuneP_pt[GeV]_muon(+)"<<";"<<"Tracker_pt[GeV]_muon(-)"<<";"<<"Tracker_pt[GeV]_muon(+)"<<";"<<"Ratio(TuneP/Tracker)_muon(-)"<<";"<<"Ratio(TuneP/Tracker)_muon(+)"<<";"<<"eta_muon(-)"<<";"<<"eta_muon(+)"<<";"<<"phi_muon(-)"<<";"<<"phi_muon(+)"<<";"<<"MET_pt[GeV]"<<";"<<"MET_phi"<<";"<<"Run"<<";"<<"Lumi"<<";"<<"Event"<<";"<<"Category"<<std::endl;
    
    
    /* std::cout<<std::setw(10)<<"No"<<std::setw(10)<<"Mass[GeV]"<<std::setw(10)<<"TuneP_pt[GeV]_muon(-)"<<std::setw(10)<<"TuneP_pt[GeV]_muon(+)"<<std::setw(10)<<"Tracker_pt[GeV]_muon(-)"<<std::setw(10)<<"Tracker_pt[GeV]_muon(+)"<<std::setw(10)<<"Ratio(TuneP/Tracker)_muon(-)"<<std::setw(10)<<"Ratio(TuneP/Tracker)_muon(+)"<<std::setw(10)<<"eta_muon(-)"<<std::setw(10)<<"eta_muon(+)"<<std::setw(10)<<"phi_muon(-)"<<std::setw(10)<<"phi_muon(+)"<<std::setw(10)<<"MET_pt[GeV]"<<std::setw(10)<<"MET_phi"<<std::setw(10)<<"Category"<<std::setw(10)<<"Run"<<std::setw(10)<<"Lumi"<<std::setw(10)<<"Event"<<std::endl;
     
     */
    for(int r=0; r<3500; r++){
        
        for(int c=0; c<20; c++){
            
            std::cout.setf(ios::fixed);
            std::cout<<setprecision(0)<<array[r][c]<<";";
            
        }
        std::cout<<std::endl;
        
    }
    
}

/*
 void Check_ratio(TH1F* DY, TH1F* ttbar, TH1F* tW, TH1F* Wantitop, TH1F* WW, TH1F* WZ, TH1F* ZZ, TH1F* Zjets, TH1F* DATA){
 
 
 TH1F* MC = (TH1F*)DY->Clone();
 MC->Add(ttbar);
 MC->Add(tW);
 MC->Add(Wantitop);
 MC->Add(WW);
 MC->Add(WZ);
 MC->Add(ZZ);
 MC->Add(Zjets);
 int mass_bins=MC->GetNbinsX();
 
 for(int i=1; i<mass_bins+1; i++){
 std::cout<<i<<"    "<<MC->GetXaxis()->GetBinLowEdge(i)<<"  "<<MC->GetBinContent(i)<<"    "<<DATA->GetBinContent(i)<<"    "<<(float)(DATA->GetBinContent(i)-MC->GetBinContent(i))/MC->GetBinContent(i)<<std::endl;
 }
 
 
 
 }
 
 */



