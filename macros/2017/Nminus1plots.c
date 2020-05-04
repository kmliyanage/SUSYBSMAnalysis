////
//  N-1 plots.c
//
//
//  Created by Kalpanie Liyanage on 10/31/19.



//  To Compare kinematic distributions of all MC and DATA
// Comment the cut that you don't need and make plots
// Also can use to make basic distributons with all cuts applied for highPt selections

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



void DrawPlot(TH1F* DY, TH1F* ttbar, TH1F* tW, TH1F* Wantitop, TH1F* WW, TH1F* WZ, TH1F* ZZ, TH1F* Wjets, TH1F* Zjets, /*Th1F* qcd, */TH1F* Zprime, TH1F* DATA, THStack* DATA_MC, TString title, TString name, bool logx);

//void PrintKine(float array[600][20]);

void Nminus1plots(){
    
    TFile *f = new TFile("out.root", "RECREATE");
    //file1->cd();
    
    const int mc = 19;
    const int data = 5;
    const double PI = 3.141592654;
    
    
    float data_kinemat[600][20]={0};
    
    /////////////////////////////////////Data//////////////////////////////////
    
    TString samp[data] =  {
        "Run2017MuonsOnly_SingleMuonRun2017B-31Mar2018-v1",
        "Run2017MuonsOnly_SingleMuonRun2017C-17Nov2017-v1",
        "Run2017MuonsOnly_SingleMuonRun2017D-17Nov2017-v1",
        "Run2017MuonsOnly_SingleMuonRun2017E-17Nov2017-v1",
        "Run2017MuonsOnly_SingleMuonRun2017F-17Nov2017-v1",
    };
    
    
    TString DATA_samples[data] =  {
        "Run2017B",
        "Run2017C",
        "Run2017D",
        "Run2017E",
        "Run2017F",
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
        2863000,    //0    dy50to120
        100000,     //1    dy120to200
        100000,     //2    dy200to400
        100000,     //3    dy400to800
        100000,     //4    dy800to1400
        100000,     //5    dy1400to2300
        100000,     //6    dy2300to3500
        100000,     //7    dy3500to4500
        100000,     //8    dy4500to6000
        100000,     //9    dy6000toInf
        960752,   //10  ttbar
        7794186,    //11    tW
        7977430,    //12    Wantitop
        7765828,    //13    WW
        3928630,    //14    WZ
        1949768,    //15    ZZ
        30008250,   //16    Wjets
        15296830,  //17 Zjets_filtered
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
    
    // using cross sections in AN_2018_011
    /* float sigma[mc] = {
     2113,     //0
     20.55,      //1
     2.886,      //2
     0.2513,     //3
     0.01707,    //4
     1.366E-3,   //5
     8.178E-5,   //6
     3.191E-6,   //7
     2.787E-7,   //8
     9.569E-9,   //9
     88.29,      //10
     34.91,      //11
     34.97,      //12
     75.8,       //13
     27.6,       //14
     12.14,      //15
     52940.0,    //16    Wjets
     6077.22,    //17
     5.515E-5,    //18
     };
     */
    
    
    
    
    
    
    
    float dil_mass;
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
    float prev_lep_pt1 = 0.0;
    float prev_lep_pt2 = 0.0;
    float prev_vertex_m = 0.0;
    float prev_dil_pt = 0.0;
    
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
    float vertex_m;
    float gen_dil_mass;
    float genWeight;
    float gen_lep_qOverPt[2];
    float lep_qOverPt[2]; //for data
    float gen_lep_eta[2];
    float gen_lep_pt[2];
    
    
    
    
    //ratios between DATA and MC in the range 60 - 120 GeV (for normalizing the background events to z peak), considering also the luminosity sections applied for all backgrounds
    float Z_peak = 1.0282;
    float Z_peak_BB = 1.0286;
    float Z_peak_BE = 1.0278;
    
    int mass_bins=200;
    
    
    //k-factor applied for DY samples only to accamodate NNPDF corrections and PI bkgs
    float gM;
    float kFactor;
    float kFactor_BB;
    float kFactor_BE;
    
    
    const double LUMINOSITY = 42135.25562; //in pb
    
    float weight[mc] = {0.0};
    float weight_BB[mc] = {0.0};
    float weight_BE[mc] = {0.0};
    
    for(int i = 0; i<mc; i++){ //Normalizing MC to the luminosity of DATA
        weight[i] = LUMINOSITY * (float) sigma[i] / (float) events[i];
        weight_BB[i] = weight[i] * (float) Z_peak_BB ;
        weight_BE[i] = weight[i] * (float) Z_peak_BE;
        weight[i] *= Z_peak;
        std::cout<<MC_samples[i]<<" "<<weight[i]<<"    "<<weight_BB[i]<<"    "<<weight_BE[i]<<std::endl;
    }
    
    
    const int    NMBINS = 50;
    const double MMIN = 60., MMAX =6000.;
    double logMbins[NMBINS+1];
    
    for (int ibin = 0; ibin <= NMBINS; ibin++){
        logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
        //       if(ibin != 0 && logMbins[ibin] <= logMbins[ibin-1]) std::cout<<"PROBLEMA"<<std::endl;
        //ss  std::cout<<logMbins[ibin]<<",";
    }
    
    
    
    double pt_bins[9] = {50.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 1000.0, 3200.0};
    Int_t  binnum_pt = sizeof(pt_bins)/sizeof(Double_t)-1;
    
    //////////////////Stack plots////////////////////////////
    
    THStack* dilepton_pt = new THStack("dil_pt ", "dil_pt");
    THStack* dilepton_rapidity = new THStack("dil_rapidity ", "dil_rapidity");
    THStack* dilepton_rapidity_BB = new THStack("dil_rapidity_BB ", "dil_rapidity_BB");
    THStack* dilepton_rapidity_BE = new THStack("dil_rapidity_BE ", "dil_rapidity_BE");
    THStack* dilepton_phi = new THStack("dil_phi ", "dil_phi");
    THStack* lepton_pt = new THStack("lep_pt", "lep_pt");
    THStack* lepton_pt_BB = new THStack("lep_pt_BB ", "lep_pt_BB");
    THStack* lepton_pt_BE = new THStack("lep_pt_BE ", "lep_pt_BE");
    THStack* lepton_pt_plus = new THStack("lep_pt_plus ", "lep_pt_plus");
    THStack* lepton_pt_minus = new THStack("lep_pt_minus ", "lep_pt_minus");
    THStack* lepton_eta = new THStack("lep_eta", "lep_eta");
    THStack* lepton_eta_BB = new THStack("lep_eta_BB", "lep_eta_BB");
    THStack* lepton_eta_BE = new THStack("lep_eta_BE", "lep_eta_BE");
    THStack* lepton_phi = new THStack("lep_phi", "lep_phi");
    THStack* dB_clear = new THStack("lep_dB", "lep_dB");
    THStack* PixelHit_clear = new THStack("lep_PixelHit", "lep_PixelHit");
    THStack* TkLayer_clear = new THStack("lep_TkLayer", "lep_TkLayer");
    THStack* Iso_clear = new THStack("dil_Iso", "dil_Iso");
    THStack* relpTErr_clear = new THStack("dil_relpTErr", "dil_relpTErr");
    THStack* Vtx_clear = new THStack("dil_Vtx", "dil_Vtx");
    THStack* ValidMu_clear = new THStack("dil_ValidMu", "dil_ValidMu");
    
    
    
    ////////////////Histograms///////////////////////////////
    
    
    
    
    //////////////////////////////DY/////////////////////////////////////////////
    TH1F* DY_pt = new TH1F("dil_pt_DY ", "dil_pt_DY", 100, 0, 2000);
    TH1F* DY_rapidity = new TH1F("dil_rapidity_DY ", "dil_rapidity_DY", 50, -3, 3);
    TH1F* DY_rapidity_BB = new TH1F("dil_rapidity_DY_BB ", "dil_rapidity_DY_BB", 50, -3, 3);
    TH1F* DY_rapidity_BE = new TH1F("dil_rapidity_DY_BE ", "dil_rapidity_DY_BE", 50, -3, 3);
    TH1F* DY_phi = new TH1F("dil_phi_DY ", "dil_phi_DY", 50, -4, 4);
    TH1F* DY_lep_pt = new TH1F("lep_pt_DY ", "lep_pt_DY", 36, 50, 1850);
    DY_lep_pt->Sumw2();
    TH1F* DY_lep_pt_BB = new TH1F("lep_pt_DY_BB ", "lep_pt_DY_BB", 36, 50, 1850);
    DY_lep_pt_BB->Sumw2();
    TH1F* DY_lep_pt_BE = new TH1F("lep_pt_DY_BE ", "lep_pt_DY_BE", 36, 50, 1850);
    DY_lep_pt_BE->Sumw2();
    TH1F* DY_lep_pt_plus = new TH1F("lep_pt_DY_plus ", "lep_pt_DY_plus", 100, 0, 2000);
    TH1F* DY_lep_pt_minus = new TH1F("lep_pt_DY_minus ", "lep_pt_DY_minus", 100, 0, 2000);
    TH1F* DY_lep_eta = new TH1F("lep_eta_DY ", "lep_eta_DY", 50, -3, 3);
    TH1F* DY_lep_eta_BB = new TH1F("lep_eta_DY_BB ", "lep_eta_DY_BB", 50, -3, 3);
    TH1F* DY_lep_eta_BE = new TH1F("lep_eta_DY_BE ", "lep_eta_DY_BE", 50, -3, 3);
    TH1F* DY_lep_phi = new TH1F("lep_phi_DY ", "lep_phi_DY", 50, -4, 4);
    TH1F* dB_DY_clear = new TH1F("DY_dB", "DY_dB",  50, 0, 0.25);
    TH1F* PixelHit_DY_clear = new TH1F("DY_PixelHit", "DY_PixelHit",  15, 0,  15);
    TH1F* TkLayer_DY_clear = new TH1F("DY_TkLayer", "DY_TkLayer", 20, 0, 20);
    TH1F* Iso_DY_clear = new TH1F("DY_Iso", "DY_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_DY_clear = new TH1F("DY_relpTErr", "DY_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_DY_clear = new TH1F("DY_Vtx", "DY_Vtx",  44, 0, 22);
    TH1F* ValidMu_DY_clear = new TH1F("DY_ValidMu", "DY_ValidMu",  55, 0, 55);
    
    
    
    
    //////////////////////////////ttbar/////////////////////////////////////////////
    TH1F* ttbar_pt = new TH1F("dil_pt_ttbar ", "dil_pt_ttbar", 100, 0, 2000);
    TH1F* ttbar_rapidity = new TH1F("dil_rapidity_ttbar ", "dil_rapidity_ttbar", 50, -3, 3);
    TH1F* ttbar_rapidity_BB = new TH1F("dil_rapidity_ttbar_BB ", "dil_rapidity_ttbar_BB", 50, -3, 3);
    TH1F* ttbar_rapidity_BE = new TH1F("dil_rapidity_ttbar_BE ", "dil_rapidity_ttbar_BE", 50, -3, 3);
    TH1F* ttbar_phi = new TH1F("dil_phi_ttbar ", "dil_phi_ttbar", 50, -4, 4);
    TH1F* ttbar_lep_pt = new TH1F("lep_pt_ttbar ", "lep_pt_ttbar", 36, 50, 1850);
    ttbar_lep_pt->Sumw2();
    TH1F* ttbar_lep_pt_BB = new TH1F("lep_pt_ttbar_BB ", "lep_pt_ttbar_BB", 36, 50, 1850);
    ttbar_lep_pt_BB->Sumw2();
    TH1F* ttbar_lep_pt_BE = new TH1F("lep_pt_ttbar_BE ", "lep_pt_ttbar_BE", 36, 50, 1850);
    ttbar_lep_pt_BE->Sumw2();
    TH1F* ttbar_lep_pt_plus = new TH1F("lep_pt_ttbar_plus ", "lep_pt_ttbar_plus", 100, 0, 2000);
    TH1F* ttbar_lep_pt_minus = new TH1F("lep_pt_ttbar_minus ", "lep_pt_ttbar_minus", 100, 0, 2000);
    TH1F* ttbar_lep_eta = new TH1F("lep_eta_ttbar ", "lep_eta_ttbar", 50, -3, 3);
    TH1F* ttbar_lep_eta_BB = new TH1F("lep_eta_ttbar_BB ", "lep_eta_ttbar_BB", 50, -3, 3);
    TH1F* ttbar_lep_eta_BE = new TH1F("lep_eta_ttbar_BE ", "lep_eta_ttbar_BE", 50, -3, 3);
    TH1F* ttbar_lep_phi = new TH1F("lep_phi_ttbar ", "lep_phi_ttbar", 50, -4, 4);
    TH1F* dB_ttbar_clear = new TH1F("ttbar_dB", "ttbar_dB",  50, 0, 0.25);
    TH1F* PixelHit_ttbar_clear = new TH1F("ttbar_PixelHit", "ttbar_PixelHit",  15, 0,  15);
    TH1F* TkLayer_ttbar_clear = new TH1F("ttbar_TkLayer", "ttbar_TkLayer", 20, 0, 20);
    TH1F* Iso_ttbar_clear = new TH1F("ttbar_Iso", "ttbar_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_ttbar_clear = new TH1F("ttbar_relpTErr", "ttbar_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_ttbar_clear = new TH1F("ttbar_Vtx", "ttbar_Vtx",  44, 0, 22);
    TH1F* ValidMu_ttbar_clear = new TH1F("ttbar_ValidMu", "ttbar_ValidMu",  55, 0, 55);
    
    
    
    
    
    //////////////////////////////tW/////////////////////////////////////////////
    TH1F* tW_pt = new TH1F("dil_pt_tW ", "dil_pt_tW", 100, 0, 2000);
    TH1F* tW_rapidity = new TH1F("dil_rapidity_tW ", "dil_rapidity_tW", 50, -3, 3);
    TH1F* tW_rapidity_BB = new TH1F("dil_rapidity_tW_BB ", "dil_rapidity_tW_BB", 50, -3, 3);
    TH1F* tW_rapidity_BE = new TH1F("dil_rapidity_tW_BE ", "dil_rapidity_tW_BE", 50, -3, 3);
    TH1F* tW_phi = new TH1F("dil_phi_tW ", "dil_phi_tW", 50, -4, 4);
    TH1F* tW_lep_pt = new TH1F("lep_pt_tW ", "lep_pt_tW", 36, 50, 1850);
    tW_lep_pt->Sumw2();
    TH1F* tW_lep_pt_BB = new TH1F("lep_pt_tW_BB ", "lep_pt_tW_BB", 36, 50, 1850);
    tW_lep_pt_BB->Sumw2();
    TH1F* tW_lep_pt_BE = new TH1F("lep_pt_tW_BE ", "lep_pt_tW_BE", 36, 50, 1850);
    tW_lep_pt_BE->Sumw2();
    TH1F* tW_lep_pt_plus = new TH1F("lep_pt_tW_plus ", "lep_pt_tW_plus", 100, 0, 2000);
    TH1F* tW_lep_pt_minus = new TH1F("lep_pt_tW_minus ", "lep_pt_tW_minus", 100, 0, 2000);
    TH1F* tW_lep_eta = new TH1F("lep_eta_tW ", "lep_eta_tW", 50, -3, 3);
    TH1F* tW_lep_eta_BB = new TH1F("lep_eta_tW_BB ", "lep_eta_tW_BB", 50, -3, 3);
    TH1F* tW_lep_eta_BE = new TH1F("lep_eta_tW_BE ", "lep_eta_tW_BE", 50, -3, 3);
    TH1F* tW_lep_phi = new TH1F("lep_phi_tW ", "lep_phi_tW", 50, -4, 4);
    TH1F* dB_tW_clear = new TH1F("tW_dB", "tW_dB",  50, 0, 0.25);
    TH1F* PixelHit_tW_clear = new TH1F("tW_PixelHit", "tW_PixelHit",  15, 0,  15);
    TH1F* TkLayer_tW_clear = new TH1F("tW_TkLayer", "tW_TkLayer", 20, 0, 20);
    TH1F* Iso_tW_clear = new TH1F("tW_Iso", "tW_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_tW_clear = new TH1F("tW_relpTErr", "tW_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_tW_clear = new TH1F("tW_Vtx", "tW_Vtx",  44, 0, 22);
    TH1F* ValidMu_tW_clear = new TH1F("tW_ValidMu", "tW_ValidMu",  55, 0, 55);
    
    
    
    
    //////////////////////////////Wantitop/////////////////////////////////////////////
    TH1F* Wantitop_pt = new TH1F("dil_pt_Wantitop ", "dil_pt_Wantitop", 100, 0, 2000);
    TH1F* Wantitop_rapidity = new TH1F("dil_rapidity_Wantitop ", "dil_rapidity_Wantitop", 50, -3, 3);
    TH1F* Wantitop_rapidity_BB = new TH1F("dil_rapidity_Wantitop_BB ", "dil_rapidity_Wantitop_BB", 50, -3, 3);
    TH1F* Wantitop_rapidity_BE = new TH1F("dil_rapidity_Wantitop_BE ", "dil_rapidity_Wantitop_BE", 50, -3, 3);
    TH1F* Wantitop_phi = new TH1F("dil_phi_Wantitop ", "dil_phi_Wantitop", 50, -4, 4);
    TH1F* Wantitop_lep_pt = new TH1F("lep_pt_Wantitop ", "lep_pt_Wantitop", 36, 50, 1850);
    Wantitop_lep_pt->Sumw2();
    TH1F* Wantitop_lep_pt_BB = new TH1F("lep_pt_Wantitop_BB ", "lep_pt_Wantitop_BB", 36, 50, 1850);
    Wantitop_lep_pt_BB->Sumw2();
    TH1F* Wantitop_lep_pt_BE = new TH1F("lep_pt_Wantitop_BE ", "lep_pt_Wantitop_BE", 36, 50, 1850);
    Wantitop_lep_pt_BE->Sumw2();
    TH1F* Wantitop_lep_pt_plus = new TH1F("lep_pt_Wantitop_plus ", "lep_pt_Wantitop_plus", 100, 0, 2000);
    TH1F* Wantitop_lep_pt_minus = new TH1F("lep_pt_Wantitop_minus ", "lep_pt_Wantitop_minus", 100, 0, 2000);
    TH1F* Wantitop_lep_eta = new TH1F("lep_eta_Wantitop ", "lep_eta_Wantitop", 50, -3, 3);
    TH1F* Wantitop_lep_eta_BB = new TH1F("lep_eta_Wantitop_BB ", "lep_eta_Wantitop_BB", 50, -3, 3);
    TH1F* Wantitop_lep_eta_BE = new TH1F("lep_eta_Wantitop_BE ", "lep_eta_Wantitop_BE", 50, -3, 3);
    TH1F* Wantitop_lep_phi = new TH1F("lep_phi_Wantitop ", "lep_phi_Wantitop", 50, -4, 4);
    TH1F* dB_Wantitop_clear = new TH1F("Wantitop_dB", "Wantitop_dB",  50, 0, 0.25);
    TH1F* PixelHit_Wantitop_clear = new TH1F("Wantitop_PixelHit", "Wantitop_PixelHit",  15, 0,  15);
    TH1F* TkLayer_Wantitop_clear = new TH1F("Wantitop_TkLayer", "Wantitop_TkLayer", 20, 0, 20);
    TH1F* Iso_Wantitop_clear = new TH1F("Wantitop_Iso", "Wantitop_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_Wantitop_clear = new TH1F("Wantitop_relpTErr", "Wantitop_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_Wantitop_clear = new TH1F("Wantitop_Vtx", "Wantitop_Vtx",  44, 0, 22);
    TH1F* ValidMu_Wantitop_clear = new TH1F("Wantitop_ValidMu", "Wantitop_ValidMu",  55, 0, 55);
    
    
    
    
    //////////////////////////////WW/////////////////////////////////////////////
    TH1F* WW_pt = new TH1F("dil_pt_WW ", "dil_pt_WW", 100, 0, 2000);
    TH1F* WW_rapidity = new TH1F("dil_rapidity_WW ", "dil_rapidity_WW", 50, -3, 3);
    TH1F* WW_rapidity_BB = new TH1F("dil_rapidity_WW_BB ", "dil_rapidity_WW_BB", 50, -3, 3);
    TH1F* WW_rapidity_BE = new TH1F("dil_rapidity_WW_BE ", "dil_rapidity_WW_BE", 50, -3, 3);
    TH1F* WW_phi = new TH1F("dil_phi_WW ", "dil_phi_WW", 50, -4, 4);
    TH1F* WW_lep_pt = new TH1F("lep_pt_WW ", "lep_pt_WW", 36, 50, 1850);
    WW_lep_pt->Sumw2();
    TH1F* WW_lep_pt_BB = new TH1F("lep_pt_WW_BB ", "lep_pt_WW_BB", 36, 50, 1850);
    WW_lep_pt_BB->Sumw2();
    TH1F* WW_lep_pt_BE = new TH1F("lep_pt_WW_BE ", "lep_pt_WW_BE", 36, 50, 1850);
    WW_lep_pt_BE->Sumw2();
    TH1F* WW_lep_pt_plus = new TH1F("lep_pt_WW_plus ", "lep_pt_WW_plus", 100, 0, 2000);
    TH1F* WW_lep_pt_minus = new TH1F("lep_pt_WW_minus ", "lep_pt_WW_minus", 100, 0, 2000);
    TH1F* WW_lep_eta = new TH1F("lep_eta_WW ", "lep_eta_WW", 50, -3, 3);
    TH1F* WW_lep_eta_BB = new TH1F("lep_eta_WW_BB ", "lep_eta_WW_BB", 50, -3, 3);
    TH1F* WW_lep_eta_BE = new TH1F("lep_eta_WW_BE ", "lep_eta_WW_BE", 50, -3, 3);
    TH1F* WW_lep_phi = new TH1F("lep_phi_WW ", "lep_phi_WW", 50, -4, 4);
    TH1F* dB_WW_clear = new TH1F("WW_dB", "WW_dB",  50, 0, 0.25);
    TH1F* PixelHit_WW_clear = new TH1F("WW_PixelHit", "WW_PixelHit",  15, 0,  15);
    TH1F* TkLayer_WW_clear = new TH1F("WW_TkLayer", "WW_TkLayer", 20, 0, 20);
    TH1F* Iso_WW_clear = new TH1F("WW_Iso", "WW_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_WW_clear = new TH1F("WW_relpTErr", "WW_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_WW_clear = new TH1F("WW_Vtx", "WW_Vtx",  44, 0, 22);
    TH1F* ValidMu_WW_clear = new TH1F("WW_ValidMu", "WW_ValidMu",  55, 0, 55);
    
    
    
    
    //////////////////////////////WZ/////////////////////////////////////////////
    TH1F* WZ_pt = new TH1F("dil_pt_WZ ", "dil_pt_WZ", 100, 0, 2000);
    TH1F* WZ_rapidity = new TH1F("dil_rapidity_WZ ", "dil_rapidity_WZ", 50, -3, 3);
    TH1F* WZ_rapidity_BB = new TH1F("dil_rapidity_WZ_BB ", "dil_rapidity_WZ_BB", 50, -3, 3);
    TH1F* WZ_rapidity_BE = new TH1F("dil_rapidity_WZ_BE ", "dil_rapidity_WZ_BE", 50, -3, 3);
    TH1F* WZ_phi = new TH1F("dil_phi_WZ ", "dil_phi_WZ", 50, -4, 4);
    TH1F* WZ_lep_pt = new TH1F("lep_pt_WZ ", "lep_pt_WZ", 36, 50, 1850);
    WZ_lep_pt->Sumw2();
    TH1F* WZ_lep_pt_BB = new TH1F("lep_pt_WZ_BB ", "lep_pt_WZ_BB", 36, 50, 1850);
    WZ_lep_pt_BB->Sumw2();
    TH1F* WZ_lep_pt_BE = new TH1F("lep_pt_WZ_BE ", "lep_pt_WZ_BE", 36, 50, 1850);
    WZ_lep_pt_BE->Sumw2();
    TH1F* WZ_lep_pt_plus = new TH1F("lep_pt_WZ_plus ", "lep_pt_WZ_plus", 100, 0, 2000);
    TH1F* WZ_lep_pt_minus = new TH1F("lep_pt_WZ_minus ", "lep_pt_WZ_minus", 100, 0, 2000);
    TH1F* WZ_lep_eta = new TH1F("lep_eta_WZ ", "lep_eta_WZ", 50, -3, 3);
    TH1F* WZ_lep_eta_BB = new TH1F("lep_eta_WZ_BB ", "lep_eta_WZ_BB", 50, -3, 3);
    TH1F* WZ_lep_eta_BE = new TH1F("lep_eta_WZ_BE ", "lep_eta_WZ_BE", 50, -3, 3);
    TH1F* WZ_lep_phi = new TH1F("lep_phi_WZ ", "lep_phi_WZ", 50, -4, 4);
    TH1F* dB_WZ_clear = new TH1F("WZ_dB", "WZ_dB",  50, 0, 0.25);
    TH1F* PixelHit_WZ_clear = new TH1F("WZ_PixelHit", "WZ_PixelHit",  15, 0,  15);
    TH1F* TkLayer_WZ_clear = new TH1F("WZ_TkLayer", "WZ_TkLayer", 20, 0, 20);
    TH1F* Iso_WZ_clear = new TH1F("WZ_Iso", "WZ_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_WZ_clear = new TH1F("WZ_relpTErr", "WZ_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_WZ_clear = new TH1F("WZ_Vtx", "WZ_Vtx",  44, 0, 22);
    TH1F* ValidMu_WZ_clear = new TH1F("WZ_ValidMu", "WZ_ValidMu",  55, 0, 55);
    
    
    
    
    
    //////////////////////////////ZZ/////////////////////////////////////////////
    TH1F* ZZ_pt = new TH1F("dil_pt_ZZ ", "dil_pt_ZZ", 100, 0, 2000);
    TH1F* ZZ_rapidity = new TH1F("dil_rapidity_ZZ ", "dil_rapidity_ZZ", 50, -3, 3);
    TH1F* ZZ_rapidity_BB = new TH1F("dil_rapidity_ZZ_BB ", "dil_rapidity_ZZ_BB", 50, -3, 3);
    TH1F* ZZ_rapidity_BE = new TH1F("dil_rapidity_ZZ_BE ", "dil_rapidity_ZZ_BE", 50, -3, 3);
    TH1F* ZZ_phi = new TH1F("dil_phi_ZZ ", "dil_phi_ZZ", 50, -4, 4);
    TH1F* ZZ_lep_pt = new TH1F("lep_pt_ZZ ", "lep_pt_ZZ", 36, 50, 1850);
    ZZ_lep_pt->Sumw2();
    TH1F* ZZ_lep_pt_BB = new TH1F("lep_pt_ZZ_BB ", "lep_pt_ZZ_BB", 36, 50, 1850);
    ZZ_lep_pt_BB->Sumw2();
    TH1F* ZZ_lep_pt_BE = new TH1F("lep_pt_ZZ_BE ", "lep_pt_ZZ_BE", 36, 50, 1850);
    ZZ_lep_pt_BE->Sumw2();
    TH1F* ZZ_lep_pt_plus = new TH1F("lep_pt_ZZ_plus ", "lep_pt_ZZ_plus", 100, 0, 2000);
    TH1F* ZZ_lep_pt_minus = new TH1F("lep_pt_ZZ_minus ", "lep_pt_ZZ_minus", 100, 0, 2000);
    TH1F* ZZ_lep_eta = new TH1F("lep_eta_ZZ ", "lep_eta_ZZ", 50, -3, 3);
    TH1F* ZZ_lep_eta_BB = new TH1F("lep_eta_ZZ_BB ", "lep_eta_ZZ_BB", 50, -3, 3);
    TH1F* ZZ_lep_eta_BE = new TH1F("lep_eta_ZZ_BE ", "lep_eta_ZZ_BE", 50, -3, 3);
    TH1F* ZZ_lep_phi = new TH1F("lep_phi_ZZ ", "lep_phi_ZZ", 50, -4, 4);
    TH1F* dB_ZZ_clear = new TH1F("ZZ_dB", "ZZ_dB",  50, 0, 0.25);
    TH1F* PixelHit_ZZ_clear = new TH1F("ZZ_PixelHit", "ZZ_PixelHit",  15, 0,  15);
    TH1F* TkLayer_ZZ_clear = new TH1F("ZZ_TkLayer", "ZZ_TkLayer", 20, 0, 20);
    TH1F* Iso_ZZ_clear = new TH1F("ZZ_Iso", "ZZ_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_ZZ_clear = new TH1F("ZZ_relpTErr", "ZZ_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_ZZ_clear = new TH1F("ZZ_Vtx", "ZZ_Vtx",  44, 0, 22);
    TH1F* ValidMu_ZZ_clear = new TH1F("ZZ_ValidMu", "ZZ_ValidMu",  55, 0, 55);
    
    
    
    
    //////////////////////////////Wjets/////////////////////////////////////////////
    TH1F* Wjets_pt = new TH1F("dil_pt_Wjets ", "dil_pt_Wjets", 100, 0, 2000);
    TH1F* Wjets_rapidity = new TH1F("dil_rapidity_Wjets ", "dil_rapidity_Wjets", 50, -3, 3);
    TH1F* Wjets_rapidity_BB = new TH1F("dil_rapidity_Wjets_BB ", "dil_rapidity_Wjets_BB", 50, -3, 3);
    TH1F* Wjets_rapidity_BE = new TH1F("dil_rapidity_Wjets_BE ", "dil_rapidity_Wjets_BE", 50, -3, 3);
    TH1F* Wjets_phi = new TH1F("dil_phi_Wjets ", "dil_phi_Wjets", 50, -4, 4);
    TH1F* Wjets_lep_pt = new TH1F("lep_pt_Wjets ", "lep_pt_Wjets", 36, 50, 1850);
    Wjets_lep_pt->Sumw2();
    TH1F* Wjets_lep_pt_BB = new TH1F("lep_pt_Wjets_BB ", "lep_pt_Wjets_BB", 36, 50, 1850);
    Wjets_lep_pt_BB->Sumw2();
    TH1F* Wjets_lep_pt_BE = new TH1F("lep_pt_Wjets_BE ", "lep_pt_Wjets_BE", 36, 50, 1850);
    Wjets_lep_pt_BE->Sumw2();
    TH1F* Wjets_lep_pt_plus = new TH1F("lep_pt_Wjets_plus ", "lep_pt_Wjets_plus", 100, 0, 2000);
    TH1F* Wjets_lep_pt_minus = new TH1F("lep_pt_Wjets_minus ", "lep_pt_Wjets_minus", 100, 0, 2000);
    TH1F* Wjets_lep_eta = new TH1F("lep_eta_Wjets ", "lep_eta_Wjets", 50, -3, 3);
    TH1F* Wjets_lep_eta_BB = new TH1F("lep_eta_Wjets_BB ", "lep_eta_Wjets_BB", 50, -3, 3);
    TH1F* Wjets_lep_eta_BE = new TH1F("lep_eta_Wjets_BE ", "lep_eta_Wjets_BE", 50, -3, 3);
    TH1F* Wjets_lep_phi = new TH1F("lep_phi_Wjets ", "lep_phi_Wjets", 50, -4, 4);
    TH1F* dB_Wjets_clear = new TH1F("Wjets_dB", "Wjets_dB",  50, 0, 0.25);
    TH1F* PixelHit_Wjets_clear = new TH1F("Wjets_PixelHit", "Wjets_PixelHit",  15, 0,  15);
    TH1F* TkLayer_Wjets_clear = new TH1F("Wjets_TkLayer", "Wjets_TkLayer", 20, 0, 20);
    TH1F* Iso_Wjets_clear = new TH1F("Wjets_Iso", "Wjets_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_Wjets_clear = new TH1F("Wjets_relpTErr", "Wjets_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_Wjets_clear = new TH1F("Wjets_Vtx", "Wjets_Vtx",  44, 0, 22);
    TH1F* ValidMu_Wjets_clear = new TH1F("Wjets_ValidMu", "Wjets_ValidMu",  55, 0, 55);
    
    
    //////////////////////////////Zjets/////////////////////////////////////////////
    TH1F* Zjets_pt = new TH1F("dil_pt_Zjets ", "dil_pt_Zjets", 100, 0, 2000);
    TH1F* Zjets_rapidity = new TH1F("dil_rapidity_Zjets ", "dil_rapidity_Zjets", 50, -3, 3);
    TH1F* Zjets_rapidity_BB = new TH1F("dil_rapidity_Zjets_BB ", "dil_rapidity_Zjets_BB", 50, -3, 3);
    TH1F* Zjets_rapidity_BE = new TH1F("dil_rapidity_Zjets_BE ", "dil_rapidity_Zjets_BE", 50, -3, 3);
    TH1F* Zjets_phi = new TH1F("dil_phi_Zjets ", "dil_phi_Zjets", 50, -4, 4);
    TH1F* Zjets_lep_pt = new TH1F("lep_pt_Zjets ", "lep_pt_Zjets", 36, 50, 1850);
    Zjets_lep_pt->Sumw2();
    TH1F* Zjets_lep_pt_BB = new TH1F("lep_pt_Zjets_BB ", "lep_pt_Zjets_BB", 36, 50, 1850);
    Zjets_lep_pt_BB->Sumw2();
    TH1F* Zjets_lep_pt_BE = new TH1F("lep_pt_Zjets_BE ", "lep_pt_Zjets_BE", 36, 50, 1850);
    Zjets_lep_pt_BE->Sumw2();
    TH1F* Zjets_lep_pt_plus = new TH1F("lep_pt_Zjets_plus ", "lep_pt_Zjets_plus", 100, 0, 2000);
    TH1F* Zjets_lep_pt_minus = new TH1F("lep_pt_Zjets_minus ", "lep_pt_Zjets_minus", 100, 0, 2000);
    TH1F* Zjets_lep_eta = new TH1F("lep_eta_Zjets ", "lep_eta_Zjets", 50, -3, 3);
    TH1F* Zjets_lep_eta_BB = new TH1F("lep_eta_Zjets_BB ", "lep_eta_Zjets_BB", 50, -3, 3);
    TH1F* Zjets_lep_eta_BE = new TH1F("lep_eta_Zjets_BE ", "lep_eta_Zjets_BE", 50, -3, 3);
    TH1F* Zjets_lep_phi = new TH1F("lep_phi_Zjets ", "lep_phi_Zjets", 50, -4, 4);
    TH1F* dB_Zjets_clear = new TH1F("Zjets_dB", "Zjets_dB",  50, 0, 0.25);
    TH1F* PixelHit_Zjets_clear = new TH1F("Zjets_PixelHit", "Zjets_PixelHit",  15, 0,  15);
    TH1F* TkLayer_Zjets_clear = new TH1F("Zjets_TkLayer", "Zjets_TkLayer", 20, 0, 20);
    TH1F* Iso_Zjets_clear = new TH1F("Zjets_Iso", "Zjets_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_Zjets_clear = new TH1F("Zjets_relpTErr", "Zjets_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_Zjets_clear = new TH1F("Zjets_Vtx", "Zjets_Vtx",  44, 0, 22);
    TH1F* ValidMu_Zjets_clear = new TH1F("Zjets_ValidMu", "Zjets_ValidMu",  55, 0, 55);
    
    
    
    
    //////////////////////////////Zprime/////////////////////////////////////////////
    TH1F* Zprime_pt = new TH1F("dil_pt_Zprime ", "dil_pt_Zprime", 100, 0, 2000);
    TH1F* Zprime_rapidity = new TH1F("dil_rapidity_Zprime ", "dil_rapidity_Zprime", 50, -3, 3);
    TH1F* Zprime_rapidity_BB = new TH1F("dil_rapidity_Zprime_BB ", "dil_rapidity_Zprime_BB", 50, -3, 3);
    TH1F* Zprime_rapidity_BE = new TH1F("dil_rapidity_Zprime_BE ", "dil_rapidity_Zprime_BE", 50, -3, 3);
    TH1F* Zprime_phi = new TH1F("dil_phi_Zprime ", "dil_phi_Zprime", 50, -4, 4);
    TH1F* Zprime_lep_pt = new TH1F("lep_pt_Zprime ", "lep_pt_Zprime", 36, 50, 1850);
    Zprime_lep_pt->Sumw2();
    TH1F* Zprime_lep_pt_BB = new TH1F("lep_pt_Zprime_BB ", "lep_pt_Zprime_BB", 36, 50, 1850);
    Zprime_lep_pt_BB->Sumw2();
    TH1F* Zprime_lep_pt_BE = new TH1F("lep_pt_Zprime_BE ", "lep_pt_Zprime_BE", 36, 50, 1850);
    Zprime_lep_pt_BE->Sumw2();
    TH1F* Zprime_lep_pt_plus = new TH1F("lep_pt_Zprime_plus ", "lep_pt_Zprime_plus", 100, 0, 2000);
    TH1F* Zprime_lep_pt_minus = new TH1F("lep_pt_Zprime_minus ", "lep_pt_Zprime_minus", 100, 0, 2000);
    TH1F* Zprime_lep_eta = new TH1F("lep_eta_Zprime ", "lep_eta_Zprime", 50, -3, 3);
    TH1F* Zprime_lep_eta_BB = new TH1F("lep_eta_Zprime_BB ", "lep_eta_Zprime_BB", 50, -3, 3);
    TH1F* Zprime_lep_eta_BE = new TH1F("lep_eta_Zprime_BE ", "lep_eta_Zprime_BE", 50, -3, 3);
    TH1F* Zprime_lep_phi = new TH1F("lep_phi_Zprime ", "lep_phi_Zprime", 50, -4, 4);
    TH1F* dB_Zprime_clear = new TH1F("Zprime_dB", "Zprime_dB",  50, 0, 0.25);
    TH1F* PixelHit_Zprime_clear = new TH1F("Zprime_PixelHit", "Zprime_PixelHit",  15, 0,  15);
    TH1F* TkLayer_Zprime_clear = new TH1F("Zprime_TkLayer", "Zprime_TkLayer", 20, 0, 20);
    TH1F* Iso_Zprime_clear = new TH1F("Zprime_Iso", "Zprime_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_Zprime_clear = new TH1F("Zprime_relpTErr", "Zprime_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_Zprime_clear = new TH1F("Zprime_Vtx", "Zprime_Vtx",  44, 0, 22);
    TH1F* ValidMu_Zprime_clear = new TH1F("Zprime_ValidMu", "Zprime_ValidMu",  55, 0, 55);
    
    
    ////////////////////////////////////Data////////////////////////////////////////////
    
    TH1F* DATA_pt = new TH1F("dil_pt_DATA ", "dil_pt_DATA", 100, 0, 2000);
    TH1F* DATA_rapidity = new TH1F("dil_rap_DATA ", "dil_rap_DATA", 50, -3, 3);
    TH1F* DATA_rapidity_BB = new TH1F("dil_rap_DATA_BB ", "dil_rap_DATA_BB", 50, -3, 3);
    TH1F* DATA_rapidity_BE = new TH1F("dil_rap_DATA_BE ", "dil_rap_DATA_BE", 50, -3, 3);
    TH1F* DATA_phi = new TH1F("dil_phi_DATA ", "dil_phi_DATA", 50, -4, 4);
    TH1F* DATA_lep_pt = new TH1F("lep_pt_DATA ", "lep_pt_DATA", 36, 50, 1850);
    DATA_lep_pt->Sumw2();
    TH1F* DATA_lep_pt_BB = new TH1F("lep_pt_DATA_BB ", "lep_pt_DATA_BB", 36, 50, 1850);
    DATA_lep_pt_BB->Sumw2();
    TH1F* DATA_lep_pt_BE = new TH1F("lep_pt_DATA_BE ", "lep_pt_DATA_BE", 36, 50, 1850);
    DATA_lep_pt_BE->Sumw2();
    TH1F* DATA_lep_pt_plus = new TH1F("lep_pt_DATA_plus ", "lep_pt_DATA_plus", 100, 0, 2000);
    TH1F* DATA_lep_pt_minus = new TH1F("lep_pt_DATA_minus ", "lep_pt_DATA_minus", 100, 0, 2000);
    TH1F* DATA_lep_eta = new TH1F("lep_eta_DATA ", "lep_eta_DATA", 50, -3, 3);
    TH1F* DATA_lep_eta_BB = new TH1F("lep_eta_DATA_BB ", "lep_eta_DATA_BB", 50, -3, 3);
    TH1F* DATA_lep_eta_BE = new TH1F("lep_eta_DATA_BE ", "lep_eta_DATA_BE", 50, -3, 3);
    TH1F* DATA_lep_phi = new TH1F("lep_phi_DATA ", "lep_phi_DATA", 50, -4, 4);
    TH1F* dB_DATA_clear = new TH1F("DATA_dB", "DATA_dB",  50, 0, 0.25);
    TH1F* PixelHit_DATA_clear = new TH1F("DATA_PixelHit", "DATA_PixelHit",  15, 0,  15);
    TH1F* TkLayer_DATA_clear = new TH1F("DATA_TkLayer", "DATA_TkLayer", 20, 0, 20);
    TH1F* Iso_DATA_clear = new TH1F("DATA_Iso", "DATA_Iso", 60, 0.0, 0.1);
    TH1F* relpTErr_DATA_clear = new TH1F("DATA_relpTErr", "DATA_relpTErr", 50, 10e-3, 0.3);
    TH1F* Vtx_DATA_clear = new TH1F("DATA_Vtx", "DATA_Vtx",  44, 0, 22);
    TH1F* ValidMu_DATA_clear = new TH1F("DATA_ValidMu", "DATA_ValidMu",  55, 0, 55);
    
    
    
    
    //strat of looping over MC samples
    for(int j=0; j<mc; j++){
        
        std::cout<<"opening.. "<<MC_samples[j]<<" --- "<<" No of evnets: "<<events[j]<<" Cross section: "<<sigma[j]<<std::endl;
        
        TChain *treeMC = new TChain("SimpleNtupler/t");
        
        treeMC->Add("/eos/user/k/kaliyana/2017_MC/FileBased/"+ MC_samples[j] +"/"+ MC_samples[j] +".root");
        
        Long64_t ne = treeMC->GetEntries();
        
        treeMC->SetBranchAddress("genWeight",&genWeight);
        treeMC->SetBranchAddress("event",&event);
        treeMC->SetBranchAddress("run",&run);
        treeMC->SetBranchAddress("lumi",&lumi);
        treeMC->SetBranchAddress("dil_mass",&dil_mass);
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
            //if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            treeMC->GetEntry(p);
            
            /*    if(j==0 && run==1 && lumi==1318 && event==1317117){
             
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[0] << std::setw(15) << lep_pt[0]  << std::setw(15) << lep_isTrackerMuon[0] << std::setw(15) << lep_isGlobalMuon[0] << std::setw(15) << fabs(lep_dB[0]) << std::setw(15) << lep_glb_numberOfValidPixelHits[0] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[0] << std::setw(15)  << lep_sumPt[0]/lep_tk_pt[0] << std::setw(15) << lep_pt_err[0]/lep_pt[0] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[0] << std::setw(15) << lep_triggerMatchPt[0] << std::setw(15)  << lep_id[0] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[1] << std::setw(15) << lep_pt[1]  << std::setw(15) << lep_isTrackerMuon[1] << std::setw(15) << lep_isGlobalMuon[1] << std::setw(15) << fabs(lep_dB[1]) << std::setw(15) << lep_glb_numberOfValidPixelHits[1] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[1] << std::setw(15)  << lep_sumPt[1]/lep_tk_pt[1] << std::setw(15) << lep_pt_err[1]/lep_pt[1] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[1] << std::setw(15) << lep_triggerMatchPt[1] << std::setw(15)  << lep_id[1] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             }
             
             if(j==2 && run==1 && lumi==40 && event==39294){
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[0] << std::setw(15) << lep_pt[0]  << std::setw(15) << lep_isTrackerMuon[0] << std::setw(15) << lep_isGlobalMuon[0] << std::setw(15) << fabs(lep_dB[0]) << std::setw(15) << lep_glb_numberOfValidPixelHits[0] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[0] << std::setw(15)  << lep_sumPt[0]/lep_tk_pt[0] << std::setw(15) << lep_pt_err[0]/lep_pt[0] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[0] << std::setw(15) << lep_triggerMatchPt[0] << std::setw(15)  << lep_id[0] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[1] << std::setw(15) << lep_pt[1]  << std::setw(15) << lep_isTrackerMuon[1] << std::setw(15) << lep_isGlobalMuon[1] << std::setw(15) << fabs(lep_dB[1]) << std::setw(15) << lep_glb_numberOfValidPixelHits[1] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[1] << std::setw(15)  << lep_sumPt[1]/lep_tk_pt[1] << std::setw(15) << lep_pt_err[1]/lep_pt[1] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[1] << std::setw(15) << lep_triggerMatchPt[1] << std::setw(15)  << lep_id[1] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             }
             
             if(j==0 && run==1 && lumi==838 && event==837137){
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[0] << std::setw(15) << lep_pt[0]  << std::setw(15) << lep_isTrackerMuon[0] << std::setw(15) << lep_isGlobalMuon[0] << std::setw(15) << fabs(lep_dB[0]) << std::setw(15) << lep_glb_numberOfValidPixelHits[0] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[0] << std::setw(15)  << lep_sumPt[0]/lep_tk_pt[0] << std::setw(15) << lep_pt_err[0]/lep_pt[0] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[0] << std::setw(15) << lep_triggerMatchPt[0] << std::setw(15)  << lep_id[0] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[1] << std::setw(15) << lep_pt[1]  << std::setw(15) << lep_isTrackerMuon[1] << std::setw(15) << lep_isGlobalMuon[1] << std::setw(15) << fabs(lep_dB[1]) << std::setw(15) << lep_glb_numberOfValidPixelHits[1] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[1] << std::setw(15)  << lep_sumPt[1]/lep_tk_pt[1] << std::setw(15) << lep_pt_err[1]/lep_pt[1] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[1] << std::setw(15) << lep_triggerMatchPt[1] << std::setw(15)  << lep_id[1] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             }
             
             if(j==1 && run==1 && lumi==64 && event==63967){
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[0] << std::setw(15) << lep_pt[0]  << std::setw(15) << lep_isTrackerMuon[0] << std::setw(15) << lep_isGlobalMuon[0] << std::setw(15) << fabs(lep_dB[0]) << std::setw(15) << lep_glb_numberOfValidPixelHits[0] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[0] << std::setw(15)  << lep_sumPt[0]/lep_tk_pt[0] << std::setw(15) << lep_pt_err[0]/lep_pt[0] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[0] << std::setw(15) << lep_triggerMatchPt[0] << std::setw(15)  << lep_id[0] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[1] << std::setw(15) << lep_pt[1]  << std::setw(15) << lep_isTrackerMuon[1] << std::setw(15) << lep_isGlobalMuon[1] << std::setw(15) << fabs(lep_dB[1]) << std::setw(15) << lep_glb_numberOfValidPixelHits[1] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[1] << std::setw(15)  << lep_sumPt[1]/lep_tk_pt[1] << std::setw(15) << lep_pt_err[1]/lep_pt[1] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[1] << std::setw(15) << lep_triggerMatchPt[1] << std::setw(15)  << lep_id[1] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             }
             
             if(j==2 && run==1 && lumi==40 && event==39651){
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[0] << std::setw(15) << lep_pt[0]  << std::setw(15) << lep_isTrackerMuon[0] << std::setw(15) << lep_isGlobalMuon[0] << std::setw(15) << fabs(lep_dB[0]) << std::setw(15) << lep_glb_numberOfValidPixelHits[0] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[0] << std::setw(15)  << lep_sumPt[0]/lep_tk_pt[0] << std::setw(15) << lep_pt_err[0]/lep_pt[0] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[0] << std::setw(15) << lep_triggerMatchPt[0] << std::setw(15)  << lep_id[0] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             
             std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[1] << std::setw(15) << lep_pt[1]  << std::setw(15) << lep_isTrackerMuon[1] << std::setw(15) << lep_isGlobalMuon[1] << std::setw(15) << fabs(lep_dB[1]) << std::setw(15) << lep_glb_numberOfValidPixelHits[1] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[1] << std::setw(15)  << lep_sumPt[1]/lep_tk_pt[1] << std::setw(15) << lep_pt_err[1]/lep_pt[1] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[1] << std::setw(15) << lep_triggerMatchPt[1] << std::setw(15)  << lep_id[1] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
             }
             
             */
            //start selection on entries
            if(
               //GoodVtx &&
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
               ){
                
                
                /*   if(j==4 && event==67760){
                 std::cout<< "event: "<<event<< " lumi: "<<lumi<< " run: "<<run<<std::endl;
                 std::cout<< "dil_pt: "<<dil_pt<<std::endl;
                 std::cout<< "vertex_m: "<<vertex_m<<std::endl;
                 std::cout<< "lep_pt[0]: "<<lep_pt[0]<<std::endl;
                 std::cout<< "lep_pt[1]: "<<lep_pt[1]<<std::endl;
                 std::cout<< "--------------------------------------"<<std::endl;
                 }
                 
                 */
                //std::cout<<dil_mass<<std::endl;
                if(prev_event == event) {
                    //std::cout<< "event: "<<event<<std::endl;
                    //std::cout<< "prev_event: "<<prev_event<<std::endl;
                    continue; //to remove double event; take the first --> they are already sorted in pt
                }
                
                prev_event = event;
                
                
                if(dil_mass < 0 || gen_dil_mass < 0){
                    std::cout<<MC_samples[j]<<"   "<<run<<"   "<<lumi<<"  "<<event<<" "<<dil_mass<<" "<<gen_dil_mass<<std::endl;
                    //continue;
                }
              
                
                weight[j] *= genWeight;
                weight_BB[j] *= genWeight;
                weight_BE[j] *= genWeight;
                
                
                
                
                if(j==10){ //for ttbar
                    
                    gM = gen_dil_mass;
                    //kFactor =0.994078695151 + gM*2.64819793287e-05 - gM*gM*3.73996461024e-08 - gM*gM*gM*1.11452866827e-11;
                    kFactor = 0.991403 + (3.05593e-05 * gM) - (2.21967e-07 * pow(gM,2)) + (6.63658e-11 * pow(gM,3));
                    kFactor_BB = 0.990973 + (6.17124e-06 * gM) - (3.31244e-07 * pow(gM,2)) + (1.2125e-10 * pow(gM,3));
                    kFactor_BE = 0.990038 + (5.44269e-05 * gM) - (2.43311e-07 * pow(gM,2)) + (5.9748e-11 * pow(gM,3));
                    
                    
                    weight[j] = weight[j] * (float) kFactor;
                    weight_BB[j] = weight_BB[j] * (float) kFactor_BB;
                    weight_BE[j] = weight_BE[j] * (float) kFactor_BE;
                    
                    ttbar_pt->Fill(dil_pt, weight[j]);
                    ttbar_rapidity->Fill(dil_rap, weight[j]);
                    ttbar_phi->Fill(dil_phi, weight[j]);
                    Vtx_ttbar_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        ttbar_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        ttbar_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        ttbar_lep_pt->Fill(lep_pt[h],  weight[j]);
                        ttbar_lep_eta->Fill(lep_eta[h],  weight[j]);
                        ttbar_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_ttbar_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_ttbar_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_ttbar_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_ttbar_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_ttbar_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_ttbar_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            ttbar_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            ttbar_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            ttbar_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            ttbar_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            ttbar_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            ttbar_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    weight[j] /= (float) kFactor;
                    weight_BB[j] /= (float) kFactor_BB;
                    weight_BE[j] /= (float) kFactor_BE;
                    
                }// end of ttbar
                
                
                
                
                else if(j==11){ //for tW
                    
                    tW_pt->Fill(dil_pt, weight[j]);
                    tW_rapidity->Fill(dil_rap, weight[j]);
                    tW_phi->Fill(dil_phi, weight[j]);
                    Vtx_tW_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        tW_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        tW_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    
                    for(int h = 0; h < 2; h++){
                        
                        tW_lep_pt->Fill(lep_pt[h],  weight[j]);
                        tW_lep_eta->Fill(lep_eta[h],  weight[j]);
                        tW_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_tW_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_tW_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_tW_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_tW_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_tW_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_tW_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            tW_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            tW_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            tW_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            tW_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            tW_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            tW_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of tW
                
                
                else if(j==12){ //for Wantitop
                    
                    Wantitop_pt->Fill(dil_pt, weight[j]);
                    Wantitop_rapidity->Fill(dil_rap, weight[j]);
                    Wantitop_phi->Fill(dil_phi, weight[j]);
                    Vtx_Wantitop_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Wantitop_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        Wantitop_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    
                    for(int h = 0; h < 2; h++){
                        
                        Wantitop_lep_pt->Fill(lep_pt[h],  weight[j]);
                        Wantitop_lep_eta->Fill(lep_eta[h],  weight[j]);
                        Wantitop_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_Wantitop_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_Wantitop_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_Wantitop_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_Wantitop_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_Wantitop_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_Wantitop_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            Wantitop_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            Wantitop_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            Wantitop_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            Wantitop_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            Wantitop_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            Wantitop_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of Wantitop
                
                
                
                
                else if(j==13){ //for WW
                    
                    WW_pt->Fill(dil_pt, weight[j]);
                    WW_rapidity->Fill(dil_rap, weight[j]);
                    WW_phi->Fill(dil_phi, weight[j]);
                    Vtx_WW_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        WW_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        WW_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        WW_lep_pt->Fill(lep_pt[h],  weight[j]);
                        WW_lep_eta->Fill(lep_eta[h],  weight[j]);
                        WW_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_WW_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_WW_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_WW_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_WW_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_WW_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_WW_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            WW_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            WW_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            WW_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            WW_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            WW_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            WW_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of WW
                
                
                
                else if(j==14){ //for WZ
                    
                    WZ_pt->Fill(dil_pt, weight[j]);
                    WZ_rapidity->Fill(dil_rap, weight[j]);
                    WZ_phi->Fill(dil_phi, weight[j]);
                    Vtx_WZ_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        WZ_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        WZ_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        WZ_lep_pt->Fill(lep_pt[h],  weight[j]);
                        WZ_lep_eta->Fill(lep_eta[h],  weight[j]);
                        WZ_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_WZ_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_WZ_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_WZ_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_WZ_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_WZ_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_WZ_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            WZ_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            WZ_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            WZ_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            WZ_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            WZ_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            WZ_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of WZ
                
                
                else if(j==15){ //for ZZ
                    
                    ZZ_pt->Fill(dil_pt, weight[j]);
                    ZZ_rapidity->Fill(dil_rap, weight[j]);
                    ZZ_phi->Fill(dil_phi, weight[j]);
                    Vtx_ZZ_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        ZZ_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        ZZ_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        ZZ_lep_pt->Fill(lep_pt[h],  weight[j]);
                        ZZ_lep_eta->Fill(lep_eta[h],  weight[j]);
                        ZZ_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_ZZ_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_ZZ_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_ZZ_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_ZZ_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_ZZ_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_ZZ_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            ZZ_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            ZZ_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            ZZ_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            ZZ_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            ZZ_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            ZZ_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of ZZ
                
                
                else if(j==16){ //for Wjets
                    
                    Wjets_pt->Fill(dil_pt, weight[j]);
                    Wjets_rapidity->Fill(dil_rap, weight[j]);
                    Wjets_phi->Fill(dil_phi, weight[j]);
                    Vtx_Wjets_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Wjets_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        Wjets_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        Wjets_lep_pt->Fill(lep_pt[h],  weight[j]);
                        Wjets_lep_eta->Fill(lep_eta[h],  weight[j]);
                        Wjets_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_Wjets_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_Wjets_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_Wjets_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_Wjets_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_Wjets_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_Wjets_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            Wjets_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            Wjets_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            Wjets_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            Wjets_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            Wjets_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            Wjets_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of Wjets
                
                
                
                else if(j==17){ //for Zjets
                    
                    Zjets_pt->Fill(dil_pt, weight[j]);
                    Zjets_rapidity->Fill(dil_rap, weight[j]);
                    Zjets_phi->Fill(dil_phi, weight[j]);
                    Vtx_Zjets_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Zjets_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        Zjets_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        Zjets_lep_pt->Fill(lep_pt[h],  weight[j]);
                        Zjets_lep_eta->Fill(lep_eta[h],  weight[j]);
                        Zjets_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_Zjets_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_Zjets_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_Zjets_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_Zjets_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_Zjets_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_Zjets_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            Zjets_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            Zjets_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            Zjets_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            Zjets_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            Zjets_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            Zjets_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                        
                        
                    } // for on muons
                    
                    
                    
                    
                    
                    
                }// end of Zjets
                
                
                
                
                
                
                else if(j==18){ //for Zprime
                    
                    Zprime_pt->Fill(dil_pt, weight[j]);
                    Zprime_rapidity->Fill(dil_rap, weight[j]);
                    Zprime_phi->Fill(dil_phi, weight[j]);
                    Vtx_Zprime_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        Zprime_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        Zprime_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        Zprime_lep_pt->Fill(lep_pt[h],  weight[j]);
                        Zprime_lep_eta->Fill(lep_eta[h],  weight[j]);
                        Zprime_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_Zprime_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_Zprime_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_Zprime_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_Zprime_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_Zprime_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_Zprime_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            Zprime_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            Zprime_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            Zprime_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            Zprime_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            Zprime_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            Zprime_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                    } // for on muons
                    
                    
                }// end of Zprime
                
                
                else{ //for DY as given in AN_2018_0
                    //std::cout<<j<<"\t"<<MC_samples[j]<<"\t"<<dil_mass<<std::endl;
                    
                    //To reject odd events
                    if(run==1 && (lumi==1318 || lumi==40 || lumi==838 || lumi==64) && (event==1317117 || event==39294 || event==837137 || event==63967 || event==39651)) continue;
                    
                    
                    gM = gen_dil_mass - 400.0;
                    kFactor = 1.067 - (0.000112 * gM) + (3.176e-08 * pow(gM,2)) - (4.068e-12 * pow(gM,3));
                    kFactor_BB = 1.036 - (0.0001441 * gM) + (5.058e-8 * pow(gM,2)) - (7.581e-12 * pow(gM,3));
                    kFactor_BE = 1.052 - (0.0001471 * gM) + (5.903e-8 * pow(gM,2)) - (9.037e-12 * pow(gM,3));
                    
                    weight[j] = weight[j] * (float) kFactor;
                    weight_BB[j] = weight_BB[j] * (float) kFactor_BB;
                    weight_BE[j] = weight_BE[j] * (float) kFactor_BE;
                    
                    DY_pt->Fill(dil_pt, weight[j]);
                    DY_rapidity->Fill(dil_rap, weight[j]);
                    DY_phi->Fill(dil_phi, weight[j]);
                    Vtx_DY_clear->Fill(vertex_chi2,  weight[j]);
                    
                    
                    if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                        DY_rapidity_BB->Fill(dil_rap,  weight_BB[j]);
                    }
                    else{
                        DY_rapidity_BE->Fill(dil_rap,  weight_BE[j]);
                    }
                    
                    for(int h = 0; h < 2; h++){
                        
                        DY_lep_pt->Fill(lep_pt[h],  weight[j]);
                        DY_lep_eta->Fill(lep_eta[h],  weight[j]);
                        DY_lep_phi->Fill(lep_phi[h],  weight[j]);
                        dB_DY_clear->Fill(lep_dB[h],  weight[j]);
                        PixelHit_DY_clear->Fill(lep_glb_numberOfValidPixelHits[h],  weight[j]);
                        TkLayer_DY_clear->Fill(lep_glb_numberOfValidTrackerLayers[h],  weight[j]);
                        Iso_DY_clear->Fill(lep_sumPt[h]/lep_tk_pt[h],  weight[j]);
                        relpTErr_DY_clear->Fill(lep_pt_err[h]/lep_pt[h],  weight[j]);
                        ValidMu_DY_clear->Fill(lep_glb_numberOfValidMuonHits[h],  weight[j]);
                        
                        if(fabs(lep_eta[h]) < 1.2){
                            DY_lep_pt_BB->Fill(lep_pt[h],  weight_BB[j]);
                            DY_lep_eta_BB->Fill(lep_eta[h],  weight_BB[j]);
                        }
                        else{
                            DY_lep_pt_BE->Fill(lep_pt[h],  weight_BE[j]);
                            DY_lep_eta_BE->Fill(lep_eta[h],  weight_BE[j]);
                        }
                        
                        if(gen_lep_qOverPt[h] > 0)
                            DY_lep_pt_plus->Fill(lep_pt[h],  weight[j]);
                        else
                            DY_lep_pt_minus->Fill(lep_pt[h],  weight[j]);
                        
                        
                        
                        ////////////////for_debugging//////////////
                        
                        //std::cout<<"Weight"<<";"<<"eta_muon"<<";"<<"phi_muon"<<";"<<"pt[GeV]"<<";"<<"Run"<<";"<<"Lumi"<<";"<<"Event"<<std::endl;
                        
                        //  if(lep_pt[h] > 700.0 && lep_pt[h] < 850.0){
                        
                        //  std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[h] << std::setw(15) << lep_pt[h]  << std::setw(15) << lep_isTrackerMuon[h] << std::setw(15) << lep_isGlobalMuon[h] << std::setw(15) << fabs(lep_dB[h]) << std::setw(15) << lep_glb_numberOfValidPixelHits[h] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[h] << std::setw(15)  << lep_sumPt[h]/lep_tk_pt[h] << std::setw(15) << lep_pt_err[h]/lep_pt[h] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[h] << std::setw(15) << lep_triggerMatchPt[h] << std::setw(15)  << lep_id[h] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
                        
                        //std::cout<< std::setw(15) << j << std::setw(15) << genWeight << std::setw(15) << weight_BB[j] << std::setw(15) << lep_eta[1] << std::setw(15) << lep_pt[1]  << std::setw(15) << lep_isTrackerMuon[1] << std::setw(15) << lep_isGlobalMuon[1] << std::setw(15) << fabs(lep_dB[1]) << std::setw(15) << lep_glb_numberOfValidPixelHits[1] <<  std::setw(15) << lep_glb_numberOfValidTrackerLayers[1] << std::setw(15)  << lep_sumPt[1]/lep_tk_pt[1] << std::setw(15) << lep_pt_err[1]/lep_pt[1] << std::setw(15) << cos_angle << std::setw(15) << vertex_chi2 << std::setw(15) << dil_mass << std::setw(15)  << gen_dil_mass << std::setw(15) << lep_numberOfMatchedStations[1] << std::setw(15) << lep_triggerMatchPt[1] << std::setw(15)  << lep_id[1] << std::setw(15)  << run << std::setw(15) << lumi << std::setw(15) << event << std::endl;
                        
                        // }
                        
                        /*792.721 96.0208 141.08  877.536
                         1        9.80134       -1.27918        792.721              1              1     0.00092082              4             12    0.000500506      0.0566727       0.966236       0.204573              1             45          40014
                         1        9.80134       -1.31698        96.0208              1              1    0.000410422              4             11              0      0.0217881       0.966236       0.204573              1             45          40014
                         */
                        
                        
                        
                    } // for on muons
                    
                    weight[j] = weight[j]/ (float) kFactor;
                    weight_BB[j] = weight_BB[j]/(float) kFactor_BB;
                    weight_BB[j] = weight_BB[j]/(float) kFactor_BE;
                    
                }// end of DY
                
                
                
                
                weight[j] /= genWeight;
                weight_BB[j] /= genWeight;
                weight_BE[j] /= genWeight;
                
                
                
            } //end selection on entries
            
        }//end looping over entries
        
    } //end looping over MC samples
    
    
    
    
    
    
    ///////////////////////////////////DATA//////////////////////////
    
    
    prev_event=-88;
    
    //strat of looping over Data samples
    for(int j=0; j<data; j++){
        
        std::cout<<"opening.. "<<DATA_samples[j]<<std::endl;
        
        TChain *treeDATA = new TChain("SimpleNtupler/t");
        
        treeDATA->Add("/eos/user/k/kaliyana/2017_data/"+ samp[j] +"/"+ DATA_samples[j] +".root");
        Long64_t ne = treeDATA->GetEntries();
        
        //treeDATA->SetBranchAddress("genWeight",&genWeight);
        treeDATA->SetBranchAddress("event",&event);
        treeDATA->SetBranchAddress("run",&run);
        treeDATA->SetBranchAddress("lumi",&lumi);
        treeDATA->SetBranchAddress("dil_mass",&dil_mass);
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
            if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            treeDATA->GetEntry(p);
            
            //start selection on entries
            if(
               //GoodVtx &&
               fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4 &&
               lep_pt[0]>53. && lep_pt[1]>53. && //offline reconstruction pt threshold
               lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1 &&
               lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1 &&
               fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2 &&
               //(lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || (lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16) && lep_numberOfMatchedRPCLayers[0] > 2))  &&
               //(lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || (lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16) && lep_numberOfMatchedRPCLayers[1] > 2)) &&
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
               ){
                
                /*  if(j==4 && event==67760){
                 std::cout<< "event: "<<event<< " lumi: "<<lumi<< " run: "<<run<<std::endl;
                 std::cout<< "dil_pt: "<<dil_pt<<std::endl;
                 std::cout<< "vertex_m: "<<vertex_m<<std::endl;
                 std::cout<< "lep_pt[0]: "<<lep_pt[0]<<std::endl;
                 std::cout<< "lep_pt[1]: "<<lep_pt[1]<<std::endl;
                 std::cout<< "--------------------------------------"<<std::endl;
                 }
                 */
                
                //std::cout<<dil_mass<<std::endl;
                if(prev_event == event) {
                    //std::cout<< "event: "<<event<<std::endl;
                    //std::cout<< "prev_event: "<<prev_event<<std::endl;
                    continue; //to remove double event; take the first --> they are already sorted in pt
                }
                
                prev_event = event;
                
                
                if(dil_mass < 0){
                    std::cout<<DATA_samples[j]<<"   "<<run<<"   "<<lumi<<"  "<<event<<" "<<dil_mass<<std::endl;
                    continue;
                }
            
                
                // std::cout<<p<<"\t"<<event<<"\t"<<dil_mass<<std::endl;
                
                DATA_pt->Fill(dil_pt);
                DATA_rapidity->Fill(dil_rap);
                DATA_phi->Fill(dil_phi);
                Vtx_DATA_clear->Fill(vertex_chi2);
                
                
                if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
                    DATA_rapidity_BB->Fill(dil_rap);
                }
                else{
                    DATA_rapidity_BE->Fill(dil_rap);
                }
                
                for(int h = 0; h < 2; h++){
                    
                    DATA_lep_pt->Fill(lep_pt[h]);
                    DATA_lep_eta->Fill(lep_eta[h]);
                    DATA_lep_phi->Fill(lep_phi[h]);
                    dB_DATA_clear->Fill(lep_dB[h]);
                    PixelHit_DATA_clear->Fill(lep_glb_numberOfValidPixelHits[h]);
                    TkLayer_DATA_clear->Fill(lep_glb_numberOfValidTrackerLayers[h]);
                    Iso_DATA_clear->Fill(lep_sumPt[h]/lep_tk_pt[h]);
                    relpTErr_DATA_clear->Fill(lep_pt_err[h]/lep_pt[h]);
                    ValidMu_DATA_clear->Fill(lep_glb_numberOfValidMuonHits[h]);
                    
                    if(fabs(lep_eta[h]) < 1.2){
                        DATA_lep_pt_BB->Fill(lep_pt[h]);
                        DATA_lep_eta_BB->Fill(lep_eta[h]);
                    }
                    else{
                        DATA_lep_pt_BE->Fill(lep_pt[h]);
                        DATA_lep_eta_BE->Fill(lep_eta[h]);
                    }
                    
                    if(gen_lep_qOverPt[h] > 0)
                        DATA_lep_pt_plus->Fill(lep_pt[h]);
                    else
                        DATA_lep_pt_minus->Fill(lep_pt[h]);
                    
                  //  if (lep_pt[h] >= 1000.0){
                     //   std::cout<< "lep_pt" << "  "<< "vertex_m"<<std::endl;
                      //  std::cout<< lep_pt[h] << "  "<< vertex_m<<" "<<lep_eta[h]<<std::endl;
                   // }
                    
                    
                } // for on muons
                
                
                ///////////////////
                
                
                
            } //end selection on entries
            
        }//end looping over entries
        
    } //end looping over Data
    
    
    
    f->Write();
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"DATA_lep_pt"<<" "<<"DATA_lep_pt_Barrel"<<" "<<"DATA_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<DATA_lep_pt->GetBinContent(i)<<" "<<DATA_lep_pt_BB->GetBinContent(i)<<" "<<DATA_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"DY_lep_pt"<<" "<<"DY_lep_pt_Barrel"<<" "<<"DY_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<DY_lep_pt->GetBinContent(i)<<" "<<DY_lep_pt_BB->GetBinContent(i)<<" "<<DY_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"ttbar_lep_pt"<<" "<<"ttbar_lep_pt_Barrel"<<" "<<"ttbar_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<ttbar_lep_pt->GetBinContent(i)<<" "<<ttbar_lep_pt_BB->GetBinContent(i)<<" "<<ttbar_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"tW_lep_pt"<<" "<<"tW_lep_pt_Barrel"<<" "<<"tW_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<tW_lep_pt->GetBinContent(i)<<" "<<tW_lep_pt_BB->GetBinContent(i)<<" "<<tW_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"Wantitop_lep_pt"<<" "<<"Wantitop_lep_pt_Barrel"<<" "<<"Wantitop_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<Wantitop_lep_pt->GetBinContent(i)<<" "<<Wantitop_lep_pt_BB->GetBinContent(i)<<" "<<Wantitop_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"WW_lep_pt"<<" "<<"WW_lep_pt_Barrel"<<" "<<"WW_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<WW_lep_pt->GetBinContent(i)<<" "<<WW_lep_pt_BB->GetBinContent(i)<<" "<<WW_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"WZ_lep_pt"<<" "<<"WZ_lep_pt_Barrel"<<" "<<"WZ_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<WZ_lep_pt->GetBinContent(i)<<" "<<WZ_lep_pt_BB->GetBinContent(i)<<" "<<WZ_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"ZZ_lep_pt"<<" "<<"ZZ_lep_pt_Barrel"<<" "<<"ZZ_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<ZZ_lep_pt->GetBinContent(i)<<" "<<ZZ_lep_pt_BB->GetBinContent(i)<<" "<<ZZ_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"Wjets_lep_pt"<<" "<<"Wjets_lep_pt_Barrel"<<" "<<"Wjets_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<Wjets_lep_pt->GetBinContent(i)<<" "<<Wjets_lep_pt_BB->GetBinContent(i)<<" "<<Wjets_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    
    std::cout<< "Pt bin values"<<std::endl;
    
    std::cout<<"Bin No"<<" "<<"Zjets_lep_pt"<<" "<<"Zjets_lep_pt_Barrel"<<" "<<"Zjets_lep_pt_Endcap"<<std::endl;
    
    for(int i=1; i<=36; i++){
        
        std::cout<<i<<" "<<Zjets_lep_pt->GetBinContent(i)<<" "<<Zjets_lep_pt_BB->GetBinContent(i)<<" "<<Zjets_lep_pt_BE->GetBinContent(i)<<std::endl;
        
    }
    
    ///////////////////////////////////////////////////////////////////////////////////
    
    
    //DrawPlot(DY_pt, ttbar_pt, tW_pt, Wantitop_pt, WW_pt, WZ_pt, ZZ_pt, Wjets_pt, Zjets_pt,  Zprime_pt, DATA_pt, dilepton_pt, "pT[GeV](#mu^{+}#mu^{-})", "DATA_MC_pt", 0);
    
    // DrawPlot(DY_rapidity_BB, ttbar_rapidity_BB, tW_rapidity_BB, Wantitop_rapidity_BB, WW_rapidity_BB, WZ_rapidity_BB, ZZ_rapidity_BB, Wjets_rapidity_BB, Zjets_rapidity_BB, Zprime_rapidity_BB, DATA_rapidity_BB, dilepton_rapidity_BB, "rapidity_BB(#mu^{+}#mu^{-})", "DATA_MC_rapidity_BB", 0);
    
    //DrawPlot(DY_rapidity_BE, ttbar_rapidity_BE, tW_rapidity_BE, Wantitop_rapidity_BE, WW_rapidity_BE, WZ_rapidity_BE, ZZ_rapidity_BE, Wjets_rapidity_BE, Zjets_rapidity_BE,  Zprime_rapidity_BE, DATA_rapidity_BE, dilepton_rapidity_BE, "rapidity_BE(#mu^{+}#mu^{-})", "DATA_MC_rapidity_BE", 0);
    
    // DrawPlot(DY_rapidity, ttbar_rapidity, tW_rapidity, Wantitop_rapidity, WW_rapidity, WZ_rapidity, ZZ_rapidity, Wjets_rapidity, Zjets_rapidity,  Zprime_rapidity,  DATA_rapidity, dilepton_rapidity, "rapidity(#mu^{+}#mu^{-})", "DATA_MC_rapidity", 0);
    
    //  DrawPlot(DY_phi, ttbar_phi, tW_phi, Wantitop_phi, WW_phi, WZ_phi, ZZ_phi, Wjets_phi, Zjets_phi,  Zprime_phi, DATA_phi, dilepton_phi, "phi", "DATA_MC_phi", 0);
    
    //  DrawPlot(Vtx_DY_clear, Vtx_ttbar_clear, Vtx_tW_clear, Vtx_Wantitop_clear, Vtx_WW_clear, Vtx_WZ_clear, Vtx_ZZ_clear, Vtx_Wjets_clear, Vtx_Zjets_clear,  Vtx_Zprime_clear, Vtx_DATA_clear, Vtx_clear, "vertex_chi2(#mu^{+}#mu^{-})", "DATA_MC_vtx", 0);
    
    
    DrawPlot(DY_lep_pt, ttbar_lep_pt, tW_lep_pt, Wantitop_lep_pt, WW_lep_pt, WZ_lep_pt, ZZ_lep_pt, Wjets_lep_pt, Zjets_lep_pt, Zprime_lep_pt, DATA_lep_pt, lepton_pt, "pT[GeV]", "DATA_MC_lep_pt", 0);
    
    DrawPlot(DY_lep_pt_BB, ttbar_lep_pt_BB, tW_lep_pt_BB, Wantitop_lep_pt_BB, WW_lep_pt_BB, WZ_lep_pt_BB, ZZ_lep_pt_BB, Wjets_lep_pt_BB, Zjets_lep_pt_BB,  Zprime_lep_pt_BB, DATA_lep_pt_BB, lepton_pt_BB, "pT_Barrel[GeV]", "DATA_MC_lep_pt_Barrel", 0);
    
    DrawPlot(DY_lep_pt_BE, ttbar_lep_pt_BE, tW_lep_pt_BE, Wantitop_lep_pt_BE, WW_lep_pt_BE, WZ_lep_pt_BE, ZZ_lep_pt_BE, Wjets_lep_pt_BE, Zjets_lep_pt_BE,   Zprime_lep_pt_BE, DATA_lep_pt_BE, lepton_pt_BE, "pT_Endcap[GeV]", "DATA_MC_lep_pt_Endcap", 0);
    
    //  DrawPlot(DY_lep_pt_plus, ttbar_lep_pt_plus, tW_lep_pt_plus, Wantitop_lep_pt_plus, WW_lep_pt_plus, WZ_lep_pt_plus, ZZ_lep_pt_plus, Wjets_lep_pt_plus, Zjets_lep_pt_plus, Zprime_lep_pt_plus, DATA_lep_pt_plus, lepton_pt_plus, "pT_(#mu^{+})[GeV]", "DATA_MC_lep_pt_plus", 0);
    
    //  DrawPlot(DY_lep_pt_minus, ttbar_lep_pt_minus, tW_lep_pt_minus, Wantitop_lep_pt_minus, WW_lep_pt_minus, WZ_lep_pt_minus, ZZ_lep_pt_minus, Wjets_lep_pt_minus, Zjets_lep_pt_minus, Zprime_lep_pt_minus, DATA_lep_pt_minus, lepton_pt_minus, "pT_(#mu^{-})[GeV]", "DATA_MC_lep_pt_minus", 0);
    
    //   DrawPlot(DY_lep_eta, ttbar_lep_eta, tW_lep_eta, Wantitop_lep_eta, WW_lep_eta, WZ_lep_eta, ZZ_lep_eta,  Wjets_lep_eta, Zjets_lep_eta, Zprime_lep_eta, DATA_lep_eta, lepton_eta, "lep_eta", "DATA_MC_lep_eta", 0);
    
    // DrawPlot(DY_lep_eta_BB, ttbar_lep_eta_BB, tW_lep_eta_BB, Wantitop_lep_eta_BB, WW_lep_eta_BB, WZ_lep_eta_BB, ZZ_lep_eta_BB,  Wjets_lep_eta_BB, Zjets_lep_eta_BB, Zprime_lep_eta_BB, DATA_lep_eta_BB, lepton_eta_BB, "lep_eta_Barrel", "DATA_MC_lep_eta_Barrel", 0);
    
    //  DrawPlot(DY_lep_eta_BE, ttbar_lep_eta_BE, tW_lep_eta_BE, Wantitop_lep_eta_BE, WW_lep_eta_BE, WZ_lep_eta_BE, ZZ_lep_eta_BE,  Wjets_lep_eta_BE, Zjets_lep_eta_BE, Zprime_lep_eta_BE, DATA_lep_eta_BE, lepton_eta_BE, "lep_eta_Endcap", "DATA_MC_lep_eta_Endcap", 0);
    
    //  DrawPlot(DY_lep_phi, ttbar_lep_phi, tW_lep_phi, Wantitop_lep_phi, WW_lep_phi, WZ_lep_phi, ZZ_lep_phi,  Wjets_lep_phi, Zjets_lep_phi,   Zprime_lep_phi, DATA_lep_phi, lepton_phi, "lep_phi", "DATA_MC_lep_phi", 0);
    
    //  DrawPlot(dB_DY_clear, dB_ttbar_clear, dB_tW_clear, dB_Wantitop_clear, dB_WW_clear, dB_WZ_clear, dB_ZZ_clear, dB_Wjets_clear, dB_Zjets_clear, dB_Zprime_clear, dB_DATA_clear, dB_clear, "dB(cm)", "DATA_MC_dB", 0);
    
    //   DrawPlot(PixelHit_DY_clear, PixelHit_ttbar_clear, PixelHit_tW_clear, PixelHit_Wantitop_clear, PixelHit_WW_clear, PixelHit_WZ_clear, PixelHit_ZZ_clear, PixelHit_Wjets_clear, PixelHit_Zjets_clear, PixelHit_Zprime_clear, PixelHit_DATA_clear, PixelHit_clear, "No.of_PixelHits", "DATA_MC_PixelHit", 0);
    
    //  DrawPlot(TkLayer_DY_clear, TkLayer_ttbar_clear, TkLayer_tW_clear, TkLayer_Wantitop_clear, TkLayer_WW_clear, TkLayer_WZ_clear, TkLayer_ZZ_clear, TkLayer_Wjets_clear, TkLayer_Zjets_clear, TkLayer_Zprime_clear, TkLayer_DATA_clear, TkLayer_clear, "No.of_TkLayers", "DATA_MC_TkLayer", 0);
    
    //  DrawPlot(Iso_DY_clear, Iso_ttbar_clear, Iso_tW_clear, Iso_Wantitop_clear, Iso_WW_clear, Iso_WZ_clear, Iso_ZZ_clear, Iso_Wjets_clear, Iso_Zjets_clear,  Iso_Zprime_clear, Iso_DATA_clear, Iso_clear, "relative tracpker-only isolation", "DATA_MC_Iso", 0);
    
    //  DrawPlot(relpTErr_DY_clear, relpTErr_ttbar_clear, relpTErr_tW_clear, relpTErr_Wantitop_clear, relpTErr_WW_clear, relpTErr_WZ_clear, relpTErr_ZZ_clear, relpTErr_Wjets_clear, relpTErr_Zjets_clear,  relpTErr_Zprime_clear, relpTErr_DATA_clear, relpTErr_clear, "relative pT Err", "DATA_MC_relpTErr", 0);
    
    //  DrawPlot(ValidMu_DY_clear, ValidMu_ttbar_clear, ValidMu_tW_clear, ValidMu_Wantitop_clear, ValidMu_WW_clear, ValidMu_WZ_clear, ValidMu_ZZ_clear, ValidMu_Wjets_clear, ValidMu_Zjets_clear, ValidMu_Zprime_clear, ValidMu_DATA_clear, ValidMu_clear, "glb_numberOfValidMuonHits", "DATA_MC_ValidMu", 0);
    
    
    
    
    f->Close();
    
    
    
    
}






void DrawPlot(TH1F* DY, TH1F* ttbar, TH1F* tW, TH1F* Wantitop, TH1F* WW, TH1F* WZ, TH1F* ZZ, TH1F* Wjets, TH1F* Zjets, /*Th1F* qcd, */TH1F* Zprime, TH1F* DATA, THStack* DATA_MC, TString title, TString name, bool logx){
    
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
    TText *t2 = tText2->AddText("CMS 2017 Data 42.4 fb^{-1} (13TeV)");
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
    ratio->Sumw2();
    TH1F* MC = (TH1F*)DY->Clone();
    MC->Sumw2();
    MC->Add(Wjets);
    MC->Add(ttbar);
    MC->Add(tW);
    MC->Add(Wantitop);
    MC->Add(WW);
    MC->Add(WZ);
    MC->Add(ZZ);
    MC->Add(Zjets);
    //ratio->Add(MC,-1);
    ratio->Divide(MC);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerSize(0.95);
    ratio->GetYaxis()->SetTitle("DATA / Bkg");
    ratio->GetXaxis()->SetTitle(title);
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
    ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
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
    
    
    
    
    TLine* line = new TLine(Xmin, 1.0, Xmax, 1.0);
    line->SetLineColor(kGreen);
    line->SetLineWidth(1);
    if(logx) line->DrawLine(pow(10,Xmin), 1.0, pow(10,Xmax), 1.0);
    else line->DrawLine(Xmin, 1.0, Xmax, 1.0);
    pad2->Update();
    
    
    
    c3->Update();
    
    c3->SaveAs(name+".pdf","pdf");
    c3->SaveAs(name+".png","png");
    
    
}

/*
 
 void PrintKine(float array[600][20]){
 
 std::cout<<"Printing muon kinematics"<<std::endl;
 
 std::cout<<"No"<<";"<<"Mass[GeV]"<<";"<<"lep_id_muon(-)"<<";"<<"lep_id_muon(+)"<<";"<<"TuneP_pt[GeV]_muon(-)"<<";"<<"TuneP_pt[GeV]_muon(+)"<<";"<<"Tracker_pt[GeV]_muon(-)"<<";"<<"Tracker_pt[GeV]_muon(+)"<<";"<<"Ratio(TuneP/Tracker)_muon(-)"<<";"<<"Ratio(TuneP/Tracker)_muon(+)"<<";"<<"eta_muon(-)"<<";"<<"eta_muon(+)"<<";"<<"phi_muon(-)"<<";"<<"phi_muon(+)"<<";"<<"MET_pt[GeV]"<<";"<<"MET_phi"<<";"<<"Category"<<";"<<"Run"<<";"<<"Lumi"<<";"<<"Event"<<std::endl;
 
 
 std::cout<<std::setw(10)<<"No"<<std::setw(10)<<"Mass[GeV]"<<std::setw(10)<<"TuneP_pt[GeV]_muon(-)"<<std::setw(10)<<"TuneP_pt[GeV]_muon(+)"<<std::setw(10)<<"Tracker_pt[GeV]_muon(-)"<<std::setw(10)<<"Tracker_pt[GeV]_muon(+)"<<std::setw(10)<<"Ratio(TuneP/Tracker)_muon(-)"<<std::setw(10)<<"Ratio(TuneP/Tracker)_muon(+)"<<std::setw(10)<<"eta_muon(-)"<<std::setw(10)<<"eta_muon(+)"<<std::setw(10)<<"phi_muon(-)"<<std::setw(10)<<"phi_muon(+)"<<std::setw(10)<<"MET_pt[GeV]"<<std::setw(10)<<"MET_phi"<<std::setw(10)<<"Category"<<std::setw(10)<<"Run"<<std::setw(10)<<"Lumi"<<std::setw(10)<<"Event"<<std::endl;
 
 
 for(int r=0; r<600; r++){
 
 for(int c=0; c<20; c++){
 std::cout<<array[r][c]<<";";
 }
 std::cout<<std::endl;
 }
 
 }
 
 */








