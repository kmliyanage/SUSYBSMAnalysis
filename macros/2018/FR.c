//
//  FR.c
//
//
//  Created by Kalpanie Liyanage on 10/15/19.
//

#include <stdio.h>
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
#include "TH1D.h"
#include "TPaveStats.h"
#include "TPad.h"
#include <TString.h>
#include "TEfficiency.h"

TFile *file1 = new TFile("out.root", "RECREATE");

void DrawPlot(TH1D* qcd, TH1D* Wjets, TH1D* Zjets, TH1D* ttbar, TH1D* WW, TH1D* WZ, TH1D* ZZ, TH1D* tW, TH1D* Wantitop, TH1D* dy, TH1D* DATA, THStack* DATA_MC, TString title, TString name, bool logx);

void RatePlot(TH1D* Wjets_den, TH1D* Zjets_den, TH1D* ttbar_den, TH1D* WW_den, TH1D* WZ_den, TH1D* ZZ_den, TH1D* tW_den, TH1D* Wantitop_den, TH1D* dy_den, TH1D* Data_den, TH1D* Wjets_num, TH1D* Zjets_num, TH1D* ttbar_num, TH1D* WW_num, TH1D* WZ_num, TH1D* ZZ_num, TH1D* tW_num, TH1D* Wantitop_num, TH1D* dy_num, TH1D* qcd_den, TH1D* qcd_num, TH1D* Data_num, TString title, TString name);

void QCDPlot(TH1D* qcd_den, TH1D* qcd_num, TString name);

//void EWKPlot(TH1D* Wjets_den, TH1D* Zjets_den, TH1D* ttbar_den, TH1D* WW_den, TH1D* WZ_den, TH1D* ZZ_den, TH1D* tW_den, TH1D* Wantitop_den, TH1D* dy_den, TH1D* Wjets_num, TH1D* Zjets_num, TH1D* ttbar_num, TH1D* WW_num, TH1D* WZ_num, TH1D* ZZ_num, TH1D* tW_num, TH1D* Wantitop_num, TH1D* dy_num, TString title, TString name);

//void calRate(TH1D* Wjets_den, TH1D* Zjets_den, TH1D* ttbar_den, TH1D* WW_den, TH1D* WZ_den, TH1D* ZZ_den, TH1D* tW_den, TH1D* Wantitop_den, TH1D* dy_den, TH1D* Data_den, TH1D* Wjets_num, TH1D* Zjets_num, TH1D* ttbar_num, TH1D* WW_num, TH1D* WZ_num, TH1D* ZZ_num, TH1D* tW_num, TH1D* Wantitop_num, TH1D* dy_num, TH1D* Data_num, TString name);

//void DrawHist(TH1D* qcd, TH1D* Wjets, TH1D* Zjets, TH1D* ttbar, TH1D* WW, TH1D* WZ, TH1D* ZZ, TH1D* tW, TH1D* Wantitop, TH1D* dy, TH1D* DATA, THStack* DATA_MC, TString title, TString name, bool logx);

//void PrintKine(double data_kinemat[400][20]);

void FR(){
    
    const int mc = 33;
    const int data = 4;
    const double PI = 3.141592654;
    
    double data_kinemat[400][20]={0};
    
    //k-factor applied for DY samples only to accamodate NNPDF corrections and PI bkgs
    double gM;
    double kFactor;
    double kFactor_BB;
    double kFactor_BE;
    
    TString MC_samples[mc] =  {
        "qcd15to30",        //0
        "qcd30to50",        //1
        "qcd50to80",        //2
        "qcd80to120",       //3
        "qcd120to170",      //4
        "qcd170to300",      //5
        "qcd300to470",      //6
        "qcd470to600",      //7
        "qcd600to800",      //8
        "qcd800to1000",     //9
        "qcd1000to1400",    //10
        "qcd1400to1800",    //11
        "qcd1800to2400",    //12
        "qcd2400to3200",    //13
        "qcd3200toInf",     //14
        "Wjets",            //15
        "Zjets_filtered",   //16
        "ttbar",            //17
        "WW",               //18
        "WZ",               //19
        "ZZ",               //20
        "tW",               //21
        "Wantitop",         //22
        "dy50to120",    //23
        "dy120to200",   //24
        "dy200to400",   //25
        "dy400to800",   //26
        "dy800to1400",  //27
        "dy1400to2300", //28
        "dy2300to3500", //29
        "dy3500to4500", //30
        "dy4500to6000", //31
        "dy6000toInf",  //32
    };
    
    double events[mc] = {
        19451000,   //0 qcd15to30
        18872000,   //1 qcd30to50
        12909000,   //2 qcd50to80
        29535000,   //3 qcd80to120
        25255000,   //4 qcd120to170
        29710000,   //5 qcd170to300
        41744000,   //6 qcd300to470
        17712000,   //7 qcd470to600
        64061000,   //8 qcd600to800
        37598000,   //9 qcd800to1000
        18485000,   //10    qcd1000to1400
        2160000,    //11    qcd1400to1800
        1445800,    //12    qcd1800to2400
        1440000,    //13    qcd2400to3200
        800000,     //14    qcd3200toInf  151000
        71026861,   //15    Wjets
        31470710,  //16    Zjets
        64310000,   //17    ttbar
        7850000,    //18    WW
        3885000,    //19    WZ
        1979000,    //20    ZZ
        9598000,    //21    tW
        7623000,    //22    Wantitop
        2982000,    //23    dy50to120
        100000,     //24    dy120to200
        100000,     //25    dy200to400
        100000,     //26    dy400to800
        100000,     //27    dy800to1400
        100000,     //28    dy1400to2300
        100000,     //29    dy2300to3500
        100000,     //30    dy3500to4500
        100000,     //31    dy4500to6000
        100000,     //32    dy6000toInf
        
        
    };
    
    //Using cross sections in XSDB
    double sigma[mc] = {
        1246000000.0,   //0 qcd15to30
        106900000.0,    //1 qcd30to50
        15710000.0,     //2 qcd50to80
        2336000.0,      //3 qcd80to120
        407300.0,       //4 qcd120to170
        103500.0,       //5 qcd170to300
        6830.0,         //6 qcd300to470
        552.1,          //7 qcd470to600
        156.5,          //8 qcd600to800
        26.28,          //9 qcd800to1000
        7.477,           //10    qcd1000to1400
        0.6484,         //11    qcd1400to1800
        0.0875,        //12    qcd1800to2400
        0.005236,       //13    qcd2400to3200
        0.0001357,      //14    qcd3200toInf
        61526.7,        //15    Wjets
        6077.22,        //16    Zjets
        88.29,          //17    ttbar
        118.7,           //18    WW
        50.2,           //19    WZ
        16.523,          //20    ZZ
        35.6,          //21    tW
        35.6,          //22    Wantitop
        2112.904,         //23    dy50to120
        20.553,          //24    dy120to200
        2.886,          //25    dy200to400
        0.2517,         //26    dy400to800
        0.01707,        //27    dy800to1400
        1.366E-3,       //28    dy1400to2300
        8.178E-5,       //29    dy2300to3500
        3.191E-6,       //30    dy3500to4500
        2.787E-7,       //31    dy4500to6000
        9.569E-9,       //32    dy6000toInf
    };
    
    //Using cross sections in XSDB and AN_2018_011
    /*  double sigma[mc] = {
     1246000000.0,   //0 qcd15to30
     106900000.0,    //1 qcd30to50
     15710000.0,     //2 qcd50to80
     2336000.0,      //3 qcd80to120
     407300.0,       //4 qcd120to170
     103500.0,       //5 qcd170to300
     6830.0,         //6 qcd300to470
     552.1,          //7 qcd470to600
     156.5,          //8 qcd600to800
     26.28,          //9 qcd800to1000
     7.47,           //10    qcd1000to1400
     0.6484,         //11    qcd1400to1800
     0.08743,        //12    qcd1800to2400
     0.005236,       //13    qcd2400to3200
     0.0001357,      //14    qcd3200toInf
     52940.0,        //15    Wjets
     5765.4,        //16    Zjets
     87.31,          //17    ttbar
     118.7,           //18    WW
     47.13,           //19    WZ
     16.523,          //20    ZZ
     35.6,          //21    tW
     35.6,          //22    Wantitop
     2112.904,         //23    dy50to120
     20.553,          //24    dy120to200
     2.886,          //25    dy200to400
     0.2517,         //26    dy400to800
     0.01707,        //27    dy800to1400
     1.366E-3,       //28    dy1400to2300
     8.178E-5,       //29    dy2300to3500
     3.191E-6,       //30    dy3500to4500
     2.787E-7,       //31    dy4500to6000
     9.569E-9,       //32    dy6000toInf
     };
     */
    
    
    
    
    
    
    
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
    
    double LUMINOSITY = 61298.77523; //in pb (14217.529374833 + 6867.781918123 + 6612.242975873 +  32717.848260849)
    
    double weight[mc] = {0};
    double weight_BB = 1.0;
    double weight_BE = 1.0;
    
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
    float leading_phi = -1.0;
    float delta_phi = -1.0;
    Int_t c_event;
    Int_t n_event;
    
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
    float lep_glb_pt[2];
    float lep_picky_pt[2];
    float lep_tpfms_pt[2];
    float lep_dB[2];
    float lep_tk_dz[2];
    float lep_tuneP_pt[2];
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
    double Z_peak = 1.0062;
    double Z_peak_BB = 1.0124;
    double Z_peak_BE = 1.0017;
    
    
    
    std::cout<<"MC_sample"<<" "<<"Weight"<<std::endl;
    for(int i = 0; i<mc; i++){
        weight[i] = LUMINOSITY * sigma[i] / events[i];
        std::cout<<MC_samples[i]<<" "<<weight[i]<<std::endl;
    }
    
    int highPt_MC[mc] = {0};
    int highPt_Data[data] = {0};
    int loose = 0;
    
    int den_MC[mc][2] = {0};
    int den_DATA[data][2] = {0};
    
    int num_MC[mc][2] = {0};
    int num_DATA[data][2] = {0};
    
    int dimu_MC[mc][2] = {0};
    int dimu_DATA[data][2] = {0};
    double pt_bins[9] = {50.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0, 1000.0, 5000.0};
    Int_t  binnum_pt = sizeof(pt_bins)/sizeof(Double_t)-1;
    //Int_t  binnum_pt = 7;
    
    const int    NMBINS = 51;
    const double MMIN = 70., MMAX =4000.;
    double logMbins[NMBINS+1];
    
    for (int ibin = 0; ibin <= NMBINS; ibin++){
        logMbins[ibin] = exp(log(MMIN) + (log(MMAX)-log(MMIN))*ibin/NMBINS);
        //       if(ibin != 0 && logMbins[ibin] <= logMbins[ibin-1]) std::cout<<"PROBLEMA"<<std::endl;
        //ss  std::cout<<logMbins[ibin]<<",";
    }
    
    
    
    ///////////////Stack plots//////////////////////////////////////
    THStack* num_pt_Barrel = new THStack("num_pt_Barrel", "num_pt_Barrel");
    THStack* den_pt_Barrel = new THStack("den_pt_Barrel", "den_pt_Barrel");
    THStack* num_pt_Endcap = new THStack("num_pt_Endcap", "num_pt_Endcap");
    THStack* den_pt_Endcap = new THStack("den_pt_Endcap", "den_pt_Endcap");
    
    THStack* num_eta = new THStack("num_eta", "num_eta");
    THStack* den_eta = new THStack("den_eta", "den_eta");
    THStack* num_eta_Barrel = new THStack("num_eta_Barrel", "num_eta_Barrel");
    THStack* den_eta_Barrel = new THStack("den_eta_Barrel", "den_eta_Barrel");
    THStack* num_eta_Endcap = new THStack("num_eta_Endcap", "num_eta_Endcap");
    THStack* den_eta_Endcap = new THStack("den_eta_Endcap", "den_eta_Endcap");
    
    THStack* num_dilmass_Barrel = new THStack("num_dilmass_Barrel", "num_dilmass_Barrel");
    THStack* den_dilmass_Barrel = new THStack("den_dilmass_Barrel", "den_dilmass_Barrel");
    THStack* num_dilmass_Endcap = new THStack("num_dilmass_Endcap", "num_dilmass_Endcap");
    THStack* den_dilmass_Endcap = new THStack("den_dilmass_Endcap", "den_dilmass_Endcap");
    
    
    THStack* num_mass_log = new THStack("num_mass_log", "num_mass_log");
    THStack* num_mass_log_BB = new THStack("num_mass_log_BB", "num_mass_log_BB");
    THStack* num_mass_log_BE = new THStack("num_mass_log_BE", "num_mass_log_BE");
    
    THStack* den_mass_log = new THStack("den_mass_log", "den_mass_log");
    THStack* den_mass_log_BB = new THStack("den_mass_log_BB", "den_mass_log_BB");
    THStack* den_mass_log_BE = new THStack("den_mass_log_BE", "den_mass_log_BE");
    
    
    
    THStack* num_transverse_mass = new THStack("num_transverse_mass", "num_transverse_mass");
    THStack* num_dz = new THStack("num_dz", "num_dz");
    THStack* num_MET_pt = new THStack("num_MET_pt", "num_MET_pt");
    THStack* num_MET_phi = new THStack("num_MET_phi", "num_MET_phi");
    
    THStack* den_transverse_mass = new THStack("den_transverse_mass", "den_transverse_mass");
    THStack* den_dz = new THStack("den_dz", "den_dz");
    THStack* den_MET_pt = new THStack("den_MET_pt", "den_MET_pt");
    THStack* den_MET_phi = new THStack("den_MET_phi", "den_MET_phi");
    
    
    ///////////////Histograms/////////////////////////////////
    
    ////////////////////qcd///////////////////////////////////
    TH1D* qcd_num_pt_Barrel = new TH1D("qcd_num_pt_Barrel", "qcd_num_pt_Barrel", binnum_pt, pt_bins);
    qcd_num_pt_Barrel->Sumw2();
    TH1D* qcd_den_pt_Barrel = new TH1D("qcd_den_pt_Barrel", "qcd_den_pt_Barrel", binnum_pt, pt_bins);
    qcd_den_pt_Barrel->Sumw2();
    TH1D* qcd_num_pt_Endcap = new TH1D("qcd_num_pt_Endcap", "qcd_num_pt_Endcap", binnum_pt, pt_bins);
    qcd_num_pt_Endcap->Sumw2();
    TH1D* qcd_den_pt_Endcap = new TH1D("qcd_den_pt_Endcap", "qcd_den_pt_Endcap", binnum_pt, pt_bins);
    qcd_den_pt_Endcap->Sumw2();
    
    TH1D* qcd_num_eta = new TH1D("qcd_num_eta", "qcd_num_eta", 50, -3, 3);
    TH1D* qcd_den_eta = new TH1D("qcd_den_eta", "qcd_den_eta", 50, -3, 3);
    TH1D* qcd_num_eta_Barrel = new TH1D("qcd_num_eta_Barrel", "qcd_num_eta_Barrel", 50, -3, 3);
    TH1D* qcd_den_eta_Barrel = new TH1D("qcd_den_eta_Barrel", "qcd_den_eta_Barrel", 50, -3, 3);
    TH1D* qcd_num_eta_Endcap = new TH1D("qcd_num_eta_Endcap", "qcd_num_eta_Endcap", 50, -3, 3);
    TH1D* qcd_den_eta_Endcap = new TH1D("qcd_den_eta_Endcap", "qcd_den_eta_Endcap", 50, -3, 3);
    
    TH1D* qcd_num_mass_log = new TH1D("qcd_num_mass_log", "qcd_num_mass_log", NMBINS, logMbins);
    TH1D* qcd_num_mass_log_BB = new TH1D("qcd_num_mass_log_BB", "qcd_num_mass_log_BB", NMBINS, logMbins);
    TH1D* qcd_num_mass_log_BE = new TH1D("qcd_num_mass_log_BE", "qcd_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* qcd_den_mass_log = new TH1D("qcd_den_mass_log", "qcd_den_mass_log", NMBINS, logMbins);
    TH1D* qcd_den_mass_log_BB = new TH1D("qcd_den_mass_log_BB", "qcd_den_mass_log_BB", NMBINS, logMbins);
    TH1D* qcd_den_mass_log_BE = new TH1D("qcd_den_mass_log_BE", "qcd_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_qcd_transverse_mass = new TH1D("num_qcd_transverse_mass", "num_qcd_transverse_mass", 50, 0, 500);
    TH1D* num_qcd_dz = new TH1D("num_qcd_dz", "num_qcd_dz", 50, 0, 10);
    TH1D* num_qcd_met_pt = new TH1D("num_qcd_met_pt", "num_qcd_met_pt", 100, 0, 1000);
    TH1D* num_qcd_met_phi = new TH1D("num_qcd_met_phi", "num_qcd_met_phi", 20, -3, 3);
    
    TH1D* den_qcd_transverse_mass = new TH1D("den_qcd_transverse_mass", "den_qcd_transverse_mass", 50, 0, 500);
    TH1D* den_qcd_dz = new TH1D("den_qcd_dz", "den_qcd_dz", 50, 0, 10);
    TH1D* den_qcd_met_pt = new TH1D("den_qcd_met_pt", "den_qcd_met_pt", 100, 0, 1000);
    TH1D* den_qcd_met_phi = new TH1D("den_qcd_met_phi", "den_qcd_met_phi", 20, -3, 3);
    
    
    /////////////////////////////////Wjets////////////////////////////////////////////
    TH1D* Wjets_num_pt_Barrel = new TH1D("Wjets_num_pt_Barrel", "Wjets_num_pt_Barrel", binnum_pt, pt_bins);
    Wjets_num_pt_Barrel->Sumw2();
    TH1D* Wjets_den_pt_Barrel = new TH1D("Wjets_den_pt_Barrel", "Wjets_den_pt_Barrel", binnum_pt, pt_bins);
    Wjets_den_pt_Barrel->Sumw2();
    TH1D* Wjets_num_pt_Endcap = new TH1D("Wjets_num_pt_Endcap", "Wjets_num_pt_Endcap", binnum_pt, pt_bins);
    Wjets_num_pt_Endcap->Sumw2();
    TH1D* Wjets_den_pt_Endcap = new TH1D("Wjets_den_pt_Endcap", "Wjets_den_pt_Endcap", binnum_pt, pt_bins);
    Wjets_den_pt_Endcap->Sumw2();
    
    TH1D* Wjets_num_eta = new TH1D("Wjets_num_eta", "Wjets_num_eta", 50, -3, 3);
    TH1D* Wjets_den_eta = new TH1D("Wjets_den_eta", "Wjets_den_eta", 50, -3, 3);
    TH1D* Wjets_num_eta_Barrel = new TH1D("Wjets_num_eta_Barrel", "Wjets_num_eta_Barrel", 50, -3, 3);
    TH1D* Wjets_den_eta_Barrel = new TH1D("Wjets_den_eta_Barrel", "Wjets_den_eta_Barrel", 50, -3, 3);
    TH1D* Wjets_num_eta_Endcap = new TH1D("Wjets_num_eta_Endcap", "Wjets_num_eta_Endcap", 50, -3, 3);
    TH1D* Wjets_den_eta_Endcap = new TH1D("Wjets_den_eta_Endcap", "Wjets_den_eta_Endcap", 50, -3, 3);
    
    TH1D* Wjets_num_mass_log = new TH1D("Wjets_num_mass_log", "Wjets_num_mass_log", NMBINS, logMbins);
    TH1D* Wjets_num_mass_log_BB = new TH1D("Wjets_num_mass_log_BB", "Wjets_num_mass_log_BB", NMBINS, logMbins);
    TH1D* Wjets_num_mass_log_BE = new TH1D("Wjets_num_mass_log_BE", "Wjets_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* Wjets_den_mass_log = new TH1D("Wjets_den_mass_log", "Wjets_den_mass_log", NMBINS, logMbins);
    TH1D* Wjets_den_mass_log_BB = new TH1D("Wjets_den_mass_log_BB", "Wjets_den_mass_log_BB", NMBINS, logMbins);
    TH1D* Wjets_den_mass_log_BE = new TH1D("Wjets_den_mass_log_BE", "Wjets_den_mass_log_BE", NMBINS, logMbins);
    
    
    TH1D* num_Wjets_transverse_mass = new TH1D("num_Wjets_transverse_mass", "num_Wjets_transverse_mass", 50, 0, 500);
    TH1D* num_Wjets_dz = new TH1D("num_Wjets_dz", "num_Wjets_dz", 50, 0, 10);
    TH1D* num_Wjets_met_pt = new TH1D("num_Wjets_met_pt", "num_Wjets_met_pt", 100, 0, 1000);
    TH1D* num_Wjets_met_phi = new TH1D("num_Wjets_met_phi", "num_Wjets_met_phi", 20, -3, 3);
    
    TH1D* den_Wjets_transverse_mass = new TH1D("den_Wjets_transverse_mass", "den_Wjets_transverse_mass", 50, 0, 500);
    TH1D* den_Wjets_dz = new TH1D("den_Wjets_dz", "den_Wjets_dz", 50, 0, 10);
    TH1D* den_Wjets_met_pt = new TH1D("den_Wjets_met_pt", "den_Wjets_met_pt", 100, 0, 1000);
    TH1D* den_Wjets_met_phi = new TH1D("den_Wjets_met_phi", "den_Wjets_met_phi", 20, -3, 3);
    
    
    /////////////////////////////////Zjets////////////////////////////////////////////
    TH1D* Zjets_num_pt_Barrel = new TH1D("Zjets_num_pt_Barrel", "Zjets_num_pt_Barrel", binnum_pt, pt_bins);
    Zjets_num_pt_Barrel->Sumw2();
    TH1D* Zjets_den_pt_Barrel = new TH1D("Zjets_den_pt_Barrel", "Zjets_den_pt_Barrel", binnum_pt, pt_bins);
    Zjets_den_pt_Barrel->Sumw2();
    TH1D* Zjets_num_pt_Endcap = new TH1D("Zjets_num_pt_Endcap", "Zjets_num_pt_Endcap", binnum_pt, pt_bins);
    Zjets_num_pt_Endcap->Sumw2();
    TH1D* Zjets_den_pt_Endcap = new TH1D("Zjets_den_pt_Endcap", "Zjets_den_pt_Endcap", binnum_pt, pt_bins);
    Zjets_den_pt_Endcap->Sumw2();
    
    TH1D* Zjets_num_eta = new TH1D("Zjets_num_eta", "Zjets_num_eta", 50, -3, 3);
    TH1D* Zjets_den_eta = new TH1D("Zjets_den_eta", "Zjets_den_eta", 50, -3, 3);
    TH1D* Zjets_num_eta_Barrel = new TH1D("Zjets_num_eta_Barrel", "Zjets_num_eta_Barrel", 50, -3, 3);
    TH1D* Zjets_den_eta_Barrel = new TH1D("Zjets_den_eta_Barrel", "Zjets_den_eta_Barrel", 50, -3, 3);
    TH1D* Zjets_num_eta_Endcap = new TH1D("Zjets_num_eta_Endcap", "Zjets_num_eta_Endcap", 50, -3, 3);
    TH1D* Zjets_den_eta_Endcap = new TH1D("Zjets_den_eta_Endcap", "Zjets_den_eta_Endcap", 50, -3, 3);
    
    TH1D* Zjets_num_mass_log = new TH1D("Zjets_num_mass_log", "Zjets_num_mass_log", NMBINS, logMbins);
    TH1D* Zjets_num_mass_log_BB = new TH1D("Zjets_num_mass_log_BB", "Zjets_num_mass_log_BB", NMBINS, logMbins);
    TH1D* Zjets_num_mass_log_BE = new TH1D("Zjets_num_mass_log_BE", "Zjets_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* Zjets_den_mass_log = new TH1D("Zjets_den_mass_log", "Zjets_den_mass_log", NMBINS, logMbins);
    TH1D* Zjets_den_mass_log_BB = new TH1D("Zjets_den_mass_log_BB", "Zjets_den_mass_log_BB", NMBINS, logMbins);
    TH1D* Zjets_den_mass_log_BE = new TH1D("Zjets_den_mass_log_BE", "Zjets_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_Zjets_transverse_mass = new TH1D("num_Zjets_transverse_mass", "num_Zjets_transverse_mass", 50, 0, 500);
    TH1D* num_Zjets_dz = new TH1D("num_Zjets_dz", "num_Zjets_dz", 50, 0, 10);
    TH1D* num_Zjets_met_pt = new TH1D("num_Zjets_met_pt", "num_Zjets_met_pt", 100, 0, 1000);
    TH1D* num_Zjets_met_phi = new TH1D("num_Zjets_met_phi", "num_Zjets_met_phi", 20, -3, 3);
    
    TH1D* den_Zjets_transverse_mass = new TH1D("den_Zjets_transverse_mass", "den_Zjets_transverse_mass", 50, 0, 500);
    TH1D* den_Zjets_dz = new TH1D("den_Zjets_dz", "den_Zjets_dz", 50, 0, 10);
    TH1D* den_Zjets_met_pt = new TH1D("den_Zjets_met_pt", "den_Zjets_met_pt", 100, 0, 1000);
    TH1D* den_Zjets_met_phi = new TH1D("den_Zjets_met_phi", "den_Zjets_met_phi", 20, -3, 3);
    
    ////////////////////////////////ttbar////////////////////////////////////////////
    TH1D* ttbar_num_pt_Barrel = new TH1D("ttbar_num_pt_Barrel", "ttbar_num_pt_Barrel", binnum_pt, pt_bins);
    ttbar_num_pt_Barrel->Sumw2();
    TH1D* ttbar_den_pt_Barrel = new TH1D("ttbar_den_pt_Barrel", "ttbar_den_pt_Barrel", binnum_pt, pt_bins);
    ttbar_den_pt_Barrel->Sumw2();
    TH1D* ttbar_num_pt_Endcap = new TH1D("ttbar_num_pt_Endcap", "ttbar_num_pt_Endcap", binnum_pt, pt_bins);
    ttbar_num_pt_Endcap->Sumw2();
    TH1D* ttbar_den_pt_Endcap = new TH1D("ttbar_den_pt_Endcap", "ttbar_den_pt_Endcap", binnum_pt, pt_bins);
    ttbar_den_pt_Endcap->Sumw2();
    
    TH1D* ttbar_num_eta = new TH1D("ttbar_num_eta", "ttbar_num_eta", 50, -3, 3);
    TH1D* ttbar_den_eta = new TH1D("ttbar_den_eta", "ttbar_den_eta", 50, -3, 3);
    TH1D* ttbar_num_eta_Barrel = new TH1D("ttbar_num_eta_Barrel", "ttbar_num_eta_Barrel", 50, -3, 3);
    TH1D* ttbar_den_eta_Barrel = new TH1D("ttbar_den_eta_Barrel", "ttbar_den_eta_Barrel", 50, -3, 3);
    TH1D* ttbar_num_eta_Endcap = new TH1D("ttbar_num_eta_Endcap", "ttbar_num_eta_Endcap", 50, -3, 3);
    TH1D* ttbar_den_eta_Endcap = new TH1D("ttbar_den_eta_Endcap", "ttbar_den_eta_Endcap", 50, -3, 3);
    
    TH1D* ttbar_num_mass_log = new TH1D("ttbar_num_mass_log", "ttbar_num_mass_log", NMBINS, logMbins);
    TH1D* ttbar_num_mass_log_BB = new TH1D("ttbar_num_mass_log_BB", "ttbar_num_mass_log_BB", NMBINS, logMbins);
    TH1D* ttbar_num_mass_log_BE = new TH1D("ttbar_num_mass_log_BE", "ttbar_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* ttbar_den_mass_log = new TH1D("ttbar_den_mass_log", "ttbar_den_mass_log", NMBINS, logMbins);
    TH1D* ttbar_den_mass_log_BB = new TH1D("ttbar_den_mass_log_BB", "ttbar_den_mass_log_BB", NMBINS, logMbins);
    TH1D* ttbar_den_mass_log_BE = new TH1D("ttbar_den_mass_log_BE", "ttbar_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_ttbar_transverse_mass = new TH1D("num_ttbar_transverse_mass", "num_ttbar_transverse_mass", 50, 0, 500);
    TH1D* num_ttbar_dz = new TH1D("num_ttbar_dz", "num_ttbar_dz", 50, 0, 10);
    TH1D* num_ttbar_met_pt = new TH1D("num_ttbar_met_pt", "num_ttbar_met_pt", 100, 0, 1000);
    TH1D* num_ttbar_met_phi = new TH1D("num_ttbar_met_phi", "num_ttbar_met_phi", 20, -3, 3);
    
    TH1D* den_ttbar_transverse_mass = new TH1D("den_ttbar_transverse_mass", "den_ttbar_transverse_mass", 50, 0, 500);
    TH1D* den_ttbar_dz = new TH1D("den_ttbar_dz", "den_ttbar_dz", 50, 0, 10);
    TH1D* den_ttbar_met_pt = new TH1D("den_ttbar_met_pt", "den_ttbar_met_pt", 100, 0, 1000);
    TH1D* den_ttbar_met_phi = new TH1D("den_ttbar_met_phi", "den_ttbar_met_phi", 20, -3, 3);
    
    
    /////////////////////////////////WW////////////////////////////////////////////
    TH1D* WW_num_pt_Barrel = new TH1D("WW_num_pt_Barrel", "WW_num_pt_Barrel", binnum_pt, pt_bins);
    WW_num_pt_Barrel->Sumw2();
    TH1D* WW_den_pt_Barrel = new TH1D("WW_den_pt_Barrel", "WW_den_pt_Barrel", binnum_pt, pt_bins);
    WW_den_pt_Barrel->Sumw2();
    TH1D* WW_num_pt_Endcap = new TH1D("WW_num_pt_Endcap", "WW_num_pt_Endcap", binnum_pt, pt_bins);
    WW_num_pt_Endcap->Sumw2();
    TH1D* WW_den_pt_Endcap = new TH1D("WW_den_pt_Endcap", "WW_den_pt_Endcap", binnum_pt, pt_bins);
    WW_den_pt_Endcap->Sumw2();
    
    TH1D* WW_num_eta = new TH1D("WW_num_eta", "WW_num_eta", 50, -3, 3);
    TH1D* WW_den_eta = new TH1D("WW_den_eta", "WW_den_eta", 50, -3, 3);
    TH1D* WW_num_eta_Barrel = new TH1D("WW_num_eta_Barrel", "WW_num_eta_Barrel", 50, -3, 3);
    TH1D* WW_den_eta_Barrel = new TH1D("WW_den_eta_Barrel", "WW_den_eta_Barrel", 50, -3, 3);
    TH1D* WW_num_eta_Endcap = new TH1D("WW_num_eta_Endcap", "WW_num_eta_Endcap", 50, -3, 3);
    TH1D* WW_den_eta_Endcap = new TH1D("WW_den_eta_Endcap", "WW_den_eta_Endcap", 50, -3, 3);
    
    TH1D* WW_num_mass_log = new TH1D("WW_num_mass_log", "WW_num_mass_log", NMBINS, logMbins);
    TH1D* WW_num_mass_log_BB = new TH1D("WW_num_mass_log_BB", "WW_num_mass_log_BB", NMBINS, logMbins);
    TH1D* WW_num_mass_log_BE = new TH1D("WW_num_mass_log_BE", "WW_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* WW_den_mass_log = new TH1D("WW_den_mass_log", "WW_den_mass_log", NMBINS, logMbins);
    TH1D* WW_den_mass_log_BB = new TH1D("WW_den_mass_log_BB", "WW_den_mass_log_BB", NMBINS, logMbins);
    TH1D* WW_den_mass_log_BE = new TH1D("WW_den_mass_log_BE", "WW_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_WW_transverse_mass = new TH1D("num_WW_transverse_mass", "num_WW_transverse_mass", 50, 0, 500);
    TH1D* num_WW_dz = new TH1D("num_WW_dz", "num_WW_dz", 50, 0, 10);
    TH1D* num_WW_met_pt = new TH1D("num_WW_met_pt", "num_WW_met_pt", 100, 0, 1000);
    TH1D* num_WW_met_phi = new TH1D("num_WW_met_phi", "num_WW_met_phi", 20, -3, 3);
    
    TH1D* den_WW_transverse_mass = new TH1D("den_WW_transverse_mass", "den_WW_transverse_mass", 50, 0, 500);
    TH1D* den_WW_dz = new TH1D("den_WW_dz", "den_WW_dz", 50, 0, 10);
    TH1D* den_WW_met_pt = new TH1D("den_WW_met_pt", "den_WW_met_pt", 100, 0, 1000);
    TH1D* den_WW_met_phi = new TH1D("den_WW_met_phi", "den_WW_met_phi", 20, -3, 3);
    
    
    /////////////////////////////////WZ////////////////////////////////////////////
    TH1D* WZ_num_pt_Barrel = new TH1D("WZ_num_pt_Barrel", "WZ_num_pt_Barrel", binnum_pt, pt_bins);
    WZ_num_pt_Barrel->Sumw2();
    TH1D* WZ_den_pt_Barrel = new TH1D("WZ_den_pt_Barrel", "WZ_den_pt_Barrel", binnum_pt, pt_bins);
    WZ_den_pt_Barrel->Sumw2();
    TH1D* WZ_num_pt_Endcap = new TH1D("WZ_num_pt_Endcap", "WZ_num_pt_Endcap", binnum_pt, pt_bins);
    WZ_num_pt_Endcap->Sumw2();
    TH1D* WZ_den_pt_Endcap = new TH1D("WZ_den_pt_Endcap", "WZ_den_pt_Endcap", binnum_pt, pt_bins);
    WZ_den_pt_Endcap->Sumw2();
    
    TH1D* WZ_num_eta = new TH1D("WZ_num_eta", "WZ_num_eta", 50, -3, 3);
    TH1D* WZ_den_eta = new TH1D("WZ_den_eta", "WZ_den_eta", 50, -3, 3);
    TH1D* WZ_num_eta_Barrel = new TH1D("WZ_num_eta_Barrel", "WZ_num_eta_Barrel", 50, -3, 3);
    TH1D* WZ_den_eta_Barrel = new TH1D("WZ_den_eta_Barrel", "WZ_den_eta_Barrel", 50, -3, 3);
    TH1D* WZ_num_eta_Endcap = new TH1D("WZ_num_eta_Endcap", "WZ_num_eta_Endcap", 50, -3, 3);
    TH1D* WZ_den_eta_Endcap = new TH1D("WZ_den_eta_Endcap", "WZ_den_eta_Endcap", 50, -3, 3);
    
    TH1D* WZ_num_mass_log = new TH1D("WZ_num_mass_log", "WZ_num_mass_log", NMBINS, logMbins);
    TH1D* WZ_num_mass_log_BB = new TH1D("WZ_num_mass_log_BB", "WZ_num_mass_log_BB", NMBINS, logMbins);
    TH1D* WZ_num_mass_log_BE = new TH1D("WZ_num_mass_log_BE", "WZ_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* WZ_den_mass_log = new TH1D("WZ_den_mass_log", "WZ_den_mass_log", NMBINS, logMbins);
    TH1D* WZ_den_mass_log_BB = new TH1D("WZ_den_mass_log_BB", "WZ_den_mass_log_BB", NMBINS, logMbins);
    TH1D* WZ_den_mass_log_BE = new TH1D("WZ_den_mass_log_BE", "WZ_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_WZ_transverse_mass = new TH1D("num_WZ_transverse_mass", "num_WZ_transverse_mass", 50, 0, 500);
    TH1D* num_WZ_dz = new TH1D("num_WZ_dz", "num_WZ_dz", 50, 0, 10);
    TH1D* num_WZ_met_pt = new TH1D("num_WZ_met_pt", "num_WZ_met_pt", 100, 0, 1000);
    TH1D* num_WZ_met_phi = new TH1D("num_WZ_met_phi", "num_WZ_met_phi", 20, -3, 3);
    
    TH1D* den_WZ_transverse_mass = new TH1D("den_WZ_transverse_mass", "den_WZ_transverse_mass", 50, 0, 500);
    TH1D* den_WZ_dz = new TH1D("den_WZ_dz", "den_WZ_dz", 50, 0, 10);
    TH1D* den_WZ_met_pt = new TH1D("den_WZ_met_pt", "den_WZ_met_pt", 100, 0, 1000);
    TH1D* den_WZ_met_phi = new TH1D("den_WZ_met_phi", "den_WZ_met_phi", 20, -3, 3);
    
    /////////////////////////////////ZZ////////////////////////////////////////////
    TH1D* ZZ_num_pt_Barrel = new TH1D("ZZ_num_pt_Barrel", "ZZ_num_pt_Barrel", binnum_pt, pt_bins);
    ZZ_num_pt_Barrel->Sumw2();
    TH1D* ZZ_den_pt_Barrel = new TH1D("ZZ_den_pt_Barrel", "ZZ_den_pt_Barrel", binnum_pt, pt_bins);
    ZZ_den_pt_Barrel->Sumw2();
    TH1D* ZZ_num_pt_Endcap = new TH1D("ZZ_num_pt_Endcap", "ZZ_num_pt_Endcap", binnum_pt, pt_bins);
    ZZ_num_pt_Endcap->Sumw2();
    TH1D* ZZ_den_pt_Endcap = new TH1D("ZZ_den_pt_Endcap", "ZZ_den_pt_Endcap", binnum_pt, pt_bins);
    ZZ_den_pt_Endcap->Sumw2();
    
    TH1D* ZZ_num_eta = new TH1D("ZZ_num_eta", "ZZ_num_eta", 50, -3, 3);
    TH1D* ZZ_den_eta = new TH1D("ZZ_den_eta", "ZZ_den_eta", 50, -3, 3);
    TH1D* ZZ_num_eta_Barrel = new TH1D("ZZ_num_eta_Barrel", "ZZ_num_eta_Barrel", 50, -3, 3);
    TH1D* ZZ_den_eta_Barrel = new TH1D("ZZ_den_eta_Barrel", "ZZ_den_eta_Barrel", 50, -3, 3);
    TH1D* ZZ_num_eta_Endcap = new TH1D("ZZ_num_eta_Endcap", "ZZ_num_eta_Endcap", 50, -3, 3);
    TH1D* ZZ_den_eta_Endcap = new TH1D("ZZ_den_eta_Endcap", "ZZ_den_eta_Endcap", 50, -3, 3);
    
    TH1D* ZZ_num_mass_log = new TH1D("ZZ_num_mass_log", "ZZ_num_mass_log", NMBINS, logMbins);
    TH1D* ZZ_num_mass_log_BB = new TH1D("ZZ_num_mass_log_BB", "ZZ_num_mass_log_BB", NMBINS, logMbins);
    TH1D* ZZ_num_mass_log_BE = new TH1D("ZZ_num_mass_log_BE", "ZZ_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* ZZ_den_mass_log = new TH1D("ZZ_den_mass_log", "ZZ_den_mass_log", NMBINS, logMbins);
    TH1D* ZZ_den_mass_log_BB = new TH1D("ZZ_den_mass_log_BB", "ZZ_den_mass_log_BB", NMBINS, logMbins);
    TH1D* ZZ_den_mass_log_BE = new TH1D("ZZ_den_mass_log_BE", "ZZ_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_ZZ_transverse_mass = new TH1D("num_ZZ_transverse_mass", "num_ZZ_transverse_mass", 50, 0, 500);
    TH1D* num_ZZ_dz = new TH1D("num_ZZ_dz", "num_ZZ_dz", 50, 0, 10);
    TH1D* num_ZZ_met_pt = new TH1D("num_ZZ_met_pt", "num_ZZ_met_pt", 100, 0, 1000);
    TH1D* num_ZZ_met_phi = new TH1D("num_ZZ_met_phi", "num_ZZ_met_phi", 20, -3, 3);
    
    TH1D* den_ZZ_transverse_mass = new TH1D("den_ZZ_transverse_mass", "den_ZZ_transverse_mass", 50, 0, 500);
    TH1D* den_ZZ_dz = new TH1D("den_ZZ_dz", "den_ZZ_dz", 50, 0, 10);
    TH1D* den_ZZ_met_pt = new TH1D("den_ZZ_met_pt", "den_ZZ_met_pt", 100, 0, 1000);
    TH1D* den_ZZ_met_phi = new TH1D("den_ZZ_met_phi", "den_ZZ_met_phi", 20, -3, 3);
    
    /////////////////////////////////tW////////////////////////////////////////////
    TH1D* tW_num_pt_Barrel = new TH1D("tW_num_pt_Barrel", "tW_num_pt_Barrel", binnum_pt, pt_bins);
    tW_num_pt_Barrel->Sumw2();
    TH1D* tW_den_pt_Barrel = new TH1D("tW_den_pt_Barrel", "tW_den_pt_Barrel", binnum_pt, pt_bins);
    tW_den_pt_Barrel->Sumw2();
    TH1D* tW_num_pt_Endcap = new TH1D("tW_num_pt_Endcap", "tW_num_pt_Endcap", binnum_pt, pt_bins);
    tW_num_pt_Endcap->Sumw2();
    TH1D* tW_den_pt_Endcap = new TH1D("tW_den_pt_Endcap", "tW_den_pt_Endcap", binnum_pt, pt_bins);
    tW_den_pt_Endcap->Sumw2();
    
    TH1D* tW_num_eta = new TH1D("tW_num_eta", "tW_num_eta", 50, -3, 3);
    TH1D* tW_den_eta = new TH1D("tW_den_eta", "tW_den_eta", 50, -3, 3);
    TH1D* tW_num_eta_Barrel = new TH1D("tW_num_eta_Barrel", "tW_num_eta_Barrel", 50, -3, 3);
    TH1D* tW_den_eta_Barrel = new TH1D("tW_den_eta_Barrel", "tW_den_eta_Barrel", 50, -3, 3);
    TH1D* tW_num_eta_Endcap = new TH1D("tW_num_eta_Endcap", "tW_num_eta_Endcap", 50, -3, 3);
    TH1D* tW_den_eta_Endcap = new TH1D("tW_den_eta_Endcap", "tW_den_eta_Endcap", 50, -3, 3);
    
    TH1D* tW_num_mass_log = new TH1D("tW_num_mass_log", "tW_num_mass_log", NMBINS, logMbins);
    TH1D* tW_num_mass_log_BB = new TH1D("tW_num_mass_log_BB", "tW_num_mass_log_BB", NMBINS, logMbins);
    TH1D* tW_num_mass_log_BE = new TH1D("tW_num_mass_log_BE", "tW_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* tW_den_mass_log = new TH1D("tW_den_mass_log", "tW_den_mass_log", NMBINS, logMbins);
    TH1D* tW_den_mass_log_BB = new TH1D("tW_den_mass_log_BB", "tW_den_mass_log_BB", NMBINS, logMbins);
    TH1D* tW_den_mass_log_BE = new TH1D("tW_den_mass_log_BE", "tW_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_tW_transverse_mass = new TH1D("num_tW_transverse_mass", "num_tW_transverse_mass", 50, 0, 500);
    TH1D* num_tW_dz = new TH1D("num_tW_dz", "num_tW_dz", 50, 0, 10);
    TH1D* num_tW_met_pt = new TH1D("num_tW_met_pt", "num_tW_met_pt", 100, 0, 1000);
    TH1D* num_tW_met_phi = new TH1D("num_tW_met_phi", "num_tW_met_phi", 20, -3, 3);
    
    TH1D* den_tW_transverse_mass = new TH1D("den_tW_transverse_mass", "den_tW_transverse_mass", 50, 0, 500);
    TH1D* den_tW_dz = new TH1D("den_tW_dz", "den_tW_dz", 50, 0, 10);
    TH1D* den_tW_met_pt = new TH1D("den_tW_met_pt", "den_tW_met_pt", 100, 0, 1000);
    TH1D* den_tW_met_phi = new TH1D("den_tW_met_phi", "den_tW_met_phi", 20, -3, 3);
    
    
    /////////////////////////////////Wantitop////////////////////////////////////////////
    TH1D* Wantitop_num_pt_Barrel = new TH1D("Wantitop_num_pt_Barrel", "Wantitop_num_pt_Barrel", binnum_pt, pt_bins);
    Wantitop_num_pt_Barrel->Sumw2();
    TH1D* Wantitop_den_pt_Barrel = new TH1D("Wantitop_den_pt_Barrel", "Wantitop_den_pt_Barrel", binnum_pt, pt_bins);
    Wantitop_den_pt_Barrel->Sumw2();
    TH1D* Wantitop_num_pt_Endcap = new TH1D("Wantitop_num_pt_Endcap", "Wantitop_num_pt_Endcap", binnum_pt, pt_bins);
    Wantitop_num_pt_Endcap->Sumw2();
    TH1D* Wantitop_den_pt_Endcap = new TH1D("Wantitop_den_pt_Endcap", "Wantitop_den_pt_Endcap", binnum_pt, pt_bins);
    Wantitop_den_pt_Endcap->Sumw2();
    
    TH1D* Wantitop_num_eta = new TH1D("Wantitop_num_eta", "Wantitop_num_eta", 50, -3, 3);
    TH1D* Wantitop_den_eta = new TH1D("Wantitop_den_eta", "Wantitop_den_eta", 50, -3, 3);
    TH1D* Wantitop_num_eta_Barrel = new TH1D("Wantitop_num_eta_Barrel", "Wantitop_num_eta_Barrel", 50, -3, 3);
    TH1D* Wantitop_den_eta_Barrel = new TH1D("Wantitop_den_eta_Barrel", "Wantitop_den_eta_Barrel", 50, -3, 3);
    TH1D* Wantitop_num_eta_Endcap = new TH1D("Wantitop_num_eta_Endcap", "Wantitop_num_eta_Endcap", 50, -3, 3);
    TH1D* Wantitop_den_eta_Endcap = new TH1D("Wantitop_den_eta_Endcap", "Wantitop_den_eta_Endcap", 50, -3, 3);
    
    TH1D* Wantitop_num_mass_log = new TH1D("Wantitop_num_mass_log", "Wantitop_num_mass_log", NMBINS, logMbins);
    TH1D* Wantitop_num_mass_log_BB = new TH1D("Wantitop_num_mass_log_BB", "Wantitop_num_mass_log_BB", NMBINS, logMbins);
    TH1D* Wantitop_num_mass_log_BE = new TH1D("Wantitop_num_mass_log_BE", "Wantitop_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* Wantitop_den_mass_log = new TH1D("Wantitop_den_mass_log", "Wantitop_den_mass_log", NMBINS, logMbins);
    TH1D* Wantitop_den_mass_log_BB = new TH1D("Wantitop_den_mass_log_BB", "Wantitop_den_mass_log_BB", NMBINS, logMbins);
    TH1D* Wantitop_den_mass_log_BE = new TH1D("Wantitop_den_mass_log_BE", "Wantitop_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_Wantitop_transverse_mass = new TH1D("num_Wantitop_transverse_mass", "num_Wantitop_transverse_mass", 50, 0, 500);
    TH1D* num_Wantitop_dz = new TH1D("num_Wantitop_dz", "num_Wantitop_dz", 50, 0, 10);
    TH1D* num_Wantitop_met_pt = new TH1D("num_Wantitop_met_pt", "num_Wantitop_met_pt", 100, 0, 1000);
    TH1D* num_Wantitop_met_phi = new TH1D("num_Wantitop_met_phi", "num_Wantitop_met_phi", 20, -3, 3);
    
    TH1D* den_Wantitop_transverse_mass = new TH1D("den_Wantitop_transverse_mass", "den_Wantitop_transverse_mass", 50, 0, 500);
    TH1D* den_Wantitop_dz = new TH1D("den_Wantitop_dz", "den_Wantitop_dz", 50, 0, 10);
    TH1D* den_Wantitop_met_pt = new TH1D("den_Wantitop_met_pt", "den_Wantitop_met_pt", 100, 0, 1000);
    TH1D* den_Wantitop_met_phi = new TH1D("den_Wantitop_met_phi", "den_Wantitop_met_phi", 20, -3, 3);
    
    
    /////////////////////////////////dy////////////////////////////////////////////
    TH1D* dy_num_pt_Barrel = new TH1D("dy_num_pt_Barrel", "dy_num_pt_Barrel", binnum_pt, pt_bins);
    dy_num_pt_Barrel->Sumw2();
    TH1D* dy_den_pt_Barrel = new TH1D("dy_den_pt_Barrel", "dy_den_pt_Barrel", binnum_pt, pt_bins);
    dy_den_pt_Barrel->Sumw2();
    TH1D* dy_num_pt_Endcap = new TH1D("dy_num_pt_Endcap", "dy_num_pt_Endcap", binnum_pt, pt_bins);
    dy_num_pt_Endcap->Sumw2();
    TH1D* dy_den_pt_Endcap = new TH1D("dy_den_pt_Endcap", "dy_den_pt_Endcap", binnum_pt, pt_bins);
    dy_den_pt_Endcap->Sumw2();
    
    TH1D* dy_num_eta = new TH1D("dy_num_eta", "dy_num_eta", 50, -3, 3);
    TH1D* dy_den_eta = new TH1D("dy_den_eta", "dy_den_eta", 50, -3, 3);
    TH1D* dy_num_eta_Barrel = new TH1D("dy_num_eta_Barrel", "dy_num_eta_Barrel", 50, -3, 3);
    TH1D* dy_den_eta_Barrel = new TH1D("dy_den_eta_Barrel", "dy_den_eta_Barrel", 50, -3, 3);
    TH1D* dy_num_eta_Endcap = new TH1D("dy_num_eta_Endcap", "dy_num_eta_Endcap", 50, -3, 3);
    TH1D* dy_den_eta_Endcap = new TH1D("dy_den_eta_Endcap", "dy_den_eta_Endcap", 50, -3, 3);
    
    TH1D* dy_num_mass_log = new TH1D("dy_num_mass_log", "dy_num_mass_log", NMBINS, logMbins);
    TH1D* dy_num_mass_log_BB = new TH1D("dy_num_mass_log_BB", "dy_num_mass_log_BB", NMBINS, logMbins);
    TH1D* dy_num_mass_log_BE = new TH1D("dy_num_mass_log_BE", "dy_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* dy_den_mass_log = new TH1D("dy_den_mass_log", "dy_den_mass_log", NMBINS, logMbins);
    TH1D* dy_den_mass_log_BB = new TH1D("dy_den_mass_log_BB", "dy_den_mass_log_BB", NMBINS, logMbins);
    TH1D* dy_den_mass_log_BE = new TH1D("dy_den_mass_log_BE", "dy_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_dy_transverse_mass = new TH1D("num_dy_transverse_mass", "num_dy_transverse_mass", 50, 0, 500);
    TH1D* num_dy_dz = new TH1D("num_dy_dz", "num_dy_dz", 50, 0, 10);
    TH1D* num_dy_met_pt = new TH1D("num_dy_met_pt", "num_dy_met_pt", 100, 0, 1000);
    TH1D* num_dy_met_phi = new TH1D("num_dy_met_phi", "num_dy_met_phi", 20, -3, 3);
    
    TH1D* den_dy_transverse_mass = new TH1D("den_dy_transverse_mass", "den_dy_transverse_mass", 50, 0, 500);
    TH1D* den_dy_dz = new TH1D("den_dy_dz", "den_dy_dz", 50, 0, 10);
    TH1D* den_dy_met_pt = new TH1D("den_dy_met_pt", "den_dy_met_pt", 100, 0, 1000);
    TH1D* den_dy_met_phi = new TH1D("den_dy_met_phi", "den_dy_met_phi", 20, -3, 3);
    
    /////////////////////////////////Data////////////////////////////////////////////
    TH1D* Data_num_pt_Barrel = new TH1D("Data_num_pt_Barrel", "Data_num_pt_Barrel", binnum_pt, pt_bins);
    Data_num_pt_Barrel->Sumw2();
    TH1D* Data_den_pt_Barrel = new TH1D("Data_den_pt_Barrel", "Data_den_pt_Barrel", binnum_pt, pt_bins);
    Data_den_pt_Barrel->Sumw2();
    TH1D* Data_num_pt_Endcap = new TH1D("Data_num_pt_Endcap", "Data_num_pt_Endcap", binnum_pt, pt_bins);
    Data_num_pt_Endcap->Sumw2();
    TH1D* Data_den_pt_Endcap = new TH1D("Data_den_pt_Endcap", "Data_den_pt_Endcap", binnum_pt, pt_bins);
    Data_den_pt_Endcap->Sumw2();
    
    TH1D* Data_num_eta = new TH1D("Data_num_eta", "Data_num_eta", 50, -3, 3);
    TH1D* Data_den_eta = new TH1D("Data_den_eta", "Data_den_eta", 50, -3, 3);
    TH1D* Data_num_eta_Barrel = new TH1D("Data_num_eta_Barrel", "Data_num_eta_Barrel", 50, -3, 3);
    TH1D* Data_den_eta_Barrel = new TH1D("Data_den_eta_Barrel", "Data_den_eta_Barrel", 50, -3, 3);
    TH1D* Data_num_eta_Endcap = new TH1D("Data_num_eta_Endcap", "Data_num_eta_Endcap", 50, -3, 3);
    TH1D* Data_den_eta_Endcap = new TH1D("Data_den_eta_Endcap", "Data_den_eta_Endcap", 50, -3, 3);
    
    TH1D* Data_num_mass_log = new TH1D("Data_num_mass_log", "Data_num_mass_log", NMBINS, logMbins);
    TH1D* Data_num_mass_log_BB = new TH1D("Data_num_mass_log_BB", "Data_num_mass_log_BB", NMBINS, logMbins);
    TH1D* Data_num_mass_log_BE = new TH1D("Data_num_mass_log_BE", "Data_num_mass_log_BE", NMBINS, logMbins);
    
    TH1D* Data_den_mass_log = new TH1D("Data_den_mass_log", "Data_den_mass_log", NMBINS, logMbins);
    TH1D* Data_den_mass_log_BB = new TH1D("Data_den_mass_log_BB", "Data_den_mass_log_BB", NMBINS, logMbins);
    TH1D* Data_den_mass_log_BE = new TH1D("Data_den_mass_log_BE", "Data_den_mass_log_BE", NMBINS, logMbins);
    
    TH1D* num_Data_transverse_mass = new TH1D("num_Data_transverse_mass", "num_Data_transverse_mass", 50, 0, 500);
    TH1D* num_Data_dz = new TH1D("num_Data_dz", "num_Data_dz", 50, 0, 10);
    TH1D* num_Data_met_pt = new TH1D("num_Data_met_pt", "num_Data_met_pt", 100, 0, 1000);
    TH1D* num_Data_met_phi = new TH1D("num_Data_met_phi", "num_Data_met_phi", 20, -3, 3);
    
    TH1D* den_Data_transverse_mass = new TH1D("den_Data_transverse_mass", "den_Data_transverse_mass", 50, 0, 500);
    TH1D* den_Data_dz = new TH1D("den_Data_dz", "den_Data_dz", 50, 0, 10);
    TH1D* den_Data_met_pt = new TH1D("den_Data_met_pt", "den_Data_met_pt", 100, 0, 1000);
    TH1D* den_Data_met_phi = new TH1D("den_Data_met_phi", "den_Data_met_phi", 20, -3, 3);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    TH1D* total_num_pt_Barrel = new TH1D("total_num_pt_Barrel", "total_num_pt_Barrel", binnum_pt, pt_bins);
    total_num_pt_Barrel->Sumw2();
    TH1D* total_den_pt_Barrel = new TH1D("total_den_pt_Barrel", "total_den_pt_Barrel", binnum_pt, pt_bins);
    total_den_pt_Barrel->Sumw2();
    TH1D* total_num_pt_Endcap = new TH1D("total_num_pt_Endcap", "total_num_pt_Endcap", binnum_pt, pt_bins);
    total_num_pt_Endcap->Sumw2();
    TH1D* total_den_pt_Endcap = new TH1D("total_den_pt_Endcap", "total_den_pt_Endcap", binnum_pt, pt_bins);
    total_den_pt_Endcap->Sumw2();
    TH1D* FR_pt_Endcap = new TH1D("FR_pt_Endcap", "FR_pt_Endcap", binnum_pt, pt_bins);
    FR_pt_Endcap->Sumw2();
    TH1D* FR_pt_Barrel = new TH1D("FR_pt_Barrel", "FR_pt_Barrel", binnum_pt, pt_bins);
    FR_pt_Barrel->Sumw2();
    
    
    
    
    //strat of looping over MC samples
    for(int j=0; j<mc; j++){
        
        
        std::cout<<"opening.. "<<MC_samples[j]<<" --- "<<" No of evnets: "<<events[j]<<" Cross section: "<<sigma[j]<<std::endl;
        
        TChain *treeMC = new TChain("SimpleNtupler/t");
        
        treeMC->Add("/eos/user/k/kaliyana/2018_MC/FileBased/"+ MC_samples[j] +"/"+ MC_samples[j] +".root");
        
        treeMC->SetBranchAddress("genWeight",&genWeight);
        treeMC->SetBranchAddress("event",&event);
        treeMC->SetBranchAddress("run",&run);
        treeMC->SetBranchAddress("lumi",&lumi);
        treeMC->SetBranchAddress("dil_mass",&dil_mass);
        treeMC->SetBranchAddress("dil_pt",&dil_pt);
        treeMC->SetBranchAddress("dil_eta",&dil_eta);
        treeMC->SetBranchAddress("dil_rap",&dil_rap);
        treeMC->SetBranchAddress("dil_phi",&dil_phi);
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
        
        Long64_t ne = treeMC->GetEntries();
        dimu_MC[j][0] = ne;
        
        std::cout<<"START: "<<ne<<" --For control region"<<std::endl;
        
        
        //looping over entries for control region (denominator)
        for(int p=0; p<ne ;p++){
            
            if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            treeMC->GetEntry(p); //takes pth event
            
            if(j==15 && run==1 && lumi==24733 && event==148397835) continue;
            
           // if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
           // prev_event = event;
            
            
            if (!(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50)) continue; //dimuon evnts should pass the trigger
            
            // c_event=event;
            
            
            
            //Consider collisions where you have only two muons. Events with dimuons having the same event number are rejected.
            /*   for(int i=p+1; i<ne; i++){
             
             treeMC->GetEntry(i); //takes the next event
             
             if(c_event==event){
             dimu_MC[j][1] ++;
             p=i+1;
             //std::cout<<p<<" "<<i<<" "<<c_event<<"   "<<event<<std::endl;
             continue;
             }
             else{
             treeMC->GetEntry(p);
             c_event=event;
             treeMC->GetEntry(p+1);
             if(c_event==event){
             p++;
             continue;
             }
             else{
             treeMC->GetEntry(p);
             break;
             }
             }
             }
             
             
             
             
             // vetoing dimuon events with two high pt muons passing the ID and isolation cuts
             if(
             GoodVtx &&
             (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4) &&
             (lep_pt[0]>53. && lep_pt[1]>53.) &&
             (lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1) &&
             (lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1) &&
             (fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2) &&
             (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || ((lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) && lep_numberOfMatchedRPCLayers[0] > 2)) &&
             (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || ((lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) && lep_numberOfMatchedRPCLayers[1] > 2)) &&
             ((lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)) &&
             (lep_glb_numberOfValidPixelHits[0]>0 && lep_glb_numberOfValidPixelHits[1]>0) &&
             (lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5) &&
             (lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10) &&
             (lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3) &&
             cos_angle>-0.9998 &&
             lep_id[0]*lep_id[1]<0 &&
             (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) &&
             //fabs(lep_tk_dz[0]) < 1.0 && fabs(lep_tk_dz[1]) < 1.0 &&
             vertex_chi2 < 20
             ){
             highPt_MC[j]++;
             continue;
             
             }
             */
            weight[j] *= genWeight;
            
            
            
            
            
            if(lep_pt[0]>lep_pt[1]){
                leading_pt = lep_pt[0];
                leading_phi = lep_phi[0];
            }
            else{
                leading_pt = lep_pt[1];
                leading_phi = lep_phi[1];
            }
            delta_phi = std::abs(met_phi-leading_phi);
            if(delta_phi>PI) delta_phi=float(2*PI)-delta_phi;
            mt = sqrt(2*leading_pt*met_pt*(1-cos(delta_phi)));
            
            
            
            
            //if(mt>35)  continue; // against W+jets
            
            // if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
            // prev_event = event;
            
            //std:: cout<< mt << std::endl;
            
            // std::cout<<"lep_eta[n]"<<"   "<<"lep_isTrackerMuon[n]"<<"  "<<"lep_isGlobalMuon[n]"<<"   "<<"lep_tk_dz[n]"<<"   "<<"lep_glb_numberOfValidPixelHits[n]"<<"  "<<"lep_glb_numberOfValidTrackerLayers[n]"<<"    "<<"lep_triggerMatchPt[n]"<<std::endl;
            
            //For muons in the dimuon event
            for(int n=0; n<2; n++){
                //check conditions for the control region (denominator)
                if(
                   //fabs(lep_eta[n])<2.4 &&
                   //lep_pt[n]>=53. &&
                   lep_isTrackerMuon[n]==1 &&
                   lep_isGlobalMuon[n]==1 &&
                   //fabs(lep_tk_dz[n]) < 1.0 &&
                   fabs(lep_dB[n]) < 0.2 &&
                   lep_glb_numberOfValidPixelHits[n] > 0.0 &&
                   lep_glb_numberOfValidTrackerLayers[n] > 5.0
                   //lep_id[0]*lep_id[1]<0 &&
                   //lep_triggerMatchPt[n] > 50.0
                   ){
                    
                    den_MC[j][0]++;
                    
                    // std::cout<<lep_eta[n]<<"   "<<lep_isTrackerMuon[n]<<"  "<<lep_isGlobalMuon[n]<<"   "<<lep_tk_dz[n]<<"   "<<lep_glb_numberOfValidPixelHits[n]<<"  "<<lep_glb_numberOfValidTrackerLayers[n]<<"    "<<lep_triggerMatchPt[n]<<std::endl;
                    
                    if(j==15 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //if((run==1 && lumi==14910 && event==37272871) || (run==1 && lumi==10525 && event==26311724) || (run==1 && lumi==6454 && event==16134369) || (run==1 && lumi==20725 && event==51812270)) continue;
                        
                        //den_Wjets_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        Wjets_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            Wjets_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            Wjets_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            Wjets_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            Wjets_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                        
                        
                        if(lep_pt[n]>400 && lep_pt[n]<1000){
                            std::cout<<"den"<<" "<<run<<"   "<<lumi<<"  "<<event<<" "<<lep_pt[n]<<" "<<lep_eta[n]<<" "<<weight[j]<<std::endl;
                        }
                        
                    }
                    
                    
                    
                    else if(j==16 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //den_Zjets_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        Zjets_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            Zjets_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            Zjets_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            Zjets_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            Zjets_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                    }
                    
                    else if(j==17 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        
                        gM = gen_dil_mass;
                        //kFactor =0.994078695151 + gM*2.64819793287e-05 - gM*gM*3.73996461024e-08 - gM*gM*gM*1.11452866827e-11;
                        kFactor = 0.991403 + 3.05593e-05 * gM - 2.21967e-07 * pow(gM,2) + 6.63658e-11 * pow(gM,3);
                        kFactor_BB = 0.990973 + 6.17124e-06 * gM - 3.31244e-07 * pow(gM,2) + 1.2125e-10 * pow(gM,3);
                        kFactor_BE = 0.990038 + 5.44269e-05 * gM - 2.43311e-07 * pow(gM,2) + 5.9748e-11 * pow(gM,3);
                        
                        
                        //den_ttbar_dz->Fill(lep_tk_dz[n], kFactor*weight[j]*Z_peak);
                        ttbar_den_eta->Fill(lep_eta[n], kFactor*weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            ttbar_den_pt_Barrel->Fill(lep_pt[n], kFactor_BB*weight[j]*Z_peak_BB);
                            ttbar_den_eta_Barrel->Fill(lep_eta[n], kFactor_BB*weight[j]*Z_peak_BB);
                        }
                        else {
                            ttbar_den_pt_Endcap->Fill(lep_pt[n], kFactor_BE*weight[j]*Z_peak_BE);
                            ttbar_den_eta_Endcap->Fill(lep_eta[n], kFactor_BE*weight[j]*Z_peak_BE);
                        }
                    }
                    
                    
                    
                    else if(j==18 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //den_WW_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        WW_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            WW_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            WW_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            WW_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            WW_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                    }
                    
                    else if(j==19 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //den_WZ_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        WZ_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            WZ_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            WZ_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            WZ_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            WZ_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                    }
                    
                    else if(j==20 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //den_ZZ_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        ZZ_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            ZZ_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            ZZ_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            ZZ_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            ZZ_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                    }
                    
                    else if(j==21 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //den_tW_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        tW_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            tW_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            tW_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            tW_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            tW_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                    }
                    
                    else if(j==22 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        //den_Wantitop_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        Wantitop_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            Wantitop_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            Wantitop_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            Wantitop_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            Wantitop_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                    }
                    else if((j==23 || j==24 || j==25 || j==26 || j==27 || j==28 || j==29 || j==30 || j==31 || j==32)/* && mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                        
                        
                        
                        gM = gen_dil_mass - 400.0;
                        kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
                        kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
                        kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
                        
                        
                        //den_dy_dz->Fill(lep_tk_dz[n], kFactor*weight[j]*Z_peak);
                        dy_den_eta->Fill(lep_eta[n], kFactor*weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            dy_den_pt_Barrel->Fill(lep_pt[n], kFactor_BB*weight[j]*Z_peak_BB);
                            dy_den_eta_Barrel->Fill(lep_eta[n], kFactor_BB*weight[j]*Z_peak_BB);
                        }
                        else {
                            dy_den_pt_Endcap->Fill(lep_pt[n], kFactor_BE*weight[j]*Z_peak_BE);
                            dy_den_eta_Endcap->Fill(lep_eta[n], kFactor_BE*weight[j]*Z_peak_BE);
                        }
                    }
                    
                    else{
                        //if(mt<35 && fabs(lep_tk_dz[n]) < 1.0){
                        
                        //den_qcd_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                        qcd_den_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            qcd_den_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                            qcd_den_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                        }
                        else {
                            qcd_den_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                            qcd_den_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                        }
                        //}
                    }
                    
                    
                    
                    
                    //dimuon events with two high pt muons passing the ID and isolation cuts for numerator
                    if(
                       //GoodVtx &&
                       fabs(lep_eta[n])<2.4 &&
                       lep_pt[n]>=53. &&
                       lep_isTrackerMuon[n]==1 &&
                       lep_isGlobalMuon[n]==1 &&
                       fabs(lep_dB[n]) < 0.2 &&
                       //(lep_numberOfMatchedStations[n] > 1|| (lep_numberOfMatchedStations[n] == 1 && !(lep_stationMask[n] == 1 || lep_stationMask[n] == 16)) || ((lep_numberOfMatchedStations[n] == 1 && (lep_stationMask[n] == 1 || lep_stationMask[n] == 16)) && lep_numberOfMatchedRPCLayers[n] > 2)) &&
                       (lep_numberOfMatchedStations[n] > 1 || (lep_numberOfMatchedStations[n] == 1 && (lep_expectedNnumberOfMatchedStations[n] < 2 || !(lep_stationMask[n] == 1 || lep_stationMask[n] == 16) || lep_numberOfMatchedRPCLayers[n] > 2))) &&
                       (lep_glb_numberOfValidMuonHits[n]>0 || lep_TuneP_numberOfValidMuonHits[n]>0)  &&
                       lep_glb_numberOfValidPixelHits[n]>0 &&
                       lep_glb_numberOfValidTrackerLayers[n]>5 &&
                       lep_sumPt[n]/lep_tk_pt[n]<0.10 &&
                       lep_pt_err[n]/lep_pt[n]<0.3 &&
                       //cos_angle>-0.9998 &&
                       //lep_id[0]*lep_id[1]<0 &&
                       lep_triggerMatchPt[n]>50
                       //vertex_chi2 < 20
                       //fabs(lep_tk_dz[n]) < 1.0
                       ){
                        
                        
                        num_MC[j][0]++;
                        
                        if(j==15 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //if((run==1 && lumi==14910 && event==37272871) || (run==1 && lumi==10525 && event==26311724) || (run==1 && lumi==6454 && event==16134369) ) continue;
                            
                            //num_Wjets_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            Wjets_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                Wjets_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                Wjets_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                Wjets_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                Wjets_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                            
                            if(lep_pt[n]>400 && lep_pt[n]<1000){
                                std::cout<<"num"<<" "<<run<<"   "<<lumi<<"  "<<event<<" "<<lep_pt[n]<<" "<<lep_eta[n]<<" "<<weight[j]<<std::endl;
                            }
                            
                        }
                        
                        
                        
                        else if(j==16 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //num_Zjets_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            Zjets_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                Zjets_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                Zjets_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                Zjets_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                Zjets_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                        }
                        
                        else if(j==17 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            
                            gM = gen_dil_mass;
                            //kFactor =0.994078695151 + gM*2.64819793287e-05 - gM*gM*3.73996461024e-08 - gM*gM*gM*1.11452866827e-11;
                            kFactor = 0.991403 + 3.05593e-05 * gM - 2.21967e-07 * pow(gM,2) + 6.63658e-11 * pow(gM,3);
                            kFactor_BB = 0.990973 + 6.17124e-06 * gM - 3.31244e-07 * pow(gM,2) + 1.2125e-10 * pow(gM,3);
                            kFactor_BE = 0.990038 + 5.44269e-05 * gM - 2.43311e-07 * pow(gM,2) + 5.9748e-11 * pow(gM,3);
                            
                            
                            //num_ttbar_dz->Fill(lep_tk_dz[n], kFactor*weight[j]*Z_peak);
                            ttbar_num_eta->Fill(lep_eta[n], kFactor*weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                ttbar_num_pt_Barrel->Fill(lep_pt[n], kFactor_BB*weight[j]*Z_peak_BB);
                                ttbar_num_eta_Barrel->Fill(lep_eta[n], kFactor_BB*weight[j]*Z_peak_BB);
                            }
                            else {
                                ttbar_num_pt_Endcap->Fill(lep_pt[n], kFactor_BE*weight[j]*Z_peak_BE);
                                ttbar_num_eta_Endcap->Fill(lep_eta[n], kFactor_BE*weight[j]*Z_peak_BE);
                            }
                        }
                        
                        
                        
                        else if(j==18 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //num_WW_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            WW_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                WW_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                WW_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                WW_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                WW_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                        }
                        
                        else if(j==19 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //num_WZ_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            WZ_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                WZ_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                WZ_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                WZ_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                WZ_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                        }
                        
                        else if(j==20 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //num_ZZ_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            ZZ_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                ZZ_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                ZZ_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                ZZ_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                ZZ_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                        }
                        
                        else if(j==21 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //num_tW_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            tW_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                tW_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                tW_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                tW_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                tW_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                        }
                        
                        else if(j==22 /*&& mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            //num_Wantitop_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            Wantitop_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                Wantitop_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                Wantitop_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                Wantitop_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                Wantitop_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                        }
                        else if((j==23 || j==24 || j==25 || j==26 || j==27 || j==28 || j==29 || j==30 || j==31 || j==32)/* && mt<35 && fabs(lep_tk_dz[n]) < 1.0*/){
                            
                            
                            
                            gM = gen_dil_mass - 400.0;
                            kFactor = 1.067 - 0.000112 * gM + 3.176e-08 * pow(gM,2) - 4.068e-12 * pow(gM,3);
                            kFactor_BB = 1.036 - 0.0001441 * gM + 5.058e-8 * pow(gM,2) - 7.581e-12 * pow(gM,3);
                            kFactor_BE = 1.052 - 0.0001471 * gM + 5.903e-8 * pow(gM,2) - 9.037e-12 * pow(gM,3);
                            
                            
                            //num_dy_dz->Fill(lep_tk_dz[n], kFactor*weight[j]*Z_peak);
                            dy_num_eta->Fill(lep_eta[n], kFactor*weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                dy_num_pt_Barrel->Fill(lep_pt[n], kFactor_BB*weight[j]*Z_peak_BB);
                                dy_num_eta_Barrel->Fill(lep_eta[n], kFactor_BB*weight[j]*Z_peak_BB);
                            }
                            else {
                                dy_num_pt_Endcap->Fill(lep_pt[n], kFactor_BE*weight[j]*Z_peak_BE);
                                dy_num_eta_Endcap->Fill(lep_eta[n], kFactor_BE*weight[j]*Z_peak_BE);
                            }
                        }
                        
                        else{
                            //if(mt<35 && fabs(lep_tk_dz[n]) < 1.0){
                            
                            //num_qcd_dz->Fill(lep_tk_dz[n], weight[j]*Z_peak);
                            qcd_num_eta->Fill(lep_eta[n], weight[j]*Z_peak);
                            
                            
                            if(fabs(lep_eta[n]) < 1.2){
                                qcd_num_pt_Barrel->Fill(lep_pt[n], weight[j]*Z_peak_BB);
                                qcd_num_eta_Barrel->Fill(lep_eta[n], weight[j]*Z_peak_BB);
                            }
                            else {
                                qcd_num_pt_Endcap->Fill(lep_pt[n], weight[j]*Z_peak_BE);
                                qcd_num_eta_Endcap->Fill(lep_eta[n], weight[j]*Z_peak_BE);
                            }
                            //}
                        }
                        
                        
                        
                    }//end conditions for numerator
                    
                    
                }//end conditions for the control region(denominator)
                
            }//end looping over muons
            
            weight[j] /= genWeight;
            
        }//end looping over entries
        
        
    }//end looping over MC samples
    
    
    
    
    
    int di =-1;
    prev_event = -88;
    
    //strat of looping over Data
    for(int j=0; j<data; j++) {
        
        
        std::cout<<"opening.. "<<DATA_samples[j]<<std::endl;
        
        
        
        TChain *treeDATA = new TChain("SimpleNtupler/t");
        
        treeDATA->Add("/eos/user/k/kaliyana/2018_data/"+ samp[j] +"/"+ DATA_samples[j] +".root");
        
        
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
        treeDATA->SetBranchAddress("GoodVtx",&GoodVtx);
        treeDATA->SetBranchAddress("lep_stationMask",&lep_stationMask);
        treeDATA->SetBranchAddress("lep_numberOfMatchedRPCLayers",lep_numberOfMatchedRPCLayers);
        //treeDATA->SetBranchAddress("gen_dil_mass", &gen_dil_mass);
        //treeDATA->SetBranchAddress("gen_lep_qOverPt", gen_lep_qOverPt);
        //treeDATA->SetBranchAddress("gen_lep_eta", gen_lep_eta);
        //treeDATA->SetBranchAddress("gen_lep_pt", gen_lep_pt);
        
        treeDATA->SetBranchAddress("met_pt",&met_pt);
        treeDATA->SetBranchAddress("met_phi",&met_phi);
        treeDATA->SetBranchAddress("lep_tuneP_pt",lep_tuneP_pt);
        treeDATA->SetBranchAddress("lep_tk_dz",lep_tk_dz);
        
        Long64_t ne = treeDATA->GetEntries();
        dimu_DATA[j][0] = ne;
        
        std::cout<<"START: "<<"<<ne<<"<<" --For control region"<<std::endl;
        
        
        
        //looping over entries
        for ( int p=0; p<ne ;p++) {
            if(p % 100000 == 0) std::cout<<p<<std::endl;
            
            
            
            treeDATA->GetEntry(p); //takes pth event
            
            //if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
            //prev_event = event;
            
            if (!(lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50)) continue; //dimuon evnts should pass the trigger
            
            
            
            //  c_event=event;
            
            
            //Consider collisions where you have only two muons. Events with dimuons having the same event number are rejected.
            /*     for(int i=p+1; i<ne; i++){
             
             treeDATA->GetEntry(i); //takes the next event
             
             if(c_event==event){
             dimu_DATA[j][1] ++;
             p=i+1;
             //std::cout<<p<<" "<<i<<" "<<c_event<<"   "<<event<<std::endl;
             continue;
             }
             else{
             treeDATA->GetEntry(p);
             c_event=event;
             treeDATA->GetEntry(p+1);
             if(c_event==event){
             p++;
             continue;
             }
             else{
             treeDATA->GetEntry(p);
             break;
             }
             }
             }
             
             
             
             // vetoing dimuon events with two high pt muons passing the ID and isolation cuts
             if(
             GoodVtx &&
             (fabs(lep_eta[0])<2.4 && fabs(lep_eta[1])<2.4) &&
             (lep_pt[0]>53. && lep_pt[1]>53.) &&
             (lep_isTrackerMuon[0]==1 && lep_isTrackerMuon[1]==1) &&
             (lep_isGlobalMuon[0]==1 && lep_isGlobalMuon[1]==1) &&
             (fabs(lep_dB[0]) < 0.2 && fabs(lep_dB[1]) < 0.2) &&
             (lep_numberOfMatchedStations[0] > 1 || (lep_numberOfMatchedStations[0] == 1 && !(lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) || ((lep_numberOfMatchedStations[0] == 1 && (lep_stationMask[0] == 1 || lep_stationMask[0] == 16)) && lep_numberOfMatchedRPCLayers[0] > 2)) &&
             (lep_numberOfMatchedStations[1] > 1 || (lep_numberOfMatchedStations[1] == 1 && !(lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) || ((lep_numberOfMatchedStations[1] == 1 && (lep_stationMask[1] == 1 || lep_stationMask[1] == 16)) && lep_numberOfMatchedRPCLayers[1] > 2)) &&
             ((lep_glb_numberOfValidMuonHits[0]>0 || lep_TuneP_numberOfValidMuonHits[0]>0) && (lep_glb_numberOfValidMuonHits[1]>0 || lep_TuneP_numberOfValidMuonHits[1]>0)) &&
             (lep_glb_numberOfValidPixelHits[0]>0 && lep_glb_numberOfValidPixelHits[1]>0) &&
             (lep_glb_numberOfValidTrackerLayers[0]>5 && lep_glb_numberOfValidTrackerLayers[1]>5) &&
             (lep_sumPt[0]/lep_tk_pt[0]<0.10 && lep_sumPt[1]/lep_tk_pt[1]<0.10) &&
             (lep_pt_err[0]/lep_pt[0]<0.3 && lep_pt_err[1]/lep_pt[1]<0.3) &&
             cos_angle>-0.9998 &&
             lep_id[0]*lep_id[1]<0 &&
             (lep_triggerMatchPt[0]>50 || lep_triggerMatchPt[1]>50) &&
             //fabs(lep_tk_dz[0]) < 1.0 && fabs(lep_tk_dz[1]) < 1.0 &&
             vertex_chi2 < 20
             ){
             highPt_Data[j]++;
             
             //For the table
             
             if(dil_mass>1000){
             
             di++;
             data_kinemat[di][0] = (float)di+1;
             data_kinemat[di][1] = dil_mass;
             data_kinemat[di][14] = met_pt;
             data_kinemat[di][15] = met_phi;
             data_kinemat[di][16] = (float)run;
             data_kinemat[di][17] = (float)lumi;
             data_kinemat[di][18] = (float)event;
             
             if(fabs(lep_eta[0]) < 1.2 && fabs(lep_eta[1]) < 1.2){
             data_kinemat[di][14] = 1.0 ;
             }
             else data_kinemat[di][14] = 2.0;
             
             if(lep_id[0]<0 && lep_id[1]>0){
             
             data_kinemat[di][2] = lep_id[0];
             data_kinemat[di][3] = lep_id[1];
             
             data_kinemat[di][4] = lep_tuneP_pt[0];
             data_kinemat[di][5] = lep_tuneP_pt[1];
             
             data_kinemat[di][6] = lep_tk_pt[0];
             data_kinemat[di][7] = lep_tk_pt[1];
             
             data_kinemat[di][8] = lep_tuneP_pt[0]/lep_tk_pt[0];
             data_kinemat[di][9] = lep_tuneP_pt[1]/lep_tk_pt[1];
             
             data_kinemat[di][10] = lep_eta[0];
             data_kinemat[di][11] = lep_phi[0];
             
             data_kinemat[di][12] = lep_eta[1];
             data_kinemat[di][13] = lep_phi[1];
             
             
             }
             
             else{
             
             data_kinemat[di][2] = lep_id[1];
             data_kinemat[di][3] = lep_id[0];
             
             data_kinemat[di][4] = lep_tuneP_pt[1];
             data_kinemat[di][5] = lep_tuneP_pt[0];
             
             data_kinemat[di][6] = lep_tk_pt[1];
             data_kinemat[di][7] = lep_tk_pt[0];
             
             data_kinemat[di][8] = lep_tuneP_pt[1]/lep_tk_pt[1];
             data_kinemat[di][9] = lep_tuneP_pt[0]/lep_tk_pt[0];
             
             data_kinemat[di][10] = lep_eta[1];
             data_kinemat[di][11] = lep_phi[1];
             
             data_kinemat[di][12] = lep_eta[0];
             data_kinemat[di][13] = lep_phi[0];
             }
             } //end for table
             
             
             continue;
             
             }
             */
            if(lep_pt[0]>lep_pt[1]){
                leading_pt = lep_pt[0];
                leading_phi = lep_phi[0];
            }
            else{
                leading_pt = lep_pt[1];
                leading_phi = lep_phi[1];
            }
            delta_phi = std::abs(met_phi-leading_phi);
            if(delta_phi>PI) delta_phi=float(2*PI)-delta_phi;
            mt = sqrt(2*leading_pt*met_pt*(1-cos(delta_phi)));
            
            //if(mt>35)  continue; // against W+jets
            
            //if(prev_event == event) continue; //to remove double event; take the first --> they are already sorted in pt
           // prev_event = event;
            
            //std:: cout<< mt << std::endl;
            
            //For muons in the dimuon event
            
            for(int n=0; n<2; n++){
                //check conditions for the control region (denominator)
                if(
                   //fabs(lep_eta[n]) < 2.4 &&
                   //lep_pt[n]>=53. &&
                   lep_isTrackerMuon[n]==1 &&
                   lep_isGlobalMuon[n]==1 &&
                   //fabs(lep_tk_dz[n]) < 1.0 &&
                   fabs(lep_dB[n]) < 0.2 &&
                   lep_glb_numberOfValidPixelHits[n] > 0.0 &&
                   lep_glb_numberOfValidTrackerLayers[n] > 5.0
                   //lep_id[0]*lep_id[1]<0 &&
                   //lep_triggerMatchPt[n] > 50.0
                   ){
                    
                    den_DATA[j][0]++;
                    
                    //if(mt<35 && fabs(lep_tk_dz[n]) < 1.0){
                    //den_Data_dz->Fill(lep_tk_dz[n]);
                    Data_den_eta->Fill(lep_eta[n]);
                    
                    
                    if(fabs(lep_eta[n]) < 1.2){
                        Data_den_pt_Barrel->Fill(lep_pt[n]);
                        Data_den_eta_Barrel->Fill(lep_eta[n]);
                    }
                    else {
                        Data_den_pt_Endcap->Fill(lep_pt[n]);
                        Data_den_eta_Endcap->Fill(lep_eta[n]);
                    }
                    //}
                    
                    
                    //dimuon events with two high pt muons passing the ID and isolation cuts for numerator
                    if(
                       //GoodVtx &&
                       fabs(lep_eta[n])<2.4 &&
                       lep_pt[n]>=53. &&
                       lep_isTrackerMuon[n]==1 &&
                       lep_isGlobalMuon[n]==1 &&
                       fabs(lep_dB[n]) < 0.2 &&
                       //(lep_numberOfMatchedStations[n] > 1 || (lep_numberOfMatchedStations[n] == 1 && !(lep_stationMask[n] == 1 || lep_stationMask[n] == 16)) || ((lep_numberOfMatchedStations[n] == 1 && (lep_stationMask[n] == 1 || lep_stationMask[n] == 16)) && lep_numberOfMatchedRPCLayers[n] > 2)) &&
                       (lep_numberOfMatchedStations[n] > 1 || (lep_numberOfMatchedStations[n] == 1 && (lep_expectedNnumberOfMatchedStations[n] < 2 || !(lep_stationMask[n] == 1 || lep_stationMask[n] == 16) || lep_numberOfMatchedRPCLayers[n] > 2))) &&
                       (lep_glb_numberOfValidMuonHits[n]>0 || lep_TuneP_numberOfValidMuonHits[n]>0)  &&
                       lep_glb_numberOfValidPixelHits[n]>0 &&
                       lep_glb_numberOfValidTrackerLayers[n]>5 &&
                       lep_sumPt[n]/lep_tk_pt[n]<0.10 &&
                       lep_pt_err[n]/lep_pt[n]<0.3 &&
                       //cos_angle>-0.9998 &&
                       //lep_id[0]*lep_id[1]<0 &&
                       lep_triggerMatchPt[n]>50
                       //vertex_chi2 < 20
                       //fabs(lep_tk_dz[n]) < 1.0
                       ){
                        
                        //std::cout<<"Kal"<<std::endl;
                        
                        num_DATA[j][0]++;
                        
                        //std::cout<<"Kal"<<std::endl;
                        //num_Data_dz->Fill(lep_tk_dz[n]);
                        Data_num_eta->Fill(lep_eta[n]);
                        
                        
                        if(fabs(lep_eta[n]) < 1.2){
                            Data_num_pt_Barrel->Fill(lep_pt[n]);
                            Data_num_eta_Barrel->Fill(lep_eta[n]);
                        }
                        else {
                            Data_num_pt_Endcap->Fill(lep_pt[n]);
                            Data_num_eta_Endcap->Fill(lep_eta[n]);
                        }
                        
                        
                    }// end conditions for numerator
                    
                }//end conditions for the control region(denominator)
                
            }//end looping over muons
            
            
        }//end looping over entries
        
        
    }//end looping over Data samples
    
    
    
    
    file1->Write();
    
    
    std::cout<<"Number"<<"  "<<"MC sample"<<"   "<<"Total entries"<<"    "<<"Entries in the same event"<<"      "<<"Ratio"<<"   "<<"# of High pt dimuons"<<"    "<<"Entries considered for the cuts"<<std::endl;
    for(int i=0; i<23; i++){
        std::cout<<i<<"    "<<MC_samples[i]<<"  "<<dimu_MC[i][0]<<"    "<<dimu_MC[i][1]<<"    "<<(float)dimu_MC[i][1]/dimu_MC[i][0]<<"  "<<highPt_MC[i]<<"  "<<dimu_MC[i][0]-dimu_MC[i][1]-highPt_MC[i]<<std::endl;
        std::cout<<num_MC[i][0]<<std::endl;
        std::cout<<den_MC[i][0]<<std::endl;
        std::cout<<"========================================"<<std::endl;
    }
    
    std::cout<<"Number"<<"  "<<"Data sample"<<"   "<<"Total entries"<<"    "<<"Entries in the same event"<<"    "<<"Ratio"<<"   "<<"# of High pt dimuons"<<"    "<<"Entries considered for the cuts"<<std::endl;
    for(int i=0; i<4; i++){
        std::cout<<i<<"    "<<DATA_samples[i]<<"  "<<dimu_DATA[i][0]<<"    "<<dimu_DATA[i][1]<<"    "<<(float)dimu_DATA[i][1]/dimu_DATA[i][0]<<"  "<<highPt_Data[i]<<"  "<<dimu_DATA[i][0]-dimu_DATA[i][1]-highPt_Data[i]<<std::endl;
        std::cout<<num_DATA[i][0]<<std::endl;
        std::cout<<den_DATA[i][0]<<std::endl;
        std::cout<<"========================================"<<std::endl;
    }
    
    
    
    
    
    
    
    
    
    DrawPlot(qcd_num_pt_Barrel, Wjets_num_pt_Barrel, Zjets_num_pt_Barrel, ttbar_num_pt_Barrel, WW_num_pt_Barrel, WZ_num_pt_Barrel, ZZ_num_pt_Barrel, tW_num_pt_Barrel, Wantitop_num_pt_Barrel,  dy_num_pt_Barrel, Data_num_pt_Barrel, num_pt_Barrel, "pT(#mu)[GeV]", "Barrel_numerator", 1);
    
    DrawPlot(qcd_num_pt_Endcap, Wjets_num_pt_Endcap, Zjets_num_pt_Endcap, ttbar_num_pt_Endcap, WW_num_pt_Endcap, WZ_num_pt_Endcap, ZZ_num_pt_Endcap, tW_num_pt_Endcap, Wantitop_num_pt_Endcap, dy_num_pt_Endcap, Data_num_pt_Endcap, num_pt_Endcap, "pT(#mu)[GeV]", "Endcap_numerator", 1);
    
    DrawPlot(qcd_den_pt_Barrel, Wjets_den_pt_Barrel, Zjets_den_pt_Barrel, ttbar_den_pt_Barrel, WW_den_pt_Barrel, WZ_den_pt_Barrel, ZZ_den_pt_Barrel, tW_den_pt_Barrel, Wantitop_den_pt_Barrel,  dy_den_pt_Barrel, Data_den_pt_Barrel, den_pt_Barrel, "pT(#mu)[GeV]", "Barrel_denominator", 1);
    
    DrawPlot(qcd_den_pt_Endcap, Wjets_den_pt_Endcap, Zjets_den_pt_Endcap, ttbar_den_pt_Endcap, WW_den_pt_Endcap, WZ_den_pt_Endcap, ZZ_den_pt_Endcap, tW_den_pt_Endcap, Wantitop_den_pt_Endcap, dy_den_pt_Endcap, Data_den_pt_Endcap, den_pt_Endcap, "pT(#mu)[GeV]", "Endcap_denominator", 1);
    
    
    
    //DrawPlot(num_qcd_transverse_mass, num_Wjets_transverse_mass, num_Zjets_transverse_mass, num_ttbar_transverse_mass, num_WW_transverse_mass, num_WZ_transverse_mass, num_ZZ_transverse_mass, num_tW_transverse_mass, num_Wantitop_transverse_mass, num_dy_transverse_mass, num_Data_transverse_mass, num_transverse_mass, "mT(#mu+MET) [GeV]", "num_Transeverse_mass", 0);
    
    //DrawPlot(den_qcd_transverse_mass, den_Wjets_transverse_mass, den_Zjets_transverse_mass, den_ttbar_transverse_mass, den_WW_transverse_mass, den_WZ_transverse_mass, den_ZZ_transverse_mass, den_tW_transverse_mass, den_Wantitop_transverse_mass, den_dy_transverse_mass, den_Data_transverse_mass, den_transverse_mass, "mT(#mu+MET) [GeV]", "den_Transeverse_mass", 0);
    
    //DrawPlot(num_qcd_met_phi, num_Wjets_met_phi, num_Zjets_met_phi, num_ttbar_met_phi, num_WW_met_phi, num_WZ_met_phi, num_ZZ_met_phi, num_tW_met_phi, num_Wantitop_met_phi, num_dy_met_phi, num_Data_met_phi, num_MET_phi, "MET_phi [rad]", "num_met_phi", 0);
    
    //DrawPlot(den_qcd_met_phi, den_Wjets_met_phi, den_Zjets_met_phi, den_ttbar_met_phi, den_WW_met_phi, den_WZ_met_phi, den_ZZ_met_phi, den_tW_met_phi, den_Wantitop_met_phi, den_dy_met_phi, den_Data_met_phi, den_MET_phi, "MET_phi [rad]", "den_met_phi", 0);
    
    //DrawPlot(num_qcd_met_pt, num_Wjets_met_pt, num_Zjets_met_pt, num_ttbar_met_pt, num_WW_met_pt, num_WZ_met_pt, num_ZZ_met_pt, num_tW_met_pt, num_Wantitop_met_pt, num_dy_met_pt, num_Data_met_pt, num_MET_pt, "MET_pT [GeV]", "num_met_pt", 0);
    
    //DrawPlot(den_qcd_met_pt, den_Wjets_met_pt, den_Zjets_met_pt, den_ttbar_met_pt, den_WW_met_pt, den_WZ_met_pt, den_ZZ_met_pt, den_tW_met_pt, den_Wantitop_met_pt, den_dy_met_pt, den_Data_met_pt, den_MET_pt, "MET_pT [GeV]", "den_met_pt", 0);
    
    //DrawPlot(num_qcd_dz, num_Wjets_dz, num_Zjets_dz, num_ttbar_dz, num_WW_dz, num_WZ_dz, num_ZZ_dz, num_tW_dz, num_Wantitop_dz, num_dy_dz, num_Data_dz, num_dz, "dz [cm]", "dz", 0);
    
    // DrawPlot(den_qcd_dz, den_Wjets_dz, den_Zjets_dz, den_ttbar_dz, den_WW_dz, den_WZ_dz, den_ZZ_dz, den_tW_dz, den_Wantitop_dz, den_dy_dz, den_Data_dz, den_dz, "dz [cm]", "dz", 0);
    
    // DrawPlot(qcd_num_eta_Barrel, Wjets_num_eta_Barrel, Zjets_num_eta_Barrel, ttbar_num_eta_Barrel, WW_num_eta_Barrel, WZ_num_eta_Barrel, ZZ_num_eta_Barrel, tW_num_eta_Barrel, Wantitop_num_eta_Barrel,  dy_num_eta_Barrel, Data_num_eta_Barrel, num_eta_Barrel, "eta_Barrel(#mu)", "Barrel_numerator_eta", 0);
    
    //  DrawPlot(qcd_num_eta_Endcap, Wjets_num_eta_Endcap, Zjets_num_eta_Endcap, ttbar_num_eta_Endcap, WW_num_eta_Endcap, WZ_num_eta_Endcap, ZZ_num_eta_Endcap, tW_num_eta_Endcap, Wantitop_num_eta_Endcap, dy_num_eta_Endcap, Data_num_eta_Endcap, num_eta_Endcap, "eta_Endcap#mu)", "Endcap_numerator_eta", 0);
    
    //   DrawPlot(qcd_den_eta_Barrel, Wjets_den_eta_Barrel, Zjets_den_eta_Barrel, ttbar_den_eta_Barrel, WW_den_eta_Barrel, WZ_den_eta_Barrel, ZZ_den_eta_Barrel, tW_den_eta_Barrel, Wantitop_den_eta_Barrel,  dy_den_eta_Barrel, Data_den_eta_Barrel, den_eta_Barrel, "eta_Barrel(#mu)", "Barrel_denominator_eta", 0);
    
    //  DrawPlot(qcd_den_eta_Endcap, Wjets_den_eta_Endcap, Zjets_den_eta_Endcap, ttbar_den_eta_Endcap, WW_den_eta_Endcap, WZ_den_eta_Endcap, ZZ_den_eta_Endcap, tW_den_eta_Endcap, Wantitop_den_eta_Endcap, dy_den_eta_Endcap, Data_den_eta_Endcap, den_eta_Endcap, "eta_Endcap(#mu)", "Endcap_denominator_eta", 0);
    
    
    //  DrawPlot(qcd_den_eta, Wjets_den_eta, Zjets_den_eta, ttbar_den_eta, WW_den_eta, WZ_den_eta, ZZ_den_eta, tW_den_eta, Wantitop_den_eta, dy_den_eta, Data_den_eta, den_eta, "eta_den(#mu)", "den_eta", 0);
    
    //  DrawPlot(qcd_num_eta, Wjets_num_eta, Zjets_num_eta, ttbar_num_eta, WW_num_eta, WZ_num_eta, ZZ_num_eta, tW_num_eta, Wantitop_num_eta, dy_num_eta, Data_num_eta, num_eta, "eta_num(#mu)", "num_eta", 0);
    
    //DrawPlot(qcd_num_mass_log, Wjets_num_mass_log, Zjets_num_mass_log, ttbar_num_mass_log, WW_num_mass_log, WZ_num_mass_log, ZZ_num_mass_log, tW_num_mass_log, Wantitop_num_mass_log, dy_num_mass_log, Data_num_mass_log, num_mass_log, "m(#mu^{+}#mu^{-})[GeV]", "num_mass_log", 1);
    
    
    //DrawPlot(qcd_den_mass_log, Wjets_den_mass_log, Zjets_den_mass_log, ttbar_den_mass_log, WW_den_mass_log, WZ_den_mass_log, ZZ_den_mass_log, tW_den_mass_log, Wantitop_den_mass_log, dy_den_mass_log, Data_den_mass_log, den_mass_log, "m(#mu^{+}#mu^{-})[GeV]", "den_mass_log", 1);
    
    //DrawPlot(qcd_num_mass_log_BB, Wjets_num_mass_log_BB, Zjets_num_mass_log_BB, ttbar_num_mass_log_BB, WW_num_mass_log_BB, WZ_num_mass_log_BB, ZZ_num_mass_log_BB, tW_num_mass_log_BB, Wantitop_num_mass_log_BB, dy_num_mass_log_BB, Data_num_mass_log_BB, num_mass_log_BB, "m(#mu^{+}#mu^{-})[GeV]", "num_mass_log_BB", 1);
    
    
    // DrawPlot(qcd_den_mass_log_BB, Wjets_den_mass_log_BB, Zjets_den_mass_log_BB, ttbar_den_mass_log_BB, WW_den_mass_log_BB, WZ_den_mass_log_BB, ZZ_den_mass_log_BB, tW_den_mass_log_BB, Wantitop_den_mass_log_BB, dy_den_mass_log_BB, Data_den_mass_log_BB, den_mass_log_BB, "m(#mu^{+}#mu^{-})[GeV]", "den_mass_log_BB", 1);
    
    
    // DrawPlot(qcd_num_mass_log_BE, Wjets_num_mass_log_BE, Zjets_num_mass_log_BE, ttbar_num_mass_log_BE, WW_num_mass_log_BE, WZ_num_mass_log_BE, ZZ_num_mass_log_BE, tW_num_mass_log_BE, Wantitop_num_mass_log_BE, dy_num_mass_log_BE, Data_num_mass_log_BE, num_mass_log_BE, "m(#mu^{+}#mu^{-})[GeV]", "num_mass_log_BE+EE", 1);
    
    
    //DrawPlot(qcd_den_mass_log_BE, Wjets_den_mass_log_BE, Zjets_den_mass_log_BE, ttbar_den_mass_log_BE, WW_den_mass_log_BE, WZ_den_mass_log_BE, ZZ_den_mass_log_BE, tW_den_mass_log_BE, Wantitop_den_mass_log_BE, dy_den_mass_log_BE, Data_den_mass_log_BE, den_mass_log_BE, "m(#mu^{+}#mu^{-})[GeV]", "den_mass_log_BE+EE", 1);
    
    
    
    
    std::cout<<"Barrel"<<std::endl;
    
    RatePlot(Wjets_den_pt_Barrel, Zjets_den_pt_Barrel, ttbar_den_pt_Barrel, WW_den_pt_Barrel,  WZ_den_pt_Barrel, ZZ_den_pt_Barrel, tW_den_pt_Barrel, Wantitop_den_pt_Barrel, dy_den_pt_Barrel,  Data_den_pt_Barrel, Wjets_num_pt_Barrel,  Zjets_num_pt_Barrel, ttbar_num_pt_Barrel,  WW_num_pt_Barrel,  WZ_num_pt_Barrel,  ZZ_num_pt_Barrel, tW_num_pt_Barrel, Wantitop_num_pt_Barrel,  dy_num_pt_Barrel,  qcd_den_pt_Barrel, qcd_num_pt_Barrel, Data_num_pt_Barrel, "Fake Rate Barrel", "Fake_Rate_Barrel");
    
    std::cout<<"Endcap"<<std::endl;
    
    RatePlot(Wjets_den_pt_Endcap, Zjets_den_pt_Endcap,  ttbar_den_pt_Endcap, WW_den_pt_Endcap,  WZ_den_pt_Endcap,  ZZ_den_pt_Endcap,  tW_den_pt_Endcap, Wantitop_den_pt_Endcap, dy_den_pt_Endcap,  Data_den_pt_Endcap,  Wjets_num_pt_Endcap,  Zjets_num_pt_Endcap,  ttbar_num_pt_Endcap,  WW_num_pt_Endcap,  WZ_num_pt_Endcap,  ZZ_num_pt_Endcap,  tW_num_pt_Endcap,  Wantitop_num_pt_Endcap, dy_num_pt_Endcap,  qcd_den_pt_Endcap, qcd_num_pt_Endcap, Data_num_pt_Endcap, "Fake Rate Endcap", "Fake_Rate_Endcap");
    
    
    QCDPlot(qcd_den_pt_Barrel, qcd_num_pt_Barrel, "FR_QCD_Barrel");
    
    QCDPlot(qcd_den_pt_Endcap, qcd_num_pt_Endcap, "FR_QCD_Endcap");
    
    //EWKPlot(Wjets_den_pt_Barrel, Zjets_den_pt_Barrel,  ttbar_den_pt_Barrel, WW_den_pt_Barrel,  WZ_den_pt_Barrel,  ZZ_den_pt_Barrel,  tW_den_pt_Barrel, Wantitop_den_pt_Barrel, dy_den_pt_Barrel,  Wjets_num_pt_Barrel,  Zjets_num_pt_Barrel,  ttbar_num_pt_Barrel,  WW_num_pt_Barrel,  WZ_num_pt_Barrel,  ZZ_num_pt_Barrel,  tW_num_pt_Barrel,  Wantitop_num_pt_Barrel, dy_num_pt_Barrel,  "FR_EWK_BB", "FR_EWK_BB");
    
    
    //EWKPlot(Wjets_den_pt_Endcap, Zjets_den_pt_Endcap,  ttbar_den_pt_Endcap, WW_den_pt_Endcap,  WZ_den_pt_Endcap,  ZZ_den_pt_Endcap,  tW_den_pt_Endcap, Wantitop_den_pt_Endcap, dy_den_pt_Endcap,  Wjets_num_pt_Endcap,  Zjets_num_pt_Endcap,  ttbar_num_pt_Endcap,  WW_num_pt_Endcap,  WZ_num_pt_Endcap,  ZZ_num_pt_Endcap,  tW_num_pt_Endcap,  Wantitop_num_pt_Endcap, dy_num_pt_Endcap,  "FR_EWK_BE+EE", "FR_EWK_BE+EE");
    
    //  calRate(Wjets_den_pt_Barrel, Zjets_den_pt_Barrel, ttbar_den_pt_Barrel, WW_den_pt_Barrel,  WZ_den_pt_Barrel, ZZ_den_pt_Barrel, tW_den_pt_Barrel, Wantitop_den_pt_Barrel, dy_den_pt_Barrel,  Data_den_pt_Barrel, Wjets_num_pt_Barrel,  Zjets_num_pt_Barrel, ttbar_num_pt_Barrel,  WW_num_pt_Barrel,  WZ_num_pt_Barrel,  ZZ_num_pt_Barrel, tW_num_pt_Barrel, Wantitop_num_pt_Barrel,  dy_num_pt_Barrel,  Data_num_pt_Barrel, "Fake Rate Barrel");
    
    //  calRate(Wjets_den_pt_Endcap, Zjets_den_pt_Endcap,  ttbar_den_pt_Endcap, WW_den_pt_Endcap,  WZ_den_pt_Endcap,  ZZ_den_pt_Endcap,  tW_den_pt_Endcap, Wantitop_den_pt_Endcap, dy_den_pt_Endcap,  Data_den_pt_Endcap,  Wjets_num_pt_Endcap,  Zjets_num_pt_Endcap,  ttbar_num_pt_Endcap,  WW_num_pt_Endcap,  WZ_num_pt_Endcap,  ZZ_num_pt_Endcap,  tW_num_pt_Endcap,  Wantitop_num_pt_Endcap, dy_num_pt_Endcap,  Data_num_pt_Endcap, "Fake Rate Endcap");
    
    
    //file1->Write();
    file1->Close();
    
}


void DrawPlot(TH1D* qcd, TH1D* Wjets, TH1D* Zjets, TH1D* ttbar, TH1D* WW, TH1D* WZ, TH1D* ZZ, TH1D* tW, TH1D* Wantitop, TH1D* dy, TH1D* DATA, THStack* DATA_MC, TString title, TString name, bool logx){
    
    
    
    gStyle->SetOptStat(0);
    //gStyle->SetOptFit(1);
    gStyle->SetLegendTextSize(0.02);
    
    
    // float r = 0.25;
    // float epsilon = 0.02;
    
    TCanvas *c3 = new TCanvas("DATA_MC", "DATA_MC", 800, 800);
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
    
    
    WW->SetLineColor(kBlack);
    WW->SetFillColor(kAzure+1);
    DATA_MC->Add(WW, "HIST");
    //WW->Draw("sameHIST");
    c3->Update();
    
    ttbar->SetLineColor(kBlack);
    ttbar->SetFillColor(kOrange-7);
    DATA_MC->Add(ttbar, "HIST");
    //ttbar->Draw("sameHIST");
    c3->Update();
    
    
    
    Wjets->SetLineColor(kBlack);
    Wjets->SetFillColor(kRed-4);
    DATA_MC->Add(Wjets, "HIST");
    //Wjets->Draw("sameHIST");
    c3->Update();
    
    
    
    qcd->SetLineColor(kBlack);
    qcd->SetFillColor(kOrange+1);
    DATA_MC->Add(qcd, "HIST");
    //qcd->Draw("sameHIST");
    c3->Update();
    
    dy->SetLineColor(kBlack);
    dy->SetFillColor(kYellow-8);
    DATA_MC->Add(dy, "HIST");
    //qcd->Draw("sameHIST");
    c3->Update();
    
    //Removing the contribution of DY from data
    //TH1D* Z = (TH1D*) dy->Clone();
    //TH1D* d = (TH1D*) DATA->Clone();
    //d->Add(dy,-1);
    
    /*  d->SetMarkerStyle(20);
     d->SetMarkerColor(kBlack);
     d->SetMarkerSize(0.95);
     d->SetMarkerSize(0.5);
     DATA_MC->Add(d, "P");
     c3->Update();
     */
    
    DATA->SetMarkerStyle(20);
    DATA->SetMarkerColor(kBlack);
    DATA->SetMarkerSize(0.95);
    // DATA_MC->Add(d, "P");
    
    
    
    DATA_MC->Draw();
    c3->Update();
    
    DATA->Draw("SAMEPE");
    c3->Update();
    
    
    DATA_MC->GetXaxis()->SetTitle(" ");
    DATA_MC->GetYaxis()->SetTitle("#Events");
    DATA_MC->GetYaxis()->SetTitleSize(18);
    DATA_MC->GetYaxis()->SetTitleFont(43);
    DATA_MC->GetYaxis()->SetTitleOffset(1.55);
    DATA_MC->SetMinimum(0.0001);
    DATA_MC->SetMaximum(10000000.0);
    DATA_MC->GetXaxis()->SetLabelSize(0);
    //DATA_MC->GetYaxis()->SetRangeUser(0.000001, 1000000);
    DATA_MC->SetTitle(name);
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
    
    
    
    
    TLegend *l1 = new TLegend(0.75,0.65,0.95,0.95);
    //l1->SetBorderSize(0);
    l1->AddEntry(DATA, "DATA", "lep");
    l1->AddEntry(dy, "#gamma^{*}/Z #rightarrow #mu^{+}#mu^{-}", "f");
    l1->AddEntry(qcd, "qcd", "f");
    l1->AddEntry(Wjets, "W + jets", "f");
    l1->AddEntry(ttbar, "t#bar{t}", "f");
    l1->AddEntry(WW, "WW", "f");
    l1->AddEntry(WZ, "WZ", "f");
    l1->AddEntry(ZZ, "ZZ", "f");
    l1->AddEntry(tW, "tW", "f");
    l1->AddEntry(Wantitop, "#bar{t}W", "f");
    l1->AddEntry(Zjets, "#gamma^{*}/Z #rightarrow #tau^{+}#tau^{-}", "f");
    
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
    
    
    
    
    TH1D* ratio = (TH1D*) DATA->Clone();
    ratio->Sumw2();
    TH1D* MC = (TH1D*) qcd->Clone();
    MC->Sumw2();
    MC->Add(dy);
    MC->Add(Wjets);
    MC->Add(Zjets);
    MC->Add(ttbar);
    MC->Add(WW);
    MC->Add(WZ);
    MC->Add(ZZ);
    MC->Add(Wantitop);
    MC->Add(tW);
    //MC->Sumw2();
    //   ratio->Add(MC,-1);
    ratio->Divide(MC);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerColor(kBlack);
    ratio->SetMarkerSize(0.95);
    
    
    std::cout<<"data/MC Ratio Info"<<std::endl;
    std::cout<<"Bin"<<" "<<"Bin content"<<"   "<<"Bin error"<<std::endl;
    for(int i=1; i<9; i++){
        double r = ratio->GetBinContent(i);
        double err = ratio->GetBinError(i);
        //double err_cal = r * sqrt((1/DATA->GetBinContent(i)) + (1/MC->GetBinContent(i)));
        //ratio->SetBinError(i,err);
        std::cout<<i<<" "<<r<<"   "<<err<<std::endl;
        
    }
    
    
    
    
    ratio->GetYaxis()->SetTitle("Data / Bkg");
    ratio->GetYaxis()->SetTitleSize(18);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(15);
    ratio->GetYaxis()->SetTitleOffset(1.55); //Distance between the axis and the title
    
    
    ratio->GetXaxis()->SetTitle(title);
    ratio->GetXaxis()->SetTitleSize(18);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetLabelFont(43);
    if(logx){
        ratio->GetXaxis()->SetLabelSize(0);
        ratio->GetXaxis()->SetLabelOffset(999);
    }
    else{
        ratio->GetXaxis()->SetLabelFont(43);
        ratio->GetXaxis()->SetLabelSize(15);
    }
    ratio->GetXaxis()->SetTitleOffset(4.);
    ratio->GetYaxis()->SetRangeUser(0, 2);
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
    
    
    //For x axis lables for log scale
    if(logx){
        const Int_t n = 9;
        ratio->GetXaxis()->SetLabelSize(15);
        char *pt[n] = {"50","100","150","200","300","400", "500", "1000", "5000"};
        Float_t x, y;
        y = Ymin - 0.2*ratio->GetYaxis()->GetBinWidth(1);
        TText t;
        //t.SetTextAngle(60);
        t.SetTextSize(0.07);
        t.SetTextAlign(22);
        for (int i=0;i<9;i++) {
            x = ratio->GetXaxis()->GetBinLowEdge(i+1);
            t.DrawText(x,y,pt[i]);
        }
    }
    
    
    
    
    
    TLine* line = new TLine(Xmin, 1, Xmax, 1);
    line->SetLineColor(kGreen);
    line->SetLineWidth(1);
    if(logx) line->DrawLine(pow(10,Xmin), 1, pow(10,Xmax), 1);
    else
        line->DrawLine(Xmin, 1, Xmax, 1);
    pad2->Update();
    
    
    c3->Update();
    
    c3->SaveAs(name+".pdf","pdf");
    c3->SaveAs(name+".png","png");
    
    
    
    
    
    
}


void RatePlot(TH1D* Wjets_den, TH1D* Zjets_den, TH1D* ttbar_den, TH1D* WW_den, TH1D* WZ_den, TH1D* ZZ_den, TH1D* tW_den, TH1D* Wantitop_den, TH1D* dy_den, TH1D* Data_den, TH1D* Wjets_num, TH1D* Zjets_num, TH1D* ttbar_num, TH1D* WW_num, TH1D* WZ_num, TH1D* ZZ_num, TH1D* tW_num, TH1D* Wantitop_num, TH1D* dy_num, TH1D* qcd_den, TH1D* qcd_num, TH1D* Data_num, TString title, TString name){
    
    TCanvas *c3 = new TCanvas("rate", "rate", 800, 800);
    c3->cd();
    
    
    
    TH1D* total_den = (TH1D*) Data_den->Clone();
    total_den->Sumw2();
    TH1D* MC_den = (TH1D*) Wjets_den->Clone();
    MC_den->Sumw2();
    MC_den->Add(dy_den);
    MC_den->Add(Zjets_den);
    MC_den->Add(ttbar_den);
    MC_den->Add(WW_den);
    MC_den->Add(WZ_den);
    MC_den->Add(ZZ_den);
    MC_den->Add(Wantitop_den);
    MC_den->Add(tW_den);
    total_den->Add(MC_den,-1);
    
    
    
    TH1D* total_num = (TH1D*) Data_num->Clone();
    total_num->Sumw2();
    TH1D* MC_num = (TH1D*) Wjets_num->Clone();
    MC_num->Sumw2();
    MC_num->Add(dy_num);
    MC_num->Add(Zjets_num);
    MC_num->Add(ttbar_num);
    MC_num->Add(WW_num);
    MC_num->Add(WZ_num);
    MC_num->Add(ZZ_num);
    MC_num->Add(Wantitop_num);
    MC_num->Add(tW_num);
    total_num->Add(MC_num,-1);
    
    
    TH1D* FR_pt= (TH1D*) total_num->Clone();
    FR_pt->Sumw2();
    //TArrayD *p = total_num->GetSumw2();
    //Int_t  n = sizeof(p)/sizeof(p[0]);
    //for(int i=0; i<n ; n++){
    //if (p && (p->fN > 0)) p->Set(0); // "destroy" fSumw2,
    // FR_pt->Sumw2();
    //FR_pt->Divide(total_den);
    //}
    std::cout<<"FR info"<<std::endl;
    
    double num=0.0;
    double num_err = 0.0;
    double den=0.0;
    double den_err = 0.0;
    double r =0.0;
    
    for(int i=1; i<9; i++){
        
        num = abs(total_num->GetBinContent(i));
        den = abs(total_den->GetBinContent(i));
        num_err = total_num->GetBinError(i);
        den_err = total_den->GetBinError(i);
        r = num/den;
        //double err = r * sqrt(((Data_num->GetBinContent(i)+MC_num->GetBinContent(i))/(total_num->GetBinContent(i)*total_num->GetBinContent(i))) + ((Data_den->GetBinContent(i)+MC_den->GetBinContent(i))/(total_den->GetBinContent(i)*total_den->GetBinContent(i))));
        double err = r * sqrt((num_err*num_err)/(num*num) + (den_err*den_err)/(den*den));
        //std::cout<<i<<" "<<total_num->GetBinContent(i)<<"   "<<total_den->GetBinContent(i)<<"   "<<num<<"   "<<den<<"   "<<r<<" "<<err<<std::endl;
        FR_pt->SetBinContent(i,r);
        FR_pt->SetBinError(i,err);
    }
    
    
    
    
    
    
    // To print numerator and denominator events in each bin in each sample
    
    /*  std::cout<<"Data"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<Data_num->GetBinContent(i)<<"    "<<Data_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"Wjets"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<Wjets_num->GetBinContent(i)<<"    "<<Wjets_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"Zjets"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<Zjets_num->GetBinContent(i)<<"    "<<Zjets_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"ttbar"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<ttbar_num->GetBinContent(i)<<"    "<<ttbar_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"WW"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<WW_num->GetBinContent(i)<<"    "<<WW_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"WZ"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<WZ_num->GetBinContent(i)<<"    "<<WZ_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"ZZ"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<ZZ_num->GetBinContent(i)<<"    "<<ZZ_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"tW"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<tW_num->GetBinContent(i)<<"    "<<tW_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"Wantitop"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<Wantitop_num->GetBinContent(i)<<"    "<<Wantitop_den->GetBinContent(i)<<std::endl;
     }
     
     std::cout<<"dy"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<dy_num->GetBinContent(i)<<"    "<<dy_den->GetBinContent(i)<<std::endl;
     }
     
     
     std::cout<<"QCD"<<std::endl;
     std::cout<<"Bin"<<" "<<"Numerator"<<"   "<<"Denominator"<<" "<<std::endl;
     for(int i=1; i<9; i++){
     std::cout<<i<<" "<<qcd_num->GetBinContent(i)<<"    "<<qcd_den->GetBinContent(i)<<std::endl;
     }
     */
    //////////////////////////////////////////////
    
    //FR_pt->Sumw2();
    
    std::cout<<"Bin"<<" "<<"Bin Content"<<"   "<<"Bin error"<<std::endl;
    for(int i=1; i<9; i++) {
        std::cout<<i<<" "<<FR_pt->GetBinContent(i)<<"   "<<FR_pt->GetBinError(i)<<std::endl;
    }
    
    FR_pt->SetMarkerStyle(20);
    FR_pt->SetMarkerColor(kBlack);
    FR_pt->SetMarkerSize(0.95);
    
    
    FR_pt->GetYaxis()->SetTitle("Fake Rate");
    FR_pt->GetXaxis()->SetTitle("pT(#mu)[GeV]");
    FR_pt->SetTitle(title);
    FR_pt->GetYaxis()->SetRangeUser(0, 1.5);
    c3->SetGrid();
    
    FR_pt->Draw("PE");
    FR_pt->SaveAs(name+".root","root");
    c3->Update();
    
    TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText2->SetBorderSize(0);
    tText2->SetFillColor(0);
    tText2->SetFillStyle(0);
    TText *t2 = tText2->AddText("CMS 2018 Data 61.3 fb^{-1} (13TeV)");
    tText2->SetTextSize(0.035);
    tText2->Draw();
    c3->Update();
    
    
    
    c3->SaveAs(name+".pdf","pdf");
    c3->SaveAs(name+".png","png");
    
}



void QCDPlot(TH1D* qcd_den, TH1D* qcd_num, TString name){
    
    TCanvas *c3 = new TCanvas("rate", "rate", 800, 800);
    c3->cd();
    
    TH1D* qcd_d = (TH1D*) qcd_den->Clone();
    qcd_d->Sumw2();
    TH1D* qcd_n = (TH1D*) qcd_num->Clone();
    qcd_n->Sumw2();
    qcd_n->Divide(qcd_d);
    
    std::cout<< "FR qcd" <<std::endl;
    std::cout<<"Bin"<<" "<<"Bin Content"<<"   "<<"Bin error"<<std::endl;
    for(int i=1; i<9; i++){
        Double_t binContent = qcd_n->GetBinContent(i);
        
        //double err = binContent * sqrt(1/qcd_den->GetBinContent(i) + 1/qcd_num->GetBinContent(i));
        std::cout<<i<<"   "<<binContent<<"  "<<qcd_n->GetBinError(i)<<std::endl;
        //qcd_n->SetBinError(i,err);
    }
    
    
    qcd_n->SetMarkerStyle(20);
    qcd_n->SetMarkerColor(kBlue);
    qcd_n->SetMarkerSize(0.95);
    
    
    qcd_n->GetYaxis()->SetTitle("QCD Fake Rate");
    qcd_n->GetXaxis()->SetTitle("pT(#mu)[GeV]");
    qcd_n->SetTitle(name);
    qcd_n->GetYaxis()->SetRangeUser(0, 1.5);
    c3->SetGrid();
    
    qcd_n->Draw("PE");
    qcd_n->SaveAs(name+".root","root");
    c3->Update();
    
    TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText2->SetBorderSize(0);
    tText2->SetFillColor(0);
    tText2->SetFillStyle(0);
    TText *t2 = tText2->AddText("CMS 2018 Data 61.3 fb^{-1} (13TeV)");
    tText2->SetTextSize(0.035);
    tText2->Draw();
    c3->Update();
    
    
    c3->SaveAs(name+".pdf","pdf");
    c3->SaveAs(name+".png","png");
    
}

void EWKPlot(TH1D* Wjets_den, TH1D* Zjets_den, TH1D* ttbar_den, TH1D* WW_den, TH1D* WZ_den, TH1D* ZZ_den, TH1D* tW_den, TH1D* Wantitop_den, TH1D* dy_den, TH1D* Wjets_num, TH1D* Zjets_num, TH1D* ttbar_num, TH1D* WW_num, TH1D* WZ_num, TH1D* ZZ_num, TH1D* tW_num, TH1D* Wantitop_num, TH1D* dy_num, TString title, TString name){
    
    TCanvas *c3 = new TCanvas("rate", "rate", 800, 800);
    c3->cd();
    
    
    
    
    TH1D* MC_den = (TH1D*) Wjets_den->Clone();
    MC_den->Add(dy_den);
    MC_den->Add(Zjets_den);
    MC_den->Add(ttbar_den);
    MC_den->Add(WW_den);
    MC_den->Add(WZ_den);
    MC_den->Add(ZZ_den);
    MC_den->Add(Wantitop_den);
    MC_den->Add(tW_den);
    
    
    
    TH1D* MC_num = (TH1D*) Wjets_num->Clone();
    MC_num->Add(dy_num);
    MC_num->Add(Zjets_num);
    MC_num->Add(ttbar_num);
    MC_num->Add(WW_num);
    MC_num->Add(WZ_num);
    MC_num->Add(ZZ_num);
    MC_num->Add(Wantitop_num);
    MC_num->Add(tW_num);
    
    
    
    
    
    MC_den->SetMarkerStyle(20);
    MC_den->SetMarkerColor(kBlack);
    MC_den->SetMarkerSize(0.95);
    
    
    MC_den->GetYaxis()->SetTitle("EWK Fake Rate");
    MC_den->GetXaxis()->SetTitle("pT(#mu)[GeV]");
    MC_den->SetTitle(title);
    MC_den->GetYaxis()->SetRangeUser(-1, 1.5);
    c3->SetGrid();
    
    //MC_den->Rebin(5);
    MC_den->Draw("PE");
    c3->Update();
    
    MC_num->SetMarkerStyle(20);
    MC_num->SetMarkerColor(kRed);
    MC_num->SetMarkerSize(0.95);
    //MC_num->Rebin(5);
    MC_num->Draw("SAMEPE");
    c3->Update();
    
    TPaveText* tText2 = new TPaveText(0.2, 0.90, 0.4, 0.95, "brNDC");
    tText2->SetBorderSize(0);
    tText2->SetFillColor(0);
    tText2->SetFillStyle(0);
    TText *t2 = tText2->AddText("CMS 2017 Data 42.4 fb^{-1} (13TeV)");
    tText2->SetTextSize(0.035);
    tText2->Draw();
    c3->Update();
    
    TLegend *l1 = new TLegend(0.75,0.65,0.95,0.95);
    //l1->SetBorderSize(0);
    l1->AddEntry(MC_den, "EWK den", "lep");
    l1->AddEntry(MC_num, "EWK num", "lep");
    l1->Draw();
    c3->Update();
    c3->cd();
    
    // for(int bin=1; bin<61; bin++){
    // Double_t binContent = FR_pt->GetBinContent(bin);
    // std::cout<<bin<<"   "<<binContent<<std::endl;
    // }
    
    c3->SaveAs(name+".pdf","pdf");
    c3->SaveAs(name+".png","png");
    
}

/*
 void PrintKine(double array[600][20]){
 
 std::cout<<"Printing muon kinematics"<<std::endl;
 
 std::cout<<"No"<<";"<<"Mass[GeV]"<<";"<<"lep_id_muon(-)"<<";"<<"lep_id_muon(+)"<<";"<<"TuneP_pt[GeV]_muon(-)"<<";"<<"TuneP_pt[GeV]_muon(+)"<<";"<<"Tracker_pt[GeV]_muon(-)"<<";"<<"Tracker_pt[GeV]_muon(+)"<<";"<<"Ratio(TuneP/Tracker)_muon(-)"<<";"<<"Ratio(TuneP/Tracker)_muon(+)"<<";"<<"eta_muon(-)"<<";"<<"eta_muon(+)"<<";"<<"phi_muon(-)"<<";"<<"phi_muon(+)"<<";"<<"MET_pt[GeV]"<<";"<<"MET_phi"<<";"<<"Category"<<";"<<"Run"<<";"<<"Lumi"<<";"<<"Event"<<std::endl;
 
 
 std::cout<<std::setw(10)<<"No"<<std::setw(10)<<"Mass[GeV]"<<std::setw(10)<<"TuneP_pt[GeV]_muon(-)"<<std::setw(10)<<"TuneP_pt[GeV]_muon(+)"<<std::setw(10)<<"Tracker_pt[GeV]_muon(-)"<<std::setw(10)<<"Tracker_pt[GeV]_muon(+)"<<std::setw(10)<<"Ratio(TuneP/Tracker)_muon(-)"<<std::setw(10)<<"Ratio(TuneP/Tracker)_muon(+)"<<std::setw(10)<<"eta_muon(-)"<<std::setw(10)<<"eta_muon(+)"<<std::setw(10)<<"phi_muon(-)"<<std::setw(10)<<"phi_muon(+)"<<std::setw(10)<<"MET_pt[GeV]"<<std::setw(10)<<"MET_phi"<<std::setw(10)<<"Category"<<std::setw(10)<<"Run"<<std::setw(10)<<"Lumi"<<std::setw(10)<<"Event"<<std::endl;
 
 
 for(int r=0; r<400; r++){
 
 for(int c=0; c<20; c++){
 std::cout<<array[r][c]<<";";
 }
 std::cout<<std::endl;
 }
 
 }
 */







