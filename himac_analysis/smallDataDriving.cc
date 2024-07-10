#include <iostream>
#include <sstream>
#include <string>

#include "GETDecoder.hh"
#include "GETPad.hh"
#include "GETSiPad.hh"
#include "dataStructure.hh"

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TClonesArray.h"
#include "TGaxis.h"
#include "TCutG.h"

using namespace TMath;
using namespace std;

// ========= ATTPC-constant ===============
const double fTPCHieghtToBeam = 89.; //[mm]
const double fTimeoffsetL = 239.45; //[mm] //5.7-245.35 //5.6-239.45 //5.5-233.60 //5.4-227.75 //5.3-221.86 //5.8-251.21
const double fTimeoffsetR = 238.08; //[mm] //5.7-243.91 //5.6-238.08 //5.5-232.25 //5.4-226.41 //5.3-220.58 //5.8-249.76
const double fTPCADCThreshold = 50.; // [ADC]
const double fPadHeight = 11.9; //[mm]
const double fPadWeith = 1.9; //[mm]
const double fPadGap = 0.1; //[mm]
const double fhalfH = 89; //[mm]
const double fDistance_TPC = 550; //[mm]
const double fDistance_Si = 715.1; //[mm]
// ========= ATTPC-constant ===============

void smallDataDriving();
void checkCSI();
void checkADC();
void numEvent();

double getTargetZ(int wing, int tracks, double s, double c){
    if(wing == 0){
        // return (-fDistance_Si * s + ( x - 37.5 )) / (Sin(40*Pi()/180) + s *Cos(40*Pi()/180)); 
        // return (-(fDistance_TPC-(fPadHeight+fPadGap)*3.5)*s + c - fhalfH - fTimeoffsetR)/(Sin(40*Pi()/180)+ s * Cos(40*Pi()/180));
        return (-(fDistance_TPC-(fPadHeight+fPadGap)*3.5)*s + c - fhalfH - fTimeoffsetR)/(Cos(50*Pi()/180));
    }else{
        // return (-fDistance_Si * s - ( x - 37.5 )) /(-Sin(40*Pi()/180) + s *Cos(40*Pi()/180)); 
        // return -((fDistance_TPC+(fPadHeight+fPadGap)*3.5)*s + c - fhalfH - fTimeoffsetL)/(-Sin(40*Pi()/180)+ s * Cos(40*Pi()/180));
        return ((fDistance_TPC+(fPadHeight+fPadGap)*3.5)*s + c - fhalfH - fTimeoffsetL)/(Cos(50*Pi()/180));
    }
    
}

int main()
{
    smallDataDriving();
    // checkCSI();
    // numEvent();

    return 0;
}


void smallDataDriving()
{
    TString TypeName[4] = {"Right-Down", "Right-Up", "Left-Down", "Left-Up"};
    const double tpcADCThreshold = 50.;

    TFile* fileIn = new TFile("../data/SmallHIMACData_v5.6.root","read");
    TTree* treeIn = (TTree*)fileIn -> Get("data");

    double trackFitData[2][2][2][3]; //[Right, Left][trkNum][xy, yz][slope, const]
    double trackData[2][2][8][32][4]; //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]

    double csiData[2][2][2]; //[Right, Left][Down, Up][Energy sumADC]
    int siHitNum;
    double siData[2][5][3][2]; //[Right, Left][hitNum][Energy, ADC, pos] Energy=(omic, junc), ADC=(omic, junc), pos=(x, y)

    double siDataBeamADCOmic[2][4][512]; // [Right, Left][omic][tb]
    double siDataBeamADCJunc[2][16][512]; // [Right, Left][junc][tb]
    double siDataBeamHitOmic[2][4];  // [Right, Left][omic][hit]
    double siDataBeamHitJunc[2][16]; // [Right, Left][junc][hit]
    double degree = 180. / Pi();
    // cout<<degree<<endl;

    treeIn -> SetBranchAddress("trackFitData", &trackFitData);
    treeIn -> SetBranchAddress("trkData", &trackData);
    treeIn -> SetBranchAddress("csiData", &csiData);
    treeIn -> SetBranchAddress("siHitNum", &siHitNum);
    treeIn -> SetBranchAddress("siData", &siData);
    treeIn -> SetBranchAddress("siBADCOmic", &siDataBeamADCOmic);
    treeIn -> SetBranchAddress("siBADCJunc", &siDataBeamADCJunc);
    treeIn -> SetBranchAddress("siBHitOmic", &siDataBeamHitOmic);
    treeIn -> SetBranchAddress("siBHitJunc", &siDataBeamHitJunc);

    auto hThetaYX = new TH1D("hThetaYX", " ;#theta_{XY} [#circ]; Counts", 20, -5, 5); hThetaYX->SetStats(0);
    auto hYZR_dz = new TH1D("hYZR_dz-DT", " ;z_{vtx} (mm); Counts", 100, -200, 100); hYZR_dz->SetLineColor(kBlue); hYZR_dz->GetYaxis()->SetRangeUser(1., 4000.);  
    auto hYZL_dz = new TH1D("hYZL_dz-DT", " ;z_{vtx} (mm); Counts", 100, -200, 100); hYZL_dz->SetLineColor(kRed); hYZL_dz->GetYaxis()->SetRangeUser(1., 4000.); 
    auto hYZR_dzp = new TH1D("hYZR_dz-DTp", " ;#Delta Z [mm]", 100, -200, 100); hYZR_dzp->SetLineColor(kBlue);
    auto hYZL_dzp = new TH1D("hYZL_dz-DTp", " ;#Delta Z [mm]", 100, -200, 100); hYZL_dzp->SetLineColor(kRed);
    auto h2dz = new TH2D("h2dz", " ;TPC-Right z_{vtx} (mm);TPC-Left z_{vtx} (mm)", 25, -200, 100, 25, -199.9, 100); 
    TLegend* sc = new TLegend(0.2, 0.57, 0.5, 0.87); sc->SetFillStyle(0); sc->SetBorderSize(0); sc -> SetTextSize(0.045);



    TH1D* hTPCdedx[2]; 
    TH1D* dedxdeut[2][2];
    TH1D* dedxAlpha[2][2];
    TH1D* dedxP[2][2];
    TH1D* dedxEle[2][2];
    for(int i=0; i<2; i++){
        hTPCdedx[i] = new TH1D(Form("hTPCdedx_%i", i),"", 100, 0., 1000);
        dedxAlpha[i][0] = new TH1D(Form("dedxAlpha_d_%i", i),"", 100, 0., 1000);
        dedxAlpha[i][1] = new TH1D(Form("dedxAlpha_u_%i", i),"", 100, 0., 1000);
        dedxdeut[i][0] = new TH1D(Form("dedxdeut_d_%i", i),"", 100, 0., 1000);
        dedxdeut[i][1] = new TH1D(Form("dedxdeut_u_%i", i),"", 100, 0., 1000);
        dedxP[i][0] = new TH1D(Form("dedxP_d_%i", i),"", 100, 0., 1000);
        dedxP[i][1] = new TH1D(Form("dedxP_u_%i", i),"", 100, 0., 1000);
        dedxEle[i][0] = new TH1D(Form("dedxElec_d_%i", i),"", 100, 0., 1000);
        dedxEle[i][1] = new TH1D(Form("dedxElec_u_%i", i),"", 100, 0., 1000);
    }

    TH2D* pidTotal = new TH2D("hPID", "", 100, 0., 210, 100, 0, 14);
    TH2D* hPID[2][2]; // wing, down up
    TH2D* htpcPID[2][2]; // wing down up
    for(int wing=0; wing<2; wing++){
        for(int du=0; du<2; du++){
            hPID[wing][du] = new TH2D(Form("hPID_%i_%i", wing, du), "", 100, 0., 210, 100, 0, 14);
            htpcPID[wing][du] = new TH2D(Form("htpcPID_%i_%i", wing, du), "", 100, 0., 210, 100, 0, 10);
        }
    }
    auto p_cut = new TCutG("p_cut", 11);
    p_cut-> SetPoint(0, 0, 5);
    p_cut-> SetPoint(1, 40, 3.5);
    p_cut-> SetPoint(2, 80, 3);
    p_cut-> SetPoint(3, 120, 2);
    p_cut-> SetPoint(4, 140, 2);
    p_cut-> SetPoint(5, 140, 3.5);
    p_cut-> SetPoint(6, 120, 3.5);
    p_cut-> SetPoint(7, 80, 4);
    p_cut-> SetPoint(8, 40, 5);
    p_cut-> SetPoint(9, 0, 7);
    p_cut-> SetPoint(10, 0, 5);

    // auto p_cut_LUsi = new TCutG("p_cut_LUsi", 11);
    // p_cut_LUsi-> SetPoint(0, 0, 3.5);
    // p_cut_LUsi-> SetPoint(1, 40, 2.5);
    // p_cut_LUsi-> SetPoint(2, 80, 2);
    // p_cut_LUsi-> SetPoint(3, 120, 1.5);
    // p_cut_LUsi-> SetPoint(4, 140, 1.5);
    // p_cut_LUsi-> SetPoint(5, 140, 2.5);
    // p_cut_LUsi-> SetPoint(6, 120, 2.5);
    // p_cut_LUsi-> SetPoint(7, 80, 3);
    // p_cut_LUsi-> SetPoint(8, 40, 4);
    // p_cut_LUsi-> SetPoint(9, 0, 5);
    // p_cut_LUsi-> SetPoint(10, 0, 3.5);

    // auto d_cut = new TCutG("d_cut", 11);
    // d_cut-> SetPoint(0, 0, 7.5);
    // d_cut-> SetPoint(1, 40, 6);
    // d_cut-> SetPoint(2, 80, 4.5);
    // d_cut-> SetPoint(3, 120, 3.5);
    // d_cut-> SetPoint(4, 140, 3.5);
    // d_cut-> SetPoint(5, 140, 5);
    // d_cut-> SetPoint(6, 120, 5);
    // d_cut-> SetPoint(7, 80, 6);
    // d_cut-> SetPoint(8, 40, 7.5);
    // d_cut-> SetPoint(9, 0, 9.5);
    // d_cut-> SetPoint(10, 0, 7.5);

    // auto t_cut = new TCutG("t_cut", 11);
    // t_cut-> SetPoint(0, 0, 10);
    // t_cut-> SetPoint(1, 40, 6.5);
    // t_cut-> SetPoint(2, 80, 5.5);
    // t_cut-> SetPoint(3, 120, 3.5);
    // t_cut-> SetPoint(4, 140, 3.5);
    // t_cut-> SetPoint(5, 140, 5.5);
    // t_cut-> SetPoint(6, 120, 5.5);
    // t_cut-> SetPoint(7, 80, 6.5);
    // t_cut-> SetPoint(8, 40, 8);
    // t_cut-> SetPoint(9, 0, 11.5);
    // t_cut-> SetPoint(10, 0, 10);
    double fit2, fit4;
    int eventNum = treeIn -> GetEntries();
    for(int event=0; event<eventNum; event++){
        if(event%1000==0){cout << "event: " << event << " / " << eventNum << endl;}
        treeIn -> GetEntry(event);

        double dEdxWing[2];
        double CsIEWing[2];
        double SiWing[2];
        memset(dEdxWing, 0, sizeof(dEdxWing));
        memset(CsIEWing, 0, sizeof(CsIEWing));
        memset(SiWing, 0, sizeof(SiWing));

        double dzR[10] = { 0 };
        double dzL[10] = { 0 };
        int numR = 0;
        int numL = 0;
        int proton[2] = { 0 };
        for(int wing=0; wing<2; wing++){ // right, left
        double CsIDE = csiData[wing][0][0];
            double CsIUE = csiData[wing][1][0];
            double CsIDADC = csiData[wing][0][1];
            double CsIUADC = csiData[wing][1][1];

            if(siHitNum < 1){continue;}

            bool isUp = false;
            double csiEnergy[2];
            double siEnergy[2];
            memset(csiEnergy, 0, sizeof(csiEnergy));
            memset(siEnergy, 0, sizeof(siEnergy));
            if(siHitNum > 1){continue;}

            for(int i=0; i<siHitNum; i++){
                double adc = siData[wing][i][1][1];
                if(adc < 10){continue;}
                double siPosY = siData[wing][i][2][1];
                double omicADC = siData[wing][i][1][0];
                double juncADC = siData[wing][i][1][1];
                double omicE = siData[wing][i][0][0];
                double juncE = siData[wing][i][0][1];

                double siADC = omicADC + juncADC;
                double siE = omicE + juncE;

                if(siPosY < 0.){
                    hPID[wing][0] -> Fill(CsIDE, siE);
                    csiEnergy[i] = CsIDE;
                    siEnergy[i] = siE;
                    pidTotal -> Fill(CsIDE, siE);
                    if( p_cut->IsInside(CsIDE, siE) ){
                        proton[wing] = 1;
                    }
                }
                else{
                    hPID[wing][1] -> Fill(CsIUE, siE);
                    csiEnergy[i] = CsIUE;
                    siEnergy[i] = siE;
                    isUp = true;
                    pidTotal -> Fill(CsIUE, siE);
                    if (wing == 1){
                        if(p_cut -> IsInside(CsIUE, siE)){
                            proton[wing] = 1;
                        }
                    }else{
                        if(p_cut -> IsInside(CsIUE, siE)){
                            proton[wing] = 1;
                        }
                    }
                }
            }
            for(int trkIdx=0; trkIdx<1; trkIdx++){     
                
                if( trackFitData[wing][trkIdx][1][0] < -990 || trackFitData[wing][trkIdx][1][1] < -990 ) { continue; }
                hThetaYX->Fill(degree * (ATan(1. / trackFitData[wing][trkIdx][0][0])));
                // cout<< (180. / Pi()) * ATan(1. / trackFitData[wing][trkIdx][0][0])<<endl;
                // double xySlope = 1./tan((angle - 180.)*(Pi()/180.));
                double dz_drift = getTargetZ(wing, trkIdx, trackFitData[wing][trkIdx][1][0], trackFitData[wing][trkIdx][1][1] );
                if(wing == 0){
                    dzR[numR] = dz_drift;
                    hYZR_dz->Fill(dz_drift);
                    numR++;
                    if(proton[0] == 1){
                        if(dz_drift > -130 && dz_drift < -116){
                            fit2++;
                        }
                        if(dz_drift > -12 && dz_drift < 12){
                            fit4++;
                        }
                    }
                    
                }else{
                    dzL[numL] = dz_drift;
                    hYZL_dz->Fill(dz_drift);
                    numL++;
                }
            }    
        }
        
        if(numR > 0 && numL > 0){
            for(int i=0; i < numR; i++){
                for(int j=0; j < numL; j++){
                    h2dz->Fill(dzR[i], dzL[j]);
                }
            }
        }
        if(proton[0] == 1 && numR > 0){
            for(int i =0; i< numR; i++){
                hYZR_dzp->Fill(dzR[i]);
            }
        }
        if(proton[1] == 1 && numL > 0){
            for(int i =0; i< numL; i++){
                hYZL_dzp->Fill(dzL[i]);
            }
        }
    }
    double pars2[3], pars4[3];
    
    auto sub1 = new TF1("sub1", "gaus", -200, -160);
    auto sub2 = new TF1("sub2", "gaus", -145, -110);
    auto sub3 = new TF1("sub3", "gaus", -100, -60);
    auto sub4 = new TF1("sub4", "gaus", -30, 30);

    auto sub1R = new TF1("sub1R", "gaus", -200, -160);
    auto sub2R = new TF1("sub2R", "gaus", -145, -110);
    auto sub3R = new TF1("sub3R", "gaus", -100, -60);
    auto sub4R = new TF1("sub4R", "gaus", -30, 30);
    sub1R->SetLineColor(kBlue);
    sub2R->SetLineColor(kBlue);
    sub3R->SetLineColor(kBlue);
    sub4R->SetLineColor(kBlue);
    hYZR_dzp->Fit(sub2R, "R");
    hYZR_dzp->Fit(sub4R, "R");
    sub2R->GetParameters(pars2);
    sub4R->GetParameters(pars4);
    TCanvas* cTarget = new TCanvas("cTarget", "", 1500, 1500);
    cTarget->SetTopMargin(0.03);
    cTarget->SetRightMargin(0.05);
    cTarget->SetBottomMargin(0.11);
    cTarget->SetLeftMargin(0.16);

    TCanvas* cTarget1 = new TCanvas("cTarget1", "", 1600, 1500);
    cTarget1->SetTopMargin(0.03);
    cTarget1->SetRightMargin(0.12);
    cTarget1->SetBottomMargin(0.135);
    cTarget1->SetLeftMargin(0.14);
    auto latex = new TLatex();
    latex->SetTextSize(0.085);

    auto latex1 = new TLatex();
    latex1->SetTextSize(0.03);
    // cTarget->Divide(2,1);

    hYZL_dz->GetYaxis()->SetTitleOffset(1.65);
    hYZL_dz->GetXaxis()->SetTitleOffset(0.94);
    hYZL_dz->GetYaxis()->SetLabelSize(0.051);
    hYZL_dz->GetXaxis()->SetLabelSize(0.051);
    hYZL_dz->GetYaxis()->SetTitleSize(0.053);
    hYZL_dz->GetXaxis()->SetTitleSize(0.053);
    hYZR_dz->SetLineWidth(2);
    hYZL_dz->SetLineWidth(2);
    hYZR_dz->SetStats(0);
    hYZL_dz->SetStats(0);
    hYZR_dzp->SetLineWidth(2);
    hYZL_dzp->SetLineWidth(2);
    hYZR_dzp->SetStats(0);
    hYZL_dzp->SetStats(0);
    hYZR_dzp->SetLineStyle(3);
    hYZL_dzp->SetLineStyle(3);

    h2dz->SetStats(0);
    h2dz->GetYaxis()->SetTitleOffset(1.25);
    h2dz->GetXaxis()->SetTitleOffset(1.2);
    h2dz->GetYaxis()->SetLabelSize(0.045);
    h2dz->GetXaxis()->SetLabelSize(0.045);
    h2dz->GetYaxis()->SetTitleSize(0.053);
    h2dz->GetXaxis()->SetTitleSize(0.053);
    sc->AddEntry(hYZR_dz,"TPC-Right","l");
    sc->AddEntry(hYZL_dz,"TPC-Left","l");
    sc->AddEntry(hYZR_dzp,"TPC-Right (proton)","l");
    sc->AddEntry(hYZL_dzp,"TPC-Left (proton)","l");
    sc->SetEntrySeparation(0.4);

    cTarget->cd();
    hYZL_dz->Draw();
    hYZR_dz->Draw("same");
    hYZR_dzp->Draw("same");
    hYZL_dzp->Draw("same");
    sub2R->Draw("same");
    sub4R->Draw("same");
    sc->Draw("same");
    latex->DrawLatexNDC(0.2, 0.9, "(a)");
    latex1->DrawLatexNDC(0.72, 0.9, "C_{2}H_{4} target");
    latex1->DrawLatexNDC(0.3, 0.4, Form("entries %.2f", fit2));
    latex1->DrawLatexNDC(0.8, 0.4, Form("entries %.2f", fit4));
    cTarget->SaveAs("../csiPicture/target_v=5.6ff.pdf");
    cTarget->SaveAs("../csiPicture/target_v=5.6ff.png");
    // latex -> SetTextSize(0.03);
    // latex -> DrawLatexNDC(0.15, 0.3, Form("#color[2]{%.1f}", par_sub1[1]));
    // latex -> DrawLatexNDC(0.27, 0.3, Form("#color[2]{%.1f}", par_sub2[1]));
    // latex -> DrawLatexNDC(0.4, 0.3, Form("#color[2]{%.1f}", par_sub3[1]));
    // latex -> DrawLatexNDC(0.75, 0.8, Form("#color[2]{mean : %.1f}", par_sub4[1]));
    // latex -> DrawLatexNDC(0.75, 0.75, Form("#color[2]{sigma : %.1f}", par_sub4[2]));

    // latex -> DrawLatexNDC(0.15, 0.25, Form("#color[4]{%.1f}", par_sub1R[1]));
    // latex -> DrawLatexNDC(0.27, 0.25, Form("#color[4]{%.1f}", par_sub2R[1]));
    // latex -> DrawLatexNDC(0.4, 0.25, Form("#color[4]{%.1f}", par_sub3R[1]));
    // latex -> DrawLatexNDC(0.75, 0.65, Form("#color[4]{mean : %.1f}", par_sub4R[1]));
    // latex -> DrawLatexNDC(0.75, 0.6, Form("#color[4]{sigma : %.1f}", par_sub4R[2]));
    // latex -> DrawLatexNDC(0.2, 0.75, Form("#color[2]{timeoffsetL %.2f}", fTimeoffsetL));
    // latex -> DrawLatexNDC(0.2, 0.7, Form("#color[4]{timeoffsetR %.2f}", fTimeoffsetR));

    cTarget1->cd();
    // h2dz->Draw("colz");
    // latex->DrawLatexNDC(0.2, 0.9, "(b)");
    // gPad->SetLogz();
    // // latex -> SetTextSize(0.06);
    // // latex -> DrawLatexNDC(0.2, 0.8, "V=5.8cm/#mus");
    // // cTarget->SaveAs("../csiPicture/target_v=5.8.png");
    // cTarget1->SaveAs("../csiPicture/2Dtarget_v=5.6.pdf");
    // cTarget1->SaveAs("../csiPicture/2Dtarget_v=5.6.png");
    hThetaYX->Draw();
    cTarget1->SaveAs("../csiPicture/angledis.png");


    // TCanvas* cPID = new TCanvas("cPID","", 2000,2000);
    // cPID -> Divide(2, 2);

    // int tmpIdx = 0;
    // for(int wing=0; wing<2; wing++){
    //     for(int du=0; du<2; du++){
    //         const int cidx = tmpIdx+1;
    //         cPID -> cd(cidx);
    //         cPID -> SetRightMargin(0.03);
    //         cPID -> SetTopMargin(0.03);
    //         gPad -> SetLogz();

    //         hPID[wing][du] -> SetStats(0);
    //         hPID[wing][du] -> SetTitle(Form("%s CsI-Si; CsI Energy [MeV]; Si Energy [MeV]", TypeName[tmpIdx].Data()));
    //         hPID[wing][du] -> Draw("colz");
    //         if(wing == 1 && du ==1 ){
    //             p_cut_LUsi->SetLineColor(kBlue);
    //             p_cut_LUsi->SetLineWidth(2);
    //             p_cut_LUsi->Draw("same");    
    //         }else{
    //             p_cut->SetLineColor(kBlue);
    //             p_cut->SetLineWidth(2);
    //             p_cut->Draw("same");

    //             d_cut->SetLineColor(kGreen);
    //             d_cut->SetLineWidth(2);
    //             d_cut->Draw("same");
    //         }
    //         tmpIdx++;
    //     }
    // }
    // cPID -> Draw();
    // cPID -> SaveAs("../csiPicture/PID_Energy_zoom.png");

    // TCanvas* cPIDtotal = new TCanvas();
    // pidTotal -> SetStats(0);
    // pidTotal -> SetTitle(Form("CsI-Si All; CsI Energy [MeV]; Si Energy [MeV]"));
    // pidTotal -> Draw("colz");
    // cPIDtotal -> Draw();
    // cPIDtotal -> SetLogz();
    // cPIDtotal -> SaveAs("../csiPicture/PID_CsiSiTotal.png");

    // TCanvas* cTPC = new TCanvas("","",1200,600);
    // cTPC -> Divide(2,1);
    // TLegend* sc1 = new TLegend(0.75, 0.55, 0.9, 0.9); sc1->SetFillStyle(0); sc1->SetBorderSize(0);
    // sc1->AddEntry(dedxAlpha[0][0],"#alpha","l");
    // sc1->AddEntry(dedxdeut[0][0],"D","l");
    // sc1->AddEntry(dedxP[0][0],"P","l");
    // sc1->AddEntry(dedxEle[0][0],"e","l");
    // sc1->SetEntrySeparation(0.4);

    // cTPC -> cd(1);
    // hTPCdedx[0] -> SetStats(0);
    // hTPCdedx[0] ->SetTitle("TPC-Right; dADCdx [ADC/mm]");
    // hTPCdedx[0] -> SetLineColor(kBlack);
    // hTPCdedx[0] -> Draw();
    // dedxP[0][0] -> SetLineColor(kBlue);
    // dedxP[0][0] -> Draw("same");
    // dedxdeut[0][0] -> SetLineColor(kGreen+1);
    // dedxdeut[0][0] -> SetLineStyle(2);
    // dedxdeut[0][0] -> Draw("same");
    // dedxAlpha[0][0] -> SetLineColor(kOrange);
    // dedxAlpha[0][0] -> SetLineStyle(2);
    // dedxAlpha[0][0] -> Draw("same");
    // dedxEle[0][0] -> SetLineColor(kRed);
    // dedxEle[0][0] -> SetLineStyle(2);
    // dedxEle[0][0] -> Draw("same");
    // sc1->Draw("same");


    // cTPC -> cd(2);
    // hTPCdedx[1] -> SetStats(0);
    // hTPCdedx[1] ->SetTitle("TPC-Left; dADCdx [ADC/mm]");
    // hTPCdedx[1] -> SetLineColor(kBlack);
    // hTPCdedx[1] -> Draw();
    // dedxP[1][0] -> SetLineColor(kBlue);
    // dedxP[1][0] -> Draw("same");
    // dedxdeut[1][0] -> SetLineColor(kGreen+1);
    // dedxdeut[1][0] -> SetLineStyle(2);
    // dedxdeut[1][0] -> Draw("same");
    // dedxAlpha[1][0] -> SetLineColor(kOrange);
    // dedxAlpha[1][0] -> SetLineStyle(2);
    // dedxAlpha[1][0] -> Draw("same");
    // dedxEle[1][0] -> SetLineColor(kRed);
    // dedxEle[1][0] -> SetLineStyle(2);
    // dedxEle[1][0] -> Draw("same");
    // sc1->Draw("same");
    

    // cTPC -> Draw();
    // cTPC -> SaveAs("../csiPicture/tpc_dedx.png");
}

void numEvent()
{
    TString TypeName[4] = {"Right-Down", "Right-Up", "Left-Down", "Left-Up"};
    const double tpcADCThreshold = 50.;

    TFile* fileIn = new TFile("../data/SmallHIMACData_v5.6.root","read");
    TTree* treeIn = (TTree*)fileIn -> Get("data");

    double trackFitData[2][2][2][3]; //[Right, Left][trkNum][xy, yz][slope, const]
    double trackData[2][2][8][32][4]; //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]

    double csiData[2][2][2]; //[Right, Left][Down, Up][Energy sumADC]
    int siHitNum;
    double siData[2][5][3][2]; //[Right, Left][hitNum][Energy, ADC, pos] Energy=(omic, junc), ADC=(omic, junc), pos=(x, y)

    double siDataBeamADCOmic[2][4][512]; // [Right, Left][omic][tb]
    double siDataBeamADCJunc[2][16][512]; // [Right, Left][junc][tb]
    double siDataBeamHitOmic[2][4];  // [Right, Left][omic][hit]
    double siDataBeamHitJunc[2][16]; // [Right, Left][junc][hit]

    treeIn -> SetBranchAddress("trackFitData", &trackFitData);
    treeIn -> SetBranchAddress("trkData", &trackData);
    treeIn -> SetBranchAddress("csiData", &csiData);
    treeIn -> SetBranchAddress("siHitNum", &siHitNum);
    treeIn -> SetBranchAddress("siData", &siData);
    treeIn -> SetBranchAddress("siBADCOmic", &siDataBeamADCOmic);
    treeIn -> SetBranchAddress("siBADCJunc", &siDataBeamADCJunc);
    treeIn -> SetBranchAddress("siBHitOmic", &siDataBeamHitOmic);
    treeIn -> SetBranchAddress("siBHitJunc", &siDataBeamHitJunc);

    auto hYZR_dz = new TH1D("hYZR_dz-DT", " ;#Delta Z [mm]; entries", 100, -200, 100); hYZR_dz->SetLineColor(kBlue);
    auto hYZL_dz = new TH1D("hYZL_dz-DT", " ;#Delta Z [mm]; entries", 100, -200, 100); hYZL_dz->SetLineColor(kRed);
    auto hYZR_dzp = new TH1D("hYZR_dz-DTp", " ;#Delta Z [mm]", 100, -200, 100); hYZR_dzp->SetLineColor(kBlue);
    auto hYZL_dzp = new TH1D("hYZL_dz-DTp", " ;#Delta Z [mm]", 100, -200, 100); hYZL_dzp->SetLineColor(kRed);
    auto h2dz = new TH2D("h2dz", " ;TPC-Right  #Delta Z [mm];TPC-Left  #Delta Z [mm]", 30, -200, 100, 30, -200, 100);
    TLegend* sc = new TLegend(0.15, 0.5, 0.5, 0.9); sc->SetFillStyle(0); sc->SetBorderSize(0);



    // TH1D* hTPCdedx[2]; 
    // TH1D* dedxdeut[2][2];
    // TH1D* dedxAlpha[2][2];
    // TH1D* dedxP[2][2];
    // TH1D* dedxEle[2][2];
    // for(int i=0; i<2; i++){
    //     hTPCdedx[i] = new TH1D(Form("hTPCdedx_%i", i),"", 100, 0., 1000);
    //     dedxAlpha[i][0] = new TH1D(Form("dedxAlpha_d_%i", i),"", 100, 0., 1000);
    //     dedxAlpha[i][1] = new TH1D(Form("dedxAlpha_u_%i", i),"", 100, 0., 1000);
    //     dedxdeut[i][0] = new TH1D(Form("dedxdeut_d_%i", i),"", 100, 0., 1000);
    //     dedxdeut[i][1] = new TH1D(Form("dedxdeut_u_%i", i),"", 100, 0., 1000);
    //     dedxP[i][0] = new TH1D(Form("dedxP_d_%i", i),"", 100, 0., 1000);
    //     dedxP[i][1] = new TH1D(Form("dedxP_u_%i", i),"", 100, 0., 1000);
    //     dedxEle[i][0] = new TH1D(Form("dedxElec_d_%i", i),"", 100, 0., 1000);
    //     dedxEle[i][1] = new TH1D(Form("dedxElec_u_%i", i),"", 100, 0., 1000);
    // }

    TH2D* pidTotal = new TH2D("hPID", "", 100, 0., 210, 100, 0, 14);
    TH2D* hPID[2][2]; // wing, down up
    TH2D* htpcPID[2][2]; // wing down up
    for(int wing=0; wing<2; wing++){
        for(int du=0; du<2; du++){
            hPID[wing][du] = new TH2D(Form("hPID_%i_%i", wing, du), "", 100, 0., 210, 100, 0, 14);
            htpcPID[wing][du] = new TH2D(Form("htpcPID_%i_%i", wing, du), "", 100, 0., 210, 100, 0, 10);
        }
    }
    TH2D* hSpecTotal = new TH2D("hSpecTotal", "", 100, 0., 210, 100, 0, 14);
    TH2D* hNumSpec[2][3];
    for(int wing =0; wing<2; wing++){
        for(int pid=0; pid<3; pid++){
            hNumSpec[wing][pid] = new TH2D(Form("hNumSpec_%i_%i", wing, pid), "", 100, 0, 210, 100, 0, 14);
        }
    }
    auto p_cut = new TCutG("p_cut", 13);
    p_cut-> SetPoint(0, 0, 5);
    p_cut-> SetPoint(1, 40, 3.5);
    p_cut-> SetPoint(2, 80, 3);
    p_cut-> SetPoint(3, 120, 2);
    p_cut-> SetPoint(4, 140, 2);
    p_cut-> SetPoint(5, 160, 2);
    p_cut-> SetPoint(6, 160, 2.5);
    p_cut-> SetPoint(7, 140, 3.5);
    p_cut-> SetPoint(8, 120, 3.7);
    p_cut-> SetPoint(9, 80, 4.3);
    p_cut-> SetPoint(10, 40, 5.5);
    p_cut-> SetPoint(11, 0, 7);
    p_cut-> SetPoint(12, 0, 5);

    

    auto d_cut = new TCutG("d_cut", 13);
    d_cut-> SetPoint(0, 0, 7.5);
    d_cut-> SetPoint(1, 40, 6);
    d_cut-> SetPoint(2, 80, 4.5);
    d_cut-> SetPoint(3, 120, 3.7);
    d_cut-> SetPoint(4, 140, 3.5);
    d_cut-> SetPoint(5, 160, 2.5);
    d_cut-> SetPoint(6, 160, 4);
    d_cut-> SetPoint(7, 140, 4.8);
    d_cut-> SetPoint(8, 120, 5);
    d_cut-> SetPoint(9, 80, 6.2);
    d_cut-> SetPoint(10, 40, 7.7);
    d_cut-> SetPoint(11, 0, 10);
    d_cut-> SetPoint(12, 0, 7.5);

    auto t_cut = new TCutG("t_cut", 13);
    t_cut-> SetPoint(0, 0, 10);
    t_cut-> SetPoint(1, 40, 7.7);
    t_cut-> SetPoint(2, 80, 6.2);
    t_cut-> SetPoint(3, 120, 5);
    t_cut-> SetPoint(4, 140, 4.8);
    t_cut-> SetPoint(5, 160, 4);
    t_cut-> SetPoint(6, 160, 6);
    t_cut-> SetPoint(7, 140, 6.5);
    t_cut-> SetPoint(8, 120, 7);
    t_cut-> SetPoint(9, 80, 8.2);
    t_cut-> SetPoint(10, 40, 9.7);
    t_cut-> SetPoint(11, 0, 11.5);
    t_cut-> SetPoint(12, 0, 10);

    auto p_cut_LUsi = new TCutG("p_cut_LUsi", 11);
    p_cut_LUsi-> SetPoint(0, 0, 3.5);
    p_cut_LUsi-> SetPoint(1, 40, 2.5);
    p_cut_LUsi-> SetPoint(2, 80, 2);
    p_cut_LUsi-> SetPoint(3, 120, 1.5);
    p_cut_LUsi-> SetPoint(4, 140, 1.5);
    p_cut_LUsi-> SetPoint(5, 160, 1.5);
    p_cut_LUsi-> SetPoint(6, 160, 2.5);
    p_cut_LUsi-> SetPoint(7, 140, 2.5);
    p_cut_LUsi-> SetPoint(8, 120, 2.5);
    p_cut_LUsi-> SetPoint(9, 80, 3.5);
    p_cut_LUsi-> SetPoint(10, 40, 4.5);
    p_cut_LUsi-> SetPoint(11, 0, 5.5);
    p_cut_LUsi-> SetPoint(12, 0, 3.5);
    auto d_cut_LUsi = new TCutG("d_cut_LUsi", 13);
    d_cut_LUsi-> SetPoint(0, 0, 5.5);
    d_cut_LUsi-> SetPoint(1, 40, 4.5);
    d_cut_LUsi-> SetPoint(2, 80, 3.5);
    d_cut_LUsi-> SetPoint(3, 120, 2.5);
    d_cut_LUsi-> SetPoint(4, 140, 2.5);
    d_cut_LUsi-> SetPoint(5, 160, 2.5);
    d_cut_LUsi-> SetPoint(6, 160, 3.5);
    d_cut_LUsi-> SetPoint(7, 140, 3.5);
    d_cut_LUsi-> SetPoint(8, 120, 3.5);
    d_cut_LUsi-> SetPoint(9, 80, 4.5);
    d_cut_LUsi-> SetPoint(10, 40, 6);
    d_cut_LUsi-> SetPoint(11, 0, 7.5);
    d_cut_LUsi-> SetPoint(12, 0, 5.5);
    auto t_cut_LUsi = new TCutG("t_cut_LUsi", 13);
    t_cut_LUsi-> SetPoint(0, 0, 7.5);
    t_cut_LUsi-> SetPoint(1, 40, 6);
    t_cut_LUsi-> SetPoint(2, 80, 4.5);
    t_cut_LUsi-> SetPoint(3, 120, 3.5);
    t_cut_LUsi-> SetPoint(4, 140, 3.5);
    t_cut_LUsi-> SetPoint(5, 160, 3.5);
    t_cut_LUsi-> SetPoint(6, 160, 5);
    t_cut_LUsi-> SetPoint(7, 140, 5);
    t_cut_LUsi-> SetPoint(8, 120, 5);
    t_cut_LUsi-> SetPoint(9, 80, 6);
    t_cut_LUsi-> SetPoint(10, 40, 7);
    t_cut_LUsi-> SetPoint(11, 0, 9);
    t_cut_LUsi-> SetPoint(12, 0, 7.5);

    // auto line = new TGraph();
    // for(int i=0; i<200; i++){
    //     line->AddPoint(i, -5./(160.)* i +12.);
    // }
    // auto line_LUsi = new TGraph();
    // for(int i=0; i<200; i++){
    //     line_LUsi->AddPoint(i, -7./(160.)* i +12.);
    // }
    // double luPIDdu = (2.-4.7)/(160.)*(csi)+7.7;
    
    int eventNum = treeIn -> GetEntries();
    for(int event=0; event<eventNum; event++){
        if(event%1000==0){cout << "event: " << event << " / " << eventNum << endl;}
        treeIn -> GetEntry(event);

        double dEdxWing[2];
        double CsIEWing[2];
        double SiWing[2];
        memset(dEdxWing, 0, sizeof(dEdxWing));
        memset(CsIEWing, 0, sizeof(CsIEWing));
        memset(SiWing, 0, sizeof(SiWing));

        double PID[2];
        double energy[2][2];
        memset(PID, 0, sizeof(PID));
        memset(energy, 0, sizeof(energy));

        for(int wing=0; wing<2; wing++){ // right, left
            double CsIDE = csiData[wing][0][0];
            double CsIUE = csiData[wing][1][0];
            double CsIDADC = csiData[wing][0][1];
            double CsIUADC = csiData[wing][1][1];

            if(siHitNum < 1){continue;}

            bool isUp = false;
            double csiEnergy[2];
            double siEnergy[2];
            memset(csiEnergy, 0, sizeof(csiEnergy));
            memset(siEnergy, 0, sizeof(siEnergy));
            if(siHitNum > 1){continue;}

            for(int i=0; i<siHitNum; i++){
                double adc = siData[wing][i][1][1];
                if(adc < 10){continue;}
                double siPosY = siData[wing][i][2][1];
                double omicADC = siData[wing][i][1][0];
                double juncADC = siData[wing][i][1][1];
                double omicE = siData[wing][i][0][0];
                double juncE = siData[wing][i][0][1];

                double siADC = omicADC + juncADC;
                double siE = omicE + juncE;

                if(siPosY < 0.){
                    hPID[wing][0] -> Fill(CsIDE, siE);
                    csiEnergy[i] = CsIDE;
                    siEnergy[i] = siE;
                    pidTotal -> Fill(CsIDE, siE);
                    if( p_cut->IsInside(CsIDE, siE) ){
                        PID[wing] = 1;
                    }
                    if( d_cut->IsInside(CsIDE, siE)){
                        PID[wing] = 2;
                    }
                    if( t_cut->IsInside(CsIDE, siE)){
                        PID[wing] = 3;
                    }
                    if(siE > -5./(160.)* CsIDE +12.){
                        PID[wing] = 4;
                    }
                    energy[wing][0] = CsIDE;
                    energy[wing][1] = siE;
                    
                }
                else{
                    hPID[wing][1] -> Fill(CsIUE, siE);
                    csiEnergy[i] = CsIUE;
                    siEnergy[i] = siE;
                    isUp = true;
                    pidTotal -> Fill(CsIUE, siE);
                    if( wing == 0){
                        if( p_cut->IsInside(CsIDE, siE) ){
                            PID[wing] = 1;
                        }
                        if( d_cut->IsInside(CsIDE, siE)){
                            PID[wing] = 2;
                        }
                        if( t_cut->IsInside(CsIDE, siE)){
                            PID[wing] = 3;
                        }
                        if(siE > -5./(160.)* CsIDE +12.){
                            PID[wing] = 4;
                        }
                        energy[wing][0] = CsIDE;
                        energy[wing][1] = siE;
                    }else{
                        if( p_cut_LUsi->IsInside(CsIDE, siE) ){
                            PID[wing] = 1;
                        }
                        if( d_cut_LUsi->IsInside(CsIDE, siE)){
                            PID[wing] = 2;
                        }
                        if( t_cut_LUsi->IsInside(CsIDE, siE)){
                            PID[wing] = 3;
                        }
                        if(siE > -7./(160.)* CsIDE +12.){
                            PID[wing] = 4;
                        }
                        energy[wing][0] = CsIDE;
                        energy[wing][1] = siE;
                    }
                    
                    
                }
            }       
        } 
        if( PID[0] == 1 ){
            if(energy[1][0] > 0 && energy[1][1] > 0){
                hNumSpec[0][0]-> Fill(energy[1][0],energy[1][1]);
            }
            
            // if( PID[1] == 1){
            //     hNumSpec[0][0] -> Fill(energy[1][0],energy[1][1]);
            // }
            // if( PID[1] == 2){
            //     hNumSpec[0][1] -> Fill(energy[1][0],energy[1][1]);
            // }
            // if( PID[1] == 3){
            //     hNumSpec[0][2] -> Fill(energy[1][0],energy[1][1]);
            // }
        }
        if( PID[1] == 1 ){
            if(energy[0][0] > 0 && energy[0][1] > 0){
                hNumSpec[1][0]-> Fill(energy[0][0],energy[0][1]);
            }
            // if( PID[0] == 1){
            //     hNumSpec[1][0] -> Fill(energy[0][0],energy[0][1]);
            // }
            // if( PID[0] == 2){
            //     hNumSpec[1][1] -> Fill(energy[0][0],energy[0][1]);
            // }
            // if( PID[0] == 3){
            //     hNumSpec[1][2] -> Fill(energy[0][0],energy[0][1]);
            // }
        }
        
        if(energy[0][0] > 0.1 && energy[0][1] > 0.1 && energy[1][1] > 0.1 && energy[1][1] > 0.1){
            hSpecTotal->Fill(energy[0][0], energy[0][1]);
            hSpecTotal->Fill(energy[1][0], energy[1][1]);
        }

        
    }
    TCanvas* cPID = new TCanvas("cPID","", 2000,2000);
    cPID -> Divide(2, 2);
    TCanvas* cPID2 = new TCanvas("cPID2","", 3000,2000);
    cPID2 -> Divide(3, 2);
    TCanvas* cPID3 = new TCanvas("cPID3","", 1000,1000);
    // int tmpIdx = 0;
    // for(int wing=0; wing<2; wing++){
    //     for(int du=0; du<2; du++){
    //         const int cidx = tmpIdx+1;
    //         cPID -> cd(cidx);
    //         cPID -> SetRightMargin(0.03);
    //         cPID -> SetTopMargin(0.03);
    //         gPad -> SetLogz();

    //         hPID[wing][du] -> SetStats(0);
    //         hPID[wing][du] -> SetTitle(Form("%s CsI-Si; CsI Energy [MeV]; Si Energy [MeV]", TypeName[tmpIdx].Data()));
    //         hPID[wing][du] -> Draw("colz");
    //         if(wing == 1 && du ==1 ){
    //             p_cut_LUsi->SetLineColor(kBlue);
    //             p_cut_LUsi->SetLineWidth(2);
    //             p_cut_LUsi->Draw("same");
    //             d_cut_LUsi->SetLineColor(kGreen);
    //             d_cut_LUsi->SetLineWidth(2);
    //             d_cut_LUsi->Draw("same");
    //             t_cut_LUsi->SetLineColor(kBlack);
    //             t_cut_LUsi->SetLineWidth(2);
    //             t_cut_LUsi->Draw("same");  
    //             // line_LUsi->Draw("same");
    //         }else{
    //             p_cut->SetLineColor(kBlue);
    //             p_cut->SetLineWidth(2);
    //             p_cut->Draw("same");

    //             d_cut->SetLineColor(kGreen);
    //             d_cut->SetLineWidth(2);
    //             d_cut->Draw("same");

    //             t_cut->SetLineColor(kBlack);
    //             t_cut->SetLineWidth(2);
    //             t_cut->Draw("same");

    //             // line->Draw("same");
    //         }
    //         tmpIdx++;
    //     }
    // }
    // cPID -> Draw();
    // cPID -> SaveAs("../csiPicture/PID_Energy_zoom.png");

    int tmpcId = 0;
    for(int wing=0; wing<2; wing++){
        for(int spec=0; spec<3; spec++){
            const int cidx = tmpcId + 1;
            cPID2->cd(cidx);
            hNumSpec[wing][spec]->Draw("colz");

            tmpcId++;
        }
        
    }
    cPID2-> Draw();
    cPID2-> SaveAs("../csiPicture/PID_spec.png");

    cPID3->cd();
    hSpecTotal->Draw("colz");
    cPID3->SaveAs("../csiPicture/PID_spec_Total.png");
}



void checkCSI()
{
    TString TypeName[4] = {"Right-Down", "Right-Up", "Left-Down", "Left-Up"};
    const double tpcADCThreshold = 50.;

    TFile* fileIn = new TFile("../data/SmallHIMACData_new.root","read");
    TTree* treeIn = (TTree*)fileIn -> Get("data");

    double trackFitData[2][2][2][2]; //[Right, Left][trkNum][xy, yz][slope, const]
    double trackData[2][2][8][32][4]; //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]

    double csiData[2][2][2]; //[Right, Left][Down, Up][Energy sumADC]
    double csiDataADC[2][2][512]; //[Right, Left][Down, Up][tb]
    int siHitNum;
    double siData[2][5][3][2]; //[Right, Left][hitNum][Energy, ADC, pos] Energy=(omic, junc), ADC=(omic, junc), pos=(x, y)

    double siDataBeamADCOmic[2][4][512]; // [Right, Left][omic][tb]
    double siDataBeamADCJunc[2][16][512]; // [Right, Left][junc][tb]
    double siDataBeamHitOmic[2][4];  // [Right, Left][omic][hit]
    double siDataBeamHitJunc[2][16]; // [Right, Left][junc][hit]

    treeIn -> SetBranchAddress("trackFitData", &trackFitData);
    treeIn -> SetBranchAddress("trkData", &trackData);
    treeIn -> SetBranchAddress("csiData", &csiData);
    treeIn -> SetBranchAddress("csiDataADC", &csiDataADC);
    treeIn -> SetBranchAddress("siHitNum", &siHitNum);
    treeIn -> SetBranchAddress("siData", &siData);
    treeIn -> SetBranchAddress("siBADCOmic", &siDataBeamADCOmic);
    treeIn -> SetBranchAddress("siBADCJunc", &siDataBeamADCJunc);
    treeIn -> SetBranchAddress("siBHitOmic", &siDataBeamHitOmic);
    treeIn -> SetBranchAddress("siBHitJunc", &siDataBeamHitJunc);


    TString csiType[4] = {"rd", "ru", "ld", "lu"};


    TCanvas* c1 = new TCanvas();
    c1 -> Divide(2,1);
    TCanvas* c3 = new TCanvas();

    TH2D* hPulse[2][2];
    TH1I* hMaxBin[2][2];
    TH1I* hMinBin[2][2];

    for(int wing=0; wing<2; wing++){
        for(int i=0; i<2; i++){
            hPulse[wing][i] = new TH2D(Form("hPulse%i_%i", wing, i), "", 500, 0.5, 500.5, 1200, -500, 700);
            hMaxBin[wing][i] = new TH1I(Form("hMaxBin%i_%i", wing, i),"",100, 150, 250);
            hMinBin[wing][i] = new TH1I(Form("hMinBin%i_%i", wing, i),"",100, 270, 370);
        } 
    }

    int eventNum = treeIn -> GetEntries();
    for(int event=0; event<eventNum; event++){
        if(event%1000==0){cout << "event: " << event << " / " << eventNum << endl;}
        treeIn -> GetEntry(event);
        
        for(int wing=0; wing<2; wing++){
            for(int i=0; i<2; i++){
                
                int maxBin = 0;
                double tmpMax = 0;
                int minBin = 0;
                double tmpMin = 999;

                for(int tb=1; tb<=500; tb++){
                    double adc = csiDataADC[wing][i][tb];   
                    // hPulse[wing][i] -> Fill(tb, adc);

                    if(tb > 170 && tb < 235){
                        if(tmpMax < adc){
                            tmpMax = adc;
                            maxBin = tb;
                        }
                    }
                    if(tb > 270 && tb < 350){
                        if(tmpMin > adc){
                            tmpMin = adc;
                            minBin = tb;
                        }
                    }
                }

                if(tmpMax < 50){continue;}
                if(tmpMin > -50){continue;}

                hMaxBin[wing][i] -> Fill(maxBin);
                hMinBin[wing][i] -> Fill(minBin);

                for(int tb=1; tb<=500; tb++){
                    double adc = csiDataADC[wing][i][tb];   
                    hPulse[wing][i] -> Fill(tb, adc);
                }

            }
        }
    }


    for(int wing=0; wing<2; wing++){
        for(int i=0; i<2; i++){
            c3 -> cd();
            hPulse[wing][i] -> SetStats(0);
            hPulse[wing][i] -> SetTitle(Form("CsI wing%i_du; tb; adc", wing, i));
            hPulse[wing][i] -> Draw("colz");
            c3 -> SetLogz();
            c3 -> SaveAs(Form("../csiPicture/MainPulse_%i_%i.png", wing, i));

            c1 -> cd(1);
            hMaxBin[wing][i] -> SetTitle(Form("CsI wing%i_ud%i;MaxBin [tb]", wing, i));
            hMaxBin[wing][i] -> Draw();

            c1 -> cd(2);
            hMinBin[wing][i] -> SetTitle(Form("CsI wing%i_ud%i;MinBin [tb]", wing, i));
            hMinBin[wing][i] -> Draw();

            c1 -> SaveAs(Form("../csiPicture/MainMaxMinBin_wing%i_ud%i.png", wing, i));

        }
    }

}

void checkADC()
{
    TString TypeName[4] = {"Right-Down", "Right-Up", "Left-Down", "Left-Up"};
    const double tpcADCThreshold = 50.;

    TFile* fileIn = new TFile("../data/SmallHIMACData.root","read");
    TTree* treeIn = (TTree*)fileIn -> Get("data");

    double trackFitData[2][2][2][2]; //[Right, Left][trkNum][xy, yz][slope, const]
    double trackData[2][2][8][32][4]; //[Right, Left][trkNum][yIdx][xIdx][x, y, z, adc]

    double csiData[2][2][2]; //[Right, Left][Down, Up][Energy sumADC]
    int siHitNum;
    double siData[2][5][3][2]; //[Right, Left][hitNum][Energy, ADC, pos] Energy=(omic, junc), ADC=(omic, junc), pos=(x, y)

    double siDataBeamADCOmic[2][4][512]; // [Right, Left][omic][tb]
    double siDataBeamADCJunc[2][16][512]; // [Right, Left][junc][tb]
    double siDataBeamHitOmic[2][4];  // [Right, Left][omic][hit]
    double siDataBeamHitJunc[2][16]; // [Right, Left][junc][hit]

    treeIn -> SetBranchAddress("trackFitData", &trackFitData);
    treeIn -> SetBranchAddress("trkData", &trackData);
    treeIn -> SetBranchAddress("csiData", &csiData);
    treeIn -> SetBranchAddress("siHitNum", &siHitNum);
    treeIn -> SetBranchAddress("siData", &siData);
    treeIn -> SetBranchAddress("siBADCOmic", &siDataBeamADCOmic);
    treeIn -> SetBranchAddress("siBADCJunc", &siDataBeamADCJunc);
    treeIn -> SetBranchAddress("siBHitOmic", &siDataBeamHitOmic);
    treeIn -> SetBranchAddress("siBHitJunc", &siDataBeamHitJunc);

    int eventNum = treeIn -> GetEntries();
    for(int event=0; event<eventNum; event++){
        if(event%1000==0){cout << "event: " << event << " / " << eventNum << endl;}
        treeIn -> GetEntry(event);

    }
    
}