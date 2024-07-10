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

using namespace TMath;
using namespace std;

// ========= ATTPC-constant ===============
const double fDriftvelocity = 55; //[mm/us]
const double fPadHeight = 11.9; //[mm]
const double fPadWeith = 1.9; //[mm]
const double fPadGap = 0.1; //[mm]
const double fDistance_TPC = 550; //[mm]
const double fDistance_Si = 715.1; //[mm]
const double fhalfH = 89; //[mm]
const double fTimeoffsetL = 259.6; //[mm]
const double fTimeoffsetR = 258.6; //[mm]
// 1timebucket = 0.6238 [mm]
const vector<double> fAngle;
// ============ Si-constant ===============
const double fSi_dy = 5.0375; // [mm]
const double fGap_ud = 7.55; //[mm]
const double fSi_y = 40.3; //[mm]
const double fcenterX_ = 9.375; //[mm]
const double fSi_Cdtocenter = fSi_y + fGap_ud; //[mm]
//=========================================
enum adcOrTb
{
    adc,
    tb
};
struct parsSi_xy{
    double x[5];
    double y[5];
}parsSi_xy; // Si-position x, y
struct TPC_hough{
    double rho[5];
    double angle[5];
}hough;
struct TPC_YZ{
    double s[5];
    double c[5];
}tpc_yz; // slope and const of the trajectory
double point_eval(double x, double resultangle, double resultrho){
    return ( - x + sin(TMath::DegToRad()*(resultangle-180))*resultrho)/tan(TMath::DegToRad()*(resultangle-90)) + cos(TMath::DegToRad()*resultangle)*resultrho;
    //returns x position.
    //hough histogram -> line component. 
}
double houghTrans(double x, double y, double theta){ //transforming function
    theta *= TMath::DegToRad();
    return x*cos(theta)+ y*sin(theta);
}
double houghTransformation(int event, vector<double> angles, int tracks, vector<double> point_x, vector<double> point_y, vector<double> point_adc){
    auto c1 = new TCanvas("c1","", 1650, 1500);
    double dtheta = 0.1;
    double low_degree = 120.;
    double high_degree = 240.;

    //hough transform
    vector<double> transformed_point_xy;
    for(int i = 0; i < point_x.size(); i++){
        for(int j=0; j < angles.size(); j++){
            transformed_point_xy.push_back( houghTrans( point_x[i], point_y[i], angles[j]) );
            // cout<<houghTrans( point_x[i], point_y[i], angles[j])<<endl;
            //hough transform
        }   
    }
    auto houghHist_xy =new TH2D("houghHist_xy","houghHist_xy projection",
                                        (high_degree-low_degree) / dtheta,
                                        low_degree,
                                        high_degree,
                                        (high_degree-low_degree) / dtheta,
                                        0,
                                        100 );
    for (int i = 0; i < transformed_point_xy.size(); i++){
        houghHist_xy->Fill( angles[i%angles.size()], -transformed_point_xy[i],point_adc[i/angles.size()] ); // (x, y, w)
        //fill hough points into histogram ( angle-rho histogram )
    }
    double hough_sigma = 100;
    double result_rho_xy = -999;
    double result_angle_xy = -999;

    // find intersection point 
    for (int i = 0; i<(high_degree-low_degree)/dtheta; ++i){
        auto find_sigma = (TH1D*) houghHist_xy->ProjectionY("",i-1, i+1);
        if(find_sigma->GetStdDev() < hough_sigma){
            hough_sigma = find_sigma->GetStdDev();
            result_rho_xy = find_sigma->GetMean();
            result_angle_xy = i;
        }
        find_sigma->Reset(); 
    }
    result_angle_xy += low_degree/dtheta;
    result_angle_xy *= dtheta;

    c1->cd();
    houghHist_xy->SetStats(0);
    houghHist_xy->SetTitle(" ;angle  [#circ]; #rho  [mm]");
    houghHist_xy->Draw("colz");
    houghHist_xy -> GetYaxis()->SetTitleSize(0.04);
    houghHist_xy -> GetYaxis()->SetTitleOffset(1.1);
    houghHist_xy -> GetXaxis()->SetTitleSize(0.04);
    houghHist_xy -> GetXaxis()->SetTitleOffset(0.9);
    houghHist_xy -> GetZaxis()->SetTitle("Entries");
    houghHist_xy -> GetZaxis() -> SetTitleOffset(1.8);
    // point->Draw("p same");

    auto latex = new TLatex();
    latex -> SetTextSize(0.085);
    // latex -> DrawLatexNDC(0.15, 0.87, "Hough-Trans Phase");
    latex->DrawLatexNDC(0.15, 0.9, "(b)");
    c1->SetTopMargin(0.03);
    c1->SetRightMargin(0.18);
    c1->Draw();
    c1->SaveAs("../result/hough_paper.png");
    c1->SaveAs("../result/hough_paper.pdf");
    hough.rho[tracks] = result_rho_xy;
    hough.angle[tracks] = 180 - result_angle_xy;
    houghHist_xy->Delete();
    // cout<<hough.angle[tracks]<< " : "<< hough.rho[tracks]<<endl;
    return 0;
}

double getTPCSlope(int event, int tracks, vector<double> point_y, vector<double> point_z, vector<double> point_adc){
    auto trackYZ = new TH2D("tracksYZ", "trackYZ", 8, -6.25, 93.75, 700, 0, 700); // YZ pad plane
    auto fTrackYZ = new TF1("fTrackYZ", "pol1", -6.25, 93.75);
    double pars[2];

    for(int i=0; i< point_y.size(); i++){

        trackYZ-> Fill(point_y[i] , point_z[i], point_adc[i]);    
    }
    trackYZ->Fit("fTrackYZ", "QNR");
    fTrackYZ->GetParameters(pars);
    tpc_yz.s[tracks] = pars[1];
    tpc_yz.c[tracks] = pars[0];
    
    trackYZ->Delete();
    return 0;
}
double getSiTrack(int event, int wing, dataStructure* data){ // wing = 0 -> Rightwing, 1 -> Leftwing
    auto c1 = new TCanvas("c1", "", 1000, 1000);
    c1->SetRightMargin(0.16);
    c1->SetTopMargin(0.02);
    GETSiPad siPad;
    int index; // Index 0-RD, 1-RU, 2-LD, 3-LU
    

    if(wing == 0){
        index = 0;
    }else{
        index = 2;
    }
    // Si
    double dis[4] = {0., 6., 6., 0,};
    auto juncPoint = new TGraph();
    int nHits_SO = 0;
    int nPos_SJ = 0;
    int adcSi_ohmic[2][4] = { 0 }; //[Up or Down][ohmicIdx]
    int adcSi_junc[2][16] = { 0 }; //[Up or Down][ohmicIdx]
    double adcSi_junc_U_even[8] = { 0 };
    double adcSi_junc_U_odd[8] = { 0 };
    double adcSi_junc_D_even[8] = { 0 };
    double adcSi_junc_D_odd[8] = { 0 };
    double centerOmicU[4] = { 0 };
    double centerOmicD[4] = { 0 };

    for(int siIndex = index; siIndex<index+2; siIndex++){ 
        //Right wing is inversed Pad
        if(siIndex == index){ //Down
            for(int ohmicIdx=0; ohmicIdx<4; ohmicIdx++){ 
                adcSi_ohmic[siIndex - index][ohmicIdx] = double(data -> getSiOmicHits(siIndex, ohmicIdx)); //adc
                siPad.ADC->Fill(fcenterX_ + 18.7 * ohmicIdx, fSi_y / 2 , adcSi_ohmic[siIndex - index][ohmicIdx]);  // For pad drawing
                if(adcSi_ohmic[siIndex - index][ohmicIdx] > 0){
                    nHits_SO ++;
                    centerOmicD[ohmicIdx] = 1;
                }
            }
            int i = 0;    
            for(int juncIdx=0; juncIdx<16; juncIdx++){ 
                adcSi_junc[siIndex - index][juncIdx] = double(data -> getSiJuncHits(siIndex, juncIdx)); //adc
                if(juncIdx % 2 == 0){
                    adcSi_junc_D_even[i] = adcSi_junc[siIndex - index][juncIdx];
                }else{
                    adcSi_junc_D_odd[i] = adcSi_junc[siIndex - index][juncIdx];
                    i++;
                }
            }
            
        }else{ //Up
            for(int ohmicIdx=0; ohmicIdx<4; ohmicIdx++){ 
                adcSi_ohmic[siIndex - index][ohmicIdx] = double(data -> getSiOmicHits(siIndex, ohmicIdx)); //adc
                siPad.ADC->Fill(fcenterX_ + 18.7 * ohmicIdx, 1.5 * fSi_y + fGap_ud, adcSi_ohmic[siIndex - index][ohmicIdx]);  // For pad drawing
                if(adcSi_ohmic[siIndex - index][ohmicIdx] > 0){
                    nHits_SO ++;
                    centerOmicU[ohmicIdx] = 1;
                }
            }
            int i = 0;    
            for(int juncIdx=0; juncIdx<16; juncIdx++){ 
                adcSi_junc[siIndex - index][juncIdx] = double(data -> getSiJuncHits(siIndex, juncIdx)); //adc
                if(juncIdx % 2 == 0){
                    adcSi_junc_U_even[i]=adcSi_junc[siIndex - index][juncIdx];
                }else{
                    adcSi_junc_U_odd[i]=adcSi_junc[siIndex - index][juncIdx];
                    i++;
                }
            }
        }
    }
    double y_positionD = fSi_dy / 2; // first strip Y position
    double y_positionU = fSi_dy / 2 + fSi_y + 2 * fGap_ud; // first strip Y position of Up Si
    double omicpositionUX[4] = { 0 }; // index, adc
    double omicpositionDX[4] = { 0 };
    for(int i=0; i<4; i++){
        if(centerOmicU[i] == 1){
            omicpositionUX[i] = fcenterX_ + 2 * i * fcenterX_; //Up-Omic channel 1-4
        }
        if(centerOmicD[i] == 1){
            omicpositionDX[i] = fcenterX_ + 2 * i * fcenterX_; //Down-Omic channel 1-4
        }
    }

    for(int i=0 ; i<8; i++){ //strip goes to up
        if(adcSi_junc_U_even[i]!= 0 && adcSi_junc_U_odd[i]!= 0){
            double juncposition = 75 * (1 - (adcSi_junc_U_even[i]/(adcSi_junc_U_even[i] + adcSi_junc_U_odd[i]))); // find the X-position of up-Si
            for(int j=0; j<4; j++){
                if(juncposition > omicpositionUX[j] - fcenterX_ - dis[j] && juncposition < omicpositionUX[j] + fcenterX_ + dis[j+1]){
                    juncPoint->SetPoint(juncPoint->GetN(), juncposition, y_positionU);
                    nPos_SJ ++;

                    adcSi_ohmic[index+1][j];
                }
            }
        }
        if(adcSi_junc_D_even[i]!= 0 && adcSi_junc_D_odd[i]!= 0){
            double juncposition = 75 * (1 - (adcSi_junc_D_even[i]/(adcSi_junc_D_even[i] + adcSi_junc_D_odd[i]))); // find the X-position of down-Si
            for(int j=0; j<4; j++){
                if(juncposition > omicpositionDX[j] - fcenterX_ - dis[j] && juncposition < omicpositionDX[j] + fcenterX_ + dis[j+1]){
                    juncPoint->SetPoint(juncPoint->GetN(), juncposition, y_positionD);
                    nPos_SJ ++;

                    adcSi_ohmic[index][j];
                }
            }
        }
        y_positionU += fSi_dy;
        y_positionD += fSi_dy;
    }
    // if(nPos_SJ > 0){ // threshold- at least 1 omic-hit and 2 junc-hit
    for(int i=0; i<juncPoint->GetN(); i++){
        parsSi_xy.x[i] = juncPoint->GetPointX(i);
        parsSi_xy.y[i] = juncPoint->GetPointY(i); 
    }
    c1->cd();                    // for the draw Si detector
    siPad.ADC->GetZaxis()->SetTitle("Ohmic ADC");
    siPad.ADC->GetZaxis()->SetTitleOffset(1.7);
    siPad.ADC->SetStats(0);
    siPad.ADC->Draw("colz");
    siPad.ADC->SetMinimum(1);
    siPad.ADC->SetMaximum(1000);
    siPad.Frame->Draw("same");
    juncPoint->SetMarkerStyle(104);
    juncPoint->SetMarkerSize(5);
    juncPoint->Draw("same p");
    c1->SaveAs(Form("../SiADC/siPad_%d.pdf", event));
    siPad.Clear();
    return juncPoint->GetN(); // return num of juncPoint = num of tracks estimated in Si  
}
double getTargetX(int wing, int tracks, double s, double x){
    if(wing == 0){
        return (- fDistance_Si * s + ( x - 37.5 )) / (Cos(40*Pi()/180) - s * Sin(40*Pi()/180));
    }else{
        return ( -fDistance_Si * s - ( x - 37.5 )) / (Cos(40*Pi()/180) + s * Sin(40*Pi()/180));
    }
    
}
double getTargetZ(int wing, int tracks, double s, double c){
    if(wing == 0){
        // return (-fDistance_Si * s + ( x - 37.5 )) / (Sin(40*Pi()/180) + s *Cos(40*Pi()/180)); 
        return (-(fDistance_TPC-(fPadHeight+fPadGap)*3.5)*s + c - fhalfH - fTimeoffsetR)/(Sin(40*Pi()/180)+ s * Cos(40*Pi()/180));
    }else{
        // return (-fDistance_Si * s - ( x - 37.5 )) /(-Sin(40*Pi()/180) + s *Cos(40*Pi()/180)); 
        return -((fDistance_TPC+(fPadHeight+fPadGap)*3.5)*s + c - fhalfH - fTimeoffsetL)/(-Sin(40*Pi()/180)+ s * Cos(40*Pi()/180));
    }
    
}
void dataAnalysis();

int main() {
    dataAnalysis();
}

void dataAnalysis() {
    GETDecoder decoder;
    GETPad pad;
    
    // auto c1 = new TCanvas("c1", "", 1000, 1000);
    // c1->Divide(2,1);
    // auto c2 = new TCanvas("c1", "", 800, 800);
    auto cPad = new TCanvas("cPad", "", 1650, 1500);
    cPad->SetTopMargin(0.03);
    cPad->SetRightMargin(0.16);
    // cPad->Divide(2,1);
    auto latex = new TLatex();
    latex->SetTextSize(0.085);

    double dtheta = 0.1;
    double low_degree = 120.;
    double high_degree = 240.;
    vector<double> angles;
    for(int i = 0; i < (high_degree-low_degree)/dtheta; i ++){
        double angle = low_degree + 0.1 * i;
        angles.push_back(angle);  
        //give calculation info to every delta degrees
    }
    auto fGaus = new TF1("fGaus", "gaus", 0, 31);               // The cluster of x_pad direction

    auto hYZR_transform_dx = new TH1D("hYZR_transform_dx-DT", " ;#it{#Delta X} [mm]", 100, -200, 100);
    auto hYZR_transform_dz = new TH1D("hYZR_transform_dz-DT", " ;#Delta Z [mm]", 100, -200, 100); hYZR_transform_dz->SetLineColor(kBlue);
    auto hYZL_transform_dx = new TH1D("hYZL_transform_dx-DT", " ;#Delta X [mm]", 100, -100, 200);
    auto hYZL_transform_dz = new TH1D("hYZL_transform_dz-DT", " ;#Delta Z [mm]", 100, -200, 100); hYZL_transform_dz->SetLineColor(kRed);
    auto h2dz = new TH2D("h2dz", " ;TPC-Right  #Delta Z [mm];TPC-Left  #Delta Z [mm]", 30, -200, 100, 30, -200, 100);
    // auto ftarget = new TF1("ftarget", "gaus", -50, 50);
    TLegend* sc = new TLegend(0.15, 0.7, 0.4, 0.9); sc->SetFillStyle(0); sc->SetBorderSize(0);

    //data taking
    TFile* fileIn = new TFile("../data/treeOfData_mainRun_temp.root", "read");
    // TFile* fileIn = new TFile("/home/shlee/workspace/forHIMAC/data/treeOfData_C12Run.root", "read");
    TTree* treeIn = (TTree*)fileIn -> Get("event");
    TClonesArray* clonesArray = new TClonesArray("dataStructure");
    treeIn -> SetBranchAddress("dataStructure", &clonesArray);

    auto trackYZ = new TH2D("tracksYZ", "time bucket; y_{pad} [mm]; timebucket [20 ns]", 8, -6.25, 93.75, 400, 0, 400);
    auto fTracksYZ = new TF1("fTrackYZ", "pol1", -6.25, 93.75);
    
    int numSi = 0;
    int eventNum = treeIn -> GetEntries();
    for(int event = 102; event < 103; event++){
        if(event%1000==0){cout << "event: " << event<< " / "<< eventNum << endl;}
        treeIn -> GetEntry(event);
        dataStructure* data = (dataStructure*)clonesArray -> At(0);
        // ============= general info ====================
        int runNumber = data -> getRunID();
        int eventNumber = data -> getEventID();
        bool isBrokenEvents = data -> getBrokenEvent();
        if(isBrokenEvents){ continue; }

        double dzR[10] = { 0 };
        double dzL[10] = { 0 };
        int numR = 0;
        int numL = 0;
        for(int wing = 0; wing < 2; wing++){ // 0 = right wing, else = left wing
            // ======================== tpc data ========================
            // ATTPC
            double adcPad[8][32] = { 0 };
            double timePad[8][32] = { 0 };
            
            // hough transform
            vector<vector<double>> point_x, point_y, point_z, point_adc;
            vector<double> point_x0, point_y0, point_z0, point_adc0;
            vector<double> point_x1, point_y1, point_z1, point_adc1;

            // ============== event-driven method ====================
            // auto hCluster = new TH1D("hCluseter","", 32, -0.5, 31.5);
            // double xPos, yPos;
            // auto trackXY = new TGraph();
            // =======================================================
            for(int yId = 1; yId < 7; yId++){
                for(int xId = 0; xId < 32; xId++){
                    if(wing == 0){
                        adcPad[yId][xId] =(double) data -> getRightTPCHits(yId, xId, 0);
                        timePad[yId][xId] =(double) data -> getRightTPCHits(yId, xId, 1);
                    }else{
                        adcPad[yId][xId] =(double) data -> getLeftTPCHits(yId, xId, 0);
                        timePad[yId][xId] =(double) data -> getLeftTPCHits(yId, xId, 1);
                    }
                    
                    if(adcPad[yId][xId] >= 50) {
                        point_adc0.push_back( adcPad[yId][xId] );
                        point_x0.push_back( 2 * xId );
                        point_y0.push_back( 12 * yId );
                        point_z0.push_back( 1.2 * timePad[yId][xId] );
                        pad.ADC->Fill( 2. * xId, 12. * yId, adcPad[yId][xId]); // For pad drawing
                    }
                }           
            }
            // =======================separate the track==========================
            if(point_x0.size() == 0) continue ;
            int standard_x = point_x0.front();
            for(int i = 0; i < point_x0.size(); i++){
                if( (point_x0[i] - standard_x) < -6 || (point_x0[i] - standard_x) > 14){
                    point_x1.push_back(point_x0[i]);
                    point_y1.push_back(point_y0[i]);
                    point_z1.push_back(point_z0[i]);
                    point_adc1.push_back(point_adc0[i]);
                    
                    point_x0.erase(point_x0.begin() + i);
                    point_y0.erase(point_y0.begin() + i);
                    point_z0.erase(point_z0.begin() + i);
                    point_adc0.erase(point_adc0.begin() + i);
                    i--;
                }
            }
            if(point_x0.size()>17){ 
                point_x.push_back(point_x0); 
                point_y.push_back(point_y0);
                point_z.push_back(point_z0);
                point_adc.push_back(point_adc0);
            }
            if(point_x1.size()>17){ 
                point_x.push_back(point_x1); 
                point_y.push_back(point_y1); 
                point_z.push_back(point_z1); 
                point_adc.push_back(point_adc1);
            }
            // =======================separate the track==========================
            // houghTransformation(event, angles, 0, point_x0, point_y0, point_adc0);
            for(int tracks=0; tracks< point_adc.size(); tracks++){
                houghTransformation(event, angles, tracks, point_x[tracks], point_y[tracks], point_adc[tracks]);
                getTPCSlope(event, tracks, point_y[tracks], point_z[tracks], point_adc[tracks]);   
            }
            TF1* fTracksXY[2];
            int siPointN = getSiTrack(event, wing, data);
            if(point_adc.size() > 0) {

                // =============================Drawing TPC-Pad ==================================
                cPad->cd(); // tpc R    // for draw the TPC-XY pad                
                pad.ADC->SetStats(0);
                pad.ADC->GetYaxis() -> SetTitleSize(0.04);
                pad.ADC->GetXaxis() -> SetTitleSize(0.04);
                pad.ADC->GetZaxis() -> SetTitle("ADC");
                pad.ADC->GetZaxis() -> SetTitleOffset(1.5);
                pad.ADC->SetMinimum(1);
                pad.ADC->SetMaximum(3500);
                pad.ADC->Draw("colz");
                pad.Frame->Draw("same");
                latex->DrawLatexNDC(0.15, 0.9, "(a)");
                // =============================Drawing TPC-Pad ==================================
                for(int i=0; i< point_adc.size(); i++){
                    fTracksXY[i] = new TF1(Form("fTrackYX_%i", i), "pol1", -6.25, 93.75);
                    fTracksXY[i]->SetParameters( - hough.rho[i] / tan((hough.angle[i] - 180) *Pi()/180) , 1 / tan((hough.angle[i] - 180)*(Pi()/180)) );
                    double dz_genCoord_drift = getTargetZ(wing, i, tpc_yz.s[i], tpc_yz.c[i]);
                    if(wing == 0){
                        dzR[numR] = dz_genCoord_drift;
                        hYZR_transform_dz->Fill(dz_genCoord_drift);
                        numR++;
                    }else{
                        dzL[numL] = dz_genCoord_drift;
                        hYZL_transform_dz->Fill(dz_genCoord_drift);
                        numL++;
                    }
                    // =============================Drawing TPC-Pad ==================================
                    fTracksXY[i]->SetLineWidth(3.);
                    fTracksXY[i]->SetLineColor(kRed);
                    fTracksXY[i]->SetLineWidth(4);
                    fTracksXY[i]->Draw("same"); // XY pad plane
                    // =============================Drawing TPC-Pad ==================================
                    // ===========================================================================
                    trackYZ->Reset("ICESM");
                    fTracksYZ->SetParameters(0, 0);
                    for(int j=0; j< point_y0.size(); j++){
                        trackYZ-> Fill(point_y[i][j] , point_z[i][j], point_adc[i][j]);    
                    }
                    // cPad->cd(2);
                    // fTracksYZ->SetParameters( tpc_yz.c[0] , tpc_yz.s[0] );
                    // trackYZ->Draw("colz");
                    // fTracksYZ->Draw("same"); 
                    // auto latex = new TLatex();
                    // latex -> SetTextSize(0.06);
                    // latex -> DrawLatexNDC(0.1, 0.4, Form("const : %f", tpc_yz.c[0]));
                    // latex -> DrawLatexNDC(0.1, 0.2, Form("slope : %f", tpc_yz.s[0]));
                    cPad->SaveAs(Form("../result/TPC_Si_%d.pdf", event)); 
                    // ===========================================================================

                    // for(int j=0; j<siPointN; j++){
                    //     double dx_genCoord_drift = getTargetX(wing, i, tpc_yz.s[i], parsSi_xy.x[j]);
                    //     double dz_genCoord_drift = getTargetZ(wing, i, tpc_yz.s[i], parsSi_xy.x[j]);
                    //     if(wing == 0){
                    //         dzR[numR] = dz_genCoord_drift;
                    //         hYZR_transform_dx->Fill(dx_genCoord_drift);
                    //         hYZR_transform_dz->Fill(dz_genCoord_drift);

                    //         numR++;
                    //     }else{
                    //         dzL[numL] = dz_genCoord_drift;
                    //         hYZL_transform_dx->Fill(dx_genCoord_drift);
                    //         hYZL_transform_dz->Fill(dz_genCoord_drift);

                    //         numL++;
                    //     }
                    // }
                }
                // cPad->SetTopMargin(0.03);
                // cPad->SetRightMargin(0.03);
                // cPad->SaveAs(Form("../result/TPC_Si_before%d.png", event));   
                // cPad->SaveAs(Form("../result/TPC_Si_before%d.pdf", event));   
                
            }
              
            pad.Clear(); 
            
        }
        if(numR > 0 && numL > 0){
            numSi ++;
            for(int i=0; i < numR; i++){
                for(int j=0; j < numL; j++){
                    h2dz->Fill(dzR[i], dzL[j]);
                }
            }
        }
    }
    // hYZR_transform_dz->SetLineWidth(2);
    // hYZL_transform_dz->SetLineWidth(2);
    // hYZR_transform_dz->SetStats(0);
    // hYZL_transform_dz->SetStats(0);
    // h2dz->SetStats(0);
    // h2dz->GetYaxis()->SetTitleOffset(1.3);
    // sc->AddEntry(hYZR_transform_dz,"TPC-Right","l");
    // sc->AddEntry(hYZL_transform_dz,"TPC-Left","l");
    // sc->SetEntrySeparation(0.4);
    // c1->cd(1);
    
    // hYZL_transform_dz->Draw();
    // hYZR_transform_dz->Draw("same");
    // sc->Draw("same");
    // c1->cd(2);
    // h2dz->Draw("colz");
    // gPad->SetLogz();
    // hYZR_transform_dz->Draw();
    // c1->cd(3);
    // hYZL_transform_dx->Draw();
    // c1->cd(4);
    // hYZL_transform_dz->Draw();
    // c1->SaveAs("../result/after_hough/hYZ_transform_dxdydz_v=5.7_onlyTPC.png");
    // c1->SaveAs("../result/after_hough/hYZ_transform_dxdydz_v=5_onlyTPC.pdf");
    // cout<<numSi<<endl;
    // c2->cd();
    // hTPCSi->Draw();
    // // hdiff->Draw();
    // c2->SaveAs("../result/after_hough/hHoughAngle_L.png");
}
