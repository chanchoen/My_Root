#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "TSystem.h"
#include "GETDecoder.hh"
#include "GETPad.hh"
#include "TObject.h"
#include "dataStructure.hh"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLatex.h"

using namespace TMath;
using namespace std;

#define nAsAd 4


const Int_t runList[8] = {0, 1, 2, 3, 4, 5, 6, 25};

const Int_t thresholdOmic[6] = {50, 40, 30, 40, 100, 90};
const Int_t thresholdJunc[6] = {40, 40, 30, 30, 60, 70};
const Int_t thresholdCsI = 14;

double arrSigFcnTR[1000], arrSigFcnTL[1000];

double sigFcnTR(double* x, double* par) {
    int idx = x[0] - (par[0] - 149);
    if (idx < 0) {
        return par[1] * arrSigFcnTR[0];
    } else if (idx >= 350) {
        return par[1] * arrSigFcnTR[349];
    } else {
        return par[1] * arrSigFcnTR[idx];
    }
}
double sigFcnTL(double* x, double* par) {
    int idx = x[0] - (par[0] - 149);
    if (idx < 0) {
        return par[1] * arrSigFcnTL[0];
    } else if (idx >= 350) {
        return par[1] * arrSigFcnTL[349];
    } else {
        return par[1] * arrSigFcnTL[idx];
    }
}
double fitFcnTR(double* x, double* par) {
    return sigFcnTR(x, par) < par[2] ? sigFcnTR(x, par) : par[2];
}
double fitFcnTL(double* x, double* par) {
    return sigFcnTL(x, par) < par[2] ? sigFcnTL(x, par) : par[2];
}


void readSigFcn() {
    TFile* fileIn = new TFile("../source/hSignalTemplate.root", "read");
    auto hSignal1DTemplateTR = (TH1D*)fileIn->Get("hSignal1DTemplateTR");
    auto hSignal1DTemplateTL = (TH1D*)fileIn->Get("hSignal1DTemplateTL");
    
    for (int buckId = 1; buckId <= 350; buckId++) {
        arrSigFcnTR[buckId - 1] = hSignal1DTemplateTR->GetBinContent(buckId);
        arrSigFcnTL[buckId - 1] = hSignal1DTemplateTL->GetBinContent(buckId);
        arrSigFcnTR[buckId - 1] /= 1000.;
        arrSigFcnTL[buckId - 1] /= 1000.;
    }
    fileIn->Close();
}

void grawToRoot();
void pulsegenerate();
std::string getLabel(int asadID, int agetID, int chanID);
int mode(int array[], int size);

int main() 
{
    // pulsegenerate();
    grawToRoot();
}

void grawToRoot()
{
    GETPad pad;

    TFile* fileOut = new TFile("../data/treeOfData_mainRun_temp.root", "recreate");
    TTree* treeOut = new TTree("event","event");

    TClonesArray* outClonesArray = new TClonesArray("dataStructure");
    treeOut -> Branch("dataStructure", &outClonesArray);
    dataStructure* data;
    

    readSigFcn();

    int signal[512];
    bool isBroken;

    TCanvas* cc = new TCanvas();
    TLatex* latexTmp = new TLatex();

    auto hSignal = new TH1D("hSignal", "", 500, 0.5, 500.5);
    auto hPedestal = new TH1D("hPedestal", "", 1000, -499.5, 500.5);
    TF1* fFitFcnTR = new TF1("fFitFcnTR", fitFcnTR, 1, 500, 3);
    TF1* fFitFcnTL = new TF1("fFitFcnTL", fitFcnTL, 1, 500, 3);
    
    for (int runID = 0; runID < 3; runID++) {
        GETDecoder decoder[nAsAd];

        for (int asadID = 0; asadID < nAsAd; asadID++){
            // gSystem->Exec(Form("ls /home/shlee/workspace/forHIMAC/data/carbon_200MeV/run%02d/CoBo0_AsAd%d* > ../data/fileList%02d.txt", runList[runID], asadID, asadID));
	    gSystem->Exec(Form("ls /data/public_data/HIMACAT-TPC/HIMACData/raw/run%02d/CoBo0_AsAd%d* > ../data/fileList%02d.txt", runList[runID], asadID, asadID));
            decoder[asadID].OpenFromList(Form("../data/fileList%02d.txt", asadID));
	    cout<<"asadID : "<<asadID<<endl;
        }

        int tmpEventNum = 0;
        while (decoder[1].Run()) {
            if(tmpEventNum%100==0){cout << "processing... (runID: "<< runList[runID] << " | EventID: " << tmpEventNum << ") " << endl;}
            tmpEventNum++;
            // if(tmpEventNum==101){break;}

            outClonesArray -> Clear("C");
            data = (dataStructure*)outClonesArray -> ConstructedAt(0);

            // Event info
            for (int asadID = 0; asadID < 4; asadID++) {
                if (asadID != 1) decoder[asadID].Run();
                if (asadID == 1) {
                    data -> setRunID(runList[runID]);
                    data -> setEventID(int(decoder[asadID].EventID));
                    data -> setEventTime(double(decoder[asadID].EventTime));
                    data -> setEventDiffTime(double(decoder[asadID].DiffTime));
                    isBroken = false;
                }
            }
            
            // ================================== CsI production ================================== 
            for (int asadID = 1; asadID <= 3; asadID += 2) {
                int agetID = 3;
                int chanIDs[2][2] = {{49, 15}, {66, 32}};  // {{RD, RU}, {LD, LU}}
                for (int chanID : chanIDs[(asadID - 1) / 2]) {
                    int FPNChanID = decoder[asadID].FPNChanID(chanID);
                    std::fill_n(signal, nBUCK, -9999);
                    std::string label = getLabel(asadID, agetID, chanID);
                    int csiIdx = -1;
                    if(label=="crd"){csiIdx = CsI_RD;}
                    if(label=="cru"){csiIdx = CsI_RU;}
                    if(label=="cld"){csiIdx = CsI_LD;}
                    if(label=="clu"){csiIdx = CsI_LU;}

                    for (int buckID = 1; buckID <= 500; buckID++) {
                        int ADC = decoder[asadID].ADC[agetID][chanID][buckID];
                        int FPN = decoder[asadID].ADC[agetID][FPNChanID][buckID];
                        signal[buckID] = ADC - FPN;
                        // For finding a broken event
                        if (buckID >= 2) {
                            int oldADC = decoder[asadID].ADC[agetID][chanID][buckID - 1];
                            if (std::fabs(ADC - oldADC) > 50) {
                                isBroken = true;
                                break;  
                            }
                        }
                        if (buckID <= 140){hPedestal->Fill(ADC - FPN);}
                    }
                    if (isBroken) break;                    

                    int ped = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());

                    int maxADC = 0;
                    for(int buckID = 1; buckID <= 500; buckID++) {
                        int adc = signal[buckID]-ped;
                        data -> setCsIADC(csiIdx, buckID, adc);

                        if(buckID > 150 && buckID < 220){
                            if(maxADC < adc){
                                maxADC = adc;
                            }
                        }
                    }
                    hPedestal->Reset("ICESM");
                    if(maxADC < thresholdCsI){maxADC = 0;}
                    data -> setCsIHits(csiIdx, maxADC);
                }
                if (isBroken) break;  // Break AsAdID loop
            }
            if (isBroken) continue;  // Skip this event

            // ================================== AT-TPC production ================================== 
            for (int asadID = 0; asadID <= 2; asadID += 2) {
                for (int agetID = 0; agetID < 4; agetID++) {
                    for (int chanID = 0; chanID < 68; chanID++) {
                        if (decoder[asadID].IsFPN(chanID) || decoder[asadID].IsDead(chanID)) continue;
                        int FPNChanID = decoder[asadID].FPNChanID(chanID);
                        for (int buckID = 1; buckID <= 500; buckID++) {
                            int ADC = decoder[asadID].ADC[agetID][chanID][buckID]-decoder[asadID].ADC[agetID][FPNChanID][buckID];

                            if (buckID <= 150) {
                                hPedestal->Fill(ADC);
                            } else {
                                hSignal->SetBinContent(buckID, ADC);
                            }
                        }

                        int pedADC = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());
                        int maxBin = hSignal->GetMaximumBin();
                        int maxADC = hSignal->GetMaximum() - pedADC;
                        // cout<<maxADC<<endl;
                        hSignal->Reset("ICESM");
                        hPedestal->Reset("ICESM");

                        for(int buckID = 1; buckID <= 500; buckID++) {
                            int adc = signal[buckID]-pedADC;
                            if(asadID==0){data->setRightTPCADC(pad.GetYId(agetID, chanID), pad.GetXId(agetID, chanID), signal[buckID], buckID);}
                            if(asadID==2){data->setLeftTPCADC(pad.GetYId(agetID, chanID), pad.GetXId(agetID, chanID),  signal[buckID], buckID);}
                        }

                        if (150 < maxBin && maxBin <= 400 && 50 < maxADC) {
                            for (int buckID = 1; buckID <= 500; buckID++) {
                                int ADC = decoder[asadID].ADC[agetID][chanID][buckID]-decoder[asadID].ADC[agetID][FPNChanID][buckID];
                                hSignal->SetBinContent(buckID, ADC - pedADC);
                            }
                            if (asadID == 0){                                
                                fFitFcnTR -> SetParameters(double(maxBin)-0.5, maxADC, double(maxADC+50.));
                                fFitFcnTR -> SetParLimits(0, double(maxBin)-20.5, double(maxBin)+20.5);
                                fFitFcnTR -> SetParLimits(1, (maxADC) - 50, (maxADC) + 50);
                                hSignal->Fit(fFitFcnTR, "QNRSM");
                                double parsTR[3];
                                fFitFcnTR->GetParameters(parsTR);
                                // if(maxADC > 3500){
                                //     cc -> cd();
                                //     hSignal -> Draw();
                                //     fFitFcnTR -> Draw("same");
                                    

                                //     latexTmp -> Clear();
                                //     latexTmp -> DrawLatexNDC(0.15, 0.8, Form("tb_hist = %.2f", double(maxBin)));
                                //     latexTmp -> DrawLatexNDC(0.15, 0.7, Form("tb_temple = %.2f", parsTR[0]));
                                //     latexTmp -> DrawLatexNDC(0.15, 0.5, Form("adc_hist = %.2f", double(maxADC)));
                                //     latexTmp -> DrawLatexNDC(0.15, 0.4, Form("adc_temple = %.2f", parsTR[1]));
                                //     latexTmp -> Draw("same");

                                //     cc -> Draw();
                                //     cc -> SaveAs(Form("../test_temp/event_%i_%i_%i-R.png", tmpEventNum, agetID, chanID));
                                // }
                                data -> setRightTPCHits(pad.GetYId(agetID, chanID), pad.GetXId(agetID, chanID), parsTR[1], parsTR[0]);
                            }
                            else{
                                fFitFcnTL -> SetParameters(double(maxBin)-0.5, maxADC, double(maxADC+50.));
                                fFitFcnTL -> SetParLimits(0, double(maxBin)-20.5, double(maxBin)+20.5);
                                fFitFcnTL -> SetParLimits(1, (maxADC) - 50, (maxADC) + 50);
                                hSignal->Fit(fFitFcnTL, "QNRSM");
                                double parsTL[3];
                                fFitFcnTL->GetParameters(parsTL);

                                // cc -> cd();
                                // hSignal -> Draw();
                                // fFitFcnTL -> Draw("same");
                                

                                // latexTmp -> Clear();
                                // latexTmp -> DrawLatexNDC(0.15, 0.8, Form("tb_hist = %.2f", double(maxBin)));
                                // latexTmp -> DrawLatexNDC(0.15, 0.7, Form("tb_temple = %.2f", parsTL[0]));
                                // latexTmp -> DrawLatexNDC(0.15, 0.5, Form("adc_hist = %.2f", double(maxADC)));
                                // latexTmp -> DrawLatexNDC(0.15, 0.4, Form("adc_temple = %.2f", parsTL[1]));
                                // latexTmp -> Draw("same");

                                // cc -> Draw();
                                // cc -> SaveAs(Form("../test_temp/event_%i_%i_%i-L.png", tmpEventNum, agetID, chanID));

                                data -> setLeftTPCHits(pad.GetYId(agetID, chanID), pad.GetXId(agetID, chanID), parsTL[1], parsTL[0]);
                            }
                        }
                        hSignal->Reset("ICESM");
                    }
                }
            }



            // ================================== Si production ================================== 
            int chanIdx[6][2];
            memset(chanIdx, 0, sizeof(chanIdx));

            for (int asadID = 1; asadID <= 3; asadID += 2) {
                for(int agetID = 0; agetID<2; agetID++){
                    int chanMaxID = (agetID == 0)? 13 : 51;

                    for(int chanID = 0; chanID < chanMaxID; chanID++){

                        if (decoder[asadID].IsFPN(chanID)) continue;

                        TString label = getLabel(asadID, agetID, chanID);
                        
                        int FPNChanID = decoder[asadID].FPNChanID(chanID);
                        std::fill_n(signal, nBUCK, -9999);

                        int avgADC = 0;
                        for (int buckID = 1; buckID <= 500; buckID++) {
                            int ADC = decoder[asadID].ADC[agetID][chanID][buckID];
                            int FPN = decoder[asadID].ADC[agetID][FPNChanID][buckID];
                            signal[buckID] = ADC - FPN;
                            if (buckID <= 140){
                                avgADC+= signal[buckID];
                            }
                        }

                        int ped = avgADC/140;

                        int SiDetIdx = -1;
                        int omicOrJunc = (label.Index("o") != -1)? 0 : 1;
                        if(label.Index("rd") != -1){SiDetIdx = Si_RD;}
                        if(label.Index("ru") != -1){SiDetIdx = Si_RU;}
                        if(label.Index("ld") != -1){SiDetIdx = Si_LD;}
                        if(label.Index("lu") != -1){SiDetIdx = Si_LU;}
                        if(label.Index("rb") != -1){SiDetIdx = Si_BR;}
                        if(label.Index("lb") != -1){SiDetIdx = Si_BL;}
                        
                        for(int buckID=1; buckID<=500; buckID++){
                            if(omicOrJunc==0){
                                data->setSiOmicADC(SiDetIdx, chanIdx[SiDetIdx][omicOrJunc], buckID, signal[buckID]-ped);
                            }
                            else if(omicOrJunc == 1){
                                data->setSiJuncADC(SiDetIdx, chanIdx[SiDetIdx][omicOrJunc], buckID, signal[buckID]-ped);
                            }
                        }
                        chanIdx[SiDetIdx][omicOrJunc]++;
                    }
                }
            }
            for(int si=0; si<6; si++){
                
                // omic-side 
                // find a peak channel
                int maxOmicADCIdx = 0;
                int tmpOmicADC = 0;
                for(int idx=0; idx<4; idx++){
                    int avgADC = 0;
                    for(int tb=185; tb<201; tb++){ // around peak
                        avgADC += data->getSiOmicADC(si, idx, tb);
                    }
                    avgADC /= 16;

                    if(tmpOmicADC < abs(avgADC)){
                        tmpOmicADC = avgADC;
                        maxOmicADCIdx = idx;
                    }
                }

                // Junction-side
                // find a peak channel
                int maxJuncADCIdx = 0;
                int tmpJuncADC = 4096;
                for(int idx=0; idx<16; idx++){
                    int avgADC = 0;
                    for(int tb=185; tb<201; tb++){ // around peak
                        avgADC += data->getSiJuncADC(si, idx, tb);
                    }
                    avgADC /= 16;

                    if(tmpJuncADC > avgADC){
                        tmpJuncADC = avgADC;
                        maxJuncADCIdx = idx;
                    }
                }

                int passJuncIdx[2];
                if(maxJuncADCIdx % 2 == 0){
                    passJuncIdx[0] = maxJuncADCIdx;
                    passJuncIdx[1] = maxJuncADCIdx+1;
                }
                else if(maxJuncADCIdx % 2 == 1){
                    passJuncIdx[0] = maxJuncADCIdx;
                    passJuncIdx[1] = maxJuncADCIdx-1;
                }

                int pedAvg[2];
                int omicPedADC[512];
                int juncPedADC[512];
                int omicSigAvg[4]; // average ADC chan by chan
                int juncSigAvg[4]; // average ADC chan by chan
                memset(pedAvg, 0, sizeof(pedAvg));
                memset(omicPedADC, 0, sizeof(omicPedADC));
                memset(juncPedADC, 0, sizeof(juncPedADC));
                memset(omicSigAvg, 0, sizeof(omicSigAvg));
                memset(juncSigAvg, 0, sizeof(juncSigAvg));
                
                // pedestal subtraction
                for(int tb=1; tb<=500; tb++){
                    for(int omicIdx=0; omicIdx<4; omicIdx++){
                        if(tb<=140){omicSigAvg[omicIdx] += data->getSiOmicADC(si, omicIdx, tb);}
                        if(omicIdx == maxOmicADCIdx){continue;}
                        omicPedADC[tb] += data->getSiOmicADC(si, omicIdx, tb);
                    }
                    omicPedADC[tb] /= 3;

                    for(int juncIdx=0; juncIdx<16; juncIdx++){
                        if(tb<=140){juncSigAvg[juncIdx] += data->getSiJuncADC(si, juncIdx, tb);}
                        if(juncIdx == passJuncIdx[0] || juncIdx == passJuncIdx[1]){continue;}
                        juncPedADC[tb] += data->getSiJuncADC(si, juncIdx, tb);
                    }
                    juncPedADC[tb] /= 14;
                    
                    if(tb<=140){
                        pedAvg[0] += omicPedADC[tb];
                        pedAvg[1] += juncPedADC[tb];
                    }
                }
                pedAvg[0] /= 140;
                pedAvg[1] /= 140;

                for(int idx=0; idx<4; idx++){      
                    int maxADC = 0;
                    omicSigAvg[idx] /= 140;
                    for(int tb=1; tb<=500; tb++){
                        int omicADCCorr = data->getSiOmicADC(si, idx, tb) - (omicPedADC[tb]+ (omicSigAvg[idx] - pedAvg[0]) );
                        data->setSiOmicADC(si, idx, tb, omicADCCorr);
                        
                        if(tb < 140 || tb > 220){continue;}
                        if(maxADC < omicADCCorr){maxADC = omicADCCorr;}
                    }
                    if(maxADC < thresholdOmic[si]){maxADC = 0;}
                    data->setSiOmicHits(si, idx, maxADC);
                }

                for(int idx=0; idx<16; idx++){
                    int maxADC = 0;
                    juncSigAvg[idx] /= 140;
                    for(int tb=1; tb<=500; tb++){
                        int juncADCCorr = data->getSiJuncADC(si, idx, tb) - (juncPedADC[tb]+ (juncSigAvg[idx] - pedAvg[1]) );
                        data->setSiJuncADC(si, idx, tb, juncADCCorr);  

                        if(tb < 140 || tb > 220){continue;}
                        if(maxADC > juncADCCorr){maxADC = juncADCCorr;}           
                    }
                    if(abs(maxADC) < thresholdJunc[si]){maxADC = 0;}
                    data -> setSiJuncHits(si, idx, abs(maxADC));
                }
            }

            treeOut->Fill();

        }
    }
    fileOut -> cd();
    treeOut -> Write();
    fileOut -> Close();

    cout << " done !!! " << endl;
}

std::string getLabel(int asadID, int agetID, int chanID) {
    std::stringstream label;
    std::string detector;
    std::string rightOrLeft = (asadID == 0 || asadID == 1) ? "r" : "l";
    std::string downOrUp;
    std::string ohmicOrJunction;
    if (chanID == 11 || chanID == 22 || chanID == 45 || chanID == 56) {
        label << "FPN";
        return label.str();
    }
    if (asadID == 0 || asadID == 2) {
        detector = "t";
    } else if (agetID == 2) {
        label << "none";
        return label.str();
    } else if (asadID == 1 && agetID == 3 && chanID == 49) {
        detector = "c";
        downOrUp = "d";
    } else if (asadID == 1 && agetID == 3 && chanID == 15) {
        detector = "c";
        downOrUp = "u";
    } else if (asadID == 3 && agetID == 3 && chanID == 66) {
        detector = "c";
        downOrUp = "d";
    } else if (asadID == 3 && agetID == 3 && chanID == 32) {
        detector = "c";
        downOrUp = "u";
    } else if (agetID == 0 && chanID <= 12) {
        detector = "s";
        ohmicOrJunction = "o";
        downOrUp = (chanID <= 3)    ? "d"
                   : (chanID <= 7)  ? "u"
                   : (chanID <= 12) ? "b"
                                    : "";
    } else if (agetID == 1 && chanID <= 50) {
        detector = "s";
        ohmicOrJunction = "j";
        downOrUp = (chanID <= 16)   ? "d"
                   : (chanID <= 33) ? "u"
                   : (chanID <= 50) ? "b"
                                    : "";
    } else {
        label << "none";
        return label.str();
    }

    label << detector << rightOrLeft << downOrUp << ohmicOrJunction;
    return label.str();
}

int mode(int array[], int size) {
    std::sort(array, array + size);
    int counter = 1;
    int max = 0;
    int mode = array[0];
    for (int idx = 0; idx < size - 1; idx++) {
        if (array[idx] == array[idx + 1]) {
            counter++;
            if (counter > max) {
                max = counter;
                mode = array[idx];
            }
        } else
            counter = 1;  // reset counter.
    }
    return mode;
}

void pulsegenerate()
{
    GETPad pad;

    // TFile* fileOut = new TFile("../data/pulsegeneration.root", "recreate");
    // TTree* treeOut = new TTree("event","event");

    TClonesArray* outClonesArray = new TClonesArray("dataStructure");
    // treeOut -> Branch("dataStructure", &outClonesArray);
    dataStructure* data;
    

    readSigFcn();

    TFile* fileOut = new TFile("../data/pulsegeneration.root", "recreate");
    TTree* treeOut = new TTree("event","event");

    double ans_time;
    double ans_ADC;
    double binADC[500];

    treeOut-> Branch("ans_time", &ans_time, "ans_time/D");
    treeOut-> Branch("ans_ADC", &ans_ADC, "ans_ADC/D");
    treeOut-> Branch("binADC", &binADC, "binADC[500]/D");

    int signal[512];
    bool isBroken;

    TCanvas* cc = new TCanvas();
    TLatex* latexTmp = new TLatex();

    auto hSignal = new TH1D("hSignal", "", 500, 0.5, 500.5);
    auto hPedestal = new TH1D("hPedestal", "", 1000, -499.5, 500.5);
    TF1* fFitFcnTR = new TF1("fFitFcnTR", fitFcnTR, 1, 500, 3);
    TF1* fFitFcnTL = new TF1("fFitFcnTL", fitFcnTL, 1, 500, 3);
    
    for (int runID = 2; runID < 3; runID++) {
        GETDecoder decoder[nAsAd];

        for (int asadID = 0; asadID < nAsAd; asadID++){
            // gSystem->Exec(Form("ls /home/shlee/workspace/forHIMAC/data/carbon_200MeV/run%02d/CoBo0_AsAd%d* > ../data/fileList%02d.txt", runList[runID], asadID, asadID));
            decoder[asadID].OpenFromList(Form("../data/fileList%02d.txt", asadID));
        }

        int tmpEventNum = 0;
        while (decoder[1].Run()) {
            if(tmpEventNum%100==0){cout << "processing... (runID: "<< runList[runID] << " | EventID: " << tmpEventNum << ") " << endl;}
            tmpEventNum++;
            // if(tmpEventNum==101){break;}

            outClonesArray -> Clear("C");
            data = (dataStructure*)outClonesArray -> ConstructedAt(0);
            ans_time = -999;
            ans_ADC = -999;
            memset(binADC, -999, sizeof(binADC));

            // Event info
            for (int asadID = 0; asadID < 1; asadID++) {
                if (asadID != 1) decoder[asadID].Run();
                // if (asadID == 1) {
                //     data -> setRunID(runList[runID]);
                //     data -> setEventID(int(decoder[asadID].EventID));
                //     data -> setEventTime(double(decoder[asadID].EventTime));
                //     data -> setEventDiffTime(double(decoder[asadID].DiffTime));
                //     isBroken = false;
                // }
            }
            
            // ================================== AT-TPC production ================================== 
            for (int asadID = 0; asadID <= 1; asadID ++) {
                for (int agetID = 0; agetID < 4; agetID++) {
                    for (int chanID = 0; chanID < 68; chanID++) {
                        if (decoder[asadID].IsFPN(chanID) || decoder[asadID].IsDead(chanID)) continue;
                        int FPNChanID = decoder[asadID].FPNChanID(chanID);
                        for (int buckID = 1; buckID <= 500; buckID++) {
                            int ADC = decoder[asadID].ADC[agetID][chanID][buckID]-decoder[asadID].ADC[agetID][FPNChanID][buckID];
                            

                            if (buckID <= 150) {
                                hPedestal->Fill(ADC);
                            } else {
                                hSignal->SetBinContent(buckID, ADC);
                                
                            }
                        }

                        int pedADC = hPedestal->GetBinCenter(hPedestal->GetMaximumBin());
                        int maxBin = hSignal->GetMaximumBin();
                        int maxADC = hSignal->GetMaximum() - pedADC;
                        // cout<<maxBin<< " : "<<maxADC<<endl;
                        
                        // cout<<maxADC<<endl;
                        hSignal->Reset("ICESM");
                        hPedestal->Reset("ICESM");

                        if (150 < maxBin && maxBin <= 400 && 2000 < maxADC && maxADC < 3500) {
                            for (int buckID = 1; buckID <= 500; buckID++) {
                                int ADC = decoder[asadID].ADC[agetID][chanID][buckID]-decoder[asadID].ADC[agetID][FPNChanID][buckID];
                                hSignal->SetBinContent(buckID, ADC - pedADC);
                                binADC[buckID] = ADC - pedADC;
                            }

                            fFitFcnTR -> SetParameters(double(maxBin)-0.5, maxADC, double(maxADC+50.));
                            fFitFcnTR -> SetParLimits(0, double(maxBin)-20.5, double(maxBin)+20.5);
                            fFitFcnTR -> SetParLimits(1, (maxADC) - 50, (maxADC) + 50);
                            hSignal->Fit(fFitFcnTR, "QNRSM");
                            double parsTR[3];
                            fFitFcnTR->GetParameters(parsTR);
                            ans_time = parsTR[0];
                            ans_ADC = parsTR[1];
  
                        }
                        hSignal->Reset("ICESM");
                    }
                }
            }
            if(ans_time > 0 && ans_ADC > 0){
                treeOut->Fill();
            }
        }
    }
    fileOut -> cd();
    treeOut -> Write();
    fileOut -> Close();

    cout << " done !!! " << endl;
}
