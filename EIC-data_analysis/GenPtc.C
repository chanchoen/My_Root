#include "TROOT.h"
#include "TH2F.h"

#include <vector>
#include <utility>
#include <map>
#include <tuple>

void GenPtc() {
    //for show the generated position, and then the particle enter at absorber or fiber
    //hit's histo has two peak, if cut the other peak's data, how show the generated position
    TH2F *Genxy = new TH2F("Genxy", "Gen vx : vy",100, 446, 452, 100, -16, -5);

    //for call the sipmData define same Struct in the file
    typedef std::pair<float,float> hitRange;
    typedef std::pair<int,int> hitXY;
    typedef std::map<hitRange, int> DRsimTimeStruct;
    typedef std::map<hitRange, int> DRsimWavlenSpectrum;
    typedef std::tuple<float,float,float> threeVector;

    struct DRsimSiPMData { 
    DRsimSiPMData() {};   
    virtual ~DRsimSiPMData() {};

    int count;
    int SiPMnum;
    int x;
    int y;
    threeVector pos;
    DRsimTimeStruct timeStruct;
    DRsimWavlenSpectrum wavlenSpectrum;
    };

    TFile* output = new TFile("./tungsten/Tungsten_ele_200GeV.root", "READ"); //call target file
    TTree* tree = (TTree*)output ->Get("DRsim");

    Float_t vx, vy, vz;//Generated particle position

    tree-> SetBranchAddress("GenPtcs.vx", &vx);
    tree-> SetBranchAddress("GenPtcs.vy", &vy);
    tree-> SetBranchAddress("GenPtcs.vz", &vz);

    Int_t num = (Int_t)tree->GetEntries();

    std::vector<DRsimSiPMData> SIPMs; //Sipm data

    tree-> SetBranchAddress("towers.SiPMs", &SIPMs);
    Float_t t_hits = 0;
    
    for(int i=0; i< SIPMs.size(); i++){
        DRsimTimeStruct timeItr = SIPMs[i].timeStruct;

        for(auto TmpItr = timeItr.begin(); TmpItr != timeItr.end(); ++TmpItr){  
            auto timeData = *TmpItr;
            t_hits += timeData.second;
            if(t_hits < 55000){ //cut hit point
                for(int j=0; j<num; j++){
                    tree->GetEntry(j);
                    Genxy->Fill(vx, vy);
                }
            } 
        }        
    }
    cout<< num<<endl;
    cout<< t_hits<< endl;

    TCanvas *c1 = new TCanvas();

    Genxy->Draw("surf1");
    c1->Draw(); 
}
