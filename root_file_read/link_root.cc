#include "/home/chanchoen/workspace/nptool/NPLib/include/TTrackInfo.h"

using namespace std;

int link_root(){
    TFile *rootFile = TFile::Open("/home/chanchoen/workspace/nptool/Outputs/Simulation/hoyle_Reaction.root");
    
    TTree *tree = (TTree*)rootFile->Get("SimulatedTree");
    TTrackInfo* trackInfo = new TTrackInfo();

    tree->SetBranchAddress("TrackInfo", &trackInfo);
    int nEntries = tree->GetEntries();

    TFile* fileOut = new TFile("./hoyle_reaction.root","recreate");
    TTree* treeOut = new TTree("data", "data");

    Int_t particleSize;
    vector<Double_t> KineticE;
    vector<Double_t> Momentum_X, Momentum_Y, Momentum_Z;
    vector<Double_t> Position_X, Position_Y, Position_Z;
    vector<Int_t> IonZ, IonA;

    treeOut->Branch("particleSize", &particleSize);
    treeOut->Branch("KineticE", &KineticE);
    treeOut->Branch("Momentum_X", &Momentum_X);
    treeOut->Branch("Momentum_Y", &Momentum_Y);
    treeOut->Branch("Momentum_Z", &Momentum_Z);
    treeOut->Branch("Position_X", &Position_X);
    treeOut->Branch("Position_Y", &Position_Y);
    treeOut->Branch("Position_Z", &Position_Z);
    treeOut->Branch("IonZ", &IonZ);
    treeOut->Branch("IonA", &IonA);
    
    for(int events = 0; events < nEntries; events++){
        if(events%100==0){cout << "event: " << events<< " / "<< nEntries << endl;}
        tree->GetEntry(events);

        KineticE.clear();
        Momentum_X.clear();
        Momentum_Y.clear();
        Momentum_Z.clear();
        Position_X.clear();
        Position_Y.clear();
        Position_Z.clear();
        IonZ.clear();
        IonA.clear();

        int Size = trackInfo->GetParticleMultiplicity();
        particleSize = trackInfo->GetParticleMultiplicity();
        // cout<<"multiplicity : "<<particleSize<<endl;
        
        for(int i=0; i<Size; i++){
            int ionZ = trackInfo->GetZ(i);
            int ionA = trackInfo->GetA(i);
            double kE = trackInfo->GetKineticEnergy(i);
            double momX = trackInfo->GetMomentumX(i);
            double momY = trackInfo->GetMomentumY(i);
            double momZ = trackInfo->GetMomentumZ(i);
            double posX = trackInfo->GetPositionX(i);
            double posY = trackInfo->GetPositionY(i);
            double posZ = trackInfo->GetPositionZ(i);

            if(ionZ == 0){ 
                particleSize = particleSize -1; 
                continue; 
            }

            KineticE.push_back(kE);
            Momentum_X.push_back(momX);
            Momentum_Y.push_back(momY);
            Momentum_Z.push_back(momZ);

            Position_X.push_back(posX);
            Position_Y.push_back(posY);
            Position_Z.push_back(posZ);
        
            IonZ.push_back(ionZ);
            IonA.push_back(ionA);
        }
        treeOut->Fill();

    }
    fileOut->cd();
    treeOut->Write();
    fileOut->Close();
    
    
    return 0;
}