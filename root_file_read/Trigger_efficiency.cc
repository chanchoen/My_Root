#include "/home/chanchoen/workspace/nptool_origin/NPLib/include/TActarData.h"

using namespace std;

int Trigger_efficiency(){
    TFile *rootFile = TFile::Open("/home/chanchoen/workspace/nptool_origin/Outputs/Simulation/hoyle_Reaction.root");
    
    TTree *tree = (TTree*)rootFile->Get("SimulatedTree");
    TActarData* actardata = new TActarData();
    
    tree->SetBranchAddress("Actar", &actardata);
    int nEntries = tree->GetEntries();

    TFile* fileOut = new TFile("./hoyle_reaction_CsI.root","recreate");
    TTree* treeOut = new TTree("data", "data");

    Int_t sizeCsI;
    vector<Int_t> fCsI_CrystalNumber;
    vector<Double_t> fCsI_Energy;

    treeOut->Branch("sizeCsI", &sizeCsI);
    treeOut->Branch("fCsI_CrystalNumber", &fCsI_CrystalNumber);
    treeOut->Branch("fCsI_Energy", &fCsI_Energy);
    
    for(int events = 0; events < nEntries; events++){
        if(events%100==0){cout << "event: " << events<< " / "<< nEntries << endl;}
        tree->GetEntry(events);
        fCsI_CrystalNumber.clear();
        fCsI_Energy.clear();
        sizeCsI = actardata->GetCsIMult();
        for(int i=0; i<sizeCsI; i++){
            int numCsI = actardata->Get_CsICrystalNumber(i);
            double eCsI = actardata->Get_CsIEnergy(i);

            fCsI_CrystalNumber.push_back(numCsI);
            fCsI_Energy.push_back(eCsI);
        }
        treeOut->Fill();

    }
    fileOut->cd();
    treeOut->Write();
    fileOut->Close();
    
    
    return 0;
}