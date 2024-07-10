
int link_root_read(){
    TFile *rootFile = TFile::Open("./hoyle_reaction.root");
    TTree *tree = (TTree*)rootFile->Get("data");

    Int_t particleSize;
    vector<Double_t> *vMom_X = new vector<Double_t>;
    vector<Double_t> *vMom_Y = new vector<Double_t>;
    vector<Double_t> *vMom_Z = new vector<Double_t>;
    vector<Double_t> *vPosX = new vector<Double_t>;
    vector<Double_t> *vPosY = new vector<Double_t>;
    vector<Double_t> *vPosZ = new vector<Double_t>;

    tree->SetBranchAddress("particleSize", &particleSize);
    tree->SetBranchAddress("Momentum_X", &vMom_X);
    tree->SetBranchAddress("Momentum_Y", &vMom_Y);
    tree->SetBranchAddress("Momentum_Z", &vMom_Z);
    tree->SetBranchAddress("Position_X", &vPosX);
    tree->SetBranchAddress("Position_Y", &vPosY);
    tree->SetBranchAddress("Position_Z", &vPosZ);
    int nEntries = tree->GetEntries();
    
    for(int event=218; event<220; event++){
        tree->GetEntry(event);
        cout<<particleSize<<endl;
        for(int i=0; i<particleSize; i++){
            double momX = vMom_X->at(i);
            double momY = vMom_Y->at(i);
            double momZ = vMom_Z->at(i);

            double posX = vPosX->at(i);
            double posY = vPosY->at(i);
            double posZ = vPosZ->at(i) + 125.;
        
            // cout<<"momX : "<<momY<<" momY : "<<momX<<" momZ : "<<momZ<<endl;
            cout<<"posX : "<<posY<<" posY : "<<posX<<" posZ : "<<posZ<<endl;
        }
    }
    
    rootFile->Close();
    
    return 0;
}