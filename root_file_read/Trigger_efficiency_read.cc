int Trigger_efficiency_read(){
    auto c1 = new TCanvas();
    auto hNumCsI = new TH1D("hNumCsI"," ", 10, 0, 10);
    auto htrigger1 = new TH1D("htrigger1"," ", 2, 0, 2);
    auto htrigger2 = new TH1D("htrigger2"," ", 2, 0, 2);
    auto htrigger3 = new TH1D("htrigger3"," ", 2, 0, 2); htrigger3->SetStats(0);
    auto htrigger4 = new TH1D("htrigger4"," ", 2, 0, 2);
    



    TFile *rootFile = TFile::Open("./hoyle_reaction_CsI.root");
    TTree *tree = (TTree*)rootFile->Get("data");

    Int_t sizeCsI;
    vector<Int_t> *fCsI_CrystalNumber = new vector<Int_t>;

    tree->SetBranchAddress("sizeCsI", &sizeCsI);
    tree->SetBranchAddress("fCsI_CrystalNumber", &fCsI_CrystalNumber);
    int nEntries = tree->GetEntries();

    int nEvent = 0;
    
    for(int event=0; event<nEntries; event++){
        tree->GetEntry(event);
        int channel1 = 0;
        int channel2 = 0;
        int channel3 = 0;
        int channel4 = 0;
        if(sizeCsI > 0){
            nEvent ++;
        }
        for(int i=0; i<sizeCsI; i++){
            int nCsI = fCsI_CrystalNumber->at(i);

            if(nCsI == 0 || nCsI == 4 || nCsI == 8 || nCsI == 12){
                channel1 = 1;
            }
            if(nCsI == 1 || nCsI == 5 || nCsI == 9 || nCsI == 13){
                channel2 = 1;
            }
            if(nCsI == 2 || nCsI == 6 || nCsI == 10 || nCsI == 14){
                channel3 = 1;
            }
            if(nCsI == 3 || nCsI == 7 || nCsI == 11 || nCsI == 15){
                channel4 = 1;
            }
        }
        // cout<<"channel1 : " <<channel1<<" channel2 : " <<channel2<<" channel3 : " <<channel3<<" channel4 : " <<channel4<<endl;
        if( channel1 == 1 && channel4 == 1 ){
            htrigger1->Fill(1);
        }
        if( channel1 == 1 && channel3 == 1 ){
            htrigger2->Fill(1);
        }
        if( channel2 == 1 && channel3 == 1 ){
            htrigger3->Fill(1);
        }
        if( channel2 == 1 && channel4 == 1 ){
            htrigger4->Fill(1);
        }
        hNumCsI->Fill(sizeCsI);
    }
    cout<<nEvent<<endl;
    rootFile->Close();
    auto latex = new TLatex();
    latex->SetTextSize(0.05);
    c1->cd();
    // hNumCsI->Draw();
    htrigger1->SetLineColor(kRed);
    htrigger2->SetLineColor(kBlue);
    htrigger3->SetLineColor(kGreen);
    htrigger4->SetLineColor(kMagenta);

    
    htrigger3->Draw();
    htrigger1->Draw("same");
    htrigger2->Draw("same");
    htrigger4->Draw("same");
    latex->DrawLatexNDC(0.12, 0.6, Form("#color[2]{logic A : %.1f %%}", (htrigger1->GetEntries() / nEntries) * 100));
    latex->DrawLatexNDC(0.12, 0.5, Form("#color[4]{logic B : %.1f %%}", (htrigger2->GetEntries() / nEntries) * 100));
    latex->DrawLatexNDC(0.12, 0.4, Form("#color[3]{logic C : %.1f %%}", (htrigger3->GetEntries() / nEntries) * 100));
    latex->DrawLatexNDC(0.12, 0.3, Form("#color[6]{logic D : %.1f %%}", (htrigger4->GetEntries() / nEntries) * 100));

    latex->SetTextSize(0.08);
    latex->DrawLatexNDC(0.12, 0.75, "Logic efficiency");
    c1->Draw();
    
    return 0;
}