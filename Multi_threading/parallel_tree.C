{
  
  int nthreads = 10;
  ROOT::EnableImplicitMT(nthreads);
  
  TFile *file = TFile::Open("filename.txt");

  TTree *tree = nullptr;
  file->GetObject<TTree>("h1", tree);
 
  
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
     tree->GetEntry(i);
  }
  
  tree->SetImplicitMT(false);
 
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
     tree->GetEntry(i); 
  }

  tree->SetImplicitMT(true);
  
  
  ROOT::DisableImplicitMT();

  
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
     tree->GetEntry(i);
  }
  return 0;
}
