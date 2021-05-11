{
  fstream file;
  file.open("<the name>.txt", ios::in);

  double x, y;

  TFile *output = new TFile("<new created file>.root", "RECREATE");
  TTree *tree = new TTree("tree", "tree");

  tree->Branch("x", &x, "x/D"); //get the branch and then add the data
  tree->Branch("y", &y, "y/D");

  while(1){
    file>>x>>y;
    
    if(file.eof()) break;
    cout << x<< y<<endl;

    tree->Fill();
  }
  output->Write();
  output->Close();


  file.close();
}
