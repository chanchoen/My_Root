{
  float alpha_energy = 3731.879;//Mev
  float alpha_Mass = 3728.43;//Mev
  float alpha_momentum = 160.408;
  float carbon_energy = 11195.636;//Mev
  float carbon_Mass = 11187.986; //Mev
  
  float resolution = 0.01;
  int nEvent = 10000;
  double a = TMath::Pi()/3;
 
  TH1D *Etotal = new TH1D("Etotal", "",1000, 10900, 11500);
  TH1D *henergy1 = new TH1D("henergy", "", 1000, 3550, 3900);
  TH1D *henergy2 = new TH1D("henergy", "", 1000, 3550, 3900);
  TH1D *henergy3 = new TH1D("henergy", "", 1000, 3550, 3900);
    
  for(int i=0; i<nEvent; i++){
    TRandom3 *rand = new TRandom3(i);
    TLorentzVector v1;
    TLorentzVector v2;
    TLorentzVector v3;
    TLorentzVector v0;
    
    v1.SetPxPyPzE(0,0,alpha_momentum, alpha_energy*(1+rand->Gaus(0,resolution)));
    v2.SetPxPyPzE(0,-cos(a)*alpha_momentum,-sin(a)*alpha_momentum, alpha_energy*(1+rand->Gaus(0,resolution)));
    v3.SetPxPyPzE(0, cos(a)*alpha_momentum,-sin(a)*alpha_momentum,alpha_energy*(1+rand->Gaus(0,resolution)));

    float px1 = v1.Px();
    float py1 = v1.Py();
    float pz1 = v1.Pz();
    float e1 = v1.E();
    float Mass1 = v1.Mag();

    float px2 = v2.Px();
    float py2 = v2.Py();
    float pz2 = v2.Pz();
    float e2 = v2.E();
    float Mass2 = v2.Mag();

    float px3 = v3.Px();
    float py3 = v3.Py();
    float pz3 = v3.Pz();
    float e3 = v3.E();
    float Mass3 = v3.Mag();

    v0.SetPxPyPzE(0,0,0, e1+e2+e3);
    float Massc =v0.Mag();

    henergy1->Fill(Mass1);
    henergy2->Fill(Mass2);
    henergy3->Fill(Mass3);
    Etotal->Fill(Massc);
  }
  
  TF1 *Fit1 = new TF1("Fit1", "gaus", 0, 20);
  TF1 *Fit2 = new TF1("Fit2", "gaus", 0, 20);
  TF1 *Fit3 = new TF1("Fit3", "gaus", 0, 20);
  TF1 *Fit0 = new TF1("Fit0", "gaus", 0, 20);
     
  TCanvas *c1 = new TCanvas("c1", "", 500, 1000);
  c1->Divide(1,2);

  c1->cd(1);
  henergy1->Draw();
  henergy2->Draw("same");
  henergy3->Draw("same");
  henergy1->Fit("Fit1");
  henergy2->Fit("Fit2");
  henergy3->Fit("Fit3");
  c1->Draw();
  float sigma1 = Fit1->GetParameter(2);
  float sigma2 = Fit2->GetParameter(2);
  float sigma3 = Fit3->GetParameter(2);
  float sigmatot = (sigma1+sigma2+sigma3)/3;
  cout<< sigmatot <<endl;

  c1->cd(2);
  Etotal->Draw();
  Etotal->Fit("Fit0");
  c1->Draw();

  float sigma0 =Fit0->GetParameter(2);

  cout<< sigma0<<endl;
  
}
