{
  TH2Poly *h2p = new TH2Poly();
  
  for(double i=0; i<5;i++){
    for(double j=0; j<5; j++){
      double x1[] ={1+1.8*j, 1+1.8*j, 1.8+1.8*j, 2.7+1.8*j, 2.7+1.8*j, 1.8+1.8*j};
      double y1[] = {1+3.5*i, 2+3.5*i, 2.5+3.5*i, 2+3.5*i, 1+3.5*i, 0.5+3.5*i};
      h2p->AddBin(6,x1,y1);
    }
  }
  for(double i=0; i<5; i++){
    for(double j=0; j<5; j++){
      double x2[] ={1.8+1.8*j, 1.8+1.8*j, 2.7+1.8*j, 3.6+1.8*j, 3.6+1.8*j, 2.7+1.8*j};
      double y2[] ={2.7+3.5*i, 3.7+3.5*i, 4.2+3.5*i, 3.7+3.5*i, 2.7+3.5*i, 2.2+3.5*i};
      h2p->AddBin(6,x2,y2);
    }
  }

  h2p->Fill(1,5,100);
  h2p->Draw("colz");
  c1->Draw();
}
