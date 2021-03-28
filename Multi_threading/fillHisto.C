const Int_t nNumbers = 20000000;

const Int_t nWorkers = 4;
Int_t mp001_fillHistos()
{
   
   auto workItem = [] (workerID){
      
      TRandom3 workerRndm(workID); 
      TFile f(Form("myFile_%u.root", workerID), "RECREATE");
      TH1F h(Form("myHisto_%u", workerID), "The Histogram", nbin, low, up);
      for (UInt_t i = 0; i < nNumbers; ++i) {
         h.Fill(workerRndm.Gaus());
      }
      h.Write();
      return 0;
   };
  
   ROOT::TProcessExecutor workers(nWorkers);
 
   workers.Map(workItem, ROOT::TSeqI(nWorkers));
   return 0;
}
