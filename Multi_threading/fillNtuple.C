const Int_t nNumbers = 20000000;

const Int_t nWorkers = 4;

const auto workSize = nNumbers / nWorkers;

void fillRandom (TNtuple & ntuple, TRandom3 & rndm, Int_t n)
{
   for (auto i : ROOT::TSeqI(n)) ntuple.Fill(rndm.Gaus());
}
Int_t mp101_fillNtuples()
{
   
   gROOT->SetBatch();

   TRandom3 rndm(1);
   TFile ofile("example.root", "RECREATE");
   TNtuple randomNumbers(" ", " ", "r");
   fillRandom(randomNumbers, rndm, nNumbers);
   randomNumbers.Write();
   ofile.Close();
  
   auto workItem = [](UInt_t workerID) {
      TRandom3 workerRndm(workerID);
     
      TFile ofile(Form("mp101_multiCore_%u.root", workerID), "RECREATE");
      TNtuple workerRandomNumbers("multiCore", "Random Numbers", "r");
      fillRandom(workerRandomNumbers, workerRndm, workSize);
      workerRandomNumbers.Write();
     
      return 0;
   };
  
   ROOT::TProcessExecutor workers(nWorkers);

   workers.Map(workItem, ROOT::TSeqI(nWorkers));
  
   return 0;
}
