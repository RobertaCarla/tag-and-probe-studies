void rebin() {

    TFile rebin("rebin.root", "RECREATE");


   //create a fix bin histogram
   TH1F *h = new TH1F("h","test rebin",100,-3,3);
   Int_t nentries = 1000;
   h->FillRandom("gaus",nentries);
   Double_t xbins[1001];
   Int_t k=0;
   TAxis *axis = h->GetXaxis();
   for (Int_t i=1;i<=100;i++) {
      Int_t y = (Int_t)h->GetBinContent(i);
      if (y <=0) continue;
      Double_t dx = axis->GetBinWidth(i)/y;
      Double_t xmin = axis->GetBinLowEdge(i);
      for (Int_t j=0;j<y;j++) {
         xbins[k] = xmin +j*dx;     //divide ogni bin in n bin dove n è inversamente proporzionale al numero di entries del bin (y)
         k++;                      //"aggiorna" il numero dell'elemento di xbins[] nel ciclo for con j (invece di mettere j così xbins[j] come avrei fatto io)
      }
   }
   xbins[k] = axis->GetXmax();    //credo stia dando all'elemento k-esimo il valore massimo dell'asse x. Ci sta nel senso che è quello che ho fatto anch'io quando cercavo di costruire un array con gli estremi dei bin
   //create a variable bin-width histogram out of fix bin histogram
   //new rebinned histogram should have about 10 entries per bin
   TH1F *hnew = new TH1F("hnew","rebinned",k,xbins);   //ora costruisce un istogramma con 10*il numero di dati che aveva prima, utilizzando il nuovo binning
   hnew->FillRandom("gaus",10*nentries);               //xbins[k+1]   k è il numero di bin
 
   //rebin hnew keeping only 50% of the bins
   Double_t xbins2[501];
   Int_t kk=0;
   for (Int_t j=0;j<k;j+=2) {    //j salta di due in due mentre kk avanza di uno alla volta: in questo modo quando arrivo a j=1000, k=500 
      xbins2[kk] = xbins[j];
      kk++;
   }
   xbins2[kk] = xbins[k];       //l'elemento kk-esimo (il 501) di xbins2 è uguale all'elemento k-esimo (il 1001) di xbins
   TH1F *hnew2 = (TH1F*)hnew->Rebin(kk,"hnew2",xbins2);   //ORA FACCIO IL REBINNING VERO E PROPRIO: infatti, come si può vedere, parto da hnew (la mia funzione di prima)
   // e la rebinno mettendo il numero di bin, il nuovo nome(altrimenti non mi fa una copia rebinnata ma mi rebinna l'istogramma iniziale) e l'array con gli estremi dei bin
   //NB: (TH1F*) per dire il tipo
 
   //draw the 3 histograms
   TCanvas *c1 = new TCanvas("c1","c1",800,1000);
   c1->Divide(1,3);
   c1->cd(1);
   h->Draw();
   c1->cd(2);
   hnew->Draw();
   c1->cd(3);
   hnew2->Draw();

    rebin.cd();
    c1->Write();

}