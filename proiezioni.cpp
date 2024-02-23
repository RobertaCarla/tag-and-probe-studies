using namespace std;

void proiezioni(){
    TFile plotTesi("plotTesi.root","UPDATE");

    TFile finput("AnalysisResults.root");
    TH2F *piRec=(TH2F*)finput.Get("efficiency-q-a/piRec");   //istogramma in 2D
    TH1F *piRecAll=(TH1F*)piRec->ProjectionY("piRecAll",1,1);   //ne prendo la proiezione lungo y, la proiezione del bin "piRecAll"--> è un istogramma 1D
    TH1F *piRecITSonly=(TH1F*)piRec->ProjectionY("piRecITSonly",2,2); 
    TH1F *piRecITSTPC=(TH1F*)piRec->ProjectionY("piRecITS+TPC",4,4); 

    TH1F *All= new TH1F(*piRecAll);
    TH1F *ITSonly= new TH1F(*piRecITSonly);
    TH1F *ITS_TPC= new TH1F(*piRecITSTPC);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    All->SetLineColor(4);      //blue
    ITSonly->SetLineColor(2);  //red
    ITS_TPC->SetLineColor(3);  //green

     TCanvas *Cprojection = new TCanvas("Cprojection","projections",200,50,600,400); 
     gStyle->SetOptStat(0);  
     Cprojection->cd(0);
     All->Draw();
     ITSonly->Draw("lSAME");
     ITS_TPC->Draw("lSAME");

     TLegend *leg= new TLegend(0.73,0.85,0.87,0.65);
     leg->SetHeader("Projections", "C");   //C -> viene messo al centro
     leg->AddEntry(All, "All","l");
     leg->AddEntry(ITSonly, "ITS only","l");
     leg->AddEntry(ITS_TPC, "ITS+TPC","l");
     leg->SetTextSize(0.038);
     leg->Draw();

   Cprojection->SaveAs("projections.pdf");


     plotTesi.cd();
     Cprojection->Write();
   /*    gDirectory->Delete("Cprojection;1");*/
     plotTesi.Close();
    
}