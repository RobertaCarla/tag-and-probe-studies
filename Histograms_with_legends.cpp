void Histograms_with_legends()
{
   TFile rebin("rebin.root", "RECREATE");
    auto canvas = new TCanvas("canvas","window title",200,10,600,400); //first entry has to match variable name
    canvas->SetGrid();
    gStyle->SetOptStat(0);
    
    //definisco i punti del primo set dati
    float Vb01[]={0.55,0.57,0.58,0.59};
    float sVb01[]={0.01,0.01,0.01,0.01};

    float Ib01[]={0.096,0.201,0.29,0.392};
    float sIb01[]={0.01,0.01,0.01,0.01};

    // definisco un secondo set di dati
    float Vb02[]={0.60,0.62,0.63,0.64};
    float sVb02[]={0.01,0.01,0.01,0.01};

    float Ib02[]={0.096,0.201,0.29,0.392};
    float sIb02[]={0.01,0.01,0.01,0.01};

    TCanvas *c1= new TCanvas("c1", "c1", 200,50,600,400);
    
    // creo una classe TH2F per mostrare piu' set di dati 
    auto histo = new TH2F("histo", "histogram title",4,0.54,0.65,4,0,0.42);
    histo->SetXTitle("V [V]");
    histo->SetYTitle("i [mA]");
      
    c1->cd();
    histo->Draw();

    // mostro i punti del primo set di dati
    auto graphErrors_set1 = new TGraphErrors(4,Vb01,Ib01,sVb01,sIb01);
    graphErrors_set1->SetMarkerSize(0.6);
    graphErrors_set1->SetMarkerStyle(21);

    c1->cd();
    graphErrors_set1->Draw("LP");

    // mostro i punti del secondo set di dati
    auto graphErrors_set2 = new TGraphErrors(4,Vb02,Ib02,sVb02,sIb02);
    graphErrors_set2->SetMarkerSize(0.6);
    graphErrors_set2->SetMarkerStyle(21);
    graphErrors_set2->SetLineColor(5);

    c1->cd();
    graphErrors_set2->Draw("LP");

    // provo ad utilizzare TLegend
    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("The Legend Title"); // option "C" allows to center the header
    legend->AddEntry(graphErrors_set1,"Primo set di dati","l");
    legend->AddEntry(graphErrors_set2,"Secondo set di dati","l");
/*     legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
    legend->AddEntry("gr","Graph with error bars","lep"); */

    c1->cd();
    legend->Draw();





    rebin.cd();
    c1->Write();
    


}