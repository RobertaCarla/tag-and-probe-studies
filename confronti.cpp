using namespace std;

void confronti(){
    TFile effWithTpcSegment("effWithTpcSegment.root");
    TFile effWithTpcSegmentMomDependent("effWithTpcSegmentMomDependent.root");
    TFile confronti("confronti.root","RECREATE");
 
    TH1F *foundOverFindable_pos= (TH1F*)effWithTpcSegment.Get("rebinned_histo_and_eff/foundOverFindableEfficiency_pos_rebin");
    TH1F *foundOverFindable_pos_MomDep= (TH1F*)effWithTpcSegmentMomDependent.Get("rebinned_histo_and_eff/foundOverFindableEfficiency_pos_rebin");
    
    TH1F *foundOverFindable_neg= (TH1F*)effWithTpcSegment.Get("rebinned_histo_and_eff/CfoundOverFindableEfficiency_neg_rebin");
    TH1F *foundOverFindable_neg_MomDep= (TH1F*)effWithTpcSegmentMomDependent.Get("rebinned_histo_and_eff/foundOverFindableEfficiency_neg_rebin");

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    cout<<"n: "<<foundOverFindable_pos->GetNbinsX()<<endl;

    TCanvas *Ccomparison_foundOverFindable_pos= new TCanvas("Ccomparison_algorithm_pos","Confronto  tra efficienze dell'algoritmo di rematching pi+",200,10,700,500);
    foundOverFindable_pos->SetLineColor(2);
    foundOverFindable_pos->GetXaxis()->SetRangeUser(0., 4.);
    foundOverFindable_pos->GetYaxis()->SetRangeUser(0., 1.1);

    foundOverFindable_pos->SetTitle("Confronto ta efficienze dell'algoritmo di rematching #pi^{+}");

    Ccomparison_foundOverFindable_pos->cd();
    foundOverFindable_pos->Draw("PE");
    foundOverFindable_pos_MomDep->Draw("PESAME");

   TLatex t1;
    t1.SetNDC();
    t1.SetTextSize(0.031);
    t1.DrawLatex(0.21, 0.265, "#scale[1.25]{Confronto tra le efficienze dell'algoritmo di rematching per #font[22]{#pi^{+}}}");  //un po’ centrato
    t1.DrawLatex(0.21, 0.215, "#scale[1.25]{#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza prima dell'ottimizzazione}");
    t1.DrawLatex(0.21, 0.165, "#scale[1.25]{#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza dopo l'ottimizzazione}");
 
    confronti.cd();
    Ccomparison_foundOverFindable_pos->Write();
    Ccomparison_foundOverFindable_pos->Print("confronto tra efficienze dell'algoritmo.pdf");

    confronti.Close();
}