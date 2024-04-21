#include "Config.h"

using namespace std;

int dimension(TH1F *hist, int nbin_hist, float &min) //n= numero di bin che avrà ciascuno dei due istogrammi
{
    int num;
    float width= hist->GetBinWidth(1);
    float bin_ledge[nbin_hist];
    for(int i=1; i<=nbin_hist; i++){
        bin_ledge[i]= hist->GetBinLowEdge(i);
        if( abs(bin_ledge[i]-0.)<width/2){
            min=i;
            num= min-1;
            break;
        }                               
    }

    return num;
}



TH1F* rebin(TH1F *g_hist)
{
    //sto riempendo il mio istogramma rebinnato
    TH1F* h_rebin= new TH1F("h_rebin", "hist rebinned", kNbinsx, k_xbins);     
    for(int i=1; i<=kNbinsx; i++){     //loop su x nuovo
        int xmin= g_hist->GetXaxis()->FindBin(k_xbins[i-1]+epsilon);
        int xmax= g_hist->GetXaxis()->FindBin(k_xbins[i]-epsilon);
        double integral= g_hist->Integral(xmin, xmax);
        h_rebin->SetBinContent(i, integral);
      } 

    return h_rebin;
}




void effWithTpcSegmentMomDependent(){

      //ATTENZIONE: d'ora in poi uso questa maro perché i dati che uso sono sempre momentum-dependent

      //NB: non è la macro che "fa" i tagli momentum-dependent ma i dati contenuti in questo file .root che sono ottenuti con tagli momentum-dependent    
      cout<<"\ntagli momentum-dependent on"<<endl;
      TFile fAnalysis("AnalysisResults_mc_20240328.root");
      TFile effWithTpcSegmentMomDependent("effWithTpcSegmentMomDependent.root", "RECREATE");
      effWithTpcSegmentMomDependent.mkdir("rebinned_histo_and_eff");

      TH2F* tpcSegment= (TH2F*)fAnalysis.Get("efficiency-q-a/tpcSegment");

      TH1F* hasTpcS= (TH1F*)tpcSegment->ProjectionY("hasTpcSegment", 1, 1);
      TH1F* foundTpcS= (TH1F*)tpcSegment->ProjectionY("foundTpcSegment", 2, 2);
      TH1F* allFoundTpcS= (TH1F*)tpcSegment->ProjectionY("allFoundTpcSegment", 3, 3);
      TH1F* foundTpcSWithFake= (TH1F*)tpcSegment->ProjectionY("foundTpcSegmentw/fake", 4, 4);

      gStyle->SetOptStat(0);

      //-------------------------------------------------------------------------------------------//
      int nbin= hasTpcS->GetNbinsX();
      float min;
      int n = dimension(hasTpcS, nbin, min);
      float max= hasTpcS->GetXaxis()->GetBinUpEdge(nbin);

      //------------------------ p+ ------------------------//
      TH1F *hasTpcS_pos= new TH1F("hasTpcS_pos", "hasTPCSegment pi+;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      hasTpcS_pos->SetName("hasTPCSegment_pi_pos");

      hasTpcS_pos->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         hasTpcS_pos->SetBinContent(i, hasTpcS->GetBinContent(min-1+i));
      }

      TCanvas *ChasTpcS_pos= new TCanvas("ChasTpcS_pos", "hasTPCSegment pi+", 200,50,600,400);
      ChasTpcS_pos->cd(0);
      hasTpcS_pos->Draw();

      effWithTpcSegmentMomDependent.cd();
      ChasTpcS_pos->Write();

      
      TH1F* hasTpcS_pos_rebin= rebin(hasTpcS_pos);
      hasTpcS_pos_rebin->SetName("hasTpcS_pos_rebin");
      hasTpcS_pos_rebin->SetTitle("hasTpcS_pos rebinned;#it{p}_{T} (GeV/#it{c});Entries");

      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      hasTpcS_pos_rebin->Write();




      TH1F *foundTpcS_pos= new TH1F("foundTpcS_pos", "foundTPCSegment pi+;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      foundTpcS_pos->SetName("foundTPCSegment_pi_pos");

      foundTpcS_pos->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         foundTpcS_pos->SetBinContent(i, foundTpcS->GetBinContent(min-1+i));
      }

      TCanvas *CfoundTpcS_pos= new TCanvas("CfoundTpcS_pos", "foundTPCSegment pi+", 200,50,600,400);
      CfoundTpcS_pos->cd(0);
      foundTpcS_pos->Draw();

      effWithTpcSegmentMomDependent.cd();   
      CfoundTpcS_pos->Write();


      TH1F* foundTpcS_pos_rebin= rebin(foundTpcS_pos);
      foundTpcS_pos_rebin->SetName("foundTpcS_pos_rebin");
      foundTpcS_pos_rebin->SetTitle("foundTpcS_pos rebinned;#it{p}_{T} (GeV/#it{c});Entries");
      
      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      foundTpcS_pos_rebin->Write();




      TH1F *allFoundTpcS_pos= new TH1F("allFoundTpcS_pos", "allFoundTPCSegment pi+;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      allFoundTpcS_pos->SetName("allFoundTPCSegment_pi_pos");

      allFoundTpcS_pos->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         allFoundTpcS_pos->SetBinContent(i, allFoundTpcS->GetBinContent(min-1+i));
      }

      TCanvas *CallFoundTpcS_pos= new TCanvas("CallFoundTpcS_pos", "allFoundTPCSegment pi+", 200,50,600,400);
      CallFoundTpcS_pos->cd(0);
      allFoundTpcS_pos->Draw();

      effWithTpcSegmentMomDependent.cd();   
      CallFoundTpcS_pos->Write();


      TH1F* allFoundTpcS_pos_rebin= rebin(allFoundTpcS_pos);
      allFoundTpcS_pos_rebin->SetName("allFoundTpcS_pos_rebin");
      allFoundTpcS_pos_rebin->SetTitle("allFoundTpcS_pos rebinned;#it{p}_{T} (GeV/#it{c});Entries");

      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      allFoundTpcS_pos_rebin->Write();




      TH1F *foundTpcSWithFake_pos= new TH1F("foundTpcSWithFake_pos", "foundTPCSegmentWithFake pi+;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      foundTpcSWithFake_pos->SetName("foundTPCSegmentWithFake_pi_pos");

      foundTpcSWithFake_pos->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         foundTpcSWithFake_pos->SetBinContent(i, foundTpcSWithFake->GetBinContent(min-1+i));
      }

      TCanvas *CfoundTpcSWithFake_pos= new TCanvas("CfoundTpcSWithFake_pos", "foundTPCSegmentWithFake pi+", 200,50,600,400);
      CfoundTpcSWithFake_pos->cd(0);
      foundTpcSWithFake_pos->Draw();

      effWithTpcSegmentMomDependent.cd();   
      CfoundTpcSWithFake_pos->Write();


      TH1F* foundTpcSWithFake_pos_rebin= rebin(foundTpcSWithFake_pos);
      foundTpcSWithFake_pos_rebin->SetName("foundTpcSWithFake_pos_rebin");
      foundTpcSWithFake_pos_rebin->SetTitle("foundTpcSWithFake_pos rebinned;#it{p}_{T} (GeV/#it{c});Entries");
      
      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      foundTpcSWithFake_pos_rebin->Write();

      //------------------------ p- ------------------------//
      TH1F *hasTpcS_neg= new TH1F("hasTpcS_neg", "hasTPCSegment pi-;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      hasTpcS_neg->SetName("hasTPCSegment_pi_neg");

      hasTpcS_neg->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         hasTpcS_neg->SetBinContent(i, hasTpcS->GetBinContent(min-i));
      }

      TCanvas *ChasTpcS_neg= new TCanvas("ChasTpcS_neg", "hasTPCSegment pi-", 200,50,600,400);
      ChasTpcS_neg->cd(0);
      hasTpcS_neg->Draw();

      effWithTpcSegmentMomDependent.cd();   
      ChasTpcS_neg->Write();


      TH1F* hasTpcS_neg_rebin= rebin(hasTpcS_neg);
      hasTpcS_neg_rebin->SetName("hasTpcS_neg_rebin");
      hasTpcS_neg_rebin->SetTitle("hasTpcS_neg rebinned;#it{p}_{T} (GeV/#it{c});Entries");
      
      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      hasTpcS_neg_rebin->Write();




      TH1F *foundTpcS_neg= new TH1F("foundTpcS_neg", "foundTPCSegment pi-;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      foundTpcS_neg->SetName("foundTPCSegment_pi_neg");

      foundTpcS_neg->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         foundTpcS_neg->SetBinContent(i, foundTpcS->GetBinContent(min-i));
      }

      TCanvas *CfoundTpcS_neg= new TCanvas("CfoundTpcS_neg", "foundTPCSegment pi-", 200,50,600,400);
      CfoundTpcS_neg->cd(0);
      foundTpcS_neg->Draw();

      effWithTpcSegmentMomDependent.cd();   
      CfoundTpcS_neg->Write();


      TH1F* foundTpcS_neg_rebin= rebin(foundTpcS_neg);
      foundTpcS_neg_rebin->SetName("foundTpcS_neg_rebin");
      foundTpcS_neg_rebin->SetTitle("foundTpcS_neg rebinned;#it{p}_{T} (GeV/#it{c});Entries");
      
      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      foundTpcS_neg_rebin->Write();




      TH1F *allFoundTpcS_neg= new TH1F("allFoundTpcS_neg", "allFoundTPCSegment pi-;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      allFoundTpcS_neg->SetName("allFoundTPCSegment_pi_neg");

      allFoundTpcS_neg->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         allFoundTpcS_neg->SetBinContent(i, allFoundTpcS->GetBinContent(min-i));
      }

      TCanvas *CallFoundTpcS_neg= new TCanvas("CallFoundTpcS_neg", "allFoundTPCSegment pi-", 200,50,600,400);
      CallFoundTpcS_neg->cd(0);
      allFoundTpcS_neg->Draw();

      effWithTpcSegmentMomDependent.cd();   
      CallFoundTpcS_neg->Write();


      TH1F* allFoundTpcS_neg_rebin= rebin(allFoundTpcS_neg);
      allFoundTpcS_neg_rebin->SetName("allFoundTpcS_neg_rebin");
      allFoundTpcS_neg_rebin->SetTitle("allFoundTpcS_neg rebinned;#it{p}_{T} (GeV/#it{c});Entries");
      
      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      allFoundTpcS_neg_rebin->Write();




      TH1F *foundTpcSWithFake_neg= new TH1F("foundTpcSWithFake_neg", "foundTPCSegmentWithFake pi-;#it{p}_{T} (GeV/#it{c});Entries",n,0.,max);
      foundTpcSWithFake_neg->SetName("foundTPCSegmentWithFake_pi_neg");

      foundTpcSWithFake_neg->SetBinContent(1, 0);
      for(int i=2; i<=n; i++){
         foundTpcSWithFake_neg->SetBinContent(i, foundTpcSWithFake->GetBinContent(min-i));
      }

      TCanvas *CfoundTpcSWithFake_neg= new TCanvas("CfoundTpcSWithFake_neg", "foundTPCSegmentWithFake pi-", 200,50,600,400);
      CfoundTpcSWithFake_neg->cd(0);
      foundTpcSWithFake_neg->Draw();

      effWithTpcSegmentMomDependent.cd();   
      CfoundTpcSWithFake_neg->Write();


      TH1F* foundTpcSWithFake_neg_rebin= rebin(foundTpcSWithFake_neg);
      foundTpcSWithFake_neg_rebin->SetName("foundTpcSWithFake_neg_rebin");
      foundTpcSWithFake_neg_rebin->SetTitle("foundTpcSWithFake_neg rebinned;#it{p}_{T} (GeV/#it{c});Entries");
      
      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");   
      foundTpcSWithFake_neg_rebin->Write();


      //------------------------ Calcolo delle efficienze ------------------------//

      //------------------------ p+ ------------------------//
      TCanvas *CfoundOverFindable_pos = new TCanvas("CfoundOverFindable_pos","Efficienza found over findable pi+",200,10,700,500);

      TH1F*bin2_1_pos= new TH1F(*foundTpcS_pos);
      bin2_1_pos->SetName("foundOverFindableEfficiency_pos");
      bin2_1_pos->Divide(foundTpcS_pos, hasTpcS_pos, 1, 1, "B");
      bin2_1_pos->SetTitle("Efficienza found over findable #pi^{+};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_1_pos->GetXaxis()->SetRangeUser(0., 5.);

      CfoundOverFindable_pos->cd();
      bin2_1_pos->Draw("PE");

      //CfoundOverFindable_pos->Print("efficienza found over findable pi+.pdf");




      TCanvas *CfakeFraction_pos = new TCanvas("CfakeFraction_pos","Frazione di fake pi+",200,10,700,500);

      TH1F*bin32_pos= new TH1F(*allFoundTpcS_pos);
      bin32_pos->Add(foundTpcS_pos, -1);
      TH1F*bin32_3_pos= new TH1F(*bin32_pos);
      bin32_3_pos->SetName("fakeFractionEfficiency_pos");
      bin32_3_pos->Divide(bin32_pos, allFoundTpcS_pos, 1, 1, "B");
      bin32_3_pos->SetTitle("Efficienza frazione di fake #pi^{+};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin32_3_pos->GetXaxis()->SetRangeUser(0., 5.);

      CfakeFraction_pos->cd();
      bin32_3_pos->Draw("PE");

      //CfakeFraction_pos->Print("efficienza frazione di fake pi+.pdf");




      TCanvas *CtrueFraction_pos = new TCanvas("CtrueFraction_pos","Frazione di true pi+",200,10,700,500);

      TH1F*bin2_3_pos= new TH1F(*foundTpcS_pos);
      bin2_3_pos->SetName("trueFractionEfficiency_pos");
      bin2_3_pos->Divide(foundTpcS_pos, allFoundTpcS_pos, 1, 1, "B");
      bin2_3_pos->SetTitle("Efficienza frazione di true #pi^{+};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_3_pos->GetXaxis()->SetRangeUser(0., 5.);

      CtrueFraction_pos->cd();
      bin2_3_pos->Draw("PE");

      //CtrueFraction_pos->Print("efficienza frazione di true pi+.pdf");




      TCanvas *Cfake_and_trueFraction_pos = new TCanvas("Cfake_and_trueFraction_pos","Confronto tra frazione di fake e frazione di true pi+",200,10,700,500);

      bin32_3_pos->SetTitle("Confronto tra frazione di fake e frazione di true #pi^{+}");
      bin32_3_pos->SetLineColor(2);
      bin2_3_pos->SetLineColor(3);
      bin32_3_pos->GetYaxis()->SetRangeUser(-0.05, 1.05);

      Cfake_and_trueFraction_pos->cd();
      bin32_3_pos->Draw("PE");
      bin2_3_pos->Draw("PESAME");

      TLatex t1;
      t1.SetNDC();
      t1.SetTextSize(0.031);
      t1.DrawLatex(0.5, 0.48, "Confronto tra frazione di fake e frazione di true per #font[22]{#pi^{+}}");  //un po’ centrato
      t1.DrawLatex(0.5, 0.43, "#font[22]{#color[2]{#scale[2.2]{-}}} Frazione di fake");
      t1.DrawLatex(0.5, 0.38, "#font[22]{#color[3]{#scale[2.2]{-}}} Frazione di true");

      effWithTpcSegmentMomDependent.cd();
      CfoundOverFindable_pos->Write();
      CfakeFraction_pos->Write();
      CtrueFraction_pos->Write();
      Cfake_and_trueFraction_pos->Write();
      bin2_1_pos->Write();
      //Cfake_and_trueFraction_pos->Print("confronto frazione di fake e frazione di true pi+.pdf");



      //------------------------ p- ------------------------//
      TCanvas *CfoundOverFindable_neg = new TCanvas("CfoundOverFindable_neg","Efficienza found over findable pi-",200,10,700,500);

      TH1F*bin2_1_neg= new TH1F(*foundTpcS_neg);
      bin2_1_neg->SetName("foundOverFoundableEfficiency_neg");
      bin2_1_neg->Divide(foundTpcS_neg, hasTpcS_neg, 1, 1, "B");
      bin2_1_neg->SetTitle("Efficienza found over findable #pi^{-};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_1_neg->GetXaxis()->SetRangeUser(0., 5.);

      CfoundOverFindable_neg->cd();
      bin2_1_neg->Draw("PE");

      //CfoundOverFindable_neg->Print("efficienza found over findable pi-.pdf");




      TCanvas *CfakeFraction_neg = new TCanvas("CfakeFraction_neg","Efficienza frazione di fake pi-",200,10,700,500);

      TH1F*bin32_neg= new TH1F(*allFoundTpcS_neg);
      bin32_neg->Add(foundTpcS_neg, -1);
      TH1F*bin32_3_neg= new TH1F(*bin32_neg);
      bin32_3_neg->SetName("fakeFractionEfficiency_neg");
      bin32_3_neg->Divide(bin32_neg, allFoundTpcS_neg, 1, 1, "B");
      bin32_3_neg->SetTitle("Efficienza frazione di fake #pi^{-};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin32_3_neg->GetXaxis()->SetRangeUser(0., 5.);

      CfakeFraction_neg->cd();
      bin32_3_neg->Draw("PE");

      //CfakeFraction_neg->Print("efficienza frazione di fake pi-.pdf");




      TCanvas *CtrueFraction_neg = new TCanvas("CtrueFraction_neg","Efficienza frazione di true pi-",200,10,700,500);

      TH1F*bin2_3_neg= new TH1F(*foundTpcS_neg);
      bin2_3_neg->SetName("trueFractionEfficiency_neg");
      bin2_3_neg->Divide(foundTpcS_neg, allFoundTpcS_neg, 1, 1, "B");
      bin2_3_neg->SetTitle("Efficienza frazione di true #pi^{-};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_3_neg->GetXaxis()->SetRangeUser(0., 5.);

      CtrueFraction_neg->cd();
      bin2_3_neg->Draw("PE");

      //CtrueFraction_neg->Print("efficienza frazione di true pi-.pdf");




      TCanvas *Cfake_and_trueFraction_neg = new TCanvas("Cfake_and_trueFraction_neg","Confronto tra frazione di fake e frazione di true pi-",200,10,700,500);

      bin32_3_neg->SetTitle("Confronto tra frazione di fake e frazione di true #pi^{-}");
      bin32_3_neg->SetLineColor(2);
      bin2_3_neg->SetLineColor(3);
      bin32_3_neg->GetYaxis()->SetRangeUser(-0.05, 1.05);

      Cfake_and_trueFraction_neg->cd();
      bin32_3_neg->Draw("PE");
      bin2_3_neg->Draw("PESAME");

      TLatex t2;
      t2.SetNDC();
      t2.SetTextSize(0.031);
      t2.DrawLatex(0.5, 0.48, "Confronto tra frazione di fake e frazione di true per #font[22]{#pi^{-}}");  //un po’ centrato
      t2.DrawLatex(0.5, 0.43, "#font[22]{#color[2]{#scale[2.2]{-}}} Frazione di fake");
      t2.DrawLatex(0.5, 0.38, "#font[22]{#color[3]{#scale[2.2]{-}}} Frazione di true");

      effWithTpcSegmentMomDependent.cd();
      CfoundOverFindable_neg->Write();
      CfakeFraction_neg->Write();
      CtrueFraction_neg->Write();
      Cfake_and_trueFraction_neg->Write();
      bin2_1_neg->Write();
      //Cfake_and_trueFraction_neg->Print("confronto frazione di fake e frazione di true pi-.pdf");



      //----------------------------------- Calcolo delle efficienze con gli istogrammi rebinnati -----------------------------------//

      //----------------------------------- p+ -----------------------------------//
      TCanvas *CfoundOverFindable_pos_rebin = new TCanvas("CfoundOverFindable_pos_rebin","Efficienza found over findable_rebin pi+",200,10,700,500);

      TH1F*bin2_1_pos_rebin= new TH1F(*foundTpcS_pos_rebin);
      bin2_1_pos_rebin->SetName("foundOverFindableEfficiency_pos_rebin");
      bin2_1_pos_rebin->Divide(foundTpcS_pos_rebin, hasTpcS_pos_rebin, 1, 1, "B");
      bin2_1_pos_rebin->SetTitle("Efficienza dell'algoritmo #pi^{+};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_1_pos_rebin->GetXaxis()->SetRangeUser(0., 4.);
      bin2_1_pos_rebin->GetYaxis()->SetRangeUser(0., 1.1);

      CfoundOverFindable_pos_rebin->cd();
      bin2_1_pos_rebin->Draw("PE");
      
      CfoundOverFindable_pos_rebin->Print("efficienza found over findable pos rebinned.pdf");




      TCanvas *CfakeFraction_pos_rebin = new TCanvas("CfakeFraction_pos_rebin","Frazione di fake rebinned pi+",200,10,700,500);

      TH1F*bin32_pos_rebin= new TH1F(*allFoundTpcS_pos_rebin);
      bin32_pos_rebin->Add(foundTpcS_pos_rebin, -1);
      TH1F*bin32_3_pos_rebin= new TH1F(*bin32_pos_rebin);
      bin32_3_pos_rebin->SetName("fakeFractionEfficiency_pos_rebin");
      bin32_3_pos_rebin->Divide(bin32_pos_rebin, allFoundTpcS_pos_rebin, 1, 1, "B");
      bin32_3_pos_rebin->SetTitle("Frazione di fake rebinned #pi^{+};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin32_3_pos_rebin->GetXaxis()->SetRangeUser(0., 4.);
      bin32_3_pos_rebin->GetYaxis()->SetRangeUser(0., 1.);

      CfakeFraction_pos_rebin->cd();
      bin32_3_pos_rebin->Draw("PE");

      CfakeFraction_pos_rebin->Print("efficienza frazione di fake rebinned pi+.pdf");




      TCanvas *CtrueFraction_pos_rebin = new TCanvas("CtrueFraction_pos_rebin","Frazione di true rebinned pi+",200,10,700,500);

      TH1F*bin2_3_pos_rebin= new TH1F(*foundTpcS_pos_rebin);
      bin2_3_pos_rebin->SetName("trueFractionEfficiency_pos_rebin");
      bin2_3_pos_rebin->Divide(foundTpcS_pos_rebin, allFoundTpcS_pos_rebin, 1, 1, "B");
      bin2_3_pos_rebin->SetTitle("Frazione di true rebinned #pi^{+};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_3_pos_rebin->GetXaxis()->SetRangeUser(0., 4.);
      bin2_3_pos_rebin->GetYaxis()->SetRangeUser(0., 1.);

      CtrueFraction_pos_rebin->cd();
      bin2_3_pos_rebin->Draw("PE");

      //CtrueFraction_pos_rebin->Print("efficienza frazione di true rebinned pi+.pdf");




      TCanvas *Cfake_and_trueFraction_pos_rebin = new TCanvas("Cfake_and_trueFraction_pos_rebin","Confronto tra frazione di fake e frazione di true rebinned pi+",200,10,700,500);

      bin32_3_pos_rebin->SetTitle("Confronto tra frazione di fake e frazione di true #pi^{+}");
      bin32_3_pos_rebin->SetLineColor(2);
      bin2_3_pos_rebin->SetLineColor(kGreen+2);
      bin32_3_pos_rebin->GetXaxis()->SetRangeUser(0, 4);
      bin32_3_pos_rebin->GetYaxis()->SetRangeUser(-0.05, 1.15);

      Cfake_and_trueFraction_pos_rebin->cd();
      bin32_3_pos_rebin->Draw("PE0");
      bin2_3_pos_rebin->Draw("PESAME");

      TLatex t1_rebin;
      t1_rebin.SetNDC();
      t1_rebin.SetTextSize(0.031);
      t1_rebin.DrawLatex(0.4, 0.48, "Confronto tra frazione di fake e frazione di true per #font[22]{#pi^{+}}");  //un po’ centrato
      t1_rebin.DrawLatex(0.4, 0.43, "#font[22]{#color[2]{#scale[2.2]{-}}} Frazione di fake");
      t1_rebin.DrawLatex(0.4, 0.38, "#font[22]{#color[8]{#scale[2.2]{-}}} Frazione di true");

      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");
      CfoundOverFindable_pos_rebin->Write();
      CfakeFraction_pos_rebin->Write();
      CtrueFraction_pos_rebin->Write();
      Cfake_and_trueFraction_pos_rebin->Write();
      bin2_1_pos_rebin->Write();
      Cfake_and_trueFraction_pos_rebin->Print("confronto frazione di fake e frazione di true rebinned pi+.pdf");



      //------------------------ p- ------------------------//
      TCanvas *CfoundOverFindable_neg_rebin = new TCanvas("CfoundOverFindable_neg_rebin","Efficienza found over findable_rebin pi-",200,10,700,500);

      TH1F*bin2_1_neg_rebin= new TH1F(*foundTpcS_neg_rebin);
      bin2_1_neg_rebin->SetName("foundOverFindableEfficiency_neg_rebin");
      bin2_1_neg_rebin->Divide(foundTpcS_neg_rebin, hasTpcS_neg_rebin, 1, 1, "B");
      bin2_1_neg_rebin->SetTitle("Efficienza dell'algoritmo #pi^{-};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_1_neg_rebin->GetXaxis()->SetRangeUser(0., 5.);
      bin2_1_neg_rebin->GetYaxis()->SetRangeUser(0., 1.);

      CfoundOverFindable_neg_rebin->cd();
      bin2_1_neg_rebin->Draw("PE");




      TCanvas *CfakeFraction_neg_rebin = new TCanvas("CfakeFraction_neg_rebin","Efficienza frazione di fake rebinned pi-",200,10,700,500);

      TH1F*bin32_neg_rebin= new TH1F(*allFoundTpcS_neg_rebin);
      bin32_neg_rebin->Add(foundTpcS_neg_rebin, -1);
      TH1F*bin32_3_neg_rebin= new TH1F(*bin32_neg_rebin);
      bin32_3_neg_rebin->SetName("fakeFractionEfficiency_neg_rebin");
      bin32_3_neg_rebin->Divide(bin32_neg_rebin, allFoundTpcS_neg_rebin, 1, 1, "B");
      bin32_3_neg_rebin->SetTitle("Frazione di fake rebinned #pi^{-};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin32_3_neg_rebin->GetXaxis()->SetRangeUser(0., 5.);
      bin32_3_neg_rebin->GetYaxis()->SetRangeUser(0., 1.);

      CfakeFraction_neg_rebin->cd();
      bin32_3_neg_rebin->Draw("PE");

      //CfakeFraction_neg_rebin->Print("efficienza frazione di fake rebinned pi-.pdf");




      TCanvas *CtrueFraction_neg_rebin = new TCanvas("CtrueFraction_neg_rebin","Efficienza frazione di true rebinned pi-",200,10,700,500);

      TH1F*bin2_3_neg_rebin= new TH1F(*foundTpcS_neg_rebin);
      bin2_3_neg_rebin->SetName("trueFractionEfficiency_neg_rebin");
      bin2_3_neg_rebin->Divide(foundTpcS_neg_rebin, allFoundTpcS_neg_rebin, 1, 1, "B");
      bin2_3_neg_rebin->SetTitle("Frazione di true rebinned #pi^{-};#it{p}_{T} (GeV/#it{c});#epsilon");
      bin2_3_neg_rebin->GetXaxis()->SetRangeUser(0., 5.);
      bin2_3_neg_rebin->GetYaxis()->SetRangeUser(0., 1.);

      CtrueFraction_neg_rebin->cd();
      bin2_3_neg_rebin->Draw("PE");

      //CtrueFraction_neg_rebin->Print("efficienza frazione di true rebinned pi-.pdf");




      TCanvas *Cfake_and_trueFraction_neg_rebin = new TCanvas("Cfake_and_trueFraction_neg_rebin","Confronto tra frazione di fake e frazione di true rebinned pi-",200,10,700,500);

      bin32_3_neg_rebin->SetTitle("Confronto tra frazione di fake e frazione di true rebinned #pi^{-}");
      bin32_3_neg_rebin->SetLineColor(2);
      bin2_3_neg_rebin->SetLineColor(3);
      bin32_3_neg_rebin->GetYaxis()->SetRangeUser(-0.05, 1.05);

      Cfake_and_trueFraction_neg_rebin->cd();
      bin32_3_neg_rebin->Draw("PE");
      bin2_3_neg_rebin->Draw("PESAME");

      TLatex t2_rebin;
      t2_rebin.SetNDC();
      t2_rebin.SetTextSize(0.031);
      t2_rebin.DrawLatex(0.4, 0.48, "Confronto tra frazione di fake e frazione di true per #font[22]{#pi^{-}} rebinned histo");  //un po’ centrato
      t2_rebin.DrawLatex(0.4, 0.43, "#font[22]{#color[2]{#scale[2.2]{-}}} Frazione di fake");
      t2_rebin.DrawLatex(0.4, 0.38, "#font[22]{#color[3]{#scale[2.2]{-}}} Frazione di true");

      effWithTpcSegmentMomDependent.cd("rebinned_histo_and_eff");
      CfoundOverFindable_neg_rebin->Write();
      CfakeFraction_neg_rebin->Write();
      CtrueFraction_neg_rebin->Write();
      Cfake_and_trueFraction_neg_rebin->Write();
      bin2_1_neg_rebin->Write();
      //Cfake_and_trueFraction_neg_rebin->Print("confronto frazione di fake e frazione di true rebinned pi-.pdf");




      //------------------------------------------ Calcolo efficienze con bin4 () al posto di bin2 ------------------------------------------//
      //se mi serve, devo copiarlo dal file effWithTpcSegment
   



      //--------------------------------- Plot con efficienze per pi+ e pi- con istogrammi rebinnati ---------------------------------//
      TCanvas *Ceff_algoritmo = new TCanvas("Ceff_algoritmo","Confronto tra efficienza found over findable pi+ e p-",200,10,700,500);

      bin2_1_pos_rebin->SetTitle("Confronto tra efficienza found over findable #pi^{+} e #pi^{-}");
      bin2_1_pos_rebin->SetLineColor(2);
      bin2_1_neg_rebin->SetLineColor(4);
      bin2_1_pos_rebin->GetYaxis()->SetRangeUser(-0.05, 1.10);

      Ceff_algoritmo->cd();
      bin2_1_pos_rebin->Draw("PE");
      bin2_1_neg_rebin->Draw("PESAME");

      TLatex t3;
      t3.SetNDC();
      t3.SetTextSize(0.031);
      t3.DrawLatex(0.5, 0.48, "Confronto tra efficienze dell'algoritmo");  //un po’ centrato
      t3.DrawLatex(0.5, 0.43, "#font[22]{#color[2]{#scale[2.2]{-}}} #font[22]{#pi^{+}}");
      t3.DrawLatex(0.5, 0.38, "#font[22]{#color[4]{#scale[2.2]{-}}} #font[22]{#pi^{-}}");

      Ceff_algoritmo->Print("confronto efficienze pound over findable per pi+ e pi-.pdf");




      TCanvas *Cfake = new TCanvas("Cfake","Confronto tra frazione di fake per pi+ e p-",200,10,700,500);

      bin32_3_pos_rebin->SetTitle("Confronto tra frazione di fake per #pi^{+} e #pi^{-}");
      bin32_3_pos_rebin->SetLineColor(2);
      bin32_3_neg_rebin->SetLineColor(4);
      bin32_3_pos_rebin->GetYaxis()->SetRangeUser(-0.05, 1.05);

      Cfake->cd();
      bin32_3_pos_rebin->Draw("PE");
      bin32_3_neg_rebin->Draw("PESAME");

      TLatex t4;
      t4.SetNDC();
      t4.SetTextSize(0.031);
      t4.DrawLatex(0.5, 0.48, "Confronto tra frazioni di fake");  //un po’ centrato
      t4.DrawLatex(0.5, 0.43, "#font[22]{#color[2]{#scale[2.2]{-}}} #font[22]{#pi^{+}}");
      t4.DrawLatex(0.5, 0.38, "#font[22]{#color[4]{#scale[2.2]{-}}} #font[22]{#pi^{-}}");

      Cfake->Print("confronto frazioni di fake per pi+ e pi-.pdf");

   effWithTpcSegmentMomDependent.Close();
}