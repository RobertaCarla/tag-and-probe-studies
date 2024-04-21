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




void efficiencyMomDependent(){

        //ATTENZIONE: d'ora in poi uso questa maro perché i dati che uso sono sempre momentum-dependent
    
        cout<<"\ntagli momentum-dependent on"<<endl;
        //per tagli momentum-depentent per pT, eta e phi
        //NB: non è la macro che "fa" i tagli momentum-dependent ma i dati contenuti in questo file .root che sono ottenuti con tagli momentum-dependent
        TFile fAnalysis("AnalysisResults_mc_20240328.root");
        TFile efficiencyMomDependent("efficiencyMomDependent.root", "RECREATE");
        efficiencyMomDependent.mkdir("eff_with_rebinned_histos");

        TH3F *piRec= (TH3F*)fAnalysis.Get("efficiency-q-a/piRec");
 
        TH2F* piRec_xy=(TH2F*)piRec->Project3D("yx");
 
        TH1F *Decays= (TH1F*)piRec->ProjectionY("Decays",1,1,1,900);
        TH1F *ITS_TPC= (TH1F*)piRec->ProjectionY("ITS_TPC",6,6,1,900);
        TH1F *ITSwithTPC= (TH1F*)piRec->ProjectionY("ITS with TPC",8,8,1,900);

        TH2F*tpcSegment= (TH2F*)fAnalysis.Get("efficiency-q-a/tpcSegment");
        TH1F*hasTpcSegment= (TH1F*)tpcSegment->ProjectionY("hasTpcSegment", 1, 1);
 
        gStyle->SetOptStat(0);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);



    //--------------------------------------------------------- decays ---------------------------------------------------------//
        int nbin= Decays->GetNbinsX();
        float min;
        int n = dimension(Decays, nbin, min);

    //------------------- decays pi+ -----------------//
        float massimox= Decays->GetXaxis()->GetBinUpEdge(nbin);
        TH1F *pi_pos_d= new TH1F("pi_pos_d", "pi+ decays;#it{p}_{T} (GeV/#it{c}); Entries",n,0.,massimox);
        pi_pos_d->SetName("decays_pi_pos");
   
        pi_pos_d->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_pos_d->SetBinContent(i, Decays->GetBinContent(min-1+i));
        }

        TCanvas *Cpi_pos_d= new TCanvas("Cpi_pos_d", "pi+ decays", 200,50,600,400);
        Cpi_pos_d->cd(0);
        pi_pos_d->Draw();

        efficiencyMomDependent.cd();   
        Cpi_pos_d->Write();
        //Cpi_pos_d->Print("efficiency_tagli_mom_dependent/decay pi+.pdf");

    //------------------- decays pi+ rebinned -----------------//
        TH1F* pi_pos_d_rebin= rebin(pi_pos_d);
        pi_pos_d_rebin->SetName("pi_pos_d_rebin");
        pi_pos_d_rebin->SetTitle("pi_pos_d rebinned");
        pi_pos_d_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_pos_d_rebin->GetYaxis()->SetTitle("Entries");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_pos_d_rebin->Write(); 
    
    //------------- decays pi- -----------------------//
        TH1F *pi_neg_d= new TH1F("pi_neg_d", "pi- decays;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_neg_d->SetName("decays_pi_neg");

        pi_neg_d->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_neg_d->SetBinContent(i, Decays->GetBinContent(min-i));
        }

        TCanvas *Cpi_neg_d= new TCanvas("Cpi_neg_d", "pi- decays", 200,50,600,400);
        Cpi_neg_d->cd(0);
        pi_neg_d->Draw(); 
 
        efficiencyMomDependent.cd();   
        Cpi_neg_d->Write();
        //Cpi_neg_d->Print("efficiency_tagli_mom_dependent/decay pi-.pdf");
   
    //------------------- decays pi- rebinned -----------------//
        TH1F* pi_neg_d_rebin= rebin(pi_neg_d);
        pi_neg_d_rebin->SetName("pi_neg_d_rebin");
        pi_neg_d_rebin->SetTitle("pi_neg_d rebinned");
        pi_neg_d_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_neg_d_rebin->GetYaxis()->SetTitle("Entries");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_neg_d_rebin->Write();  


    //--------------------------------------------------------- ITS+TPC ---------------------------------------------------------//
    
    //------------------- ITS+TPC pi+ -----------------//
        TH1F *pi_pos_IT= new TH1F("pi_pos_ITS_TPC", "pi+ ITS+TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_pos_IT->SetName("ITS+TPC pi_pos");

        pi_pos_IT->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_pos_IT->SetBinContent(i, ITS_TPC->GetBinContent(min-1+i));
        }

        TCanvas *Cpi_pos_IT= new TCanvas("Cpi_pos_IT", "pi+ ITS+TPC", 200,50,600,400);
        Cpi_pos_IT->cd(0);
        pi_pos_IT->Draw(); 

        efficiencyMomDependent.cd();   
        Cpi_pos_IT->Write();
        //Cpi_pos_IT->Print("efficiency_tagli_mom_dependent/ITS+TPC pi+.pdf");
        
    //------------------- ITS+TPC pi+ rebinned -----------------//
        TH1F* pi_pos_IT_rebin= rebin(pi_pos_IT);
        pi_pos_IT_rebin->SetName("pi_pos_IT_rebin");
        pi_pos_IT_rebin->SetTitle("pi_pos_IT rebinned");
        pi_pos_IT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_pos_IT_rebin->GetYaxis()->SetTitle("Entries");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_pos_IT_rebin->Write();

    //------------- ITS+TPC pi- -----------------------//
        TH1F *pi_neg_IT= new TH1F("pi_neg_ITS_TPC", "pi- ITS+TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_neg_IT->SetName("ITS+TPC pi_neg");

        pi_neg_IT->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_neg_IT->SetBinContent(i, ITS_TPC->GetBinContent(min-i));
        }

        TCanvas *Cpi_neg_IT= new TCanvas("Cpi_neg_IT", "pi- ITS+TPC", 200,50,600,400);
        Cpi_neg_IT->cd(0);
        pi_neg_IT->Draw(); 

        efficiencyMomDependent.cd();   
        Cpi_neg_IT->Write(); 
        //Cpi_neg_IT->Print("efficiency_tagli_mom_dependent/ITS+TPC pi-.pdf");  

    //------------------- ITS+TPC pi- rebinned -----------------//
        TH1F* pi_neg_IT_rebin= rebin(pi_neg_IT);
        pi_neg_IT_rebin->SetName("pi_neg_IT_rebin");
        pi_neg_IT_rebin->SetTitle("pi_neg_IT rebinned");
        pi_neg_IT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_neg_IT_rebin->GetYaxis()->SetTitle("Entries");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_neg_IT_rebin->Write();


    //--------------------------------------------------------- ITS with TPC ---------------------------------------------------------//
    
    //------------------- ITS with TPC pi+ -----------------//
        TH1F *pi_pos_W= new TH1F("pi_pos_ITSwithTPC", "pi+ ITS with TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_pos_W->SetName("ITS with TPC pi_pos");

        pi_pos_W->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_pos_W->SetBinContent(i, ITSwithTPC->GetBinContent(min-1+i));
        }

        TCanvas *Cpi_pos_W= new TCanvas("Cpi_pos_W", "pi+ ITS with TPC", 200,50,600,400);
        Cpi_pos_W->cd(0);
        pi_pos_W->Draw();

        efficiencyMomDependent.cd();   
        Cpi_pos_W->Write();
        //Cpi_pos_W->Print("efficiency_tagli_mom_dependent/ITS with TPC pi+.pdf");

    //------------------- ITS with TPC pi+ rebinned -----------------//
        TH1F* pi_pos_W_rebin= rebin(pi_pos_W);
        pi_pos_W_rebin->SetName("pi_pos_W_rebin");
        pi_pos_W_rebin->SetTitle("pi_pos_W rebinned");
        pi_pos_W_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_pos_W_rebin->GetYaxis()->SetTitle("Entries");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_pos_W_rebin->Write();  
            
    //------------- ITS with TPC pi- -----------------------//
        TH1F *pi_neg_W= new TH1F("pi_neg_ITSwithTPC", "pi- ITS with TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_neg_W->SetName("ITS with TPC pi_neg");

        pi_neg_W->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_neg_W->SetBinContent(i, ITSwithTPC->GetBinContent(min-i));
        }

        TCanvas *Cpi_neg_W= new TCanvas("Cpi_neg_W", "pi- ITS with TPC", 200,50,600,400);
        Cpi_neg_W->cd(0);
        pi_neg_W->Draw();

        efficiencyMomDependent.cd();   
        Cpi_neg_W->Write(); 
        //Cpi_neg_W->Print("efficiency_tagli_mom_dependent/ITS with TPC pi-.pdf");
        
    //------------------- ITS with TPC pi- rebinned -----------------//
        TH1F* pi_neg_W_rebin= rebin(pi_neg_W);
        pi_neg_W_rebin->SetName("pi_neg_W_rebin");
        pi_neg_W_rebin->SetTitle("pi_neg_W rebinned");
        pi_neg_W_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_neg_W_rebin->GetYaxis()->SetTitle("Entries");
 
        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_neg_W_rebin->Write();
        

    //--------------------------------------------------------- TPC segment with MC label ---------------------------------------------------------//

    //------------------- hasTPCSegment pi+ -----------------//
        TH1F *pi_pos_TpcSegment= new TH1F("pi_pos_TpcSegment", "pi+ hasTpcSegment MC label;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_pos_TpcSegment->SetName("hasTpcSegment_pi_pos");
   
        pi_pos_TpcSegment->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_pos_TpcSegment->SetBinContent(i, hasTpcSegment->GetBinContent(min-1+i));
        }

        TCanvas *Cpi_pos_TpcSegment= new TCanvas("Cpi_pos_TpcSegment", "pi+ hasTpcSegment MC label", 200,50,600,400);
        Cpi_pos_TpcSegment->cd(0);
        pi_pos_TpcSegment->Draw();

        efficiencyMomDependent.cd();   
        Cpi_pos_TpcSegment->Write();
        //Cpi_pos_TpcSegment->Print("efficiency_tagli_mom_dependent/hasTpcSegment pi+.pdf");

    //------------------- hasTPCSegment pi+ rebinned -----------------//
        TH1F* pi_pos_TpcSegment_rebin= rebin(pi_pos_TpcSegment);
        pi_pos_TpcSegment_rebin->SetName("pi_pos_TpcSegment_rebin");
        pi_pos_TpcSegment_rebin->SetTitle("pi_pos_TpcSegment rebinned");
        pi_pos_TpcSegment_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_pos_TpcSegment_rebin->GetYaxis()->SetTitle("Entries");
 
        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_pos_TpcSegment_rebin->Write();
    
    //------------- hasTPCSegment pi- -----------------------//
        TH1F *pi_neg_TpcSegment= new TH1F("pi_neg_TpcSegment", "pi- hasTpsSegment MC label;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        pi_neg_TpcSegment->SetName("hasTpcSegment_pi_neg");

        pi_neg_TpcSegment->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        pi_neg_TpcSegment->SetBinContent(i, hasTpcSegment->GetBinContent(min-i));
        }

        TCanvas *Cpi_neg_TpcSegment= new TCanvas("Cpi_neg_TpcSegment", "pi- hasTpcSegment MC label", 200,50,600,400);
        Cpi_neg_TpcSegment->cd(0);
        pi_neg_TpcSegment->Draw();
 
        efficiencyMomDependent.cd();   
        Cpi_neg_TpcSegment->Write();
        //Cpi_neg_TpcSegment->Print("efficiency_tagli_mom_dependent/hasTpcSegment pi-.pdf");

    //------------------- hasTPCSegment pi- rebinned -----------------//
        TH1F* pi_neg_TpcSegment_rebin= rebin(pi_neg_TpcSegment);
        pi_neg_TpcSegment_rebin->SetName("pi_neg_TpcSegment_rebin");
        pi_neg_TpcSegment_rebin->SetTitle("pi_neg_TpcSegment rebinned");
        pi_neg_TpcSegment_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        pi_neg_TpcSegment_rebin->GetYaxis()->SetTitle("Entries");
 
        efficiencyMomDependent.cd("eff_with_rebinned_histos");   
        pi_neg_TpcSegment_rebin->Write(); 



    //------------------------- calcolo efficienze -------------------------//
    //------------------------- calcoli per p+ -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC --------------------//
        TCanvas *CE_glob_pos = new TCanvas("CE_glob_pos","Efficienza globale per traccia ITS+TPC pi+",200,10,700,500);

        TH1F *gE_glob_pos = new TH1F(*pi_pos_IT);
        gE_glob_pos->SetName("efficiency_pi_pos_glob");
        gE_glob_pos->Divide(pi_pos_IT, pi_pos_d, 1., 1., "B");
        gE_glob_pos->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{+}");
        gE_glob_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_glob_pos->GetYaxis()->SetTitle("#epsilon");
    
        gE_glob_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gE_glob_pos->GetXaxis()->SetRangeUser(0., 3.);

        CE_glob_pos->cd();
        gE_glob_pos->Draw("PE");

        efficiencyMomDependent.cd();
        CE_glob_pos->Write();
        CE_glob_pos->Print("efficiency_tagli_mom_dependent/Efficienza globale per traccia ITS+TPC pi+.pdf");

    //------------------- Efficienza di tracciamento TPC SA --------------------//
        TCanvas *CE_tpc_pos = new TCanvas("CE_tpc_pos","Efficienza di tracciamento TPC SA pi+",200,10,700,500);

        TH1F *Copy_pi_pos_IT= new TH1F(*pi_pos_IT);
        Copy_pi_pos_IT->SetName("copia_pi_pos_IT");
        Copy_pi_pos_IT->Add(pi_pos_W);

        TH1F *gE_tpc_pos= new TH1F(*Copy_pi_pos_IT);
        gE_tpc_pos->SetName("tpc_efficiency_pi_pos");
        gE_tpc_pos->Divide(Copy_pi_pos_IT, pi_pos_d, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_pos->SetTitle("Efficienza di tracciamento TPC SA #pi^{+}");
        gE_tpc_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_pos->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_pos->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così
    

        CE_tpc_pos->cd();
        gE_tpc_pos->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori

        efficiencyMomDependent.cd();
        CE_tpc_pos->Write();
        CE_tpc_pos->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA pi+.pdf");

    //------------------- efficienza matching its-tpc --------------------//
        TCanvas *CE_its_tpc_pos = new TCanvas("CE_its_tpc_pos","Efficienza matching ITS TPC pi+",200,10,700,500);

        TH1F *Copy_pi_pos_W= new TH1F(*pi_pos_W);
        Copy_pi_pos_W->SetName("copia_pi_pos_W");
        Copy_pi_pos_W->Add(pi_pos_IT);

        TH1F *gE_its_tpc_pos= new TH1F(*Copy_pi_pos_W);
        gE_its_tpc_pos->SetName("its_tpc_efficiency_pi_pos");
        gE_its_tpc_pos->Divide(pi_pos_IT, Copy_pi_pos_W, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_pos->SetName("Efficienza_matching_ITS_TPC_pi_os");
        gE_its_tpc_pos->SetTitle("Efficienza matching ITS TPC #pi^{+}");
        gE_its_tpc_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_pos->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_pos->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così
    
        CE_its_tpc_pos->cd();
        gE_its_tpc_pos->Draw("PE");

        efficiencyMomDependent.cd();
        CE_its_tpc_pos->Write();
        CE_its_tpc_pos->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC pi+.pdf");


    //------------------------- Calcoli per pi- -------------------------//

    //------------------- efficienza gloabale --------------------//
        TCanvas *CE_glob_neg = new TCanvas("CE_glob_neg","Efficienza globale per traccia ITS+TPC pi-",200,10,700,500);

        TH1F *gE_glob_neg = new TH1F(*pi_neg_IT);
        gE_glob_neg->SetName("glob_efficiency_pi_neg");
        gE_glob_neg->Divide(pi_neg_IT, pi_neg_d, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_glob_neg->SetName("Efficienza_globale_pi_neg");
        gE_glob_neg->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{-}");
        gE_glob_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_glob_neg->GetYaxis()->SetTitle("#epsilon");
        
        gE_glob_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gE_glob_neg->GetXaxis()->SetRangeUser(0., 3.); 

        CE_glob_neg->cd();
        gE_glob_neg->Draw("PE");

        efficiencyMomDependent.cd();
        CE_glob_neg->Write();
        CE_glob_neg->Print("efficiency_tagli_mom_dependent/Efficienza globale per traccia ITS+TPC pi-.pdf");

    //------------------- Efficienza di tracciamento TPC SA --------------------//
        TCanvas *CE_tpc_neg = new TCanvas("CE_tpc_neg","Efficienza di tracciamento TPC SA pi-",200,10,700,500);

        TH1F *Copy_pi_neg_IT = new TH1F(*pi_neg_IT);
        Copy_pi_neg_IT ->SetName("copia_pi_neg_IT");
        Copy_pi_neg_IT ->Add(pi_neg_W);
        TH1F *gE_tpc_neg = new TH1F(*Copy_pi_neg_IT);
        gE_tpc_neg->SetName("Efficienza_tpc_pi_neg");
        gE_tpc_neg->Divide(Copy_pi_neg_IT, pi_neg_d, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_neg->SetTitle("Efficienza di tracciamento TPC SA #pi^{-}");
        gE_tpc_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_neg->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg->GetXaxis()->SetRangeUser(0., 3.); 

        CE_tpc_neg->cd();
        gE_tpc_neg->Draw("PE");

        efficiencyMomDependent.cd();
        CE_tpc_neg->Write();
        CE_tpc_neg->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA pi-.pdf");

    //------------------- efficienza matching its-tpc --------------------//
        TCanvas *CE_its_tpc_neg = new TCanvas("CE_its_tpc_neg","Efficienza matching ITS TPC pi-",200,10,700,500);
    
        TH1F *Copy_pi_neg_W = new TH1F(*pi_neg_W);
        Copy_pi_neg_W ->SetName("copia_pi_neg_W");
        Copy_pi_neg_W ->Add(pi_neg_IT);
        TH1F *gE_its_tpc_neg = new TH1F(*Copy_pi_neg_W);
        gE_its_tpc_neg->SetName("Efficienza_its_tpc_pi_neg");
        gE_its_tpc_neg->Divide(pi_neg_IT, Copy_pi_neg_W, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_neg->SetTitle("Efficienza matching ITS TPC #pi^{-}");
        gE_its_tpc_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_neg->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg->GetXaxis()->SetRangeUser(0., 3.); 

        CE_its_tpc_neg->cd();
        gE_its_tpc_neg->Draw("PE");

        efficiencyMomDependent.cd();
        CE_its_tpc_neg->Write();
        CE_its_tpc_neg->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC pi-.pdf");



    //-------------------------------------------- efficienze per pi+ e pi- -------------------------------------//

    //--------------------------------------------------- Efficienza globale per traccia ITS+TPC ------------------------------------------------------------//
        TCanvas *CE_glob = new TCanvas("CE_glob","Efficienza globale per traccia ITS+TPC pi+ e pi-", 400, 600);

        gE_glob_neg->SetLineColor(2);

        gE_glob_pos->GetXaxis()->SetRangeUser(0., 4);
        gE_glob_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gE_glob_neg->GetXaxis()->SetRangeUser(0., 4.); 

        gE_glob_neg->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{+} e #pi^{-}");

        CE_glob->cd();
        gE_glob_neg->Draw("PE");
        gE_glob_pos->Draw("PESAME");

        TLine *l1_glob= new TLine(0.5, 0, 0.5, 1.1);

        l1_glob->SetLineStyle(3);  //anche 7 (e 5) non era male

        l1_glob->Draw();

        TBox* box1_glob= new TBox(0., 0., 0.5, 1.1);
        box1_glob->SetFillColorAlpha(17, 0.50);   //17: colore grigino
        box1_glob->Draw();

        TLatex t1;
        t1.SetNDC();
        t1.SetTextSize(0.031);
        t1.DrawLatex(0.50, 0.83, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza globale per traccia ITS+TPC #font[22]{#pi^{+}}");
        t1.DrawLatex(0.50, 0.78, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza globale per traccia ITS+TPC #font[22]{#pi^{-}}");

        efficiencyMomDependent.cd();
        CE_glob->Write();
        CE_glob->Print("efficiency_tagli_mom_dependent/Efficienza globale per traccia ITS+TPC pi+ e pi-.pdf");

    //--------------------------------------------------- Efficienza di tracciamento TPC SA ------------------------------------------------------------//
        TCanvas *CE_TPC = new TCanvas("CE_TPC","Efficienza di tracciamento TPC SA pi+ e pi-",200,10,700,500);

        gE_tpc_neg->SetLineColor(2);

        gE_tpc_pos->GetXaxis()->SetRangeUser(0., 4);
        gE_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg->GetXaxis()->SetRangeUser(0., 4.); 

        gE_tpc_neg->SetTitle("Efficienza di tracciamento TPC SA #pi^{+} e #pi^{-}");

        CE_TPC->cd();
        gE_tpc_neg->Draw("PE");
        gE_tpc_pos->Draw("PESAME");

        TLine *l1_TPC= new TLine(0.5, 0, 0.5, 1.1);

        l1_TPC->SetLineStyle(3);  //anche 7 (e 5) non era male

        l1_TPC->Draw();

        TBox* box1_TPC= new TBox(0., 0., 0.5, 1.1);
        box1_TPC->SetFillColorAlpha(17, 0.50);   //17: colore grigino
        box1_TPC->Draw();

        TLatex t2;
        t2.SetNDC();
        t2.SetTextSize(0.031);
        t2.DrawLatex(0.50, 0.83, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA #font[22]{#pi^{+}}");
        t2.DrawLatex(0.50, 0.78, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA #font[22]{#pi^{-}}"); 

        efficiencyMomDependent.cd();
        CE_TPC->Write();
        CE_TPC->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA pi+ e pi-.pdf");

    //--------------------------------------------------- Efficienza di matching ITS TPC ------------------------------------------------------------//
        TCanvas *CE_ITS_TPC = new TCanvas("CE_ITS_TPC","Efficienza di matching ITS TPC pi+ e pi-",200,10,700,500);

        gE_its_tpc_neg->SetLineColor(2);

        gE_its_tpc_pos->GetXaxis()->SetRangeUser(0., 4);
        gE_its_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg->GetXaxis()->SetRangeUser(0., 4.); 

        gE_its_tpc_neg->SetTitle("Efficienza di matching ITS TPC #pi^{+} e #pi^{-}");

        CE_ITS_TPC->cd();
        gE_its_tpc_neg->Draw("PE");
        gE_its_tpc_pos->Draw("PESAME");

        TLine *l1_ITS_TPC= new TLine(0.5, 0, 0.5, 1.1);
        l1_ITS_TPC->SetLineStyle(3);  //anche 7 (e 5) non era male
        l1_ITS_TPC->Draw();
        TBox* box1_ITS_TPC= new TBox(0., 0., 0.5, 1.1);
        box1_ITS_TPC->SetFillColorAlpha(17, 0.50);   //17: colore grigino
        box1_ITS_TPC->Draw();

        TLatex t3;
        t3.SetNDC();
        t3.SetTextSize(0.031);
        t3.DrawLatex(0.50, 0.63, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di matching ITS TPC #font[22]{#pi^{+}}");
        t3.DrawLatex(0.50, 0.58, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di matching ITS TPC #font[22]{#pi^{-}}");

        efficiencyMomDependent.cd();
        CE_ITS_TPC->Write();
        CE_ITS_TPC->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC pi+ e pi-.pdf");



    //------------------------------- efficienze con hasTpcSegment -------------------------//
    //------------------------- Calcoli per pi+ -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC --------------------//  questa rimane uguale perché hasTpcSegment è l’equivalente di ITSwithTPC ma fatto con le mc label
  
    //------------------- Efficienza di tracciamento TPC SA --------------------//
        TCanvas *CE_tpc_pos_MClabel = new TCanvas("CE_tpc_pos_MClabel ","Efficienza di tracciamento TPC SA con MClabel  pi+",200,10,700,500);

        TH1F *Copy_pi_pos_IT_MClabel = new TH1F(*pi_pos_IT);
        Copy_pi_pos_IT_MClabel->SetName("copia_pi_pos_IT_MClabel ");
        Copy_pi_pos_IT_MClabel->Add(pi_pos_TpcSegment);

        TH1F *gE_tpc_pos_MClabel = new TH1F(*Copy_pi_pos_IT_MClabel );
        gE_tpc_pos_MClabel->SetName("tpc_efficiency_pi_pos_w/MClabel ");
        gE_tpc_pos_MClabel->Divide(Copy_pi_pos_IT_MClabel, pi_pos_d, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_pos_MClabel->SetTitle("Efficienza di tracciamento TPC SA con MClabel  #pi^{+}");
        gE_tpc_pos_MClabel->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_pos_MClabel->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_pos_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_pos_MClabel->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così
    
        CE_tpc_pos_MClabel->cd();
        gE_tpc_pos_MClabel->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori

        efficiencyMomDependent.cd();
        CE_tpc_pos_MClabel->Write();
        CE_tpc_pos_MClabel->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA con MClabel pi+.pdf");

    //------------------- efficienza matching its-tpc --------------------//
        TCanvas *CE_its_tpc_pos_MClabel = new TCanvas("CE_its_tpc_pos_MClabel ","Efficienza matching ITS TPC con MClabel pi+",200,10,700,500);
        
        TH1F *Copy_pi_pos_MClabel = new TH1F(*pi_pos_TpcSegment);
        Copy_pi_pos_MClabel->SetName("copia_pi_pos_TpcSegment");
        Copy_pi_pos_MClabel->Add(pi_pos_IT);
        TH1F *gE_its_tpc_pos_MClabel = new TH1F(*Copy_pi_pos_MClabel);
        gE_its_tpc_pos_MClabel->SetName("Efficienza_its_tpc_pi_pos_MClabel");
        gE_its_tpc_pos_MClabel->Divide(pi_pos_IT, Copy_pi_pos_MClabel, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_pos_MClabel->SetTitle("Efficienza matching ITS TPC con MClabel #pi^{+}");
        gE_its_tpc_pos_MClabel->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_pos_MClabel->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_pos_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_pos_MClabel->GetXaxis()->SetRangeUser(0., 3.); 

        CE_its_tpc_pos_MClabel->cd();
        gE_its_tpc_pos_MClabel->Draw("PE");

        efficiencyMomDependent.cd();
        CE_its_tpc_pos_MClabel->Write();
        CE_its_tpc_pos_MClabel->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC con MClabel pi+.pdf");


    //------------------------ Calcoli per pi- -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC --------------------//questa rimane uguale

    //------------------- Efficienza di tracciamento TPC SA --------------------//
        TCanvas *CE_tpc_neg_MClabel = new TCanvas("CE_tpc_neg_MClabel","Efficienza di tracciamento TPC SA con MClabel pi-",200,10,700,500);

        TH1F *Copy_pi_neg_IT_MClabel= new TH1F(*pi_neg_IT);
        Copy_pi_neg_IT_MClabel->SetName("copia_pi_neg_IT_MClabel");
        Copy_pi_neg_IT_MClabel->Add(pi_neg_TpcSegment);
        TH1F *gE_tpc_neg_MClabel = new TH1F(*Copy_pi_neg_IT_MClabel);
        gE_tpc_neg_MClabel->SetName("Efficienza_tpc_pi_neg_MClabel");
        gE_tpc_neg_MClabel->Divide(Copy_pi_neg_IT_MClabel, pi_neg_d, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_neg_MClabel->SetTitle("Efficienza di tracciamento TPC SA con MClabel #pi^{-}");
        gE_tpc_neg_MClabel->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_neg_MClabel->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_neg_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg_MClabel->GetXaxis()->SetRangeUser(0., 3.); 

        CE_tpc_neg_MClabel->cd();
        gE_tpc_neg_MClabel->Draw("PE");

        efficiencyMomDependent.cd();
        CE_tpc_neg_MClabel->Write();
        CE_tpc_neg_MClabel->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA con MClabel pi-.pdf");




        

    //------------------- efficienza matching its-tpc --------------------//
        TCanvas *CE_its_tpc_neg_MClabel = new TCanvas("CE_its_tpc_neg_MClabel","Efficienza matching ITS TPC con MClabel pi-",200,10,700,500);
    
        TH1F *Copy_pi_neg_MClabel = new TH1F(*pi_neg_TpcSegment);
        Copy_pi_neg_MClabel->SetName("copia_pi_neg_MClabel");
        Copy_pi_neg_MClabel->Add(pi_neg_IT);
        TH1F *gE_its_tpc_neg_MClabel = new TH1F(*Copy_pi_neg_MClabel);
        gE_its_tpc_neg_MClabel->SetName("Efficienza_its_tpc_pi_neg_MClabel");
        gE_its_tpc_neg_MClabel->Divide(pi_neg_IT, Copy_pi_neg_MClabel, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_neg_MClabel->SetTitle("Efficienza matching ITS TPC con MClabel #pi^{-}");
        gE_its_tpc_neg_MClabel->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_neg_MClabel->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_neg_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg_MClabel->GetXaxis()->SetRangeUser(0., 3.); 

        CE_its_tpc_neg_MClabel->cd();
        gE_its_tpc_neg_MClabel->Draw("PE");

        efficiencyMomDependent.cd();
        CE_its_tpc_neg_MClabel->Write();
        CE_its_tpc_neg_MClabel->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC con MClabel pi-.pdf");




    //----------------------------- confronto tra efficienze ottenute con rematching e quelle ottenute con Mclabel -----------------------------//

    //--------------------------------- Confronto tra efficienze TPC per p+ ---------------------------------//
        TCanvas *Ccomparison_tpc_pi_pos= new TCanvas("Ccomparison_tpc_pi_pos","Confronto Efficienza di tracciamento TPC SA pi+",200,10,700,500);

        gE_tpc_pos->SetLineColor(2);

        gE_tpc_pos->GetXaxis()->SetRangeUser(0., 3);
        gE_tpc_pos_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_pos_MClabel->GetXaxis()->SetRangeUser(0., 3);

        gE_tpc_pos->SetTitle("Confronto Efficienza di tracciamento TPC SA #pi^{+}");

        Ccomparison_tpc_pi_pos->cd();
        gE_tpc_pos->Draw("PE");
        gE_tpc_pos_MClabel->Draw("PESAME");
        
        TLatex t4;
        t4.SetNDC();
        t4.SetTextSize(0.031);
        t4.DrawLatex(0.5, 0.45, "Confronto tra efficienze TPC per #font[22]{#pi^{+}}");  //un po’ centrato
        t4.DrawLatex(0.5, 0.40, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con MC label");
        t4.DrawLatex(0.5, 0.35, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA ottenuta da rematching");

        efficiencyMomDependent.cd();
        Ccomparison_tpc_pi_pos->Write();
        Ccomparison_tpc_pi_pos->Print("efficiency_tagli_mom_dependent/confronto Efficienza di tracciamento TPC SA pi+.pdf");


    //--------------------------------- Confronto tra efficienze TPC per p- ---------------------------------//
        TCanvas *Ccomparison_tpc_pi_neg= new TCanvas("Ccomparison_tpc_pi_neg","Confronto Efficienza di tracciamento TPC SA pi-",200,10,700,500);

        gE_tpc_neg->SetLineColor(2);

        gE_tpc_neg->GetXaxis()->SetRangeUser(0., 3);
        gE_tpc_neg_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg_MClabel->GetXaxis()->SetRangeUser(0., 3);

        gE_tpc_neg->SetTitle("Confronto Efficienza di tracciamento TPC SA #pi^{-}");

        Ccomparison_tpc_pi_neg->cd();
        gE_tpc_neg->Draw("PE");
        gE_tpc_neg_MClabel->Draw("PESAME");
        
        TLatex t5;
        t5.SetNDC();
        t5.SetTextSize(0.031);
        t5.DrawLatex(0.5, 0.45, "Confronto tra efficienze TPC per #font[22]{#pi^{-}}");  //un po’ centrato
        t5.DrawLatex(0.5, 0.40, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con MC label");
        t5.DrawLatex(0.5, 0.35, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA ottenuta da rematching");

        efficiencyMomDependent.cd();
        Ccomparison_tpc_pi_neg->Write();
        Ccomparison_tpc_pi_neg->Print("efficiency_tagli_mom_dependent/confronto Efficienza di tracciamento TPC SA pi-.pdf");


    //---------------------------------- Confronto tra efficienze ITS-TPC per p+ ---------------------------------//
        TCanvas *Ccomparison_its_tpc_pi_pos= new TCanvas("Ccomparison_its_tpc_pi_pos","Confronto efficienza ITS-TPC pi+",200,10,700,500);

        gE_its_tpc_pos->SetLineColor(2);

        gE_its_tpc_pos->GetXaxis()->SetRangeUser(0., 3);
        gE_its_tpc_pos_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_pos_MClabel->GetXaxis()->SetRangeUser(0., 3);

        gE_its_tpc_pos->SetTitle("Confronto efficienza ITS-TPC #pi^{+}");

        Ccomparison_its_tpc_pi_pos->cd();
        gE_its_tpc_pos->Draw("PE");
        gE_its_tpc_pos_MClabel->Draw("PESAME");
        
        TLatex t6;
        t6.SetNDC();
        t6.SetTextSize(0.031);
        t6.DrawLatex(0.5, 0.48, "Confronto tra efficienze ITS-TPC per #font[22]{#pi^{+}}");  //un po’ centrato
        t6.DrawLatex(0.5, 0.43, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con MC label");
        t6.DrawLatex(0.5, 0.38, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA ottenuta da rematching");

        efficiencyMomDependent.cd();
        Ccomparison_its_tpc_pi_pos->Write();
        Ccomparison_its_tpc_pi_pos->Print("efficiency_tagli_mom_dependent/confronto efficienza its-tpc pi+.pdf");


    //--------------------------------- Confronto tra efficienze ITS-TPC per p- ---------------------------------//
        TCanvas *Ccomparison_its_tpc_pi_neg= new TCanvas("Ccomparison_its_tpc_pi_neg","Confronto efficienza ITS-TPC pi-",200,10,700,500);

        gE_its_tpc_neg->SetLineColor(2);

        gE_its_tpc_neg->GetXaxis()->SetRangeUser(0., 3);
        gE_its_tpc_neg_MClabel->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg_MClabel->GetXaxis()->SetRangeUser(0., 3);

        gE_its_tpc_neg->SetTitle("Confronto efficienza ITS-TPC #pi^{-}");

        Ccomparison_its_tpc_pi_neg->cd();
        gE_its_tpc_neg->Draw("PE");
        gE_its_tpc_neg_MClabel->Draw("PESAME");
        
        TLatex t7;
        t7.SetNDC();
        t7.SetTextSize(0.031);
        t7.DrawLatex(0.5, 0.48, "Confronto tra efficienze ITS-TPC per #font[22]{#pi^{-}}");  //un po’ centrato
        t7.DrawLatex(0.5, 0.43, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con MC label");
        t7.DrawLatex(0.5, 0.38, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA ottenuta da rematching");

        efficiencyMomDependent.cd();
        Ccomparison_its_tpc_pi_neg->Write();
        Ccomparison_its_tpc_pi_neg->Print("efficiency_tagli_mom_dependent/confronto efficienza its-tpc pi-.pdf");



    //------------------------- CALCOLO DELLE EFFICIENZE CON ISTOGRAMMI REBINNATI -------------------------//
    //------------------------- Calcoli per pi+ rebinned -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC rebinned --------------------//
        TCanvas *CE_glob_pos_rebin = new TCanvas("CE_glob_pos_rebin","Efficienza globale per traccia ITS+TPC rebinned pi+", 400, 600);

        TH1F *gE_glob_pos_rebin = new TH1F(*pi_pos_IT_rebin);
        gE_glob_pos_rebin->SetName("efficiency_pi_pos_glob_rebin");
        gE_glob_pos_rebin->Divide(pi_pos_IT_rebin, pi_pos_d_rebin, 1., 1., "B");
        gE_glob_pos_rebin->SetTitle("Efficienza globale per traccia ITS+TPC rebinned #pi^{+}");
        gE_glob_pos_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_glob_pos_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_glob_pos_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_glob_pos_rebin->GetXaxis()->SetRangeUser(0., 4.);

        CE_glob_pos_rebin->cd();
        gE_glob_pos_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_glob_pos_rebin->Write();
        CE_glob_pos_rebin->Print("efficiency_tagli_mom_dependent/Efficienza globale per traccia ITS+TPC rebinned pi+.pdf");

    //------------------- Efficienza di tracciamento TPC SA rebinned --------------------//
        TCanvas *CE_tpc_pos_rebin = new TCanvas("CE_tpc_pos_rebin","Efficienza di tracciamento TPC SA rebinned pi+",200,10,700,500);

        TH1F *Copy_pi_pos_IT_rebin= new TH1F(*pi_pos_IT_rebin);
        Copy_pi_pos_IT_rebin->SetName("copia_pi_pos_IT_rebin");
        Copy_pi_pos_IT_rebin->Add(pi_pos_W_rebin);

        TH1F *gE_tpc_pos_rebin= new TH1F(*Copy_pi_pos_IT_rebin);
        gE_tpc_pos_rebin->SetName("tpc_efficiency_pi_pos_rebin");
        gE_tpc_pos_rebin->Divide(Copy_pi_pos_IT_rebin, pi_pos_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_pos_rebin->SetTitle("Efficienza di tracciamento TPC SA rebinned #pi^{+}");
        gE_tpc_pos_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_pos_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_pos_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_pos_rebin->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così

        CE_tpc_pos_rebin->cd();
        gE_tpc_pos_rebin->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_tpc_pos_rebin->Write();
        CE_tpc_pos_rebin->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA rebinned pi+.pdf");

    //------------------- efficienza matching its-tpc rebinned --------------------//
        TCanvas *CE_its_tpc_pos_rebin = new TCanvas("CE_its_tpc_pos_rebin","Efficienza matching ITS TPC rebinned pi+",200,10,700,500);
    
        TH1F *Copy_pi_pos_W_rebin= new TH1F(*pi_pos_W_rebin);
        Copy_pi_pos_W_rebin->SetName("copia_pi_pos_W_rebin");
        Copy_pi_pos_W_rebin->Add(pi_pos_IT_rebin);

        TH1F *gE_its_tpc_pos_rebin= new TH1F(*Copy_pi_pos_W_rebin);
        gE_its_tpc_pos_rebin->SetName("its_tpc_efficiency_pi_pos_rebin");
        gE_its_tpc_pos_rebin->Divide(pi_pos_IT_rebin, Copy_pi_pos_W_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_pos_rebin->SetName("Efficienza_matching_ITS_TPC_pi_pos_rebin");
        gE_its_tpc_pos_rebin->SetTitle("Efficienza matching ITS TPC rebinned #pi^{+}");
        gE_its_tpc_pos_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_pos_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_pos_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_pos_rebin->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così
    
        CE_its_tpc_pos_rebin->cd();
        gE_its_tpc_pos_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_its_tpc_pos_rebin->Write();
        CE_its_tpc_pos_rebin->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC rebinned pi+.pdf");



    //------------------------- Calcoli per pi- rebinned -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC rebinned --------------------//
        TCanvas *CE_glob_neg_rebin = new TCanvas("CE_glob_neg_rebin","Efficienza globale per traccia ITS+TPC rebinned pi-",200,10,700,500);

        TH1F *gE_glob_neg_rebin = new TH1F(*pi_neg_IT_rebin);
        gE_glob_neg_rebin->SetName("glob_efficiency_pi_neg_rebin");
        gE_glob_neg_rebin->Divide(pi_neg_IT_rebin, pi_neg_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_glob_neg_rebin->SetName("Efficienza_globale_pi_neg_rebin");
        gE_glob_neg_rebin->SetTitle("Efficienza globale per traccia ITS+TPC rebinned #pi^{-}");
        gE_glob_neg_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_glob_neg_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_glob_neg_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_glob_neg_rebin->GetXaxis()->SetRangeUser(0., 3.); 

        CE_glob_neg_rebin->cd();
        gE_glob_neg_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_glob_neg_rebin->Write();
        CE_glob_neg_rebin->Print("efficiency_tagli_mom_dependent/Efficienza globale per traccia ITS+TPC rebinned pi-.pdf");

    //------------------- Efficienza di tracciamento TPC SA rebinned --------------------//
        TCanvas *CE_tpc_neg_rebin = new TCanvas("CE_tpc_neg_rebin","Efficienza di tracciamento TPC SA rebinned pi-",200,10,700,500);

        TH1F *Copy_pi_neg_IT_rebin = new TH1F(*pi_neg_IT_rebin);
        Copy_pi_neg_IT_rebin ->SetName("copia_pi_neg_IT_rebin");
        Copy_pi_neg_IT_rebin ->Add(pi_neg_W_rebin);
        TH1F *gE_tpc_neg_rebin = new TH1F(*Copy_pi_neg_IT_rebin);
        gE_tpc_neg_rebin->SetName("Efficienza_tpc_pi_neg_rebin");
        gE_tpc_neg_rebin->Divide(Copy_pi_neg_IT_rebin, pi_neg_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_neg_rebin->SetTitle("Efficienza di tracciamento TPC SA rebinned #pi^{-}");
        gE_tpc_neg_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_neg_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_neg_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg_rebin->GetXaxis()->SetRangeUser(0., 3.); 

        CE_tpc_neg_rebin->cd();
        gE_tpc_neg_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_tpc_neg_rebin->Write();
        CE_tpc_neg_rebin->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA rebinned pi-.pdf");

    //------------------- efficienza matching its-tpc rebinned --------------------//
        TCanvas *CE_its_tpc_neg_rebin = new TCanvas("CE_its_tpc_neg_rebin","Efficienza matching ITS TPC rebinned pi-",200,10,700,500);

        TH1F *Copy_pi_neg_W_rebin = new TH1F(*pi_neg_W_rebin);
        Copy_pi_neg_W_rebin ->SetName("copia_pi_neg_W_rebin");
        Copy_pi_neg_W_rebin ->Add(pi_neg_IT_rebin);
        TH1F *gE_its_tpc_neg_rebin = new TH1F(*Copy_pi_neg_W_rebin);
        gE_its_tpc_neg_rebin->SetName("Efficienza_its_tpc_pi_neg_rebin");
        gE_its_tpc_neg_rebin->Divide(pi_neg_IT_rebin, Copy_pi_neg_W_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_neg_rebin->SetTitle("Efficienza matching ITS TPC rebinned #pi^{-}");
        gE_its_tpc_neg_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_neg_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_neg_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg_rebin->GetXaxis()->SetRangeUser(0., 3.); 

        CE_its_tpc_neg_rebin->cd();
        gE_its_tpc_neg_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_its_tpc_neg_rebin->Write();
        CE_its_tpc_neg_rebin->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC rebinned pi-.pdf");




    //-------------------------------------------- efficienze per pi+ e pi- rebinned -------------------------------------//

    //--------------------------------------------------- Efficienza globale per traccia ITS+TPC rebinned ------------------------------------------------------------//
        TCanvas *CE_glob_rebin = new TCanvas("CE_glob_rebin","Efficienza globale per traccia ITS+TPC rebinned pi+ e pi-", 600, 600);

        gE_glob_neg_rebin->SetLineColor(2);

        gE_glob_pos_rebin->GetXaxis()->SetRangeUser(0., 4);
        gE_glob_neg_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_glob_neg_rebin->GetXaxis()->SetRangeUser(0., 4.); 

        gE_glob_neg_rebin->SetTitle("Efficienza globale ITS+TPC per pioni carichi");

        CE_glob_rebin->cd();
        gE_glob_neg_rebin->Draw("PE");
        gE_glob_pos_rebin->Draw("PESAME");

        TLine *l1_glob_pos_rebin= new TLine(0.5, 0, 0.5, 1.1);
        l1_glob_pos_rebin->SetLineStyle(3);  //anche 7 (e 5) non era male
        l1_glob_pos_rebin->Draw();
        TBox* box1_glob_pos_rebin= new TBox(0., 0., 0.5, 1.1);
        box1_glob_pos_rebin->SetFillColorAlpha(17, 0.50);   //17: colore grigino
        box1_glob_pos_rebin->Draw();
        
        TLatex t1_rebin;
        t1_rebin.SetNDC();
        t1_rebin.SetTextSize(0.031);
        t1_rebin.DrawLatex(0.3, 0.83, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza globale per traccia ITS+TPC #font[22]{#pi^{+}}");
        t1_rebin.DrawLatex(0.3, 0.78, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza globale per traccia ITS+TPC #font[22]{#pi^{-}}");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_glob_rebin->Write();
        CE_glob_rebin->Print("efficiency_tagli_mom_dependent/Efficienza globale per traccia ITS+TPC rebinned pi+ e pi-.pdf");

    //--------------------------------------------------- Efficienza di tracciamento TPC SA rebinned ------------------------------------------------------------//
        TCanvas *CE_TPC_rebin = new TCanvas("CE_TPC_rebin","Efficienza di tracciamento TPC SA rebinned pi+ e pi-", 600, 600);

        gE_tpc_neg_rebin->SetLineColor(2);

        gE_tpc_pos_rebin->GetXaxis()->SetRangeUser(0., 4);
        gE_tpc_neg_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg_rebin->GetXaxis()->SetRangeUser(0., 4.); 

        gE_tpc_neg_rebin->SetTitle("Efficienza di tracciamento TPC SA per #pi^{+} e #pi^{-}");

        CE_TPC_rebin->cd();
        gE_tpc_neg_rebin->Draw("PE");
        gE_tpc_pos_rebin->Draw("PESAME");

        TLine *l1_TPC_rebin= new TLine(0.6, 0, 0.6, 1.1);
        l1_TPC_rebin->SetLineStyle(3);  //anche 7 (e 5) non era male
        l1_TPC_rebin->Draw();
        TBox* box1_TPC_rebin= new TBox(0., 0., 0.6, 1.1);
        box1_TPC_rebin->SetFillColorAlpha(17, 0.50);   //17: colore grigino
        box1_TPC_rebin->Draw();

        TLatex t2_rebin;
        t2_rebin.SetNDC();
        t2_rebin.SetTextSize(0.031);
        t2_rebin.DrawLatex(0.3, 0.83, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA #font[22]{#pi^{+}}");
        t2_rebin.DrawLatex(0.3, 0.78, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA #font[22]{#pi^{-}}"); 

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_TPC_rebin->Write();
        CE_TPC_rebin->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA rebinned pi+ e pi-.pdf");

    //--------------------------------------------------- Efficienza di matching ITS TPC rebinned ------------------------------------------------------------//
        TCanvas *CE_ITS_TPC_rebin = new TCanvas("CE_ITS_TPC_rebin","Efficienza di matching ITS TPC rebinned pi+ e pi-", 600, 600);

        gE_its_tpc_neg_rebin->SetLineColor(2);

        gE_its_tpc_pos_rebin->GetXaxis()->SetRangeUser(0., 4);
        gE_its_tpc_neg_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg_rebin->GetXaxis()->SetRangeUser(0., 4.); 

        gE_its_tpc_neg_rebin->SetTitle("Efficienza di matching ITS TPC #pi^{+} e #pi^{-}");

        CE_ITS_TPC_rebin->cd();
        gE_its_tpc_neg_rebin->Draw("PE");
        gE_its_tpc_pos_rebin->Draw("PESAME");

        TLine *l1_ITS_TPC_rebin= new TLine(0.6, 0, 0.6, 1.1);
        l1_ITS_TPC_rebin->SetLineStyle(3);  //anche 7 (e 5) non era male
        l1_ITS_TPC_rebin->Draw();
        TBox* box1_ITS_TPC_rebin= new TBox(0., 0., 0.6, 1.1);
        box1_ITS_TPC_rebin->SetFillColorAlpha(17, 0.50);   //17: colore grigino
        box1_ITS_TPC_rebin->Draw();

        TLatex t3_rebin;
        t3_rebin.SetNDC();
        t3_rebin.SetTextSize(0.031);
        t3_rebin.DrawLatex(0.3, 0.63, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di matching ITS TPC #font[22]{#pi^{+}}");
        t3_rebin.DrawLatex(0.3, 0.58, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di matching ITS TPC #font[22]{#pi^{-}}");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_ITS_TPC_rebin->Write();
        CE_ITS_TPC_rebin->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC rebinned pi+ e pi-.pdf");




    //------------------------------- efficienze con hasTpcSegment rebinned -------------------------//
    //------------------------- Calcoli per pi+ rebinned -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC rebinned --------------------//  questa rimane uguale perché hasTpcSegment rebinned è l’equivalente di ITSwithTPC ma fatto con le mc label
    
    //------------------- Efficienza di tracciamento TPC SA rebinned --------------------//
        TCanvas *CE_tpc_pos_MClabel_rebin = new TCanvas("CE_tpc_pos_MClabel_rebin ","Efficienza di tracciamento TPC SA con MClabel rebinned pi+",200,10,700,500);

        TH1F *Copy_pi_pos_IT_MClabel_rebin = new TH1F(*pi_pos_IT_rebin);
        Copy_pi_pos_IT_MClabel_rebin->SetName("copia_pi_pos_IT_MClabel_rebin");
        Copy_pi_pos_IT_MClabel_rebin->Add(pi_pos_TpcSegment_rebin);

        TH1F *gE_tpc_pos_MClabel_rebin = new TH1F(*Copy_pi_pos_IT_MClabel_rebin );
        gE_tpc_pos_MClabel_rebin->SetName("tpc_efficiency_pi_pos_w/MClabel_rebin");
        gE_tpc_pos_MClabel_rebin->Divide(Copy_pi_pos_IT_MClabel_rebin, pi_pos_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_pos_MClabel_rebin->SetTitle("Efficienza di tracciamento TPC SA con MClabel rebinned  #pi^{+}");
        gE_tpc_pos_MClabel_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_pos_MClabel_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_pos_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_pos_MClabel_rebin->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così

        CE_tpc_pos_MClabel_rebin->cd();
        gE_tpc_pos_MClabel_rebin->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_tpc_pos_MClabel_rebin->Write();
        CE_tpc_pos_MClabel->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA con MClabel rebinned pi+.pdf");

    //------------------- efficienza matching its-tpc rebinned --------------------//
        TCanvas *CE_its_tpc_pos_MClabel_rebin = new TCanvas("CE_its_tpc_pos_MClabel_rebin ","Efficienza matching ITS TPC con MClabel rebinned pi+",200,10,700,500);

        TH1F *Copy_pi_pos_MClabel_rebin = new TH1F(*pi_pos_TpcSegment_rebin);
        Copy_pi_pos_MClabel_rebin->SetName("copia_pi_pos_TpcSegment_rebin");
        Copy_pi_pos_MClabel_rebin->Add(pi_pos_IT_rebin);
        TH1F *gE_its_tpc_pos_MClabel_rebin = new TH1F(*Copy_pi_pos_MClabel_rebin);
        gE_its_tpc_pos_MClabel_rebin->SetName("Efficienza_its_tpc_pi_pos_MClabel_rebin");
        gE_its_tpc_pos_MClabel_rebin->Divide(pi_pos_IT_rebin, Copy_pi_pos_MClabel_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_pos_MClabel_rebin->SetTitle("Efficienza matching ITS TPC con MClabel rebinned #pi^{+}");
        gE_its_tpc_pos_MClabel_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_pos_MClabel_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_pos_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_pos_MClabel_rebin->GetXaxis()->SetRangeUser(0., 3.); 

        CE_its_tpc_pos_MClabel_rebin->cd();
        gE_its_tpc_pos_MClabel_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_its_tpc_pos_MClabel_rebin->Write();
        CE_its_tpc_pos_MClabel_rebin->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC con MClabel rebinned pi+.pdf");


    //------------------------- Calcoli per pi- rebinned -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC rebinned --------------------//questa rimane uguale

    //------------------- Efficienza di tracciamento TPC SA rebinned --------------------//
        TCanvas *CE_tpc_neg_MClabel_rebin = new TCanvas("CE_tpc_neg_MClabel_rebin","Efficienza di tracciamento TPC SA con MClabel rebinned pi-",200,10,700,500);

        TH1F *Copy_pi_neg_IT_MClabel_rebin= new TH1F(*pi_neg_IT_rebin);
        Copy_pi_neg_IT_MClabel_rebin->SetName("copia_pi_neg_IT_MClabel_rebin");
        Copy_pi_neg_IT_MClabel_rebin->Add(pi_neg_TpcSegment_rebin);
        TH1F *gE_tpc_neg_MClabel_rebin = new TH1F(*Copy_pi_neg_IT_MClabel_rebin);
        gE_tpc_neg_MClabel_rebin->SetName("Efficienza_tpc_pi_neg_MClabel_rebin");
        gE_tpc_neg_MClabel_rebin->Divide(Copy_pi_neg_IT_MClabel_rebin, pi_neg_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_tpc_neg_MClabel_rebin->SetTitle("Efficienza di tracciamento TPC SA con MClabel rebinned #pi^{-}");
        gE_tpc_neg_MClabel_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_tpc_neg_MClabel_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_tpc_neg_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg_MClabel_rebin->GetXaxis()->SetRangeUser(0., 3.); 

        CE_tpc_neg_MClabel_rebin->cd();
        gE_tpc_neg_MClabel_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_tpc_neg_MClabel_rebin->Write();
        CE_tpc_neg_MClabel_rebin->Print("efficiency_tagli_mom_dependent/Efficienza di tracciamento TPC SA con MClabel rebinned pi-.pdf");

    //------------------- efficienza matching its-tpc rebinned --------------------//
        TCanvas *CE_its_tpc_neg_MClabel_rebin = new TCanvas("CE_its_tpc_neg_MClabel_rebin","Efficienza matching ITS TPC con MClabel rebinned pi-",200,10,700,500);

        TH1F *Copy_pi_neg_MClabel_rebin = new TH1F(*pi_neg_TpcSegment_rebin);
        Copy_pi_neg_MClabel_rebin->SetName("copia_pi_neg_MClabel_rebin");
        Copy_pi_neg_MClabel_rebin->Add(pi_neg_IT_rebin);
        TH1F *gE_its_tpc_neg_MClabel_rebin = new TH1F(*Copy_pi_neg_MClabel_rebin);
        gE_its_tpc_neg_MClabel_rebin->SetName("Efficienza_its_tpc_pi_neg_MClabel_rebin");
        gE_its_tpc_neg_MClabel_rebin->Divide(pi_neg_IT_rebin, Copy_pi_neg_MClabel_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gE_its_tpc_neg_MClabel_rebin->SetTitle("Efficienza matching ITS TPC con MClabel rebinned #pi^{-}");
        gE_its_tpc_neg_MClabel_rebin->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gE_its_tpc_neg_MClabel_rebin->GetYaxis()->SetTitle("#epsilon");
        
        gE_its_tpc_neg_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg_MClabel_rebin->GetXaxis()->SetRangeUser(0., 3.); 

        CE_its_tpc_neg_MClabel_rebin->cd();
        gE_its_tpc_neg_MClabel_rebin->Draw("PE");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        CE_its_tpc_neg_MClabel_rebin->Write();
        CE_its_tpc_neg_MClabel_rebin->Print("efficiency_tagli_mom_dependent/efficienza matching ITS TPC con MClabel rebinned pi-.pdf");



    //----------------------------- confronto tra efficienze ottenute con rematching e quelle ottenute con Mclabel rebinned -----------------------------//

    //--------------------------------- Confronto tra efficienze TPC rebinned per p+ ---------------------------------//
        TCanvas *Ccomparison_tpc_pi_pos_rebin= new TCanvas("Ccomparison_tpc_pi_pos_rebin","Confronto Efficienza di tracciamento TPC SA rebinned pi+",200,10,700,500);

        gE_tpc_pos_rebin->SetLineColor(2);

        gE_tpc_pos_rebin->GetXaxis()->SetRangeUser(0., 4);
        gE_tpc_pos_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_pos_MClabel_rebin->GetXaxis()->SetRangeUser(0., 4);

        gE_tpc_pos_rebin->SetTitle("");

        Ccomparison_tpc_pi_pos_rebin->cd();
        gE_tpc_pos_rebin->Draw("PE");
        gE_tpc_pos_MClabel_rebin->Draw("PESAME");

        TLatex t4_rebin;
        t4_rebin.SetNDC();
        t4_rebin.SetTextSize(0.031);
        t4_rebin.DrawLatex(0.25, 0.45, "Confronto tra efficienze TPC per #font[22]{#pi^{+}}");  //un po’ centrato
        t4_rebin.DrawLatex(0.25, 0.40, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con verit#grave{a} MC");
        t4_rebin.DrawLatex(0.25, 0.35, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con tag-and-probe+rematching");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        Ccomparison_tpc_pi_pos_rebin->Write();
        Ccomparison_tpc_pi_pos_rebin->Print("efficiency_tagli_mom_dependent/confronto Efficienza di tracciamento TPC SA rebinned pi+.pdf");





    //--------------------------------- Confronto tra efficienze TPC rebinned per p- ---------------------------------//
        TCanvas *Ccomparison_tpc_pi_neg_rebin= new TCanvas("Ccomparison_tpc_pi_neg_rebin","Confronto Efficienza di tracciamento TPC SA rebinned pi-",200,10,700,500);

        gE_tpc_neg_rebin->SetLineColor(2);

        gE_tpc_neg_rebin->GetXaxis()->SetRangeUser(0., 3);
        gE_tpc_neg_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_tpc_neg_MClabel_rebin->GetXaxis()->SetRangeUser(0., 3);

        gE_tpc_neg_rebin->SetTitle("Confronto Efficienza di tracciamento TPC SA rebinned #pi^{-}");

        Ccomparison_tpc_pi_neg_rebin->cd();
        gE_tpc_neg_rebin->Draw("PE");
        gE_tpc_neg_MClabel_rebin->Draw("PESAME");
        
        TLatex t5_rebin;
        t5_rebin.SetNDC();
        t5_rebin.SetTextSize(0.031);
        t5_rebin.DrawLatex(0.3, 0.45, "Confronto tra efficienze TPC per #font[22]{#pi^{-}}");  //un po’ centrato
        t5_rebin.DrawLatex(0.3, 0.40, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA con MC label");
        t5_rebin.DrawLatex(0.3, 0.35, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di tracciamento TPC SA ottenuta da rematching");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        Ccomparison_tpc_pi_neg_rebin->Write();
        Ccomparison_tpc_pi_neg_rebin->Print("efficiency_tagli_mom_dependent/confronto Efficienza di tracciamento TPC SA rebinned pi-.pdf");




    //--------------------------------- Confronto tra efficienze ITS-TPC rebinned per p+ ---------------------------------//
        TCanvas *Ccomparison_its_tpc_pi_pos_rebin= new TCanvas("Ccomparison_its_tpc_pi_pos_rebin","Confronto efficienza ITS-TPC rebinned pi+",200,10,700,500);

        gE_its_tpc_pos_rebin->SetLineColor(2);

        gE_its_tpc_pos_rebin->GetXaxis()->SetRangeUser(0., 4);
        gE_its_tpc_pos_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_pos_MClabel_rebin->GetXaxis()->SetRangeUser(0., 4);

        gE_its_tpc_pos_rebin->SetTitle("");

        Ccomparison_its_tpc_pi_pos_rebin->cd();
        gE_its_tpc_pos_rebin->Draw("PE");
        gE_its_tpc_pos_MClabel_rebin->Draw("PESAME");
        
        TLatex t6_rebin;
        t6_rebin.SetNDC();
        t6_rebin.SetTextSize(0.031);
        t6_rebin.DrawLatex(0.25, 0.48, "Confronto tra efficienze di matching ITS-TPC per #font[22]{#pi^{+}}");  //un po’ centrato
        t6_rebin.DrawLatex(0.25, 0.43, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di matching ITS-TPC con verit#grave{a} MC");
        t6_rebin.DrawLatex(0.25, 0.38, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di matching ITS-TPC con tag-and-probe+rematching");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        Ccomparison_its_tpc_pi_pos_rebin->Write();
        Ccomparison_its_tpc_pi_pos_rebin->Print("efficiency_tagli_mom_dependent/confronto efficienza its-tpc rebinned pi+.pdf");


    //--------------------------------- Confronto tra efficienze ITS-TPC rebinned per p- ---------------------------------//
        TCanvas *Ccomparison_its_tpc_pi_neg_rebin= new TCanvas("Ccomparison_its_tpc_pi_neg_rebin","Confronto efficienza ITS-TPC rebinned pi-",200,10,700,500);

        gE_its_tpc_neg_rebin->SetLineColor(2);

        gE_its_tpc_neg_rebin->GetXaxis()->SetRangeUser(0., 3);
        gE_its_tpc_neg_MClabel_rebin->GetYaxis()->SetRangeUser(0., 1.1);
        gE_its_tpc_neg_MClabel_rebin->GetXaxis()->SetRangeUser(0., 3);

        gE_its_tpc_neg_rebin->SetTitle("Confronto efficienza ITS-TPC rebinned #pi^{-}");

        Ccomparison_its_tpc_pi_neg_rebin->cd();
        gE_its_tpc_neg_rebin->Draw("PE");
        gE_its_tpc_neg_MClabel_rebin->Draw("PESAME");

        TLatex t7_rebin;
        t7_rebin.SetNDC();
        t7_rebin.SetTextSize(0.031);
        t7_rebin.DrawLatex(0.3, 0.48, "Confronto tra efficienze di matching ITS-TPC per #font[22]{#pi^{-}}");  //un po’ centrato
        t7_rebin.DrawLatex(0.3, 0.43, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza di matching ITS-TPC con MC label");
        t7_rebin.DrawLatex(0.3, 0.38, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza di matching ITS-TPC ottenuta da rematching");

        efficiencyMomDependent.cd("eff_with_rebinned_histos");
        Ccomparison_its_tpc_pi_neg_rebin->Write();
        Ccomparison_its_tpc_pi_neg_rebin->Print("efficiency_tagli_mom_dependent/confronto efficienza its-tpc rebinned pi-.pdf");


        efficiencyMomDependent.Close();

}