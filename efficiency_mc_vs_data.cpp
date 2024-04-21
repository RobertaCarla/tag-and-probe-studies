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




void efficiency_mc_vs_data(){
    TFile mcFile("AnalysisResults_mc_20240328.root");
    TFile dataFile("AnalysisResults_data_20240328.root");
    TFile efficiency_mc_vs_data("efficiency_mc_vs_data.root", "RECREATE");

    TH3F *mc_piRec= (TH3F*)mcFile.Get("efficiency-q-a/piRec");
    TH3F *data_piRec = (TH3F*)dataFile.Get("efficiency-q-a/piRec");

    TH1F *mcDecays= (TH1F*)mc_piRec->ProjectionY("Decays",1,1,1,900);
    mcDecays->SetName("mcDecays");
    TH1F *mcITS_TPC= (TH1F*)mc_piRec->ProjectionY("ITS_TPC",6,6,1,900);
    mcITS_TPC->SetName("mcITS_TPC");
    TH1F *mcITSwithTPC= (TH1F*)mc_piRec->ProjectionY("ITS with TPC",8,8,1,900);
    mcITSwithTPC->SetName("mcITSwithTPC");

    TH1F *dataDecays= (TH1F*)data_piRec->ProjectionY("Decays",1,1,1,900);
    dataDecays->SetName("dataDecays");
    TH1F *dataITS_TPC= (TH1F*)data_piRec->ProjectionY("ITS_TPC",6,6,1,900);
    dataITS_TPC->SetName("dataITS_TPC");
    TH1F *dataITSwithTPC= (TH1F*)data_piRec->ProjectionY("ITS with TPC",8,8,1,900);
    dataITSwithTPC->SetName("dataITSwithTPC");

    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);

    efficiency_mc_vs_data.cd();
    mcDecays->Write();
    dataDecays->Write();
    mcITS_TPC->Write();
    dataITS_TPC->Write();

    //--------------------------------------------------------- decays ---------------------------------------------------------//
        int nbin= mcDecays->GetNbinsX();
        float min;
        int n = dimension(mcDecays, nbin, min);

    //------------------- decays pi+ -----------------//
        float massimox= mcDecays->GetXaxis()->GetBinUpEdge(nbin);
        TH1F *mc_pi_pos_d= new TH1F("mc_pi_pos_d", "pi+ decays;#it{p}_{T} (GeV/#it{c}); Entries",n,0.,massimox);
        mc_pi_pos_d->SetName("mc_decays_pi_pos");
   
        mc_pi_pos_d->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        mc_pi_pos_d->SetBinContent(i, mcDecays->GetBinContent(min-1+i));
        }

    //------------------- decays pi+ rebinned -----------------//
        TH1F* mc_pi_pos_d_rebin= rebin(mc_pi_pos_d);
        mc_pi_pos_d_rebin->SetName("mc_pi_pos_d_rebin");
        mc_pi_pos_d_rebin->SetTitle("mc_pi_pos_d rebinned");
        mc_pi_pos_d_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        mc_pi_pos_d_rebin->GetYaxis()->SetTitle("Entries");

    
    //------------- decays pi- -----------------------//
        TH1F *mc_pi_neg_d= new TH1F("mc_pi_neg_d", "pi- decays;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        mc_pi_neg_d->SetName("mc_decays_pi_neg");

        mc_pi_neg_d->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        mc_pi_neg_d->SetBinContent(i, mcDecays->GetBinContent(min-i));
        }
   
    //------------------- decays pi- rebinned -----------------//
        TH1F* mc_pi_neg_d_rebin= rebin(mc_pi_neg_d);
        mc_pi_neg_d_rebin->SetName("mc_pi_neg_d_rebin");
        mc_pi_neg_d_rebin->SetTitle("mc_pi_neg_d rebinned");
        mc_pi_neg_d_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        mc_pi_neg_d_rebin->GetYaxis()->SetTitle("Entries");



    //--------------------------------------------------------- ITS+TPC ---------------------------------------------------------//
    
    //------------------- ITS+TPC pi+ -----------------//
        TH1F *mc_pi_pos_IT= new TH1F("mc_pi_pos_ITS_TPC", "pi+ ITS+TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        mc_pi_pos_IT->SetName("mc_ITS+TPC pi_pos");

        mc_pi_pos_IT->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        mc_pi_pos_IT->SetBinContent(i, mcITS_TPC->GetBinContent(min-1+i));
        }
        
    //------------------- ITS+TPC pi+ rebinned -----------------//
        TH1F* mc_pi_pos_IT_rebin= rebin(mc_pi_pos_IT);
        mc_pi_pos_IT_rebin->SetName("mc_pi_pos_IT_rebin");
        mc_pi_pos_IT_rebin->SetTitle("mc_pi_pos_IT rebinned");
        mc_pi_pos_IT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        mc_pi_pos_IT_rebin->GetYaxis()->SetTitle("Entries");


    //------------- ITS+TPC pi- -----------------------//
        TH1F *mc_pi_neg_IT= new TH1F("mc_pi_neg_ITS_TPC", "pi- ITS+TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        mc_pi_neg_IT->SetName("mc_ITS+TPC pi_neg");

        mc_pi_neg_IT->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        mc_pi_neg_IT->SetBinContent(i, mcITS_TPC->GetBinContent(min-i));
        }

    //------------------- ITS+TPC pi- rebinned -----------------//
        TH1F* mc_pi_neg_IT_rebin= rebin(mc_pi_neg_IT);
        mc_pi_neg_IT_rebin->SetName("mc_pi_neg_IT_rebin");
        mc_pi_neg_IT_rebin->SetTitle("mc_pi_neg_IT rebinned");
        mc_pi_neg_IT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        mc_pi_neg_IT_rebin->GetYaxis()->SetTitle("Entries");


    //--------------------------------------------------------- ITS with TPC ---------------------------------------------------------//
    
    //------------------- ITS with TPC pi+ -----------------//
        TH1F *mc_pi_pos_W= new TH1F("mc_pi_pos_ITSwithTPC", "pi+ ITS with TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        mc_pi_pos_W->SetName("mc_ITS with TPC pi_pos");

        mc_pi_pos_W->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        mc_pi_pos_W->SetBinContent(i, mcITSwithTPC->GetBinContent(min-1+i));
        }

    //------------------- ITS with TPC pi+ rebinned -----------------//
        TH1F* mc_pi_pos_W_rebin= rebin(mc_pi_pos_W);
        mc_pi_pos_W_rebin->SetName("mc_pi_pos_W_rebin");
        mc_pi_pos_W_rebin->SetTitle("mc_pi_pos_W rebinned");
        mc_pi_pos_W_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        mc_pi_pos_W_rebin->GetYaxis()->SetTitle("Entries");

    //------------- ITS with TPC pi- -----------------------//
        TH1F *mc_pi_neg_W= new TH1F("mc_pi_neg_ITSwithTPC", "pi- ITS with TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        mc_pi_neg_W->SetName("mc_ITS with TPC pi_neg");

        mc_pi_neg_W->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        mc_pi_neg_W->SetBinContent(i, mcITSwithTPC->GetBinContent(min-i));
        }
        
    //------------------- ITS with TPC pi- rebinned -----------------//
        TH1F* mc_pi_neg_W_rebin= rebin(mc_pi_neg_W);
        mc_pi_neg_W_rebin->SetName("mc_pi_neg_W_rebin");
        mc_pi_neg_W_rebin->SetTitle("mc_pi_neg_W rebinned");
        mc_pi_neg_W_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        mc_pi_neg_W_rebin->GetYaxis()->SetTitle("Entries");




    //--------------------------------- con data  ---------------------------------//
            //------------------- decays pi+ -----------------//
        TH1F *data_pi_pos_d= new TH1F("data_pi_pos_d", "pi+ decays;#it{p}_{T} (GeV/#it{c}); Entries",n,0.,massimox);
        data_pi_pos_d->SetName("data_decays_pi_pos");
   
        data_pi_pos_d->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        data_pi_pos_d->SetBinContent(i, dataDecays->GetBinContent(min-1+i));
        }

    //------------------- decays pi+ rebinned -----------------//
        TH1F* data_pi_pos_d_rebin= rebin(data_pi_pos_d);
        data_pi_pos_d_rebin->SetName("data_pi_pos_d_rebin");
        data_pi_pos_d_rebin->SetTitle("data_pi_pos_d rebinned");
        data_pi_pos_d_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        data_pi_pos_d_rebin->GetYaxis()->SetTitle("Entries");

    
    //------------- decays pi- -----------------------//
        TH1F *data_pi_neg_d= new TH1F("data_pi_neg_d", "pi- decays;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        data_pi_neg_d->SetName("data_decays_pi_neg");

        data_pi_neg_d->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        data_pi_neg_d->SetBinContent(i, dataDecays->GetBinContent(min-i));
        }
   
    //------------------- decays pi- rebinned -----------------//
        TH1F* data_pi_neg_d_rebin= rebin(data_pi_neg_d);
        data_pi_neg_d_rebin->SetName("data_pi_neg_d_rebin");
        data_pi_neg_d_rebin->SetTitle("data_pi_neg_d rebinned");
        data_pi_neg_d_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        data_pi_neg_d_rebin->GetYaxis()->SetTitle("Entries");


    //--------------------------------------------------------- ITS+TPC ---------------------------------------------------------//
    
    //------------------- ITS+TPC pi+ -----------------//
        TH1F *data_pi_pos_IT= new TH1F("data_pi_pos_ITS_TPC", "pi+ ITS+TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        data_pi_pos_IT->SetName("data_ITS+TPC pi_pos");

        data_pi_pos_IT->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        data_pi_pos_IT->SetBinContent(i, dataITS_TPC->GetBinContent(min-1+i));
        }
        
    //------------------- ITS+TPC pi+ rebinned -----------------//
        TH1F* data_pi_pos_IT_rebin= rebin(data_pi_pos_IT);
        data_pi_pos_IT_rebin->SetName("data_pi_pos_IT_rebin");
        data_pi_pos_IT_rebin->SetTitle("data_pi_pos_IT rebinned");
        data_pi_pos_IT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        data_pi_pos_IT_rebin->GetYaxis()->SetTitle("Entries");


    //------------- ITS+TPC pi- -----------------------//
        TH1F *data_pi_neg_IT= new TH1F("data_pi_neg_ITS_TPC", "pi- ITS+TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        data_pi_neg_IT->SetName("data_ITS+TPC pi_neg");

        data_pi_neg_IT->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        data_pi_neg_IT->SetBinContent(i, dataITS_TPC->GetBinContent(min-i));
        }

    //------------------- ITS+TPC pi- rebinned -----------------//
        TH1F* data_pi_neg_IT_rebin= rebin(data_pi_neg_IT);
        data_pi_neg_IT_rebin->SetName("data_pi_neg_IT_rebin");
        data_pi_neg_IT_rebin->SetTitle("data_pi_neg_IT rebinned");
        data_pi_neg_IT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        data_pi_neg_IT_rebin->GetYaxis()->SetTitle("Entries");


    //--------------------------------------------------------- ITS with TPC ---------------------------------------------------------//
    
    //------------------- ITS with TPC pi+ -----------------//
        TH1F *data_pi_pos_W= new TH1F("data_pi_pos_ITSwithTPC", "pi+ ITS with TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        data_pi_pos_W->SetName("data_ITS with TPC pi_pos");

        data_pi_pos_W->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        data_pi_pos_W->SetBinContent(i, dataITSwithTPC->GetBinContent(min-1+i));
        }

    //------------------- ITS with TPC pi+ rebinned -----------------//
        TH1F* data_pi_pos_W_rebin= rebin(data_pi_pos_W);
        data_pi_pos_W_rebin->SetName("data_pi_pos_W_rebin");
        data_pi_pos_W_rebin->SetTitle("data_pi_pos_W rebinned");
        data_pi_pos_W_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        data_pi_pos_W_rebin->GetYaxis()->SetTitle("Entries");

    //------------- ITS with TPC pi- -----------------------//
        TH1F *data_pi_neg_W= new TH1F("data_pi_neg_ITSwithTPC", "pi- ITS with TPC;#it{p}_{T} (GeV/#it{c});Entries",n,0.,massimox);
        data_pi_neg_W->SetName("data_ITS with TPC pi_neg");

        data_pi_neg_W->SetBinContent(1, 0);
        for(int i=2; i<=n; i++){
        data_pi_neg_W->SetBinContent(i, dataITSwithTPC->GetBinContent(min-i));
        }
        
    //------------------- ITS with TPC pi- rebinned -----------------//
        TH1F* data_pi_neg_W_rebin= rebin(data_pi_neg_W);
        data_pi_neg_W_rebin->SetName("data_pi_neg_W_rebin");
        data_pi_neg_W_rebin->SetTitle("data_pi_neg_W rebinned");
        data_pi_neg_W_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
        data_pi_neg_W_rebin->GetYaxis()->SetTitle("Entries");

        



        //------------------------- calcolo efficienze -------------------------//
    //------------------------- calcoli per p+ -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC --------------------//
        TCanvas *CE_glob_pos = new TCanvas("CE_glob_pos","Efficienza globale per traccia ITS+TPC pi+ - confronto tra dati veri e mc",500, 500);

        TCanvas *Cmc_E_glob_pos = new TCanvas("Cmc_E_glob_pos","Efficienza globale per traccia ITS+TPC pi+ mc",200,10,700,500);

        TH1F *gmc_E_glob_pos = new TH1F(*mc_pi_pos_IT_rebin);
        gmc_E_glob_pos->SetName("mc_efficiency_pi_pos_glob");
        gmc_E_glob_pos->Divide(mc_pi_pos_IT_rebin, mc_pi_pos_d_rebin, 1., 1., "B");
        gmc_E_glob_pos->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{+} con mc");
        gmc_E_glob_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gmc_E_glob_pos->GetYaxis()->SetTitle("#epsilon");
        Cmc_E_glob_pos->cd();
        gmc_E_glob_pos->Draw("PE");
    
        TCanvas *Cdata_E_glob_pos = new TCanvas("Cdata_E_glob_pos","Efficienza globale per traccia ITS+TPC pi+ dati veri",200,10,700,500);

        TH1F *gdata_E_glob_pos = new TH1F(*data_pi_pos_IT_rebin);
        gdata_E_glob_pos->SetName("data_efficiency_pi_pos_glob");
        gdata_E_glob_pos->Divide(data_pi_pos_IT_rebin, data_pi_pos_d_rebin, 1., 1., "B");
        gdata_E_glob_pos->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{+} con dati veri");
        gdata_E_glob_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gdata_E_glob_pos->GetYaxis()->SetTitle("#epsilon");
        Cdata_E_glob_pos->cd();
        gdata_E_glob_pos->Draw("PE");

        gdata_E_glob_pos->SetLineColor(2);
        gmc_E_glob_pos->SetLineStyle(7);

        gmc_E_glob_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gmc_E_glob_pos->GetXaxis()->SetRangeUser(0., 4.); 
        gdata_E_glob_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gdata_E_glob_pos->GetXaxis()->SetRangeUser(0., 4.); 

        gmc_E_glob_pos->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{+}");
        gStyle->SetTitleW(1.5);

        CE_glob_pos->cd();
        gmc_E_glob_pos->Draw("PE");
        gdata_E_glob_pos->Draw("PESAME");

        TLatex t1;
        t1.SetNDC();
        t1.SetTextSize(0.031);
        t1.DrawLatex(0.25, 0.42, "#scale[1.5]{Efficienze globali per traccia}");  //un po’ centrato
        t1.DrawLatex(0.25, 0.37, "#scale[1.5]{ ITS+TPC per #font[22]{#pi^{+}}}");
        t1.DrawLatex(0.25, 0.32, "#font[22]{#color[4]{#scale[2.2]{--}}} #scale[1.5]{MC}");
        t1.DrawLatex(0.25, 0.27, "#font[22]{#color[2]{#scale[2.2]{-}}} #scale[1.5]{Dati veri}");

        efficiency_mc_vs_data.cd();
        CE_glob_pos->Write();
        CE_glob_pos->Print("Efficienza globale per traccia ITS+TPC pi+ - confronto tra dati veri e mc.pdf");


    //------------------------------------- Efficienza di tracciamento TPC SA  -------------------------------------//
        TCanvas *CE_tpc_pos = new TCanvas("CE_tpc_pos","Efficienza di tracciamento TPC SA pi+ - confronto tra dati veri e mc",200,10,500,500);

        TH1F *mc_Copy_pi_pos_IT= new TH1F(*mc_pi_pos_IT_rebin);
        mc_Copy_pi_pos_IT->SetName("mc_copia_pi_pos_IT");
        mc_Copy_pi_pos_IT->Add(mc_pi_pos_W_rebin);

        TH1F *gmc_E_tpc_pos= new TH1F(*mc_Copy_pi_pos_IT);
        gmc_E_tpc_pos->SetName("mc_tpc_efficiency_pi_pos");
        gmc_E_tpc_pos->Divide(mc_Copy_pi_pos_IT, mc_pi_pos_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gmc_E_tpc_pos->SetTitle("Efficienza di tracciamento TPC SA #pi^{+} con mc");
        gmc_E_tpc_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gmc_E_tpc_pos->GetYaxis()->SetTitle("#epsilon");


        TH1F *data_Copy_pi_pos_IT= new TH1F(*data_pi_pos_IT_rebin);
        data_Copy_pi_pos_IT->SetName("data_copia_pi_pos_IT");
        data_Copy_pi_pos_IT->Add(data_pi_pos_W_rebin);

        TH1F *gdata_E_tpc_pos= new TH1F(*data_Copy_pi_pos_IT);
        gdata_E_tpc_pos->SetName("data_tpc_efficiency_pi_pos");
        gdata_E_tpc_pos->Divide(data_Copy_pi_pos_IT, data_pi_pos_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gdata_E_tpc_pos->SetTitle("Efficienza di tracciamento TPC SA #pi^{+} con mc");
        gdata_E_tpc_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gdata_E_tpc_pos->GetYaxis()->SetTitle("#epsilon");

        gdata_E_tpc_pos->SetLineColor(2);
        gmc_E_tpc_pos->SetLineStyle(7);

        gmc_E_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gmc_E_tpc_pos->GetXaxis()->SetRangeUser(0., 4.);
        gdata_E_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gdata_E_tpc_pos->GetXaxis()->SetRangeUser(0., 4.);

        gmc_E_tpc_pos->SetTitle("Efficienza di tracciamento TPC SA #pi^{+}");
        gStyle->SetTitleW(1.5);    

        CE_tpc_pos->cd();
        gmc_E_tpc_pos->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori
        gdata_E_tpc_pos->Draw("PESAME");

        TLatex t2;
        t2.SetNDC();
        t2.SetTextSize(0.031);
        t2.DrawLatex(0.25, 0.42, "#scale[1.5]{Efficienze di tracciamento}");  //un po’ centrato
        t2.DrawLatex(0.25, 0.37, "#scale[1.5]{TPC SA per #font[22]{#pi^{+}}}");  //un po’ centrato
        t2.DrawLatex(0.25, 0.32, "#font[22]{#color[4]{#scale[2.2]{--}}} #scale[1.5]{MC}");
        t2.DrawLatex(0.25, 0.27, "#font[22]{#color[2]{#scale[2.2]{-}}} #scale[1.5]{Dati veri}");

        efficiency_mc_vs_data.cd();
        CE_tpc_pos->Write();
        CE_tpc_pos->Print("Efficienza di tracciamento TPC SA pi+ - confronto tra dati veri e mc.pdf");


    //------------------- efficienza matching its-tpc --------------------//
        TCanvas *CE_its_tpc_pos = new TCanvas("CE_its_tpc_pos","Efficienza matching ITS TPC pi+ - confronto tra dati veri e mc",200,10,500,500);

        TH1F *mc_Copy_pi_pos_W= new TH1F(*mc_pi_pos_W_rebin);
        mc_Copy_pi_pos_W->SetName("mc_copia_pi_pos_W");
        mc_Copy_pi_pos_W->Add(mc_pi_pos_IT_rebin);

        TH1F *gmc_E_its_tpc_pos= new TH1F(*mc_Copy_pi_pos_W);
        gmc_E_its_tpc_pos->SetName("mc_its_tpc_efficiency_pi_pos");
        gmc_E_its_tpc_pos->Divide(mc_pi_pos_IT_rebin, mc_Copy_pi_pos_W, 1., 1., "B");   //B mi calcolal'errore binomiale
        gmc_E_its_tpc_pos->SetName("mc_Efficienza_matching_ITS_TPC_pi_pos");
        gmc_E_its_tpc_pos->SetTitle("Efficienza matching ITS TPC #pi^{+} con mc");
        gmc_E_its_tpc_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gmc_E_its_tpc_pos->GetYaxis()->SetTitle("#epsilon");


        TH1F *data_Copy_pi_pos_W= new TH1F(*data_pi_pos_W_rebin);
        data_Copy_pi_pos_W->SetName("data_copia_pi_pos_W");
        data_Copy_pi_pos_W->Add(data_pi_pos_IT_rebin);

        TH1F *gdata_E_its_tpc_pos= new TH1F(*data_Copy_pi_pos_W);
        gdata_E_its_tpc_pos->SetName("data_its_tpc_efficiency_pi_pos");
        gdata_E_its_tpc_pos->Divide(data_pi_pos_IT_rebin, data_Copy_pi_pos_W, 1., 1., "B");   //B mi calcolal'errore binomiale
        gdata_E_its_tpc_pos->SetName("data_Efficienza_matching_ITS_TPC_pi_pos");
        gdata_E_its_tpc_pos->SetTitle("Efficienza matching ITS TPC #pi^{+} con data");
        gdata_E_its_tpc_pos->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gdata_E_its_tpc_pos->GetYaxis()->SetTitle("#epsilon");

        gdata_E_its_tpc_pos->SetLineColor(2);
        gmc_E_its_tpc_pos->SetLineStyle(7);
        
        gmc_E_its_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gmc_E_its_tpc_pos->GetXaxis()->SetRangeUser(0., 4.);
        gdata_E_its_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
        gdata_E_its_tpc_pos->GetXaxis()->SetRangeUser(0., 4.);

        gmc_E_its_tpc_pos->SetTitle("Efficienza di matching ITS TPC #pi^{+}");
        gStyle->SetTitleW(1.5);  
    
        CE_its_tpc_pos->cd();
        gmc_E_its_tpc_pos->Draw("PE");
        gdata_E_its_tpc_pos->Draw("PESAME");

        TLatex t3;
        t3.SetNDC();
        t3.SetTextSize(0.031);
        t3.DrawLatex(0.2, 0.42, "#scale[1.5]{Efficienze di matching ITS TPC per #font[22]{#pi^{+}}}");  //un po’ centrato
        t3.DrawLatex(0.2, 0.37, "#font[22]{#color[4]{#scale[2.2]{--}}} #scale[1.5]{MC}");
        t3.DrawLatex(0.2, 0.32, "#font[22]{#color[2]{#scale[2.2]{-}}} #scale[1.5]{Dati veri}");


        efficiency_mc_vs_data.cd();
        CE_its_tpc_pos->Write();
        CE_its_tpc_pos->Print("efficienza matching ITS TPC pi+ - confronto tra dati veri e mc.pdf");





          //------------------------- calcoli per p- -------------------------//

    //------------------- Efficienza globale per traccia ITS+TPC --------------------//
        TCanvas *Cmc_E_glob_neg = new TCanvas("Cmc_E_glob_neg","Efficienza globale per traccia ITS+TPC pi- mc",200,10,700,500);

        TH1F *gmc_E_glob_neg = new TH1F(*mc_pi_neg_IT_rebin);
        gmc_E_glob_neg->SetName("mc_efficiency_pi_neg_glob");
        gmc_E_glob_neg->Divide(mc_pi_neg_IT_rebin, mc_pi_neg_d_rebin, 1., 1., "B");
        gmc_E_glob_neg->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{-} con mc");
        gmc_E_glob_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gmc_E_glob_neg->GetYaxis()->SetTitle("#epsilon");
        Cmc_E_glob_neg->cd();
        gmc_E_glob_neg->Draw("PE");
    
        TCanvas *Cdata_E_glob_neg = new TCanvas("Cdata_E_glob_neg","Efficienza globale per traccia ITS+TPC pi- dati veri",200,10,700,500);
        gStyle->SetTitleW(1.5);

        TH1F *gdata_E_glob_neg = new TH1F(*data_pi_neg_IT_rebin);
        gdata_E_glob_neg->SetName("data_efficiency_pi_neg_glob");
        gdata_E_glob_neg->Divide(data_pi_neg_IT_rebin, data_pi_neg_d_rebin, 1., 1., "B");
        gdata_E_glob_neg->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{-} con dati veri");
        gdata_E_glob_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gdata_E_glob_neg->GetYaxis()->SetTitle("#epsilon");
        Cdata_E_glob_neg->cd();
        gdata_E_glob_neg->Draw("PE");

        gdata_E_glob_neg->SetLineColor(2);

        gmc_E_glob_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gmc_E_glob_neg->GetXaxis()->SetRangeUser(0., 4.); 
        gdata_E_glob_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gdata_E_glob_neg->GetXaxis()->SetRangeUser(0., 4.); 

        gmc_E_glob_neg->SetTitle("Efficienza globale per traccia ITS+TPC #pi^{-} - confronto tra dati veri e mc");

        gmc_E_glob_neg->SetLineStyle(7);

        TCanvas *CE_glob_neg = new TCanvas("CE_glob_neg","Efficienza globale per traccia ITS+TPC pi- - confronto tra dati veri e mc",700, 500);

        CE_glob_neg->cd();
        gmc_E_glob_neg->Draw("PE");
        gdata_E_glob_neg->Draw("PESAME");

        TLatex t4;
        t4.SetNDC();
        t4.SetTextSize(0.031);
        t4.DrawLatex(0.25, 0.42, "#scale[1.5]{Efficienze globali per traccia ITS+TPC per #font[22]{#pi^{-}}}");  //un po’ centrato
        t4.DrawLatex(0.25, 0.37, "#font[22]{#color[4]{#scale[2.2]{--}}} #scale[1.5]{MC}");
        t4.DrawLatex(0.25, 0.32, "#font[22]{#color[2]{#scale[2.2]{-}}} #scale[1.5]{Dati veri}");

        efficiency_mc_vs_data.cd();
        CE_glob_neg->Write();
        CE_glob_neg->Print("Efficienza globale per traccia ITS+TPC pi- - confronto tra dati veri e mc.pdf");


    //------------------------------------- Efficienza di tracciamento TPC SA  -------------------------------------//
        TCanvas *CE_tpc_neg = new TCanvas("CE_tpc_neg","Efficienza di tracciamento TPC SA pi- - confronto tra dati veri e mc",200,10,700,500);

        TH1F *mc_Copy_pi_neg_IT= new TH1F(*mc_pi_neg_IT_rebin);
        mc_Copy_pi_neg_IT->SetName("mc_copia_pi_neg_IT");
        mc_Copy_pi_neg_IT->Add(mc_pi_neg_W_rebin);

        TH1F *gmc_E_tpc_neg= new TH1F(*mc_Copy_pi_neg_IT);
        gmc_E_tpc_neg->SetName("mc_tpc_efficiency_pi_neg");
        gmc_E_tpc_neg->Divide(mc_Copy_pi_neg_IT, mc_pi_neg_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gmc_E_tpc_neg->SetTitle("Efficienza di tracciamento TPC SA #pi^{-} con mc");
        gmc_E_tpc_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gmc_E_tpc_neg->GetYaxis()->SetTitle("#epsilon");


        TH1F *data_Copy_pi_neg_IT= new TH1F(*data_pi_neg_IT_rebin);
        data_Copy_pi_neg_IT->SetName("data_copia_pi_neg_IT");
        data_Copy_pi_neg_IT->Add(data_pi_neg_W_rebin);

        TH1F *gdata_E_tpc_neg= new TH1F(*data_Copy_pi_pos_IT);
        gdata_E_tpc_neg->SetName("data_tpc_efficiency_pi_neg");
        gdata_E_tpc_neg->Divide(data_Copy_pi_neg_IT, data_pi_neg_d_rebin, 1., 1., "B");   //B mi calcolal'errore binomiale
        gdata_E_tpc_neg->SetTitle("Efficienza di tracciamento TPC SA #pi^{-} con mc");
        gdata_E_tpc_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gdata_E_tpc_neg->GetYaxis()->SetTitle("#epsilon");

        gdata_E_tpc_neg->SetLineColor(2);
        gmc_E_tpc_neg->SetLineStyle(7);

        gmc_E_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gmc_E_tpc_neg->GetXaxis()->SetRangeUser(0., 4.);
        gdata_E_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gdata_E_tpc_neg->GetXaxis()->SetRangeUser(0., 4.);

        gmc_E_tpc_neg->SetTitle("Efficienza di tracciamento TPC SA #pi^{-} - confronto tra dati veri e mc");
        gStyle->SetTitleW(1.5);    

        CE_tpc_neg->cd();
        gmc_E_tpc_neg->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori
        gdata_E_tpc_neg->Draw("PESAME");

        TLatex t5;
        t5.SetNDC();
        t5.SetTextSize(0.031);
        t5.DrawLatex(0.25, 0.42, "#scale[1.5]{Efficienze di tracciamento TPC SA per #font[22]{#pi^{-}}}");  //un po’ centrato
        t5.DrawLatex(0.25, 0.37, "#font[22]{#color[4]{#scale[2.2]{--}}} #scale[1.5]{MC}");
        t5.DrawLatex(0.25, 0.32, "#font[22]{#color[2]{#scale[2.2]{-}}} #scale[1.5]{Dati veri}");

        efficiency_mc_vs_data.cd();
        CE_tpc_neg->Write();
        CE_tpc_neg->Print("Efficienza di tracciamento TPC SA pi- - confronto tra dati veri e mc.pdf");


    //------------------- efficienza matching its-tpc --------------------//
        TCanvas *CE_its_tpc_neg = new TCanvas("CE_its_tpc_neg","Efficienza matching ITS TPC pi- - confronto tra dati veri e mc",200,10,700,500);

        TH1F *mc_Copy_pi_neg_W= new TH1F(*mc_pi_neg_W_rebin);
        mc_Copy_pi_neg_W->SetName("mc_copia_pi_neg_W");
        mc_Copy_pi_neg_W->Add(mc_pi_neg_IT_rebin);

        TH1F *gmc_E_its_tpc_neg= new TH1F(*mc_Copy_pi_neg_W);
        gmc_E_its_tpc_neg->SetName("mc_its_tpc_efficiency_pi_neg");
        gmc_E_its_tpc_neg->Divide(mc_pi_neg_IT_rebin, mc_Copy_pi_neg_W, 1., 1., "B");   //B mi calcolal'errore binomiale
        gmc_E_its_tpc_neg->SetName("mc_Efficienza_matching_ITS_TPC_pi_neg");
        gmc_E_its_tpc_neg->SetTitle("Efficienza matching ITS TPC #pi^{-} con mc");
        gmc_E_its_tpc_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gmc_E_its_tpc_neg->GetYaxis()->SetTitle("#epsilon");


        TH1F *data_Copy_pi_neg_W= new TH1F(*data_pi_neg_W_rebin);
        data_Copy_pi_neg_W->SetName("data_copia_pi_neg_W");
        data_Copy_pi_neg_W->Add(data_pi_neg_IT_rebin);

        TH1F *gdata_E_its_tpc_neg= new TH1F(*data_Copy_pi_neg_W);
        gdata_E_its_tpc_neg->SetName("data_its_tpc_efficiency_pi_neg");
        gdata_E_its_tpc_neg->Divide(data_pi_neg_IT_rebin, data_Copy_pi_neg_W, 1., 1., "B");   //B mi calcolal'errore binomiale
        gdata_E_its_tpc_neg->SetName("data_Efficienza_matching_ITS_TPC_pi_neg");
        gdata_E_its_tpc_neg->SetTitle("Efficienza matching ITS TPC #pi^{-} con data");
        gdata_E_its_tpc_neg->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        gdata_E_its_tpc_neg->GetYaxis()->SetTitle("#epsilon");

        gdata_E_its_tpc_neg->SetLineColor(2);
        gmc_E_its_tpc_neg->SetLineStyle(7);
        
        gmc_E_its_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gmc_E_its_tpc_neg->GetXaxis()->SetRangeUser(0., 4.);
        gdata_E_its_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
        gdata_E_its_tpc_neg->GetXaxis()->SetRangeUser(0., 4.);

        gmc_E_its_tpc_neg->SetTitle("Efficienza di matching ITS TPC #pi^{-} - confronto tra dati veri e mc");
        gStyle->SetTitleW(1.5);  
    
        CE_its_tpc_neg->cd();
        gmc_E_its_tpc_neg->Draw("PE");
        gdata_E_its_tpc_neg->Draw("PESAME");

        TLatex t6;
        t6.SetNDC();
        t6.SetTextSize(0.031);
        t6.DrawLatex(0.25, 0.42, "#scale[1.5]{Efficienze di matching ITS TPC per #font[22]{#pi^{-}}}");  //un po’ centrato
        t6.DrawLatex(0.25, 0.37, "#font[22]{#color[4]{#scale[2.2]{--}}} #scale[1.5]{MC}");
        t6.DrawLatex(0.25, 0.32, "#font[22]{#color[2]{#scale[2.2]{-}}} #scale[1.5]{Dati veri}");


        efficiency_mc_vs_data.cd();
        CE_its_tpc_neg->Write();
        CE_its_tpc_neg->Print("efficienza matching ITS TPC pi- - confronto tra dati veri e mc.pdf");


}