#include "Config.h"

using namespace std;


void estremi(TH1D *hist, float a, float b, float &m, float &M)
{
    float aa= hist->GetXaxis()->FindBin(a);
    m= hist->GetXaxis()->GetBinLowEdge(aa);

    float bb= hist->GetXaxis()->FindBin(b);
    M= hist->GetXaxis()->GetBinUpEdge(bb);
}






void double_fit(TH1D *hist, int i, TH1D *h_cov_matrix_status, TH1D *h_chi2, TH1D *h_prob, TH1D *h_cov_matrix_status2, TH1D *h_chi2_2, TH1D *h_prob2, float &a, float &sa, TF1 *f1, TF1 *f2)
{
    
        float mean, dev_std, mean2, dev_std2;
        float low_edge, up_edge, low_edge2, up_edge2;


        float estr_inf, estr_sup;
        mean= hist->GetMean(1);   //int axis=1 mi dà il valore medio sull'asse X (2<->Y,  3<->Z)
        dev_std= hist->GetStdDev(1);
        low_edge= mean-dev_std;
        up_edge= mean+dev_std;
        estremi(hist, low_edge, up_edge, estr_inf, estr_sup);
        cout<<endl;
        cout<<"\nlow edge: "<<low_edge<<endl;
        cout<<"up edge: "<<up_edge<<endl;
        cout<<"estr_inf: "<<estr_inf<<endl;
        cout<<"estr_sup: "<<estr_sup<<endl;
        cout<<endl;

        
        TFitResultPtr FitResult=hist->Fit("gaus", "RMLS", "", estr_inf, estr_sup);    //specifico il range su cui effettuare il fit
        cout<<"\nStato della matrice di covarianza: "<<FitResult->CovMatrixStatus()<<endl;
        cout<<"Chi2: "<<FitResult->Chi2()<<endl;
        cout<<"Probabilità: "<<FitResult->Prob()<<endl;

        //grafico stato del fit, stato della matrice di covarianza, chi2 e probabilità
        h_cov_matrix_status->SetBinContent(i, FitResult->CovMatrixStatus());
        h_chi2->SetBinContent(i, FitResult->Chi2());
        h_prob->SetBinContent(i, FitResult->Prob());
       
         


        hist->GetFunction("gaus")->SetName("gaus1");
        f1 = (TF1*)hist->GetListOfFunctions()->FindObject("gaus1");
 
        cout<<"\n-------------- Secondo fit --------------"<<endl;
        float estr_inf2, estr_sup2;
        mean2=hist->GetFunction("gaus1")->GetParameter(1);     //f(x) = p0*exp(-0.5*((x-p1)/p2)^2) 
        dev_std2=hist->GetFunction("gaus1")->GetParameter(2); 
        low_edge2= mean2-dev_std2;
        up_edge2= mean2+dev_std2;
        estremi(hist, low_edge2, up_edge2, estr_inf2, estr_sup2);
        cout<<"\nlow edge2: "<<low_edge2<<endl;
        cout<<"up edge2: "<<up_edge2<<endl;
        cout<<"estr_inf2: "<<estr_inf2<<endl;
        cout<<"estr_sup2: "<<estr_sup2<<endl;
        cout<<endl;

        TFitResultPtr FitResult2 = hist->Fit("gaus", "RMLS+", "",estr_inf2, estr_sup2);     //+ mi fitta questa gaussiana sul "grafico" in cui c'è già l'altra e non la sovrascrive
        cout<<"\nStato della matrice di covarianza: "<<FitResult2->CovMatrixStatus()<<endl;
        cout<<"Chi2: "<<FitResult2->Chi2()<<endl;
        cout<<"Probabilità: "<<FitResult2->Prob()<<endl;

        h_cov_matrix_status2->SetBinContent(i, FitResult2->CovMatrixStatus());
        h_chi2_2->SetBinContent(i, FitResult2->Chi2());
        h_prob2->SetBinContent(i, FitResult2->Prob());



        hist->GetFunction("gaus")->SetName("gaus2");
        a= hist->GetFunction("gaus2")->GetParameter(2);
        sa= hist->GetFunction("gaus2")->GetParError(2);
   

        
        f2 = (TF1*)hist->GetListOfFunctions()->FindObject("gaus2");
        f2->SetLineColor(4);

        
/*         if(i==14){
        TCanvas* Cproiezione= new TCanvas("Cproiezione", "proiezione - (1.3, 1.4) Gev/#it{c}", 200,50,600,400); 
        gStyle->SetOptStat(0);

        Cproiezione->cd(0);
        hist->GetXaxis()->SetRangeUser(estr_inf, estr_sup);
        hist->Draw();

        TLegend *leg= new TLegend(0.73,0.85,0.88,0.65);
        leg->SetHeader("Fit gaussiani", "C");  
        leg->AddEntry(f1, "primo fit","l");
        leg->AddEntry(f2, "secondo fit","l");
        leg->SetTextSize(0.038);
        leg->Draw();

        Cproiezione->Print("proiezione_1.3_1.4.pdf");
        }


        if(i==3){
        TCanvas* Cproiezione1= new TCanvas("Cproiezione1", "proiezione - (0.2, 0.3) Gev/#it{c}", 200,50,600,400); 
        gStyle->SetOptStat(0);

        Cproiezione1->cd(0);
        hist->GetXaxis()->SetRangeUser(-0.12, 0.12);
        hist->Draw();

        TLegend *leg1= new TLegend(0.73,0.85,0.88,0.65);
        leg1->SetHeader("Fit gaussiani", "C");   
        leg1->AddEntry(f1, "primo fit","l");
        leg1->AddEntry(f2, "secondo fit","l");
        leg1->SetTextSize(0.038);
        leg1->Draw();

        Cproiezione1->Print("proiezione_0.2_0.3.pdf");
       }
 */


}





TH2F* add(TH2 *start_hist, int ny)
{
    
    //prendo l'istogramma 2D e cerco il binning del suo asse x e del suo asse y-> mi servono per sommare la parte positiva e quella negativa
    int nbinx= start_hist->GetXaxis()->GetLast();
    float width_x= start_hist->GetXaxis()->GetBinWidth(1);
    float *bin_ledgex = new float[nbinx];
    float minx;
    int nx;     
    for(int i=1; i<=nbinx; i++){
        bin_ledgex[i]= start_hist->GetXaxis()->GetBinLowEdge(i);
        if( abs(bin_ledgex[i]-0.)<width_x/2){
            minx=i;
            cout<<"il bin minimo è: "<<minx<<endl;
            nx= minx-1;
            break;
        }                               
    }

    cout<<"nx: "<<nx<<endl;
    float max_x= start_hist->GetXaxis()->GetBinUpEdge(nbinx);

    float width_y= start_hist->GetYaxis()->GetBinWidth(1);
    float max_y= start_hist->GetYaxis()->GetBinUpEdge(ny);
    float min_y= start_hist->GetYaxis()->GetBinLowEdge(1); 


    // prendo l'istogramma 2D solo per pt positivi e sommo i bin speculari rispetto allo zero-> quello che ottengo è un TH2F i cui conteggi contengono sia pi+ che pi-

    TH2F* hist_add= new TH2F("hist_add", "#hist(#it{p}_{T}) - somma dati #pi+ e #pi-", nx, 0., max_x, ny, min_y, max_y);
    hist_add->SetName("somma dati #pi+ e #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=ny; j++){
            hist_add->SetBinContent(i, j, start_hist->GetBinContent(nx+i, j));
        }
    }


    TH2F* hist_pi_neg= new TH2F("pi-", "#eta(#it{p}_{T}) - #pi-", nx, 0., max_x, ny, min_y, max_y);
    hist_pi_neg->SetName("eta - #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=ny; j++){
            hist_pi_neg->SetBinContent(i, j, start_hist->GetBinContent(nx+1-i, j));
        }
    }


 
    TH2F *g_hist_add= new TH2F(*hist_add);   //ci metto sopra un grafico   sto passando l'oggetto puntato, ne sto facendo una copia
    gStyle->SetOptStat(0);
    g_hist_add->Add(hist_pi_neg);


    return g_hist_add;
}






TH2F* rebin(TH2F *g_hist_add, int ny, const double ybins[])
{
    //sto riempendo il mio istogramma rebinnato
    TH2F* h_rebin= new TH2F("h_rebin", "hist rebinned", kNbinsx, k_xbins, ny, ybins);    //istogramma g_eta_add rebinnato   
    for(int i=1; i<=kNbinsx; i++){     //loop su x nuovo
        int xmin= g_hist_add->GetXaxis()->FindBin(k_xbins[i-1]+epsilon);
        int xmax= g_hist_add->GetXaxis()->FindBin(k_xbins[i]-epsilon);
        //printf("xmin = %d, xmax = %d, k_xbins[j-1] = %f, k_xbins[j] = %f\n", xmin, xmax, k_xbins[i-1], k_xbins[i]);
        for(int j=1; j<=ny; j++){   //loop su y nuovo
            double integral= g_hist_add->Integral(xmin, xmax, j, j);
            h_rebin->SetBinContent(i, j, integral);
        }
    } 

    return h_rebin;
}














void eta(const double nsigma=1.){

    TFile fAnalysis("AnalysisResultsLHC23k2f.root");
    TFile parameterisation_study("parameterisation_study.root", "RECREATE");
    parameterisation_study.mkdir("delta_eta");
    parameterisation_study.mkdir("delta_phi");
    parameterisation_study.mkdir("delta_pT");

    parameterisation_study.mkdir("delta_eta/Projections_eta");
    parameterisation_study.mkdir("delta_phi/Projections_phi");
    parameterisation_study.mkdir("delta_pT/Projections_pT");

    

//------------------------------------------------------------------------------ eta ------------------------------------------------------------------------------
    TH2F* etaTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/etaTpcIts");

    TCanvas *Ceta_add= new TCanvas("Ceta_add", "eta(p_T) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    gStyle->SetOptStat(0);

    TH2F *g_eta_add = add(etaTpcIts, kNbinsy_eta);
    g_eta_add->SetName("#Delta#eta(#it{p}_{T}) - somma dati #pi+ e #pi-");
    g_eta_add->SetTitle("#Delta#eta(#it{p}_{T})");
    g_eta_add->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    g_eta_add->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");


    TH2F* h_deltaeta_rebin= rebin(g_eta_add, kNbinsy_eta, k_ybins_eta);  
    h_deltaeta_rebin->SetName("h_deltaeta_rebin");
    h_deltaeta_rebin->SetTitle("#Delta#eta(#it{p}_{T}) rebinned");
    h_deltaeta_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
    h_deltaeta_rebin->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");
    


    Ceta_add->cd();
    g_eta_add->Draw();

    parameterisation_study.cd("delta_eta");   
    Ceta_add->Write();
    h_deltaeta_rebin->Write();


 
//------------------------------------------------------------------------------ pT ------------------------------------------------------------------------------
    TH2F* ptTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/ptTpcIts");

    TCanvas *CpT_add= new TCanvas("CpT_add", "pT(p_T) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    gStyle->SetOptStat(0);

    TH2F *g_pT_add = add(ptTpcIts, kNbinsy_pT);
    g_pT_add->SetName("#Delta#it{p}_{T}(#it{p}_{T}) - somma dati #pi+ e #pi-");
    g_pT_add->SetTitle("#Delta#it{p}_{T}(#it{p}_{T})");
    g_pT_add->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    g_pT_add->GetYaxis()->SetTitle("#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (Gev/#it{c})");


    TH2F* h_deltapT_rebin= rebin(g_pT_add, kNbinsy_pT, k_ybins_pT);
    h_deltapT_rebin->SetName("h_deltapT_rebin"); 
    h_deltapT_rebin->SetTitle("#Delta#it{p}_{T}(#it{p}_{T}) rebinned");
    h_deltapT_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
    h_deltapT_rebin->GetYaxis()->SetTitle("#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (Gev/#it{c})");
 


    CpT_add->cd();
    g_pT_add->Draw();

    parameterisation_study.cd("delta_pT");   
    CpT_add->Write();
    h_deltapT_rebin->Write();


 
//------------------------------------------------------------------------------ phi ------------------------------------------------------------------------------
    TH2F* phiTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/phiTpcIts");

    TCanvas *Cphi_add= new TCanvas("Cphi_add", "#phi(p_T) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    gStyle->SetOptStat(0);

    TH2F *g_phi_add = add(phiTpcIts, kNbinsy_phi);
    g_phi_add->SetName("#Delta#phi(#it{p}_{T}) - somma dati #pi+ e #pi-");
    g_phi_add->SetTitle("#Delta#phi(#it{p}_{T})");
    g_phi_add->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    g_phi_add->GetYaxis()->SetTitle("#phi^{TPC} - #phi^{ITS} (rad)");


    TH2F* h_deltaphi_rebin= rebin(g_phi_add, kNbinsy_phi, k_ybins_phi);  
    h_deltaphi_rebin->SetName("h_deltaphi_rebin");
    h_deltaphi_rebin->SetTitle("#Delta#phi(#it{p}_{T}) rebinned");
    h_deltaphi_rebin->GetXaxis()->SetTitle("#it{p}_{T} (Gev/#it{c})");
    h_deltaphi_rebin->GetYaxis()->SetTitle("#phi^{TPC} - #phi^{ITS} (rad)");
 


    Cphi_add->cd();
    g_phi_add->Draw();

    parameterisation_study.cd("delta_phi");   
    Cphi_add->Write();
    h_deltaphi_rebin->Write();
 




   
//ora devo fare un ciclo for in cui prendo proiezione sull'asse y dell'istogramma 2D, faccio due fit gaussiani e grafico dev std al variare di pT

//------------------------------------------------------------------------------ eta ------------------------------------------------------------------------------
    TH1F* hist_full_eta= new TH1F("delta_eta", "#sigma_{#Delta#eta}(#it{p}_{T})", kNbinsx, k_xbins);                 
    hist_full_eta->GetYaxis()->SetTitle("#sigma_{#Delta#eta}");
    hist_full_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    
    TH1D** h_proj_eta= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status_eta= new TH1D("CovMatrixStatus_eta", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_eta= new TH1D("chi2_eta", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob_eta= new TH1D("prob_eta", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2_eta= new TH1D("CovMatrixStatus2_eta", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status2_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_2_eta= new TH1D("chi2_2_eta", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_2_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob2_eta= new TH1D("prob2_eta", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob2_eta->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    int entries_eta;

    TF1 *f1_eta;
    TF1 *f2_eta;

    float delta_eta, sdelta_eta;



    //------------------------------------------------------------------------------ pT ------------------------------------------------------------------------------
    TH1F* hist_full_pT= new TH1F("delta_pT", "#sigma_{#Delta#it{p}_{T}}(#it{p}_{T})", kNbinsx, k_xbins);                 
    hist_full_pT->GetYaxis()->SetTitle("#sigma_{#Delta#it{p}_{T}}");
    hist_full_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    
    TH1D** h_proj_pT= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status_pT= new TH1D("CovMatrixStatus_pT", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_pT= new TH1D("chi2_pT", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob_pT= new TH1D("prob_pT", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2_pT= new TH1D("CovMatrixStatus2_pT", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status2_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_2_pT= new TH1D("chi2_2_pT", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_2_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob2_pT= new TH1D("prob2_pT", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob2_pT->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    int entries_pT;

    TF1 *f1_pT;
    TF1 *f2_pT;

    float delta_pT, sdelta_pT;



//------------------------------------------------------------------------------ phi ------------------------------------------------------------------------------
    TH1F* hist_full_phi= new TH1F("delta_phi", "#sigma_{#Delta#phi}(#it{p}_{T})", kNbinsx, k_xbins);                 
    hist_full_phi->GetYaxis()->SetTitle("#sigma_{#Delta#phi}");
    hist_full_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    
    TH1D** h_proj_phi= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status_phi= new TH1D("CovMatrixStatus_phi", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_phi= new TH1D("chi2_phi", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob_phi= new TH1D("prob_phi", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2_phi= new TH1D("CovMatrixStatus2_phi", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status2_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_2_phi= new TH1D("chi2_2_phi", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_2_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob2_phi= new TH1D("prob2_phi", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob2_phi->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    int entries_phi;

    TF1 *f1_phi;
    TF1 *f2_phi;

    float delta_phi, sdelta_phi;





    for(int i=1; i<=kNbinsx; i++){

        cout<<"-----------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"------------------------------------------------------------- eta -------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in pt: ("<<h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

        h_proj_eta[i-1]=h_deltaeta_rebin->ProjectionY(Form("ProjectionY_%.1f_%.1f", h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)), i, i);  
        //questa parte mi dava problemi perché, quando ci sono pochi o nessun dato, il fit, naturalmente, non viene ma mi blocca tutto quanto e mi dà cose brutte anche per ricavare lo stato del fit e della matrice di covarianza
        entries_eta=h_proj_eta[i-1]->GetEntries();
        if(entries_eta>=5){

        double_fit(h_proj_eta[i-1], i, h_cov_matrix_status_eta, h_chi2_eta, h_prob_eta, h_cov_matrix_status2_eta, h_chi2_2_eta, h_prob2_eta, delta_eta, sdelta_eta, f1_eta, f2_eta);
        
        cout<<"\nDelta eta: "<<delta_eta<<" +- "<<sdelta_eta<<endl;

        hist_full_eta->SetBinContent(i, delta_eta);
        hist_full_eta->SetBinError(i, sdelta_eta);

        h_proj_eta[i-1]->GetXaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");
        h_proj_eta[i-1]->GetYaxis()->SetTitle("Entries");
        h_proj_eta[i-1]->SetTitle(Form("ProjectionY #Delta#eta - Intervallo in #it{p}_{T} = (%.1f, %.1f) Gev/#it{c}", h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)));


        parameterisation_study.cd("delta_eta/Projections_eta");
        h_proj_eta[i-1]->Write();
        }


        
        
        cout<<"------------------------------------------------------------- pT -------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in pt: ("<<h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltapT_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

        h_proj_pT[i-1]=h_deltapT_rebin->ProjectionY(Form("ProjectionY_%.1f_%.1f", h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i), h_deltapT_rebin->GetXaxis()->GetBinUpEdge(i)), i, i);  
        //questa parte mi dava problemi perché, quando ci sono pochi o nessun dato, il fit, naturalmente, non viene ma mi blocca tutto quanto e mi dà cose brutte anche per ricavare lo stato del fit e della matrice di covarianza
        entries_pT=h_proj_pT[i-1]->GetEntries();
        if(entries_pT>=5){

        double_fit(h_proj_pT[i-1], i, h_cov_matrix_status_pT, h_chi2_pT, h_prob_pT, h_cov_matrix_status2_pT, h_chi2_2_pT, h_prob2_pT, delta_pT, sdelta_pT, f1_pT, f2_pT);
        
        cout<<"\nDelta pT: "<<delta_pT<<" +- "<<sdelta_pT<<endl;

        hist_full_pT->SetBinContent(i, delta_pT);
        hist_full_pT->SetBinError(i, sdelta_pT);

        h_proj_pT[i-1]->GetXaxis()->SetTitle("#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (Gev/#it{c})");
        h_proj_pT[i-1]->GetYaxis()->SetTitle("Entries");
        h_proj_pT[i-1]->SetTitle(Form("ProjectionY #Delta#it{p}_{T} - Intervallo in #it{p}_{T} = (%.1f, %.1f) Gev/#it{c}", h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i), h_deltapT_rebin->GetXaxis()->GetBinUpEdge(i)));


        parameterisation_study.cd("delta_pT/Projections_pT");
        h_proj_pT[i-1]->Write();
        }




        cout<<"------------------------------------------------------------- phi -------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in pt: ("<<h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltaphi_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

        h_proj_phi[i-1]=h_deltaphi_rebin->ProjectionY(Form("ProjectionY_%.1f_%.1f", h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaphi_rebin->GetXaxis()->GetBinUpEdge(i)), i, i);  
        //questa parte mi dava problemi perché, quando ci sono pochi o nessun dato, il fit, naturalmente, non viene ma mi blocca tutto quanto e mi dà cose brutte anche per ricavare lo stato del fit e della matrice di covarianza
        entries_phi=h_proj_phi[i-1]->GetEntries();
        if(entries_phi>=5){

        double_fit(h_proj_phi[i-1], i, h_cov_matrix_status_phi, h_chi2_phi, h_prob_phi, h_cov_matrix_status2_phi, h_chi2_2_phi, h_prob2_phi, delta_phi, sdelta_phi, f1_phi, f2_phi);
        
        cout<<"\nDelta phi: "<<delta_phi<<" +- "<<sdelta_phi<<endl;

        hist_full_phi->SetBinContent(i, delta_phi);
        hist_full_phi->SetBinError(i, sdelta_phi);

        h_proj_phi[i-1]->GetXaxis()->SetTitle("#phi^{TPC} - #phi^{ITS} (rad)");
        h_proj_phi[i-1]->GetYaxis()->SetTitle("Entries");
        h_proj_phi[i-1]->SetTitle(Form("ProjectionY #Delta#phi - Intervallo in #it{p}_{T} = (%.1f, %.1f) Gev/#it{c}", h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaphi_rebin->GetXaxis()->GetBinUpEdge(i)));


        parameterisation_study.cd("delta_phi/Projections_phi");
        h_proj_phi[i-1]->Write();
        }
       

    }




    cout<<endl;
    cout<<endl;
    cout<<"-------------------------- eta - Fit con power law --------------------------"<<endl;
    TF1* power_law= new TF1("power_law", "[0]*TMath::Power(x, [1])");
    TFitResultPtr FitResult_eta= hist_full_eta->Fit(power_law, "RMS", "0", 0.2, 5.);    //R: in the specified range
    cout<<"\nCHi2: "<<FitResult_eta->Chi2()<<endl;
    cout<<"Prob: "<<FitResult_eta->Prob()<<endl;


    TCanvas *plot_eta= new TCanvas("Cplot_eta", "plot_eta", 200,50,600,400);
    hist_full_eta->SetStats(1);   //activates the statistics box

    //plot_eta->Update();
    plot_eta->cd();
    hist_full_eta->Draw();
    gStyle->SetOptFit(1111);
    power_law->SetRange(0., 5.);
    power_law->Draw("SAME");

    TLatex t_eta;
    t_eta.SetNDC();
    t_eta.SetTextSize(0.071);
    t_eta.DrawLatex(0.41, 0.50, "#it{f}(x) = p0 x^{p1}");


    parameterisation_study.cd("delta_eta");
    plot_eta->Write();
    plot_eta->Print("plot_eta.pdf");
    

 




    cout<<endl;
    cout<<"-------------------------- pT - Fit con parabola --------------------------"<<endl;
    TF1* parabola= new TF1("parabola", "pol2");
    TFitResultPtr FitResult_pT= hist_full_pT->Fit(parabola, "RMS", "0", 0.1, 2.8);
    cout<<"\nCHi2: "<<FitResult_pT->Chi2()<<endl;
    cout<<"Prob: "<<FitResult_pT->Prob()<<endl;


    TCanvas *plot_pT= new TCanvas("Cplot_pT", "plot_pT", 200,50,600,400);
    hist_full_pT->SetStats(1);   //activates the statistics box

    //plot_pT->Update();
    plot_pT->cd();
    hist_full_pT->Draw();
    gStyle->SetOptFit(1111);
    parabola->SetRange(0., 5.);
    parabola->Draw("SAME");

    TLatex t_pT;
    t_pT.SetNDC();
    t_pT.SetTextSize(0.071);
    t_pT.DrawLatex(0.28, 0.50, "#it{f}(x) = p0 + p1 x + p2 x^{2}");


    parameterisation_study.cd("delta_pT");
    plot_pT->Write();
    plot_pT->Print("plot_pT.pdf");


 




    cout<<endl;
    cout<<"-------------------------- phi - Fit con esponenziale negativo --------------------------"<<endl;
    TF1* exp_neg= new TF1("exp_neg", "[0]+[1]*exp(-([2]*x))");
    TFitResultPtr FitResult_phi= hist_full_phi->Fit(exp_neg, "RMS", "", 0.1, 1.5);
    cout<<"\nCHi2: "<<FitResult_phi->Chi2()<<endl;
    cout<<"Prob: "<<FitResult_phi->Prob()<<endl;
 

    TCanvas *plot_phi= new TCanvas("Cplot_phi", "plot_phi", 200,50,600,400);
    hist_full_phi->SetStats(1);   //activates the statistics box

    plot_phi->Update();
    hist_full_phi->Draw();
    plot_phi->cd();
    gStyle->SetOptFit(1111);
    exp_neg->SetRange(0., 5.);
    exp_neg->Draw("SAME");

    TLatex t_phi;
    t_phi.SetNDC();
    t_phi.SetTextSize(0.071);
    t_phi.DrawLatex(0.36, 0.50, "#it{f}(x) = p0 + p1 e^{p2 x}");


    parameterisation_study.cd("delta_phi");
    plot_phi->Write();
    plot_phi->Print("plot_phi.pdf");







    
 
 


    


    parameterisation_study.cd("delta_eta");
    h_cov_matrix_status_eta->Write();
    h_chi2_eta->Write();
    h_prob_eta->Write();
    h_cov_matrix_status2_eta->Write();
    h_chi2_2_eta->Write();
    h_prob2_eta->Write();

    

    parameterisation_study.cd("delta_pT");
    h_cov_matrix_status_pT->Write();
    h_chi2_pT->Write();
    h_prob_pT->Write();
    h_cov_matrix_status2_pT->Write();
    h_chi2_2_pT->Write();
    h_prob2_pT->Write();



    parameterisation_study.cd("delta_phi");
    h_cov_matrix_status_phi->Write();
    h_chi2_phi->Write();
    h_prob_phi->Write();
    h_cov_matrix_status2_phi->Write();
    h_chi2_2_phi->Write();
    h_prob2_phi->Write();

    parameterisation_study.Close();
 
      
}