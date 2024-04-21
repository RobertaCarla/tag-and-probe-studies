#include "Config.h"
//è una copia di eta.cpp, direi più aggiornata--> mi conviene usare questa
using namespace std;


void estremi(TH1D *hist, float a, float b, float &m, float &M)
{
    float aa= hist->GetXaxis()->FindBin(a);
    m= hist->GetXaxis()->GetBinLowEdge(aa);

    float bb= hist->GetXaxis()->FindBin(b);
    M= hist->GetXaxis()->GetBinUpEdge(bb);
}






void double_fit(TH1D *hist, int i, TH1D *h_cov_matrix_status, TH1D *h_chi2, TH1D *h_prob, TH1D *h_cov_matrix_status2, TH1D *h_chi2_2, TH1D *h_prob2, float &a, float &sa, TF1 *f1, TF1 *f2, float nsigma=1.f)
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

        
        TFitResultPtr FitResult=hist->Fit("gaus", "RMLSI", "", estr_inf, estr_sup);    //specifico il range su cui effettuare il fit
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
        low_edge2= mean2-nsigma*dev_std2;
        up_edge2= mean2+nsigma*dev_std2;
        estremi(hist, low_edge2, up_edge2, estr_inf2, estr_sup2);
        cout<<"\nlow edge2: "<<low_edge2<<endl;
        cout<<"up edge2: "<<up_edge2<<endl;
        cout<<"estr_inf2: "<<estr_inf2<<endl;
        cout<<"estr_sup2: "<<estr_sup2<<endl;
        cout<<endl;

        TFitResultPtr FitResult2 = hist->Fit("gaus", "RMLSI+", "",estr_inf2, estr_sup2);     //+ mi fitta questa gaussiana sul "grafico" in cui c'è già l'altra e non la sovrascrive
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

    TH2F* hist_add= new TH2F("hist_add", "#hist(#it{p}^{ITS}) - somma dati #pi+ e #pi-", nx, 0., max_x, ny, min_y, max_y);
    hist_add->SetName("somma dati #pi+ e #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=ny; j++){
            hist_add->SetBinContent(i, j, start_hist->GetBinContent(nx+i, j));
        }
    }


    TH2F* hist_pi_neg= new TH2F("pi-", "#eta(#it{p}^{ITS}) - #pi-", nx, 0., max_x, ny, min_y, max_y);
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














void parameter(){

    TFile fAnalysis("AnalysisResults_LHC23k2f_20240327.root");
    TFile parameterisation_study("parameterisation_study.root", "RECREATE");
    parameterisation_study.mkdir("delta_eta");
    parameterisation_study.mkdir("delta_phi");
    parameterisation_study.mkdir("delta_pT");

    parameterisation_study.mkdir("delta_eta/Projections_eta");
    parameterisation_study.mkdir("delta_phi/Projections_phi");
    parameterisation_study.mkdir("delta_pT/Projections_pT");



//------------------------------------------------------------------------------ eta ------------------------------------------------------------------------------
    TH2F* etaTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/etaTpcIts");

    TCanvas *Ceta_add= new TCanvas("Ceta_add", "eta(#it{p}^{ITS}) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    gStyle->SetOptStat(0);

    TH2F *g_eta_add = add(etaTpcIts, kNbinsy_eta);
    g_eta_add->SetName("#Delta#eta(#it{p}^{ITS}) - somma dati #pi+ e #pi-");
    g_eta_add->SetTitle("#Delta#eta(#it{p}^{ITS})");
    g_eta_add->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    g_eta_add->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");


    TH2F* h_deltaeta_rebin= rebin(g_eta_add, kNbinsy_eta, k_ybins_eta);  
    h_deltaeta_rebin->SetName("h_deltaeta_rebin");
    h_deltaeta_rebin->SetTitle("#Delta#eta(#it{p}^{ITS}) rebinned");
    h_deltaeta_rebin->GetXaxis()->SetTitle("#it{p}^{ITS} (Gev/#it{c})");
    h_deltaeta_rebin->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");
    


    Ceta_add->cd();
    g_eta_add->Draw();

    parameterisation_study.cd("delta_eta");   
    Ceta_add->Write();
    h_deltaeta_rebin->Write();


 
//------------------------------------------------------------------------------ pT ------------------------------------------------------------------------------
    TH2F* ptTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/ptTpcIts");

    TCanvas *CpT_add= new TCanvas("CpT_add", "pT(p) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    gStyle->SetOptStat(0);

    TH2F *g_pT_add = add(ptTpcIts, kNbinsy_pT);
    g_pT_add->SetName("#Delta#it{p}_{T}(#it{p}^{ITS}) - somma dati #pi+ e #pi-");
    g_pT_add->SetTitle("#Delta#it{p}_{T}(#it{p}^{ITS})");
    g_pT_add->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    g_pT_add->GetYaxis()->SetTitle("#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (Gev/#it{c})");


    TH2F* h_deltapT_rebin= rebin(g_pT_add, kNbinsy_pT, k_ybins_pT);
    h_deltapT_rebin->SetName("h_deltapT_rebin"); 
    h_deltapT_rebin->SetTitle("#Delta#it{p}_{T}(#it{p}^{ITS}) rebinned");
    h_deltapT_rebin->GetXaxis()->SetTitle("#it{p}^{ITS} (Gev/#it{c})");
    h_deltapT_rebin->GetYaxis()->SetTitle("#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (Gev/#it{c})");



    CpT_add->cd();
    g_pT_add->Draw();

    parameterisation_study.cd("delta_pT");   
    CpT_add->Write();
    h_deltapT_rebin->Write();


 
//------------------------------------------------------------------------------ phi ------------------------------------------------------------------------------
    TH2F* phiTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/phiTpcIts");

    TCanvas *Cphi_add= new TCanvas("Cphi_add", "#phi(#it{p}^{ITS}) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    gStyle->SetOptStat(0);

    TH2F *g_phi_add = add(phiTpcIts, kNbinsy_phi);
    g_phi_add->SetName("#Delta#phi(#it{p}^{ITS}) - somma dati #pi+ e #pi-");
    g_phi_add->SetTitle("#Delta#phi(#it{p}^{ITS})");
    g_phi_add->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    g_phi_add->GetYaxis()->SetTitle("#phi^{TPC} - #phi^{ITS} (rad)");


    TH2F* h_deltaphi_rebin= rebin(g_phi_add, kNbinsy_phi, k_ybins_phi);  
    h_deltaphi_rebin->SetName("h_deltaphi_rebin");
    h_deltaphi_rebin->SetTitle("#Delta#phi(#it{p}^{ITS}) rebinned");
    h_deltaphi_rebin->GetXaxis()->SetTitle("#it{p}^{ITS} (Gev/#it{c})");
    h_deltaphi_rebin->GetYaxis()->SetTitle("#phi^{TPC} - #phi^{ITS} (rad)");
 


    Cphi_add->cd();
    g_phi_add->Draw();

    parameterisation_study.cd("delta_phi");   
    Cphi_add->Write();
    h_deltaphi_rebin->Write();
 




  
//ora devo fare un ciclo for in cui prendo proiezione sull'asse y dell'istogramma 2D, faccio due fit gaussiani e grafico dev std al variare di pT

//------------------------------------------------------------------------------ eta ------------------------------------------------------------------------------
    TH1F* hist_full_eta= new TH1F("delta_eta", "#sigma_{#Delta#eta}(#it{p}^{ITS})", kNbinsx, k_xbins);                 
    hist_full_eta->GetYaxis()->SetTitle("#sigma_{#Delta#eta}");
    hist_full_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    hist_full_eta->SetBinContent(1, 0.);
    hist_full_eta->SetBinError(1, 0.);
  
    TH1D** h_proj_eta= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status_eta= new TH1D("CovMatrixStatus_eta", "stato della matrice di covarianza in funzione di #it{p}^{ITS}", kNbinsx, k_xbins);
    h_cov_matrix_status_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_chi2_eta= new TH1D("chi2_eta", "chi2(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_chi2_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_prob_eta= new TH1D("prob_eta", "prob(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_prob_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2_eta= new TH1D("CovMatrixStatus2_eta", "stato della matrice di covarianza in funzione di #it{p}^{ITS}", kNbinsx, k_xbins);
    h_cov_matrix_status2_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_chi2_2_eta= new TH1D("chi2_2_eta", "chi2(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_chi2_2_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_prob2_eta= new TH1D("prob2_eta", "prob(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_prob2_eta->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    int entries_eta;

    TF1 *f1_eta;
    TF1 *f2_eta;

    float delta_eta, sdelta_eta;



    //------------------------------------------------------------------------------ pT ------------------------------------------------------------------------------
    TH1F* hist_full_pT= new TH1F("delta_pT", "#sigma_{#Delta#it{p}_{T}}(#it{p}^{ITS})", kNbinsx, k_xbins);                 
    hist_full_pT->GetYaxis()->SetTitle("#sigma_{#Delta#it{p}_{T}} (GeV/#it{c})");
    hist_full_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    hist_full_pT->SetBinContent(1, 0.);
    hist_full_pT->SetBinError(1, 0.);
   
    TH1D** h_proj_pT= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status_pT= new TH1D("CovMatrixStatus_pT", "stato della matrice di covarianza in funzione di #it{p}^{ITS}", kNbinsx, k_xbins);
    h_cov_matrix_status_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_chi2_pT= new TH1D("chi2_pT", "chi2(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_chi2_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_prob_pT= new TH1D("prob_pT", "prob(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_prob_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2_pT= new TH1D("CovMatrixStatus2_pT", "stato della matrice di covarianza in funzione di #it{p}^{ITS}", kNbinsx, k_xbins);
    h_cov_matrix_status2_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_chi2_2_pT= new TH1D("chi2_2_pT", "chi2(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_chi2_2_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_prob2_pT= new TH1D("prob2_pT", "prob(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_prob2_pT->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    int entries_pT;

    TF1 *f1_pT;
    TF1 *f2_pT;

    float delta_pT, sdelta_pT;



//------------------------------------------------------------------------------ phi ------------------------------------------------------------------------------
    TH1F* hist_full_phi= new TH1F("delta_phi", "#sigma_{#Delta#phi}(#it{p}^{ITS})", kNbinsx, k_xbins);                 
    hist_full_phi->GetYaxis()->SetTitle("#sigma_{#Delta#phi} (rad)");
    hist_full_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    hist_full_phi->SetBinContent(1, 0.);
    hist_full_phi->SetBinError(1, 0.);
  
    TH1D** h_proj_phi= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status_phi= new TH1D("CovMatrixStatus_phi", "stato della matrice di covarianza in funzione di #it{p}^{ITS}", kNbinsx, k_xbins);
    h_cov_matrix_status_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_chi2_phi= new TH1D("chi2_phi", "chi2(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_chi2_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_prob_phi= new TH1D("prob_phi", "prob(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_prob_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2_phi= new TH1D("CovMatrixStatus2_phi", "stato della matrice di covarianza in funzione di #it{p}^{ITS}", kNbinsx, k_xbins);
    h_cov_matrix_status2_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_chi2_2_phi= new TH1D("chi2_2_phi", "chi2(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_chi2_2_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    TH1D* h_prob2_phi= new TH1D("prob2_phi", "prob(#it{p}^{ITS})", kNbinsx, k_xbins);
    h_prob2_phi->GetXaxis()->SetTitle("#it{p}^{ITS} (GeV/#it{c})");
    int entries_phi;

    TF1 *f1_phi;
    TF1 *f2_phi;

    float delta_phi, sdelta_phi;





    for(int i=2; i<=kNbinsx; i++){

        cout<<"-----------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
        cout<<"------------------------------------------------------------- eta -------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in p: ("<<h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

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
        h_proj_eta[i-1]->SetTitle(Form("ProjectionY #Delta#eta - Intervallo in #it{p}^{ITS} = (%.1f, %.1f) Gev/#it{c}", h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)));


        parameterisation_study.cd("delta_eta/Projections_eta");
        h_proj_eta[i-1]->Write();
        }


        
        
        cout<<"------------------------------------------------------------- pT -------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in p: ("<<h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltapT_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

        h_proj_pT[i-1]=h_deltapT_rebin->ProjectionY(Form("ProjectionY_%.1f_%.1f", h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i), h_deltapT_rebin->GetXaxis()->GetBinUpEdge(i)), i, i);  
        //questa parte mi dava problemi perché, quando ci sono pochi o nessun dato, il fit, naturalmente, non viene ma mi blocca tutto quanto e mi dà cose brutte anche per ricavare lo stato del fit e della matrice di covarianza
        entries_pT=h_proj_pT[i-1]->GetEntries();
        if(entries_pT>=5){
        float nsigma_pT=3;
        if(h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i)<1.99){
        nsigma_pT=1.5;
        } else if(h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i)>3.59){
            nsigma_pT=20.;
        }
        double_fit(h_proj_pT[i-1], i, h_cov_matrix_status_pT, h_chi2_pT, h_prob_pT, h_cov_matrix_status2_pT, h_chi2_2_pT, h_prob2_pT, delta_pT, sdelta_pT, f1_pT, f2_pT, nsigma_pT);
        
        cout<<"\nDelta pT: "<<delta_pT<<" +- "<<sdelta_pT<<endl;

        hist_full_pT->SetBinContent(i, delta_pT);
        hist_full_pT->SetBinError(i, sdelta_pT);

        h_proj_pT[i-1]->GetXaxis()->SetTitle("#it{p}_{T}^{TPC} - #it{p}_{T}^{ITS} (Gev/#it{c})");
        h_proj_pT[i-1]->GetYaxis()->SetTitle("Entries");
        h_proj_pT[i-1]->SetTitle(Form("ProjectionY #Delta#it{p}_{T} - Intervallo in #it{p}^{ITS} = (%.1f, %.1f) Gev/#it{c}", h_deltapT_rebin->GetXaxis()->GetBinLowEdge(i), h_deltapT_rebin->GetXaxis()->GetBinUpEdge(i)));


        parameterisation_study.cd("delta_pT/Projections_pT");
        h_proj_pT[i-1]->Write();
        }




        cout<<"------------------------------------------------------------- phi -------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in p: ("<<h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltaphi_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

        h_proj_phi[i-1]=h_deltaphi_rebin->ProjectionY(Form("ProjectionY_%.1f_%.1f", h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaphi_rebin->GetXaxis()->GetBinUpEdge(i)), i, i);  
        //questa parte mi dava problemi perché, quando ci sono pochi o nessun dato, il fit, naturalmente, non viene ma mi blocca tutto quanto e mi dà cose brutte anche per ricavare lo stato del fit e della matrice di covarianza
        entries_phi=h_proj_phi[i-1]->GetEntries();
        if(entries_phi>=5){
        float nsigma_phi=3;
        if(h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i)<1.49){
        nsigma_phi=1.;
        } else if(h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i)<3.99){
            nsigma_phi=2.;
        }

        double_fit(h_proj_phi[i-1], i, h_cov_matrix_status_phi, h_chi2_phi, h_prob_phi, h_cov_matrix_status2_phi, h_chi2_2_phi, h_prob2_phi, delta_phi, sdelta_phi, f1_phi, f2_phi, nsigma_phi);
        
        cout<<"\nDelta phi: "<<delta_phi<<" +- "<<sdelta_phi<<endl;

        hist_full_phi->SetBinContent(i, delta_phi);
        hist_full_phi->SetBinError(i, sdelta_phi);

        h_proj_phi[i-1]->GetXaxis()->SetTitle("#phi^{TPC} - #phi^{ITS} (rad)");
        h_proj_phi[i-1]->GetYaxis()->SetTitle("Entries");
        h_proj_phi[i-1]->SetTitle(Form("ProjectionY #Delta#phi - Intervallo in #it{p}^{ITS} = (%.1f, %.1f) Gev/#it{c}", h_deltaphi_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaphi_rebin->GetXaxis()->GetBinUpEdge(i)));


        parameterisation_study.cd("delta_phi/Projections_phi");
        h_proj_phi[i-1]->Write();
        }
       

    }

//--------------------------------------------- pdf priezioni ---------------------------------------------//
    TCanvas *Cprojection_eta = new TCanvas("Cprojection_eta", "projection_eta_1.3_1.4", 500, 600);
    gStyle->SetOptFit(0);
    Cprojection_eta->cd();
    h_proj_eta[13]->SetTitle("Intervallo #it{p}^{ITS} = (1.3, 1.4) GeV/#it{c}");
    h_proj_eta[13]->GetXaxis()->SetRangeUser(-0.04, 0.04);
    h_proj_eta[13]->Draw("PE");

    TLatex t1;
    t1.SetNDC();
    t1.SetTextSize(0.071);
    t1.DrawLatex(0.63, 0.70, "#font[42]{#color[2]{#scale[1.2]{-}} #scale[0.7]{Primo fit}}");
    t1.DrawLatex(0.63, 0.65, "#font[42]{#color[4]{#scale[1.2]{-}} #scale[0.7]{Secondo fit}}");
    t1.DrawLatex(0.61, 0.55, Form("#scale[0.4]{#sigma_{#Delta#eta} = (%.2f #pm %.2f) #times 10^{-3}}", hist_full_eta->GetBinContent(13)/pow(10., -3), hist_full_eta->GetBinError(13)/pow(10., -3)));

    Cprojection_eta->Print("projection_eta_1.3_1.4.pdf");



    TCanvas *Cprojection_pT = new TCanvas("Cprojection_pT", "projection_pT_1.3_1.4", 500, 600);
    gStyle->SetOptFit(0);
    Cprojection_pT->cd();    
    h_proj_pT[13]->SetTitle("Intervallo #it{p}^{ITS} = (1.3, 1.4) GeV/#it{c}");
    h_proj_pT[13]->GetXaxis()->SetRangeUser(-0.3, 0.3);
    h_proj_pT[13]->Draw("PE");

    TLatex t2;
    t2.SetNDC();
    t2.SetTextSize(0.071);
    t2.DrawLatex(0.63, 0.70, "#font[42]{#color[2]{#scale[1.2]{-}} #scale[0.7]{Primo fit}}");
    t2.DrawLatex(0.63, 0.65, "#font[42]{#color[4]{#scale[1.2]{-}} #scale[0.7]{Secondo fit}}");
    t1.DrawLatex(0.61, 0.55, Form("#scale[0.31]{#sigma_{#Delta#it{p}_{T}} = (%.1f #pm %.1f) #times 10^{-2} GeV/#it{c}}", hist_full_pT->GetBinContent(13)/pow(10., -2), hist_full_pT->GetBinError(13)/pow(10., -2)));

    Cprojection_pT->Print("projection_pT_1.3_1.4.pdf");



    TCanvas *Cprojection_phi = new TCanvas("Cprojection_phi", "projection_phi_1.3_1.4", 500, 600);
    gStyle->SetOptFit(0);
    Cprojection_phi->cd();
    h_proj_phi[13]->SetTitle("Intervallo #it{p}^{ITS} = (1.3, 1.4) GeV/#it{c}");
    h_proj_phi[13]->GetXaxis()->SetRangeUser(-0.05, 0.05);
    h_proj_phi[13]->Draw("PE");

    TLatex t3;
    t3.SetNDC();
    t3.SetTextSize(0.071);
    t3.DrawLatex(0.63, 0.70, "#font[42]{#color[2]{#scale[1.2]{-}} #scale[0.7]{Primo fit}}");
    t3.DrawLatex(0.63, 0.65, "#font[42]{#color[4]{#scale[1.2]{-}} #scale[0.7]{Secondo fit}}");
    t1.DrawLatex(0.62, 0.55, Form("#scale[0.35]{#sigma_{#Delta#phi} = (%.1f #pm %.1f) #times 10^{-3} rad}", hist_full_phi->GetBinContent(13)/pow(10., -3), hist_full_phi->GetBinError(13)/pow(10., -3)));

    Cprojection_phi->Print("projection_phi_1.3_1.4.pdf");



    cout<<endl;
    cout<<endl;
    cout<<"-------------------------- eta - Fit con power law --------------------------"<<endl;
    TF1* power_law= new TF1("power_law", "[0]*TMath::Power(x, [1])");
    TFitResultPtr FitResult_eta= hist_full_eta->Fit(power_law, "RMSI", "0", 0.3, 4);    //R: in the specified range    I: Uses the integral of function in the bin instead of the default bin center value.
    cout<<"\nCHi2: "<<FitResult_eta->Chi2()<<endl;
    cout<<"Prob: "<<FitResult_eta->Prob()<<endl;


    TCanvas *plot_eta= new TCanvas("Cplot_eta", "plot_eta", 560,440);
    hist_full_eta->SetStats(1);   //activates the statistics box
    //hist_full_eta->GetYaxis()->SetTitleSize(0.04);
    //hist_full_eta->GetYaxis()->SetLabelSize(0.02);
    gStyle->SetOptStat(1100);   //I: mean, II: dev std, III: entries, IV: title of the statbox
    gStyle->SetOptFit(1111);

    //plot_eta->Update();
    plot_eta->cd();
    hist_full_eta->GetXaxis()->SetRangeUser(0., 4.);
    hist_full_eta->Draw();
    power_law->SetRange(0., 4.);
    power_law->Draw("SAME");

    TLine *l1_eta= new TLine(0.3, 0, 0.3, 0.035);
    //TLine *l2_eta= new TLine(4., 0, 4., 0.0245);
    l1_eta->SetLineStyle(3);
    //l2_eta->SetLineStyle(3);
    l1_eta->Draw();
    //l2_eta->Draw();

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
    TFitResultPtr FitResult_pT= hist_full_pT->Fit(parabola, "RMSI", "0", 0.3, 4);
    cout<<"\nCHi2: "<<FitResult_pT->Chi2()<<endl;
    cout<<"Prob: "<<FitResult_pT->Prob()<<endl;


    TCanvas *plot_pT= new TCanvas("Cplot_pT", "plot_pT", 500,440);
    hist_full_pT->SetStats(1);   //activates the statistics box

    //plot_pT->Update();
    plot_pT->cd();
    hist_full_pT->GetXaxis()->SetRangeUser(0., 4.);
    hist_full_pT->Draw();
    gStyle->SetOptFit(1111);
    parabola->SetRange(0., 4.);
    parabola->Draw("SAME");
    auto ps = (TPaveStats *)hist_full_pT->GetListOfFunctions()->FindObject("stats");
    ps->SetY1NDC(0.15);
    ps->SetY2NDC(0.45);
    plot_pT->Modified();
    plot_pT->Update();

    TLine *l1_pT= new TLine(0.3, 0, 0.3, 0.39);
    //TLine *l2_pT= new TLine(4., 0.165, 4., 0.39);
    l1_pT->SetLineStyle(3);
    //l2_pT->SetLineStyle(3);
    l1_pT->Draw();
    //l2_pT->Draw();

    TLatex t_pT;
    t_pT.SetNDC();
    t_pT.SetTextSize(0.071);
    t_pT.DrawLatex(0.20, 0.75, "#scale[0.9]{#it{f}(x) = p0 + p1 x + p2 x^{2}}");


    parameterisation_study.cd("delta_pT");
    plot_pT->Write();
    plot_pT->Print("plot_pT.pdf");


 




    cout<<endl;
    cout<<"-------------------------- phi - Fit con esponenziale negativo --------------------------"<<endl;
    TF1* exp_neg= new TF1("exp_neg", "[0]+[1]*exp(-([2]*x))");
    exp_neg->SetParameter(0, 7.21022e-03);
    exp_neg->SetParameter(1, 5.62625e-02);
    exp_neg->SetParameter(2, 2.97602e+00);
    TFitResultPtr FitResult_phi= hist_full_phi->Fit(exp_neg, "RMSI", "", 0.3, 4);
    cout<<"\nCHi2: "<<FitResult_phi->Chi2()<<endl;
    cout<<"Prob: "<<FitResult_phi->Prob()<<endl;
 

    TCanvas *plot_phi= new TCanvas("Cplot_phi", "plot_phi", 500,440);
    hist_full_phi->SetStats(1);   //activates the statistics box

    plot_phi->Update();
    hist_full_phi->GetXaxis()->SetRangeUser(0., 4.);
    hist_full_phi->GetYaxis()->SetRangeUser(0., 0.065);
    hist_full_phi->Draw();
    plot_phi->cd();
    gStyle->SetOptFit(1111);
    exp_neg->SetRange(0., 4.);
    exp_neg->Draw("SAME");

    TLine *l1_phi= new TLine(0.3, 0, 0.3, 0.065);
    //TLine *l2_phi= new TLine(4., 0, 4., 0.0405);
    l1_phi->SetLineStyle(3);
    //l2_phi->SetLineStyle(3);
    l1_phi->Draw();
    //l2_phi->Draw();

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







    //--------------------------------------------- timeTpcIts ---------------------------------------------//
    TH1F *timeTpcIts= (TH1F*)fAnalysis.Get("efficiency-q-a/timeTpcItsNoNorm");
    TCanvas *CtimeTpcIts= new TCanvas ("CtimeTpcIts", "timeTpcIts", 700, 600); 
    gStyle->SetOptStat(0);
    //CtimeTpcIts->SetLogy();
    timeTpcIts->GetXaxis()->SetRangeUser(-70e3, 70e3);
    CtimeTpcIts->cd();
    timeTpcIts->Draw();

    parameterisation_study.cd();
    CtimeTpcIts->Write();
    CtimeTpcIts->Print("timeTpcIts.pdf");

 
    parameterisation_study.Close();
    
}