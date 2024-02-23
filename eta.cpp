#include "Config.h"

using namespace std;


void estremi(TH1D *hist, float a, float b, float &m, float &M)
{
    float aa= hist->GetXaxis()->FindBin(a);
    m= hist->GetXaxis()->GetBinLowEdge(aa);

    float bb= hist->GetXaxis()->FindBin(b);
    M= hist->GetXaxis()->GetBinUpEdge(bb);
}




void eta(const double nsigma=1.){

    TFile fAnalysis("AnalysisResults-mc.root");
    TFile delta_eta_study("delta_eta_study.root", "RECREATE");


    cout<<"------------------------------------- studio su delta eta --------------------------------------------"<<endl;
    cout<<endl;

    TH2F* etaTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/etaTpcIts");



    //prendo l'istogramma 2D e cerco il binning del suo asse x e del suo asse y-> mi servono per sommare la parte positiva e quella negativa

    int nbinx= etaTpcIts->GetXaxis()->GetLast();
    float width_x= etaTpcIts->GetXaxis()->GetBinWidth(1);
    float *bin_ledgex = new float[nbinx];
    float minx;
    int nx;     
    for(int i=1; i<=nbinx; i++){
        bin_ledgex[i]= etaTpcIts->GetXaxis()->GetBinLowEdge(i);
        if( abs(bin_ledgex[i]-0.)<width_x/2){
            minx=i;
            cout<<"il bin minimo è: "<<minx<<endl;
            nx= minx-1;
            break;
        }                               
    }

    cout<<"nx: "<<nx<<endl;

    float max_x= etaTpcIts->GetXaxis()->GetBinUpEdge(nbinx);





     
    int nbiny = etaTpcIts->GetYaxis()->GetLast();
    float *bin_ledgey = new float[nbiny];
    float miny;
    float width_y= etaTpcIts->GetYaxis()->GetBinWidth(1);
    float max_y= etaTpcIts->GetYaxis()->GetBinUpEdge(nbiny);
    float min_y= etaTpcIts->GetYaxis()->GetBinLowEdge(1); 



    // prendo l'istogramma 2D solo per pt positivi e sommo i bin speculari rispetto allo zero-> quello che ottengo è un TH2F i cui conteggi contengono sia pi+ che pi-

    TH2F* eta_add= new TH2F("eta_add", "#eta(#it{p}_{T}) - somma dati #pi+ e #pi-", nx, 0., max_x, nbiny, min_y, max_y);
    eta_add->SetName("eta - somma dati #pi+ e #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=nbiny; j++){
            eta_add->SetBinContent(i, j, etaTpcIts->GetBinContent(nx+i, j));
        }
    }



    TH2F* eta_pi_neg= new TH2F("eta_pi-", "#eta(#it{p}_{T}) - #pi-", nx, 0., max_x, nbiny, min_y, max_y);
    eta_pi_neg->SetName("eta - #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=nbiny; j++){
            eta_pi_neg->SetBinContent(i, j, etaTpcIts->GetBinContent(nx+1-i, j));
        }
    }


    
  
    TCanvas *Ceta_add= new TCanvas("Ceta_add", "eta(p_T) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    TH2F *g_eta_add = new TH2F(*eta_add);   //ci metto sopra un grafico   sto passando l'oggetto puntato, ne sto facendo una copia
    gStyle->SetOptStat(0);

    g_eta_add->SetName("#eta(#it{p}_{T}) - somma dati #pi+ e #pi-");
    g_eta_add->Add(eta_pi_neg);
    g_eta_add->SetTitle("#eta(#it{p}_{T})");
    g_eta_add->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    g_eta_add->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");


    double* ybins= new double[nbiny+1];
    for(int i=1; i<=nbiny; i++){
        ybins[i-1]= g_eta_add->GetYaxis()->GetBinLowEdge(i);
    }

    ybins[nbiny]=g_eta_add->GetYaxis()->GetBinUpEdge(nbiny);


  
  //sto riempendo il mio istogramma rebinnato
    TH2F* h_deltaeta_rebin= new TH2F("h_deltaeta_rebin", "#eta(#it{p_{T}}) rebinned", kNbinsx, k_xbins, nbiny, ybins);    //istogramma g_eta_add rebinnato
    for(int i=1; i<=kNbinsx; i++){     //loop su x nuovo
        int xmin= g_eta_add->GetXaxis()->FindBin(k_xbins[i-1]+epsilon);
        int xmax= g_eta_add->GetXaxis()->FindBin(k_xbins[i]-epsilon);
        printf("xmin = %d, xmax = %d, k_xbins[j-1] = %f, k_xbins[j] = %f\n", xmin, xmax, k_xbins[i-1], k_xbins[i]);
        for(int j=1; j<=nbiny; j++){   //loop su y nuovo
            double integral= g_eta_add->Integral(xmin, xmax, j, j);
            h_deltaeta_rebin->SetBinContent(i, j, integral);
        }
    } 

    h_deltaeta_rebin->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");
    h_deltaeta_rebin->GetXaxis()->SetTitle("#it{p_{T}} (Gev/#it{c})");









    

    Ceta_add->cd();
    g_eta_add->Draw();

    delta_eta_study.cd();   
    Ceta_add->Write();
    h_deltaeta_rebin->Write();
    //Ceta_add->Print("eta_add.pdf");




//ora devo fare un ciclo for in cui prendo proiezione sull'asse y dell'istogramma 2D, faccio due fit gaussiani e grafico dev std al variare di pT
    TH1F* hist_full= new TH1F("delta_eta", "#sigma_{#Delta#eta}(#it{p_{T}})", kNbinsx, k_xbins);                 
    hist_full->GetYaxis()->SetTitle("#sigma_{#Delta#eta}");
    hist_full->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    
    TH1D** h_proj= new TH1D*[kNbinsx];
    TH1D* h_cov_matrix_status= new TH1D("CovMatrixStatus", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2= new TH1D("chi2", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob= new TH1D("prob", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");

    TH1D* h_cov_matrix_status2= new TH1D("CovMatrixStatus2", "stato della matrice di covarianza in funzione di #it{p}_{T}", kNbinsx, k_xbins);
    h_cov_matrix_status2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_chi2_2= new TH1D("chi2_2", "chi2(#it{p}_{T})", kNbinsx, k_xbins);
    h_chi2_2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    TH1D* h_prob2= new TH1D("prob2", "prob(#it{p}_{T})", kNbinsx, k_xbins);
    h_prob2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    int entries;




    for(int i=1; i<=kNbinsx; i++){

        
        float mean, dev_std, mean2, dev_std2;
        float low_edge, up_edge, low_edge2, up_edge2, low_edge3, up_edge3;
        float delta, sdelta;


        cout<<"-----------------------------------------------------------------------------------------"<<endl;
        cout<<"\nNUMERO: "<<i<<"              Intervallo in pt: ("<<h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i)<<","<< h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)<<")"<<endl;

        h_proj[i]=h_deltaeta_rebin->ProjectionY(Form("ProjectionY_%.1f_%.1f", h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)), i, i);   //prendo la proiezione i   Attenzione a qual è l'istogramma che voglio utilizzare
        //questa parte mi dava problemi perché, quando ci sono pochi o nessun dato, il fit, naturalmente, non viene ma mi blocca tutto quanto e mi dà cose brutte anche per ricavare lo stato del fit e della matrice di covarianza
        entries=h_proj[i]->GetEntries();
        cout<<"\nEntries: "<<entries<<endl;
        if(entries<5) continue;
        
        
        float estr_inf, estr_sup;
        mean= h_proj[i]->GetMean(1);   //int axis=1 mi dà il valore medio sull'asse X (2<->Y,  3<->Z)
        dev_std= h_proj[i]->GetStdDev(1);
        low_edge= mean-dev_std;
        up_edge= mean+dev_std;
        estremi(h_proj[i], low_edge, up_edge, estr_inf, estr_sup);
        cout<<endl;
        cout<<"\nlow edge: "<<low_edge<<endl;
        cout<<"up edge: "<<up_edge<<endl;
        cout<<"estr_inf: "<<estr_inf<<endl;
        cout<<"estr_sup: "<<estr_sup<<endl;
        cout<<endl;

        
        TFitResultPtr FitResult=h_proj[i]->Fit("gaus", "RMLS", "", estr_inf, estr_sup);    //specifico il range su cui effettuare il fit
        //cout<<"\nStato del fit: "<<FitResult->Status()<<endl;
        cout<<"\nStato della matrice di covarianza: "<<FitResult->CovMatrixStatus()<<endl;
        cout<<"Chi2: "<<FitResult->Chi2()<<endl;
        cout<<"Probabilità: "<<FitResult->Prob()<<endl;

        //grafico stato del fit, stato della matrice di covarianza, chi2 e probabilità
        h_cov_matrix_status->SetBinContent(i, FitResult->CovMatrixStatus());
        h_chi2->SetBinContent(i, FitResult->Chi2());
        h_prob->SetBinContent(i, FitResult->Prob());
       
         


        h_proj[i]->GetFunction("gaus")->SetName("gaus1");
        TF1 *f1 = (TF1*)h_proj[i]->GetListOfFunctions()->FindObject("gaus1");
 
        cout<<"\n-------------- Secondo fit --------------"<<endl;
        float estr_inf2, estr_sup2;
        mean2=h_proj[i]->GetFunction("gaus1")->GetParameter(1);     //f(x) = p0*exp(-0.5*((x-p1)/p2)^2)  la gaussiana ha tre parametri: a me interessa p2, dev std
        dev_std2=h_proj[i]->GetFunction("gaus1")->GetParameter(2); 
        low_edge2= mean2-dev_std2;
        up_edge2= mean2+dev_std2;
        estremi(h_proj[i], low_edge2, up_edge2, estr_inf2, estr_sup2);
        cout<<"\nlow edge2: "<<low_edge2<<endl;
        cout<<"up edge2: "<<up_edge2<<endl;
        cout<<"estr_inf2: "<<estr_inf2<<endl;
        cout<<"estr_sup2: "<<estr_sup2<<endl;
        cout<<endl;

        TFitResultPtr FitResult2 = h_proj[i]->Fit("gaus", "RMLS+", "",estr_inf2, estr_sup2);     //+ mi fitta questa gaussiana sul "grafico" in cui c'è già l'altra e non la sovrascrive
        cout<<"\nStato della matrice di covarianza: "<<FitResult2->CovMatrixStatus()<<endl;
        cout<<"Chi2: "<<FitResult2->Chi2()<<endl;
        cout<<"Probabilità: "<<FitResult2->Prob()<<endl;

        h_cov_matrix_status2->SetBinContent(i, FitResult2->CovMatrixStatus());
        h_chi2_2->SetBinContent(i, FitResult2->Chi2());
        h_prob2->SetBinContent(i, FitResult2->Prob());



        h_proj[i]->GetFunction("gaus")->SetName("gaus2");
        delta= h_proj[i]->GetFunction("gaus2")->GetParameter(2);
        sdelta= h_proj[i]->GetFunction("gaus2")->GetParError(2);
        cout<<"\nDelta eta: "<<delta<<" +- "<<sdelta<<endl;

        hist_full->SetBinContent(i, delta);
        hist_full->SetBinError(i, sdelta);



        TF1 *f2 = (TF1*)h_proj[i]->GetListOfFunctions()->FindObject("gaus2");
        f2->SetLineColor(4);


       
        h_proj[i]->GetXaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");
        h_proj[i]->GetYaxis()->SetTitle("Entries");
        h_proj[i]->SetTitle(Form("ProjectionY - Intervallo in #it{p_{T}} = (%.1f, %.1f) Gev/#it{c}", h_deltaeta_rebin->GetXaxis()->GetBinLowEdge(i), h_deltaeta_rebin->GetXaxis()->GetBinUpEdge(i)));


        if(i==14){
        TCanvas* Cproiezione= new TCanvas("Cproiezione", "proiezione - (1.3, 1.4) Gev/#it{c}", 200,50,600,400); 
        gStyle->SetOptStat(0);

        Cproiezione->cd(0);
        h_proj[i]->GetXaxis()->SetRangeUser(estr_inf, estr_sup);
        h_proj[i]->Draw();

        TLegend *leg= new TLegend(0.73,0.85,0.87,0.65);
        leg->SetHeader("Fit gaussiani", "C");   //C -> viene messo al centro
        leg->AddEntry(f1, "primo fit","l");
        leg->AddEntry(f2, "secondo fit","l");
        leg->SetTextSize(0.038);
        leg->Draw();




        TCanvas* Cproiezione1= new TCanvas("Cproiezione1", "proiezione - (0.2, 0.3) Gev/#it{c}", 200,50,600,400); 
        gStyle->SetOptStat(0);

        Cproiezione1->cd(0);
        h_proj[i-11]->GetXaxis()->SetRangeUser(-0.12, 0.12);
        h_proj[i-11]->Draw();

        TLegend *leg1= new TLegend(0.73,0.85,0.87,0.65);
        leg1->SetHeader("Fit gaussiani", "C");   //C -> viene messo al centro
        leg1->AddEntry(f1, "primo fit","l");
        leg1->AddEntry(f2, "secondo fit","l");
        leg1->SetTextSize(0.038);
        leg1->Draw();


        delta_eta_study.cd();
        Cproiezione->Print("proiezione_1.3_1.4.pdf");
        Cproiezione1->Print("proiezione_0.2_0.3.pdf");
       }




        delta_eta_study.cd();
        h_proj[i]->Write();

       

    }


    cout<<endl;
    TF1* power_law= new TF1("power_law", "[0]*TMath::Power(x, [1])",k_xbins[2], k_xbins[kNbinsx]);
    TFitResultPtr FitResult3= hist_full->Fit(power_law, "RMS");
    cout<<"\n CHi2: "<<FitResult3->Chi2()<<endl;
    cout<<"Prob: "<<FitResult3->Prob()<<endl;

    TCanvas* Cdelta_eta= new TCanvas("Cdelta_eta", "#sigma_{#Delta#eta}(#it{p_{T}})", 200,50,600,400); 
    gStyle->SetOptStat(0);

    Cdelta_eta->cd(0); 
    hist_full->GetYaxis()->SetRangeUser(0., 0.06);
    hist_full->SetTitle("#sigma_{#Delta#eta}(#it{p_{T}})");
    hist_full->GetYaxis()->SetTitle("#sigma_{#Delta#eta}(#it{p_{T}})");
    hist_full->Draw();
    

    


    delta_eta_study.cd();
    h_cov_matrix_status->Write();
    h_chi2->Write();
    h_prob->Write();
    h_cov_matrix_status2->Write();
    h_chi2_2->Write();
    h_prob2->Write();
    hist_full->Write();
    Cdelta_eta->Print("delta_eta.pdf");
    















































    cout<<"------------------------------------------------------------- studio su pT -------------------------------------------"<<endl;
    cout<<endl;

    TH2F* ptTpcIts= (TH2F*)fAnalysis.Get("efficiency-q-a/ptTpcIts");



    //prendo l'istogramma 2D e cerco il binning del suo asse x e del suo asse y-> mi servono per sommare la parte positiva e quella negativa

    int nbinx= etaTpcIts->GetXaxis()->GetLast();
    float width_x= etaTpcIts->GetXaxis()->GetBinWidth(1);
    float *bin_ledgex = new float[nbinx];
    float minx;
    int nx;     
    for(int i=1; i<=nbinx; i++){
        bin_ledgex[i]= etaTpcIts->GetXaxis()->GetBinLowEdge(i);
        if( abs(bin_ledgex[i]-0.)<width_x/2){
            minx=i;
            cout<<"il bin minimo è: "<<minx<<endl;
            nx= minx-1;
            break;
        }                               
    }

    cout<<"nx: "<<nx<<endl;

    float max_x= etaTpcIts->GetXaxis()->GetBinUpEdge(nbinx);





     
    int nbiny = etaTpcIts->GetYaxis()->GetLast();
    float *bin_ledgey = new float[nbiny];
    float miny;
    float width_y= etaTpcIts->GetYaxis()->GetBinWidth(1);
    float max_y= etaTpcIts->GetYaxis()->GetBinUpEdge(nbiny);
    float min_y= etaTpcIts->GetYaxis()->GetBinLowEdge(1); 



    // prendo l'istogramma 2D solo per pt positivi e sommo i bin speculari rispetto allo zero-> quello che ottengo è un TH2F i cui conteggi contengono sia pi+ che pi-

    TH2F* eta_add= new TH2F("eta_add", "#eta(#it{p}_{T}) - somma dati #pi+ e #pi-", nx, 0., max_x, nbiny, min_y, max_y);
    eta_add->SetName("eta - somma dati #pi+ e #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=nbiny; j++){
            eta_add->SetBinContent(i, j, etaTpcIts->GetBinContent(nx+i, j));
        }
    }



    TH2F* eta_pi_neg= new TH2F("eta_pi-", "#eta(#it{p}_{T}) - #pi-", nx, 0., max_x, nbiny, min_y, max_y);
    eta_pi_neg->SetName("eta - #pi-");

    for(int i=1; i<=nx; i++){
        for(int j=1; j<=nbiny; j++){
            eta_pi_neg->SetBinContent(i, j, etaTpcIts->GetBinContent(nx+1-i, j));
        }
    }


    
  
    TCanvas *Ceta_add= new TCanvas("Ceta_add", "eta(p_T) - somma dati pi+ e pi-", 200,50,600,400);   //disegno il canvas
    TH2F *g_eta_add = new TH2F(*eta_add);   //ci metto sopra un grafico   sto passando l'oggetto puntato, ne sto facendo una copia
    gStyle->SetOptStat(0);

    g_eta_add->SetName("#eta(#it{p}_{T}) - somma dati #pi+ e #pi-");
    g_eta_add->Add(eta_pi_neg);
    g_eta_add->SetTitle("#eta(#it{p}_{T})");
    g_eta_add->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    g_eta_add->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");


    double* ybins= new double[nbiny+1];
    for(int i=1; i<=nbiny; i++){
        ybins[i-1]= g_eta_add->GetYaxis()->GetBinLowEdge(i);
    }

    ybins[nbiny]=g_eta_add->GetYaxis()->GetBinUpEdge(nbiny);


  
  //sto riempendo il mio istogramma rebinnato
    TH2F* h_deltaeta_rebin= new TH2F("h_deltaeta_rebin", "#eta(#it{p_{T}}) rebinned", kNbinsx, k_xbins, nbiny, ybins);    //istogramma g_eta_add rebinnato
    for(int i=1; i<=kNbinsx; i++){     //loop su x nuovo
        int xmin= g_eta_add->GetXaxis()->FindBin(k_xbins[i-1]+epsilon);
        int xmax= g_eta_add->GetXaxis()->FindBin(k_xbins[i]-epsilon);
        printf("xmin = %d, xmax = %d, k_xbins[j-1] = %f, k_xbins[j] = %f\n", xmin, xmax, k_xbins[i-1], k_xbins[i]);
        for(int j=1; j<=nbiny; j++){   //loop su y nuovo
            double integral= g_eta_add->Integral(xmin, xmax, j, j);
            h_deltaeta_rebin->SetBinContent(i, j, integral);
        }
    } 

    h_deltaeta_rebin->GetYaxis()->SetTitle("#eta^{TPC} - #eta^{ITS}");
    h_deltaeta_rebin->GetXaxis()->SetTitle("#it{p_{T}} (Gev/#it{c})");









    

    Ceta_add->cd();
    g_eta_add->Draw();

    delta_eta_study.cd();   
    Ceta_add->Write();
    h_deltaeta_rebin->Write();
    //Ceta_add->Print("eta_add.pdf");

}  