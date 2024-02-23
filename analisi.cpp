/* 
#include <TH1F.h>
#include <TH2F.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TF1.h> 
root interpreta riga per riga il codice e si prende le cose quando gli servono quindi non c'è bisogno che io metta gli include. Gli include servono se devo compilare il codice
*/

#include "Config.h"



//ho dovuto scriverla fuori perchè C does not allow nested function definitions. Move these functions out of main 
    Double_t myfunction(Double_t *x, Double_t *par)
    {
        Float_t xx =x[0];
        Double_t f;
        if(xx<par[2] || xx>par[3])       //par[2] è il limite inferiore dell'intervallo e par[3] è il limite superiore
            f = exp(par[0]+par[1]*xx);   //exp è il numero di nepero
        else
            f=0;
        return f;
    }





using namespace std;


void analisi(){
    TFile fmass("AnalysisResults.root");                      //Opens or creates a local ROOT file. ("nome del file")
    TFile plotTesi("plotTesi.root", "recreate");
    TH1F *MassV0=(TH1F*)fmass.Get("efficiency-q-a/massV0");   //lo prende dalla cartella efficiency-q-a e il nome del file è massV0
    TH1F *MassV0_exp= new TH1F(*MassV0);                      //copia usata per il fit con solo esponenziale
    TH1F *MassV0_plot= new TH1F(*MassV0);                     //copia usata per il plot
    

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

//------------------------------------------- Fit gaussiana+esponenziale --------------------------------------------------------------//


    //ora che sono risalita ai dati (che ho immagazzinati ma che non vedo) dall'immagine, creo un canvas, un ambiente in cui disegnare l'istogramma e poi lo disegno per verificare che sia andato tutto bene
    TCanvas *CMassV0 = new TCanvas("CMassV0 ","MassV0",200,50,600,400); 
     gStyle->SetOptStat(0);    //mi toglie il riquadro con mean, dev std, etc. Nel file .root si vede ancora ma una volta salvato in pdf o altro non si vede più
     

    //ora devo fittare con una gaussiana più un esponenziale decrescente
    TF1 *gauss_exp = new TF1("gauss_exp","gausn(0)+expo(3)",0.,1.);  //gausn: gaussiana normalizzata

    //sono più importanti i SetParLimits dei SetParameters
    gauss_exp->SetParLimits(0,0,kMaxCounts);   //ho una decina di bin (forse un po' di meno) con 100 conteggi-> il parametro [0] è il parametro di normalizzazione, è dell'ordine del numero di conteggi. È un valore positivo
    gauss_exp->SetParLimits(1,0.47,0.52);   //è la media quindi so che deve essere vicino a 490 MeV/c2 ma il grafico è in GeV
    gauss_exp->SetParLimits(2,0,0.1);   //è la dev std. Per stimare la dev std devo prendere la metà larghezza a metà altezza. È dell'ordine di 0.01. È un valore positivo
    gauss_exp->SetParLimits(3,0,10);   //l'esponensiale è definito come exp^(a+bx)  a =[0,10] circa
    gauss_exp->SetParLimits(4,-10,0);     // b deve essere negativo -> [-10,0]

    gauss_exp->SetLineColor(2); 


    auto fitResult=MassV0->Fit(gauss_exp,"RMLS");   //L: likelihood  -> usa il principio di massima verosomiglianza   S: get fit result pointer teoricamente
    CMassV0->cd(0);

    MassV0->DrawCopy();

    cout<<"Probability: "<<gauss_exp->GetProb()<<endl;


      

    float mu = gauss_exp->GetParameter(1);
    float sigma = gauss_exp->GetParameter(2);
    float m = mu - kMassWindow * sigma;
    float M = mu + kMassWindow * sigma;
    int nbin= MassV0->GetNbinsX();     //funzionalità dell'istogramma quindi devo mettere nome dell'istogramma
    const int n=400;                   //valore messo alla luce di quanto mi viene restituito dalla riga predecente
    cout<<endl;
    cout<<"numero di bin su x: "<<nbin<<endl;
    cout<<endl;


    cout<<"min: "<<m<<" GeV/c2"<<endl;
    cout<<"max: "<<M<<" GeV/c2"<<endl;
    cout<<endl;
    float minimo, massimo;
    int bin_min, bin_max;
    float limsx[n];
    for(int i=1;i<=nbin;i++){                                 //il binning inizia da 1
        limsx[i]=MassV0->GetXaxis()->GetBinLowEdge(i);        //nome del TH1F->GetXaxis()->GetBinLowEdge(i) mi crea un array con tutti i limiti inferiori dei bin
        if(limsx[i]>m){
            minimo=limsx[i];
            bin_min=i;
            cout<<"Questo è il limite inferiore dell'intervallo di integrazione [GeV/c2]: "<<minimo<<endl;   //se do più di un comando nell' if devo mettere {}
            cout<<"Il numero del bin contenente il minimo è: "<<bin_min<<endl;
            break;
        }
    }

    cout<<endl;
   cout<<endl;
    float limdx[n];   //l'argomento dell'array deve essere una static const int ma poi nel for la static const int non funziona. Devo mettere l'altra cosa (solo int)
    for(int i=1;i<=nbin;i++){
        limdx[i]=MassV0->GetXaxis()->GetBinUpEdge(i);    
        if(limdx[i]>M){
            massimo=limdx[i-1];
            bin_max=i-1;
            cout<<"Questo è il limite superiore dell'intervallo di integrazione [GeV/c2]: "<<massimo<<endl;   //se do più di un comando nell' if devo mettere {}
            cout<<"Il numero del bin contenente il massimo è: "<<bin_max<<endl;
            break;
        }
    }


    cout<<"\n"<<endl;
    float width= MassV0->GetBinWidth(1);  //il binning parte da 1

    float integrale=0;
    for(int i=bin_min; i<=bin_max; i++){
        integrale += MassV0->GetBinContent(i);
    }

        cout<<"Il valore dell'integrale totale è: "<<integrale<<endl;







//----------------------------- Fit con solo esponenziale sulle side bands --------------------------------------------//

cout<<"\n----------------------------- Fit con solo esponenziale sulle side bands --------------------------------------------\n"<<endl;
MassV0_exp->SetName("fit_exp_MassV0");
TCanvas *CMassV0_exp = new TCanvas("CMassV0_exp ","MassV0_exp",200,50,600,400);

    float e_inferiore= mu - kSBandsWindow * sigma;
    float e_superiore= mu + kSBandsWindow * sigma;


    float min, max;
    float exp_bin_min, exp_bin_max;
    float LowEdge[n];
    for(int i=1;i<=nbin;i++){                                 
        LowEdge[i]=MassV0->GetXaxis()->GetBinLowEdge(i);    
        if(LowEdge[i]>e_inferiore){
            min=LowEdge[i];
            exp_bin_min=i;
            cout<<"Questo è il limite inferiore dell'intervallo di integrazione [GeV/c2]: "<<min<<endl;   
            cout<<"Il numero del bin contenente il minimo è: "<<exp_bin_min<<endl;
            break;
        }
    }

    cout<<endl;
   cout<<endl;
    float UpEdge[n];   
    for(int i=1;i<=nbin;i++){
        UpEdge[i]=MassV0->GetXaxis()->GetBinUpEdge(i);    
        if(UpEdge[i]>e_superiore){
            max=UpEdge[i-1];
            exp_bin_max=i-1;
            cout<<"Questo è il limite superiore dell'intervallo di integrazione [GeV/c2]: "<<max<<endl;   
            cout<<"Il numero del bin contenente il massimo è: "<<exp_bin_max<<endl;
            break;
        }
    }
    cout<<endl;







 for(int i=exp_bin_min; i<=exp_bin_max; i++){
        MassV0_exp->SetBinContent(i, 0);   //SetBinContent(nbin, content)
    }



   auto f1 = new TF1("myfunc",myfunction,0,1,4);   //myfunction è la funzione che ho definito ha 4 parametri ma i parametri 2 e 3 sono fissati
   f1->FixParameter(2, min);
   f1->FixParameter(3, max);
   f1->SetNpx(1000);
   f1->SetParLimits(0,0,10);   //l'esponensiale è definito come exp^(a+bx)  a =[0,10] circa
   f1->SetParLimits(1,-10,0);     // b deve essere negativo -> [-10,0]
   CMassV0_exp->cd(0);
   TFitResultPtr r = MassV0_exp->Fit(f1, "RMLS"); 
   MassV0_exp->Draw();

  cout<<"Probability: "<<f1->GetProb()<<endl;


    plotTesi.cd();    //per come ho aperto il file, questo non è un puntatore quindi metto il punto
    CMassV0_exp->Write();
    MassV0_exp->Write();
    CMassV0_exp->Print("exponential_fit.pdf");






cout<<endl;
cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
cout<<endl;
cout<<"Ora verifico compatibilità tra parametri dell'esponenziale decrescente ottenuti dal fit con funzione gausn(0)+expo(3) e quelli ottenuti con fit sul solo fondo (vedi sezione precedenete o macro cont_bin_analysis.cpp)"<<endl;

float a= gauss_exp->GetParameter(3);
float b= gauss_exp->GetParameter(4);

float p0= f1->GetParameter(0);
float sp0= f1->GetParError(0);
float p1= f1->GetParameter(1);
float sp1= f1->GetParError(1);


cout<<"\na= "<<a<<endl;
cout<<"p0= "<<p0<<" +- "<<sp0<<endl;
cout<<"b= "<<b<<endl;
cout<<"p1= "<<p1<<" +- "<<sp1<<endl;


float Za= (a-p0)/sp0;
float Zb= (b-p1)/sp1;
cout<<"\nZa= "<<Za<<endl;
cout<<"Zb= "<<Zb<<endl;
cout<<endl; 


if(abs(Za)>1.96 || abs(Zb)>1.96)
    cout<<"ATTENZIONE!!! I parametri calcolati nei due modi non sono compatibili. Verificare se sono stati commessi errori o valutare quale set di parametri utilizzare di seguito"<<endl;
    


if(abs(Za)<1.96 && abs(Zb)<1.96) 
    cout<<"I parametri calcolati nei due modi sono compatibili. D'ora in poi utilizzaremo p0 e p1 perchè ha più senso usare parametri di un fit effettuato senza il picco. \nOra calcolo le quantità di cui ho bisogno utilizzando p0 e p1"<<endl;
cout<<endl;



    TF1 *esponenziale = new TF1("esponenziale","exp([0]+[1]*x)",0,1);
    esponenziale->SetLineColor(2);  
    esponenziale->SetParameters(p0,p1);


    float integrale_exp;
    integrale_exp= esponenziale->Integral(minimo, massimo)/width;
    float s_integrale_exp= esponenziale->IntegralError(minimo, massimo, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/width;

    cout<<"\n Il valore dell'integrale del solo esponenziale decrescente è: "<<integrale_exp<<"+-"<<s_integrale_exp<<endl;



    //ora definisco le quantità che devo calcolare
    float S,B, purezza, significatività;
    S= integrale-integrale_exp;   //segnale
    B= integrale_exp;             //background(fondo)
    purezza= S/(S+B);
    significatività= S/sqrt(S+B); //è l'inverso dell'incertezza relativa sul segnale 

    float sS,sB, s_purezza, s_significatività;
    sS= sqrt(integrale + pow(s_integrale_exp, 2));  //l'errore sul numero totale di conteggi è la radice dei conteggi quindi s_integrale^2=integrale
    sB= s_integrale_exp;
    s_purezza= sqrt( pow(B*sS/pow(S+B,2),2) + pow(S*sB/pow(S+B,2),2) );  
    s_significatività= sqrt( pow((S/2+B)*sS/pow(S+B,3/2),2) + pow((S/2)*sB/pow(S+B,3/2),2) );   



 cout<<"\n\nS= "<<S<<" +- "<<sS<<endl;
 cout<<"B= "<<B<<" +- "<<sB<<endl;
cout<<"purezza= "<<purezza<<" +- "<<s_purezza<<endl;
cout<<"significatività= "<<significatività<<" +- "<<s_significatività<<endl;
cout<<endl;






cout<<"\n----------------------------------------------------------------------------------------------- Plot -----------------------------------------------------------------------------------------------\n"<<endl;

    MassV0_plot->SetName("MassV0_plot");

    TCanvas *CMassV0plot = new TCanvas("CMassV0plot ","MassV0_plot",200,50,600,400); 
     gStyle->SetOptStat(0);



    CMassV0plot->cd(0);  //disegna sul canvas
    //CMassV0plot->SetLogy();
    //CMassV0plot->Draw();  
    MassV0_plot->Draw();   //disegna l'istogramma

    esponenziale->Draw("lSAME");  //superimpose on top of existing picture



float MaxBin= MassV0_plot->GetMaximumBin();   //bin con il numero massimo di entries
float yf= 0.9*(MassV0_plot->GetBinContent(MaxBin));

    TLine *l1= new TLine(minimo, 0, minimo, yf);
    TLine *l2= new TLine(massimo, 0, massimo, yf);

    l1->SetLineStyle(3);  //anche 7 (e 5) non era male
    l2->SetLineStyle(3);

    l1->Draw();
    l2->Draw();



   TLatex t;
   t.SetNDC();
   t.SetTextSize(0.031);
   t.DrawLatex(0.59, 0.83, Form("#mu = (%f #pm %f) GeV/#it{c}^{2}", gauss_exp->GetParameter(1), gauss_exp->GetParError(1)));
   t.DrawLatex(0.59, 0.78, Form("#sigma = (%f #pm %f) GeV/#it{c}^{2}", gauss_exp->GetParameter(2), gauss_exp->GetParError(2)));
   t.DrawLatex(0.59, 0.73, Form("S = %f #pm %f", S, sS));
   t.DrawLatex(0.59, 0.68, Form("B = %f #pm %f", B, sB));
   t.DrawLatex(0.59, 0.63, Form("#frac{S}{S+B} =  %f #pm %f", purezza, s_purezza));             //per scrivere il prodotto vettoriale:   #times 10^{-2}
   //t.DrawLatex(0.2, 0.35, Form("#frac{S}{#sqrt{S+B}} = %f #pm %f", significatività, s_significatività));    posso anche non mettere la significatività
    
 
 


    plotTesi.cd();
    CMassV0plot->Write(); 
    MassV0_plot->Write();
    CMassV0plot->Print("plot.pdf");  //mi salva il plot in pdf e aggiorna lo stesso pdf se faccio rigirare il codice
    plotTesi.Close();
 

/* questo codice mi permette di cancellare un istogramma (o un altro tipo di file) da un root file
 TFile *plotTesi( TFile::Open("plotTesi.root", "UPDATE") ); 
    plotTesi->cd();
    gDirectory->Delete("MassV0new;1");
    plotTesi->Close();
 */

  
    }



    
/* per fare gli assi anche a destra e in alto. In realtà, c'è un modo più semplice (gStyle->SetPadTickX(1))
float exp_xmin= MassV0_exp->GetXaxis()->GetBinLowEdge(1);
float exp_xmax= MassV0_exp->GetXaxis()->GetBinUpEdge(nbin);
float exp_ymin= 0;
float exp_ymax= 6.3;

      TGaxis *exp_Yaxis = new TGaxis(exp_xmax,exp_ymin,exp_xmax,exp_ymax,exp_ymin,exp_ymax,610,"+L");  
    exp_Yaxis->SetLabelSize(0.035);
    exp_Yaxis->SetLabelFont(40);
    exp_Yaxis->Draw();

    TGaxis *exp_Xaxis = new TGaxis(exp_xmin,exp_ymax,exp_xmax,exp_ymax,exp_xmin,exp_xmax,615,"-L");  
    exp_Xaxis->SetLabelSize(0.035);
    exp_Xaxis->SetLabelFont(40);
    exp_Xaxis->Draw();
 */
