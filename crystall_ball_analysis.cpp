
//per il momento non sono riuscita a farci niente perchè il mio root non mi prende la funzione crystal ball

using namespace std;

void crystall_ball_analysis(){

    TFile fmass("AnalysisResults.root");   
    TH1F *MassV0=(TH1F*)fmass.Get("efficiency-q-a/massV0");   
    

    TCanvas *CMassV0 = new TCanvas("CMassV0 ","MassV0",200,50,600,400); 
     gStyle->SetOptStat(0); 
     

    //ora devo fittare con una crystal ball più un esponenziale decrescente
    TF1 *crystal_exp = new TF1("crystal_exp","[p0]*ROOT::Math::crystalball_function(x,[p1],[p2],[p3],[p4])+expo(5)",7,0.,1.);   //la crystal ball ha 5 parametri

    //sono più importanti i SetParLimits dei SetParameters
    crystal_exp->SetParLimits(0,0,1000);   //ho una decina di bin (forse un po' di meno) con 100 conteggi-> il parametro [0] è il parametro di normalizzazione, è dell'ordine del numero di conteggi. È un valore positivo
    crystal_exp->SetParLimits(1,0.47,0.52);   //è la media quindi so che deve essere vicino a 490 MeV/c2 ma il grafico è in GeV
    crystal_exp->SetParLimits(2,0,0.1);   //è la dev std. Per stimare la dev std devo prendere la metà larghezza a metà altezza. È dell'ordine di 0.01. È un valore positivo
    crystal_exp->SetParLimits(3,0,10);   
    crystal_exp->SetParLimits(4,-10,0); 
    crystal_exp->SetParLimits(5,0,10);    //l'esponensiale è definito come exp^(a+bx)  a =[0,10] circa
    crystal_exp->SetParLimits(6,-10,0);    // b deve essere negativo -> [-10,0]

    crystal_exp->SetLineColor(2); 


    auto fitResult=MassV0->Fit(crystal_exp,"RMLS");   //L: likelihood  -> usa il principio di massima verosomiglianza   S: get fit result pointer teoricamente
    CMassV0->cd(0);

    MassV0->DrawCopy();

    cout<<"Probability: "<<crystal_exp->GetProb()<<endl;


    auto covMatrix = fitResult->GetCovarianceMatrix();
    cout<<"\n\nMatrice di covarianza \n";
    covMatrix.Print();

     
    //mu= 4.94632e-001 +- 1.50858e-004
    //sigma= 6.59342e-003  +- 1.42971e-004

    float m= 0.494632-0.00659346*3;
    float M= 0.494632+0.00659346*3;
    int nbin= MassV0->GetNbinsX();     //funzionalità dell'istogramma quindi devo mettere nome dell'istogramma
    const int n=400;                   //valore messo alla luce di quanto mi viene restituito dalla riga predecente
    cout<<endl;
    cout<<"numero di bin su x: "<<nbin<<endl;
    cout<<endl;


    cout<<"min: "<<m<<" GeV/c2"<<endl;
    cout<<"max: "<<M<<" GeV/c2"<<endl;
    cout<<endl;
    float limsx[n];
    for(int i=0;i<nbin;i++){
        limsx[i]=MassV0->GetXaxis()->GetBinLowEdge(i);        //nome del TH1F->GetXaxis()->GetBinLowEdge(i) mi crea un array con tutti i limiti inferiori dei bin
        if(limsx[i]>m){
            double min=limsx[i];
            cout<<"Questo è il limite inferiore dell'intervallo di integrazione [GeV/c2]: "<<min;   
            break;
        }
    }

    cout<<endl;
   cout<<endl;
    float limdx[n];   //l'argomento dell'array deve essere una static const int ma poi nel for la static const int non funziona. Devo mettere l'altra cosa (solo int)
    for(int i=0;i<nbin;i++){
        limdx[i]=MassV0->GetXaxis()->GetBinUpEdge(i);    
        if(limdx[i]>M){
            double max=limdx[i-1];
            cout<<"Questo è il limite superiore dell'intervallo di integrazione [GeV/c2]: "<<max;   
            break;
        }
    }


    cout<<"\n"<<endl;
    float integrale=crystal_exp->Integral(min, max);   //metto gli estremi di integrazione
    float err_integrale=crystal_exp->IntegralError(min, max); 
    cout<<"Il valore dell'integrale è: "<<integrale<<" +- "<<err_integrale<<endl;
    cout<<"\n"<<endl;




    float integrale_exp, err_integrale_exp;
    //la funzionalità Integral() e IntegralError() valgono solo per funzioni con cui è stato effettuato un fit -> non posso usarlo per expo soltanto
    //tuttavia posso calcolarmi l'integrale a mano  -> integrale di e^(a+bx)= (1/b)*e^(a+bx) da calcolare tra gli estremi dati
    //per gli errori, si calcolano le derivate parziali dell'integrale rispetto a tutte le quantità con incertezza, ossia a e b e come incertezze si usano
    //i valori riportati nella matrice di covarianza del fit precedente
    //i valori della matrice di covarianza che mi interessano sono quelli di parametro 3 con parametro 3, parametro 3 con parametro 4 e parametro 4 con parametro 4


    float a=crystal_exp->GetParameter(3);
    float b=crystal_exp->GetParameter(4);

    integrale_exp= (1/b)*(exp(a+b*max) - exp(a+b*min));
    float der_int_a=(1/b)*(exp(a+b*max) - exp(a+b*min));  //derivata dell'integrale rispetto ad a
    float der_int_b=(-1/pow(b,2))*(exp(a+b*max) - exp(a+b*min))+(1/b)*(max*exp(a+b*max) - min*exp(a+b*min));  //derivata dell'integrale rispetto ad b
    
    err_integrale_exp=sqrt(pow(der_int_a,2)*covMatrix(3,3) + pow(der_int_b,2)*covMatrix(4,4) + 2*der_int_a*der_int_b*covMatrix(3,4));  //regola di propagazione dell'errore   covMatrix non sono alla seconda perchè nella formula ho la DEV STD alla seconda, che sarebbe la VARIANZA, ossia quello che ho, già fatto, nella matrice di covarianza

    cout<<"\n Il valore dell'integrale del solo esponenziale decrescente è: "<<integrale_exp<<" +- "<<err_integrale_exp<<endl;





    //ora definisco le quantità che devo calcolare
    float S,B, purezza, significatività;
    S= integrale-integrale_exp;   //segnale
    B= integrale_exp;             //background(fondo)
    purezza= S/(S+B);
    significatività= S/sqrt(S+B);

    float sS,sB, s_purezza, s_significatività;
    sS= sqrt(pow(err_integrale,2) + pow(err_integrale_exp,2));
    sB= err_integrale_exp;
    s_purezza= sqrt( pow(B*sS/pow(S+B,2),2) + pow(S*sB/pow(S+B,2),2) );
    s_significatività= sqrt( pow((S/2+B)*sS/pow(S+B,3/2),2) + pow((S/2)*sB/pow(S+B,3/2),2) );



 cout<<"\n\nS= "<<S<<" +- "<<sS<<endl;
 cout<<"B= "<<B<<" +- "<<sB<<endl;
cout<<"purezza= "<<purezza<<" +- "<<s_purezza<<endl;
cout<<"significatività= "<<significatività<<" +- "<<s_significatività<<endl;




   TFile fmass("AnalysisResults.root"); 
    TH1F *MassV0new=(TH1F*)fmass.Get("efficiency-q-a/massV0");  
    

    TCanvas *CMassV0new = new TCanvas("CMassV0new ","MassV0new",200,50,600,400); 
     gStyle->SetOptStat(0); 
 
    MassV0new->Fit(crystal_exp,"RMLS"); 

    CMassV0new->cd(0);

    MassV0new->DrawCopy();

      TLegend *leg1 = new TLegend(0.57,0.88,0.89,0.45);    //posizione del lato verticale sx, pos. del lato orizzontale in alto, pos. lato verticale dx, pos. lato orizzontale in basso
  TLegendEntry *frase1 = leg1->AddEntry(crystal_exp,"S = (1.15 #pm 0.06) entry","");
  frase1->SetTextSize(0.03);
  TLegendEntry *frase2 = leg1->AddEntry(crystal_exp,"B =  (0.05 #pm 0.04) entry","");
  frase2->SetTextSize(0.03);
  TLegendEntry *frase3 = leg1->AddEntry(crystal_exp,"#frac{S}{S+B} =  0.96 #pm 0.03","");
  frase3->SetTextSize(0.03);
  TLegendEntry *frase4 = leg1->AddEntry(crystal_exp,"#frac{S}{#sqrt{S+B}} = 1.05 #pm 0.04","");
  frase4->SetTextSize(0.03);
  leg1->SetFillColor(0);
  leg1->SetTextSize(0.038);
  leg1->Draw();





}