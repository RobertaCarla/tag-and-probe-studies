using namespace std;


//ho dovuto scriverla fuori perchè C does not allow nested function definitions. Move these functions out of main 
    Double_t myfunction(Double_t *x, Double_t *par)
    {
        Float_t xx =x[0];
        Double_t f;
        if(xx<156 || xx>233)   //invece dei numeri metto par[2] e par[3]
            f = exp(par[0]+par[1]*xx);   //exp è il numero di nepero
        else
            f=0;
        return f;
    }





void cont_bin_analysis(){   //calcolo gli integrali non a partire da fit ma da conteggi di bin

    TFile fmass("AnalysisResults.root");    
    TH1F* MassV0=(TH1F*)fmass.Get("efficiency-q-a/massV0"); 
    TH1F *copyMassV0 = new TH1F(*MassV0);
    //TH1F* alias = MassV0;  //sta chiamando il puntatore vecchio in un nuovo modo, da non fare -->memory leak: un puntatore punta ancora ad un'informazione che è stata cancellata
    copyMassV0->SetName("V0_massa");   //il cambio nome è una cosa interna all'oggetto. Io da fuori continuo a chiamarlo "copyMassV0"
    const int nbin= 400;   //ottenuto da MassV0->GetNbinsX() ma, non so perché, se provo a lasciare int nbin=MassV0->GetNbinsX() non funziona per il ciclo for
    

    
    int bin_min= 156;
    int bin_max= 233;

    for(int i=bin_min; i<=bin_max; i++){
        copyMassV0->SetBinContent(i, 0);   //SetBinContent(nbin, content)
    }

 /*  serve per verificare che il ciclo for precedente faccia veramente quello che deve fare
     for(int i=1; i<=nbin; i++){
        cout<<i<<") "<<copyMassV0->GetBinContent(i)<<endl;  
    }
  */

    float min= 0.475148;
    float max= 0.514148;

 

  
   auto f1 = new TF1("myfunc",myfunction,0,1,2);   //myfunction è la funzione che ho definito
   //f1->SetParameters(1,2); se metto parameterS con la s, i valori tra parentesi sono le inizializzazioni di tutti i parametri: [0]=1, [1]=2
   //[2] e [3] SetParameter() e SetConstant(# parametro che devo settare costante) altrimenti li fitta

    f1->SetParLimits(0,0,10);   //l'esponensiale è definito come exp^(a+bx)  a =[0,10] circa
    f1->SetParLimits(1,-10,0);     // b deve essere negativo -> [-10,0]
   copyMassV0->Fit(f1, "RMLS"); 
   copyMassV0->Draw();

 cout<<"Probability: "<<f1->GetProb()<<endl;


   
    //unique_ptr<TFile> histo_exp( TFile::Open("histo_exp.root", "CREATE") );    //mi crea un file .root
    TFile *plotTesi( TFile::Open("plotTesi.root", "RECREATE") );     
    plotTesi->cd();
    plotTesi->WriteObject(copyMassV0, "copyMassV0");          //mi salva il mio oggetto all'interno del file .root   devo passare alla funzione plotTesi.root un puntatore quindi metto solo il nome del puntatore (se mettessi *copyMassV0 passerei l'oggetto puntato, non il puntatore)
    plotTesi->Close();
}