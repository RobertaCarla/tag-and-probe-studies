using namespace std;


void efficiency(){
    TFile fAnalysis("AnalysisResults-mc.root");
    TFile efficiency("efficiency.root", "RECREATE");

    TH2F *piRec= (TH2F*)fAnalysis.Get("efficiency-q-a/piRec");

    TH1F *Decays= (TH1F*)piRec->ProjectionY("Decays",1,1);
    TH1F *ITS_TPC= (TH1F*)piRec->ProjectionY("ITS_TPC",6,6);
    TH1F *ITSwithTPC= (TH1F*)piRec->ProjectionY("ITS with TPC",8,8);


    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);



//--------------------------------------------------------- decays ---------------------------------------------------------//
    int nbin= Decays->GetNbinsX();
    float width= Decays->GetBinWidth(1);

    float bin_ledge[nbin];
    float min;
    int n;     //n= numero di bin che avrà ciascuno dei due istogrammi
    for(int i=1; i<=nbin; i++){
        bin_ledge[i]= Decays->GetBinLowEdge(i);
        if( abs(bin_ledge[i]-0.)<width/2){
            min=i;
            cout<<"il bin minimo è: "<<min<<endl;
            n= min-1;
            break;
        }                               
    }

    cout<<"n: "<<n<<endl;



//------------------- decays pi+ -----------------//

    float massimox= Decays->GetXaxis()->GetBinUpEdge(nbin);
    TH1F *pi_pos_d= new TH1F("pi_pos_d", "pi+ decays;#it{p_{T}} (GeV/#it{c});Entries",n,0.,massimox);
    pi_pos_d->SetName("decays_pi_pos");
    //auto *axisX1 = new TGaxis(0,0.,1000,0.,"",510,"");
    //axisX1->SetTitle("p_{T} GeV");
    //axisX1->Draw();
   
    for(int i=1; i<=n; i++){
       pi_pos_d->SetBinContent(i, Decays->GetBinContent(min-1+i));
    }


    TCanvas *Cpi_pos_d= new TCanvas("Cpi_pos_d", "pi+ decays", 200,50,600,400);
    Cpi_pos_d->cd(0);
    pi_pos_d->Draw();

    efficiency.cd();   
    Cpi_pos_d->Write();
    Cpi_pos_d->Print("decay pi+.pdf");



    
//------------- decays pi- -----------------------//
    TH1F *pi_neg_d= new TH1F("pi_neg_d", "pi- decays;#it{p_{T}} (GeV/#it{c});Entries",n,0.,massimox);
    pi_neg_d->SetName("decays_pi_neg");


    for(int i=1; i<=n; i++){
       pi_neg_d->SetBinContent(i, Decays->GetBinContent(min-i));
    }


    TCanvas *Cpi_neg_d= new TCanvas("Cpi_neg_d", "pi- decays", 200,50,600,400);
    Cpi_neg_d->cd(0);
    pi_neg_d->Draw();

    efficiency.cd();   
    Cpi_neg_d->Write();
    Cpi_neg_d->Print("decay pi-.pdf");













//--------------------------------------------------------- ITS+TPC ---------------------------------------------------------//


//------------------- ITS+TPC pi+ -----------------//

    TH1F *pi_pos_IT= new TH1F("pi_pos_ITS_TPC", "pi+ ITS+TPC;#it{p_{T}} (GeV/#it{c});Entries",n,0.,massimox);
    pi_pos_IT->SetName("ITS+TPC pi_pos");
    //pi_pos_IT->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    //pi_pos_IT->GetYaxis()->SetTitle("Entries");

   
    for(int i=1; i<=n; i++){
       pi_pos_IT->SetBinContent(i, ITS_TPC->GetBinContent(min-1+i));
    }


    TCanvas *Cpi_pos_IT= new TCanvas("Cpi_pos_IT", "pi+ ITS+TPC", 200,50,600,400);
    Cpi_pos_IT->cd(0);
    pi_pos_IT->Draw();

    efficiency.cd();   
    Cpi_pos_IT->Write();
    Cpi_pos_IT->Print("ITS+TPC pi+.pdf");



    
//------------- ITS+TPC pi- -----------------------//
    TH1F *pi_neg_IT= new TH1F("pi_neg_ITS_TPC", "pi- ITS+TPC;#it{p_{T}} (GeV/#it{c});Entries",n,0.,massimox);
    pi_neg_IT->SetName("ITS+TPC pi_neg");


    for(int i=1; i<=n; i++){
       pi_neg_IT->SetBinContent(i, ITS_TPC->GetBinContent(min-i));
    }


    TCanvas *Cpi_neg_IT= new TCanvas("Cpi_neg_IT", "pi- ITS+TPC", 200,50,600,400);
    Cpi_neg_IT->cd(0);
    pi_neg_IT->Draw();

    efficiency.cd();   
    Cpi_neg_IT->Write();
    Cpi_neg_IT->Print("ITS+TPC pi-.pdf");  
    












//--------------------------------------------------------- ITS with TPC ---------------------------------------------------------//


//------------------- ITS with TPC pi+ -----------------//

    TH1F *pi_pos_W= new TH1F("pi_pos_ITSwithTPC", "pi+ ITS with TPC;#it{p_{T}} (GeV/#it{c});Entries",n,0.,massimox);
    pi_pos_W->SetName("ITS with TPC pi_pos");

   
    for(int i=1; i<=n; i++){
       pi_pos_W->SetBinContent(i, ITSwithTPC->GetBinContent(min-1+i));
    }


    TCanvas *Cpi_pos_W= new TCanvas("Cpi_pos_W", "pi+ ITS with TPC", 200,50,600,400);
    Cpi_pos_W->cd(0);
    pi_pos_W->Draw();

    efficiency.cd();   
    Cpi_pos_W->Write();
    Cpi_pos_W->Print("ITS with TPC pi+.pdf");



    
//------------- ITS with TPC pi- -----------------------//
    TH1F *pi_neg_W= new TH1F("pi_neg_ITSwithTPC", "pi- ITS with TPC;#it{p_{T}} (GeV/#it{c});Entries",n,0.,massimox);
    pi_neg_W->SetName("ITS with TPC pi_neg");


    for(int i=1; i<=n; i++){
       pi_neg_W->SetBinContent(i, ITSwithTPC->GetBinContent(min-i));
    }


    TCanvas *Cpi_neg_W= new TCanvas("Cpi_neg_W", "pi- ITS with TPC", 200,50,600,400);
    Cpi_neg_W->cd(0);
    pi_neg_W->Draw();

    efficiency.cd();   
    Cpi_neg_W->Write();
    Cpi_neg_W->Print("ITS with TPC pi-.pdf");  










    //------------------------- calcolo efficienze -------------------------//
    cout<<"\n\n-------------------------------------------------- CALCOLO DELLE EFFICIENZE --------------------------------------------------"<<endl;
    cout<<"\n\n------------------------- Calcoli per pi+ -------------------------"<<endl;



    //------------------- efficienza globale --------------------//
    TCanvas *CE_glob_pos = new TCanvas("CE_glob_pos","Efficienza globale pi+",200,10,700,500);


    TH1F *gE_glob_pos = new TH1F(*pi_pos_IT);
    gE_glob_pos->SetName("efficiency_pi_pos_glob");
    gE_glob_pos->Divide(pi_pos_IT, pi_pos_d, 1., 1., "B");
    gE_glob_pos->SetTitle("Efficienza globale #pi^{+}");
    gE_glob_pos->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    gE_glob_pos->GetYaxis()->SetTitle("#epsilon #times A");
    
    gE_glob_pos->GetYaxis()->SetRangeUser(0., 1.1);
    gE_glob_pos->GetXaxis()->SetRangeUser(0., 3.);

    CE_glob_pos->cd();
    gE_glob_pos->Draw("PE");

    efficiency.cd();
    CE_glob_pos->Write();
    CE_glob_pos->Print("efficienza globale pi+.pdf");




    //------------------- efficienza tpc --------------------//

        TCanvas *CE_tpc_pos = new TCanvas("CE_tpc_pos","Efficienza TPC pi+",200,10,700,500);

   TH1F *Copy_pi_pos_IT= new TH1F(*pi_pos_IT);
   Copy_pi_pos_IT->SetName("copia_pi_pos_IT");
   Copy_pi_pos_IT->Add(pi_pos_W);

   TH1F *gE_tpc_pos= new TH1F(*Copy_pi_pos_IT);
   gE_tpc_pos->SetName("tpc_efficiency_pi_pos");
   gE_tpc_pos->Divide(Copy_pi_pos_IT, pi_pos_d, 1., 1., "B");   //B mi calcolal'errore binomiale
    gE_tpc_pos->SetTitle("Efficienza tpc #pi^{+}");
    gE_tpc_pos->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    gE_tpc_pos->GetYaxis()->SetTitle("#epsilon #times A");
    
    gE_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
    gE_tpc_pos->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così
 

    CE_tpc_pos->cd();
    gE_tpc_pos->Draw("PE");    //P: mi disegna i punti    E: mi disegna gli errori

    efficiency.cd();
    CE_tpc_pos->Write();
    CE_tpc_pos->Print("efficienza TPC pi+.pdf");




    

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
    gE_its_tpc_pos->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    gE_its_tpc_pos->GetYaxis()->SetTitle("#epsilon #times A");
    
    gE_its_tpc_pos->GetYaxis()->SetRangeUser(0., 1.1);
    gE_its_tpc_pos->GetXaxis()->SetRangeUser(0., 3.);     //mi prende fino a 3 GeV ma è ok così
 
    CE_its_tpc_pos->cd();
    gE_its_tpc_pos->Draw("PE");

    efficiency.cd();
    CE_its_tpc_pos->Write();
    CE_its_tpc_pos->Print("efficienza matching ITS TPC pi+.pdf");







 

















    cout<<"\n\n------------------------- Calcoli per pi- -------------------------"<<endl;


    //------------------- efficienza gloabale --------------------//
    TCanvas *CE_glob_neg = new TCanvas("CE_glob_neg","Efficienza globale pi-",200,10,700,500);


   TH1F *gE_glob_neg = new TH1F(*pi_neg_IT);
   gE_glob_neg->SetName("glob_efficiency_pi_neg");
   gE_glob_neg->Divide(pi_neg_IT, pi_neg_d, 1., 1., "B");   //B mi calcolal'errore binomiale
   gE_glob_neg->SetName("Efficienza_globale_pi_neg");
    gE_glob_neg->SetTitle("Efficienza globale #pi^{-}");
    gE_glob_neg->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    gE_glob_neg->GetYaxis()->SetTitle("#epsilon #times A");
    
    gE_glob_neg->GetYaxis()->SetRangeUser(0., 1.1);
    gE_glob_neg->GetXaxis()->SetRangeUser(0., 3.); 

    CE_glob_neg->cd();
    gE_glob_neg->Draw("PE");

    efficiency.cd();
    CE_glob_neg->Write();
    CE_glob_neg->Print("efficienza globale pi-.pdf");




    //------------------- efficienza tpc --------------------//

        TCanvas *CE_tpc_neg = new TCanvas("CE_tpc_neg","Efficienza TPC pi-",200,10,700,500);

  TH1F *Copy_pi_neg_IT = new TH1F(*pi_neg_IT);
   Copy_pi_neg_IT ->SetName("copia_pi_neg_IT");
   Copy_pi_neg_IT ->Add(pi_neg_W);
    TH1F *gE_tpc_neg = new TH1F(*Copy_pi_neg_IT);
    gE_tpc_neg->SetName("Efficienza_tpc_pi_neg");
   gE_tpc_neg->Divide(Copy_pi_neg_IT, pi_neg_d, 1., 1., "B");   //B mi calcolal'errore binomiale
    gE_tpc_neg->SetTitle("Efficienza TPC #pi^{-}");
    gE_tpc_neg->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    gE_tpc_neg->GetYaxis()->SetTitle("#epsilon #times A");
    
    gE_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
    gE_tpc_neg->GetXaxis()->SetRangeUser(0., 3.); 

    CE_tpc_neg->cd();
    gE_tpc_neg->Draw("PE");

    efficiency.cd();
    CE_tpc_neg->Write();
    CE_tpc_neg->Print("efficienza TPC pi-.pdf");




    

    //------------------- efficienza matching its-tpc --------------------//

        TCanvas *CE_its_tpc_neg = new TCanvas("CE_its_tpc_neg","Efficienza matching ITS TPC pi-",200,10,700,500);

 
  TH1F *Copy_pi_neg_W = new TH1F(*pi_neg_W);
   Copy_pi_neg_W ->SetName("copia_pi_neg_W");
   Copy_pi_neg_W ->Add(pi_neg_IT);
    TH1F *gE_its_tpc_neg = new TH1F(*Copy_pi_neg_W);
    gE_its_tpc_neg->SetName("Efficienza_its_tpc_pi_neg");
   gE_its_tpc_neg->Divide(pi_neg_IT, Copy_pi_neg_W, 1., 1., "B");   //B mi calcolal'errore binomiale
    gE_its_tpc_neg->SetTitle("Efficienza matching ITS TPC #pi^{-}");
    gE_its_tpc_neg->GetXaxis()->SetTitle("#it{p_{T}} (GeV/#it{c})");
    gE_its_tpc_neg->GetYaxis()->SetTitle("#epsilon #times A");
    
    gE_its_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
    gE_its_tpc_neg->GetXaxis()->SetRangeUser(0., 3.); 

    CE_its_tpc_neg->cd();
    gE_its_tpc_neg->Draw("PE");

    efficiency.cd();
    CE_its_tpc_neg->Write();
    CE_its_tpc_neg->Print("efficienza matching ITS TPC pi-.pdf");














    //-------------------------------------------- efficienze per pi+ e pi- -------------------------------------//
    cout<<endl;
    cout<<endl;
    cout<<"--------------------------------------- efficienze per pi+ e pi- --------------------------------------"<<endl;

    //--------------------------------------------------- efficienza globale ------------------------------------------------------------//

    TCanvas *CE_glob = new TCanvas("CE_glob","Efficienza globale pi+ e pi-",200,10,700,500);

    gE_glob_neg->SetLineColor(2);

    gE_glob_pos->GetXaxis()->SetRangeUser(0., 3);
    gE_glob_neg->GetYaxis()->SetRangeUser(0., 1.1);
    gE_glob_neg->GetXaxis()->SetRangeUser(0., 3.); 

    gE_glob_neg->SetTitle("Efficienza globale #pi^{+} e #pi^{-}");

    CE_glob->cd();
    gE_glob_neg->Draw("PE");
    gE_glob_pos->Draw("PESAME");


    
   TLatex t1;
   t1.SetNDC();
   t1.SetTextSize(0.031);
   t1.DrawLatex(0.65, 0.83, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza globale #font[22]{#pi^{+}}");
   t1.DrawLatex(0.65, 0.78, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza globale #font[22]{#pi^{-}}");



    efficiency.cd();
    CE_glob->Write();
    CE_glob->Print("efficienza globale pi+ e pi-.pdf");





    
    //--------------------------------------------------- efficienza TPC ------------------------------------------------------------//

    TCanvas *CE_TPC = new TCanvas("CE_TPC","Efficienza TPC pi+ e pi-",200,10,700,500);

    gE_tpc_neg->SetLineColor(2);

    gE_tpc_pos->GetXaxis()->SetRangeUser(0., 3);
    gE_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
    gE_tpc_neg->GetXaxis()->SetRangeUser(0., 3.); 

    gE_tpc_neg->SetTitle("Efficienza TPC #pi^{+} e #pi^{-}");

    CE_TPC->cd();
    gE_tpc_neg->Draw("PE");
    gE_tpc_pos->Draw("PESAME");


    TLatex t2;
    t2.SetNDC();
    t2.SetTextSize(0.031);
    t2.DrawLatex(0.65, 0.83, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza TPC #font[22]{#pi^{+}}");
    t2.DrawLatex(0.65, 0.78, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza TPC #font[22]{#pi^{-}}");


    efficiency.cd();
    CE_TPC->Write();
    CE_TPC->Print("efficienza TPC pi+ e pi-.pdf");







    //--------------------------------------------------- efficienza ITS TPC ------------------------------------------------------------//

    TCanvas *CE_ITS_TPC = new TCanvas("CE_ITS_TPC","Efficienza ITS TPC pi+ e pi-",200,10,700,500);

    gE_its_tpc_neg->SetLineColor(2);

    gE_its_tpc_pos->GetXaxis()->SetRangeUser(0., 3);
    gE_its_tpc_neg->GetYaxis()->SetRangeUser(0., 1.1);
    gE_its_tpc_neg->GetXaxis()->SetRangeUser(0., 3.); 

    gE_its_tpc_neg->SetTitle("Efficienza ITS TPC #pi^{+} e #pi^{-}");

    CE_ITS_TPC->cd();
    gE_its_tpc_neg->Draw("PE");
    gE_its_tpc_pos->Draw("PESAME");


    TLatex t3;
    t3.SetNDC();
    t3.SetTextSize(0.031);
    t3.DrawLatex(0.65, 0.68, "#font[22]{#color[4]{#scale[2.2]{-}}} Efficienza ITS TPC #font[22]{#pi^{+}}");
    t3.DrawLatex(0.65, 0.63, "#font[22]{#color[2]{#scale[2.2]{-}}} Efficienza ITS TPC #font[22]{#pi^{-}}");

    efficiency.cd();
    CE_ITS_TPC->Write();
    CE_ITS_TPC->Print("efficienza matching ITS TPC pi+ e pi-.pdf");








 
    efficiency.Close(); 
}