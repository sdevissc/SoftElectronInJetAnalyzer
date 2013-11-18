#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1F.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <./EventElecComm.C>
#include <TPaveStats.h>
using namespace std;


class HistoMaker
{
  public:
	HistoMaker(const TString &,const TString &);
	~HistoMaker();
	void SetStyle(TH1F *,int, int, const TString&,const TString&, float, float);
        void SetStyle(TProfile *,int, int, const TString&,const TString&, float, float);
	void Combine();

	TH1F *Pt_gen_all;
	TH1F *Pt_gen_MatchedWithASeed;
	TH1F *Pt_gen_MatchedWithAGedGsfElec;
        TH1F *Eta_gen_all;
	TH1F *Eta_gen_MatchedWithASeed;
	TH1F *Eta_gen_MatchedWithAGedGsfElec;

	TH1F *EfficiencyVsPt_Seed;
        TH1F *EfficiencyVsEta_Seed;
        TH1F *EfficiencyVsPt_GedGsfElec;
        TH1F *EfficiencyVsEta_GedGsfElec;

	TH1F* ForIntegral;
	TGraph *eEffvspiEff;
	float vecXmva[40],vecYmva[40];

};

HistoMaker::HistoMaker(const TString& inputname,const TString & tag)
{
	 TH1::SetDefaultSumw2(kTRUE);
	cout<<"define the binning"<<endl;
	float ptx[17]= {0.0 , 2.0 , 3.0 , 4.0 , 6.0 , 8.0 , 10.0 , 13.0 , 16.0  , 20.0 , 24.0 , 30.0 , 40.0 , 55.0 , 70.0 , 90.0 , 120.0};
	float etax[14]={-2.5 , -2.1  , -1.7 , -1.3 , -0.9 , -0.5 , -0.1 , 0.1 , 0.5 , 0.9 , 1.3 , 1.7 , 2.1 , 2.5 };
	 cout<<"define pt and eta basic histos"<<endl;
	ForIntegral=new TH1F(inputname+tag+"forIntegral","",40,-100000,100000);	

	int nbinpt=16,nbineta=13;
        Pt_gen_all=new TH1F(inputname+tag+"Pt_gen_all","",nbinpt,ptx);
	cout<<"ok"<<endl;
        Eta_gen_all=new TH1F(inputname+tag+"Eta_gen_all","",nbineta,etax);
	
	Pt_gen_MatchedWithASeed=new TH1F(inputname+tag+"Pt_gen_MatchedWithSeed","",nbinpt,ptx);
        Eta_gen_MatchedWithASeed=new TH1F(inputname+tag+"Eta_gen_MatchedWithSeed","",nbineta,etax);
	
	Pt_gen_MatchedWithAGedGsfElec=new TH1F(inputname+tag+"Pt_gen_MatchedWithGedGsfElec","",nbinpt,ptx);
        Eta_gen_MatchedWithAGedGsfElec=new TH1F(inputname+tag+"Eta_gen_MatchedWithGedGsfElec","",nbineta,etax);
	
	 cout<<"define the efficiency histos"<<endl;
	EfficiencyVsPt_Seed=new TH1F(inputname+tag+"effptseed","",nbinpt,ptx);
        EfficiencyVsEta_Seed=new TH1F(inputname+tag+"effetaseed","",nbineta,etax);
        EfficiencyVsPt_GedGsfElec=new TH1F(inputname+tag+"effptged","",nbinpt,ptx);
        EfficiencyVsEta_GedGsfElec=new TH1F(inputname+tag+"effgedeta","",nbineta,etax);


	SetStyle(Pt_gen_all,1,22,"p_{T} of the generated electron","",0.001,1.2);
	SetStyle(Eta_gen_all,1,22,"eta_{T} of the generated electron","",0.001,1.2);	

        SetStyle(Pt_gen_MatchedWithASeed,1,22,"p_{T} of the generated electron","",0.001,1.2);
        SetStyle(Pt_gen_MatchedWithAGedGsfElec,2,22,"p_{T} of the generated electron","",0.001,1.2);

        SetStyle(Eta_gen_MatchedWithASeed,1,22,"#eta of the generated electron","",0.001,1.2);
        SetStyle(Eta_gen_MatchedWithAGedGsfElec,2,22,"#eta of the generated electron","",0.001,1.2);

	SetStyle(EfficiencyVsPt_Seed,1,22,"p_{T} of the generated electron","",0.001,1.2);
        SetStyle(EfficiencyVsPt_GedGsfElec,2,22,"p_{T} of the generated electron","",0.001,1.2);
	
	SetStyle(EfficiencyVsEta_Seed,1,22,"#eta of the generated electron","",0.001,1.2);
	SetStyle(EfficiencyVsEta_GedGsfElec,2,22,"#eta of the generated electron","",0.001,1.2);
	
	

}


void HistoMaker::SetStyle(TH1F *histo,int color, int style, const TString& xtitle,const TString& ytitle, float min, float max){
	histo->SetLineColor(color);
        histo->SetMarkerColor(color);
        histo->SetMarkerStyle(style);
        histo->GetXaxis()->SetTitle(xtitle);
	histo->GetYaxis()->SetTitle(ytitle);
		histo->SetMinimum(min);
      	  	histo->SetMaximum(max);
}

void HistoMaker::SetStyle(TProfile *histo,int color, int style, const TString& xtitle,const TString& ytitle, float min, float max){
        histo->SetLineColor(color);
        histo->SetMarkerColor(color);
        histo->SetMarkerStyle(style);
        histo->GetXaxis()->SetTitle(xtitle);
        histo->GetYaxis()->SetTitle(ytitle);
        if(min!=-777.0 && max!=-777.0 ){
                histo->SetMinimum(min);
                histo->SetMaximum(max);
        }
}

HistoMaker::~HistoMaker()                 // destructor, just an example
{
}





void HistoMaker::Combine(){
	TH1::SetDefaultSumw2(kTRUE);
	EfficiencyVsPt_Seed->Add(Pt_gen_MatchedWithASeed);
        EfficiencyVsPt_Seed->Divide(Pt_gen_all);
        EfficiencyVsEta_Seed->Add(Eta_gen_MatchedWithASeed);
        EfficiencyVsEta_Seed->Divide(Eta_gen_all);

	EfficiencyVsPt_GedGsfElec->Add(Pt_gen_MatchedWithAGedGsfElec);
        EfficiencyVsPt_GedGsfElec->Divide(Pt_gen_all);
        EfficiencyVsEta_GedGsfElec->Add(Eta_gen_MatchedWithAGedGsfElec);
        EfficiencyVsEta_GedGsfElec->Divide(Eta_gen_all);
	
	cout<<EfficiencyVsPt_Seed->Integral()<<" "<<EfficiencyVsEta_Seed->Integral()<<" "<<EfficiencyVsPt_GedGsfElec->Integral()<<" "<<EfficiencyVsEta_GedGsfElec->Integral()<<endl;
}


class SLPlotter                   // begin declaration of the class
{
  public:                    // begin public section
    	SLPlotter(const TString &,const TString &);     // constructor
    	~SLPlotter();                  // destructor
	void initialize();
	void setTDRStyle();
	EventElecComm *Evt;
	int nentries;
	TFile *file;
        TTree *tp,*te;
	void getPlot(int,TString&); 
	void GetEfficiencies(int,TString &);
//	void GetEscOverPGen(TString &);
	void FillHistos();
	void TheFill(HistoMaker *,EventElecComm *);
//	void MakeEffvsEffPlot(TString &);
	int nmvasteps;
//	TH1F* GetEffPvsPt();
//	TH1F* GetEffPvsEta();

	TString DefElecFromB;
	TString DefElecFromV;
	TString DefPiKaInJet;
	TString LimitBarrelEndcap;
	TString ProcessFile;
        TString ProcName;
	TString TheTag;
	HistoMaker *electron;
	HistoMaker *pion;
	HistoMaker *kaon;
	HistoMaker *electron_lowPU;
        HistoMaker *electron_midPU;
        HistoMaker *electron_highPU;
        HistoMaker *pion_lowPU;
        HistoMaker *pion_midPU;
        HistoMaker *pion_highPU;
	HistoMaker *kaon_lowPU;
        HistoMaker *kaon_midPU;
        HistoMaker *kaon_highPU;	
	

		
	TString ProcessName;
	
}
;

SLPlotter::SLPlotter(const TString &input,const TString &tag)
{
	ProcessFile = input;
	TheTag=tag;
	electron=new HistoMaker("elec",TheTag);
	pion=new HistoMaker("pion",TheTag);
	kaon=new HistoMaker("kaon",TheTag);
	electron_lowPU=new HistoMaker("elec_lpu",TheTag);
	electron_midPU=new HistoMaker("elec_mpu",TheTag);
	electron_highPU=new HistoMaker("elec_hpu",TheTag);
        pion_lowPU=new HistoMaker("pion_lpu",TheTag);
	pion_midPU=new HistoMaker("pion_mpu",TheTag);
	pion_highPU=new HistoMaker("pion_hpu",TheTag);
	kaon_lowPU=new HistoMaker("kaon_lpu",TheTag);
        kaon_midPU=new HistoMaker("kaon_mpu",TheTag);
        kaon_highPU=new HistoMaker("kaon_hpu",TheTag);

	
}

void SLPlotter::initialize(){
        ifstream Processes;
        Processes.open(ProcessFile);
        string processline;
        Processes >>   ProcName;
        Evt = new EventElecComm(ProcName);
        nentries =/* 1500000;*/Evt->fChain->GetEntriesFast();
}


SLPlotter::~SLPlotter()                 // destructor, just an example
{
}

void SLPlotter::getPlot(int log,TString& process){
	ProcessName = process;
//	GetEscOverPGen(ProcessName);
	FillHistos();
	GetEfficiencies(log,ProcessName);	
}


void SLPlotter::GetEfficiencies(int log,TString &process){
	TH1::SetDefaultSumw2(kTRUE);
	ProcessName = process;
	setTDRStyle();
        TH1::SetDefaultSumw2(kTRUE);
        TCanvas *can = new TCanvas("h", "bla",600,600);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.05);
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);
	TLegend *legend=new TLegend(0.65,0.8,0.95,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
	

//-----------------------------------

	gPad->SetGridx();
        gPad->SetGridy();
	if(log==1)gPad->SetLogy(1);
	//Plotting Pt plot


	// GSFTrack_all vs GSFTracl_ecalSeeded
	legend->AddEntry(electron->EfficiencyVsPt_Seed,"Seedk","p");
	legend->AddEntry(electron->EfficiencyVsPt_GedGsfElec,"GED Gsf elec","p");	
        legend->AddEntry(electron_lowPU->EfficiencyVsPt_Seed,"Seedk lowPU","p");
        legend->AddEntry(electron_lowPU->EfficiencyVsPt_GedGsfElec,"GED Gsf elec low PU","p");
        legend->AddEntry(electron_highPU->EfficiencyVsPt_Seed,"Seedk high PU","p");
        legend->AddEntry(electron_highPU->EfficiencyVsPt_GedGsfElec,"GED Gsf elec high PU","p");

	electron->EfficiencyVsPt_Seed->SetMaximum(1.2);
	electron->EfficiencyVsPt_Seed->Draw("pe1");
        electron->EfficiencyVsPt_GedGsfElec->Draw("pe1same");

	electron_lowPU->EfficiencyVsPt_Seed->SetMarkerStyle(23);
	electron_lowPU->EfficiencyVsPt_GedGsfElec->SetMarkerStyle(23);	
	electron_highPU->EfficiencyVsPt_Seed->SetMarkerStyle(24);
        electron_highPU->EfficiencyVsPt_GedGsfElec->SetMarkerStyle(24);

	electron_lowPU->EfficiencyVsPt_Seed->Draw("pe1same");
        electron_lowPU->EfficiencyVsPt_GedGsfElec->Draw("pe1same");
	electron_highPU->EfficiencyVsPt_Seed->Draw("pe1same");
        electron_highPU->EfficiencyVsPt_GedGsfElec->Draw("pe1same");	
	legend->Draw();	

	can->Print(ProcessName+"_Electron_RecoEff_vs_Pt.pdf");	

	electron->EfficiencyVsEta_Seed->SetMaximum(1.2);
	electron->EfficiencyVsEta_Seed->SetMinimum(0.0001);	
	electron->EfficiencyVsEta_Seed->Draw("pe1");
        electron->EfficiencyVsEta_GedGsfElec->Draw("pe1same");

	electron_lowPU->EfficiencyVsEta_Seed->SetMarkerStyle(23);
        electron_lowPU->EfficiencyVsEta_GedGsfElec->SetMarkerStyle(23);
        electron_highPU->EfficiencyVsEta_Seed->SetMarkerStyle(23);
        electron_highPU->EfficiencyVsEta_GedGsfElec->SetMarkerStyle(23);

        electron_lowPU->EfficiencyVsEta_Seed->Draw("pe1same");
        electron_lowPU->EfficiencyVsEta_GedGsfElec->Draw("pe1same");
        electron_highPU->EfficiencyVsEta_Seed->Draw("pe1same");
        electron_highPU->EfficiencyVsEta_GedGsfElec->Draw("pe1same");

	legend->Draw();
	
        can->Print(ProcessName+"_Electron_RecoEff_vs_Eta.pdf");
	//-------------------------------------------------------------------

	gPad->SetLogy(1);	
	pion->EfficiencyVsPt_Seed->SetMaximum(1.2);
        pion->EfficiencyVsPt_Seed->Draw("pe1");
        pion->EfficiencyVsPt_GedGsfElec->Draw("pe1same");
        legend->Draw();

        can->Print(ProcessName+"_Pion_RecoEff_vs_Pt.pdf");

        pion->EfficiencyVsEta_Seed->SetMaximum(1.2);
//        pion->EfficiencyVsEta_Seed->SetMinimum(0.0001);
        pion->EfficiencyVsEta_Seed->Draw("pe1");
        pion->EfficiencyVsEta_GedGsfElec->Draw("pe1same");	

	can->Print(ProcessName+"_Pion_RecoEff_vs_Eta.pdf");

	legend->Clear();

	delete legend;
	delete can;
}

void SLPlotter::TheFill(HistoMaker *obj,EventElecComm *Evt){
        int flagEta=fabs(Evt->etaGen)<1.4?0:1;

        obj->Pt_gen_all->Fill(Evt->ptGen);
        obj->Eta_gen_all->Fill(Evt->etaGen);
        if(Evt->isMatchedWithASeed==1){
	      	obj->Pt_gen_MatchedWithASeed->Fill(Evt->ptGen);
                obj->Eta_gen_MatchedWithASeed->Fill(Evt->etaGen);
        }
	if(Evt->isMatchedWithAGedGsfElec==1){
        	obj->Pt_gen_MatchedWithAGedGsfElec->Fill(Evt->ptGen);
                obj->Eta_gen_MatchedWithAGedGsfElec->Fill(Evt->etaGen);
	}	
	//cout<<"histo integral "<<obj->Pt_gen_MatchedWithASeed->Integral()<<" "<<obj->Eta_gen_MatchedWithASeed->Integral()<<endl;	
}


void SLPlotter::FillHistos(){
	for(int ientry=0; ientry<nentries; ientry++) {
		if(ientry % 20000==0)cout<<"Event number "<<ientry<<endl;
                Evt->GetEntry(ientry);
		bool basicselection= Evt->ptGen>2 && fabs(Evt->etaGen)<2.4 ;
        	bool isfromB=(Evt->origin==4 || Evt->origin==6);
        	bool isfromD=(Evt->origin==2);
        	bool isfromX=(Evt->origin==0);
	//	if(Evt->ptGen>2 && fabs(Evt->etaGen)<2.4)cout<<Evt->ptGen<<" "<<Evt->etaGen<<" "<<isfromB<<" "<<abs(Evt->pdgId)<<endl;
		if(basicselection && isfromB && abs(Evt->pdgId)==11){
	//		cout<<"an electron spotted"<<endl;
			TheFill(electron,Evt);
			if(Evt->nPV<10){
				TheFill(electron_lowPU,Evt);
			}
			if(Evt->nPV>=10 && Evt->nPV<25){
				TheFill(electron_midPU,Evt);
			}
			if(Evt->nPV>=25){
				TheFill(electron_highPU,Evt);
			}
                }
		if(basicselection && abs(Evt->pdgId)==211){
			TheFill(pion,Evt);
			if(Evt->nPV<10){
				TheFill(pion_lowPU,Evt);
			}
                        if(Evt->nPV>=10 && Evt->nPV<25){
				TheFill(pion_midPU,Evt);
			}
                        if(Evt->nPV>=25){
				TheFill(pion_highPU,Evt);
			}
                }
		if(basicselection && abs(Evt->pdgId)==311){
                        TheFill(pion,Evt);
                        if(Evt->nPV<10){
                                TheFill(kaon_lowPU,Evt);
                        }
                        if(Evt->nPV>=10 && Evt->nPV<25){
                                TheFill(kaon_midPU,Evt);
                        }
                        if(Evt->nPV>=25){
                                TheFill(kaon_highPU,Evt);
                        }
                }

        }
	electron->Combine();
	electron_lowPU->Combine();
	electron_highPU->Combine();
	pion->Combine();
	kaon->Combine();


}


void SLPlotter::setTDRStyle() {
        gStyle->SetCanvasBorderMode(0);
        gStyle->SetCanvasColor(kWhite);
        gStyle->SetCanvasDefH(1500); //Height of canvas
        gStyle->SetCanvasDefW(1500); //Width of canvas
        gStyle->SetCanvasDefX(0);   //POsition on screen
        gStyle->SetCanvasDefY(0);

        gStyle->SetPadBorderMode(0);
        gStyle->SetPadColor(kWhite);
        gStyle->SetPadGridX(false);
        gStyle->SetPadGridY(false);
        gStyle->SetGridColor(0);
        gStyle->SetGridStyle(3);
        gStyle->SetGridWidth(1);

        gStyle->SetFrameBorderMode(0);
        gStyle->SetFrameBorderSize(0.1);
        gStyle->SetFrameFillColor(0);
        gStyle->SetFrameFillStyle(0);
        gStyle->SetFrameLineColor(1);
        gStyle->SetFrameLineStyle(1);
        gStyle->SetFrameLineWidth(0.1);

        gStyle->SetHistLineColor(1);
        gStyle->SetHistLineStyle(0);
        gStyle->SetHistLineWidth(0.1);

        gStyle->SetEndErrorSize(2);
//        gStyle->SetErrorX(0.);
        gStyle->SetMarkerStyle(20);

        gStyle->SetOptFit(1);
        gStyle->SetFitFormat("5.4g");
        gStyle->SetFuncColor(2);
        gStyle->SetFuncStyle(1);
        gStyle->SetFuncWidth(1);

        gStyle->SetOptDate(0);
        gStyle->SetOptStat(0);

        // Margins:
        gStyle->SetPadTopMargin(0.05);
        gStyle->SetPadBottomMargin(0.13);
        gStyle->SetPadLeftMargin(0.16);
        gStyle->SetPadRightMargin(0.02);

        // For the Global title:

        gStyle->SetOptTitle(0);
        gStyle->SetTitleFont(42);
        gStyle->SetTitleColor(1);
        gStyle->SetTitleTextColor(1);
        gStyle->SetTitleFillColor(10);
        gStyle->SetTitleFontSize(0.05);

        // For the axis titles:

        gStyle->SetTitleColor(1, "XYZ");
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetTitleSize(0.06, "XYZ");
        // gStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
        // gStyle->SetTitleYSize(Float_t size = 0.02);
        gStyle->SetTitleXOffset(0.9);
        gStyle->SetTitleYOffset(1.25);
        // gStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

        // For the axis labels:

        gStyle->SetLabelColor(1, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetLabelOffset(0.007, "XYZ");
        gStyle->SetLabelSize(0.05, "XYZ");

        // For the axis:

        gStyle->SetAxisColor(1, "XYZ");
        gStyle->SetStripDecimals(kTRUE);
        gStyle->SetTickLength(0.03, "XYZ");
        gStyle->SetNdivisions(510, "XYZ");
        gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
        gStyle->SetPadTickY(1);

        // Change for log plots:
        gStyle->SetOptLogx(0);
        gStyle->SetOptLogy(0);
        gStyle->SetOptLogz(0);

        gROOT->ForceStyle();

}
/*
TH1F* SLPlotter::GetEffPvsPt(){
        return electron->EfficiencyVsPt_PF_all;
}
TH1F* SLPlotter::GetEffPvsEta(){
        return electron->EfficiencyVsEta_PF_all;
}

class PlotTogether
{
        public:
                PlotTogether();
                ~PlotTogether();
		void EffPlotCombination();
};

PlotTogether::PlotTogether()
{

}

PlotTogether::~PlotTogether()
{
}

void PlotTogether::EffPlotCombination(){
	SLPlotter a("","");
	a.setTDRStyle();
        TCanvas *can = new TCanvas("h", "bla",600,600);
        can->SetLeftMargin(0.17);
        can->SetTopMargin(0.05);
        can->SetRightMargin(0.05);
        can->SetBottomMargin(0.15);

        SLPlotter *TT=new SLPlotter((TString)"input_tt.txt",(TString)"tt");
        TT->initialize();
        TT->FillHistos();
	cout<<"getting the histos"<<endl;
	TH1F *EffPt_tt=TT->GetEffPvsPt();
        TH1F *EffEta_tt=TT->GetEffPvsEta();
	EffPt_tt->SetMarkerColor(kBlue);
        EffEta_tt->SetMarkerColor(kBlue);
	cout<<"getting the plots"<<endl;
	TLegend *legend=new TLegend(0.65,0.8,0.95,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
        legend->AddEntry(EffPt_qcd,"QCD","p");
        legend->AddEntry(EffPt_tt,"t#bar{t}","p");
        legend->AddEntry(EffPt_wbb,"Wbb","p");

	EffPt_tt->SetMaximum(1.2);
	EffPt_tt->SetMinimum(0.0);
	EffPt_tt->Draw("p");
	legend->Draw();
	can->Print("UltimateTest_pt.pdf");
	
	EffEta_tt->SetMaximum(1.2);
        EffEta_tt->SetMinimum(0.0);

	EffEta_tt->Draw("p");
	legend->Draw();
        can->Print("UltimateTest_eta.pdf");

}
*/
