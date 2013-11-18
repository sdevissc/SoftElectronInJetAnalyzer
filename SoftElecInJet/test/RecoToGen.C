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
#include <THStack.h>
#include <./EventPurity.C>
#include <TPaveStats.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
using namespace std;


class HistoMaker
{
  public:
	HistoMaker(const TString &,const TString &);
	~HistoMaker();

        TH1F *var[29];
	TH1F* ForIntegral;
	int nbin[29];
	float xmin[29],xmax[29],ymin[29],ymax[29];
	TString VarName[29];
};

HistoMaker::HistoMaker(const TString& inputname,const TString & tag)
{

//	cout<<"define the binning"<<endl;
	xmin[0]=0;	xmax[0]=6;	ymin[0]=0;	ymax[0]=0.2;
        xmin[1]=0;	xmax[1]=4;	ymin[1]=0;	ymax[1]=0.2;
        xmin[2]=0;      xmax[2]=1;      ymin[2]=0;      ymax[2]=0.3;
        xmin[3]=-0.5;	xmax[3]=1.5;	ymin[3]=0;	ymax[3]=0.4;
        xmin[4]=-6;	xmax[4]=-2;	ymin[4]=0;	ymax[4]=0.4;
        xmin[5]=0;	xmax[5]=0.06;	ymin[5]=0;	ymax[5]=0.4;
        xmin[6]=-0.5;	xmax[6]=1.0;	ymin[6]=0;	ymax[6]=0.4;
        xmin[7]=0;	xmax[7]=5;	ymin[7]=0;	ymax[7]=0.4;
        xmin[8]=0;	xmax[8]=5;	ymin[8]=0;	ymax[8]=0.4;
        xmin[9]=0;	xmax[9]=20;	ymin[9]=0;	ymax[9]=0.4;
        xmin[10]=0;	xmax[10]=0.2;	ymin[10]=0;	ymax[10]=0.4;
        xmin[11]=-10;	xmax[11]=10;	ymin[11]=0;	ymax[11]=0.4;
        xmin[12]=-0.1;     xmax[12]=0.1;      ymin[12]=0;      ymax[12]=0.2;
        xmin[13]=0;     xmax[13]=0.3;      ymin[13]=0;      ymax[13]=0.4;
        xmin[14]=-0.1;     xmax[14]=0.1;      ymin[14]=0;      ymax[14]=0.4;
        xmin[15]=0;    xmax[15]=0.1;     ymin[15]=0;      ymax[15]=0.4;
        xmin[16]=0;     xmax[16]=0.3;    ymin[16]=0;      ymax[16]=0.4;
        xmin[17]=0;     xmax[17]=0.4;     ymin[17]=0;      ymax[17]=0.4;
        xmin[18]=0;     xmax[18]=1;     ymin[18]=0;      ymax[18]=0.4;
        xmin[19]=0;     xmax[19]=5;      ymin[19]=0;      ymax[19]=0.4;
        xmin[20]=0;     xmax[20]=4;     ymin[20]=0;      ymax[20]=0.4;
        xmin[21]=0;     xmax[21]=2;    ymin[21]=0;      ymax[21]=0.4;
        xmin[22]=0;     xmax[22]=0.06;      ymin[22]=0;      ymax[22]=0.2;
        xmin[23]=0;     xmax[23]=2;      ymin[23]=0;      ymax[23]=0.4;
        xmin[24]=0;     xmax[24]=1;      ymin[24]=0;      ymax[24]=0.4;
        xmin[25]=-2;  xmax[25]=2;     ymin[25]=0;      ymax[25]=0.4;
        xmin[26]=-100;     xmax[26]=100;    ymin[26]=0;      ymax[26]=0.4;
        xmin[27]=0; 	 xmax[27]=100;     ymin[27]=0;      ymax[27]=1;
        xmin[28]=-2.5;     xmax[28]=2.5;    ymin[28]=0;      ymax[28]=1;

//	 cout<<"define pt and eta basic histos"<<endl;
	ForIntegral=new TH1F(inputname+tag+"forIntegral","",40,-100000,100000);	
	VarName[0]="E_{tot}/P_{in}";
	VarName[1]="E_{e}/P_{out}";
	VarName[2]="f_{brem}";
	VarName[3]="E_{brem}/Delta P";
	VarName[4]="log#sigma(#eta#eta)";
	VarName[5]="#Delta#eta(T,seed)";
	VarName[6]="H/E";
	VarName[7]="#chi^{2}_{GSF}";
	VarName[8]="#chi^{2}_{KF}";
	VarName[9]="n Hits (KF)";
	VarName[10]="#sigma(P_{T})/P_{T}";
	VarName[11]= "ln(P_{T})";
        VarName[12]="deta";
        VarName[13]="dphi";
        VarName[14]="detacalo";
        VarName[15]="see";
        VarName[16]="etawidth";
        VarName[17]="phiwidth";
        VarName[18]="e1x5e5x5";
        VarName[19]="R9";
        VarName[20]="HoE";
        VarName[21]="EoP";
        VarName[22]="IoEmIoP";
        VarName[23]="eleEoPout";
        VarName[24]="d0";
        VarName[25]="ip3d";
        VarName[26]="ip3dSig";
	VarName[27]="gedgsfElecPt";
	VarName[28]="gedgsfElecEta";

	
	for(int u=0;u<29;u++){
		nbin[u]=90;	
		if(u==9)nbin[u]=20;	
		var[u]=new TH1F(Form(inputname+"var_%d",u),";"+VarName[u]+";Normalized scale",nbin[u],xmin[u],xmax[u]);;
	}
}


HistoMaker::~HistoMaker()                 // destructor, just an example
{
}


class SLPlotter                   // begin declaration of the class
{
  public:                    // begin public section
    	SLPlotter(const TString &,const TString &,int);     // constructor
    	~SLPlotter();                  // destructor
	void initialize();
	void setTDRStyle();
	EventPurity *Evt;
	int nentries;
        TTree *tp,*te;
	void SetStyle(TH1F *,int, int, float, float);
	void getDetailedPlot(int,TString&); 
        void getPurityPlots(int log,TString& process);
//        void getSummaryPlot(int,TString&);
//	void GetEscOverPGen(TString &);
	void FillHistos();
	void RecordHistos(HistoMaker *);
	void TheFill(HistoMaker *,EventPurity *);
	bool Condition(int ,EventPurity *);
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
	HistoMaker *Object[4][2];
	TString ProcessName;
	TMVA::Reader *reader;	
	int condInput;
}
;

SLPlotter::SLPlotter(const TString &input,const TString &tag, int thecond)
{
	ProcessFile = input;
	TheTag=tag;
	for(int u=0;u<4;u++)for(int v=0;v<2;v++){
		Object[u][v]=new HistoMaker(Form("obj_%d_%d",u,v),TheTag);
	}
	condInput=thecond;

	
}

void SLPlotter::SetStyle(TH1F *histo,int color, int style, float min, float max){
        histo->SetLineColor(color);
        histo->SetMarkerColor(color);
        if(style!=0){
		histo->SetFillStyle(style);
		histo->SetFillColor(color);
	}
	
	histo->GetYaxis()->SetTitleOffset(1.5);
        histo->SetMinimum(min);
        histo->SetMaximum(max);
}


void SLPlotter::initialize(){
        ifstream Processes;
        Processes.open(ProcessFile);
        string processline;
        Processes >>   ProcName;
        Evt = new EventPurity(ProcName);
        nentries =/* 1500000;*/Evt->fChain->GetEntriesFast();
	reader = new TMVA::Reader( "!Color:!Silent" );
}


SLPlotter::~SLPlotter()                 // destructor, just an example
{
}



void SLPlotter::getDetailedPlot(int log,TString& process){
	ProcessName = process;
//	GetEscOverPGen(ProcessName);
	FillHistos();
	setTDRStyle();	
	TLegend *legend=new TLegend(0.5,0.8,0.95,0.94);
        legend->SetBorderSize(1);
        legend->SetFillColor(0);
        legend->SetTextSize(0.025);
	for(int u=0;u<4;u++)for(int v=0;v<2;v++)for(int k=0;k<29;k++){
		if(Object[u][v]->var[k]->Integral()!=0)Object[u][v]->var[k]->Scale(1.0/Object[u][v]->var[k]->Integral());
		int style=0;	
		if(u==0 && v==0)style=3004;
		int color=3*u+v+1;
		SetStyle(Object[u][v]->var[k],color,style, 0,0.2);
	}
	legend->AddEntry(Object[0][0]->var[0],"gen elec from BC decay","f");
        legend->AddEntry(Object[1][0]->var[0],"gen pi/K","l");
        legend->AddEntry(Object[2][0]->var[0],"gen X","l");
        legend->AddEntry(Object[0][1]->var[0],"GEANT elec","l");
        legend->AddEntry(Object[1][1]->var[0],"GEANT pi/K","l");                           
        legend->AddEntry(Object[2][1]->var[0],"GEANT X","l");
	TLegend *legend2=new TLegend(0.5,0.8,0.95,0.94);
        legend2->SetBorderSize(1);
        legend2->SetFillColor(0);
        legend2->SetTextSize(0.025);	
	legend2->AddEntry(Object[0][0]->var[0],"gen elec from BC decay","f");
        legend2->AddEntry(Object[3][0]->var[0],"background","l");	
	for(int k=0;k<27;k++){
		TCanvas *can = new TCanvas("h", "bla",600,600);
        	can->SetLeftMargin(0.17);
        	can->SetTopMargin(0.05);
        	can->SetRightMargin(0.05);
        	can->SetBottomMargin(0.15);
	
//		cout<<"Checking the case of var "<<k<<endl;
		for(int u=0;u<3;u++)for(int v=0;v<2;v++){
//			cout<<"-->Drawing histo "<<u<<" "<<v<<endl;
        	        if(u==0 && v==0) {
//				cout<<"--> test1 "<<Object[u][v]->var[k]->Integral()<<endl;
				Object[u][v]->var[k]->Draw("");
			}
        	        else{
//				cout<<"--> test2 "<<Object[u][v]->var[k]->Integral()<<endl;
				Object[u][v]->var[k]->Draw("same");
			}	
//			cout<<"--> out"<<endl;
        	}
//		cout<<"out of double loop "<<k<<endl;
		legend->Draw();
//		cout<<"legend drawn"<<endl;
	
        	can->Print(Form(ProcessName+"_var_%d_cond_%d.pdf",k,condInput));	
//		cout<<"ok"<<endl;
		SetStyle(Object[0][0]->var[k],1,3004 ,  0,0.2);
		SetStyle(Object[3][0]->var[k],2,0, 0,0.2);
                Object[0][0]->var[k]->Draw("");
                Object[3][0]->var[k]->Draw("same");
                legend2->Draw();
                can->Print(Form("Summary_"+ProcessName+"_var_%d_cond_%d.pdf",k,condInput));
		delete can;
	}
	delete legend;
	delete legend2;
}

void SLPlotter::getPurityPlots(int log,TString& process){
//	cout<<"in Detailed"<<endl;
        ProcessName = process;
//      GetEscOverPGen(ProcessName);
        FillHistos();
        setTDRStyle();
	TH1F *h[3][2];
        for(int k=0;k<29;k++){
		TH1F *sum=(TH1F*)Object[0][0]->var[k]->Clone();
		sum->Reset();
//		cout<<"After reset sum for variable k="<<k<<" "<<sum->Integral()<<endl;
		TLegend *legend=new TLegend(0.7,0.8,0.95,0.94);
        	legend->SetBorderSize(1);
        	legend->SetFillColor(0);
        	legend->SetTextSize(0.025);
//		cout<<"----------> PLOT "<<k<<" <-----"<<endl;
		for(int u=0;u<3;u++)for(int v=0;v<2;v++){
                        h[u][v]=(TH1F*)Object[u][v]->var[k]->Clone();
			h[u][v]->SetMaximum(1.3);
			sum->Add(h[u][v]);
//			cout<<"bins after addition of u="<<u<<" and v="<<v<<" "<<h[u][v]->GetBinContent(1)<<" "<<sum->GetBinContent(1)<<endl;
			int style=0;
                        style=3004+2*u;
                        int color=1+3*v+u;
                        SetStyle(h[u][v],color,style, 0,1000);
		}
		legend->AddEntry(h[0][0],"gen elec from BC decay","f");
                legend->AddEntry(h[1][0],"gen pi/K","f");
                legend->AddEntry(h[2][0],"gen X","f");
                legend->AddEntry(h[0][1],"GEANT elec","f");
                legend->AddEntry(h[1][1],"GEANT pi/K","f");
                legend->AddEntry(h[2][1],"GEANT X","f");

//		cout<<"Integral of sum for variable k="<<k<<" "<<sum->Integral()<<endl;
		for(int u=0;u<3;u++)for(int v=0;v<2;v++){
//			cout<<"Integral of histo  u="<<u<<" and v="<<v<<" before division "<<h[u][v]->Integral()<<endl;
                        h[u][v]->Divide(sum);	
//			cout<<" after division, category u="<<u<<" and v="<<v<<" has an value of "<<h[u][v]->Integral()<<endl;
                }
//		cout<<"legend drawn"<<endl;
		THStack *hs = new THStack("hs","test stacked histograms");
                TCanvas *can = new TCanvas("h", "bla",600,600);
                can->SetLeftMargin(0.17);
                can->SetTopMargin(0.05);
                can->SetRightMargin(0.05);
                can->SetBottomMargin(0.15);

                for(int u=0;u<3;u++)for(int v=0;v<2;v++){
                        hs->Add(h[u][v]);
                }
                legend->Draw();
		hs->Draw("");
		legend->Draw();
		hs->GetXaxis()->SetTitle(Object[0][0]->VarName[k]);
		hs->GetYaxis()->SetTitle("Purity");
		hs->SetMaximum(1.2);
		hs->SetMinimum(0.0);
                can->Print(Form("Purity"+ProcessName+"_var_%d_cond_%d.pdf",k,condInput));
                delete can;
		delete hs;	
		delete legend;
		delete sum;
        }
	//for(int u=0;u<4;u++)
	//for(int v=0;v<2;v++){
       // 	delete h[u][v];
       // }
}



void SLPlotter::TheFill(HistoMaker *obj,EventPurity *Evt){
	obj->var[0]->Fill(Evt->EtotOvePin);
	obj->var[1]->Fill(Evt->EClusOverPout);
        obj->var[2]->Fill(Evt->fbrem);
        obj->var[3]->Fill(Evt->EBremOverDeltaP);
        obj->var[4]->Fill(Evt->logSigmaEtaEta);
        obj->var[5]->Fill(Evt->DeltaEtaTrackEcalSeed);
        obj->var[6]->Fill(Evt->HOverE);
        obj->var[7]->Fill(Evt->gsfchi2);
        obj->var[8]->Fill(Evt->kfchi2);
        obj->var[9]->Fill(Evt->kfhits);
        obj->var[10]->Fill(Evt->SigmaPtOverPt);
        obj->var[11]->Fill(Evt->lnPt);
        obj->var[12]->Fill(Evt->deta);
        obj->var[13]->Fill(Evt->dphi);
        obj->var[14]->Fill(Evt->detacalo);
        obj->var[15]->Fill(Evt->see);
        obj->var[16]->Fill(Evt->etawidth);
        obj->var[17]->Fill(Evt->phiwidth);
        obj->var[18]->Fill(Evt->e1x5e5x5);
        obj->var[19]->Fill(Evt->R9);
        obj->var[20]->Fill(Evt->HoE);
        obj->var[21]->Fill(Evt->EoP);
        obj->var[22]->Fill(Evt->IoEmIoP);
        obj->var[23]->Fill(Evt->eleEoPout);
        obj->var[24]->Fill(Evt->d0);
        obj->var[25]->Fill(Evt->ip3d);
        obj->var[26]->Fill(Evt->ip3dSig);
	obj->var[27]->Fill(Evt->pt);
	obj->var[28]->Fill(Evt->eta);	


/*

	std::cout<<"Filled histograms "<<Evt->EtotOvePin<<" "<<
					Evt->EClusOverPout<<" "<<
					Evt->fbrem<<" "<<
					Evt->EBremOverDeltaP<<" "<<
					Evt->logSigmaEtaEta<<" "<<
					Evt->DeltaEtaTrackEcalSeed<<" "<<
					Evt->HOverE<<" "<<
					Evt->Chi2GSF<<" "<<
					Evt->CHi2KF<<" "<<
					Evt->nHits<<" "<<
					Evt->SigmaPtOverPt<<" "<<
					Evt->lnPt<<std::endl;
*/		
}

bool  SLPlotter::Condition(int cond,EventPurity *Evt){
	bool ptcond[3],etacond[3];
	ptcond[0]=Evt->pt>2 && Evt->pt<5;	
	ptcond[1]= Evt->pt>5 && Evt->pt<10;
	ptcond[2]=Evt->pt>10;
	etacond[0]=fabs(Evt->eta)<0.8;
	etacond[1]=fabs(Evt->eta)>0.8 && fabs(Evt->eta)<1.4;
	etacond[2]=fabs(Evt->eta)>1.4;
	if(cond==0 && ptcond[0] && etacond[0])return true;
	if(cond==1 && ptcond[1] && etacond[0])return true;
	if(cond==2 && ptcond[2] && etacond[0])return true;
	if(cond==3 && ptcond[0] && etacond[1])return true;
	if(cond==4 && ptcond[1] && etacond[1])return true;
	if(cond==5 && ptcond[2] && etacond[1])return true;
        if(cond==6 && ptcond[0] && etacond[2])return true;
        if(cond==7 && ptcond[1] && etacond[2])return true;
        if(cond==8 && ptcond[2] && etacond[2])return true;
	return false;
}

void SLPlotter::FillHistos(){
	for(int ientry=0; ientry<nentries; ientry++) {
		if(ientry % 200==0)cout<<"Event number "<<ientry<<endl;
                Evt->GetEntry(ientry);
		if(Condition(condInput,Evt)){
	        	bool isfromB=(Evt->origin==4 || Evt->origin==6);
	       // 	bool isfromD=(Evt->origin==2);
	   //     	bool isfromX=(Evt->origin==0);	
			bool isAGenElec=abs(Evt->pdgId)==11;
			bool isAGenPionKaon=abs(Evt->pdgId)==211 || abs(Evt->pdgId)==321;	  
			bool isAnElec=abs(Evt->tppdgId)==11;
	                bool isAPionKaon=abs(Evt->tppdgId)==211 || abs(Evt->tppdgId)==321;
	
			if(isfromB && isAGenElec){
	//			cout<<"A gen elec"<<endl;
				TheFill(Object[0][0],Evt);
			}
			if(isAGenPionKaon){
	//			cout<<"A gen kapion"<<endl;
				TheFill(Object[1][0],Evt);
			}	
			if(Evt->Vtx==1 && !isAGenPionKaon && !isAGenElec){
	//			cout<<"A gen X"<<endl;
				TheFill(Object[2][0],Evt);
			}
	
			if(Evt->Vtx!=1 && isAnElec){
	//			cout<<"A not gen elec"<<endl;
				TheFill(Object[0][1],Evt);
			}
	                if(Evt->Vtx!=1 && isAPionKaon){
	//			cout<<"A not gen kapion"<<endl;
				TheFill(Object[1][1],Evt);
			}
	                if(Evt->Vtx!=1 && !isAPionKaon && !isAnElec){
	//			cout<<"A not gen X"<<endl;
				TheFill(Object[2][1],Evt);
			}
	//               if(!(isfromB && isAGenElec)){
			if(isAPionKaon){
				TheFill(Object[3][0],Evt);
			}	
		}
        }
	TFile *file = new TFile(Form("TH1FList_%d.root",condInput),"RECREATE");
	RecordHistos(Object[0][0]);
        RecordHistos(Object[1][0]);
        RecordHistos(Object[2][0]);
        RecordHistos(Object[0][1]);
        RecordHistos(Object[1][1]);
        RecordHistos(Object[2][1]);	
	file->Write();
        file->Close();
	
}

void SLPlotter::RecordHistos(HistoMaker *obj){
	for(int i=0;i<29;i++){
		obj->var[i]->Write();
	}
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
        gStyle->SetErrorX(0.);
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
        gStyle->SetTitleYOffset(1.5);
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
