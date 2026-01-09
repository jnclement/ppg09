#include <dlUtility.h>

int plot_raw(int zsel = 0, int cutsel = 0, int radind = 4, int it = 2)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  const int ncut = 3;
  string cuttype[ncut] = {"both","dijet","frac"};
  const int nzt = 3;
  string ztype[nzt] = {"all","zvertex30","zvertex60"};
  
  TFile* rawdatf = TFile::Open(("input/output_data_r0"+to_string(radind)+".root").c_str(),"READ");
  TFile* epcdatf = TFile::Open(("output/output_data_"+ztype[zsel]+"_"+cuttype[cutsel]+"_purityefficiency_r0"+to_string(radind)+".root").c_str(),"READ");
  TFile* unfdatf = TFile::Open(("output/output_unfolded_"+ztype[zsel]+"_"+cuttype[cutsel]+"_r0"+to_string(radind)+".root").c_str(), "READ");
  TFile* hanpuep = TFile::Open(("/sphenix/user/hanpuj/CaloDataAna24_ppg09/offline/analysis_unfold/output_data_puritycorr_r0"+to_string(radind)+".root").c_str(),"READ");
  TFile* hanpunf = TFile::Open(("/sphenix/user/hanpuj/CaloDataAna24_ppg09/offline/analysis_unfold/output_unfolded_r0"+to_string(radind)+".root").c_str(), "READ");

  TH1D* h_meas = (TH1D*)rawdatf->Get(("h_calibjet_pt_eff_"+cuttype[cutsel]+"_"+ztype[zsel]).c_str());
  TH1D* h_corr[2];
  TH1D* h_ufld[2];
  
  h_corr[0] = (TH1D*)epcdatf->Get(("h_calibjet_pt_puritycorr_"+cuttype[cutsel]+"_"+ztype[zsel]).c_str());
  h_ufld[0] = (TH1D*)unfdatf->Get(("h_unfold_"+cuttype[cutsel]+"_"+ztype[zsel]+"_"+to_string(it)).c_str());

  h_corr[1] = (TH1D*)hanpuep->Get(("h_calibjet_pt_puritycorr_"+cuttype[cutsel]+"_"+ztype[zsel]).c_str());
  h_ufld[1] = (TH1D*)hanpunf->Get(("h_unfold_"+cuttype[cutsel]+"_"+ztype[zsel]).c_str());

  TCanvas* c = new TCanvas("","",1000,1000);
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  h_meas->SetMarkerColor(kBlack);
  h_meas->SetMarkerSize(1.5);
  h_meas->SetLineColor(kBlack);
  h_meas->SetMarkerStyle(20);
  h_meas->GetYaxis()->SetTitle("Raw Counts");
  h_meas->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  for(int i=0; i<2; ++i)
    {
      h_corr[i]->SetMarkerColor((i==0?kRed:kCyan)+1);
      h_corr[i]->SetMarkerSize(1.5);
      h_corr[i]->SetLineColor((i==0?kRed:kCyan)+1);
      h_corr[i]->SetMarkerStyle(i==0?21:20);
      h_ufld[i]->SetMarkerColor((i==0?kRed:kCyan)+1);
      h_ufld[i]->SetMarkerSize(1.5);
      h_ufld[i]->SetLineColor((i==0?kRed:kCyan)+1);
      h_ufld[i]->SetMarkerStyle(i==0?21:20);
      h_ufld[i]->GetYaxis()->SetTitle("Unfolded Counts");
      h_ufld[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
      h_corr[i]->GetYaxis()->SetTitle("Purity Corrected Counts");
      h_corr[i]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
    }

  int hibinufld = h_ufld[0]->GetNbinsX()+1;
  int hibincorr = h_corr[0]->GetNbinsX()+1;

  cout << "pT Bin Hanpu Cts   Joey Cts" << endl;
  for(int i=1; i<hibinufld; ++i)
    {
      cout << h_ufld[1]->GetBinLowEdge(i) << " " <<h_ufld[1]->GetBinLowEdge(i+1) << " " << h_ufld[1]->GetBinContent(i) << " " << h_ufld[0]->GetBinContent(i) << endl;
    }

  cout << endl << endl << endl << "pT BLE Hanpu Cts   Joey Cts" << endl;
  
 for(int i=1; i<hibincorr; ++i)
    {
      cout << h_corr[1]->GetBinLowEdge(i) << " " << h_corr[1]->GetBinLowEdge(i+1) << " " << h_corr[1]->GetBinContent(i) << " " << h_corr[0]->GetBinContent(i) << endl;
    }

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);

  TLegend* leg = new TLegend(0.5,0.7,0.8,0.8);
  leg->AddEntry(h_corr[0],"Joey","p");
  leg->AddEntry(h_corr[1],"Hanpu","p");

  h_corr[0]->Draw("PE");
  h_corr[1]->Draw("SAME PE");
  leg->Draw();
  sphenixtext(0.65,0.97);
  sqrt_s_text(0.65,0.925);
  c->SaveAs("output/drawpecorr.png");

  h_ufld[0]->Draw("PE");
  h_ufld[1]->Draw("SAME PE");
  leg->Draw();
  sphenixtext(0.65,0.97);
  sqrt_s_text(0.65,0.925);
  c->SaveAs("output/drawunf.png");

  h_meas->Draw("PE");
  sphenixtext(0.65,0.97);
  sqrt_s_text(0.65,0.925);
  c->SaveAs("output/drawraw.png");
  

  return 0;
}
