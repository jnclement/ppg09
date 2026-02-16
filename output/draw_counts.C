#include "/sphenix/user/jocl/projects/chi2checker/src/dlUtility.h"

int draw_counts()
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile* file = TFile::Open("/sphenix/user/hanpuj/CaloDataAna24_ppg09/offline/analysis_unfold_test/output_data_r04.root","READ");

  TH1D* hist = (TH1D*) file->Get("h_recojet_pt_record_all");
  
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.5);
  hist->Rebin(10);
  hist->GetXaxis()->SetRangeUser(15,75);
  hist->GetXaxis()->SetTitle("p_{T}^{uncalib} [GeV]");
  hist->GetYaxis()->SetTitle("Counts");
  TCanvas* c = new TCanvas("","",1000,1000);

  c->SetLogy();
  c->SetTopMargin(0.15);
  hist->Draw("PE");
  maintexts();
  c->SaveAs("rawcounts.pdf");

  file->Close();

  return 0;

}
