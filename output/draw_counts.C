#include "/sphenix/user/jocl/projects/chi2checker/src/dlUtility.h"

int draw_counts(int radius_index = 4)
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile* file = TFile::Open("../ana_output/output_data_r04.root","READ");//("/sphenix/user/hanpuj/CaloDataAna24_ppg09/offline/analysis_unfold_test/output_data_r0"+to_string(radius_index)+".root").c_str(),"READ");

  TH1D* hist = (TH1D*) file->Get("h_calibjet_pt_record_all_noweight");
  
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1.5);
  hist->Rebin(10);
  hist->GetXaxis()->SetRangeUser(15,75);
  hist->GetXaxis()->SetTitle("p_{T}^{calib} [GeV]");
  hist->GetYaxis()->SetTitle("Counts");
  float max = hist->GetMaximum()*2;
  hist->GetYaxis()->SetRangeUser(0.5,max);
  TCanvas* c = new TCanvas("","",1000,1000);



  c->SetLogy();
  c->SetTopMargin(0.15);
  hist->Draw("PE");
  double calibptbins[8] = {19, 24, 29, 35, 41, 48, 56, 65};
  TLine* lines[8];
  for(int i=0; i<8; ++i)
    {
      lines[i] = new TLine(calibptbins[i],0.5,calibptbins[i],max);
      lines[i]->SetLineColor(kRed);
      lines[i]->SetLineWidth(2);
      lines[i]->SetLineStyle(0);
      lines[i]->Draw();
    }
  
  maintexts(0.96,0.7,0,0.03,1,0,radius_index);
  drawText("No z_{vtx} required",0.7,0.86,0,kBlack,0.03);
  c->SaveAs("rawcounts.pdf");

  file->Close();
  /*

  TH1D* hist2 = (TH1D*) file->Get("h_calibjet_pt_record_all");
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerSize(1.5);
  hist2->Rebin(10);
  hist2->GetXaxis()->SetRangeUser(15,75);
  hist2->GetXaxis()->SetTitle("p_{T}^{calib} [GeV]");
  hist2->GetYaxis()->SetTitle("Counts");
  */
  return 0;

}
