int check_hists()
{
  gStyle->SetOptStat(0);
  TFile* jf = TFile::Open("output_jet70_1_21.root","READ");
  TFile* hf = TFile::Open("/sphenix/user/hanpuj/CaloDataAna24_ppg09/offline/analysis_unfold_test/output_sim/output_r04_Jet70GeV_0_20.root");

  TH2D* jh = (TH2D*)jf->Get("h_respmatrix_all_nosmear");
  TH2D* hh = (TH2D*)hf->Get("h_respmatrix_all");

  cout << jh->GetEntries() << endl << hh->GetEntries() << endl;
  cout << jh->Integral(1,10,1,13) << endl << hh->Integral(1,10,1,13) << endl;
  TH2D* rh = (TH2D*)jh->Clone("rathist");
  rh->GetZaxis()->SetTitle("Ratio");
  rh->Divide(jh,hh);

  TCanvas* c = new TCanvas("","",1000,1000);
  c->SetRightMargin(0.2);
  rh->GetYaxis()->SetRangeUser(7,75);

  rh->Draw("COLZ");

  c->SaveAs("respmatrix_all.pdf");

  jf->Close();
  hf->Close();
  
  return 0;

}
