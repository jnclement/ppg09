#include "/sphenix/user/jocl/projects/chi2checker/src/dlUtility.h"

int draw_respmatrix(string type = "all", int radius_index = 4)
{

  double calibptbins[] = {7, 11, 15, 19, 24, 29, 35, 41, 48, 56, 65, 75};
  double truthptbins[] = {7, 11, 15, 19, 24, 29, 35, 41, 48, 56, 65, 75, 86};

  double calibptbins_big[] = {7, 11, 15, 19, 24, 29, 35, 41, 48, 56, 65, 75, 86};
  double truthptbins_big[] = {7, 11, 15, 19, 24, 29, 35, 41, 48, 56, 65, 75, 86, 95};
  
  int calibnpt = sizeof(calibptbins) / sizeof(calibptbins[0]) - 1;
  int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;
  
  TFile* file = TFile::Open(("../output_comb_r0"+to_string(radius_index)+".root").c_str(),"READ");//("/sphenix/user/hanpuj/CaloDataAna24_ppg09/offline/analysis_unfold_test/output_sim_r0"+to_string(radius_index)+".root").c_str());

  TH2D* respmatrix = (TH2D*)file->Get(("h_respmatrix_"+type).c_str());
  TH1D* hmiss = (TH1D*)file->Get(("h_miss_"+type).c_str());
  TH1D* hfake = (TH1D*)file->Get(("h_fake_"+type).c_str());

  TH2D* h_normcal_withmiss = new TH2D("hncwm",";p_{T}^{calib} [GeV];p_{T}^{truth} [GeV];Row Normalized Counts",calibnpt+1,calibptbins_big,truthnpt,truthptbins);

  TH2D* h_normtruth_withfake = new TH2D("hntwf",";p_{T}^{calib} [GeV];p_{T}^{truth} [GeV];Column Normalized Counts",calibnpt,calibptbins,truthnpt+1,truthptbins_big);

  TPaveText* miss = new TPaveText(0.7,0.065,0.85,0.099,"NDC");
  miss->SetFillColor(kWhite);
  miss->SetBorderSize(0);
  miss->AddText("Miss");

  TPaveText* fake = new TPaveText(0.05,0.75,0.099,0.87,"NDC");
  fake->SetFillColor(kWhite);
  fake->SetBorderSize(0);
  fake->AddText("Fake");
  
  
  float rowsums[12] = {0};
  float colsums[11] = {0};

  for(int i=1; i<calibnpt+1; ++i)
    {
      colsums[i-1] += abs(respmatrix->Integral(i,i,1,13));
      cout << "before: " << colsums[i-1] << endl;
      colsums[i-1] += abs(hfake->GetBinContent(i));
      cout << "after:" << colsums[i-1] << endl;
    }

  cout << endl << endl;
  for(int i=1; i<truthnpt+1; ++i)
    {
      rowsums[i-1] += abs(respmatrix->Integral(1,13,i,i));
      cout << "before: " << rowsums[i-1] << endl;
      rowsums[i-1] += abs(hmiss->GetBinContent(i));
      cout << "after:" << rowsums[i-1] << endl;
    }


  for(int i=1; i<truthnpt+1; ++i)
    {
      for(int j=0; j<calibnpt+2; ++j)
	{
	  if(j<calibnpt+1) h_normcal_withmiss->SetBinContent(j,i,abs(respmatrix->GetBinContent(j,i))/rowsums[i-1]);
	  else h_normcal_withmiss->SetBinContent(j,i,abs(hmiss->GetBinContent(i))/rowsums[i-1]);
	}
    }

  for(int i=1; i<calibnpt+1; ++i)
    {
      for(int j=0; j<truthnpt+2; ++j)
	{
	  if(j<truthnpt+1) h_normtruth_withfake->SetBinContent(i,j,abs(respmatrix->GetBinContent(i,j))/colsums[i-1]);
	  else h_normtruth_withfake->SetBinContent(i,j,abs(hfake->GetBinContent(i))/colsums[i-1]);
	}
    }
  

  


  respmatrix->GetYaxis()->SetRangeUser(7,86);
  respmatrix->GetXaxis()->SetRangeUser(7,86);
  respmatrix->GetZaxis()->SetTitle("Reweighted Counts");
  respmatrix->GetZaxis()->SetTitleOffset(1.57);
  gStyle->SetOptStat(0);

  gStyle->SetPalette(kRainbow);
  
  TCanvas* c = new TCanvas("","",1000,1000);
  c->SetTopMargin(0.15);
  c->SetRightMargin(0.175);
  //c->SetLogz();

  respmatrix->Draw("COLZ");
  maintexts(0.96,0.7,0,0.03,0,0,radius_index);
  
  c->SaveAs("respmatrix_nosys.pdf");

  cout << h_normcal_withmiss->Integral(1,13,2,2)<< " " << h_normtruth_withfake->Integral(2,2,1,13) << endl;

  h_normcal_withmiss->GetZaxis()->SetTitleOffset(1.57);
  h_normcal_withmiss->GetZaxis()->SetRangeUser(h_normcal_withmiss->GetMinimum(),h_normcal_withmiss->GetMaximum());
  h_normcal_withmiss->Draw("COLZ");
  maintexts(0.96,0.7,0,0.03,0,0,radius_index);
  if(type=="all")drawText("No z_{vtx} required",0.7,0.86,0,kBlack,0.03);
  else drawText("|z_{vtx}|<60 cm required",0.7,0.86,0,kBlack,0.03);
  miss->Draw();
  c->SaveAs(("h_respmatrix_rownormed_"+type+"_withmiss_r0"+to_string(radius_index)+".pdf").c_str());

  h_normtruth_withfake->GetZaxis()->SetRangeUser(h_normtruth_withfake->GetMinimum(),h_normtruth_withfake->GetMaximum());
  h_normtruth_withfake->GetZaxis()->SetTitleOffset(1.57);
  h_normtruth_withfake->Draw("COLZ");
  maintexts(0.96,0.7,0,0.03,0,0,radius_index);
  if(type=="all") drawText("No z_{vtx} required",0.7,0.86,0,kBlack,0.03);
  else drawText("|z_{vtx}|<60 cm required",0.7,0.86,0,kBlack,0.03);
  fake->Draw();
  TText* tt = fake->GetLineWith("Fake");
  tt->SetTextAngle(90);
  gPad->Modified();
  c->SaveAs(("h_respmatrix_colnormed_"+type+"_withfake_r0"+to_string(radius_index)+".pdf").c_str());
  

  file->Close();

  return 0;
}
