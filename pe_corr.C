int pe_corr(int zsel = 0, int cutsel = 0, int radind = 4)
{
  const int ncut = 3;
  string cuttype[ncut] = {"both","dijet","frac"};
  const int nzt = 3;
  string ztype[nzt] = {"all","zvertex30","zvertex60"};

  TFile* insimf = TFile::Open(("input/output_sim_r0"+to_string(radind)+".root").c_str(),"READ");
  TFile* indatf = TFile::Open(("input/output_data_r0"+to_string(radind)+".root").c_str(),"READ");

  TFile* outsimf = TFile::Open(("output/output_sim_"+ztype[zsel]+"_"+cuttype[cutsel]+"_purityefficiency_r0"+to_string(radind)+".root").c_str(),"RECREATE");
  TFile* outdatf = TFile::Open(("output/output_data_"+ztype[zsel]+"_"+cuttype[cutsel]+"_purityefficiency_r0"+to_string(radind)+".root").c_str(),"RECREATE");

  const int nsys = 9;
  string systs[nsys] = {"","_jesup","_jesdown","_jerup","_jerdown","_jetup","_jetdown","_mbdup","_mbddown"};
  
  TH1D* h_pur[nsys];
  TH1D* h_eff[nsys];
  TH1D* h_res[nsys];
  
  
  for(int i=0; i<nsys; ++i)
    {
      string dsuf = "_"+cuttype[cutsel]+"_"+ztype[zsel];
      string suff = "_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i];
      if(i<5)
	{
	  TH2D* h_resp = (TH2D*)insimf->Get(("h_respmatrix"+suff).c_str());
	  TH1D* h_truth = (TH1D*)insimf->Get(("h_truth"+suff).c_str());
	  TH1D* h_meas = (TH1D*)insimf->Get(("h_measure"+suff).c_str());
	  TH1D* h_fake = (TH1D*)insimf->Get(("h_fake"+suff).c_str());
	  TH1D* h_miss = (TH1D*)insimf->Get(("h_miss"+suff).c_str());
	  h_pur[i] = (TH1D*)h_resp->ProjectionX(("h_purity"+suff).c_str());
	  h_pur[i]->Divide(h_pur[i],h_meas,1,1,"B");
	  h_eff[i] = (TH1D*)h_resp->ProjectionY(("h_efficiency"+suff).c_str());
	  h_eff[i]->Divide(h_eff[i],h_truth,1,1,"B");
	}
      else
	{
	  h_pur[i] = (TH1D*)h_pur[0]->Clone(("h_purity"+suff).c_str());
	  h_eff[i] = (TH1D*)h_eff[0]->Clone(("h_efficiency"+suff).c_str());
	}
      outsimf->cd();
      h_eff[i]->Write();
      h_pur[i]->Write();

      h_res[i] = (TH1D*)((TH1D*)indatf->Get(("h_calibjet_pt_eff"+dsuf).c_str()))->Clone(("h_calibjet_pt_puritycorr"+suff).c_str());
      h_res[i]->Multiply(h_pur[i]);
      outdatf->cd();
      h_res[i]->Write();
    }
  return 0;
}
