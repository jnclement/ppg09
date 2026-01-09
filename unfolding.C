#include "../analysis_unfold/unfold_Def.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"


int toys(RooUnfoldBayes& unfold, TH1D* h_unfold, int ntoy)
{
  vector<TVectorD> val;
  vector<TVectorD> err;
  vector<double> chi2;

  unfold.SetNToys(ntoy);
  unfold.RunToys(ntoy, val, err, chi2);

  int nbin = val[0].GetNoElements();

  TVectorD avg(nbin);
  TVectorD var(nbin);
  
  for(int i=0; i<nbin; ++i)
    {
      double sum = 0;
      for(int j=0; j<ntoy; ++j)
	{
	  sum += val[j][i];
	}
      avg[i] = sum/ntoy;
    }
  
  for(int i=0; i<nbin; ++i)
    {
      double sdif = 0;
      for(int j=0; j<ntoy; ++j)
	{
	  sdif += pow(val[j][i] - avg[i],2);
	}
      var[i] = sdif/(ntoy-1);
    }
  int nx = h_unfold->GetNbinsX();

  for(int i=1; i<nx+1; ++i)
    {
      int vindex = (i-1);
      int binerr2 = pow(h_unfold->GetBinError(i),2);
      h_unfold->SetBinError(i,sqrt(binerr2+var[vindex]));
    }

  return 0;
}

int get_unfolded_spectrum(TH2D* h_resp, TH1D* h_spec, TH1D*& h_unfolded, int nit, string histname, TH1D* h_eff, int ntoy = 1000)
{

  TH1D* h_meas = (TH1D*)h_resp->ProjectionX("h_meas");
  TH1D* h_truth = (TH1D*)h_resp->ProjectionY("h_truth");

  RooUnfoldResponse* resp = new RooUnfoldResponse(h_meas, h_truth, h_resp, "response", "");
  RooUnfoldBayes unfold(resp, h_spec, nit);
  cout << "test0" <<h_spec << " " << resp << endl;
  h_unfolded = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
    cout << "test1.5" << endl;
  h_unfolded->SetName(histname.c_str());

  cout << "test2" << endl;
  toys(unfold,h_unfolded,ntoy);
  h_unfolded->Divide(h_eff);

  return 0;
}

int unfolding(int cutsel = 0, int zsel = 0, int radind = 4)
{

  const int ncut = 3;
  string cuttype[ncut] = {"both","dijet","frac"};
  const int nz = 3;
  string ztype[nz] = {"all","zvertex30","zvertex60"};

  TFile *f_resp = new TFile(("input/output_sim_r0"+to_string(radind)+".root").c_str());
  TFile *f_sim = new TFile(("output/output_sim_"+ztype[zsel]+"_"+cuttype[cutsel]+"_purityefficiency_r0"+to_string(radind)+".root").c_str(), "READ");
  TFile *f_data = new TFile(("output/output_data_"+ztype[zsel]+"_"+cuttype[cutsel]+"_purityefficiency_r0"+to_string(radind)+".root").c_str(), "READ");
  TFile *f_out = new TFile(("output/output_unfolded_"+ztype[zsel]+"_"+cuttype[cutsel]+"_r0"+to_string(radind)+".root").c_str(), "RECREATE");
  
  
  const int nit = 20;
  const int ntoy = 1000;
  const int nsys = 9;

  string systs[nsys] = {"","_jesup","_jesdown","_jerup","_jerdown","_jetup","_jetdown","_mbdup","_mbddown"};
  
  TH1D* h_unfold[nsys][nit];
  TH2D* h_resp[nsys];
  TH1D* h_spec[nsys];
  TH1D* h_eff[nsys];

  TFile* f_truth = TFile::Open("input/output_sim_r04.root","READ");
  TH1D* h_truth = (TH1D*)f_truth->Get(("h_truth_"+cuttype[cutsel]+"_"+ztype[zsel]).c_str());


  f_out->cd();
  
  for(int i=0; i<nsys; ++i)
    {
      cout << "test" << endl;
      if(i<5)
	{
	  h_resp[i] = (TH2D*)f_resp->Get(("h_respmatrix_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]).c_str());
	  h_spec[i] = (TH1D*)f_data->Get(("h_calibjet_pt_puritycorr_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]).c_str());
	  h_eff[i] = (TH1D*)f_sim->Get(("h_efficiency_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]).c_str());
	}
      else
	{
	  h_resp[i] = (TH2D*)f_resp->Get(("h_respmatrix_"+cuttype[cutsel]+"_"+ztype[zsel]).c_str());
	  h_resp[i]->SetName(("h_respmatrix_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]).c_str());
	  h_spec[i] = (TH1D*)f_data->Get(("h_calibjet_pt_puritycorr_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]).c_str());
	  cout << h_spec[i] << endl;
	  h_eff[i] = (TH1D*)f_sim->Get(("h_efficiency_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]).c_str());
	}
      cout << 2 << endl;
      for(int j=0; j<nit; ++j)
	{
	  cout << 3 << endl;
	  string thename = "h_unfold_"+cuttype[cutsel]+"_"+ztype[zsel]+systs[i]+"_"+to_string(j);
	  get_unfolded_spectrum(h_resp[i], h_spec[i], h_unfold[i][j], j, thename, h_eff[i], ntoy);
	  cout << h_unfold[i][j] << endl;
	  h_unfold[i][j]->Write();
	}
    }

  h_truth->Write();

  f_out->Close();
  
  return 0;
}
