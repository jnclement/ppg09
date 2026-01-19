#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include "unfold_Def.h"

void get_reweight_hist(TFile* f_out, TH1D*& h_reweight, std::string h_reweightname, TH1D* h_data, TH1D* h_sim) { 
  h_data->Rebin(10); h_data->Scale(1.0 / h_data->Integral());
  h_sim->Rebin(10); h_sim->Scale(1.0 / h_sim->Integral());
  h_reweight = (TH1D*)h_data->Clone(h_reweightname.c_str());
  h_reweight->Divide(h_sim);
  h_reweight->SetName(h_reweightname.c_str());

  std::string prefix_to_remove = "h_reweight_";
  std::string new_prefix = "reweightfunc_";
  std::string plot_name = new_prefix + h_reweightname.substr(prefix_to_remove.length());

  TF1* func_reweight = new TF1(plot_name.c_str(), "[0]*TMath::Exp(-[1]*x) + [2]", 0, 200);
  func_reweight->SetParameters(1., 0.1, 0.5);
  h_reweight->Fit(func_reweight, "", "", 15, 70);

  f_out->cd();
  h_reweight->Write();
  func_reweight->Write();
}

void get_reweighthist(int radius_index = 4) {

  SetAtlasStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // Read Files
  TFile* f_data = new TFile("output_data.root", "READ");
  TFile* f_sim = new TFile("output_sim.root", "READ");
  TFile* f_out = new TFile("output_reweightfunction.root", "RECREATE");

  // Read data histograms for reweighting
  TH1D* h_calibjet_pt_all = (TH1D*)f_data->Get("h_calibjet_pt_record_all");
  TH1D* h_calibjet_pt_all_jesup = (TH1D*)h_calibjet_pt_all->Clone("h_calibjet_pt_record_all_jesup");
  TH1D* h_calibjet_pt_all_jesdown = (TH1D*)h_calibjet_pt_all->Clone("h_calibjet_pt_record_all_jesdown");
  TH1D* h_calibjet_pt_all_jerup = (TH1D*)h_calibjet_pt_all->Clone("h_calibjet_pt_record_all_jerup");
  TH1D* h_calibjet_pt_all_jerdown = (TH1D*)h_calibjet_pt_all->Clone("h_calibjet_pt_record_all_jerdown");
  TH1D* h_calibjet_pt_all_jetup = (TH1D*)f_data->Get("h_calibjet_pt_record_all_jetup");
  TH1D* h_calibjet_pt_all_jetdown = (TH1D*)f_data->Get("h_calibjet_pt_record_all_jetdown");

  TH1D* h_calibjet_pt_zvertex30 = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex30");
  TH1D* h_calibjet_pt_zvertex30_jesup = (TH1D*)h_calibjet_pt_zvertex30->Clone("h_calibjet_pt_record_zvertex30_jesup");
  TH1D* h_calibjet_pt_zvertex30_jesdown = (TH1D*)h_calibjet_pt_zvertex30->Clone("h_calibjet_pt_record_zvertex30_jesdown");
  TH1D* h_calibjet_pt_zvertex30_jerup = (TH1D*)h_calibjet_pt_zvertex30->Clone("h_calibjet_pt_record_zvertex30_jerup");
  TH1D* h_calibjet_pt_zvertex30_jerdown = (TH1D*)h_calibjet_pt_zvertex30->Clone("h_calibjet_pt_record_zvertex30_jerdown");
  TH1D* h_calibjet_pt_zvertex30_jetup = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex30_jetup");
  TH1D* h_calibjet_pt_zvertex30_jetdown = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex30_jetdown");
  TH1D* h_calibjet_pt_zvertex30_mbdup = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex30_mbdup");
  TH1D* h_calibjet_pt_zvertex30_mbddown = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex30_mbddown");

  TH1D* h_calibjet_pt_zvertex60 = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex60");
  TH1D* h_calibjet_pt_zvertex60_jesup = (TH1D*)h_calibjet_pt_zvertex60->Clone("h_calibjet_pt_record_zvertex60_jesup");
  TH1D* h_calibjet_pt_zvertex60_jesdown = (TH1D*)h_calibjet_pt_zvertex60->Clone("h_calibjet_pt_record_zvertex60_jesdown");
  TH1D* h_calibjet_pt_zvertex60_jerup = (TH1D*)h_calibjet_pt_zvertex60->Clone("h_calibjet_pt_record_zvertex60_jerup");
  TH1D* h_calibjet_pt_zvertex60_jerdown = (TH1D*)h_calibjet_pt_zvertex60->Clone("h_calibjet_pt_record_zvertex60_jerdown");
  TH1D* h_calibjet_pt_zvertex60_jetup = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex60_jetup");
  TH1D* h_calibjet_pt_zvertex60_jetdown = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex60_jetdown");
  TH1D* h_calibjet_pt_zvertex60_mbdup = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex60_mbdup");
  TH1D* h_calibjet_pt_zvertex60_mbddown = (TH1D*)f_data->Get("h_calibjet_pt_record_zvertex60_mbddown");

  // Read sim histograms for reweighting
  TH1D* h_measure_unweighted_all = (TH1D*)f_sim->Get("h_measure_unweighted_all");
  TH1D* h_measure_unweighted_all_jesup = (TH1D*)f_sim->Get("h_measure_unweighted_all_jesup");
  TH1D* h_measure_unweighted_all_jesdown = (TH1D*)f_sim->Get("h_measure_unweighted_all_jesdown");
  TH1D* h_measure_unweighted_all_jerup = (TH1D*)f_sim->Get("h_measure_unweighted_all_jerup");
  TH1D* h_measure_unweighted_all_jerdown = (TH1D*)f_sim->Get("h_measure_unweighted_all_jerdown");
  TH1D* h_measure_unweighted_all_jetup = (TH1D*)f_sim->Get("h_measure_unweighted_all_jetup");
  TH1D* h_measure_unweighted_all_jetdown = (TH1D*)f_sim->Get("h_measure_unweighted_all_jetdown");

  TH1D* h_measure_unweighted_zvertex30 = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30");
  TH1D* h_measure_unweighted_zvertex30_jesup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_jesup");
  TH1D* h_measure_unweighted_zvertex30_jesdown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_jesdown");
  TH1D* h_measure_unweighted_zvertex30_jerup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_jerup");
  TH1D* h_measure_unweighted_zvertex30_jerdown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_jerdown");
  TH1D* h_measure_unweighted_zvertex30_jetup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_jetup");
  TH1D* h_measure_unweighted_zvertex30_jetdown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_jetdown");
  TH1D* h_measure_unweighted_zvertex30_mbdup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_mbdup");
  TH1D* h_measure_unweighted_zvertex30_mbddown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex30_mbddown");

  TH1D* h_measure_unweighted_zvertex60 = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60");
  TH1D* h_measure_unweighted_zvertex60_jesup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_jesup");
  TH1D* h_measure_unweighted_zvertex60_jesdown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_jesdown");
  TH1D* h_measure_unweighted_zvertex60_jerup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_jerup");
  TH1D* h_measure_unweighted_zvertex60_jerdown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_jerdown");
  TH1D* h_measure_unweighted_zvertex60_jetup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_jetup");
  TH1D* h_measure_unweighted_zvertex60_jetdown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_jetdown");
  TH1D* h_measure_unweighted_zvertex60_mbdup = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_mbdup");
  TH1D* h_measure_unweighted_zvertex60_mbddown = (TH1D*)f_sim->Get("h_measure_unweighted_zvertex60_mbddown");

  // Form reweighting histograms
  TH1D* h_reweight_all; get_reweight_hist(f_out, h_reweight_all, "h_reweight_all", h_calibjet_pt_all, h_measure_unweighted_all);
  TH1D* h_reweight_all_jesup; get_reweight_hist(f_out, h_reweight_all_jesup, "h_reweight_all_jesup", h_calibjet_pt_all_jesup, h_measure_unweighted_all_jesup);
  TH1D* h_reweight_all_jesdown; get_reweight_hist(f_out, h_reweight_all_jesdown, "h_reweight_all_jesdown", h_calibjet_pt_all_jesdown, h_measure_unweighted_all_jesdown);
  TH1D* h_reweight_all_jerup; get_reweight_hist(f_out, h_reweight_all_jerup, "h_reweight_all_jerup", h_calibjet_pt_all_jerup, h_measure_unweighted_all_jerup);
  TH1D* h_reweight_all_jerdown; get_reweight_hist(f_out, h_reweight_all_jerdown, "h_reweight_all_jerdown", h_calibjet_pt_all_jerdown, h_measure_unweighted_all_jerdown);
  TH1D* h_reweight_all_jetup; get_reweight_hist(f_out, h_reweight_all_jetup, "h_reweight_all_jetup", h_calibjet_pt_all_jetup, h_measure_unweighted_all_jetup);
  TH1D* h_reweight_all_jetdown; get_reweight_hist(f_out, h_reweight_all_jetdown, "h_reweight_all_jetdown", h_calibjet_pt_all_jetdown, h_measure_unweighted_all_jetdown);

  TH1D* h_reweight_zvertex30; get_reweight_hist(f_out, h_reweight_zvertex30, "h_reweight_zvertex30", h_calibjet_pt_zvertex30, h_measure_unweighted_zvertex30);
  TH1D* h_reweight_zvertex30_jesup; get_reweight_hist(f_out, h_reweight_zvertex30_jesup, "h_reweight_zvertex30_jesup", h_calibjet_pt_zvertex30_jesup, h_measure_unweighted_zvertex30_jesup);
  TH1D* h_reweight_zvertex30_jesdown; get_reweight_hist(f_out, h_reweight_zvertex30_jesdown, "h_reweight_zvertex30_jesdown", h_calibjet_pt_zvertex30_jesdown, h_measure_unweighted_zvertex30_jesdown);
  TH1D* h_reweight_zvertex30_jerup; get_reweight_hist(f_out, h_reweight_zvertex30_jerup, "h_reweight_zvertex30_jerup", h_calibjet_pt_zvertex30_jerup, h_measure_unweighted_zvertex30_jerup);
  TH1D* h_reweight_zvertex30_jerdown; get_reweight_hist(f_out, h_reweight_zvertex30_jerdown, "h_reweight_zvertex30_jerdown", h_calibjet_pt_zvertex30_jerdown, h_measure_unweighted_zvertex30_jerdown);
  TH1D* h_reweight_zvertex30_jetup; get_reweight_hist(f_out, h_reweight_zvertex30_jetup, "h_reweight_zvertex30_jetup", h_calibjet_pt_zvertex30_jetup, h_measure_unweighted_zvertex30_jetup);
  TH1D* h_reweight_zvertex30_jetdown; get_reweight_hist(f_out, h_reweight_zvertex30_jetdown, "h_reweight_zvertex30_jetdown", h_calibjet_pt_zvertex30_jetdown, h_measure_unweighted_zvertex30_jetdown);
  TH1D* h_reweight_zvertex30_mbdup; get_reweight_hist(f_out, h_reweight_zvertex30_mbdup, "h_reweight_zvertex30_mbdup", h_calibjet_pt_zvertex30_mbdup, h_measure_unweighted_zvertex30_mbdup);
  TH1D* h_reweight_zvertex30_mbddown; get_reweight_hist(f_out, h_reweight_zvertex30_mbddown, "h_reweight_zvertex30_mbddown", h_calibjet_pt_zvertex30_mbddown, h_measure_unweighted_zvertex30_mbddown);

  TH1D* h_reweight_zvertex60; get_reweight_hist(f_out, h_reweight_zvertex60, "h_reweight_zvertex60", h_calibjet_pt_zvertex60, h_measure_unweighted_zvertex60);
  TH1D* h_reweight_zvertex60_jesup; get_reweight_hist(f_out, h_reweight_zvertex60_jesup, "h_reweight_zvertex60_jesup", h_calibjet_pt_zvertex60_jesup, h_measure_unweighted_zvertex60_jesup);
  TH1D* h_reweight_zvertex60_jesdown; get_reweight_hist(f_out, h_reweight_zvertex60_jesdown, "h_reweight_zvertex60_jesdown", h_calibjet_pt_zvertex60_jesdown, h_measure_unweighted_zvertex60_jesdown);
  TH1D* h_reweight_zvertex60_jerup; get_reweight_hist(f_out, h_reweight_zvertex60_jerup, "h_reweight_zvertex60_jerup", h_calibjet_pt_zvertex60_jerup, h_measure_unweighted_zvertex60_jerup);
  TH1D* h_reweight_zvertex60_jerdown; get_reweight_hist(f_out, h_reweight_zvertex60_jerdown, "h_reweight_zvertex60_jerdown", h_calibjet_pt_zvertex60_jerdown, h_measure_unweighted_zvertex60_jerdown);
  TH1D* h_reweight_zvertex60_jetup; get_reweight_hist(f_out, h_reweight_zvertex60_jetup, "h_reweight_zvertex60_jetup", h_calibjet_pt_zvertex60_jetup, h_measure_unweighted_zvertex60_jetup);
  TH1D* h_reweight_zvertex60_jetdown; get_reweight_hist(f_out, h_reweight_zvertex60_jetdown, "h_reweight_zvertex60_jetdown", h_calibjet_pt_zvertex60_jetdown, h_measure_unweighted_zvertex60_jetdown);
  TH1D* h_reweight_zvertex60_mbdup; get_reweight_hist(f_out, h_reweight_zvertex60_mbdup, "h_reweight_zvertex60_mbdup", h_calibjet_pt_zvertex60_mbdup, h_measure_unweighted_zvertex60_mbdup);
  TH1D* h_reweight_zvertex60_mbddown; get_reweight_hist(f_out, h_reweight_zvertex60_mbddown, "h_reweight_zvertex60_mbddown", h_calibjet_pt_zvertex60_mbddown, h_measure_unweighted_zvertex60_mbddown);
}
