TH1D* get_h1(TFile* f, const std::string& name)
{
  TH1D* h = (TH1D*)f->Get(name.c_str());
  if (!h) std::cerr << "Missing: " << name << std::endl;
  return h;
}

TH2D* get_h2(TFile* f, const std::string& name)
{
  TH2D* h = (TH2D*)f->Get(name.c_str());
  if (!h) std::cerr << "Missing: " << name << std::endl;
  return h;
}

void do_unfolding_iteration_check(TH2D* h_respmatrix, TH1D* h_spectrum, string figure_name) {
  TH1D* h_meas = (TH1D*)h_respmatrix->ProjectionX("h_meas");
  TH1D* h_truth = (TH1D*)h_respmatrix->ProjectionY("h_truth");
  const int niter = 20;
  TH1D* h_unfolded[niter];
  for (int i = 0; i < niter; ++i) {
    RooUnfoldResponse* response = new RooUnfoldResponse(h_meas, h_truth, h_respmatrix, "response", "");
    RooUnfoldBayes unfold(response, h_spectrum, i + 1);
    h_unfolded[i] = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
    h_unfolded[i]->SetName(Form("h_unfolded_iter_%d", i + 1));
    delete response;
  }

  TGraph* g_var_iter = new TGraph(niter-1);
  TGraph* g_var_error = new TGraph(niter-1);
  TGraph* g_var_total = new TGraph(niter-1);
  for (int i = 1; i < niter; ++i) {
    double var_iter = 0;
    double var_error = 0;
    double var_total = 0;
    for (int j = 3; j <= h_unfolded[i]->GetNbinsX()-1; ++j) {
      double diff = h_unfolded[i]->GetBinContent(j) - h_unfolded[i-1]->GetBinContent(j);
      double cont = (h_unfolded[i]->GetBinContent(j) + h_unfolded[i-1]->GetBinContent(j)) / 2.0;
      var_iter += diff * diff / (double)(cont * cont);
      var_error += h_unfolded[i]->GetBinError(j) * h_unfolded[i]->GetBinError(j) / (double)(h_unfolded[i]->GetBinContent(j) * h_unfolded[i]->GetBinContent(j));
    }
    var_iter = std::sqrt(var_iter / (double)(h_unfolded[i]->GetNbinsX()-3));
    var_error = std::sqrt(var_error / (double)(h_unfolded[i]->GetNbinsX()-3));
    var_total = std::sqrt(var_iter * var_iter + var_error * var_error);

    g_var_iter->SetPoint(i-1, i+1, var_iter);
    g_var_error->SetPoint(i-1, i+1, var_error);
    g_var_total->SetPoint(i-1, i+1, var_total);
  }
  TCanvas* can = new TCanvas("can", "Unfolding Iteration Check", 800, 600);
  g_var_iter->SetMarkerStyle(20);
  g_var_iter->SetMarkerColor(kBlue);
  g_var_iter->SetLineColor(kBlue);
  g_var_error->SetMarkerStyle(20);
  g_var_error->SetMarkerColor(kRed);
  g_var_error->SetLineColor(kRed);
  g_var_total->SetMarkerStyle(20);
  g_var_total->SetMarkerColor(kBlack);
  g_var_total->SetLineColor(kBlack);
  g_var_total->SetMinimum(0);
  g_var_total->Draw("AP");
  g_var_iter->Draw("P same");
  g_var_error->Draw("P same");
  g_var_iter->GetXaxis()->SetTitle("Iteration");
  can->SaveAs((figure_name + ".png").c_str());
}

void get_unfolded_spectrum(TH2D* h_respmatrix, TH1D* h_spectrum, TH1D*& h_unfolded, int niter, std::string hist_name, TH1D* h_eff) {
  TH1D* h_meas = (TH1D*)h_respmatrix->ProjectionX("h_meas");
  TH1D* h_truth = (TH1D*)h_respmatrix->ProjectionY("h_truth");
  RooUnfoldResponse* response = new RooUnfoldResponse(h_meas, h_truth, h_respmatrix, "response", "");
  RooUnfoldBayes unfold(response, h_spectrum, niter);
  h_unfolded = (TH1D*)unfold.Hunfold(RooUnfold::kErrors);
  h_unfolded->SetName(hist_name.c_str());
  h_unfolded->Divide(h_eff);
  delete response;
}


void do_unfolding_iter(int radius_index = 4)
{
  float jet_radius = 0.1 * radius_index;

  TFile *f_out = new TFile(Form("output_unfolded_r0%d.root", radius_index), "RECREATE");
  TFile *f_data = new TFile(Form("output_data_puritycorr_r0%d.root", radius_index), "READ");
  TFile *f_rm   = new TFile(Form("output_sim_r0%d.root", radius_index), "READ");
  TFile *f_eff  = new TFile(Form("output_purityefficiency_r0%d.root", radius_index), "READ");

  // -----------------------------
  // Central categories
  // -----------------------------
  std::vector<std::string> regions = {
    "all",
    "zvertex30",
    "zvertex60"
  };

  std::vector<std::string> variations = {
    "",
    "_jesup", "_jesdown",
    "_jerup", "_jerdown",
    "_jetup", "_jetdown",
    "_mbdup", "_mbddown"
  };

  // Store outputs
  std::vector<TH1D*> outputs;

  // -----------------------------
  // Loop over regions
  // -----------------------------
  for (const auto& region : regions)
  {
    // --- iteration check only for nominal
    TH2D* h_rm_nom = get_h2(f_rm, Form("h_respmatrix_%s", region.c_str()));
    TH1D* h_data_nom = get_h1(f_data, Form("h_calibjet_pt_puritycorr_%s", region.c_str()));

    do_unfolding_iteration_check(
        h_rm_nom,
        h_data_nom,
        Form("unfolding_iteration_check_%s", region.c_str()));

    // -----------------------------
    // Loop over variations
    // -----------------------------
    for (const auto& var : variations)
    {
      std::string name_suffix = region + var;

      TH1D* h_data = get_h1(f_data,
          Form("h_calibjet_pt_puritycorr_%s%s", region.c_str(), var.c_str()));

      TH2D* h_rm = nullptr;

      // special case: some variations don't exist → clone nominal
      if (var == "_jetup" || var == "_jetdown" ||
          var == "_mbdup" || var == "_mbddown")
      {
        h_rm = (TH2D*)h_rm_nom->Clone(Form("h_respmatrix_%s%s", region.c_str(), var.c_str()));
      }
      else
      {
        h_rm = get_h2(f_rm,
            Form("h_respmatrix_%s%s", region.c_str(), var.c_str()));
      }

      TH1D* h_eff = get_h1(f_eff,
          Form("h_efficiency_%s%s", region.c_str(), var.c_str()));

      if (!h_data || !h_rm || !h_eff) continue;

      // --- nominal unfolding
      TH1D* h_unfold = nullptr;
      get_unfolded_spectrum(
          h_rm,
          h_data,
          h_unfold,
          1,
          Form("h_unfold_%s", name_suffix.c_str()),
          h_eff);

      outputs.push_back(h_unfold);

      // --- unfolding uncertainty (niter=2, only for nominal)
      if (var == "")
      {
        TH1D* h_unfold_unc = nullptr;
        get_unfolded_spectrum(
            h_rm,
            h_data,
            h_unfold_unc,
            2,
            Form("h_unfold_%s_unfoldunc", region.c_str()),
            h_eff);

        outputs.push_back(h_unfold_unc);
      }
    }
  }

  // -----------------------------
  // Write output
  // -----------------------------
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();

  for (auto* h : outputs)
    if (h) h->Write();

  f_out->Close();
  std::cout << "All done!" << std::endl;
}
