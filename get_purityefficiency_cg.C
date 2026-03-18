#include <iostream>
#include <vector>
#include <string>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"
#include "unfold_Def.h"

// ---------------------------------------------------------
// Utility: Fetch histogram safely
// ---------------------------------------------------------
template <typename T>
T* get_histogram(TFile* file, const std::string& name)
{
  T* h = dynamic_cast<T*>(file->Get(name.c_str()));
  if (!h)
  {
    std::cerr << "Error: missing histogram " << name << std::endl;
  }
  return h;
}

// ---------------------------------------------------------
// Compute purity and efficiency from response matrix
// ---------------------------------------------------------
void compute_purity_efficiency(
    TH1D*& h_purity,
    TH1D*& h_efficiency,
    const std::string& suffix,
    TFile* f_in,
    TFile* f_out,
    int jet_radius_index)
{
  // --- Load input histograms ---
  TH1D* h_truth       = get_histogram<TH1D>(f_in, "h_truth" + suffix);
  TH1D* h_measure     = get_histogram<TH1D>(f_in, "h_measure" + suffix);
  TH2D* h_response    = get_histogram<TH2D>(f_in, "h_respmatrix" + suffix);
  TH1D* h_fake        = get_histogram<TH1D>(f_in, "h_fake" + suffix);
  TH1D* h_miss        = get_histogram<TH1D>(f_in, "h_miss" + suffix);

  if (!h_truth || !h_measure || !h_response || !h_fake || !h_miss)
  {
    std::cerr << "Error: missing inputs for suffix " << suffix << std::endl;
    return;
  }

  // -------------------------------------------------------
  // Purity = matched reco / total reco
  // -------------------------------------------------------
  h_purity = (TH1D*)h_response->ProjectionX(("h_purity" + suffix).c_str());
  h_purity->Divide(h_purity, h_measure, 1, 1, "B");

  // -------------------------------------------------------
  // Efficiency = (truth - miss) / truth
  // -------------------------------------------------------
  h_efficiency = (TH1D*)h_truth->Clone(("h_efficiency" + suffix).c_str());
  h_efficiency->Add(h_miss, -1);
  h_efficiency->Divide(h_efficiency, h_truth, 1, 1, "B");

  // -------------------------------------------------------
  // Sanity check: purity consistency
  // -------------------------------------------------------
  TH1D* h_check = (TH1D*)h_response->ProjectionX(("h_purity" + suffix + "_check").c_str());
  h_check->Add(h_fake);

  for (int i = 1; i <= h_check->GetNbinsX(); ++i)
  {
    double meas = h_measure->GetBinContent(i);
    if (meas == 0) continue;

    double ratio = h_check->GetBinContent(i) / meas;
    if (ratio != 1.0)
    {
      std::cout << "[Purity check failed] "
                << suffix << " bin " << i
                << " ratio=" << ratio
                << " check=" << h_check->GetBinContent(i)
                << " meas=" << meas << std::endl;
    }
  }

  // -------------------------------------------------------
  // Write outputs
  // -------------------------------------------------------
  f_out->cd();
  h_purity->Write();
  h_efficiency->Write();

  // -------------------------------------------------------
  // Plotting
  // -------------------------------------------------------
  TH1D* h_purity_proj = (TH1D*)h_response->ProjectionX("h_purity_proj");
  TH1D* h_eff_proj    = (TH1D*)h_truth->Clone("h_efficiency_proj");
  h_eff_proj->Add(h_miss, -1);

  std::vector<TH1F*> h_input;
  std::vector<int> color = {kBlack, kRed};
  std::vector<int> marker = {20, 20};

  std::vector<std::string> text;
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} R = 0." + std::to_string(jet_radius_index));
  text.push_back("|#eta^{jet}| < 0." + std::to_string(11 - jet_radius_index));

  if (suffix.find("zvertex30") != std::string::npos)
    text.push_back("|z-vertex| < 30 cm");
  else if (suffix.find("zvertex60") != std::string::npos)
    text.push_back("|z-vertex| < 60 cm");
  else if (suffix.find("all") != std::string::npos)
    text.push_back("No z-vertex cut");

  // --- Purity plot ---
  h_input = {(TH1F*)h_measure, (TH1F*)h_purity_proj};
  std::vector<std::string> legend = {"Reco jet spectrum", "Matched reco jet spectrum"};

  draw_1D_multiple_plot_ratio(
      h_input, color, marker,
      false, 10, true,
      true, calibptbins[0], calibptbins[calibnpt], false,
      false, 0, 0.5, true,
      false, 0.75, 1.1,
      true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Purity", 1,
      true, text, 0.45, 0.9, 0.05,
      true, legend, 0.25, 0.2, 0.05,
      Form("figure_purityefficiency/purity%s_r0%d.png", suffix.c_str(), jet_radius_index));

  // --- Efficiency plot ---
  h_input = {(TH1F*)h_truth, (TH1F*)h_eff_proj};
  legend  = {"Truth jet spectrum", "Matched truth jet spectrum"};

  draw_1D_multiple_plot_ratio(
      h_input, color, marker,
      false, 10, true,
      true, truthptbins[0], truthptbins[truthnpt-1], false,
      false, 1e-9, 1e-1, true,
      false, 0., 1.2,
      true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
      true, text, 0.45, 0.9, 0.05,
      true, legend, 0.25, 0.2, 0.05,
      Form("figure_purityefficiency/efficiency%s_r0%d.png", suffix.c_str(), jet_radius_index));
}

// ---------------------------------------------------------
// Clone helper
// ---------------------------------------------------------
void clone_purity_efficiency(
    TH1D*& h_purity,
    TH1D*& h_efficiency,
    const std::string& suffix,
    TH1D* h_purity_src,
    TH1D* h_eff_src,
    TFile* f_out)
{
  h_purity    = (TH1D*)h_purity_src->Clone(("h_purity" + suffix).c_str());
  h_efficiency = (TH1D*)h_eff_src->Clone(("h_efficiency" + suffix).c_str());

  f_out->cd();
  h_purity->Write();
  h_efficiency->Write();
}

// ---------------------------------------------------------
// Apply purity correction
// ---------------------------------------------------------
void apply_purity_correction(
    TH1D*& h_result,
    TH1D* h_purity,
    const std::string& result_suffix,
    const std::string& data_suffix,
    TFile* f_data,
    TFile* f_out)
{
  TH1D* h_input = (TH1D*)f_data->Get(("h_calibjet_pt" + data_suffix).c_str());
  if (!h_input)
  {
    std::cerr << "Missing data histogram for " << data_suffix << std::endl;
    return;
  }

  h_result = (TH1D*)h_input->Clone(("h_calibjet_pt_puritycorr" + result_suffix).c_str());
  h_result->Multiply(h_purity);

  f_out->cd();
  h_result->Write();
}
