#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>

#include "unfold_Def.h"
#include "/sphenix/user/hanpuj/plotstyle/AtlasStyle.C"
#include "/sphenix/user/hanpuj/plotstyle/AtlasUtils.C"

// ---------------------------------------------------------
// Utility: safe histogram fetch
// ---------------------------------------------------------
template <typename T>
T* get_hist(TFile* f, const std::string& name)
{
  T* h = dynamic_cast<T*>(f->Get(name.c_str()));
  if (!h)
    std::cerr << "Missing histogram: " << name << std::endl;
  return h;
}

// ---------------------------------------------------------
// Normalize histogram by bin width (used for spectra)
// ---------------------------------------------------------
void normalize_by_binwidth(TH1* h)
{
  for (int i = 1; i <= h->GetNbinsX(); ++i)
  {
    double width = h->GetBinWidth(i);
    if (width == 0) continue;

    h->SetBinContent(i, h->GetBinContent(i) / width);
    h->SetBinError(i, h->GetBinError(i) / width);
  }
}

// ---------------------------------------------------------
// Plot multiple histograms + ratio
// ---------------------------------------------------------
void draw_ratio_plot(
    const std::vector<TH1F*>& input_hists,
    const std::vector<int>& colors,
    const std::vector<int>& markers,
    bool do_rebin, int rebin_factor,
    bool do_normalize,
    bool set_xrange, float xlow, float xhigh, bool logx,
    bool set_yrange, float ylow, float yhigh, bool logy,
    bool set_ratio_range, float rlow, float rhigh,
    TF1* fit_func,
    const std::string& xtitle,
    const std::string& ytitle,
    const std::string& ratio_title,
    bool use_binomial,
    const std::vector<std::string>& text,
    const std::vector<std::string>& legend,
    const std::string& output_name)
{
  if (input_hists.size() <= 1)
  {
    std::cerr << "Need >=2 histograms for ratio plot\n";
    return;
  }

  // Clone input histograms
  std::vector<TH1F*> h;
  for (size_t i = 0; i < input_hists.size(); ++i)
  {
    h.push_back((TH1F*)input_hists[i]->Clone(Form("h_clone_%zu", i)));

    h.back()->SetMarkerStyle(markers[i]);
    h.back()->SetMarkerColor(colors[i]);
    h.back()->SetLineColor(colors[i]);

    if (do_rebin) h.back()->Rebin(rebin_factor);
    if (do_normalize) normalize_by_binwidth(h.back());
  }

  // Canvas
  TCanvas* c = new TCanvas("c_ratio", "", 800, 963);
  c->Divide(1, 2);

  // ------------------ TOP PAD ------------------
  TPad* pad1 = (TPad*)c->cd(1);
  pad1->SetPad(0, 0.4, 1, 1);
  pad1->SetLeftMargin(0.15);
  pad1->SetBottomMargin(0.035);
  if (logx) pad1->SetLogx();
  if (logy) pad1->SetLogy();

  if (set_xrange) h[0]->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_yrange) h[0]->GetYaxis()->SetRangeUser(ylow, yhigh);

  h[0]->GetXaxis()->SetTitle(xtitle.c_str());
  h[0]->GetYaxis()->SetTitle(ytitle.c_str());
  h[0]->Draw();

  for (size_t i = 1; i < h.size(); ++i)
    h[i]->Draw("same");

  // Text
  for (size_t i = 0; i < text.size(); ++i)
    myText(0.4, 0.9 - i * 0.06, 1, text[i].c_str(), 0.05);

  // Legend
  for (size_t i = 0; i < legend.size(); ++i)
    myMarkerLineText(0.25, 0.2 + i * 0.06, 1,
                     colors[i], markers[i], colors[i],
                     1, legend[i].c_str(), 0.05, true);

  // ------------------ RATIO PAD ------------------
  TPad* pad2 = (TPad*)c->cd(2);
  pad2->SetPad(0, 0, 1, 0.4);
  pad2->SetLeftMargin(0.15);
  pad2->SetBottomMargin(0.25);

  std::vector<TH1F*> h_ratio;

  for (size_t i = 1; i < h.size(); ++i)
  {
    TH1F* r = (TH1F*)h[i]->Clone(Form("h_ratio_%zu", i));
    if (use_binomial)
      r->Divide(h[i], h[0], 1, 1, "B");
    else
      r->Divide(h[0]);

    r->SetMarkerStyle(markers[i]);
    r->SetMarkerColor(colors[i]);
    r->SetLineColor(colors[i]);

    h_ratio.push_back(r);
  }

  h_ratio[0]->GetXaxis()->SetTitle(xtitle.c_str());
  h_ratio[0]->GetYaxis()->SetTitle(ratio_title.c_str());

  if (set_xrange) h_ratio[0]->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_ratio_range) h_ratio[0]->GetYaxis()->SetRangeUser(rlow, rhigh);

  h_ratio[0]->Draw();
  for (size_t i = 1; i < h_ratio.size(); ++i)
    h_ratio[i]->Draw("same");

  if (fit_func)
  {
    fit_func->SetLineColor(kRed);
    fit_func->Draw("same");
  }

  // unity line
  double xmin = set_xrange ? xlow : h[0]->GetXaxis()->GetXmin();
  double xmax = set_xrange ? xhigh : h[0]->GetXaxis()->GetXmax();

  TLine* line = new TLine(xmin, 1, xmax, 1);
  line->SetLineStyle(3);
  line->Draw();

  c->SaveAs(output_name.c_str());

  delete line;
  delete c;
}

void compute_reweight_hist(
    TFile* f_out,
    TH1D*& h_reweight,
    const std::string& name,
    TH1D* h_data,
    TH1D* h_sim)
{
  // --- Normalize inputs ---
  h_data->Rebin(10);
  h_sim->Rebin(10);

  h_data->Scale(1.0 / h_data->Integral());
  h_sim->Scale(1.0 / h_sim->Integral());

  // --- Ratio ---
  h_reweight = (TH1D*)h_data->Clone(name.c_str());
  h_reweight->Divide(h_sim);

  // --- Fit ---
  std::string func_name = "reweightfunc_" + name.substr(std::string("h_reweight_").length());

  TF1* f = new TF1(func_name.c_str(), "[0]*exp(-[1]*x)+[2]", 0, 200);
  f->SetParameters(1., 0.1, 0.5);

  h_reweight->Fit(f, "", "", 15, 70);

  // --- Plot ---
  std::vector<TH1F*> h_input = {(TH1F*)h_sim, (TH1F*)h_data};
  std::vector<int> colors = {kRed, kBlack};
  std::vector<int> markers = {24, 24};

  std::vector<std::string> text = {
      "#bf{#it{sPHENIX}} Internal",
      "Data & PYTHIA8 p+p#sqrt{s} = 200 GeV",
      "anti-k_{t} R = 0.4"};

  std::vector<std::string> legend = {
      "Simulation jet spectrum",
      "Data jet spectrum"};

  draw_ratio_plot(
      h_input, colors, markers,
      false, 10, true,
      true, calibptbins[0], calibptbins[calibnpt], false,
      false, 0, 0.5, true,
      true, 0., 2.,
      f,
      "p_{T}^{jet} [GeV]",
      "Arbitrary Unit",
      "Reweight factor",
      false,
      text,
      legend,
      Form("figure_reweight/%s.png", func_name.c_str()));

  // --- Output ---
  f_out->cd();
  h_reweight->Write();
  f->Write();
}
