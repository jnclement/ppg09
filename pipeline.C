// ========================= helpers.C =========================
// Shared utilities across all macros

#include <map>
#include <vector>
#include <string>
#include <iostream>

void do_normalization(TH1D* h, double lumi, double R)
{
  if (!h) return;
  for (int i = 1; i <= h->GetNbinsX(); ++i)
  {
    double val = h->GetBinContent(i);
    double width = h->GetBinWidth(i);
    if (width > 0 && corr > 0)
      h->SetBinContent(i, val / (lumi * width * (2*(1.1-R))));
  }
}

TH1D* get_statistical_uncertainty(TH1D* h, const std::string& name)
{
  if (!h) return nullptr;
  TH1D* out = (TH1D*)h->Clone(name.c_str());

  for (int i = 1; i <= h->GetNbinsX(); ++i)
  {
    double val = h->GetBinContent(i);
    double err = h->GetBinError(i);
    out->SetBinContent(i, val != 0 ? err / val : 0);
    out->SetBinContent(i, 0);
  }
  return out;
}

void get_uncertainty(
    TH1D*& up,
    TH1D*& down,
    TH1D* nominal,
    TH1D* hup,
    TH1D* hdown,
    const std::string& name)
{
  if (!nominal || !hup || !hdown) return;

  up = (TH1D*)nominal->Clone((name + "_up").c_str());
  down = (TH1D*)nominal->Clone((name + "_down").c_str());
  up->Reset();
  down->Reset();
  for (int i = 1; i <= nominal->GetNbinsX(); ++i)
  {
    double n = nominal->GetBinContent(i);
    double u = hup->GetBinContent(i);
    double d = hdown->GetBinContent(i);

    if (n != 0)
    {
      up->SetBinContent(i, (u - n) / n);
      down->SetBinContent(i, (n - d) / n);
    }
    else
    {
      up->SetBinContent(i, 0);
      down->SetBinContent(i, 0);
    }
  }
}


void get_total_uncertainty(TH1D* &h_total_up, TH1D* &h_total_dn, TH1D* h_stat, std::map<std::string, TH1D*> unc_up, std::map<std::string, TH1D*> unc_dn)
{
  h_total_up = (TH1D*)h_stat->Clone("h_total_up");
  h_total_dn = (TH1D*)h_stat->Clone("h_total_dn");
  h_total_up->Reset();
  h_total_dn->Reset();

  for(int i=1; i<=h_stat->GetNbinsX(); ++i)
    {
      double stat_err = h_stat->GetBinContent(i);
      double sys_up = 0;
      double sys_dn = 0;
      for(auto& [k,v] : unc_up) sys_up += v->GetBinContent(i) * v->GetBinContent(i);
      for(auto& [k,v] : unc_dn) sys_dn += v->GetBinContent(i) * v->GetBinContent(i);
      h_total_up->SetBinContent(i, sqrt(sys_up));
      h_total_dn->SetBinContent(i, sqrt(sys_dn));
      h_total_up->SetBinError(i, 0);
      h_total_dn->SetBinError(i, 0);
    }
}

// ========================= get_finalspectrum_refactored.C =========================

void get_finalspectrum_refactored(TFile* f_spectrum, TFile* f_out, string ztype = "all")
{
  std::vector<std::string> variations = {
    "", "_jesup", "_jesdown", "_jerup", "_jerdown", "_jetup", "_jetdown", "_mbdup", "_mbddown", "_unfoldunc"
  };

  std::map<std::string, TH1D*> h_unfold;

  for (const auto& var : variations)
  {
    if(ztype=="all" && (var=="_mbdup" || var=="_mbddown")) continue;
    std::string name = "h_unfold_" + ztype + var;
    h_unfold[var] = (TH1D*)f_spectrum->Get(name.c_str());

    if (!h_unfold[var])
    {
      std::cerr << "Missing histogram: " << name << std::endl;
      continue;
    }

    do_normalization(h_unfold[var], 1.0, nullptr, 1.0);
  }

  // Statistical uncertainty
  TH1D* h_stat = get_statistical_uncertainty(
      h_unfold[""], "h_unc_stat");

  // Systematics
  struct Sys { std::string name, up, down; };

  std::vector<Sys> systs = {
    {"jes", "_jesup", "_jesdown"},
    {"jer", "_jerup", "_jerdown"},
    {"jet", "_jetup", "_jetdown"},
    {"mbd", "_mbdup", "_mbddown"},
    {"unf", "_unfoldunc", "_unfoldunc"}
  };

  std::map<std::string, TH1D*> unc_up, unc_down;

  for (auto& s : systs)
  {
    if(s.name=="mbd" && ztype=="all") continue;
    get_uncertainty(
      unc_up[s.name], unc_down[s.name],
      h_unfold[""], h_unfold[s.up], h_unfold[s.down],
      "unc_" + s.name);
  }

  f_out->cd();
  h_unfold[""]->Write("final_spectrum");
  h_stat->Write();

  for (auto& [k,v] : unc_up) v->Write();
  for (auto& [k,v] : unc_down) v->Write();

  TH1D* h_total_up, * h_total_dn;

  get_total_uncertainty(h_total_up, h_total_dn, h_stat, unc_up, unc_down);

  h_total_up->Write();
  h_total_dn->Write();
  
}
