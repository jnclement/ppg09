#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
R__LOAD_LIBRARY(libBootstrapGenerator.so)
#include "BootstrapGenerator/BootstrapGenerator.h"
#include "BootstrapGenerator/TH1DBootstrap.h"
#include "BootstrapGenerator/TH2DBootstrap.h"
void trim_1Dhist_meas(TH1D* h, std::string dataset);
void trim_1Dhist_truth(TH1D* h, std::string dataset);
void trim_2Dhist(TH2D* h, std::string dataset);
void combine_truthjet(TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out);
void combine_1Dhist(string histname, string surfix, std::string datatype, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out, bool hasall = true);
void combine_2Dhist(string histname, string surfix, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out, bool hasall = true);
void combine_1Dhist_noall(string histname, string surfix, std::string datatype, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out);
void combine_2Dhist_noall(string histname, string surfix, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out);
void combine_1Dhist_closure(string surfix, string tag, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out);
void combine_2Dhist_closure(string surfix, string tag, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out);

const double mb_cross_section = 4.197e+10;
const double jet5_cross_section = 1.3878e+08;
const double jet12_cross_section = 1.4903e+06;
const double jet40_cross_section = 1.3553e+02;
const double jet20_cross_section = 6.2623e+04;
const double jet30_cross_section = 2.5298e+03;
const double jet50_cross_section = 7.3113;
const double jet60_cross_section = 3.3261e-01;

int get_combinedoutput(int radius_index = 4) {

  TFile* f_MB = new TFile(Form("ana_output/output_mb_r0%d.root", radius_index), "READ");
  TFile* f_Jet5GeV = new TFile(Form("ana_output/output_jet5_r0%d.root", radius_index), "READ");
  TFile* f_Jet12GeV = new TFile(Form("ana_output/output_jet12_r0%d.root", radius_index), "READ");
  TFile* f_Jet40GeV = new TFile(Form("ana_output/output_jet40_r0%d.root", radius_index), "READ");
  TFile* f_Jet20GeV = new TFile(Form("ana_output/output_jet20_r0%d.root", radius_index), "READ");
  TFile* f_Jet30GeV = new TFile(Form("ana_output/output_jet30_r0%d.root", radius_index), "READ");
  TFile* f_Jet50GeV = new TFile(Form("ana_output/output_jet50_r0%d.root", radius_index), "READ");
  TFile* f_Jet60GeV = new TFile(Form("ana_output/output_jet60_r0%d.root", radius_index), "READ");
  if (!f_MB || !f_Jet5GeV || !f_Jet12GeV || !f_Jet40GeV || !f_Jet20GeV || !f_Jet30GeV || !f_Jet50GeV || !f_Jet60GeV) {
    std::cout << "Error: cannot open one or more input files." << std::endl;
    return 1;
  }
  TFile* f_combined = new TFile(Form("output_comb_r0%d.root", radius_index), "RECREATE");

  combine_truthjet(f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  //combine_1Dhist("h_recojet_pt_record_nocut", "", "reco", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  //combine_1Dhist("h_recojet_pt_record", "", "reco", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine nominal histograms
  combine_1Dhist("h_truth", "", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine JES up variation histograms
  combine_1Dhist("h_truth", "_jesup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "_jesup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "_jesup", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "_jesup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "_jesup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "_jesup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "_jesup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "_jesup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine JES down variation histograms
  combine_1Dhist("h_truth", "_jesdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "_jesdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "_jesdown", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "_jesdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "_jesdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "_jesdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "_jesdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "_jesdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine JER up variation histograms
  combine_1Dhist("h_truth", "_jerup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "_jerup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "_jerup", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "_jerup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "_jerup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "_jerup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "_jerup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "_jerup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine JER down variation histograms
  combine_1Dhist("h_truth", "_jerdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "_jerdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "_jerdown", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "_jerdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "_jerdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "_jerdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "_jerdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "_jerdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine jet trigger up variation histograms
  combine_1Dhist("h_truth", "_jetup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "_jetup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "_jetup", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "_jetup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "_jetup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "_jetup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "_jetup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "_jetup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine jet trigger down variation histograms
  combine_1Dhist("h_truth", "_jetdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure", "_jetdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist("h_respmatrix", "_jetdown", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_fake", "_jetdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_miss", "_jetdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_weighted", "_jetdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_matchedtruth_unweighted", "_jetdown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist("h_measure_unweighted", "_jetdown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine mbd trigger up variation histograms
  combine_1Dhist_noall("h_truth", "_mbdup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_measure", "_mbdup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist_noall("h_respmatrix", "_mbdup", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_fake", "_mbdup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_miss", "_mbdup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_matchedtruth_weighted", "_mbdup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_matchedtruth_unweighted", "_mbdup", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_measure_unweighted", "_mbdup", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // Combine mbd trigger down variation histograms
  combine_1Dhist_noall("h_truth", "_mbddown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_measure", "_mbddown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist_noall("h_respmatrix", "_mbddown", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_fake", "_mbddown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_miss", "_mbddown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_matchedtruth_weighted", "_mbddown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_matchedtruth_unweighted", "_mbddown", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_noall("h_measure_unweighted", "_mbddown", "meas", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  // full closure test
  combine_1Dhist_closure("h_fullclosure_", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_fullclosure_", "measure", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist_closure("h_fullclosure_", "respmatrix", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_fullclosure_", "fake", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_fullclosure_", "miss", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  combine_1Dhist_closure("h_halfclosure_", "inputmeasure", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_halfclosure_", "truth", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_halfclosure_", "measure", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_2Dhist_closure("h_halfclosure_", "respmatrix", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_halfclosure_", "fake", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);
  combine_1Dhist_closure("h_halfclosure_", "miss", f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_combined);

  return 0;
}

void combine_truthjet(TFile* f_MB,  TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out) {
  TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = "h_truthjet_pt_record_all";

  TH1D* h_MB_all_forcombine = (TH1D*)f_MB->Get(histname_all.c_str());
  TH1D* h_Jet5GeV_all_forcombine = (TH1D*)f_Jet5GeV->Get(histname_all.c_str());
  TH1D* h_Jet12GeV_all_forcombine = (TH1D*)f_Jet12GeV->Get(histname_all.c_str());
  TH1D* h_Jet40GeV_all_forcombine = (TH1D*)f_Jet40GeV->Get(histname_all.c_str());
  TH1D* h_Jet20GeV_all_forcombine = (TH1D*)f_Jet20GeV->Get(histname_all.c_str());
  TH1D* h_Jet30GeV_all_forcombine = (TH1D*)f_Jet30GeV->Get(histname_all.c_str());
  TH1D* h_Jet50GeV_all_forcombine = (TH1D*)f_Jet50GeV->Get(histname_all.c_str());
  TH1D* h_Jet60GeV_all_forcombine = (TH1D*)f_Jet60GeV->Get(histname_all.c_str());

  TH1D* h_all_combined = (TH1D*)h_MB_all_forcombine->Clone(histname_all.c_str());
  
  h_all_combined->Scale(mb_scale_all);
  h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
  h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
  h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
  h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
  h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
  h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
  h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);

  f_out->cd();
  h_all_combined->Write();
}

void combine_1Dhist(string histname, string surfix, std::string datatype, TFile* f_MB,  TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out, bool hasall) {
  TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname_all = histname + "_all" + surfix;
  string histname_zvertex30 = histname + "_zvertex30" + surfix;
  string histname_zvertex60 = histname + "_zvertex60" + surfix;
  string histname_nozvtx = histname + "_nozvtx_nosmear";
  TH1DBootstrap* h_MB_all_forcombine = (TH1DBootstrap*)f_MB->Get(histname_all.c_str());
  TH1DBootstrap* h_MB_zvertex30_forcombine = (TH1DBootstrap*)f_MB->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_MB_zvertex60_forcombine = (TH1DBootstrap*)f_MB->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_MB_nozvtx_forcombine = (TH1DBootstrap*)f_MB->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet5GeV_all_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet12GeV_all_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet40GeV_all_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet20GeV_all_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet30GeV_all_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet50GeV_all_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname_all.c_str());
  TH1DBootstrap* h_Jet60GeV_all_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname_all.c_str());
  TAxis* ax = h_MB_zvertex60_forcombine->GetNominal()->GetXaxis();
  int nbins = ax->GetNbins();
  const TArrayD* bins = ax->GetXbins();
  const double* edges;
  double xmin, xmax;
  TH1DBootstrap* h_all_combined;
  TH1DBootstrap* h_nozvtx_combined;
  TH1DBootstrap* h_zvertex30_combined;
  TH1DBootstrap* h_zvertex60_combined;
  if(bins->GetSize() > 0)
    {
      edges = bins->GetArray();
      if(hasall) h_all_combined = new TH1DBootstrap(histname_all.c_str(),h_MB_all_forcombine->GetNominal()->GetTitle(), nbins, edges, h_MB_all_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_all_forcombine->GetGenerator());
      if(h_MB_nozvtx_forcombine) h_nozvtx_combined = new TH1DBootstrap(histname_nozvtx.c_str(),h_MB_nozvtx_forcombine->GetTitle(), nbins, edges, h_MB_nozvtx_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_nozvtx_forcombine->GetGenerator());
      h_zvertex30_combined = new TH1DBootstrap(histname_zvertex30.c_str(),h_MB_zvertex30_forcombine->GetTitle(), nbins, edges, h_MB_zvertex30_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_zvertex30_forcombine->GetGenerator());
      h_zvertex60_combined = new TH1DBootstrap(histname_zvertex60.c_str(),h_MB_zvertex60_forcombine->GetTitle(), nbins, edges, h_MB_zvertex60_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }
  else
    {
      xmin = ax->GetXmin();
      xmax = ax->GetXmax();
      if(hasall) h_all_combined = new TH1DBootstrap(histname_all.c_str(),h_MB_all_forcombine->GetTitle(),nbins, xmin, xmax, h_MB_all_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_all_forcombine->GetGenerator());
      if(h_MB_nozvtx_forcombine) h_nozvtx_combined = new TH1DBootstrap(histname_nozvtx.c_str(),h_MB_nozvtx_forcombine->GetTitle(),nbins, xmin, xmax, h_MB_nozvtx_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_nozvtx_forcombine->GetGenerator());
      h_zvertex30_combined = new TH1DBootstrap(histname_zvertex30.c_str(),h_MB_zvertex30_forcombine->GetTitle(),nbins, xmin, xmax, h_MB_zvertex30_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_zvertex30_forcombine->GetGenerator());
      h_zvertex60_combined = new TH1DBootstrap(histname_zvertex60.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nbins, xmin, xmax, h_MB_zvertex60_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }

  if(hasall)
    {
      h_all_combined->Add(h_MB_all_forcombine, mb_scale_all);
      h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
      h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
      h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
      h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
      h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
      h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
      h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);
    }


  TH1DBootstrap* h_Jet5GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_Jet12GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_Jet40GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_Jet20GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_Jet30GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_Jet50GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname_zvertex30.c_str());
  TH1DBootstrap* h_Jet60GeV_zvertex30_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname_zvertex30.c_str());

  h_zvertex30_combined->Add(h_MB_zvertex30_forcombine, mb_scale_all);
  h_zvertex30_combined->Add(h_Jet5GeV_zvertex30_forcombine, jet5_scale_all);
  h_zvertex30_combined->Add(h_Jet12GeV_zvertex30_forcombine, jet12_scale_all);
  h_zvertex30_combined->Add(h_Jet40GeV_zvertex30_forcombine, jet40_scale_all);
  h_zvertex30_combined->Add(h_Jet20GeV_zvertex30_forcombine, jet20_scale_all);
  h_zvertex30_combined->Add(h_Jet30GeV_zvertex30_forcombine, jet30_scale_all);
  h_zvertex30_combined->Add(h_Jet50GeV_zvertex30_forcombine, jet50_scale_all);
  h_zvertex30_combined->Add(h_Jet60GeV_zvertex30_forcombine, jet60_scale_all);


  TH1DBootstrap* h_Jet5GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet12GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet40GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet20GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet30GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet50GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname_zvertex60.c_str());
  TH1DBootstrap* h_Jet60GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname_zvertex60.c_str());

  h_zvertex60_combined->Add(h_MB_zvertex30_forcombine, mb_scale_all);
  h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);



  TH1DBootstrap* h_Jet5GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet12GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet40GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet20GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet30GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet50GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname_nozvtx.c_str());
  TH1DBootstrap* h_Jet60GeV_nozvtx_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname_nozvtx.c_str());
  if(surfix=="" && (histname=="h_fake" || histname=="h_miss" || histname=="h_measure" || histname=="h_truth"))
    {
      h_nozvtx_combined->Add(h_MB_nozvtx_forcombine, mb_scale_all);
      h_nozvtx_combined->Add(h_Jet5GeV_nozvtx_forcombine, jet5_scale_all);
      h_nozvtx_combined->Add(h_Jet12GeV_nozvtx_forcombine, jet12_scale_all);
      h_nozvtx_combined->Add(h_Jet40GeV_nozvtx_forcombine, jet40_scale_all);
      h_nozvtx_combined->Add(h_Jet20GeV_nozvtx_forcombine, jet20_scale_all);
      h_nozvtx_combined->Add(h_Jet30GeV_nozvtx_forcombine, jet30_scale_all);
      h_nozvtx_combined->Add(h_Jet50GeV_nozvtx_forcombine, jet50_scale_all);
      h_nozvtx_combined->Add(h_Jet60GeV_nozvtx_forcombine, jet60_scale_all);
    }

  f_out->cd();
  if(hasall) h_all_combined->Write();
  h_zvertex30_combined->Write();
  h_zvertex60_combined->Write();
  if(surfix=="" && (histname=="h_fake" || histname=="h_miss" || histname=="h_measure" || histname=="h_truth")) h_nozvtx_combined->Write();
}

void combine_2Dhist(string histname, string surfix, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out, bool hasall) {
  TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;
  string histname_all = histname + "_all" + surfix;
  string histname_zvertex30 = histname + "_zvertex30" + surfix;
  string histname_zvertex60 = histname + "_zvertex60" + surfix;
  string histname_nozvtx = histname + "_nozvtx" + surfix + "_nosmear";

  string histname_nentry_all = histname + "_nentry_all" + surfix;
  string histname_nentry_zvertex30 = histname + "_nentry_zvertex30" + surfix;
  string histname_nentry_zvertex60 = histname + "_nentry_zvertex60" + surfix;

  TH2DBootstrap* h_MB_all_forcombine = (TH2DBootstrap*)f_MB->Get(histname_all.c_str());
  TH2DBootstrap* h_MB_zvertex30_forcombine = (TH2DBootstrap*)f_MB->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_MB_zvertex60_forcombine = (TH2DBootstrap*)f_MB->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_MB_nozvtx_forcombine = (TH2DBootstrap*)f_MB->Get(histname_nozvtx.c_str());

  TAxis* xax = h_MB_zvertex60_forcombine->GetNominal()->GetXaxis();
  TAxis* yax = h_MB_zvertex60_forcombine->GetNominal()->GetYaxis();
  int nxbins = xax->GetNbins();
  int nybins = yax->GetNbins();

  const TArrayD* xbins = xax->GetXbins();
  const TArrayD* ybins = yax->GetXbins();
  const double* xedges;
  const double* yedges;
  double xmin, xmax, ymin, ymax;
  TH2DBootstrap* h_all_combined, *h_nozvtx_combined, *h_zvertex30_combined, *h_zvertex60_combined, *h_nentry_all_combined, *h_nentry_zvertex30_combined, *h_nentry_zvertex60_combined, *h_nentry_nozvtx_combined;

  if(xbins->GetSize() > 0)
    {
      xedges = xbins->GetArray();
      yedges = ybins->GetArray();
      if(hasall) h_all_combined = new TH2DBootstrap(histname_all.c_str(),h_MB_all_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_all_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_all_forcombine->GetGenerator());
      if(h_MB_nozvtx_forcombine) h_nozvtx_combined = new TH2DBootstrap(histname_nozvtx.c_str(),h_MB_nozvtx_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_nozvtx_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_nozvtx_forcombine->GetGenerator());
      h_zvertex30_combined = new TH2DBootstrap(histname_zvertex30.c_str(),h_MB_zvertex30_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_zvertex30_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex30_forcombine->GetGenerator());
      h_zvertex60_combined = new TH2DBootstrap(histname_zvertex60.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_zvertex60_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());


      if(hasall) h_nentry_all_combined = new TH2DBootstrap(histname_nentry_all.c_str(),h_MB_all_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_all_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_all_forcombine->GetGenerator());
      if(h_MB_nozvtx_forcombine) h_nentry_nozvtx_combined = new TH2DBootstrap(histname_nozvtx.c_str(),h_MB_nozvtx_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_nozvtx_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_nozvtx_forcombine->GetGenerator());
      h_nentry_zvertex30_combined = new TH2DBootstrap(histname_nentry_zvertex30.c_str(),h_MB_zvertex30_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_zvertex30_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex30_forcombine->GetGenerator());
      h_nentry_zvertex60_combined = new TH2DBootstrap(histname_nentry_zvertex60.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_zvertex60_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }
  else
    {
      xmin = xax->GetXmin();
      xmax = xax->GetXmax();
      ymin = yax->GetXmin();
      ymax = yax->GetXmax();

      if(hasall) h_all_combined = new TH2DBootstrap(histname_all.c_str(),h_MB_all_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_all_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_all_forcombine->GetGenerator());
      if(h_MB_nozvtx_forcombine) h_nozvtx_combined = new TH2DBootstrap(histname_nozvtx.c_str(),h_MB_nozvtx_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_nozvtx_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_nozvtx_forcombine->GetGenerator());
      h_zvertex30_combined = new TH2DBootstrap(histname_zvertex30.c_str(),h_MB_zvertex30_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_zvertex30_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex30_forcombine->GetGenerator());
      h_zvertex60_combined = new TH2DBootstrap(histname_zvertex60.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_zvertex60_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());


      if(hasall) h_nentry_all_combined = new TH2DBootstrap(histname_nentry_all.c_str(),h_MB_all_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_all_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_all_forcombine->GetGenerator());
      h_nentry_zvertex30_combined = new TH2DBootstrap(histname_nentry_zvertex30.c_str(),h_MB_zvertex30_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_zvertex30_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex30_forcombine->GetGenerator());
      h_nentry_zvertex60_combined = new TH2DBootstrap(histname_nentry_zvertex60.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_zvertex60_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }
  
  
  
  TH2DBootstrap* h_Jet5GeV_all_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet12GeV_all_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet40GeV_all_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet20GeV_all_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet30GeV_all_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet50GeV_all_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname_all.c_str());
  TH2DBootstrap* h_Jet60GeV_all_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname_all.c_str());
  if(hasall)
    {
      h_all_combined->Add(h_MB_all_forcombine,mb_scale_all);
      h_all_combined->Add(h_Jet5GeV_all_forcombine, jet5_scale_all);
      h_all_combined->Add(h_Jet12GeV_all_forcombine, jet12_scale_all);
      h_all_combined->Add(h_Jet40GeV_all_forcombine, jet40_scale_all);
      h_all_combined->Add(h_Jet20GeV_all_forcombine, jet20_scale_all);
      h_all_combined->Add(h_Jet30GeV_all_forcombine, jet30_scale_all);
      h_all_combined->Add(h_Jet50GeV_all_forcombine, jet50_scale_all);
      h_all_combined->Add(h_Jet60GeV_all_forcombine, jet60_scale_all);
      h_nentry_all_combined->Add(h_MB_all_forcombine);
      h_nentry_all_combined->Add(h_Jet5GeV_all_forcombine);
      h_nentry_all_combined->Add(h_Jet12GeV_all_forcombine);
      h_nentry_all_combined->Add(h_Jet40GeV_all_forcombine);
      h_nentry_all_combined->Add(h_Jet20GeV_all_forcombine);
      h_nentry_all_combined->Add(h_Jet30GeV_all_forcombine);
      h_nentry_all_combined->Add(h_Jet50GeV_all_forcombine);
      h_nentry_all_combined->Add(h_Jet60GeV_all_forcombine);
    }


  TH2DBootstrap* h_Jet5GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_Jet12GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_Jet40GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_Jet20GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_Jet30GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_Jet50GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname_zvertex30.c_str());
  TH2DBootstrap* h_Jet60GeV_zvertex30_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname_zvertex30.c_str());

  h_zvertex30_combined->Add(h_MB_zvertex30_forcombine,mb_scale_all);
  h_zvertex30_combined->Add(h_Jet5GeV_zvertex30_forcombine, jet5_scale_all);
  h_zvertex30_combined->Add(h_Jet12GeV_zvertex30_forcombine, jet12_scale_all);
  h_zvertex30_combined->Add(h_Jet40GeV_zvertex30_forcombine, jet40_scale_all);
  h_zvertex30_combined->Add(h_Jet20GeV_zvertex30_forcombine, jet20_scale_all);
  h_zvertex30_combined->Add(h_Jet30GeV_zvertex30_forcombine, jet30_scale_all);
  h_zvertex30_combined->Add(h_Jet50GeV_zvertex30_forcombine, jet50_scale_all);
  h_zvertex30_combined->Add(h_Jet60GeV_zvertex30_forcombine, jet60_scale_all);
  h_nentry_zvertex30_combined->Add(h_MB_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet5GeV_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet12GeV_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet40GeV_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet20GeV_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet30GeV_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet50GeV_zvertex30_forcombine);
  h_nentry_zvertex30_combined->Add(h_Jet60GeV_zvertex30_forcombine);


  TH2DBootstrap* h_Jet5GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet12GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet40GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet20GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet30GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet50GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname_zvertex60.c_str());
  TH2DBootstrap* h_Jet60GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname_zvertex60.c_str());

  h_zvertex60_combined->Add(h_MB_zvertex60_forcombine,mb_scale_all);
  h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);
  h_nentry_zvertex60_combined->Add(h_MB_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine);
  h_nentry_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine);


  TH2DBootstrap* h_Jet5GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname_nozvtx.c_str());
  TH2DBootstrap* h_Jet12GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname_nozvtx.c_str());
  TH2DBootstrap* h_Jet40GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname_nozvtx.c_str());
  TH2DBootstrap* h_Jet20GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname_nozvtx.c_str());
  TH2DBootstrap* h_Jet30GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname_nozvtx.c_str());
  TH2DBootstrap* h_Jet50GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname_nozvtx.c_str());
  TH2DBootstrap* h_Jet60GeV_nozvtx_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname_nozvtx.c_str());
  if(surfix=="")
    {
      
      h_nozvtx_combined->Add(h_MB_nozvtx_forcombine, mb_scale_all);
      h_nozvtx_combined->Add(h_Jet5GeV_nozvtx_forcombine, jet5_scale_all);
      h_nozvtx_combined->Add(h_Jet12GeV_nozvtx_forcombine, jet12_scale_all);
      h_nozvtx_combined->Add(h_Jet40GeV_nozvtx_forcombine, jet40_scale_all);
      h_nozvtx_combined->Add(h_Jet20GeV_nozvtx_forcombine, jet20_scale_all);
      h_nozvtx_combined->Add(h_Jet30GeV_nozvtx_forcombine, jet30_scale_all);
      h_nozvtx_combined->Add(h_Jet50GeV_nozvtx_forcombine, jet50_scale_all);
      h_nozvtx_combined->Add(h_Jet60GeV_nozvtx_forcombine, jet60_scale_all);
    }


  f_out->cd();
  if(hasall) h_all_combined->Write();
  h_zvertex30_combined->Write();
  h_zvertex60_combined->Write();
  if(surfix=="") h_nozvtx_combined->Write();
  if(hasall) h_nentry_all_combined->Write();
  h_nentry_zvertex30_combined->Write();
  h_nentry_zvertex60_combined->Write();
}

void combine_1Dhist_noall(string histname, string surfix, std::string datatype, TFile* f_MB,  TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out) {

  combine_1Dhist(histname,surfix,datatype,f_MB,f_Jet5GeV,f_Jet12GeV,f_Jet40GeV,f_Jet20GeV,f_Jet30GeV,f_Jet50GeV,f_Jet60GeV,f_out,false);
}

void combine_2Dhist_noall(string histname, string surfix, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out)
{
  combine_2Dhist(histname, surfix, f_MB, f_Jet5GeV, f_Jet12GeV, f_Jet40GeV, f_Jet20GeV, f_Jet30GeV, f_Jet50GeV, f_Jet60GeV, f_out, false);
}


void combine_1Dhist_closure(string surfix, string tag, TFile* f_MB, TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out) {
  TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname = surfix + tag + "_zvertex60";

  TH1DBootstrap* h_MB_zvertex60_forcombine = (TH1DBootstrap*)f_MB->Get(histname.c_str());

  TAxis* ax = h_MB_zvertex60_forcombine->GetNominal()->GetXaxis();
  int nbins = ax->GetNbins();
  const TArrayD* bins = ax->GetXbins();
  const double* edges;
  double xmin, xmax;
  TH1DBootstrap* h_zvertex60_combined;
  if(bins->GetSize() > 0)
    {
      edges = bins->GetArray();
      h_zvertex60_combined = new TH1DBootstrap(histname.c_str(),h_MB_zvertex60_forcombine->GetTitle(), nbins, edges, h_MB_zvertex60_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }
  else
    {
      xmin = ax->GetXmin();
      xmax = ax->GetXmax();
      h_zvertex60_combined = new TH1DBootstrap(histname.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nbins, xmin, xmax, h_MB_zvertex60_forcombine->GetNReplica(), (BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }

  TH1DBootstrap* h_Jet5GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet5GeV->Get(histname.c_str());
  TH1DBootstrap* h_Jet12GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet12GeV->Get(histname.c_str());
  TH1DBootstrap* h_Jet40GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet40GeV->Get(histname.c_str());
  TH1DBootstrap* h_Jet20GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet20GeV->Get(histname.c_str());
  TH1DBootstrap* h_Jet30GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet30GeV->Get(histname.c_str());
  TH1DBootstrap* h_Jet50GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet50GeV->Get(histname.c_str());
  TH1DBootstrap* h_Jet60GeV_zvertex60_forcombine = (TH1DBootstrap*)f_Jet60GeV->Get(histname.c_str());
  h_zvertex60_combined->Add(h_MB_zvertex60_forcombine, mb_scale_all);
  h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_zvertex60_combined->Write();
}

void combine_2Dhist_closure(string surfix, string tag, TFile* f_MB,  TFile* f_Jet5GeV, TFile* f_Jet12GeV, TFile* f_Jet40GeV, TFile* f_Jet20GeV, TFile* f_Jet30GeV, TFile* f_Jet50GeV, TFile* f_Jet60GeV, TFile* f_out) {
  TH1D *h_MB_event_all = (TH1D*)f_MB->Get("h_event_all"); int mb_event_all = h_MB_event_all->GetBinContent(1); double mb_scale_all = mb_cross_section / (double)mb_event_all;
  TH1D *h_Jet5GeV_event_all = (TH1D*)f_Jet5GeV->Get("h_event_all"); int jet5_event_all = h_Jet5GeV_event_all->GetBinContent(1); double jet5_scale_all = jet5_cross_section / (double)jet5_event_all;
  TH1D *h_Jet12GeV_event_all = (TH1D*)f_Jet12GeV->Get("h_event_all"); int jet12_event_all = h_Jet12GeV_event_all->GetBinContent(1); double jet12_scale_all = jet12_cross_section / (double)jet12_event_all;
  TH1D *h_Jet40GeV_event_all = (TH1D*)f_Jet40GeV->Get("h_event_all"); int jet40_event_all = h_Jet40GeV_event_all->GetBinContent(1); double jet40_scale_all = jet40_cross_section / (double)jet40_event_all;
  TH1D *h_Jet20GeV_event_all = (TH1D*)f_Jet20GeV->Get("h_event_all"); int jet20_event_all = h_Jet20GeV_event_all->GetBinContent(1); double jet20_scale_all = jet20_cross_section / (double)jet20_event_all;
  TH1D *h_Jet30GeV_event_all = (TH1D*)f_Jet30GeV->Get("h_event_all"); int jet30_event_all = h_Jet30GeV_event_all->GetBinContent(1); double jet30_scale_all = jet30_cross_section / (double)jet30_event_all;
  TH1D *h_Jet50GeV_event_all = (TH1D*)f_Jet50GeV->Get("h_event_all"); int jet50_event_all = h_Jet50GeV_event_all->GetBinContent(1); double jet50_scale_all = jet50_cross_section / (double)jet50_event_all;
  TH1D *h_Jet60GeV_event_all = (TH1D*)f_Jet60GeV->Get("h_event_all"); int jet60_event_all = h_Jet60GeV_event_all->GetBinContent(1); double jet60_scale_all = jet60_cross_section / (double)jet60_event_all;

  string histname = surfix + tag + "_zvertex60";

  TH2DBootstrap* h_MB_zvertex60_forcombine = (TH2DBootstrap*)f_MB->Get(histname.c_str());

  TAxis* xax = h_MB_zvertex60_forcombine->GetNominal()->GetXaxis();
  TAxis* yax = h_MB_zvertex60_forcombine->GetNominal()->GetYaxis();
  int nxbins = xax->GetNbins();
  int nybins = yax->GetNbins();

  const TArrayD* xbins = xax->GetXbins();
  const TArrayD* ybins = yax->GetXbins();
  const double* xedges;
  const double* yedges;
  double xmin, xmax, ymin, ymax;
  TH2DBootstrap* h_zvertex60_combined;

  if(xbins->GetSize() > 0)
    {
      xedges = xbins->GetArray();
      yedges = ybins->GetArray();
      h_zvertex60_combined = new TH2DBootstrap(histname.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nxbins,xedges,nybins,yedges,h_MB_zvertex60_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());
    }
  else
    {
      xmin = xax->GetXmin();
      xmax = xax->GetXmax();
      ymin = yax->GetXmin();
      ymax = yax->GetXmax();
      h_zvertex60_combined = new TH2DBootstrap(histname.c_str(),h_MB_zvertex60_forcombine->GetTitle(),nxbins,xmin,xmax,nybins,ymin,ymax,h_MB_zvertex60_forcombine->GetNReplica(),(BootstrapGenerator*)h_MB_zvertex60_forcombine->GetGenerator());


    }
  

  
  TH2DBootstrap* h_Jet5GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet5GeV->Get(histname.c_str());
  TH2DBootstrap* h_Jet12GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet12GeV->Get(histname.c_str());
  TH2DBootstrap* h_Jet40GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet40GeV->Get(histname.c_str());
  TH2DBootstrap* h_Jet20GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet20GeV->Get(histname.c_str());
  TH2DBootstrap* h_Jet30GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet30GeV->Get(histname.c_str());
  TH2DBootstrap* h_Jet50GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet50GeV->Get(histname.c_str());
  TH2DBootstrap* h_Jet60GeV_zvertex60_forcombine = (TH2DBootstrap*)f_Jet60GeV->Get(histname.c_str());
  h_zvertex60_combined->Add(h_MB_zvertex60_forcombine, mb_scale_all);
  h_zvertex60_combined->Add(h_Jet5GeV_zvertex60_forcombine, jet5_scale_all);
  h_zvertex60_combined->Add(h_Jet12GeV_zvertex60_forcombine, jet12_scale_all);
  h_zvertex60_combined->Add(h_Jet40GeV_zvertex60_forcombine, jet40_scale_all);
  h_zvertex60_combined->Add(h_Jet20GeV_zvertex60_forcombine, jet20_scale_all);
  h_zvertex60_combined->Add(h_Jet30GeV_zvertex60_forcombine, jet30_scale_all);
  h_zvertex60_combined->Add(h_Jet50GeV_zvertex60_forcombine, jet50_scale_all);
  h_zvertex60_combined->Add(h_Jet60GeV_zvertex60_forcombine, jet60_scale_all);

  f_out->cd();
  h_zvertex60_combined->Write();
}
