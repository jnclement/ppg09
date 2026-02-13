#include "unfold_Def.h"
#include "TRandom.h"
#include "TRandom3.h"
TRandom3 randGen(1234);

float min_dphi = 3*M_PI/4;
float jet_rad = 0.4;
void get_l_sl_jet(int& lindex, int& sindex, float* jet_pt, int jet_n)
{
  float lpt = 0;
  float spt = 0;
  lindex = -1;
  sindex = -1;

  for(int i=0; i<jet_n; ++i)
    {
      if(jet_pt[i] > lpt)
	{
	  sindex = lindex;
	  spt = lpt;
	  lindex = i;
	  lpt = jet_pt[i];
	}
      else if(jet_pt[i] > spt)
	{
	  sindex = i;
	  spt = jet_pt[i];
	}
    }
}

float get_dphi(float phi1, float phi2)
{
  float dphi = phi1-phi2;
  if(dphi > M_PI) dphi = 2*M_PI-dphi;
  if(dphi < -M_PI) dphi = 2*M_PI+dphi;
  return dphi;
}

float get_deta(float eta1, float eta2)
{
  return eta1-eta2;
}

float get_dR(float eta1, float phi1, float eta2, float phi2)
{
  float deta = get_deta(eta1, eta2);
  float dphi = get_dphi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

float get_dR_ellipsoid(float eta1, float phi1, float eta2, float phi2) {
  float deta = get_deta(eta1, eta2) / 3.;
  float dphi = get_dphi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

float get_emcal_mineta_zcorrected(float zvertex) {
  float minz_EM = -130.23;
  float radius_EM = 93.5;
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float maxz_EM = 130.23;
  float radius_EM = 93.5;
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float minz_IH = -170.299;
  float radius_IH = 127.503;
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float maxz_IH = 170.299;
  float radius_IH = 127.503;
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float minz_OH = -301.683;
  float radius_OH = 225.87;
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float maxz_OH = 301.683;
  float radius_OH = 225.87;
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zvertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zvertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zvertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zvertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zvertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zvertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zvertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

bool check_dphicut(float lphi, float sphi)
{
  return abs(get_dphi(lphi, sphi)) > min_dphi;
}

void filter_jets(vector<bool>& filter, float* jet_e, float* jet_pt, float* jet_eta, float zvertex, int jet_n)
{
  filter.clear();
  for(int i=0; i<jet_n; ++i)
    {
      filter.push_back(jet_e[i]<0 || check_bad_jet_eta(jet_eta[i], zvertex, jet_rad) || jet_pt[i]<1 || abs(jet_eta[i])>1.1-jet_rad);
    }
}

void filter_tjets(vector<bool>& filter, float* jet_e, float* jet_pt, float* jet_eta, float zvertex, int jet_n)
{
  filter.clear();
  for(int i=0; i<jet_n; ++i)
    {
      filter.push_back(jet_e[i]<0 || jet_pt[i]<1 || abs(jet_eta[i])>1.1-jet_rad);
    }
}
int sw = 0;
void match_meas_truth(std::vector<float>& meas_eta, std::vector<float>& meas_phi, std::vector<float>& meas_matched, std::vector<float>& truth_eta, std::vector<float>& truth_phi, std::vector<float>& truth_matched, float jet_radius, bool has_zvertex) {
  meas_matched.assign(meas_eta.size(), -1);
  truth_matched.assign(truth_eta.size(), -1);
  float max_match_dR = jet_radius * 0.75;
  for (int im = 0; im < meas_eta.size(); ++im) {
    float min_dR = 100;
    int match_index = -9999;
    for (int it = 0; it < truth_eta.size(); ++it) {
      float dR;
      if (has_zvertex) dR = get_dR(meas_eta[im], meas_phi[im], truth_eta[it], truth_phi[it]);
      else dR = get_dR_ellipsoid(meas_eta[im], meas_phi[im], truth_eta[it], truth_phi[it]);
      if (dR < min_dR) {
        match_index = it;
        min_dR = dR;
      }
    }
    if (min_dR < max_match_dR) {
      if (truth_matched[match_index] == -1) {
        meas_matched[im] = match_index;
        truth_matched[match_index] = im;
      } else {
        float dR1;
        if (has_zvertex) dR1 = get_dR(meas_eta[truth_matched[match_index]], meas_phi[truth_matched[match_index]], truth_eta[match_index], truth_phi[match_index]);
        else dR1 = get_dR_ellipsoid(meas_eta[truth_matched[match_index]], meas_phi[truth_matched[match_index]], truth_eta[match_index], truth_phi[match_index]);
        if (min_dR < dR1) {
          meas_matched[truth_matched[match_index]] = -1;
          meas_matched[im] = match_index;
          truth_matched[match_index] = im;
        }
      }
    }
  }
}

void fill_response_matrix(TH1D*& h_truth, TH1D*& h_meas, TH2D*& h_resp, TH1D*& h_fake, TH1D*& h_miss, TH1D*& h_matchtruth_rec, TH1D*& h_matchtruth_unw, TH1D*& h_meas_unw, double scale, TF1* f_reweight, std::vector<float>& meas_pt, std::vector<float>& meas_matched, std::vector<float>& truth_pt, std::vector<float>& truth_matched) {
  
  for (int im = 0; im < meas_pt.size(); ++im) {
    double reweight_factor = f_reweight->Eval(meas_pt[im]);
    
    if (meas_matched[im] < 0) {
      
      h_meas->Fill(meas_pt[im], scale*reweight_factor);
      h_fake->Fill(meas_pt[im], scale*reweight_factor);
      // For reweighting
      h_meas_unw->Fill(meas_pt[im], scale);
      
    } else {
      
      h_meas->Fill(meas_pt[im], scale*reweight_factor);
      h_resp->Fill(meas_pt[im], truth_pt[meas_matched[im]], scale*reweight_factor);
      // For reweighting
      h_meas_unw->Fill(meas_pt[im], scale);
      h_matchtruth_unw->Fill(truth_pt[meas_matched[im]], scale);
      h_matchtruth_rec->Fill(truth_pt[meas_matched[im]], scale*reweight_factor);
      
    }
  }
  
  for (int it = 0; it < truth_pt.size(); ++it) {
    
    if (truth_matched[it] < 0) {
      h_truth->Fill(truth_pt[it], scale);
      h_miss->Fill(truth_pt[it], scale);
    } else {
      h_truth->Fill(truth_pt[it], scale);
    }
  }
}

void build_hist(std::string tag, TH1D*& h_truth, TH1D*& h_measure, TH2D*& h_respmatrix, TH1D*& h_fake, TH1D*& h_miss, TH1D*& h_matchedtruth_weighted, TH1D*& h_matchedtruth_unweighted, TH1D*& h_measure_unweighted) {
  h_truth = new TH1D(Form("h_truth_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  h_measure = new TH1D(Form("h_measure_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  h_respmatrix = new TH2D(Form("h_respmatrix_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  h_fake = new TH1D(Form("h_fake_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  h_miss = new TH1D(Form("h_miss_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  h_matchedtruth_weighted = new TH1D(Form("h_matchedtruth_weighted_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", 1000, 0, 100);
  h_matchedtruth_unweighted = new TH1D(Form("h_matchedtruth_unweighted_%s", tag.c_str()), ";p_{T}^{Truth jet} [GeV]", 1000, 0, 100);
  h_measure_unweighted = new TH1D(Form("h_measure_unweighted_%s", tag.c_str()), ";p_{T}^{Calib jet} [GeV]", 1000, 0, 100);
}


bool isbg;
void get_calibjet(std::vector<float>& calibjet_pt, std::vector<float>& calibjet_eta, std::vector<float>& calibjet_phi, std::vector<bool>& jet_filter, int jet_n, float* jet_pt, float* jet_eta, float* jet_phi, float jes, float jer)
{
  calibjet_pt.clear();
  calibjet_eta.clear();
  calibjet_phi.clear();

  for(int i=0; i<jet_n; ++i)
    {
      if(jet_filter.at(i) || isbg) continue;
      double calib_pt;
      if(jer != 0) calib_pt = jet_pt[i]*(1+randGen.Gaus(0.0,jer))*jes;
      else calib_pt = jet_pt[i];
      if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) continue;
      calibjet_pt.push_back(calib_pt);
      calibjet_eta.push_back(jet_eta[i]);
      calibjet_phi.push_back(jet_phi[i]);
    }
}

void get_truthjet(std::vector<float>& goodtruthjet_pt, std::vector<float>& goodtruthjet_eta, std::vector<float>& goodtruthjet_phi, std::vector<bool>& truthjet_filter, int tjet_n, float* tjet_pt, float* tjet_eta, float* tjet_phi)
{
  goodtruthjet_pt.clear();
  goodtruthjet_eta.clear();
  goodtruthjet_phi.clear();
  for(int i=0; i<tjet_n; ++i)
    {
      if(truthjet_filter.at(i)) continue;
      if(tjet_pt[i] < truthptbins[0] || tjet_pt[i] > truthptbins[truthnpt]) continue;
      goodtruthjet_pt.push_back(tjet_pt[i]);
      goodtruthjet_eta.push_back(tjet_eta[i]);
      goodtruthjet_phi.push_back(tjet_phi[i]);
    }
}
  
int analyze_segment_sim(string runtype, int iseg, int nseg)
{

  double truthjet_pt_min = 0, truthjet_pt_max = 1000, recojet_pt_max = 1000;
  if (runtype == "mb") {
    if (jet_rad == 0.2) {truthjet_pt_min = 0; truthjet_pt_max = 5; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 0; truthjet_pt_max = 6; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 0; truthjet_pt_max = 7; recojet_pt_max = 14;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 0; truthjet_pt_max = 10; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 0; truthjet_pt_max = 11; recojet_pt_max = 1000;}
    truthjet_pt_min = 0; truthjet_pt_max = 7; recojet_pt_max = 14;
  } else if (runtype == "jet5") {
    if (jet_rad == 0.2) {truthjet_pt_min = 5; truthjet_pt_max = 12; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 6; truthjet_pt_max = 13; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 7; truthjet_pt_max = 14; recojet_pt_max = 24;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 10; truthjet_pt_max = 15; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 11; truthjet_pt_max = 17; recojet_pt_max = 1000;}
    truthjet_pt_min = 7; truthjet_pt_max = 14; recojet_pt_max = 24;
  } else if (runtype == "jet10") {
    if (jet_rad == 0.2) {truthjet_pt_min = 12; truthjet_pt_max = 15; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 13; truthjet_pt_max = 16; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 14; truthjet_pt_max = 17; recojet_pt_max = 33;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 15; truthjet_pt_max = 24; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 17; truthjet_pt_max = 26; recojet_pt_max = 1000;}
    truthjet_pt_min = 14; truthjet_pt_max = 17; recojet_pt_max = 33;
  } else if (runtype == "jet15") {
    if (jet_rad == 0.2) {truthjet_pt_min = 15; truthjet_pt_max = 20; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 16; truthjet_pt_max = 21; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 17; truthjet_pt_max = 22; recojet_pt_max = 45;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 24; truthjet_pt_max = 30; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 26; truthjet_pt_max = 35; recojet_pt_max = 1000;}
    truthjet_pt_min = 17; truthjet_pt_max = 22; recojet_pt_max = 45;
  } else if (runtype == "jet20") {
    if (jet_rad == 0.2) {truthjet_pt_min = 20; truthjet_pt_max = 31; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 21; truthjet_pt_max = 33; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 22; truthjet_pt_max = 35; recojet_pt_max = 59;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 30; truthjet_pt_max = 40; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 35; truthjet_pt_max = 45; recojet_pt_max = 1000;}
    truthjet_pt_min = 22; truthjet_pt_max = 35; recojet_pt_max = 59;
  } else if (runtype == "jet30") {
    if (jet_rad == 0.2) {truthjet_pt_min = 31; truthjet_pt_max = 50; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 33; truthjet_pt_max = 51; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 35; truthjet_pt_max = 52; recojet_pt_max = 72;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 40; truthjet_pt_max = 60; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 45; truthjet_pt_max = 63; recojet_pt_max = 1000;}
    truthjet_pt_min = 35; truthjet_pt_max = 52; recojet_pt_max = 72;
  } else if (runtype == "jet50") {
    if (jet_rad == 0.2) {truthjet_pt_min = 50; truthjet_pt_max = 70; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 51; truthjet_pt_max = 70; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 52; truthjet_pt_max = 71; recojet_pt_max = 1000;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 60; truthjet_pt_max = 75; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 63; truthjet_pt_max = 79; recojet_pt_max = 1000;}
    truthjet_pt_min = 52; truthjet_pt_max = 71; recojet_pt_max = 1000;
  } else if (runtype == "jet70") {
    if (jet_rad == 0.2) {truthjet_pt_min = 70; truthjet_pt_max = 3000; recojet_pt_max = 1000;}
    else if (jet_rad == 0.3) {truthjet_pt_min = 70; truthjet_pt_max = 3000; recojet_pt_max = 1000;}
    else if (jet_rad == 0.4) {truthjet_pt_min = 71; truthjet_pt_max = 3000; recojet_pt_max = 1000;}
    else if (jet_rad == 0.5) {truthjet_pt_min = 75; truthjet_pt_max = 3000; recojet_pt_max = 1000;}
    else if (jet_rad == 0.6) {truthjet_pt_min = 79; truthjet_pt_max = 3000; recojet_pt_max = 1000;}
    truthjet_pt_min = 71; truthjet_pt_max = 3000; recojet_pt_max = 1000;
  }

  
  TChain chain("jet_tree");
  
  for(int i=0; i<nseg; ++i)
    {
      chain.Add(("input/inputfile_"+to_string(i)+".root").c_str());
    }

  float zvtx;
  float jet_e[100];
  float jet_pt[100];
  float jet_et[100];
  float jet_eta[100];
  float jet_phi[100];
  float jet_t[100];
  float jet_pt_calib[100];
  float tjet_e[100];
  float tjet_pt[100];
  float tjet_et[100];
  float tjet_eta[100];
  float tjet_phi[100];
  float tjet_t[100];
  
  int jet_n, tjet_n;
  int calib_jet_n;

  chain.SetBranchAddress("zvtx",&zvtx);
  chain.SetBranchAddress("jet_et",jet_e);
  chain.SetBranchAddress("jet_pt",jet_pt);
  chain.SetBranchAddress("jet_etrans",jet_et);
  chain.SetBranchAddress("jet_eta",jet_eta);
  chain.SetBranchAddress("jet_phi",jet_phi);
  chain.SetBranchAddress("jet_t",jet_t);
  chain.SetBranchAddress("jet_pt_calib",jet_pt_calib);
  chain.SetBranchAddress("jet_n",&jet_n);
  chain.SetBranchAddress("calib_jet_n",&calib_jet_n);

  chain.SetBranchAddress("tjet_e",tjet_e);
  chain.SetBranchAddress("tjet_pt",tjet_pt);
  //chain.SetBranchAddress("tjet_et",tjet_et);
  chain.SetBranchAddress("tjet_eta",tjet_eta);
  chain.SetBranchAddress("tjet_phi",tjet_phi);
  //chain.SetBranchAddress("tjet_t",tjet_t);
  chain.SetBranchAddress("tjet_n",&tjet_n);


  TFile *f_zvertexreweight = new TFile("input_zvertexreweight.root", "READ");
  if (!f_zvertexreweight || f_zvertexreweight->IsZombie()) {
    std::cout << "Error: cannot open input_zvertexreweight.root" << std::endl;
    return -1;
  }
  TH1D *h_zvertex_reweight = (TH1D*)f_zvertexreweight->Get("h_zvertex_reweight");
  if (!h_zvertex_reweight) {
    std::cout << "Error: cannot open h_zvertex_reweight" << std::endl;
    return -1;
  }
  
  TFile* fmb = TFile::Open("output_mbdefficiency.root","READ");
  TF1* mbt[3];
  mbt[0] = (TF1*)fmb->Get("mbdtrig04_nominal");
  mbt[1] = (TF1*)fmb->Get("mbdtrig04_up");
  mbt[2] = (TF1*)fmb->Get("mbdtrig04_down");
  int jet_radius_index = 4;
  TFile *f_reweight = new TFile(Form("output_reweightfunction_r0%d.root", jet_radius_index), "READ");
  TF1 *f_reweightfunc_all = (TF1*)f_reweight->Get("reweightfunc_all");
  TF1 *f_reweightfunc_all_jesup = (TF1*)f_reweight->Get("reweightfunc_all_jesup");
  TF1 *f_reweightfunc_all_jesdown = (TF1*)f_reweight->Get("reweightfunc_all_jesdown");
  TF1 *f_reweightfunc_all_jerup = (TF1*)f_reweight->Get("reweightfunc_all_jerup");
  TF1 *f_reweightfunc_all_jerdown = (TF1*)f_reweight->Get("reweightfunc_all_jerdown");
  TF1 *f_reweightfunc_all_jetup = (TF1*)f_reweight->Get("reweightfunc_all_jetup");
  TF1 *f_reweightfunc_all_jetdown = (TF1*)f_reweight->Get("reweightfunc_all_jetdown");
  TF1 *f_reweightfunc_zvertex30 = (TF1*)f_reweight->Get("reweightfunc_zvertex30");
  TF1 *f_reweightfunc_zvertex30_jesup = (TF1*)f_reweight->Get("reweightfunc_zvertex30_jesup");
  TF1 *f_reweightfunc_zvertex30_jesdown = (TF1*)f_reweight->Get("reweightfunc_zvertex30_jesdown");
  TF1 *f_reweightfunc_zvertex30_jerup = (TF1*)f_reweight->Get("reweightfunc_zvertex30_jerup");
  TF1 *f_reweightfunc_zvertex30_jerdown = (TF1*)f_reweight->Get("reweightfunc_zvertex30_jerdown");
  TF1 *f_reweightfunc_zvertex30_jetup = (TF1*)f_reweight->Get("reweightfunc_zvertex30_jetup");
  TF1 *f_reweightfunc_zvertex30_jetdown = (TF1*)f_reweight->Get("reweightfunc_zvertex30_jetdown");
  TF1 *f_reweightfunc_zvertex30_mbdup = (TF1*)f_reweight->Get("reweightfunc_zvertex30_mbdup");
  TF1 *f_reweightfunc_zvertex30_mbddown = (TF1*)f_reweight->Get("reweightfunc_zvertex30_mbddown");
  TF1 *f_reweightfunc_zvertex60 = (TF1*)f_reweight->Get("reweightfunc_zvertex60");
  TF1 *f_reweightfunc_zvertex60_jesup = (TF1*)f_reweight->Get("reweightfunc_zvertex60_jesup");
  TF1 *f_reweightfunc_zvertex60_jesdown = (TF1*)f_reweight->Get("reweightfunc_zvertex60_jesdown");
  TF1 *f_reweightfunc_zvertex60_jerup = (TF1*)f_reweight->Get("reweightfunc_zvertex60_jerup");
  TF1 *f_reweightfunc_zvertex60_jerdown = (TF1*)f_reweight->Get("reweightfunc_zvertex60_jerdown");
  TF1 *f_reweightfunc_zvertex60_jetup = (TF1*)f_reweight->Get("reweightfunc_zvertex60_jetup");
  TF1 *f_reweightfunc_zvertex60_jetdown = (TF1*)f_reweight->Get("reweightfunc_zvertex60_jetdown");
  TF1 *f_reweightfunc_zvertex60_mbdup = (TF1*)f_reweight->Get("reweightfunc_zvertex60_mbdup");
  TF1 *f_reweightfunc_zvertex60_mbddown = (TF1*)f_reweight->Get("reweightfunc_zvertex60_mbddown");

  ////////// Histograms //////////
  TH1D *h_event_all = new TH1D("h_event_all", ";Event Number", 1, 0, 1);
  TH1D *h_event_beforecut = new TH1D("h_event_beforecut", ";Event Number", 1, 0, 1);
  TH1D *h_event_passed = new TH1D("h_event_passed", ";Event Number", 1, 0, 1);
  TH1D *h_recojet_pt_record_nocut_all = new TH1D("h_recojet_pt_record_nocut_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_all = new TH1D("h_recojet_pt_record_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_nocut_zvertex30 = new TH1D("h_recojet_pt_record_nocut_zvertex30", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_zvertex30 = new TH1D("h_recojet_pt_record_zvertex30", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_nocut_zvertex60 = new TH1D("h_recojet_pt_record_nocut_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_zvertex60 = new TH1D("h_recojet_pt_record_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_truthjet_pt_record_all = new TH1D("h_truthjet_pt_record_all", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  // Nominal histograms
  TH2D *h_respmatrix_all; TH1D *h_truth_all, *h_measure_all, *h_fake_all, *h_miss_all, *h_matchedtruth_weighted_all, *h_matchedtruth_unweighted_all, *h_measure_unweighted_all;
  build_hist("all", h_truth_all, h_measure_all, h_respmatrix_all, h_fake_all, h_miss_all, h_matchedtruth_weighted_all, h_matchedtruth_unweighted_all, h_measure_unweighted_all);

  TH2D *h_respmatrix_all_nosmear; TH1D *h_truth_all_nosmear, *h_measure_all_nosmear, *h_fake_all_nosmear, *h_miss_all_nosmear, *h_matchedtruth_weighted_all_nosmear, *h_matchedtruth_unweighted_all_nosmear, *h_measure_unweighted_all_nosmear;
  build_hist("all_nosmear", h_truth_all_nosmear, h_measure_all_nosmear, h_respmatrix_all_nosmear, h_fake_all_nosmear, h_miss_all_nosmear, h_matchedtruth_weighted_all_nosmear, h_matchedtruth_unweighted_all_nosmear, h_measure_unweighted_all_nosmear);

  TH2D *h_respmatrix_zvertex60_nosmear; TH1D *h_truth_zvertex60_nosmear, *h_measure_zvertex60_nosmear, *h_fake_zvertex60_nosmear, *h_miss_zvertex60_nosmear, *h_matchedtruth_weighted_zvertex60_nosmear, *h_matchedtruth_unweighted_zvertex60_nosmear, *h_measure_unweighted_zvertex60_nosmear;
  build_hist("zvertex60_nosmear", h_truth_zvertex60_nosmear, h_measure_zvertex60_nosmear, h_respmatrix_zvertex60_nosmear, h_fake_zvertex60_nosmear, h_miss_zvertex60_nosmear, h_matchedtruth_weighted_zvertex60_nosmear, h_matchedtruth_unweighted_zvertex60_nosmear, h_measure_unweighted_zvertex60_nosmear);

  
  TH2D *h_respmatrix_zvertex30; TH1D *h_truth_zvertex30, *h_measure_zvertex30, *h_fake_zvertex30, *h_miss_zvertex30, *h_matchedtruth_weighted_zvertex30, *h_matchedtruth_unweighted_zvertex30, *h_measure_unweighted_zvertex30;
  build_hist("zvertex30", h_truth_zvertex30, h_measure_zvertex30, h_respmatrix_zvertex30, h_fake_zvertex30, h_miss_zvertex30, h_matchedtruth_weighted_zvertex30, h_matchedtruth_unweighted_zvertex30, h_measure_unweighted_zvertex30);
  TH2D *h_respmatrix_zvertex60; TH1D *h_truth_zvertex60, *h_measure_zvertex60, *h_fake_zvertex60, *h_miss_zvertex60, *h_matchedtruth_weighted_zvertex60, *h_matchedtruth_unweighted_zvertex60, *h_measure_unweighted_zvertex60;
  build_hist("zvertex60", h_truth_zvertex60, h_measure_zvertex60, h_respmatrix_zvertex60, h_fake_zvertex60, h_miss_zvertex60, h_matchedtruth_weighted_zvertex60, h_matchedtruth_unweighted_zvertex60, h_measure_unweighted_zvertex60);

  // JES up/down histograms
  TH2D *h_respmatrix_all_jesup; TH1D *h_truth_all_jesup, *h_measure_all_jesup, *h_fake_all_jesup, *h_miss_all_jesup, *h_matchedtruth_weighted_all_jesup, *h_matchedtruth_unweighted_all_jesup, *h_measure_unweighted_all_jesup;
  build_hist("all_jesup", h_truth_all_jesup, h_measure_all_jesup, h_respmatrix_all_jesup, h_fake_all_jesup, h_miss_all_jesup, h_matchedtruth_weighted_all_jesup, h_matchedtruth_unweighted_all_jesup, h_measure_unweighted_all_jesup);
  TH2D *h_respmatrix_all_jesdown; TH1D *h_truth_all_jesdown, *h_measure_all_jesdown, *h_fake_all_jesdown, *h_miss_all_jesdown, *h_matchedtruth_weighted_all_jesdown, *h_matchedtruth_unweighted_all_jesdown, *h_measure_unweighted_all_jesdown;
  build_hist("all_jesdown", h_truth_all_jesdown, h_measure_all_jesdown, h_respmatrix_all_jesdown, h_fake_all_jesdown, h_miss_all_jesdown, h_matchedtruth_weighted_all_jesdown, h_matchedtruth_unweighted_all_jesdown, h_measure_unweighted_all_jesdown);
  TH2D *h_respmatrix_zvertex30_jesup; TH1D *h_truth_zvertex30_jesup, *h_measure_zvertex30_jesup, *h_fake_zvertex30_jesup, *h_miss_zvertex30_jesup, *h_matchedtruth_weighted_zvertex30_jesup, *h_matchedtruth_unweighted_zvertex30_jesup, *h_measure_unweighted_zvertex30_jesup;
  build_hist("zvertex30_jesup", h_truth_zvertex30_jesup, h_measure_zvertex30_jesup, h_respmatrix_zvertex30_jesup, h_fake_zvertex30_jesup, h_miss_zvertex30_jesup, h_matchedtruth_weighted_zvertex30_jesup, h_matchedtruth_unweighted_zvertex30_jesup, h_measure_unweighted_zvertex30_jesup);
  TH2D *h_respmatrix_zvertex30_jesdown; TH1D *h_truth_zvertex30_jesdown, *h_measure_zvertex30_jesdown, *h_fake_zvertex30_jesdown, *h_miss_zvertex30_jesdown, *h_matchedtruth_weighted_zvertex30_jesdown, *h_matchedtruth_unweighted_zvertex30_jesdown, *h_measure_unweighted_zvertex30_jesdown;
  build_hist("zvertex30_jesdown", h_truth_zvertex30_jesdown, h_measure_zvertex30_jesdown, h_respmatrix_zvertex30_jesdown, h_fake_zvertex30_jesdown, h_miss_zvertex30_jesdown, h_matchedtruth_weighted_zvertex30_jesdown, h_matchedtruth_unweighted_zvertex30_jesdown, h_measure_unweighted_zvertex30_jesdown);
  TH2D *h_respmatrix_zvertex60_jesup; TH1D *h_truth_zvertex60_jesup, *h_measure_zvertex60_jesup, *h_fake_zvertex60_jesup, *h_miss_zvertex60_jesup, *h_matchedtruth_weighted_zvertex60_jesup, *h_matchedtruth_unweighted_zvertex60_jesup, *h_measure_unweighted_zvertex60_jesup;
  build_hist("zvertex60_jesup", h_truth_zvertex60_jesup, h_measure_zvertex60_jesup, h_respmatrix_zvertex60_jesup, h_fake_zvertex60_jesup, h_miss_zvertex60_jesup, h_matchedtruth_weighted_zvertex60_jesup, h_matchedtruth_unweighted_zvertex60_jesup, h_measure_unweighted_zvertex60_jesup);
  TH2D *h_respmatrix_zvertex60_jesdown; TH1D *h_truth_zvertex60_jesdown, *h_measure_zvertex60_jesdown, *h_fake_zvertex60_jesdown, *h_miss_zvertex60_jesdown, *h_matchedtruth_weighted_zvertex60_jesdown, *h_matchedtruth_unweighted_zvertex60_jesdown, *h_measure_unweighted_zvertex60_jesdown;
  build_hist("zvertex60_jesdown", h_truth_zvertex60_jesdown, h_measure_zvertex60_jesdown, h_respmatrix_zvertex60_jesdown, h_fake_zvertex60_jesdown, h_miss_zvertex60_jesdown, h_matchedtruth_weighted_zvertex60_jesdown, h_matchedtruth_unweighted_zvertex60_jesdown, h_measure_unweighted_zvertex60_jesdown);

  // JER up/down histograms
  TH2D *h_respmatrix_all_jerup; TH1D *h_truth_all_jerup, *h_measure_all_jerup, *h_fake_all_jerup, *h_miss_all_jerup, *h_matchedtruth_weighted_all_jerup, *h_matchedtruth_unweighted_all_jerup, *h_measure_unweighted_all_jerup;
  build_hist("all_jerup", h_truth_all_jerup, h_measure_all_jerup, h_respmatrix_all_jerup, h_fake_all_jerup, h_miss_all_jerup, h_matchedtruth_weighted_all_jerup, h_matchedtruth_unweighted_all_jerup, h_measure_unweighted_all_jerup);
  TH2D *h_respmatrix_all_jerdown; TH1D *h_truth_all_jerdown, *h_measure_all_jerdown, *h_fake_all_jerdown, *h_miss_all_jerdown, *h_matchedtruth_weighted_all_jerdown, *h_matchedtruth_unweighted_all_jerdown, *h_measure_unweighted_all_jerdown;
  build_hist("all_jerdown", h_truth_all_jerdown, h_measure_all_jerdown, h_respmatrix_all_jerdown, h_fake_all_jerdown, h_miss_all_jerdown, h_matchedtruth_weighted_all_jerdown, h_matchedtruth_unweighted_all_jerdown, h_measure_unweighted_all_jerdown);
  TH2D *h_respmatrix_zvertex30_jerup; TH1D *h_truth_zvertex30_jerup, *h_measure_zvertex30_jerup, *h_fake_zvertex30_jerup, *h_miss_zvertex30_jerup, *h_matchedtruth_weighted_zvertex30_jerup, *h_matchedtruth_unweighted_zvertex30_jerup, *h_measure_unweighted_zvertex30_jerup;
  build_hist("zvertex30_jerup", h_truth_zvertex30_jerup, h_measure_zvertex30_jerup, h_respmatrix_zvertex30_jerup, h_fake_zvertex30_jerup, h_miss_zvertex30_jerup, h_matchedtruth_weighted_zvertex30_jerup, h_matchedtruth_unweighted_zvertex30_jerup, h_measure_unweighted_zvertex30_jerup);
  TH2D *h_respmatrix_zvertex30_jerdown; TH1D *h_truth_zvertex30_jerdown, *h_measure_zvertex30_jerdown, *h_fake_zvertex30_jerdown, *h_miss_zvertex30_jerdown, *h_matchedtruth_weighted_zvertex30_jerdown, *h_matchedtruth_unweighted_zvertex30_jerdown, *h_measure_unweighted_zvertex30_jerdown;
  build_hist("zvertex30_jerdown", h_truth_zvertex30_jerdown, h_measure_zvertex30_jerdown, h_respmatrix_zvertex30_jerdown, h_fake_zvertex30_jerdown, h_miss_zvertex30_jerdown, h_matchedtruth_weighted_zvertex30_jerdown, h_matchedtruth_unweighted_zvertex30_jerdown, h_measure_unweighted_zvertex30_jerdown);
  TH2D *h_respmatrix_zvertex60_jerup; TH1D *h_truth_zvertex60_jerup, *h_measure_zvertex60_jerup, *h_fake_zvertex60_jerup, *h_miss_zvertex60_jerup, *h_matchedtruth_weighted_zvertex60_jerup, *h_matchedtruth_unweighted_zvertex60_jerup, *h_measure_unweighted_zvertex60_jerup;
  build_hist("zvertex60_jerup", h_truth_zvertex60_jerup, h_measure_zvertex60_jerup, h_respmatrix_zvertex60_jerup, h_fake_zvertex60_jerup, h_miss_zvertex60_jerup, h_matchedtruth_weighted_zvertex60_jerup, h_matchedtruth_unweighted_zvertex60_jerup, h_measure_unweighted_zvertex60_jerup);
  TH2D *h_respmatrix_zvertex60_jerdown; TH1D *h_truth_zvertex60_jerdown, *h_measure_zvertex60_jerdown, *h_fake_zvertex60_jerdown, *h_miss_zvertex60_jerdown, *h_matchedtruth_weighted_zvertex60_jerdown, *h_matchedtruth_unweighted_zvertex60_jerdown, *h_measure_unweighted_zvertex60_jerdown;
  build_hist("zvertex60_jerdown", h_truth_zvertex60_jerdown, h_measure_zvertex60_jerdown, h_respmatrix_zvertex60_jerdown, h_fake_zvertex60_jerdown, h_miss_zvertex60_jerdown, h_matchedtruth_weighted_zvertex60_jerdown, h_matchedtruth_unweighted_zvertex60_jerdown, h_measure_unweighted_zvertex60_jerdown);

  // Jet trigger up/down histograms
  TH2D *h_respmatrix_all_jetup; TH1D *h_truth_all_jetup, *h_measure_all_jetup, *h_fake_all_jetup, *h_miss_all_jetup, *h_matchedtruth_weighted_all_jetup, *h_matchedtruth_unweighted_all_jetup, *h_measure_unweighted_all_jetup;
  build_hist("all_jetup", h_truth_all_jetup, h_measure_all_jetup, h_respmatrix_all_jetup, h_fake_all_jetup, h_miss_all_jetup, h_matchedtruth_weighted_all_jetup, h_matchedtruth_unweighted_all_jetup, h_measure_unweighted_all_jetup);
  TH2D *h_respmatrix_all_jetdown; TH1D *h_truth_all_jetdown, *h_measure_all_jetdown, *h_fake_all_jetdown, *h_miss_all_jetdown, *h_matchedtruth_weighted_all_jetdown, *h_matchedtruth_unweighted_all_jetdown, *h_measure_unweighted_all_jetdown;
  build_hist("all_jetdown", h_truth_all_jetdown, h_measure_all_jetdown, h_respmatrix_all_jetdown, h_fake_all_jetdown, h_miss_all_jetdown, h_matchedtruth_weighted_all_jetdown, h_matchedtruth_unweighted_all_jetdown, h_measure_unweighted_all_jetdown);
  TH2D *h_respmatrix_zvertex30_jetup; TH1D *h_truth_zvertex30_jetup, *h_measure_zvertex30_jetup, *h_fake_zvertex30_jetup, *h_miss_zvertex30_jetup, *h_matchedtruth_weighted_zvertex30_jetup, *h_matchedtruth_unweighted_zvertex30_jetup, *h_measure_unweighted_zvertex30_jetup;
  build_hist("zvertex30_jetup", h_truth_zvertex30_jetup, h_measure_zvertex30_jetup, h_respmatrix_zvertex30_jetup, h_fake_zvertex30_jetup, h_miss_zvertex30_jetup, h_matchedtruth_weighted_zvertex30_jetup, h_matchedtruth_unweighted_zvertex30_jetup, h_measure_unweighted_zvertex30_jetup);
  TH2D *h_respmatrix_zvertex30_jetdown; TH1D *h_truth_zvertex30_jetdown, *h_measure_zvertex30_jetdown, *h_fake_zvertex30_jetdown, *h_miss_zvertex30_jetdown, *h_matchedtruth_weighted_zvertex30_jetdown, *h_matchedtruth_unweighted_zvertex30_jetdown, *h_measure_unweighted_zvertex30_jetdown;
  build_hist("zvertex30_jetdown", h_truth_zvertex30_jetdown, h_measure_zvertex30_jetdown, h_respmatrix_zvertex30_jetdown, h_fake_zvertex30_jetdown, h_miss_zvertex30_jetdown, h_matchedtruth_weighted_zvertex30_jetdown, h_matchedtruth_unweighted_zvertex30_jetdown, h_measure_unweighted_zvertex30_jetdown);
  TH2D *h_respmatrix_zvertex60_jetup; TH1D *h_truth_zvertex60_jetup, *h_measure_zvertex60_jetup, *h_fake_zvertex60_jetup, *h_miss_zvertex60_jetup, *h_matchedtruth_weighted_zvertex60_jetup, *h_matchedtruth_unweighted_zvertex60_jetup, *h_measure_unweighted_zvertex60_jetup;
  build_hist("zvertex60_jetup", h_truth_zvertex60_jetup, h_measure_zvertex60_jetup, h_respmatrix_zvertex60_jetup, h_fake_zvertex60_jetup, h_miss_zvertex60_jetup, h_matchedtruth_weighted_zvertex60_jetup, h_matchedtruth_unweighted_zvertex60_jetup, h_measure_unweighted_zvertex60_jetup);
  TH2D *h_respmatrix_zvertex60_jetdown; TH1D *h_truth_zvertex60_jetdown, *h_measure_zvertex60_jetdown, *h_fake_zvertex60_jetdown, *h_miss_zvertex60_jetdown, *h_matchedtruth_weighted_zvertex60_jetdown, *h_matchedtruth_unweighted_zvertex60_jetdown, *h_measure_unweighted_zvertex60_jetdown;
  build_hist("zvertex60_jetdown", h_truth_zvertex60_jetdown, h_measure_zvertex60_jetdown, h_respmatrix_zvertex60_jetdown, h_fake_zvertex60_jetdown, h_miss_zvertex60_jetdown, h_matchedtruth_weighted_zvertex60_jetdown, h_matchedtruth_unweighted_zvertex60_jetdown, h_measure_unweighted_zvertex60_jetdown);

  // MBD trigger up/down histograms
  TH2D *h_respmatrix_zvertex30_mbdup; TH1D *h_truth_zvertex30_mbdup, *h_measure_zvertex30_mbdup, *h_fake_zvertex30_mbdup, *h_miss_zvertex30_mbdup, *h_matchedtruth_weighted_zvertex30_mbdup, *h_matchedtruth_unweighted_zvertex30_mbdup, *h_measure_unweighted_zvertex30_mbdup;
  build_hist("zvertex30_mbdup", h_truth_zvertex30_mbdup, h_measure_zvertex30_mbdup, h_respmatrix_zvertex30_mbdup, h_fake_zvertex30_mbdup, h_miss_zvertex30_mbdup, h_matchedtruth_weighted_zvertex30_mbdup, h_matchedtruth_unweighted_zvertex30_mbdup, h_measure_unweighted_zvertex30_mbdup);
  TH2D *h_respmatrix_zvertex30_mbddown; TH1D *h_truth_zvertex30_mbddown, *h_measure_zvertex30_mbddown, *h_fake_zvertex30_mbddown, *h_miss_zvertex30_mbddown, *h_matchedtruth_weighted_zvertex30_mbddown, *h_matchedtruth_unweighted_zvertex30_mbddown, *h_measure_unweighted_zvertex30_mbddown;
  build_hist("zvertex30_mbddown", h_truth_zvertex30_mbddown, h_measure_zvertex30_mbddown, h_respmatrix_zvertex30_mbddown, h_fake_zvertex30_mbddown, h_miss_zvertex30_mbddown, h_matchedtruth_weighted_zvertex30_mbddown, h_matchedtruth_unweighted_zvertex30_mbddown, h_measure_unweighted_zvertex30_mbddown);
  TH2D *h_respmatrix_zvertex60_mbdup; TH1D *h_truth_zvertex60_mbdup, *h_measure_zvertex60_mbdup, *h_fake_zvertex60_mbdup, *h_miss_zvertex60_mbdup, *h_matchedtruth_weighted_zvertex60_mbdup, *h_matchedtruth_unweighted_zvertex60_mbdup, *h_measure_unweighted_zvertex60_mbdup;
  build_hist("zvertex60_mbdup", h_truth_zvertex60_mbdup, h_measure_zvertex60_mbdup, h_respmatrix_zvertex60_mbdup, h_fake_zvertex60_mbdup, h_miss_zvertex60_mbdup, h_matchedtruth_weighted_zvertex60_mbdup, h_matchedtruth_unweighted_zvertex60_mbdup, h_measure_unweighted_zvertex60_mbdup);
  TH2D *h_respmatrix_zvertex60_mbddown; TH1D *h_truth_zvertex60_mbddown, *h_measure_zvertex60_mbddown, *h_fake_zvertex60_mbddown, *h_miss_zvertex60_mbddown, *h_matchedtruth_weighted_zvertex60_mbddown, *h_matchedtruth_unweighted_zvertex60_mbddown, *h_measure_unweighted_zvertex60_mbddown;
  build_hist("zvertex60_mbddown", h_truth_zvertex60_mbddown, h_measure_zvertex60_mbddown, h_respmatrix_zvertex60_mbddown, h_fake_zvertex60_mbddown, h_miss_zvertex60_mbddown, h_matchedtruth_weighted_zvertex60_mbddown, h_matchedtruth_unweighted_zvertex60_mbddown, h_measure_unweighted_zvertex60_mbddown);

  // Closure test histograms
  TH1D *h_fullclosure_measure_zvertex60 = new TH1D("h_fullclosure_measure_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_fullclosure_truth_zvertex60 = new TH1D("h_fullclosure_truth_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH2D *h_fullclosure_respmatrix_zvertex60 = new TH2D("h_fullclosure_respmatrix_zvertex60", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH1D *h_fullclosure_fake_zvertex60 = new TH1D("h_fullclosure_fake_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_fullclosure_miss_zvertex60 = new TH1D("h_fullclosure_miss_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);

  TH1D *h_halfclosure_inputmeasure_zvertex60 = new TH1D("h_halfclosure_inputmeasure_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_halfclosure_measure_zvertex60 = new TH1D("h_halfclosure_measure_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_halfclosure_truth_zvertex60 = new TH1D("h_halfclosure_truth_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH2D *h_halfclosure_respmatrix_zvertex60 = new TH2D("h_halfclosure_respmatrix_zvertex60", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH1D *h_halfclosure_fake_zvertex60 = new TH1D("h_halfclosure_fake_zvertex60", ";p_{T}^{Calib jet} [GeV]", calibnpt, calibptbins);
  TH1D *h_halfclosure_miss_zvertex60 = new TH1D("h_halfclosure_miss_zvertex60", ";p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);


  
  
  long long unsigned int nevt = chain.GetEntries();
  bool bgdj, bgnj, z30, z60, is22, is18;
  vector<bool> jet_filter, tjet_filter;
  std::vector<float> goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched;
  std::vector<float> calibjet_pt, calibjet_eta, calibjet_phi, calibjet_matched;

  std::vector<float> calibjet_pt_nosmear, calibjet_eta_nosmear, calibjet_phi_nosmear, calibjet_matched_nosmear;
    
  std::vector<float> calibjet_pt_jesup, calibjet_eta_jesup, calibjet_phi_jesup, calibjet_matched_jesup;
  std::vector<float> calibjet_pt_jesdown, calibjet_eta_jesdown, calibjet_phi_jesdown, calibjet_matched_jesdown;
  std::vector<float> calibjet_pt_jerup, calibjet_eta_jerup, calibjet_phi_jerup, calibjet_matched_jerup;
  std::vector<float> calibjet_pt_jerdown, calibjet_eta_jerdown, calibjet_phi_jerdown, calibjet_matched_jerdown;

  cout << nevt << endl;
  
  for(long long unsigned int i=0; i<nevt; ++i)
    {
            
      chain.GetEntry(i);
      z30 = abs(zvtx) < 30 && zvtx!=0;
      z60 = abs(zvtx) < 60 && zvtx!=0;
      bool has_zvtx = true;
      if(abs(zvtx)>990 || zvtx==0)
	{
	  zvtx = 0;
	  has_zvtx = false;
	}
      double scale_zvertexreweight = h_zvertex_reweight->GetBinContent(h_zvertex_reweight->GetXaxis()->FindBin(zvtx));

      bgdj = false;
      bgnj = false;

      
      int tlji = -1;
      float ltjpt = -9999;
      for(int j=0; j<tjet_n; ++j)
	{
	  if(abs(tjet_eta[j]) > 1.5-jet_rad) continue;
	  if(tjet_pt[j] > ltjpt)
	    {
	      ltjpt = tjet_pt[j];
	      tlji = j;
	    }
	}
      
      if(tlji < 0) continue;
      if(ltjpt < truthjet_pt_min || ltjpt > truthjet_pt_max) continue;
      for(int j=0; j<tjet_n; ++j)
	{
	  if(abs(tjet_eta[j]) > 1.1-jet_rad) continue;
	  h_truthjet_pt_record_all->Fill(tjet_pt[j]);
	}


      int lji = -1;
      int sji = -1;
      get_l_sl_jet(lji,sji,jet_e,jet_n);
      if(sji < 0) continue;
      if(jet_pt[lji] > recojet_pt_max) continue;

      bool dijetcut = !check_dphicut(jet_phi[lji],jet_phi[sji]) || jet_e[sji]/jet_e[lji]<0.3;
      
      bgdj = dijetcut;

      int njg5 = 0;
      for(int j=0; j<jet_n; ++j)
	{
	  if(jet_e[j]>5) ++njg5;
	}
      if(njg5>9) bgnj = true;

      

      double ljpt = jet_pt[lji];

      double meff[3];

      for(int j=0; j<3; ++j)
	{
	  meff[j] = mbt[j]->Eval(ljpt);
	  if(meff[j]<0.01) meff[j]=0.01;
	}
      double mbdtrig_scale_nominal = 1.0 / meff[0];

      
      filter_jets(jet_filter, jet_e, jet_pt, jet_eta, zvtx, jet_n);
      filter_tjets(tjet_filter, tjet_e, tjet_pt, tjet_eta, zvtx, tjet_n);
      
      for(int j=0; j<jet_filter.size(); ++j)
	{
	  if(jet_filter.at(j)) continue;

	  h_recojet_pt_record_nocut_all->Fill(jet_pt[j]);
	  if(!bgnj && !bgdj) h_recojet_pt_record_all->Fill(jet_pt[j]);

	  if(z60)
	    {
	      h_recojet_pt_record_nocut_zvertex60->Fill(jet_pt[j]);
	      if(!bgnj && !bgdj) h_recojet_pt_record_zvertex60->Fill(jet_pt[j]);
	      if(z30)
		{
		  h_recojet_pt_record_nocut_zvertex30->Fill(jet_pt[j]);
		  if(!bgnj && !bgdj) h_recojet_pt_record_zvertex30->Fill(jet_pt[j]);
		}
	    }
	}
      
      if(bgnj || bgdj) isbg = true;
      else isbg = false;
      h_event_passed->Fill(0.5);
      get_truthjet(goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, tjet_filter, tjet_n, tjet_pt, tjet_eta, tjet_phi);
      get_calibjet(calibjet_pt, calibjet_eta, calibjet_phi, jet_filter, calib_jet_n, jet_pt_calib, jet_eta, jet_phi, 1, 0.1);

      get_calibjet(calibjet_pt_nosmear, calibjet_eta_nosmear, calibjet_phi_nosmear, jet_filter, calib_jet_n, jet_pt_calib, jet_eta, jet_phi, 1, 0);
      
      get_calibjet(calibjet_pt_jesup, calibjet_eta_jesup, calibjet_phi_jesup, jet_filter, calib_jet_n, jet_pt_calib, jet_eta, jet_phi, 1.06, 0.1);
      get_calibjet(calibjet_pt_jesdown, calibjet_eta_jesdown, calibjet_phi_jesdown, jet_filter, calib_jet_n, jet_pt_calib, jet_eta, jet_phi, 0.94, 0.1);
      get_calibjet(calibjet_pt_jerup, calibjet_eta_jerup, calibjet_phi_jerup, jet_filter, calib_jet_n, jet_pt_calib, jet_eta, jet_phi, 1, 0.15);
      get_calibjet(calibjet_pt_jerdown, calibjet_eta_jerdown, calibjet_phi_jerdown, jet_filter, calib_jet_n, jet_pt_calib, jet_eta, jet_phi, 1, 0.05);
      
      sw = 1;
      match_meas_truth(calibjet_eta, calibjet_phi, calibjet_matched, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_rad, has_zvtx);
      sw = 0;
      
      fill_response_matrix(h_truth_all, h_measure_all, h_respmatrix_all, h_fake_all, h_miss_all,
                         h_matchedtruth_weighted_all, h_matchedtruth_unweighted_all, h_measure_unweighted_all,
                         scale_zvertexreweight, f_reweightfunc_all,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
      
      
      
      fill_response_matrix(h_truth_all_jetup, h_measure_all_jetup, h_respmatrix_all_jetup, h_fake_all_jetup, h_miss_all_jetup,
                         h_matchedtruth_weighted_all_jetup, h_matchedtruth_unweighted_all_jetup, h_measure_unweighted_all_jetup,
                         scale_zvertexreweight, f_reweightfunc_all_jetup,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
      
      fill_response_matrix(h_truth_all_jetdown, h_measure_all_jetdown, h_respmatrix_all_jetdown, h_fake_all_jetdown, h_miss_all_jetdown,
                         h_matchedtruth_weighted_all_jetdown, h_matchedtruth_unweighted_all_jetdown, h_measure_unweighted_all_jetdown,
                         scale_zvertexreweight, f_reweightfunc_all_jetdown,
                         calibjet_pt, calibjet_matched,
                         goodtruthjet_pt, goodtruthjet_matched);
    
    if (z30) {
      fill_response_matrix(h_truth_zvertex30, h_measure_zvertex30, h_respmatrix_zvertex30, h_fake_zvertex30, h_miss_zvertex30,
                           h_matchedtruth_weighted_zvertex30, h_matchedtruth_unweighted_zvertex30, h_measure_unweighted_zvertex30,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex30_jetup, h_measure_zvertex30_jetup, h_respmatrix_zvertex30_jetup, h_fake_zvertex30_jetup, h_miss_zvertex30_jetup,
                           h_matchedtruth_weighted_zvertex30_jetup, h_matchedtruth_unweighted_zvertex30_jetup, h_measure_unweighted_zvertex30_jetup,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_jetup,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex30_jetdown, h_measure_zvertex30_jetdown, h_respmatrix_zvertex30_jetdown, h_fake_zvertex30_jetdown, h_miss_zvertex30_jetdown,
                           h_matchedtruth_weighted_zvertex30_jetdown, h_matchedtruth_unweighted_zvertex30_jetdown, h_measure_unweighted_zvertex30_jetdown,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_jetdown,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex30_mbdup, h_measure_zvertex30_mbdup, h_respmatrix_zvertex30_mbdup, h_fake_zvertex30_mbdup, h_miss_zvertex30_mbdup,
                           h_matchedtruth_weighted_zvertex30_mbdup, h_matchedtruth_unweighted_zvertex30_mbdup, h_measure_unweighted_zvertex30_mbdup,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_mbdup,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex30_mbddown, h_measure_zvertex30_mbddown, h_respmatrix_zvertex30_mbddown, h_fake_zvertex30_mbddown, h_miss_zvertex30_mbddown,
                           h_matchedtruth_weighted_zvertex30_mbddown, h_matchedtruth_unweighted_zvertex30_mbddown, h_measure_unweighted_zvertex30_mbddown,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_mbddown,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
    }
    if (z60) {


      fill_response_matrix(h_truth_zvertex60, h_measure_zvertex60, h_respmatrix_zvertex60, h_fake_zvertex60, h_miss_zvertex60,
                           h_matchedtruth_weighted_zvertex60, h_matchedtruth_unweighted_zvertex60, h_measure_unweighted_zvertex60,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex60_jetup, h_measure_zvertex60_jetup, h_respmatrix_zvertex60_jetup, h_fake_zvertex60_jetup, h_miss_zvertex60_jetup,
                           h_matchedtruth_weighted_zvertex60_jetup, h_matchedtruth_unweighted_zvertex60_jetup, h_measure_unweighted_zvertex60_jetup,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_jetup,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex60_jetdown, h_measure_zvertex60_jetdown, h_respmatrix_zvertex60_jetdown, h_fake_zvertex60_jetdown, h_miss_zvertex60_jetdown,
                           h_matchedtruth_weighted_zvertex60_jetdown, h_matchedtruth_unweighted_zvertex60_jetdown, h_measure_unweighted_zvertex60_jetdown,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_jetdown,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex60_mbdup, h_measure_zvertex60_mbdup, h_respmatrix_zvertex60_mbdup, h_fake_zvertex60_mbdup, h_miss_zvertex60_mbdup,
                           h_matchedtruth_weighted_zvertex60_mbdup, h_matchedtruth_unweighted_zvertex60_mbdup, h_measure_unweighted_zvertex60_mbdup,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_mbdup,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
      fill_response_matrix(h_truth_zvertex60_mbddown, h_measure_zvertex60_mbddown, h_respmatrix_zvertex60_mbddown, h_fake_zvertex60_mbddown, h_miss_zvertex60_mbddown,
                           h_matchedtruth_weighted_zvertex60_mbddown, h_matchedtruth_unweighted_zvertex60_mbddown, h_measure_unweighted_zvertex60_mbddown,
                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_mbddown,
                           calibjet_pt, calibjet_matched,
                           goodtruthjet_pt, goodtruthjet_matched);
    }
    if (z60) {
      for (int im = 0; im < calibjet_pt.size(); ++im) {
        h_fullclosure_measure_zvertex60->Fill(calibjet_pt[im]);
        if ((i/5) % 2 == 0) h_halfclosure_measure_zvertex60->Fill(calibjet_pt[im]);
        else h_halfclosure_inputmeasure_zvertex60->Fill(calibjet_pt[im]);
        if (calibjet_matched[im] < 0) {
          h_fullclosure_fake_zvertex60->Fill(calibjet_pt[im]);
          if ((i/5) % 2 == 0) h_halfclosure_fake_zvertex60->Fill(calibjet_pt[im]);
        } else {
          h_fullclosure_respmatrix_zvertex60->Fill(calibjet_pt[im], goodtruthjet_pt[calibjet_matched[im]]);
          if ((i/5) % 2 == 0) h_halfclosure_respmatrix_zvertex60->Fill(calibjet_pt[im], goodtruthjet_pt[calibjet_matched[im]]);
        }
      }
      for (int it = 0; it < goodtruthjet_pt.size(); ++it) {
        h_fullclosure_truth_zvertex60->Fill(goodtruthjet_pt[it]);
        if ((i/5) % 2 == 0) h_halfclosure_truth_zvertex60->Fill(goodtruthjet_pt[it]);
        if (goodtruthjet_matched[it] < 0) {
          h_fullclosure_miss_zvertex60->Fill(goodtruthjet_pt[it]);
          if ((i/5) % 2 == 0) h_halfclosure_miss_zvertex60->Fill(goodtruthjet_pt[it]);
        }
      }
    }

    match_meas_truth(calibjet_eta_jesup, calibjet_phi_jesup, calibjet_matched_jesup, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_rad, has_zvtx);
    fill_response_matrix(h_truth_all_jesup, h_measure_all_jesup, h_respmatrix_all_jesup, h_fake_all_jesup, h_miss_all_jesup,
                         h_matchedtruth_weighted_all_jesup, h_matchedtruth_unweighted_all_jesup, h_measure_unweighted_all_jesup,
                         scale_zvertexreweight, f_reweightfunc_all_jesup,
                         calibjet_pt_jesup, calibjet_matched_jesup,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (z30) fill_response_matrix(h_truth_zvertex30_jesup, h_measure_zvertex30_jesup, h_respmatrix_zvertex30_jesup, h_fake_zvertex30_jesup, h_miss_zvertex30_jesup,
                                           h_matchedtruth_weighted_zvertex30_jesup, h_matchedtruth_unweighted_zvertex30_jesup, h_measure_unweighted_zvertex30_jesup,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_jesup,
                                           calibjet_pt_jesup, calibjet_matched_jesup,
                                           goodtruthjet_pt, goodtruthjet_matched);
    if (z60) fill_response_matrix(h_truth_zvertex60_jesup, h_measure_zvertex60_jesup, h_respmatrix_zvertex60_jesup, h_fake_zvertex60_jesup, h_miss_zvertex60_jesup,
                                           h_matchedtruth_weighted_zvertex60_jesup, h_matchedtruth_unweighted_zvertex60_jesup, h_measure_unweighted_zvertex60_jesup,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_jesup,
                                           calibjet_pt_jesup, calibjet_matched_jesup,
                                           goodtruthjet_pt, goodtruthjet_matched);
                                           
    match_meas_truth(calibjet_eta_jesdown, calibjet_phi_jesdown, calibjet_matched_jesdown, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_rad, has_zvtx);
    fill_response_matrix(h_truth_all_jesdown, h_measure_all_jesdown, h_respmatrix_all_jesdown, h_fake_all_jesdown, h_miss_all_jesdown,
                         h_matchedtruth_weighted_all_jesdown, h_matchedtruth_unweighted_all_jesdown, h_measure_unweighted_all_jesdown,
                         scale_zvertexreweight, f_reweightfunc_all_jesdown,
                         calibjet_pt_jesdown, calibjet_matched_jesdown,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (z30) fill_response_matrix(h_truth_zvertex30_jesdown, h_measure_zvertex30_jesdown, h_respmatrix_zvertex30_jesdown, h_fake_zvertex30_jesdown, h_miss_zvertex30_jesdown,
                                           h_matchedtruth_weighted_zvertex30_jesdown, h_matchedtruth_unweighted_zvertex30_jesdown, h_measure_unweighted_zvertex30_jesdown,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_jesdown,
                                           calibjet_pt_jesdown, calibjet_matched_jesdown,
                                           goodtruthjet_pt, goodtruthjet_matched);
    if (z60) fill_response_matrix(h_truth_zvertex60_jesdown, h_measure_zvertex60_jesdown, h_respmatrix_zvertex60_jesdown, h_fake_zvertex60_jesdown, h_miss_zvertex60_jesdown,
                                           h_matchedtruth_weighted_zvertex60_jesdown, h_matchedtruth_unweighted_zvertex60_jesdown, h_measure_unweighted_zvertex60_jesdown,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_jesdown,
                                           calibjet_pt_jesdown, calibjet_matched_jesdown,
                                           goodtruthjet_pt, goodtruthjet_matched);

    match_meas_truth(calibjet_eta_jerup, calibjet_phi_jerup, calibjet_matched_jerup, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_rad, has_zvtx);
    fill_response_matrix(h_truth_all_jerup, h_measure_all_jerup, h_respmatrix_all_jerup, h_fake_all_jerup, h_miss_all_jerup,
                         h_matchedtruth_weighted_all_jerup, h_matchedtruth_unweighted_all_jerup, h_measure_unweighted_all_jerup,
                         scale_zvertexreweight, f_reweightfunc_all_jerup,
                         calibjet_pt_jerup, calibjet_matched_jerup,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (z30) fill_response_matrix(h_truth_zvertex30_jerup, h_measure_zvertex30_jerup, h_respmatrix_zvertex30_jerup, h_fake_zvertex30_jerup, h_miss_zvertex30_jerup,
                                           h_matchedtruth_weighted_zvertex30_jerup, h_matchedtruth_unweighted_zvertex30_jerup, h_measure_unweighted_zvertex30_jerup,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_jerup,
                                           calibjet_pt_jerup, calibjet_matched_jerup,
                                           goodtruthjet_pt, goodtruthjet_matched);
    if (z60) fill_response_matrix(h_truth_zvertex60_jerup, h_measure_zvertex60_jerup, h_respmatrix_zvertex60_jerup, h_fake_zvertex60_jerup, h_miss_zvertex60_jerup,
                                           h_matchedtruth_weighted_zvertex60_jerup, h_matchedtruth_unweighted_zvertex60_jerup, h_measure_unweighted_zvertex60_jerup,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_jerup,
                                           calibjet_pt_jerup, calibjet_matched_jerup,
                                           goodtruthjet_pt, goodtruthjet_matched);

    match_meas_truth(calibjet_eta_jerdown, calibjet_phi_jerdown, calibjet_matched_jerdown, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_rad, has_zvtx);
    fill_response_matrix(h_truth_all_jerdown, h_measure_all_jerdown, h_respmatrix_all_jerdown, h_fake_all_jerdown, h_miss_all_jerdown,
                         h_matchedtruth_weighted_all_jerdown, h_matchedtruth_unweighted_all_jerdown, h_measure_unweighted_all_jerdown,
                         scale_zvertexreweight, f_reweightfunc_all_jerdown,
                         calibjet_pt_jerdown, calibjet_matched_jerdown,
                         goodtruthjet_pt, goodtruthjet_matched);
    if (z30) fill_response_matrix(h_truth_zvertex30_jerdown, h_measure_zvertex30_jerdown, h_respmatrix_zvertex30_jerdown, h_fake_zvertex30_jerdown, h_miss_zvertex30_jerdown,
                                           h_matchedtruth_weighted_zvertex30_jerdown, h_matchedtruth_unweighted_zvertex30_jerdown, h_measure_unweighted_zvertex30_jerdown,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex30_jerdown,
                                           calibjet_pt_jerdown, calibjet_matched_jerdown,
                                           goodtruthjet_pt, goodtruthjet_matched);
    if (z60) fill_response_matrix(h_truth_zvertex60_jerdown, h_measure_zvertex60_jerdown, h_respmatrix_zvertex60_jerdown, h_fake_zvertex60_jerdown, h_miss_zvertex60_jerdown,
                                           h_matchedtruth_weighted_zvertex60_jerdown, h_matchedtruth_unweighted_zvertex60_jerdown, h_measure_unweighted_zvertex60_jerdown,
                                           scale_zvertexreweight*mbdtrig_scale_nominal, f_reweightfunc_zvertex60_jerdown,
                                           calibjet_pt_jerdown, calibjet_matched_jerdown,
                                           goodtruthjet_pt, goodtruthjet_matched);
    match_meas_truth(calibjet_eta_nosmear, calibjet_phi_nosmear, calibjet_matched_nosmear, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_rad, has_zvtx);
    fill_response_matrix(h_truth_all_nosmear, h_measure_all_nosmear, h_respmatrix_all_nosmear, h_fake_all_nosmear, h_miss_all_nosmear,
                         h_matchedtruth_weighted_all_nosmear, h_matchedtruth_unweighted_all_nosmear, h_measure_unweighted_all_nosmear,
                         scale_zvertexreweight, f_reweightfunc_all,
                         calibjet_pt_nosmear, calibjet_matched_nosmear,
                         goodtruthjet_pt, goodtruthjet_matched);
    if(z60)       fill_response_matrix(h_truth_zvertex60_nosmear, h_measure_zvertex60_nosmear, h_respmatrix_zvertex60_nosmear, h_fake_zvertex60_nosmear, h_miss_zvertex60_nosmear,
			   h_matchedtruth_weighted_zvertex60_nosmear, h_matchedtruth_unweighted_zvertex60_nosmear, h_measure_unweighted_zvertex60_nosmear,
			   scale_zvertexreweight, f_reweightfunc_zvertex60,
			   calibjet_pt_nosmear, calibjet_matched_nosmear,
			   goodtruthjet_pt, goodtruthjet_matched);


  } // event loop end
  // Fill event histograms.
  h_event_all->Fill(0.5, chain.GetEntries());

  TFile* outf = TFile::Open(("output/output_"+runtype+"_"+to_string(iseg)+"_"+to_string(iseg+nseg)+".root").c_str(),"RECREATE");
  outf->cd();
  h_event_all->Write();
  h_event_beforecut->Write();
  h_event_passed->Write();
  h_recojet_pt_record_nocut_all->Write();
  h_recojet_pt_record_all->Write();
  h_recojet_pt_record_nocut_zvertex30->Write();
  h_recojet_pt_record_zvertex30->Write();
  h_recojet_pt_record_nocut_zvertex60->Write();
  h_recojet_pt_record_zvertex60->Write();
  h_truthjet_pt_record_all->Write();

  // Nominal histograms
  h_truth_all->Write(); h_measure_all->Write(); h_respmatrix_all->Write(); h_fake_all->Write(); h_miss_all->Write(); h_matchedtruth_weighted_all->Write(); h_matchedtruth_unweighted_all->Write(); h_measure_unweighted_all->Write();

  h_truth_all_nosmear->Write(); h_measure_all_nosmear->Write(); h_respmatrix_all_nosmear->Write(); h_fake_all_nosmear->Write(); h_miss_all_nosmear->Write(); h_matchedtruth_weighted_all_nosmear->Write(); h_matchedtruth_unweighted_all_nosmear->Write(); h_measure_unweighted_all_nosmear->Write();


    h_truth_zvertex60_nosmear->Write(); h_measure_zvertex60_nosmear->Write(); h_respmatrix_zvertex60_nosmear->Write(); h_fake_zvertex60_nosmear->Write(); h_miss_zvertex60_nosmear->Write(); h_matchedtruth_weighted_zvertex60_nosmear->Write(); h_matchedtruth_unweighted_zvertex60_nosmear->Write(); h_measure_unweighted_zvertex60_nosmear->Write();

  
  h_truth_zvertex30->Write(); h_measure_zvertex30->Write(); h_respmatrix_zvertex30->Write(); h_fake_zvertex30->Write(); h_miss_zvertex30->Write(); h_matchedtruth_weighted_zvertex30->Write(); h_matchedtruth_unweighted_zvertex30->Write(); h_measure_unweighted_zvertex30->Write();
  h_truth_zvertex60->Write(); h_measure_zvertex60->Write(); h_respmatrix_zvertex60->Write(); h_fake_zvertex60->Write(); h_miss_zvertex60->Write(); h_matchedtruth_weighted_zvertex60->Write(); h_matchedtruth_unweighted_zvertex60->Write(); h_measure_unweighted_zvertex60->Write();
  // JES up/down histograms
  h_truth_all_jesup->Write(); h_measure_all_jesup->Write(); h_respmatrix_all_jesup->Write(); h_fake_all_jesup->Write(); h_miss_all_jesup->Write(); h_matchedtruth_weighted_all_jesup->Write(); h_matchedtruth_unweighted_all_jesup->Write(); h_measure_unweighted_all_jesup->Write();
  h_truth_all_jesdown->Write(); h_measure_all_jesdown->Write(); h_respmatrix_all_jesdown->Write(); h_fake_all_jesdown->Write(); h_miss_all_jesdown->Write(); h_matchedtruth_weighted_all_jesdown->Write(); h_matchedtruth_unweighted_all_jesdown->Write(); h_measure_unweighted_all_jesdown->Write();
  h_truth_zvertex30_jesup->Write(); h_measure_zvertex30_jesup->Write(); h_respmatrix_zvertex30_jesup->Write(); h_fake_zvertex30_jesup->Write(); h_miss_zvertex30_jesup->Write(); h_matchedtruth_weighted_zvertex30_jesup->Write(); h_matchedtruth_unweighted_zvertex30_jesup->Write(); h_measure_unweighted_zvertex30_jesup->Write();
  h_truth_zvertex30_jesdown->Write(); h_measure_zvertex30_jesdown->Write(); h_respmatrix_zvertex30_jesdown->Write(); h_fake_zvertex30_jesdown->Write(); h_miss_zvertex30_jesdown->Write(); h_matchedtruth_weighted_zvertex30_jesdown->Write(); h_matchedtruth_unweighted_zvertex30_jesdown->Write(); h_measure_unweighted_zvertex30_jesdown->Write();
  h_truth_zvertex60_jesup->Write(); h_measure_zvertex60_jesup->Write(); h_respmatrix_zvertex60_jesup->Write(); h_fake_zvertex60_jesup->Write(); h_miss_zvertex60_jesup->Write(); h_matchedtruth_weighted_zvertex60_jesup->Write(); h_matchedtruth_unweighted_zvertex60_jesup->Write(); h_measure_unweighted_zvertex60_jesup->Write();
  h_truth_zvertex60_jesdown->Write(); h_measure_zvertex60_jesdown->Write(); h_respmatrix_zvertex60_jesdown->Write(); h_fake_zvertex60_jesdown->Write(); h_miss_zvertex60_jesdown->Write(); h_matchedtruth_weighted_zvertex60_jesdown->Write(); h_matchedtruth_unweighted_zvertex60_jesdown->Write(); h_measure_unweighted_zvertex60_jesdown->Write();
  // JER up/down histograms
  h_truth_all_jerup->Write(); h_measure_all_jerup->Write(); h_respmatrix_all_jerup->Write(); h_fake_all_jerup->Write(); h_miss_all_jerup->Write(); h_matchedtruth_weighted_all_jerup->Write(); h_matchedtruth_unweighted_all_jerup->Write(); h_measure_unweighted_all_jerup->Write();
  h_truth_all_jerdown->Write(); h_measure_all_jerdown->Write(); h_respmatrix_all_jerdown->Write(); h_fake_all_jerdown->Write(); h_miss_all_jerdown->Write(); h_matchedtruth_weighted_all_jerdown->Write(); h_matchedtruth_unweighted_all_jerdown->Write(); h_measure_unweighted_all_jerdown->Write();
  h_truth_zvertex30_jerup->Write(); h_measure_zvertex30_jerup->Write(); h_respmatrix_zvertex30_jerup->Write(); h_fake_zvertex30_jerup->Write(); h_miss_zvertex30_jerup->Write(); h_matchedtruth_weighted_zvertex30_jerup->Write(); h_matchedtruth_unweighted_zvertex30_jerup->Write(); h_measure_unweighted_zvertex30_jerup->Write();
  h_truth_zvertex30_jerdown->Write(); h_measure_zvertex30_jerdown->Write(); h_respmatrix_zvertex30_jerdown->Write(); h_fake_zvertex30_jerdown->Write(); h_miss_zvertex30_jerdown->Write(); h_matchedtruth_weighted_zvertex30_jerdown->Write(); h_matchedtruth_unweighted_zvertex30_jerdown->Write(); h_measure_unweighted_zvertex30_jerdown->Write();
  h_truth_zvertex60_jerup->Write(); h_measure_zvertex60_jerup->Write(); h_respmatrix_zvertex60_jerup->Write(); h_fake_zvertex60_jerup->Write(); h_miss_zvertex60_jerup->Write(); h_matchedtruth_weighted_zvertex60_jerup->Write(); h_matchedtruth_unweighted_zvertex60_jerup->Write(); h_measure_unweighted_zvertex60_jerup->Write();
  h_truth_zvertex60_jerdown->Write(); h_measure_zvertex60_jerdown->Write(); h_respmatrix_zvertex60_jerdown->Write(); h_fake_zvertex60_jerdown->Write(); h_miss_zvertex60_jerdown->Write(); h_matchedtruth_weighted_zvertex60_jerdown->Write(); h_matchedtruth_unweighted_zvertex60_jerdown->Write(); h_measure_unweighted_zvertex60_jerdown->Write();
  // Jet trigger up/down histograms
  h_truth_all_jetup->Write(); h_measure_all_jetup->Write(); h_respmatrix_all_jetup->Write(); h_fake_all_jetup->Write(); h_miss_all_jetup->Write(); h_matchedtruth_weighted_all_jetup->Write(); h_matchedtruth_unweighted_all_jetup->Write(); h_measure_unweighted_all_jetup->Write();
  h_truth_all_jetdown->Write(); h_measure_all_jetdown->Write(); h_respmatrix_all_jetdown->Write(); h_fake_all_jetdown->Write(); h_miss_all_jetdown->Write(); h_matchedtruth_weighted_all_jetdown->Write(); h_matchedtruth_unweighted_all_jetdown->Write(); h_measure_unweighted_all_jetdown->Write();
  h_truth_zvertex30_jetup->Write(); h_measure_zvertex30_jetup->Write(); h_respmatrix_zvertex30_jetup->Write(); h_fake_zvertex30_jetup->Write(); h_miss_zvertex30_jetup->Write(); h_matchedtruth_weighted_zvertex30_jetup->Write(); h_matchedtruth_unweighted_zvertex30_jetup->Write(); h_measure_unweighted_zvertex30_jetup->Write();
  h_truth_zvertex30_jetdown->Write(); h_measure_zvertex30_jetdown->Write(); h_respmatrix_zvertex30_jetdown->Write(); h_fake_zvertex30_jetdown->Write(); h_miss_zvertex30_jetdown->Write(); h_matchedtruth_weighted_zvertex30_jetdown->Write(); h_matchedtruth_unweighted_zvertex30_jetdown->Write(); h_measure_unweighted_zvertex30_jetdown->Write();
  h_truth_zvertex60_jetup->Write(); h_measure_zvertex60_jetup->Write(); h_respmatrix_zvertex60_jetup->Write(); h_fake_zvertex60_jetup->Write(); h_miss_zvertex60_jetup->Write(); h_matchedtruth_weighted_zvertex60_jetup->Write(); h_matchedtruth_unweighted_zvertex60_jetup->Write(); h_measure_unweighted_zvertex60_jetup->Write();
  h_truth_zvertex60_jetdown->Write(); h_measure_zvertex60_jetdown->Write(); h_respmatrix_zvertex60_jetdown->Write(); h_fake_zvertex60_jetdown->Write(); h_miss_zvertex60_jetdown->Write(); h_matchedtruth_weighted_zvertex60_jetdown->Write(); h_matchedtruth_unweighted_zvertex60_jetdown->Write(); h_measure_unweighted_zvertex60_jetdown->Write();
  // MBD trigger up/down histograms
  h_truth_zvertex30_mbdup->Write(); h_measure_zvertex30_mbdup->Write(); h_respmatrix_zvertex30_mbdup->Write(); h_fake_zvertex30_mbdup->Write(); h_miss_zvertex30_mbdup->Write(); h_matchedtruth_weighted_zvertex30_mbdup->Write(); h_matchedtruth_unweighted_zvertex30_mbdup->Write(); h_measure_unweighted_zvertex30_mbdup->Write();
  h_truth_zvertex30_mbddown->Write(); h_measure_zvertex30_mbddown->Write(); h_respmatrix_zvertex30_mbddown->Write(); h_fake_zvertex30_mbddown->Write(); h_miss_zvertex30_mbddown->Write(); h_matchedtruth_weighted_zvertex30_mbddown->Write(); h_matchedtruth_unweighted_zvertex30_mbddown->Write(); h_measure_unweighted_zvertex30_mbddown->Write();
  h_truth_zvertex60_mbdup->Write(); h_measure_zvertex60_mbdup->Write(); h_respmatrix_zvertex60_mbdup->Write(); h_fake_zvertex60_mbdup->Write(); h_miss_zvertex60_mbdup->Write(); h_matchedtruth_weighted_zvertex60_mbdup->Write(); h_matchedtruth_unweighted_zvertex60_mbdup->Write(); h_measure_unweighted_zvertex60_mbdup->Write();
  h_truth_zvertex60_mbddown->Write(); h_measure_zvertex60_mbddown->Write(); h_respmatrix_zvertex60_mbddown->Write(); h_fake_zvertex60_mbddown->Write(); h_miss_zvertex60_mbddown->Write(); h_matchedtruth_weighted_zvertex60_mbddown->Write(); h_matchedtruth_unweighted_zvertex60_mbddown->Write(); h_measure_unweighted_zvertex60_mbddown->Write();

  // Full closure test histograms
  h_fullclosure_truth_zvertex60->Write();
  h_fullclosure_measure_zvertex60->Write();
  h_fullclosure_respmatrix_zvertex60->Write();
  h_fullclosure_fake_zvertex60->Write();
  h_fullclosure_miss_zvertex60->Write();

  h_halfclosure_truth_zvertex60->Write();
  h_halfclosure_inputmeasure_zvertex60->Write();
  h_halfclosure_measure_zvertex60->Write();
  h_halfclosure_respmatrix_zvertex60->Write();
  h_halfclosure_fake_zvertex60->Write();
  h_halfclosure_miss_zvertex60->Write();
  outf->Close();
  return 0;
}
