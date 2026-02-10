#include "unfold_Def.h"
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
  float dphi = abs(phi1-phi2);
  if(dphi > M_PI) dphi = 2*M_PI-dphi;

  return dphi;
}

float get_deta(float eta1, float eta2)
{
  return eta1-eta2;
}

float get_dR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = get_deta(eta1, eta2);
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

int analyze_segment_data(int iseg, int nseg)
{
  TChain chain("jet_tree");
  
  for(int i=0; i<nseg; ++i)
    {
      chain.Add(("input/inputfile_"+to_string(i)+".root").c_str());
    }

  float zvtx;
  long long unsigned int trigvec = 0;
  float jet_e[100];
  float jet_pt[100];
  float jet_et[100];
  float jet_eta[100];
  float jet_phi[100];
  float jet_t[100];
  float jet_pt_calib[100];
  int jet_n;
  int calib_jet_n;

  chain.SetBranchAddress("zvtx",&zvtx);
  chain.SetBranchAddress("triggervec",&trigvec);
  chain.SetBranchAddress("jet_et",jet_e);
  chain.SetBranchAddress("jet_pt",jet_pt);
  chain.SetBranchAddress("jet_etrans",jet_et);
  chain.SetBranchAddress("jet_eta",jet_eta);
  chain.SetBranchAddress("jet_phi",jet_phi);
  chain.SetBranchAddress("jet_t",jet_t);
  chain.SetBranchAddress("jet_pt_calib",jet_pt_calib);
  chain.SetBranchAddress("jet_n",&jet_n);
  chain.SetBranchAddress("calib_jet_n",&calib_jet_n);

  TFile* fjt = TFile::Open("output_jetefficiency.root","READ");
  TF1* jt[3];
  jt[0] = (TF1*)fjt->Get("jettrig_04_pt_nominal");
  jt[1] = (TF1*)fjt->Get("jettrig_04_pt_up");
  jt[2] = (TF1*)fjt->Get("jettrig_04_pt_down");

  TFile* fmb = TFile::Open("output_mbdefficiency.root","READ");
  TF1* mbt[3];
  mbt[0] = (TF1*)fmb->Get("mbdtrig04_nominal");
  mbt[1] = (TF1*)fmb->Get("mbdtrig04_up");
  mbt[2] = (TF1*)fmb->Get("mbdtrig04_down");

  TH1D *h_event_all = new TH1D("h_event_all", ";Event Number", 1, 0, 1);
  TH1D *h_event_beforecut = new TH1D("h_event_beforecut", ";Event Number", 1, 0, 1);
  TH1D *h_event_passed = new TH1D("h_event_passed", ";Event Number", 1, 0, 1);

  TH1D *h_recojet_pt_record_nocut_all = new TH1D("h_recojet_pt_record_nocut_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_all = new TH1D("h_recojet_pt_record_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_nocut_zvertex30 = new TH1D("h_recojet_pt_record_nocut_zvertex30", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_zvertex30 = new TH1D("h_recojet_pt_record_zvertex30", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_nocut_zvertex60 = new TH1D("h_recojet_pt_record_nocut_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_zvertex60 = new TH1D("h_recojet_pt_record_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);

  TH1D *h_calibjet_pt_all = new TH1D("h_calibjet_pt_all", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_all_jetup = new TH1D("h_calibjet_pt_all_jetup", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_all_jetdown = new TH1D("h_calibjet_pt_all_jetdown", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_record_all = new TH1D("h_calibjet_pt_record_all", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_all_jetup = new TH1D("h_calibjet_pt_record_all_jetup", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_all_jetdown = new TH1D("h_calibjet_pt_record_all_jetdown", ";p_{T} [GeV]", 1000, 0, 100);
  
  TH1D *h_calibjet_pt_zvertex30 = new TH1D("h_calibjet_pt_zvertex30", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex30_jetup = new TH1D("h_calibjet_pt_zvertex30_jetup", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex30_jetdown = new TH1D("h_calibjet_pt_zvertex30_jetdown", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex30_mbdup = new TH1D("h_calibjet_pt_zvertex30_mbdup", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex30_mbddown = new TH1D("h_calibjet_pt_zvertex30_mbddown", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_record_zvertex30 = new TH1D("h_calibjet_pt_record_zvertex30", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex30_jetup = new TH1D("h_calibjet_pt_record_zvertex30_jetup", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex30_jetdown = new TH1D("h_calibjet_pt_record_zvertex30_jetdown", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex30_mbdup = new TH1D("h_calibjet_pt_record_zvertex30_mbdup", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex30_mbddown = new TH1D("h_calibjet_pt_record_zvertex30_mbddown", ";p_{T} [GeV]", 1000, 0, 100);  

  TH1D *h_calibjet_pt_zvertex60 = new TH1D("h_calibjet_pt_zvertex60", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex60_jetup = new TH1D("h_calibjet_pt_zvertex60_jetup", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex60_jetdown = new TH1D("h_calibjet_pt_zvertex60_jetdown", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex60_mbdup = new TH1D("h_calibjet_pt_zvertex60_mbdup", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_zvertex60_mbddown = new TH1D("h_calibjet_pt_zvertex60_mbddown", ";p_{T} [GeV]", calibnpt, calibptbins);
  TH1D *h_calibjet_pt_record_zvertex60 = new TH1D("h_calibjet_pt_record_zvertex60", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex60_jetup = new TH1D("h_calibjet_pt_record_zvertex60_jetup", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex60_jetdown = new TH1D("h_calibjet_pt_record_zvertex60_jetdown", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex60_mbdup = new TH1D("h_calibjet_pt_record_zvertex60_mbdup", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record_zvertex60_mbddown = new TH1D("h_calibjet_pt_record_zvertex60_mbddown", ";p_{T} [GeV]", 1000, 0, 100);

  long long unsigned int nevt = chain.GetEntries();
  bool bgdj, bgnj, z30, z60, is22, is18;
  vector<bool> jet_filter = {};

  for(long long unsigned int i=0; i<nevt; ++i)
    {
      chain.GetEntry(i);
      is22 = (trigvec >> 22) & 1;
      //is18 = (trigvec >> 18) & 1;
      if(!is22) continue;
      z30 = abs(zvtx) < 30 && zvtx!=0;
      z60 = abs(zvtx) < 60 && zvtx!=0;
      if(abs(zvtx)>990) zvtx = 0;

      bgdj = false;
      bgnj = false;

      int lji = -1;
      int sji = -1;
      get_l_sl_jet(lji,sji,jet_e,jet_n);
      if(lji < 0) continue;

      bool dijetcut = !check_dphicut(jet_phi[lji],jet_phi[sji]) || jet_e[sji]/jet_e[lji]<0.3 || sji < 0;
      bool tcut = abs(jet_t[lji]*17.6+2)>6 || abs(jet_t[lji]*17.6-jet_t[sji]*17.6)>3;

      bgdj = tcut || dijetcut;

      int njg5 = 0;
      for(int j=0; j<jet_n; ++j)
	{
	  if(jet_e[j]>5) ++njg5;
	}
      if(njg5>9) bgnj = true;

      

      double ljpt = jet_pt[lji];

      double jeff[3];
      double meff[3];

      for(int j=0; j<3; ++j)
	{
	  jeff[j] = jt[j]->Eval(ljpt);
	  if(jeff[j]<0.01) jeff[j]=0.01;
	  meff[j] = mbt[j]->Eval(ljpt);
	  if(meff[j]<0.01) meff[j]=0.01;
	}

      double maxeff = 0.95;
      filter_jets(jet_filter, jet_e, jet_pt, jet_eta, zvtx, jet_n);

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
      if(bgnj || bgdj) continue;
      h_event_passed->Fill(0.5);
      for(int j=0; j<jet_n; ++j)
	{
	  if(jet_filter.at(j)) continue;
	  h_calibjet_pt_all->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/maxeff);
          h_calibjet_pt_all_jetup->Fill(jet_pt_calib[j], 1.0/jeff[1]*1.0/maxeff);
          h_calibjet_pt_all_jetdown->Fill(jet_pt_calib[j], 1.0/jeff[2]*1.0/maxeff);
          h_calibjet_pt_record_all->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/maxeff);
          h_calibjet_pt_record_all_jetup->Fill(jet_pt_calib[j], 1.0/jeff[1]*1.0/maxeff);
          h_calibjet_pt_record_all_jetdown->Fill(jet_pt_calib[j], 1.0/jeff[2]*1.0/maxeff);

	  if(z60)
	    {
	      h_calibjet_pt_zvertex60->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[0]*1.0/maxeff);
	      h_calibjet_pt_zvertex60_jetup->Fill(jet_pt_calib[j], 1.0/jeff[1]*1.0/meff[0]*1.0/maxeff);
	      h_calibjet_pt_zvertex60_jetdown->Fill(jet_pt_calib[j], 1.0/jeff[2]*1.0/meff[0]*1.0/maxeff);
	      h_calibjet_pt_zvertex60_mbdup->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[1]*1.0/maxeff);
	      h_calibjet_pt_zvertex60_mbddown->Fill(jet_pt_calib[j], 1.0/jeff[0]* 1.0/meff[2]*1.0/maxeff);
	      h_calibjet_pt_record_zvertex60->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[0]*1.0/maxeff);
	      h_calibjet_pt_record_zvertex60_jetup->Fill(jet_pt_calib[j], 1.0/jeff[1]*1.0/meff[0]*1.0/maxeff);
	      h_calibjet_pt_record_zvertex60_jetdown->Fill(jet_pt_calib[j], 1.0/jeff[2]*1.0/meff[0]*1.0/maxeff);
	      h_calibjet_pt_record_zvertex60_mbdup->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[1]*1.0/maxeff);
	      h_calibjet_pt_record_zvertex60_mbddown->Fill(jet_pt_calib[j], 1.0/jeff[0]* 1.0/meff[2]*1.0/maxeff);
	      if(z30)
		{
		  h_calibjet_pt_zvertex30->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[0]*1.0/maxeff);
		  h_calibjet_pt_zvertex30_jetup->Fill(jet_pt_calib[j], 1.0/jeff[1]*1.0/meff[0]*1.0/maxeff);
		  h_calibjet_pt_zvertex30_jetdown->Fill(jet_pt_calib[j], 1.0/jeff[2]*1.0/meff[0]*1.0/maxeff);
		  h_calibjet_pt_zvertex30_mbdup->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[1]*1.0/maxeff);
		  h_calibjet_pt_zvertex30_mbddown->Fill(jet_pt_calib[j], 1.0/jeff[0]* 1.0/meff[2]*1.0/maxeff);
		  h_calibjet_pt_record_zvertex30->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[0]*1.0/maxeff);
		  h_calibjet_pt_record_zvertex30_jetup->Fill(jet_pt_calib[j], 1.0/jeff[1]*1.0/meff[0]*1.0/maxeff);
		  h_calibjet_pt_record_zvertex30_jetdown->Fill(jet_pt_calib[j], 1.0/jeff[2]*1.0/meff[0]*1.0/maxeff);
		  h_calibjet_pt_record_zvertex30_mbdup->Fill(jet_pt_calib[j], 1.0/jeff[0]*1.0/meff[1]*1.0/maxeff);
		  h_calibjet_pt_record_zvertex30_mbddown->Fill(jet_pt_calib[j], 1.0/jeff[0]* 1.0/meff[2]*1.0/maxeff);
		}
	    }
	}
    }
  h_event_all->Fill(0.5,chain.GetEntries());
  TFile* outf = TFile::Open(("output/output_"+to_string(iseg)+".root").c_str(),"RECREATE");
  outf->cd();

  h_recojet_pt_record_nocut_all->Write();
  h_recojet_pt_record_all->Write();
  h_recojet_pt_record_nocut_zvertex30->Write();
  h_recojet_pt_record_zvertex30->Write();
  h_recojet_pt_record_nocut_zvertex60->Write();
  h_recojet_pt_record_zvertex60->Write();

  h_calibjet_pt_all->Write();
  h_calibjet_pt_all_jetup->Write();
  h_calibjet_pt_all_jetdown->Write();
  h_calibjet_pt_record_all->Write();
  h_calibjet_pt_record_all_jetup->Write();
  h_calibjet_pt_record_all_jetdown->Write();

  h_calibjet_pt_zvertex30->Write();
  h_calibjet_pt_zvertex30_jetup->Write();
  h_calibjet_pt_zvertex30_jetdown->Write();
  h_calibjet_pt_zvertex30_mbdup->Write();
  h_calibjet_pt_zvertex30_mbddown->Write();
  h_calibjet_pt_record_zvertex30->Write();
  h_calibjet_pt_record_zvertex30_jetup->Write();
  h_calibjet_pt_record_zvertex30_jetdown->Write();
  h_calibjet_pt_record_zvertex30_mbdup->Write();
  h_calibjet_pt_record_zvertex30_mbddown->Write();

  h_calibjet_pt_zvertex60->Write();
  h_calibjet_pt_zvertex60_jetup->Write();
  h_calibjet_pt_zvertex60_jetdown->Write();
  h_calibjet_pt_zvertex60_mbdup->Write();
  h_calibjet_pt_zvertex60_mbddown->Write();
  h_calibjet_pt_record_zvertex60->Write();
  h_calibjet_pt_record_zvertex60_jetup->Write();
  h_calibjet_pt_record_zvertex60_jetdown->Write();
  h_calibjet_pt_record_zvertex60_mbdup->Write();
  h_calibjet_pt_record_zvertex60_mbddown->Write();

  outf->Close();

  return 0;
  
}
