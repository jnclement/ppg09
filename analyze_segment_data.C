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
  float dphi = phi1-phi2;
  float dphim = phi1-phi2-2*M_PI;
  float dphip = phi1-phi2+2*M_PI;

  if(abs(dphim) < abs(dphi)) dphi = dphim;
  if(abs(dphip) < abs(dphi)) dphi = dphip;

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

void check_dphicut(float lphi, float sphi)
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

int analyze_segment_data(int nunnumber, int iseg, int nseg)
{
  TChain chain("jet_tree");
  
