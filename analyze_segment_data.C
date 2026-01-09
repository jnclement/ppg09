#include "unfold_Def.h"

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

