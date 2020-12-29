// layout_minsurface.cpp
// Andrei Alexandru
// Jan 2010

#include "layout.h"
#include "layout_minsurface.h"
#include <vector>
#include <assert.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "terminate.h"

namespace qcd {

layout_minsurface::layout_minsurface(
    unsigned short nx, unsigned short ny,
    unsigned short nz, unsigned short nt,
    unsigned int np,
    unsigned short sub_x,unsigned short sub_y,unsigned short sub_z) :
  surface_minimizer(np, dimensions(nx,ny,nz,nt),sub_x,sub_y,sub_z),
  layout_subdividing(surface_minimizer::_hyp_lat, surface_minimizer::_sub_lat) { }

layout_minsurface_eo::layout_minsurface_eo(
    unsigned short nx, unsigned short ny,
    unsigned short nz, unsigned short nt,
    unsigned int np,
    unsigned short sub_x,unsigned short sub_y,unsigned short sub_z
    ) :
  surface_minimizer(np, dimensions(nx,ny,nz,nt),sub_x,sub_y,sub_z),
  layout_subdividing_eo(surface_minimizer::_hyp_lat, surface_minimizer::_sub_lat)
{
  assert(volume%(2*np)==0);
  if( ((nx%2) != 0)  || ((ny%2) !=0) || ((nz%2) != 0) || ((nt%2) != 0) )
  {
    if(get_node_rank() == 0) printf("For e-o layout all dimensions have to be even ... exiting.\n");
    terminate_error(ERROR);
  }
}
unsigned int layout_subdividing::get_site_index(const position& p) const  
{
  position lp(p);
  lp.normalize(sub_lat);
  return lp.index(sub_lat);
}

rank_t layout_subdividing::get_site_rank(const position& p) const 
{
  position lp;
  for(int i=0; i<4; ++i) lp[i] = p[i] / sub_lat[i];
  return lp.index(hyp_lat);
}

position layout_subdividing::get_position(local_idx_t idx, rank_t rnk) const
{
  position p(rnk, hyp_lat);
  for(int i=0; i<4; ++i) p[i] *= sub_lat[i];
  position dp(idx, sub_lat);
  return p + dp;
}


local_idx_t layout_subdividing_eo::get_site_index(const position& p) const
{
  position lp(p);
  lp.normalize(sub_lat);
  int par = p.parity();
  int vol2 = sites_on_node/2;
  return par*vol2 + lp.index(sub_lat)/2;
}

rank_t layout_subdividing_eo::get_site_rank(const position& p) const
{
  position lp;
  for(int i=0; i<4; ++i) lp[i] = p[i] / sub_lat[i];
  return lp.index(hyp_lat);
}

position layout_subdividing_eo::get_position(local_idx_t idx, rank_t rnk) const
{
  position p(rnk, hyp_lat);
  for(int i=0; i<4; ++i) p[i] *= sub_lat[i];
  int parity = idx/(sites_on_node/2);
  position dp((idx%(sites_on_node/2))*2, sub_lat);
  p = p + dp;
  if(parity != p.parity())
  {
    position dp2((idx%(sites_on_node/2))*2+1, sub_lat);
    p = p + dp2 - dp;
  }
  return p; 
}

static int get_first_prime(int num)
{
  for(int i=2; i*i<=num; ++i) if((num%i)==0) return i;
  return num;
}

static void factor(int num, std::vector<int>& fact)
{
  assert(num>0);
  while(num>1)
  {
    int prime = get_first_prime(num);
    fact.push_back(prime);
    num /= prime;
  }
}

static void div(int* digits, std::vector<int>& fact, dimensions & divs)
{
  divs[0] = divs[1] = divs[2] = divs[3] = 1;
  int len = fact.size();
  for(int i=0; i<len; ++i) divs[digits[i]] *= fact[i];
}

static bool isdiv(dimensions dims, dimensions divs)
{
  for(int i=0; i<4; ++i) if (dims[i]%divs[i] !=0) return false;
  return true;
}

static unsigned int surface(const dimensions &dims, const dimensions &divs)
{
  int sdims[4];
  int vol = 1;
  for(int i=0; i<4; ++i) 
  {
    sdims[i] = dims[i]/divs[i];
    vol *= sdims[i];
  } 

  int res = 0;
  for(int i=0; i<4; ++i)
    if(divs[i]>1) res += 2*vol/sdims[i];

  return res;
}

// returns true when cut1 is better than cut2
// we assume that the surfaces are equal
// The better cut is the one with larger x(0) dimension
static bool compare_cuts(const dimensions &cut1, const dimensions &cut2, const dimensions &dims)
{
  for(int i=0; i<3; ++i)
  {
    int dim1 = dims[i]/cut1[i];
    int dim2 = dims[i]/cut2[i];
    if(dim1 > dim2) return true;
    if(dim1 < dim2) return false;
  }
  return true; // we should never get here
}

surface_minimizer::surface_minimizer(unsigned int np, dimensions lat_size,
            unsigned short sub_x,unsigned short sub_y,unsigned short sub_z) : lat_sz(lat_size), 
  _hyp_lat(1,1,1,1), _sub_lat(lat_size) 
{
  if(sub_x>0&&sub_y>0&&sub_z>0&&(np%(sub_x*sub_y*sub_z)==0))
  {
      _hyp_lat[0]=sub_x;
      _hyp_lat[1]=sub_y;
      _hyp_lat[2]=sub_z;
      _hyp_lat[3]=np/(sub_x*sub_y*sub_z);
      for(int i=0; i<4; ++i) _sub_lat[i] = lat_size[i] / _hyp_lat[i];
      return;
  }

  std::vector<int> fact;
  factor(np, fact);

  int len = fact.size();
  if(len == 0) return;
  std::vector<int> digits(len);
  for(int i=0; i<len; ++i) digits[i] = 0;
  unsigned int min_surf = 2*lat_sz.compute_volume();
  dimensions divs, dims = lat_sz;

  for(int i=(int)pow(4,len); i>0; --i)
  {
    div(&digits[0], fact, divs);
    if(isdiv(dims, divs))
    {
      unsigned int tmp = surface(dims, divs);
      if ((min_surf>tmp) || ((min_surf==tmp) && compare_cuts(divs, _hyp_lat, dims)))
      {
	min_surf = tmp;
	_hyp_lat=divs;
	for(int i=0; i<4; ++i) _sub_lat[i] = lat_size[i] / _hyp_lat[i];
      }
    }
    int j = 0; while((j < len-1) && (digits[j] == 3)) digits[j++] = 0;
    ++digits[j];
  }

  if(min_surf == 2*lat_sz.compute_volume()) 
  {
    printf("Impossible to layout this lattice.\n");
    terminate_error(ERROR);
  }    
}

/*
// A simple test code

int main(int argc, char** argv)
{
int np = atoi(argv[1]);
int nx = atoi(argv[2]);
int ny = atoi(argv[3]);
int nz = atoi(argv[4]);
int nt = atoi(argv[5]);

lattice_desc lat(nx, ny, nz, nt), hyp(lat);

generate_min_surface(np, lat, hyp);

printf("To generate the minimal surface when dividing a lattice\n");
printf("of size %d %d %d %d in %d sublattices you need to break\n", 
nx, ny, nz, nt, np);
printf("it into %d %d %d %d\n", hyp.dims[0], hyp.dims[1], hyp.dims[2],
hyp.dims[3]);

return 0;
}

 */

} // qcd namespace
