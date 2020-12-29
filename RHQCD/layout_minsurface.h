// layout_minsurface.h
// Andrei Alexandru
// Jan 2010

#pragma once

#include "layout.h"
#include "comm_low.h"

namespace qcd {

struct surface_minimizer 
{
  dimensions lat_sz, _hyp_lat, _sub_lat;
  surface_minimizer(unsigned int np, dimensions lat_size, unsigned short sub_x=0,unsigned short sub_y=0,unsigned short sub_z=0);
};

struct layout_subdividing : lattice_desc
{
  dimensions sub_lat, hyp_lat;
  layout_subdividing(dimensions hyp_lat, dimensions sub_lat) :
    lattice_desc(hyp_lat*sub_lat, sub_lat.compute_volume()),
    sub_lat(sub_lat), hyp_lat(hyp_lat) { }
  unsigned int get_site_index(const position& p) const;
  rank_t get_site_rank(const position& p) const;
  position get_position(rank_t rnk, unsigned int idx) const;
  using lattice_desc::get_site_index;
  using lattice_desc::get_site_rank;
};

struct layout_subdividing_eo : lattice_desc
{
  dimensions sub_lat, hyp_lat;
  layout_subdividing_eo(dimensions hyp_lat, dimensions sub_lat) :
    lattice_desc(hyp_lat*sub_lat, sub_lat.compute_volume()),
    sub_lat(sub_lat), hyp_lat(hyp_lat) { }
  unsigned int get_site_index(const position& p) const;
  rank_t get_site_rank(const position& p) const;
  position get_position(rank_t rnk, unsigned int idx) const;
  using lattice_desc::get_site_index;
  using lattice_desc::get_site_rank;
};

struct layout_minsurface : private surface_minimizer, public layout_subdividing
{

  layout_minsurface(
      unsigned short nx, unsigned short ny, 
      unsigned short nz, unsigned short nt, unsigned int np=get_num_nodes(),
      unsigned short sub_x=0,unsigned short sub_y=0,unsigned short sub_z=0);
};

struct layout_minsurface_eo : private surface_minimizer, public layout_subdividing_eo
{

  layout_minsurface_eo(
      unsigned short nx, unsigned short ny, 
      unsigned short nz, unsigned short nt, unsigned int np=get_num_nodes(),
      unsigned short sub_x=0,unsigned short sub_y=0,unsigned short sub_z=0);

};

} // qcd namespace
