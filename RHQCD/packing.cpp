#include "cpu_vector.h"
#include "layout.h"

namespace qcd {

/*
 * Standard packing
 */
void pack_links(su3_field& lnk, double* packed_links)
{
  int volume = lnk.desc->volume;
  for(int dir=0; dir<4; ++dir)
  for(int r=0; r<3; ++r)
  for(int c=0; c<3; ++c)
  for(int i=0; i<volume; ++i)
  {
    packed_links[dir*3*3*2*volume + r*3*2*volume + c*2*volume + 0*volume + i] = lnk.data[i+volume*dir].e[r][c].real;
    packed_links[dir*3*3*2*volume + r*3*2*volume + c*2*volume + 1*volume + i] = lnk.data[i+volume*dir].e[r][c].imag;
  }
}

void pack_spinor(vector& vec, double* packed_vec)
{
  int volume = vec.desc->volume;
  for(int c=0; c<3; ++c)
  for(int s=0; s<4; ++s) 
  for(int i=0; i<volume; ++i)
  {
    packed_vec[c*4*2*volume+s*2*volume+0*volume+i] = vec.data[i].d[s].c[c].real;
    packed_vec[c*4*2*volume+s*2*volume+1*volume+i] = vec.data[i].d[s].c[c].imag;
  }
}

void unpack_spinor(double* packed_vec, vector& vec)
{
  int volume = vec.desc->volume;
  for(int c=0; c<3; ++c)
  for(int s=0; s<4; ++s) 
  for(int i=0; i<volume; ++i)
  {
    vec.data[i].d[s].c[c] = double_complex(
		    packed_vec[c*4*2*volume+s*2*volume+0*volume+i],
		    packed_vec[c*4*2*volume+s*2*volume+1*volume+i]);
  }
}


/*
 * Even/odd packing
 */
void pack_links_eo(su3_field& lnk, double* packed_links, enum parity parity)
{
  int volume_2 = lnk.desc->volume/2;
  double* start = packed_links;

  for (int i=0; i<lnk.desc->volume; ++i)
  {
    if (get_index_parity(i, lnk.desc) == parity) {
      double *off = start;

      for(int dir=0; dir<4; ++dir)      // direction
      for(int r=0; r<3; ++r)            // row
      for(int c=0; c<3; ++c)            // col
      {
	*off = lnk.data[i + lnk.desc->volume*dir].e[r][c].real;
	off += volume_2;
	*off = lnk.data[i + lnk.desc->volume*dir].e[r][c].imag;
	off += volume_2;
      }

      start++;
    }
  }
}

void pack_spinor_eo(vector& vec, double* packed_vec, enum parity parity)
{
  int volume_2 = vec.desc->volume/2; 
  double* start = packed_vec;

  for (int i=0; i<vec.desc->volume; ++i)
  {
    if (get_index_parity(i, vec.desc) == parity) {
      double *off = start;

      for(int c=0; c<3; ++c)		// color
      for(int s=0; s<4; ++s)		// spinor
      {
	*off = vec.data[i].d[s].c[c].real;
	off += volume_2;
	*off = vec.data[i].d[s].c[c].imag;
	off += volume_2;
      }

      start++;
    }
  }
}

void unpack_spinor_eo(double* packed_vec, vector& vec, enum parity parity)
{
  int volume_2 = vec.desc->volume/2; 
  double* start = packed_vec;

  for (int i=0; i<vec.desc->volume; ++i)
  {
    if (get_index_parity(i, vec.desc) == parity) {
      double *off = start;

      for(int c=0; c<3; ++c)            // color
      for(int s=0; s<4; ++s)            // spinor
      {
	double real, imag;
	real = *off;
	off += volume_2;
	imag = *off;
	off += volume_2;

	vec.data[i].d[s].c[c] = double_complex(real, imag);
      }

      start++;
    }
  }
}

} // qcd namespace

