// layout.h
// Ben Gamari
// August 2009

#ifndef _LAYOUT_H
#define _LAYOUT_H

#include <assert.h>

namespace qcd {


typedef unsigned int rank_t;
typedef unsigned int local_idx_t;
typedef unsigned int global_idx_t;

enum parity { EVEN=0, ODD=1 };
enum direction { X=0, Y=1, Z=2, T=3 };

struct position;

struct dimensions {
	unsigned short nx, ny, nz, nt;

	dimensions() {}
	dimensions(unsigned short nx, unsigned short ny, unsigned short nz, unsigned short nt) :
		nx(nx), ny(ny), nz(nz), nt(nt) { }

	unsigned int compute_volume() const;
	unsigned short& operator[](unsigned int dim);
	unsigned short  operator[](unsigned int dim) const;
	dimensions operator*(dimensions& b) const;
};

struct position {
	short pos[4];

	position();


	position(short x, short y, short z, short t);
	position(unsigned short* pos);
	position(unsigned int idx, const dimensions& dims);
  void set_index(unsigned int idx, const dimensions& dims);

	short& operator[](unsigned int i);
	const short& operator[](unsigned int i) const;
	bool operator==(position& p) const;
	bool operator!=(position& p) const;

	short& x();
	short& y();
	short& z();
	short& t();

	const short& x() const;
	const short& y() const;
	const short& z() const;
	const short& t() const;

	position operator+ (const position& p) const;
	position operator- (const position& p) const;
	void normalize(const dimensions& dims);
	void normalize_dir(int dir, int dim);
	enum parity parity() const;
	unsigned int index(const dimensions& dims) const;
};

// Unary minus for position
position operator-(const position&);

/* lattice_desc: Lattice description
 * This describes the shape of a lattice. The volume is also included to avoid
 * recomputing.
 */
struct lattice_desc : dimensions {
	const unsigned int volume;
	const unsigned int sites_on_node;

	lattice_desc(dimensions dims) :
		dimensions(dims), volume(dims.compute_volume()), sites_on_node(dims.compute_volume()) { }

	lattice_desc(dimensions dims, unsigned int sites_on_node) :
		dimensions(dims), volume(dims.compute_volume()), sites_on_node(sites_on_node) { }

	lattice_desc(unsigned short nx, unsigned short ny, unsigned short nz, unsigned short nt) :
		dimensions(nx,ny,nz,nt), volume(nx*ny*nz*nt), sites_on_node(volume) { }

	virtual ~lattice_desc() {}

	// TODO: These should be pure virtual
	virtual rank_t get_site_rank(const position &pos) const;
	virtual unsigned int get_site_index(const position &pos) const ;
	virtual position get_position(local_idx_t index, rank_t rank) const;

	virtual local_idx_t get_site_index(global_idx_t idx) const
	{
	  position p(idx, *this);
	  return get_site_index(p);
	}


	virtual rank_t get_site_rank(global_idx_t idx) const 
	{
	  position p(idx, *this);
	  return get_site_rank(p);
	}

        virtual global_idx_t get_global_index(local_idx_t index, rank_t rank) const
        {
          position p = get_position(index, rank);
          return p.index(*this);
        }

};

/*
 * Old layout functions (TODO: Deprecate these)
 */
int unsafe_index(int x, int y, int z, int t, int* dims);
int get_index(int x, int y, int z, int t, int* dims);
void get_pos(int index, int *p, int* dims);

/*
 * Even/odd layout
 */
int get_pos_parity(int x, int y, int z, int t);
int get_index_parity(int i, const lattice_desc* desc);
int unsafe_index_eo(int x, int y, int z, int t, int* dims);
int get_index_eo(int x, int y, int z, int t, int* dims);
void get_pos_eo(int index, int parity, int *p, int* dims);

} // qcd namespace

#endif


