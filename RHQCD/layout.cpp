// layout.c
// Andrei Alexandru
// Jan 2006

#include "layout.h"
#include <stdexcept>

namespace qcd {


unsigned short& dimensions::operator[](unsigned int dim)
{
	if (dim==0)
		return nx;
	else if (dim==1)
		return ny;
	else if (dim==2)
		return nz;
	else if (dim==3)
		return nt;
	else
		throw std::logic_error("Invalid dimension index");
}

unsigned short dimensions::operator[](unsigned int dim) const
{
	if (dim==0)
		return nx;
	else if (dim==1)
		return ny;
	else if (dim==2)
		return nz;
	else if (dim==3)
		return nt;
	else
		throw std::logic_error("Invalid dimension index");
}

unsigned int dimensions::compute_volume() const
{
	return nx*ny*nz*nt;
}

dimensions dimensions::operator*(dimensions& b) const
{
	return dimensions(nx*b.nx, ny*b.ny, nz*b.nz, nt*b.nt);
}

position::position()
{
  pos[0] = pos[1] = pos[2] = pos[3] = 0;
}

position::position(short x, short y, short z, short t)
{
	pos[0] = x; pos[1] = y;
	pos[2] = z; pos[3] = t;
}

position::position(unsigned short* pos)
{
	for (int i=0; i<4; i++)
		this->pos[i] = pos[i];
}

position::position(unsigned int idx, const dimensions& dims)
{
  set_index(idx, dims);
}
void position::set_index(unsigned int idx, const dimensions& dims)
{
	assert(idx < dims.compute_volume());
	x() = idx % dims.nx;
	idx /= dims.nx;
	y() = idx % dims.ny;
	idx /= dims.ny;
	z() = idx % dims.nz;
	t() = idx / dims.nz;
}

short& position::operator[](unsigned int i)
{
	assert(i < 4);
	return pos[i];
}

const short& position::operator[](unsigned int i) const
{
	assert(i < 4);
	return pos[i];
}

bool position::operator==(position& p) const
{
	if (x() == p.x() &&
			y() == p.y() &&
			z() == p.z() &&
			t() == p.t())
		return true;
	else
		return false;
}

	bool position::operator!=(position& p) const
{
	return !(*this == p);
}

short& position::x() { return pos[0]; }
short& position::y() { return pos[1]; }
short& position::z() { return pos[2]; }
short& position::t() { return pos[3]; }

const short& position::x() const { return pos[0]; }
const short& position::y() const { return pos[1]; }
const short& position::z() const { return pos[2]; }
const short& position::t() const { return pos[3]; }

position position::operator+ (const position& p) const
{
	position res = *this;
	for(int i=0; i<4; ++i) 
		res[i] += p[i];
	return res; 
}

position position::operator- (const position& p) const
{
	position res = *this;
	for(int i=0; i<4; ++i) 
		res[i] -= p[i];
	return res; 
}

void position::normalize(const dimensions& dims)
{
	for(int i=0; i<4; ++i) normalize_dir(i, dims[i]);
}

void position::normalize_dir(int dir, int dim)
{
	assert(dir < 4);
	pos[dir] = pos[dir]%dim;
	while(pos[dir] < 0) pos[dir] += dim;
}

enum parity position::parity() const
{
	return (enum parity) ((x()+y()+z()+t()) % 2);
}

unsigned int position::index(const dimensions& dims) const
{
	return x() + dims.nx*(y() + dims.ny*(z() + dims.nz*t()));
}

position operator-(const position& p) { return position()-p; }

/* 
 * TODO: These shouldn't exist and should be pure virtual
 */
rank_t lattice_desc::get_site_rank(const position &pos) const {
	throw new std::runtime_error("Lattice_desc doesn't implement layout");
}
unsigned int lattice_desc::get_site_index(const position &pos) const {
	throw new std::runtime_error("Lattice_desc doesn't implement layout");
}
position lattice_desc::get_position(local_idx_t index, rank_t rank) const {
	throw new std::runtime_error("Lattice_desc doesn't implement layout");
}


int unsafe_index(int x, int y, int z, int t, int* dims)
{
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];
 
  return x + nx*(y + ny*(z + nz*t));
}

int get_index(int x, int y, int z, int t, int* dims)
{
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];
  int nt = dims[3];

  x = (x+nx)%nx;
  y = (y+ny)%ny;
  z = (z+nz)%nz;
  t = (t+nt)%nt;
  return x + nx*(y + ny*(z + nz*t));
}

void get_pos(int index, int *p, int* dims)
{
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];

  p[0] = index % nx;
  index /= nx;
  p[1] = index % ny;
  index /= ny;
  p[2] = index % nz;
  p[3] = index/nz;
}

/*
 * Even/odd layout
 */
int get_pos_parity(int x, int y, int z, int t)
{
  return (x+y+z+t)%2;
}

int get_index_parity(int i, const lattice_desc* desc)
{
  int temp = desc->nx;
  int parity = i;

  parity += i/temp;
  temp *= desc->ny;
  parity += i/temp;
  temp *= desc->nz;
  parity += i/temp;

  return parity%2;
}

int unsafe_index_eo(int x, int y, int z, int t, int* dims)
{
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];

  int i = x + nx*(y + ny*(z + nz*t));
  return i/2;
}

int get_index_eo(int x, int y, int z, int t, int* dims)
{
  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];
  int nt = dims[3];

  x = (x+nx)%nx;
  y = (y+ny)%ny;
  z = (z+nz)%nz;
  t = (t+nt)%nt;
  return unsafe_index_eo(x, y, z, t, dims);
}

void get_pos_eo(int i, int parity, int *p, int* dims)
{
  i *= 2;

  p[0] = i % dims[0];
  i /= dims[0];
  p[1] = i % dims[1];
  i /= dims[1];
  p[2] = i % dims[2];
  i /= dims[2];
  p[3] = i % dims[3];
  p[0] += (p[1]+p[2]+p[3]+parity) % 2;
}


} // qcd namespace
