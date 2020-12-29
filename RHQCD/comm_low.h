#ifndef _COMM_COMMON_H
#define _COMM_COMMON_H

#include <unistd.h>
#include "layout.h"
#include "complex.h"

namespace qcd {

void init_machine(int &argc, char **&argv, bool single_node = false);
void shutdown_machine();
int get_num_nodes();
rank_t get_node_rank();
void synchronize();

void broadcast(char* buffer, size_t size);
void send(char* buffer, size_t size, rank_t dest, int tag=1);
void recv(char* buffer, size_t size, rank_t src, int tag=1);

struct comm_link
{
  comm_link();
  ~comm_link();
  void start_send(char* buffer, size_t size, rank_t dest, int tag=1);
  void start_recv(char* buffer, size_t size, rank_t dest, int tag=1);
  void finish_comm();
  comm_link(const comm_link&);
private:
  void* state;
};

template <typename T> void global_sum(T&);
template <typename R> void global_sum(complex<R>& v)
{
  global_sum(v.real);
  global_sum(v.imag);
}

template <typename T> void global_sum(T* val, int len)
{
  for(int i=0; i<len; ++i) global_sum(val[i]);
}

} // qcd namespace

#endif

