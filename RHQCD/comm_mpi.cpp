#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "comm_low.h"
//#include "timer.h"

namespace qcd {

#ifdef CUDA
void initialize_cuda_devices();
#endif

template<typename T> MPI_Datatype get_mpi_type() { return 0; }
template<> MPI_Datatype get_mpi_type<char>() { return MPI_CHAR; }
template<> MPI_Datatype get_mpi_type<int>() { return MPI_INT; }
template<> MPI_Datatype get_mpi_type<short>() { return MPI_SHORT; }
template<> MPI_Datatype get_mpi_type<float>() { return MPI_FLOAT; }
template<> MPI_Datatype get_mpi_type<double>() { return MPI_DOUBLE; }

static void mpi_error(MPI_Comm *comm, int *stat, ...)
{
	int len;
	char err_string[MPI_MAX_ERROR_STRING];

	fprintf(stderr, "MPI error: %d\n", *stat);
	MPI_Error_string(*stat, err_string, &len);
	fprintf(stderr, "%s\n", err_string);
	abort();
}

void init_machine(int &argc, char **&argv, bool single_node)
{
  int flag;
  MPI_Comm comm;
  MPI_Errhandler errhandler;


  flag = MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  if (flag) mpi_error(&comm, &flag);
  //flag = MPI_Errhandler_create(mpi_error, &errhandler);
  flag = MPI_Comm_create_errhandler(mpi_error, &errhandler);
  if (flag) mpi_error(&comm, &flag);
  //flag = MPI_Errhandler_set(MPI_COMM_WORLD, errhandler);
  flag = MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
  if (flag) mpi_error(&comm, &flag);

#ifdef CUDA
  initialize_cuda_devices();
#endif

  if(get_node_rank() == 0) fprintf(stderr, "MPI comm started with %d process.\n", get_num_nodes());
  //time_stamp("Start");
  
  if(single_node) //check is single node
  {
    if(get_num_nodes() !=1)
    {
      if(get_node_rank() == 0)
      fprintf(stderr, "This code is only ment to run serially ... exiting\n");
      shutdown_machine();
      abort();
    }
  }
  //global_timer.start("MPI");
}

void shutdown_machine()
{
  //time_stamp("Exit");
  //global_timer.stop();
  //synchronize();
  MPI_Finalize();
}

int get_num_nodes()
{
	int n;
	MPI_Comm_size(MPI_COMM_WORLD, &n);
	return n;
}

rank_t get_node_rank()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

void synchronize()
{
	MPI_Barrier(MPI_COMM_WORLD);
}

void broadcast(char* buffer, size_t size)
{
	MPI_Bcast(buffer, size, get_mpi_type<char>(), 0, MPI_COMM_WORLD);
}


void send(char* buffer, size_t size, rank_t dest, int tag)
{
	MPI_Send(buffer, size, get_mpi_type<char>(), dest, tag, MPI_COMM_WORLD);
}


void recv(char* buffer, size_t size, rank_t src, int tag)
{
	MPI_Status status;
	MPI_Recv(buffer, size, get_mpi_type<char>(), src, tag, MPI_COMM_WORLD, &status);
}


template<typename T>
void global_sum(T &value)
{
	T work;
	MPI_Allreduce(&value, &work, 1, get_mpi_type<T>(), MPI_SUM, MPI_COMM_WORLD);
	value = work;
}

template void global_sum(double& val); 
template void global_sum(float&  val); 
template void global_sum(int&    val); 

comm_link::comm_link()
{
  state = (void*) new MPI_Request;
}

comm_link::~comm_link()
{
  delete (MPI_Request*)state; 
}

comm_link::comm_link(const comm_link& l)
{
  MPI_Request* tmp = new MPI_Request;
  *tmp = *(MPI_Request*) l.state;
  state = (void*) tmp;
}

void comm_link::start_send(char* buffer, size_t size, rank_t dest, int tag)
{
  MPI_Issend(buffer, size, MPI_CHAR, dest, tag, MPI_COMM_WORLD,
      (MPI_Request*)state);
}

void comm_link::start_recv(char* buffer, size_t size, rank_t dest, int tag)
{
  MPI_Irecv(buffer, size, MPI_CHAR, dest, tag, MPI_COMM_WORLD,
      (MPI_Request*)state);
}

void comm_link::finish_comm()
{
  MPI_Wait((MPI_Request*)state, MPI_STATUS_IGNORE);
  assert(MPI_REQUEST_NULL == *((MPI_Request*)state));
}

} // qcd namespace
