// terminate.c
// Terminates execution gracefully
// Copied from com_mpi.c

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "terminate.h"
#include "comm_low.h"

namespace qcd {


int terminate_success()
{
  time_t time_stamp;

  if(get_node_rank() == 0)
  {
    time(&time_stamp);
    printf("exit: %s\n", ctime(&time_stamp));
  }
  return GWU_SUCCESS;
}

int terminate_error(int iError,bool all_node)
{
  time_t time_stamp;

  int tn = get_node_rank();
  time(&time_stamp);
 if(get_node_rank()==0||all_node==true)
 {
  printf("node%04d: termination: %s\n", tn, ctime(&time_stamp));
  printf("node%04d: Termination: status = %d\n", tn, iError);
 }
  fflush(stdout);
  exit(iError);
}

void check_func(bool cond, const char* msg, int error_code, const char* filename, int line)
{ 
  if(!cond)
  { 
    if(get_node_rank() == 0)
    { 
      printf("CHECK failed in file %s at line %d (error code %d):\n", filename, line, error_code);
      printf("MSG:  %s\n", msg);
    }
    terminate_error(error_code);
  }
}

} // qcd namespace
