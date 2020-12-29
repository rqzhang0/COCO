// terminate.h
// Error messages and exit control for the code.
// Andrei Alexandru
// Oct 2003


#ifndef _TERMINATE_H
#define _TERMINATE_H

// Error messages
#define GWU_SUCCESS 0
#define ERROR   1
#define ERROR_READING_INITIAL_SET 2
#define ERROR_READING_EXTENDED_SET 3
#define ERROR_MEMORY_ALLOCATION 4
#define ERROR_OVERLAP_EIGENVALUE_PROBLEM 5
#define ERROR_IO_ACCESS 6
#define ERROR_UNKNOWN 7
#define ERROR_IN_ARGUMENTS 8

#define CHECK(a, b, c) check_func(a, b, c, __FILE__, __LINE__)

#ifdef __cplusplus
namespace qcd {
extern "C" {
#endif

int terminate_error(int iError,bool all_node=true);
int terminate_success();

void check_func(bool cond, const char* msg, int error_code, const char* filename, int line);

#ifdef __cplusplus
}
} // qcd namespace
#endif

#endif

