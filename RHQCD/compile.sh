#mpicxx -o test test.cpp
#mpicxx -o test test.cpp comm_mpi.cpp packing.cpp terminate.cpp layout.cpp su3.cpp randomSU3.cpp memalign.cpp  vector_util.cpp
mpicxx -o test test.cpp comm_mpi.cpp packing.cpp terminate.cpp layout.cpp su3.cpp memalign.cpp vector_util.cpp layout_minsurface.cpp
