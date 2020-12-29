#include <iostream>
#include <stdio.h>
#include "layout.h"
#include "su3.h"
#include "random.h"
#include "layout_minsurface.h"
#include "lattice.h"
#include "comm_low.h"

using namespace qcd;

int main(int argc,char* argv[]){
    init_machine(argc, argv);
    //lattice_desc* desc;
    position* position;
    std::cout << "Hello!" <<std::endl;
    int a=1;
    std::cout<<"a="<<a<<std::endl;

    layout_minsurface_eo desc(16, 16, 16, 128);
    
    std::cout<<"number of nodes="<<get_num_nodes()<<std::endl;
    shutdown_machine();

    return 0;
}
