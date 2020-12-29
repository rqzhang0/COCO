#include "layout.h"
#include "lattice.h"
#include "comm_low.h"
#include <stdio.h>
#include <iostream>
#include "layout_minsurface.h"

using namespace qcd;

int main(int argc,char* argv[]){
    init_machine(argc, argv);
    //lattice_desc* desc;
    position* position;
    std::cout << "Hello!" <<std::endl;
    int a=1;
    std::cout<<"a="<<a<<std::endl;
    //layout_minsurface_eo desc(16, 16, 16, 128);
    return 0;
}
