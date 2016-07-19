#include <iostream>

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "ChemDriver.H"

typedef ChemDriver::Edge Edge;
typedef std::list<Edge> EdgeList;

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    ParmParse pp;

    ChemDriver cd;
    std::string QPDatom="H"; pp.query("QPDatom",QPDatom);
    EdgeList edges = cd.getEdges(QPDatom);

    for (EdgeList::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
      std::cout << *it << std::endl;
    }

    BoxLib::Finalize();
    return 0;
}
