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
    std::string QPDatom="C"; pp.query("QPDatom",QPDatom);

    int PrintVerbose = 0; pp.query("PrintVerbose",PrintVerbose);
    int HackSplitting = 0; pp.query("HackSplitting",HackSplitting);
    EdgeList edges = cd.getEdges(QPDatom,PrintVerbose,HackSplitting);

    for (EdgeList::const_iterator it = edges.begin(); it!=edges.end(); ++it) {
      //std::cout << *it << std::endl;
    }

    BoxLib::Finalize();
    return 0;
}
