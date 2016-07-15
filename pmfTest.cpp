#include <fstream>
#include <sstream>

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "Utility.H"
#include "ChemDriver.H"

void ToFab(const std::vector<std::vector<Real> >& data_vals,
	   FArrayBox& fab);

void read_file(const std::string&               infile,
	       std::vector<std::vector<Real> >& data_vals);
int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    ParmParse pp;
    std::string infile = "ALZ_F_pmf.dat";
    pp.query("infile",infile);

    ChemDriver cd;

    std::vector<std::vector<Real> > data_vals;
    read_file(infile,data_vals);

    FArrayBox fab;
    ToFab(data_vals,fab);

    BoxLib::Finalize();
    return 0;
}

void ToFab(const std::vector<std::vector<Real> >& data_vals,
		 FArrayBox& fab) 
{
  int nVals = data_vals.size();
  if (nVals == 0) return;
  int nComp = data_vals[0].size();
  if (nComp == 0) return;
  Box box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(nVals-1,0,0)));
  fab.resize(box,nComp);
  for (int i=0; i<nVals; ++i) {
    IntVect iv(D_DECL(i,0,0));
    for (int n=0; n<nComp; ++n) {
      fab(iv,n) = data_vals[i][n];
    }
  }
}

void read_file(const std::string&               infile,
	       std::vector<std::vector<Real> >& data_vals)
{
    std::ifstream ifs(infile.c_str());

    std::string line;
    std::getline(ifs,line);
    std::vector<std::string> variable_names = BoxLib::Tokenize(line,"=\", ");
    variable_names.erase(variable_names.begin());

    std::getline(ifs,line);
    std::vector<std::string> tokens = BoxLib::Tokenize(line,"=\", ");
    tokens.erase(tokens.begin());

    int nPts = 0;
    for (int i=0; i<tokens.size()/2; ++i) {
      if (tokens[2*i] == "I") {
	std::istringstream iss(tokens[2*i+1]);
	if (!(iss >> nPts)) { break; } // error
      }
    }
    if (nPts == 0) {
      BoxLib::Abort("No data in file");
    }

    int nComp = variable_names.size();
    int num_lines = 0;
    while (std::getline(ifs, line))
    {
      data_vals.push_back(std::vector<Real>(nComp));
      std::istringstream iss(line);
      for (int i=0; i<nComp; ++i) {
	if (!(iss >> data_vals[num_lines][i])) { break; } // error
      }
      num_lines++;
    }
    
    if (nPts != num_lines) {
      BoxLib::Abort("Trouble reading all of the records");
    }
}
