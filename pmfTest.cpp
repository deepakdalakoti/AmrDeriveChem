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

void write_tv(const std::vector<Real>& vals,
	      const std::string&       file_name,
	      const ChemDriver&        cd);

void read_tv(std::vector<Real>& vals,
	     const ChemDriver&  cd);

void read_tv_fab_file(Array<Real>& typical_values);

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);
    
    ParmParse pp;
    std::string infile = "ALZ_F_pmf.dat";
    pp.query("infile",infile);

    ChemDriver cd;
    int Nspec = cd.numSpecies();

    std::vector<std::vector<Real> > data_vals;
    read_file(infile,data_vals);

    FArrayBox fab;
    ToFab(data_vals,fab);
    const Box& box = fab.box();
    //Box box(IntVect(D_DECL(0,0,0)),
//	    IntVect(D_DECL(1,0,0)));

    Real Patm = 1.0;
    pp.query("Patm",Patm);
    FArrayBox& T = fab;
    int sCompT = 2;
    FArrayBox& X = fab;
    int sCompX = 3;
    FArrayBox Y(box,Nspec);
    int sCompY = 0;
    FArrayBox Rho(box,1);
    int sCompR = 0;
    FArrayBox RhoH(box,1);
    int sCompRH = 0;
    
    cd.moleFracToMassFrac(Y,X,box,sCompX,sCompY);
    cd.getRhoGivenPTY(Rho,Patm,T,Y,box,sCompT,sCompY,sCompR);
    cd.getHmixGivenTY(RhoH,T,Y,box,sCompT,sCompY,sCompRH);
    RhoH.mult(Rho,sCompR,sCompRH,1);

    FArrayBox RhoY(box,Nspec);
    int sCompRY = 0;
    RhoY.copy(Y,sCompY,sCompRY,Nspec);
    for (int i=0; i<Nspec; ++i) {
      RhoY.mult(Rho,sCompR,sCompRY+i,1);
    }

    FArrayBox Src(box,Nspec+1);
    Src.setVal(0);

    int first_spec = 3;
    int nComp = first_spec + Nspec + 4;
    Array<Real> typical_values(nComp,1);
    if (pp.countval("tv_fab_file") > 0) {
      read_tv_fab_file(typical_values);
      write_tv(typical_values,"TYPICAL_VALUES.txt",cd);
    }
    read_tv(typical_values,cd);


    Real dt = 1.e-12;
    pp.query("dt",dt);

    int nchemdiag = 0; pp.query("nchemdiag",nchemdiag);
    int nReac = cd.numReactions();
    FArrayBox chemDiag(box,nReac); chemDiag.setVal(0);
    FArrayBox Snew(box,fab.nComp());
    FArrayBox& RhoYnew = Snew;
    FArrayBox& RhoHnew = Snew;
    FArrayBox& Tnew = Snew;
    FArrayBox FC(box,1); FC.setVal(0);
    cd.solveTransient_sdc(RhoYnew,RhoHnew,Tnew,RhoY,RhoH,T,Src,FC,
			  box,sCompRY,sCompRH,sCompT,dt,Patm,typical_values,&chemDiag,nchemdiag);

    for (int i=0; i<nReac; ++i) {
      //std::cout << i << " " << chemDiag.min(i) << " " << chemDiag.max(i) << std::endl;
    }
    BoxLib::Finalize();
    return 0;
}

void read_tv_fab_file(Array<Real>& typical_values)
{
  std::string tvfile;
  ParmParse pp;
  pp.get("tv_fab_file",tvfile);
  std::ifstream tvis;
  tvis.open(tvfile.c_str(),std::ios::in|std::ios::binary);
  if (tvis.good())
  {
    FArrayBox tvfab;
    tvfab.readFrom(tvis);
    if (tvfab.nComp() != typical_values.size())
      BoxLib::Abort("Typical values file has wrong number of components");
    for (int i=0; i<typical_values.size(); ++i) {
      typical_values[i] = tvfab.dataPtr()[i];
    }
  }
}

void write_tv(const std::vector<Real>& vals,
	      const std::string&       file_name,
	      const ChemDriver&        cd)
{
  ParmParse pp;
  int Nspec = cd.numSpecies();
  std::ofstream ofs(file_name.c_str());
  ofs << "typVal_Temp = " << vals[Nspec+5] << '\n';
  ofs << "typVal_RhoH = " << vals[Nspec+3] << '\n';
  ofs << "typVal_Density = " << vals[2] << '\n';
  for (int i=0; i<Nspec; ++i) {
    ofs << "typValY_" << cd.speciesNames()[i] << " = " << vals[i+3] << '\n';
  }
}

void read_tv(std::vector<Real>& vals,
	     const ChemDriver&  cd)
{
  int Nspec = cd.numSpecies();
  const Array<std::string>& speciesNames = cd.speciesNames();
  ParmParse pp;
  int first_spec = 3;
  for (int i=0; i<Nspec; ++i) {
    const std::string ppStr = std::string("typValY_") + speciesNames[i];
    if (pp.countval(ppStr.c_str())>0) {
      pp.get(ppStr.c_str(),vals[first_spec + i]);
    }
    std::string otherKeys[5] = {"Temp", "RhoH", "Vel", "Trac", "Density"};
    int otherIdx[5] = {Nspec+5, Nspec+3, 0, Nspec+4, 2};
    for (int i=0; i<5; ++i) {
      const std::string ppStr(std::string("typVal_")+otherKeys[i]);
      if (pp.countval(ppStr.c_str())>0) {
        pp.get(ppStr.c_str(),vals[otherIdx[i]]);
      }
    }
  }
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
