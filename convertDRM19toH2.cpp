
#include <iostream>
#include <fstream>
#include <cstdlib>
#ifndef WIN32
#include <unistd.h>
#endif
#include <map>

using std::cout;
using std::cerr;
using std::endl;
using std::map;

#include "REAL.H"
#include "FArrayBox.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"
#include "ChemDriver_F.H"

#include "convertDRM19toH2_F.H"
#include "ckfuncs.H" /* include to make Fortran-looking ChemKin-based functions visible */

static
void
PrintUsage (char* progName)
{
    std::cout << "\nUsage:\n"
              << progName
              << "\n\tifile = ifile"
              << "\n\tofile = ofile"
              << "\n\t[-help]"
              << "\n\n";
    exit(1);
}

#define MEIO ParallelDescriptor::IOProcessor()

void
writePlotfile(const PArray<MultiFab>&    data,
              Real                       time,
              const Array<Real>&         probLo,
              const Array<Real>&         probHi,
              const Array<int>&          refRatio,
              const Array<Box>&          probDomain,
              const Array<Array<Real> >& dxLevel,
              int                        coordSys,
              std::string&               oFile,
              const Array<std::string>&  names,
              bool                       verbose)
{
    // This is the version of plotfile that will be written
    std::string plotFileVersion = "NavierStokes-V1.1";

    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(oFile,0755))
            BoxLib::CreateDirectoryFailed(oFile);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    
    std::string oFileHeader(oFile);
    oFileHeader += "/Header";
    
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ofstream os;
    
    //os.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Opening file = " << oFileHeader << '\n';
    
    os.open(oFileHeader.c_str(), std::ios::out|std::ios::binary);
    
    if (os.fail())
        BoxLib::FileOpenFailed(oFileHeader);
    //
    // Start writing plotfile.
    //
    os << plotFileVersion << '\n';
    int n_var = data[0].nComp();
    os << n_var << '\n';
    for (int n = 0; n < n_var; n++) os << names[n] << '\n';
    os << BL_SPACEDIM << '\n';
    os << time << '\n';
    const int finestLevel = data.size() - 1;
    os << finestLevel << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probLo[i] << ' ';
    os << '\n';
    for (int i = 0; i < BL_SPACEDIM; i++) os << probHi[i] << ' ';
    os << '\n';
    for (int i = 0; i < finestLevel; i++) os << refRatio[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) os << probDomain[i] << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++) os << 0 << ' ';
    os << '\n';
    for (int i = 0; i <= finestLevel; i++)
    {
        for (int k = 0; k < BL_SPACEDIM; k++)
            os << dxLevel[i][k] << ' ';
        os << '\n';
    }
    os << coordSys << '\n';
    os << "0\n"; // The bndry data width.
    //
    // Write out level by level.
    //
    for (int iLevel = 0; iLevel <= finestLevel; ++iLevel)
    {
        //
        // Write state data.
        //
        const BoxArray& ba = data[iLevel].boxArray();
        int nGrids = ba.size();
        char buf[64];
        sprintf(buf, "Level_%d", iLevel);
        
        if (ParallelDescriptor::IOProcessor())
        {
            os << iLevel << ' ' << nGrids << ' ' << time << '\n';
            os << 0 << '\n';
            
            for (int i = 0; i < nGrids; ++i)
            {
                const Box& b = ba[i];
                for (int n = 0; n < BL_SPACEDIM; n++)
                {
                    Real glo = b.smallEnd()[n]*dxLevel[iLevel][n];
                    Real ghi = (b.bigEnd()[n]+1)*dxLevel[iLevel][n];
                    os << glo << ' ' << ghi << '\n';
                }
            }
            //
            // Build the directory to hold the MultiFabs at this level.
            //
            std::string Level(oFile);
            Level += '/';
            Level += buf;
            
            if (!BoxLib::UtilCreateDirectory(Level, 0755))
                BoxLib::CreateDirectoryFailed(Level);
        }
        //
        // Force other processors to wait till directory is built.
        //
        ParallelDescriptor::Barrier();
        //
        // Now build the full relative pathname of the MultiFab.
        //
        static const std::string MultiFabBaseName("/MultiFab");
        
        std::string PathName(oFile);
        PathName += '/';
        PathName += buf;
        PathName += MultiFabBaseName;
        
        if (ParallelDescriptor::IOProcessor())
        {
            //
            // The full name relative to the Header file.
            //
            std::string RelativePathName(buf);
            RelativePathName += '/';
            RelativePathName += MultiFabBaseName;
            os << RelativePathName << '\n';
        }
        VisMF::Write(data[iLevel], PathName, VisMF::OneFilePerCPU);
    }
    
    os.close();
}

// Set of ChemKin funcs I need to see for the "old" mech:
//   CKINDX, CKXTY, fesymnum_, fesymname_, speciesEnthalpy, CKHBMS, cp_R, CKCPBS
extern "C" {

    void CKHBMS(double * T, double * y, int * iwrk, double * rwrk, double * hbms);
};


// Set of ChemKin funcs I need for the "new" mech
void CKINDX_new(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 27;
    *nfit = -1; /*Why do you need this anyway ?  */
}

void CKXTY_new(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    y[0] = x[0]*2.015940/XW; 
    y[1] = x[1]*1.007970/XW; 
    y[2] = x[2]*15.999400/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*18.015340/XW; 
    y[6] = x[6]*33.006770/XW; 
    y[7] = x[7]*34.014740/XW; 
    y[8] = x[8]*28.013400/XW; 

    return;
}
int fesymnum_new_(const char* s1)
{
    if (strcmp(s1, "H2")==0) return 0; 
    if (strcmp(s1, "H")==0) return 1; 
    if (strcmp(s1, "O")==0) return 2; 
    if (strcmp(s1, "O2")==0) return 3; 
    if (strcmp(s1, "OH")==0) return 4; 
    if (strcmp(s1, "H2O")==0) return 5; 
    if (strcmp(s1, "HO2")==0) return 6; 
    if (strcmp(s1, "H2O2")==0) return 7; 
    if (strcmp(s1, "N2")==0) return 8; 
    /*species name not found */
    return -1;
}
char* fesymname_new_(int sn)
{
    if (sn==0) return "H2"; 
    if (sn==1) return "H"; 
    if (sn==2) return "O"; 
    if (sn==3) return "O2"; 
    if (sn==4) return "OH"; 
    if (sn==5) return "H2O"; 
    if (sn==6) return "HO2"; 
    if (sn==7) return "H2O2"; 
    if (sn==8) return "N2"; 
    /*species name not found */
    return "NOTFOUND";
}
/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
void speciesEnthalpy_new(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 / tc[1];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 / tc[1];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 / tc[1];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 / tc[1];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 / tc[1];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 / tc[1];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 / tc[1];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 / tc[1];
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 / tc[1];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
    }
    return;
}
void CKHBMS_new(double *T, double *y, int * iwrk, double * rwrk, double * hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy_new(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/2.015940; /*H2 */
    result += y[1]*hml[1]/1.007970; /*H */
    result += y[2]*hml[2]/15.999400; /*O */
    result += y[3]*hml[3]/31.998800; /*O2 */
    result += y[4]*hml[4]/17.007370; /*OH */
    result += y[5]*hml[5]/18.015340; /*H2O */
    result += y[6]*hml[6]/33.006770; /*HO2 */
    result += y[7]*hml[7]/34.014740; /*H2O2 */
    result += y[8]*hml[8]/28.013400; /*N2 */

    *hbms = result * RT;
}

void cp_R_new(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}

void CKCPBS_new(double *T, double *y, int * iwrk, double * rwrk, double * cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9]; /* temporary storage */
    cp_R_new(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/2.01594; /*H2 */
    result += cpor[1]*y[1]/1.00797; /*H */
    result += cpor[2]*y[2]/15.9994; /*O */
    result += cpor[3]*y[3]/31.9988; /*O2 */
    result += cpor[4]*y[4]/17.0074; /*OH */
    result += cpor[5]*y[5]/18.0153; /*H2O */
    result += cpor[6]*y[6]/33.0068; /*HO2 */
    result += cpor[7]*y[7]/34.0147; /*H2O2 */
    result += cpor[8]*y[8]/28.0134; /*N2 */

    *cpbs = result * 8.314e+07;
}

static
void
ConvertDRM19toH2 (const FArrayBox&   T_old,
                  int                TComp_old,
                  const FArrayBox&   X_old,
                  int                xComp_old,
                  FArrayBox&         Y_new,
                  int                yComp_new,
                  int                nSpec_new,
                  FArrayBox&         T_new,
                  int                TComp_new,
                  ChemDriver&        cd,
                  string&            atom,
                  const Array<Real>& Xold1,
                  Real               Told1,
                  const Array<Real>& Xold2)
{
    BL_ASSERT(T_old.box() == X_old.box());
    BL_ASSERT(T_old.box() == Y_new.box());
    BL_ASSERT(T_old.box() == T_new.box());
    Real rwrk;
    int iwrk;

    const int nSpec_old = cd.numSpecies();
    const Array<string>& spNames = cd.speciesNames();
    Array<int> nu(nSpec_old);
    Real W2=0, f2=0;
    // const Real M = cd.elementAtomicWt()[cd.indexElt(atom)]; // function not yet implemented
    const Real M = 12.01115;
    for (int i=0; i<nSpec_old; ++i)
    {
        nu[i] = cd.numberOfElementXinSpeciesY(atom,spNames[i]);
        W2 += Xold2[i]*cd.speciesMolecWt()[i];
        f2 += M * nu[i] * Xold2[i];
    }
    f2 *= 1./W2;

    int idO2n = fesymnum_new_("O2");
    int idN2n = fesymnum_new_("N2");
    int idO2o = cd.index("O2");
    int idN2o = cd.index("N2");

    Array<Real> Yold1 = cd.moleFracToMassFrac(Xold1);
    Array<Real> Yold2 = cd.moleFracToMassFrac(Xold2);

    // Set up stream1 in the new species (here just gonna assume N2,O2 are the only relevant species)
    Array<Real> Xnew1(nSpec_new,0.), Ynew1(nSpec_new), Ynew2(nSpec_new);
    Xnew1[idO2n] = Xold1[idO2o];
    Xnew1[idN2n] = Xold1[idN2o];
    CKXTY_new(Xnew1.dataPtr(),&iwrk,&rwrk,Ynew1.dataPtr());

    Real errMAX=1.e-8;
    int NiterMAX = 20;
    Array<Real> res(NiterMAX);

    Array<Real> YtmpO(nSpec_old), XtmpO(nSpec_old);
    Array<Real> YtmpN(nSpec_new), XtmpN(nSpec_new);
    Array<Real> Xnew2(nSpec_new+1); // Tnew is going to be passed in the last slot
    Array<Real> Ynew(nSpec_new);
    Real Tnew2;
    for (IntVect iv = T_old.box().smallEnd(); iv <= T_old.box().bigEnd(); T_old.box().next(iv))
    {
        Real Ttmp = T_old(iv,TComp_old);

        Real f=0, Wold=0;
        for (int i=0; i<nSpec_old; ++i)
        {
            XtmpO[i] = X_old(iv,xComp_old+i);
            Wold += XtmpO[i] * cd.speciesMolecWt()[i];
            f += M * nu[i] * XtmpO[i];
        }
        f *= 1./Wold;
        Real alpha1 = std::min(1.0, std::max(0.0, ( f2 - f )/f2) );  // fraction of mass from stream1
        Real alpha2 = 1.0 - alpha1;

        YtmpO = cd.moleFracToMassFrac(XtmpO);

        Real Tnew;

        if (alpha2>0.001)
        {
            // Set enthalpy of mixture from stream2
            Real hmix, hmix1;
            CKHBMS(&Told1,Yold1.dataPtr(),&iwrk,&rwrk,&hmix1); hmix1 *= 1.e-4;
            CKHBMS(&Ttmp,YtmpO.dataPtr(),&iwrk,&rwrk,&hmix);  hmix  *= 1.e-4;
            Real hmix2 = (hmix - alpha1*hmix1)/alpha2;

            // Set mixture of stream2 fluid
            for (int i=0; i<nSpec_old; ++i)
                YtmpO[i] = (YtmpO[i] - alpha1*Yold1[i])/alpha2;
            XtmpO = cd.massFracToMoleFrac(YtmpO);
            
            // Find temperature of stream2 fluid
            Real Told2=Ttmp;
            int Niter = FORT_TfromHYpt(&Told2,&hmix2,YtmpO.dataPtr(),&errMAX,&NiterMAX,res.dataPtr());
            FORT_CONVERT_DRM19_TO_H2(&Told2, Xnew2.dataPtr());
            Tnew2 = Xnew2[Xnew2.size()-1];
            CKXTY_new(Xnew2.dataPtr(),&iwrk,&rwrk,Ynew2.dataPtr());
            CKHBMS_new(&Tnew2,Ynew2.dataPtr(),&iwrk,&rwrk,&hmix2);  hmix2  *= 1.e-4;

            // Now, blend new stream2 fluid with stream1, set total enthalpy and compute T
            for (int i=0; i<nSpec_new; ++i)
                Ynew[i] = alpha2*Ynew2[i] + alpha1*Ynew1[i];
            hmix = alpha1*hmix1 + alpha2*hmix2;
            Tnew = Tnew2; // Initial guess
            Niter = FORT_TfromHYpt_new(&Tnew,&hmix,Ynew.dataPtr(),&errMAX,&NiterMAX,res.dataPtr());
        }
        else
        {
            Tnew = Told1;
            for (int i=0; i<nSpec_new; ++i)
                Ynew[i] = Ynew1[i];
        }

        for (int i = 0; i < nSpec_new; i++)
            Y_new(iv, yComp_new + i) = Ynew[i];
        T_new(iv,TComp_new) = Tnew;
    }
}

int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    if (argc == 1) PrintUsage(argv[0]);

    ParmParse pp;
    
    if (pp.contains("help")) PrintUsage(argv[0]);
    
    FArrayBox::setFormat(FABio::FAB_NATIVE);
    
    std::string ifile, ofile;
    
    pp.query("ifile", ifile);
    if (ifile.empty()) BoxLib::Abort("You must specify `ifile'");
    
    pp.query("ofile", ofile);
    if (ofile.empty()) BoxLib::Abort("You must specify `ofile'");
    
    int verbose=0; pp.query("verbose",verbose);
    if (verbose>1) AmrData::SetVerbose(true);
    
    string tranfile="tran.asc.drm19"; pp.query("tranfile",tranfile);
    ChemDriver cd(tranfile);

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    
    DataServices dataServices(ifile, fileType);

    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    AmrData&  amrData = dataServices.AmrDataRef();

    const int idT     = amrData.StateNumber("temp");
    const int idX     = amrData.StateNumber("x_velocity");
    const int idY     = amrData.StateNumber("y_velocity");
    const int idZ     = amrData.StateNumber("z_velocity");
    const int idSpec  = amrData.StateNumber(string("X(" + cd.speciesNames()[0] + ")"));
    if (idT < 0 || idX < 0 || idY < 0 || idZ < 0 || idSpec < 0)
        BoxLib::Abort("Not all states are available");

    // Get info about new mech
    Real rwrk;
    int iwrk,nelt_new,nspec_new,nreac_new,nfit_new;
    CKINDX_new(&iwrk,&rwrk,&nelt_new,&nspec_new,&nreac_new,&nfit_new);

    int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);
    int nspec_old = cd.numSpecies();
    int ncomps = nspec_new;           // mass fractions
    ncomps    += 1;                   // temp
    ncomps    += BL_SPACEDIM;         // velocities

    // Define 2 streams, one with "atom" and one without
    // hc + a.O2 -> b.CO2 + c.H2O
    //   for hc = CH4, a=2, b=1, c=2

    string atom="C"; pp.query("atom",atom);
    Real phi=0.8; pp.query("phi",phi);
    string fuel="CH4"; pp.query("fuel",fuel);
    Real aphi=-1;
    if (fuel=="CH4")
    {
        aphi = 2.0;
    }
    else if (fuel=="C3H8")
    {
        aphi = 5.0;
    }
    else
    {
        pp.get("aphi",aphi);
    }

    Array<Real> Xold1(nspec_old,0);
    Xold1[cd.index("O2")] = 0.21;
    Xold1[cd.index("N2")] = 1.0 - Xold1[cd.index("O2")];
    Real Told1 = 300.;
    
    Array<Real> Xold2(nspec_old,0);
    Xold2[cd.index("O2")] = 1.0/(1.0 + phi/aphi  + Xold1[cd.index("N2")]/Xold1[cd.index("O2")]);
    Xold2[cd.index(fuel)] = phi * Xold2[cd.index("O2")] / aphi;
    Xold2[cd.index("N2")] = Xold1[cd.index("N2")] / Xold1[cd.index("O2")] * Xold2[cd.index("O2")];

    PArray<MultiFab>   mfout(finestLevel + 1);
    Array<std::string> newPlotNames(ncomps);

    int idx = 0;
    if (idX > -1) newPlotNames[idx++] = "x_velocity";
    if (idY > -1) newPlotNames[idx++] = "y_velocity";
    if (idZ > -1) newPlotNames[idx++] = "z_velocity";
    for (int i = 0; i < nspec_new; i++)
        newPlotNames[idx + i] = "Y(" + string(fesymname_new_(i)) + ")";
    newPlotNames[idx+nspec_new] = "temp";

    for(int lev = 0; lev <= finestLevel; ++lev)
    {
        mfout.set(lev,new MultiFab(amrData.boxArray(lev), ncomps, 0));

        MultiFab molefrac(amrData.boxArray(lev), nspec_old, 0);

        mfout[lev].setVal(0);
        molefrac.setVal(0);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Reading velocity data at level " << lev << " ...\n";
        if (idX > -1) mfout[lev].copy(amrData.GetGrids(lev,idX), 0, idX, 1);
        if (idY > -1) mfout[lev].copy(amrData.GetGrids(lev,idY), 0, idY, 1);
        if (idZ > -1) mfout[lev].copy(amrData.GetGrids(lev,idZ), 0, idZ, 1);

        amrData.FlushGrids(idX);
        amrData.FlushGrids(idY);
        amrData.FlushGrids(idZ);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Reading temperature data at level " << lev << " ...\n";
        const MultiFab& mfT = amrData.GetGrids(lev,idT);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Reading mole fractions at level " << lev << " ...\n";
        for (int i = 0; i < nspec_old; i++)
        {
            molefrac.copy(amrData.GetGrids(lev,idSpec+i), 0, i, 1);
            amrData.FlushGrids(idSpec+i);
        }

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Converting level " << lev << " ...\n";
        int yComp_new = BL_SPACEDIM;
        int TComp_new = ncomps-1;
        int xComp_old = 0;
        int TComp_old = 0;

        //Hack
#if 0
        for (int i=0; i<nspec_old; ++i)
            molefrac.setVal(Xold2[i],i,1,0);        
        MultiFab inT(amrData.boxArray(lev), 1, 0);
        inT.setVal(300);
        TComp_old = 0;
#endif

        for (MFIter mfi(molefrac); mfi.isValid(); ++mfi)
        {
            ConvertDRM19toH2(mfT[mfi], TComp_old, molefrac[mfi], xComp_old, mfout[lev][mfi],
                             yComp_new, nspec_new, mfout[lev][mfi], TComp_new, cd,
                             atom, Xold1, Told1, Xold2);
        }

        amrData.FlushGrids(idT);
    }

    const AmrData& a = amrData;
    bool verb = true;
    writePlotfile(mfout,a.Time(),a.ProbLo(),a.ProbHi(),a.RefRatio(),a.ProbDomain(),
                  a.DxLevel(),a.CoordSys(),ofile,newPlotNames,verb);

    BoxLib::Finalize();
    return 0;
}
