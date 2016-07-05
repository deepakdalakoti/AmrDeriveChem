
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

#include "convertDRM19toC3H8_F.H"
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

// Set of ChemKin funcs I need to see for the "old" mech
extern "C" {

    void CKHBMS(double * T, double * y, int * iwrk, double * rwrk, double * hbms);
};


// Set of ChemKin funcs I need for the "new" mech:
//   CKINDX, CKXTY, fesymnum_, fesymname_, speciesEnthalpy, CKHBMS, cp_R, CKCPBS

void CKINDX_new(int * iwrk, double * rwrk, int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 5;
    *kk = 39;
    *ii = 175;
    *nfit = -1; /*Why do you need this anyway ?  */
}
void CKXTY_new(double * x, int * iwrk, double * rwrk, double * y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*39.948000; /*AR */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*2.015940; /*H2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*34.014740; /*H2O2 */
    XW += x[10]*28.010400; /*CO */
    XW += x[11]*44.009800; /*CO2 */
    XW += x[12]*29.018370; /*HCO */
    XW += x[13]*30.026340; /*CH2O */
    XW += x[14]*16.042880; /*CH4 */
    XW += x[15]*15.034910; /*CH3 */
    XW += x[16]*14.026940; /*T-CH2 */
    XW += x[17]*14.026940; /*S-CH2 */
    XW += x[18]*28.053880; /*C2H4 */
    XW += x[19]*31.034310; /*CH3O */
    XW += x[20]*29.061850; /*C2H5 */
    XW += x[21]*30.069820; /*C2H6 */
    XW += x[22]*13.018970; /*CH */
    XW += x[23]*26.037940; /*C2H2 */
    XW += x[24]*27.045910; /*C2H3 */
    XW += x[25]*43.045310; /*CH2CHO */
    XW += x[26]*44.053280; /*C2H4O */
    XW += x[27]*42.037340; /*CH2CO */
    XW += x[28]*41.029370; /*HCCO */
    XW += x[29]*25.029970; /*C2H */
    XW += x[30]*31.034310; /*CH2OH */
    XW += x[31]*32.042280; /*CH3OH */
    XW += x[32]*40.064880; /*C3H4 */
    XW += x[33]*39.056910; /*C3H3 */
    XW += x[34]*41.072850; /*C3H5 */
    XW += x[35]*42.080820; /*C3H6 */
    XW += x[36]*44.096760; /*C3H8 */
    XW += x[37]*43.088790; /*I-C3H7 */
    XW += x[38]*43.088790; /*N-C3H7 */
    /*Now compute conversion */
    y[0] = x[0]*28.013400/XW; 
    y[1] = x[1]*39.948000/XW; 
    y[2] = x[2]*1.007970/XW; 
    y[3] = x[3]*31.998800/XW; 
    y[4] = x[4]*17.007370/XW; 
    y[5] = x[5]*15.999400/XW; 
    y[6] = x[6]*2.015940/XW; 
    y[7] = x[7]*18.015340/XW; 
    y[8] = x[8]*33.006770/XW; 
    y[9] = x[9]*34.014740/XW; 
    y[10] = x[10]*28.010400/XW; 
    y[11] = x[11]*44.009800/XW; 
    y[12] = x[12]*29.018370/XW; 
    y[13] = x[13]*30.026340/XW; 
    y[14] = x[14]*16.042880/XW; 
    y[15] = x[15]*15.034910/XW; 
    y[16] = x[16]*14.026940/XW; 
    y[17] = x[17]*14.026940/XW; 
    y[18] = x[18]*28.053880/XW; 
    y[19] = x[19]*31.034310/XW; 
    y[20] = x[20]*29.061850/XW; 
    y[21] = x[21]*30.069820/XW; 
    y[22] = x[22]*13.018970/XW; 
    y[23] = x[23]*26.037940/XW; 
    y[24] = x[24]*27.045910/XW; 
    y[25] = x[25]*43.045310/XW; 
    y[26] = x[26]*44.053280/XW; 
    y[27] = x[27]*42.037340/XW; 
    y[28] = x[28]*41.029370/XW; 
    y[29] = x[29]*25.029970/XW; 
    y[30] = x[30]*31.034310/XW; 
    y[31] = x[31]*32.042280/XW; 
    y[32] = x[32]*40.064880/XW; 
    y[33] = x[33]*39.056910/XW; 
    y[34] = x[34]*41.072850/XW; 
    y[35] = x[35]*42.080820/XW; 
    y[36] = x[36]*44.096760/XW; 
    y[37] = x[37]*43.088790/XW; 
    y[38] = x[38]*43.088790/XW; 

    return;
}
int fesymnum_new_(const char* s1)
{
    if (strcmp(s1, "N2")==0) return 0; 
    if (strcmp(s1, "AR")==0) return 1; 
    if (strcmp(s1, "H")==0) return 2; 
    if (strcmp(s1, "O2")==0) return 3; 
    if (strcmp(s1, "OH")==0) return 4; 
    if (strcmp(s1, "O")==0) return 5; 
    if (strcmp(s1, "H2")==0) return 6; 
    if (strcmp(s1, "H2O")==0) return 7; 
    if (strcmp(s1, "HO2")==0) return 8; 
    if (strcmp(s1, "H2O2")==0) return 9; 
    if (strcmp(s1, "CO")==0) return 10; 
    if (strcmp(s1, "CO2")==0) return 11; 
    if (strcmp(s1, "HCO")==0) return 12; 
    if (strcmp(s1, "CH2O")==0) return 13; 
    if (strcmp(s1, "CH4")==0) return 14; 
    if (strcmp(s1, "CH3")==0) return 15; 
    if (strcmp(s1, "T-CH2")==0) return 16; 
    if (strcmp(s1, "S-CH2")==0) return 17; 
    if (strcmp(s1, "C2H4")==0) return 18; 
    if (strcmp(s1, "CH3O")==0) return 19; 
    if (strcmp(s1, "C2H5")==0) return 20; 
    if (strcmp(s1, "C2H6")==0) return 21; 
    if (strcmp(s1, "CH")==0) return 22; 
    if (strcmp(s1, "C2H2")==0) return 23; 
    if (strcmp(s1, "C2H3")==0) return 24; 
    if (strcmp(s1, "CH2CHO")==0) return 25; 
    if (strcmp(s1, "C2H4O")==0) return 26; 
    if (strcmp(s1, "CH2CO")==0) return 27; 
    if (strcmp(s1, "HCCO")==0) return 28; 
    if (strcmp(s1, "C2H")==0) return 29; 
    if (strcmp(s1, "CH2OH")==0) return 30; 
    if (strcmp(s1, "CH3OH")==0) return 31; 
    if (strcmp(s1, "C3H4")==0) return 32; 
    if (strcmp(s1, "C3H3")==0) return 33; 
    if (strcmp(s1, "C3H5")==0) return 34; 
    if (strcmp(s1, "C3H6")==0) return 35; 
    if (strcmp(s1, "C3H8")==0) return 36; 
    if (strcmp(s1, "I-C3H7")==0) return 37; 
    if (strcmp(s1, "N-C3H7")==0) return 38; 
    /*species name not found */
    return -1;
}
char* fesymname_new_(int sn)
{
    if (sn==0) return "N2"; 
    if (sn==1) return "AR"; 
    if (sn==2) return "H"; 
    if (sn==3) return "O2"; 
    if (sn==4) return "OH"; 
    if (sn==5) return "O"; 
    if (sn==6) return "H2"; 
    if (sn==7) return "H2O"; 
    if (sn==8) return "HO2"; 
    if (sn==9) return "H2O2"; 
    if (sn==10) return "CO"; 
    if (sn==11) return "CO2"; 
    if (sn==12) return "HCO"; 
    if (sn==13) return "CH2O"; 
    if (sn==14) return "CH4"; 
    if (sn==15) return "CH3"; 
    if (sn==16) return "T-CH2"; 
    if (sn==17) return "S-CH2"; 
    if (sn==18) return "C2H4"; 
    if (sn==19) return "CH3O"; 
    if (sn==20) return "C2H5"; 
    if (sn==21) return "C2H6"; 
    if (sn==22) return "CH"; 
    if (sn==23) return "C2H2"; 
    if (sn==24) return "C2H3"; 
    if (sn==25) return "CH2CHO"; 
    if (sn==26) return "C2H4O"; 
    if (sn==27) return "CH2CO"; 
    if (sn==28) return "HCCO"; 
    if (sn==29) return "C2H"; 
    if (sn==30) return "CH2OH"; 
    if (sn==31) return "CH3OH"; 
    if (sn==32) return "C3H4"; 
    if (sn==33) return "C3H3"; 
    if (sn==34) return "C3H5"; 
    if (sn==35) return "C3H6"; 
    if (sn==36) return "C3H8"; 
    if (sn==37) return "I-C3H7"; 
    if (sn==38) return "N-C3H7"; 
    /*species name not found */
    return "NOTFOUND";
}
void speciesEnthalpy_new(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 / tc[1];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 2: H */
        species[2] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 / tc[1];
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
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 / tc[1];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 / tc[1];
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 / tc[1];
        /*species 7: H2O */
        species[7] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 / tc[1];
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 / tc[1];
        /*species 9: H2O2 */
        species[9] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 / tc[1];
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 / tc[1];
        /*species 11: CO2 */
        species[11] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 / tc[1];
        /*species 12: HCO */
        species[12] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 / tc[1];
        /*species 13: CH2O */
        species[13] =
            +4.79372315e+00
            -4.95416684e-03 * tc[1]
            +1.24406669e-05 * tc[2]
            -9.48213152e-09 * tc[3]
            +2.63545304e-12 * tc[4]
            -1.43089567e+04 / tc[1];
        /*species 14: CH4 */
        species[14] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 / tc[1];
        /*species 15: CH3 */
        species[15] =
            +3.67359040e+00
            +1.00547588e-03 * tc[1]
            +1.91007285e-06 * tc[2]
            -1.71779356e-09 * tc[3]
            +5.08771468e-13 * tc[4]
            +1.64449988e+04 / tc[1];
        /*species 16: T-CH2 */
        species[16] =
            +3.76267867e+00
            +4.84436072e-04 * tc[1]
            +9.31632803e-07 * tc[2]
            -9.62727883e-10 * tc[3]
            +3.37483438e-13 * tc[4]
            +4.60040401e+04 / tc[1];
        /*species 17: S-CH2 */
        species[17] =
            +4.19860411e+00
            -1.18330710e-03 * tc[1]
            +2.74432073e-06 * tc[2]
            -1.67203995e-09 * tc[3]
            +3.88629474e-13 * tc[4]
            +5.04968163e+04 / tc[1];
        /*species 18: C2H4 */
        species[18] =
            +3.95920148e+00
            -3.78526124e-03 * tc[1]
            +1.90330097e-05 * tc[2]
            -1.72897188e-08 * tc[3]
            +5.39768746e-12 * tc[4]
            +5.08977593e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +3.71180502e+00
            -1.40231653e-03 * tc[1]
            +1.25516990e-05 * tc[2]
            -1.18268022e-08 * tc[3]
            +3.73176840e-12 * tc[4]
            +1.30772484e+03 / tc[1];
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00
            -2.09329446e-03 * tc[1]
            +1.65714269e-05 * tc[2]
            -1.49781651e-08 * tc[3]
            +4.61018008e-12 * tc[4]
            +1.28416265e+04 / tc[1];
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -2.75077135e-03 * tc[1]
            +1.99812763e-05 * tc[2]
            -1.77116571e-08 * tc[3]
            +5.37371542e-12 * tc[4]
            -1.15222055e+04 / tc[1];
        /*species 22: CH */
        species[22] =
            +3.48981665e+00
            +1.61917771e-04 * tc[1]
            -5.62996883e-07 * tc[2]
            +7.90543317e-10 * tc[3]
            -2.81218134e-13 * tc[4]
            +7.07972934e+04 / tc[1];
        /*species 23: C2H2 */
        species[23] =
            +8.08681094e-01
            +1.16807815e-02 * tc[1]
            -1.18390605e-05 * tc[2]
            +7.00381092e-09 * tc[3]
            -1.70014595e-12 * tc[4]
            +2.64289807e+04 / tc[1];
        /*species 24: C2H3 */
        species[24] =
            +3.21246645e+00
            +7.57395810e-04 * tc[1]
            +8.64031373e-06 * tc[2]
            -8.94144617e-09 * tc[3]
            +2.94301746e-12 * tc[4]
            +3.48598468e+04 / tc[1];
        /*species 25: CH2CHO */
        species[25] =
            +1.01340010e+00
            +1.13407335e-02 * tc[1]
            -5.24464800e-06 * tc[2]
            +1.01228758e-09 * tc[3]
            +5.91980240e-14 * tc[4]
            +3.80428530e+02 / tc[1];
        /*species 26: C2H4O */
        species[26] =
            +4.72945950e+00
            -1.59664290e-03 * tc[1]
            +1.58449737e-05 * tc[2]
            -1.43646527e-08 * tc[3]
            +4.38622240e-12 * tc[4]
            -2.15728780e+04 / tc[1];
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00
            +9.05943605e-03 * tc[1]
            -5.79824913e-06 * tc[2]
            +2.33599392e-09 * tc[3]
            -4.02915230e-13 * tc[4]
            -7.04291804e+03 / tc[1];
        /*species 28: HCCO */
        species[28] =
            +2.25172140e+00
            +8.82751050e-03 * tc[1]
            -7.90970033e-06 * tc[2]
            +4.31893975e-09 * tc[3]
            -1.01329622e-12 * tc[4]
            +2.00594490e+04 / tc[1];
        /*species 29: C2H */
        species[29] =
            +2.88965733e+00
            +6.70498055e-03 * tc[1]
            -9.49231670e-06 * tc[2]
            +7.36977613e-09 * tc[3]
            -2.18663022e-12 * tc[4]
            +6.68393932e+04 / tc[1];
        /*species 30: CH2OH */
        species[30] =
            +4.47832317e+00
            -6.75348435e-04 * tc[1]
            +9.28279023e-06 * tc[2]
            -9.12168492e-09 * tc[3]
            +2.95813550e-12 * tc[4]
            -3.52476728e+03 / tc[1];
        /*species 31: CH3OH */
        species[31] =
            +5.71539582e+00
            -7.61545645e-03 * tc[1]
            +2.17480385e-05 * tc[2]
            -1.77701722e-08 * tc[3]
            +5.22705396e-12 * tc[4]
            -2.56427656e+04 / tc[1];
        /*species 32: C3H4 */
        species[32] =
            +2.61304450e+00
            +6.06128750e-03 * tc[1]
            +6.17996000e-06 * tc[2]
            -8.63128725e-09 * tc[3]
            +3.06701580e-12 * tc[4]
            +2.15415670e+04 / tc[1];
        /*species 33: C3H3 */
        species[33] =
            +1.35110927e+00
            +1.63705612e-02 * tc[1]
            -1.57942378e-05 * tc[2]
            +9.40774520e-09 * tc[3]
            -2.37081846e-12 * tc[4]
            +4.01057783e+04 / tc[1];
        /*species 34: C3H5 */
        species[34] =
            +1.36318350e+00
            +9.90691050e-03 * tc[1]
            +4.16568667e-06 * tc[2]
            -8.33888875e-09 * tc[3]
            +3.16931420e-12 * tc[4]
            +1.92456290e+04 / tc[1];
        /*species 35: C3H6 */
        species[35] =
            +1.49330700e+00
            +1.04625900e-02 * tc[1]
            +1.49559800e-06 * tc[2]
            -4.17228000e-09 * tc[3]
            +1.43162920e-12 * tc[4]
            +1.07482600e+03 / tc[1];
        /*species 36: C3H8 */
        species[36] =
            +9.28510930e-01
            +1.32302830e-02 * tc[1]
            +2.01108153e-06 * tc[2]
            -5.47873825e-09 * tc[3]
            +1.89923088e-12 * tc[4]
            -1.40579070e+04 / tc[1];
        /*species 37: I-C3H7 */
        species[37] =
            +1.44491990e+00
            +1.04995560e-02 * tc[1]
            +2.56787407e-06 * tc[2]
            -4.61906325e-09 * tc[3]
            +1.42565924e-12 * tc[4]
            +9.42237240e+03 / tc[1];
        /*species 38: N-C3H7 */
        species[38] =
            +1.04911730e+00
            +1.30044865e-02 * tc[1]
            +7.84750533e-07 * tc[2]
            -4.89878300e-09 * tc[3]
            +1.87440414e-12 * tc[4]
            +1.03123460e+04 / tc[1];
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 / tc[1];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 / tc[1];
        /*species 2: H */
        species[2] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 / tc[1];
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
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 / tc[1];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 / tc[1];
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 / tc[1];
        /*species 7: H2O */
        species[7] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 / tc[1];
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 / tc[1];
        /*species 9: H2O2 */
        species[9] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 / tc[1];
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 / tc[1];
        /*species 11: CO2 */
        species[11] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 / tc[1];
        /*species 12: HCO */
        species[12] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 / tc[1];
        /*species 13: CH2O */
        species[13] =
            +1.76069008e+00
            +4.60000041e-03 * tc[1]
            -1.47419604e-06 * tc[2]
            +2.51603030e-10 * tc[3]
            -1.76771128e-14 * tc[4]
            -1.39958323e+04 / tc[1];
        /*species 14: CH4 */
        species[14] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 / tc[1];
        /*species 15: CH3 */
        species[15] =
            +2.28571772e+00
            +3.61995018e-03 * tc[1]
            -9.95714493e-07 * tc[2]
            +1.48921161e-10 * tc[3]
            -9.34308788e-15 * tc[4]
            +1.67755843e+04 / tc[1];
        /*species 16: T-CH2 */
        species[16] =
            +2.87410113e+00
            +1.82819646e-03 * tc[1]
            -4.69648657e-07 * tc[2]
            +6.50448872e-11 * tc[3]
            -3.75455134e-15 * tc[4]
            +4.62636040e+04 / tc[1];
        /*species 17: S-CH2 */
        species[17] =
            +2.29203842e+00
            +2.32794318e-03 * tc[1]
            -6.70639823e-07 * tc[2]
            +1.04476500e-10 * tc[3]
            -6.79432730e-15 * tc[4]
            +5.09259997e+04 / tc[1];
        /*species 18: C2H4 */
        species[18] =
            +2.03611116e+00
            +7.32270755e-03 * tc[1]
            -2.23692638e-06 * tc[2]
            +3.68057308e-10 * tc[3]
            -2.51412122e-14 * tc[4]
            +4.93988614e+03 / tc[1];
        /*species 19: CH3O */
        species[19] =
            +4.75779238e+00
            +3.72071237e-03 * tc[1]
            -8.99017253e-07 * tc[2]
            +1.09522626e-10 * tc[3]
            -5.27074196e-15 * tc[4]
            +3.90139164e+02 / tc[1];
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00
            +8.69863610e-03 * tc[1]
            -2.66068889e-06 * tc[2]
            +4.38044223e-10 * tc[3]
            -2.99283152e-14 * tc[4]
            +1.28575200e+04 / tc[1];
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +1.08426339e-02 * tc[1]
            -3.34186890e-06 * tc[2]
            +5.53530003e-10 * tc[3]
            -3.80005780e-14 * tc[4]
            -1.14263932e+04 / tc[1];
        /*species 22: CH */
        species[22] =
            +2.87846473e+00
            +4.85456840e-04 * tc[1]
            +4.81485517e-08 * tc[2]
            -3.26719623e-11 * tc[3]
            +3.52158766e-15 * tc[4]
            +7.10124364e+04 / tc[1];
        /*species 23: C2H2 */
        species[23] =
            +4.14756964e+00
            +2.98083332e-03 * tc[1]
            -7.90982840e-07 * tc[2]
            +1.16853043e-10 * tc[3]
            -7.22470426e-15 * tc[4]
            +2.59359992e+04 / tc[1];
        /*species 24: C2H3 */
        species[24] =
            +3.01672400e+00
            +5.16511460e-03 * tc[1]
            -1.56027450e-06 * tc[2]
            +2.54408220e-10 * tc[3]
            -1.72521408e-14 * tc[4]
            +3.46128739e+04 / tc[1];
        /*species 25: CH2CHO */
        species[25] =
            +5.16620060e+00
            +5.42391300e-03 * tc[1]
            -1.48861227e-06 * tc[2]
            +2.01571370e-10 * tc[3]
            -9.68203860e-15 * tc[4]
            -7.31993470e+02 / tc[1];
        /*species 26: C2H4O */
        species[26] =
            +5.40411080e+00
            +5.86152950e-03 * tc[1]
            -1.40877123e-06 * tc[2]
            +1.70931128e-10 * tc[3]
            -8.19697260e-15 * tc[4]
            -2.25931220e+04 / tc[1];
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00
            +4.50179872e-03 * tc[1]
            -1.38979878e-06 * tc[2]
            +2.30836470e-10 * tc[3]
            -1.58967640e-14 * tc[4]
            -7.55105311e+03 / tc[1];
        /*species 28: HCCO */
        species[28] =
            +5.62820580e+00
            +2.04267005e-03 * tc[1]
            -5.31151567e-07 * tc[2]
            +7.15651300e-11 * tc[3]
            -3.88156640e-15 * tc[4]
            +1.93272150e+04 / tc[1];
        /*species 29: C2H */
        species[29] =
            +3.16780652e+00
            +2.37610951e-03 * tc[1]
            -6.12623590e-07 * tc[2]
            +7.60475630e-11 * tc[3]
            -3.54465540e-15 * tc[4]
            +6.71210650e+04 / tc[1];
        /*species 30: CH2OH */
        species[30] =
            +5.09312037e+00
            +2.97379275e-03 * tc[1]
            -6.88321747e-07 * tc[2]
            +8.07516758e-11 * tc[3]
            -3.76250104e-15 * tc[4]
            -4.05813228e+03 / tc[1];
        /*species 31: CH3OH */
        species[31] =
            +1.78970791e+00
            +7.04691460e-03 * tc[1]
            -2.12166945e-06 * tc[2]
            +3.45427713e-10 * tc[3]
            -2.34120440e-14 * tc[4]
            -2.53748747e+04 / tc[1];
        /*species 32: C3H4 */
        species[32] =
            +6.31687220e+00
            +5.56686400e-03 * tc[1]
            -1.32097927e-06 * tc[2]
            +1.58910595e-10 * tc[3]
            -7.57510800e-15 * tc[4]
            +2.01174950e+04 / tc[1];
        /*species 33: C3H3 */
        species[33] =
            +7.14221880e+00
            +3.80951002e-03 * tc[1]
            -8.91533167e-07 * tc[2]
            +1.06228700e-10 * tc[3]
            -5.02950830e-15 * tc[4]
            +3.89087427e+04 / tc[1];
        /*species 34: C3H5 */
        species[34] =
            +6.50078770e+00
            +7.16236550e-03 * tc[1]
            -1.89272107e-06 * tc[2]
            +2.77020025e-10 * tc[3]
            -1.80727774e-14 * tc[4]
            +1.74824490e+04 / tc[1];
        /*species 35: C3H6 */
        species[35] =
            +6.73225700e+00
            +7.45417000e-03 * tc[1]
            -1.64996633e-06 * tc[2]
            +1.80300550e-10 * tc[3]
            -7.53240800e-15 * tc[4]
            -9.23570300e+02 / tc[1];
        /*species 36: C3H8 */
        species[36] =
            +7.52441520e+00
            +9.44914100e-03 * tc[1]
            -2.09736803e-06 * tc[2]
            +2.30403642e-10 * tc[3]
            -9.73689560e-15 * tc[4]
            -1.65643940e+04 / tc[1];
        /*species 37: I-C3H7 */
        species[37] =
            +6.51927410e+00
            +8.61005200e-03 * tc[1]
            -1.91214057e-06 * tc[2]
            +2.10326830e-10 * tc[3]
            -8.91318260e-15 * tc[4]
            +7.32271930e+03 / tc[1];
        /*species 38: N-C3H7 */
        species[38] =
            +7.70974790e+00
            +8.01574250e-03 * tc[1]
            -1.75734127e-06 * tc[2]
            +1.89720880e-10 * tc[3]
            -7.77254380e-15 * tc[4]
            +7.97622360e+03 / tc[1];
    }
    return;
}
void CKHBMS_new(double *T, double *y, int * iwrk, double * rwrk, double * hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[39]; /* temporary storage */
    double RT = 8.314e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*hml[0]/28.013400; /*N2 */
    result += y[1]*hml[1]/39.948000; /*AR */
    result += y[2]*hml[2]/1.007970; /*H */
    result += y[3]*hml[3]/31.998800; /*O2 */
    result += y[4]*hml[4]/17.007370; /*OH */
    result += y[5]*hml[5]/15.999400; /*O */
    result += y[6]*hml[6]/2.015940; /*H2 */
    result += y[7]*hml[7]/18.015340; /*H2O */
    result += y[8]*hml[8]/33.006770; /*HO2 */
    result += y[9]*hml[9]/34.014740; /*H2O2 */
    result += y[10]*hml[10]/28.010400; /*CO */
    result += y[11]*hml[11]/44.009800; /*CO2 */
    result += y[12]*hml[12]/29.018370; /*HCO */
    result += y[13]*hml[13]/30.026340; /*CH2O */
    result += y[14]*hml[14]/16.042880; /*CH4 */
    result += y[15]*hml[15]/15.034910; /*CH3 */
    result += y[16]*hml[16]/14.026940; /*T-CH2 */
    result += y[17]*hml[17]/14.026940; /*S-CH2 */
    result += y[18]*hml[18]/28.053880; /*C2H4 */
    result += y[19]*hml[19]/31.034310; /*CH3O */
    result += y[20]*hml[20]/29.061850; /*C2H5 */
    result += y[21]*hml[21]/30.069820; /*C2H6 */
    result += y[22]*hml[22]/13.018970; /*CH */
    result += y[23]*hml[23]/26.037940; /*C2H2 */
    result += y[24]*hml[24]/27.045910; /*C2H3 */
    result += y[25]*hml[25]/43.045310; /*CH2CHO */
    result += y[26]*hml[26]/44.053280; /*C2H4O */
    result += y[27]*hml[27]/42.037340; /*CH2CO */
    result += y[28]*hml[28]/41.029370; /*HCCO */
    result += y[29]*hml[29]/25.029970; /*C2H */
    result += y[30]*hml[30]/31.034310; /*CH2OH */
    result += y[31]*hml[31]/32.042280; /*CH3OH */
    result += y[32]*hml[32]/40.064880; /*C3H4 */
    result += y[33]*hml[33]/39.056910; /*C3H3 */
    result += y[34]*hml[34]/41.072850; /*C3H5 */
    result += y[35]*hml[35]/42.080820; /*C3H6 */
    result += y[36]*hml[36]/44.096760; /*C3H8 */
    result += y[37]*hml[37]/43.088790; /*I-C3H7 */
    result += y[38]*hml[38]/43.088790; /*N-C3H7 */

    *hbms = result * RT;
}
void cp_R_new(double * species, double * tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: N2 */
        species[0] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: H2 */
        species[6] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 7: H2O */
        species[7] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: CO2 */
        species[11] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 12: HCO */
        species[12] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 13: CH2O */
        species[13] =
            +4.79372315e+00
            -9.90833369e-03 * tc[1]
            +3.73220008e-05 * tc[2]
            -3.79285261e-08 * tc[3]
            +1.31772652e-11 * tc[4];
        /*species 14: CH4 */
        species[14] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +3.67359040e+00
            +2.01095175e-03 * tc[1]
            +5.73021856e-06 * tc[2]
            -6.87117425e-09 * tc[3]
            +2.54385734e-12 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +3.76267867e+00
            +9.68872143e-04 * tc[1]
            +2.79489841e-06 * tc[2]
            -3.85091153e-09 * tc[3]
            +1.68741719e-12 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +4.19860411e+00
            -2.36661419e-03 * tc[1]
            +8.23296220e-06 * tc[2]
            -6.68815981e-09 * tc[3]
            +1.94314737e-12 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +3.95920148e+00
            -7.57052247e-03 * tc[1]
            +5.70990292e-05 * tc[2]
            -6.91588753e-08 * tc[3]
            +2.69884373e-11 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +3.71180502e+00
            -2.80463306e-03 * tc[1]
            +3.76550971e-05 * tc[2]
            -4.73072089e-08 * tc[3]
            +1.86588420e-11 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +4.30646568e+00
            -4.18658892e-03 * tc[1]
            +4.97142807e-05 * tc[2]
            -5.99126606e-08 * tc[3]
            +2.30509004e-11 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +4.29142492e+00
            -5.50154270e-03 * tc[1]
            +5.99438288e-05 * tc[2]
            -7.08466285e-08 * tc[3]
            +2.68685771e-11 * tc[4];
        /*species 22: CH */
        species[22] =
            +3.48981665e+00
            +3.23835541e-04 * tc[1]
            -1.68899065e-06 * tc[2]
            +3.16217327e-09 * tc[3]
            -1.40609067e-12 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +8.08681094e-01
            +2.33615629e-02 * tc[1]
            -3.55171815e-05 * tc[2]
            +2.80152437e-08 * tc[3]
            -8.50072974e-12 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.21246645e+00
            +1.51479162e-03 * tc[1]
            +2.59209412e-05 * tc[2]
            -3.57657847e-08 * tc[3]
            +1.47150873e-11 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +1.01340010e+00
            +2.26814670e-02 * tc[1]
            -1.57339440e-05 * tc[2]
            +4.04915030e-09 * tc[3]
            +2.95990120e-13 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            +4.72945950e+00
            -3.19328580e-03 * tc[1]
            +4.75349210e-05 * tc[2]
            -5.74586110e-08 * tc[3]
            +2.19311120e-11 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +2.13583630e+00
            +1.81188721e-02 * tc[1]
            -1.73947474e-05 * tc[2]
            +9.34397568e-09 * tc[3]
            -2.01457615e-12 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +2.25172140e+00
            +1.76550210e-02 * tc[1]
            -2.37291010e-05 * tc[2]
            +1.72757590e-08 * tc[3]
            -5.06648110e-12 * tc[4];
        /*species 29: C2H */
        species[29] =
            +2.88965733e+00
            +1.34099611e-02 * tc[1]
            -2.84769501e-05 * tc[2]
            +2.94791045e-08 * tc[3]
            -1.09331511e-11 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            +4.47832317e+00
            -1.35069687e-03 * tc[1]
            +2.78483707e-05 * tc[2]
            -3.64867397e-08 * tc[3]
            +1.47906775e-11 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            +5.71539582e+00
            -1.52309129e-02 * tc[1]
            +6.52441155e-05 * tc[2]
            -7.10806889e-08 * tc[3]
            +2.61352698e-11 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +2.61304450e+00
            +1.21225750e-02 * tc[1]
            +1.85398800e-05 * tc[2]
            -3.45251490e-08 * tc[3]
            +1.53350790e-11 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +1.35110927e+00
            +3.27411223e-02 * tc[1]
            -4.73827135e-05 * tc[2]
            +3.76309808e-08 * tc[3]
            -1.18540923e-11 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +1.36318350e+00
            +1.98138210e-02 * tc[1]
            +1.24970600e-05 * tc[2]
            -3.33555550e-08 * tc[3]
            +1.58465710e-11 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +1.49330700e+00
            +2.09251800e-02 * tc[1]
            +4.48679400e-06 * tc[2]
            -1.66891200e-08 * tc[3]
            +7.15814600e-12 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            +9.28510930e-01
            +2.64605660e-02 * tc[1]
            +6.03324460e-06 * tc[2]
            -2.19149530e-08 * tc[3]
            +9.49615440e-12 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +1.44491990e+00
            +2.09991120e-02 * tc[1]
            +7.70362220e-06 * tc[2]
            -1.84762530e-08 * tc[3]
            +7.12829620e-12 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +1.04911730e+00
            +2.60089730e-02 * tc[1]
            +2.35425160e-06 * tc[2]
            -1.95951320e-08 * tc[3]
            +9.37202070e-12 * tc[4];
    } else {
        /*species 0: N2 */
        species[0] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 1: AR */
        species[1] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 2: H */
        species[2] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: H2 */
        species[6] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 7: H2O */
        species[7] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O2 */
        species[9] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: CO2 */
        species[11] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 12: HCO */
        species[12] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 13: CH2O */
        species[13] =
            +1.76069008e+00
            +9.20000082e-03 * tc[1]
            -4.42258813e-06 * tc[2]
            +1.00641212e-09 * tc[3]
            -8.83855640e-14 * tc[4];
        /*species 14: CH4 */
        species[14] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 15: CH3 */
        species[15] =
            +2.28571772e+00
            +7.23990037e-03 * tc[1]
            -2.98714348e-06 * tc[2]
            +5.95684644e-10 * tc[3]
            -4.67154394e-14 * tc[4];
        /*species 16: T-CH2 */
        species[16] =
            +2.87410113e+00
            +3.65639292e-03 * tc[1]
            -1.40894597e-06 * tc[2]
            +2.60179549e-10 * tc[3]
            -1.87727567e-14 * tc[4];
        /*species 17: S-CH2 */
        species[17] =
            +2.29203842e+00
            +4.65588637e-03 * tc[1]
            -2.01191947e-06 * tc[2]
            +4.17906000e-10 * tc[3]
            -3.39716365e-14 * tc[4];
        /*species 18: C2H4 */
        species[18] =
            +2.03611116e+00
            +1.46454151e-02 * tc[1]
            -6.71077915e-06 * tc[2]
            +1.47222923e-09 * tc[3]
            -1.25706061e-13 * tc[4];
        /*species 19: CH3O */
        species[19] =
            +4.75779238e+00
            +7.44142474e-03 * tc[1]
            -2.69705176e-06 * tc[2]
            +4.38090504e-10 * tc[3]
            -2.63537098e-14 * tc[4];
        /*species 20: C2H5 */
        species[20] =
            +1.95465642e+00
            +1.73972722e-02 * tc[1]
            -7.98206668e-06 * tc[2]
            +1.75217689e-09 * tc[3]
            -1.49641576e-13 * tc[4];
        /*species 21: C2H6 */
        species[21] =
            +1.07188150e+00
            +2.16852677e-02 * tc[1]
            -1.00256067e-05 * tc[2]
            +2.21412001e-09 * tc[3]
            -1.90002890e-13 * tc[4];
        /*species 22: CH */
        species[22] =
            +2.87846473e+00
            +9.70913681e-04 * tc[1]
            +1.44445655e-07 * tc[2]
            -1.30687849e-10 * tc[3]
            +1.76079383e-14 * tc[4];
        /*species 23: C2H2 */
        species[23] =
            +4.14756964e+00
            +5.96166664e-03 * tc[1]
            -2.37294852e-06 * tc[2]
            +4.67412171e-10 * tc[3]
            -3.61235213e-14 * tc[4];
        /*species 24: C2H3 */
        species[24] =
            +3.01672400e+00
            +1.03302292e-02 * tc[1]
            -4.68082349e-06 * tc[2]
            +1.01763288e-09 * tc[3]
            -8.62607041e-14 * tc[4];
        /*species 25: CH2CHO */
        species[25] =
            +5.16620060e+00
            +1.08478260e-02 * tc[1]
            -4.46583680e-06 * tc[2]
            +8.06285480e-10 * tc[3]
            -4.84101930e-14 * tc[4];
        /*species 26: C2H4O */
        species[26] =
            +5.40411080e+00
            +1.17230590e-02 * tc[1]
            -4.22631370e-06 * tc[2]
            +6.83724510e-10 * tc[3]
            -4.09848630e-14 * tc[4];
        /*species 27: CH2CO */
        species[27] =
            +4.51129732e+00
            +9.00359745e-03 * tc[1]
            -4.16939635e-06 * tc[2]
            +9.23345882e-10 * tc[3]
            -7.94838201e-14 * tc[4];
        /*species 28: HCCO */
        species[28] =
            +5.62820580e+00
            +4.08534010e-03 * tc[1]
            -1.59345470e-06 * tc[2]
            +2.86260520e-10 * tc[3]
            -1.94078320e-14 * tc[4];
        /*species 29: C2H */
        species[29] =
            +3.16780652e+00
            +4.75221902e-03 * tc[1]
            -1.83787077e-06 * tc[2]
            +3.04190252e-10 * tc[3]
            -1.77232770e-14 * tc[4];
        /*species 30: CH2OH */
        species[30] =
            +5.09312037e+00
            +5.94758550e-03 * tc[1]
            -2.06496524e-06 * tc[2]
            +3.23006703e-10 * tc[3]
            -1.88125052e-14 * tc[4];
        /*species 31: CH3OH */
        species[31] =
            +1.78970791e+00
            +1.40938292e-02 * tc[1]
            -6.36500835e-06 * tc[2]
            +1.38171085e-09 * tc[3]
            -1.17060220e-13 * tc[4];
        /*species 32: C3H4 */
        species[32] =
            +6.31687220e+00
            +1.11337280e-02 * tc[1]
            -3.96293780e-06 * tc[2]
            +6.35642380e-10 * tc[3]
            -3.78755400e-14 * tc[4];
        /*species 33: C3H3 */
        species[33] =
            +7.14221880e+00
            +7.61902005e-03 * tc[1]
            -2.67459950e-06 * tc[2]
            +4.24914801e-10 * tc[3]
            -2.51475415e-14 * tc[4];
        /*species 34: C3H5 */
        species[34] =
            +6.50078770e+00
            +1.43247310e-02 * tc[1]
            -5.67816320e-06 * tc[2]
            +1.10808010e-09 * tc[3]
            -9.03638870e-14 * tc[4];
        /*species 35: C3H6 */
        species[35] =
            +6.73225700e+00
            +1.49083400e-02 * tc[1]
            -4.94989900e-06 * tc[2]
            +7.21202200e-10 * tc[3]
            -3.76620400e-14 * tc[4];
        /*species 36: C3H8 */
        species[36] =
            +7.52441520e+00
            +1.88982820e-02 * tc[1]
            -6.29210410e-06 * tc[2]
            +9.21614570e-10 * tc[3]
            -4.86844780e-14 * tc[4];
        /*species 37: I-C3H7 */
        species[37] =
            +6.51927410e+00
            +1.72201040e-02 * tc[1]
            -5.73642170e-06 * tc[2]
            +8.41307320e-10 * tc[3]
            -4.45659130e-14 * tc[4];
        /*species 38: N-C3H7 */
        species[38] =
            +7.70974790e+00
            +1.60314850e-02 * tc[1]
            -5.27202380e-06 * tc[2]
            +7.58883520e-10 * tc[3]
            -3.88627190e-14 * tc[4];
    }
    return;
}
void CKCPBS_new(double *T, double *y, int * iwrk, double * rwrk, double * cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[39]; /* temporary storage */
    cp_R(cpor, tc);
    /*multiply by y/molecularweight */
    result += cpor[0]*y[0]/28.0134; /*N2 */
    result += cpor[1]*y[1]/39.948; /*AR */
    result += cpor[2]*y[2]/1.00797; /*H */
    result += cpor[3]*y[3]/31.9988; /*O2 */
    result += cpor[4]*y[4]/17.0074; /*OH */
    result += cpor[5]*y[5]/15.9994; /*O */
    result += cpor[6]*y[6]/2.01594; /*H2 */
    result += cpor[7]*y[7]/18.0153; /*H2O */
    result += cpor[8]*y[8]/33.0068; /*HO2 */
    result += cpor[9]*y[9]/34.0147; /*H2O2 */
    result += cpor[10]*y[10]/28.0104; /*CO */
    result += cpor[11]*y[11]/44.0098; /*CO2 */
    result += cpor[12]*y[12]/29.0184; /*HCO */
    result += cpor[13]*y[13]/30.0263; /*CH2O */
    result += cpor[14]*y[14]/16.0429; /*CH4 */
    result += cpor[15]*y[15]/15.0349; /*CH3 */
    result += cpor[16]*y[16]/14.0269; /*T-CH2 */
    result += cpor[17]*y[17]/14.0269; /*S-CH2 */
    result += cpor[18]*y[18]/28.0539; /*C2H4 */
    result += cpor[19]*y[19]/31.0343; /*CH3O */
    result += cpor[20]*y[20]/29.0618; /*C2H5 */
    result += cpor[21]*y[21]/30.0698; /*C2H6 */
    result += cpor[22]*y[22]/13.019; /*CH */
    result += cpor[23]*y[23]/26.0379; /*C2H2 */
    result += cpor[24]*y[24]/27.0459; /*C2H3 */
    result += cpor[25]*y[25]/43.0453; /*CH2CHO */
    result += cpor[26]*y[26]/44.0533; /*C2H4O */
    result += cpor[27]*y[27]/42.0373; /*CH2CO */
    result += cpor[28]*y[28]/41.0294; /*HCCO */
    result += cpor[29]*y[29]/25.03; /*C2H */
    result += cpor[30]*y[30]/31.0343; /*CH2OH */
    result += cpor[31]*y[31]/32.0423; /*CH3OH */
    result += cpor[32]*y[32]/40.0649; /*C3H4 */
    result += cpor[33]*y[33]/39.0569; /*C3H3 */
    result += cpor[34]*y[34]/41.0729; /*C3H5 */
    result += cpor[35]*y[35]/42.0808; /*C3H6 */
    result += cpor[36]*y[36]/44.0968; /*C3H8 */
    result += cpor[37]*y[37]/43.0888; /*I-C3H7 */
    result += cpor[38]*y[38]/43.0888; /*N-C3H7 */

    *cpbs = result * 8.314e+07;
}

static
void
Convert (const FArrayBox&   T_old,
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
            FORT_CONVERT_DRM19_TO_C3H8(&Told2, Xnew2.dataPtr());
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
    
    string tranfile="tran.asc.propane"; pp.query("tranfile",tranfile);
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

        for (MFIter mfi(molefrac); mfi.isValid(); ++mfi)
        {
            Convert(mfT[mfi], TComp_old, molefrac[mfi], xComp_old, mfout[lev][mfi],
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
