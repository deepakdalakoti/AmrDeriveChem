
#include "REAL.H"
#include "Box.H"
#include "PArray.H"
#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "ChemDriver.H"

int
main (int   argc,
      char* argv[])
{
  BoxLib::Initialize(argc,argv);
    
  ParmParse pp;

  FArrayBox::setFormat(FABio::FAB_IEEE_32);
    
  int verbose = 0; pp.query("verbose",verbose);
  if (verbose > 2)
    AmrData::SetVerbose(true);

  std::string infile; pp.get("infile",infile);

  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(infile, fileType);
  if (!dataServices.AmrDataOk())
    DataServices::Dispatch(DataServices::ExitRequest, NULL);    
  AmrData& amrData = dataServices.AmrDataRef();
  int finestLevel = amrData.FinestLevel(); pp.query("finestLevel",finestLevel);

  Array<std::string> varNames;
  int nc = pp.countval("comps");
  pp.queryarr("comps",varNames,0,nc);

  Array<int> comps(varNames.size());
  for (int i=0; i<comps.size(); ++i) {
    comps[i] = amrData.StateNumber(varNames[i]);
  }

  Array<int> destFillComps(varNames.size());
  for (int i=0; i<varNames.size(); ++i)
    destFillComps[i] = i;

  if (verbose && ParallelDescriptor::IOProcessor())
    for (int i=0; i<varNames.size(); ++i)
      std::cout << "Getting component: " << varNames[i] << std::endl;

  const Array<Real>& plo = amrData.ProbLo();
  Array<Array<Real> > vals(nc+BL_SPACEDIM,Array<Real>());
  long cnt = 0;
  for (int iLevel=finestLevel; iLevel>=0; --iLevel)
  {
    int cr = 1;
    const BoxArray& ba = amrData.boxArray(iLevel);
    MultiFab mf(ba,nc,0,Fab_allocate);
    amrData.FillVar(mf,iLevel,varNames,destFillComps);
    const Array<Real>& dx = amrData.DxLevel()[iLevel];

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
      const FArrayBox& fab = mf[mfi];
      BoxArray ba(mfi.validbox());

      if (iLevel < finestLevel) {
        cr = amrData.RefRatio()[iLevel];
        const BoxArray baf = BoxArray(amrData.boxArray(iLevel+1)).coarsen(cr);
        ba = BoxLib::complementIn(mfi.validbox(),baf);
      }
      cnt = vals[0].size();
      for (int n=0; n<nc+BL_SPACEDIM; ++n) {
        vals[n].resize(vals[n].size() + ba.numPts());
      }
      for (int i=0; i<ba.size(); ++i) {
        const Box& box = ba[i];
        for (IntVect iv=box.smallEnd(), bend=box.bigEnd(); iv<=bend; box.next(iv)) {
          for (int d=0; d<BL_SPACEDIM; ++d) {
            vals[d][cnt] = plo[d] + (iv[d] + 0.5)*dx[d];
          }
          for (int n=0; n<nc; ++n) {
            vals[n+BL_SPACEDIM][cnt] = fab(iv,n);
          }
          cnt++;
        }
      }
    }
  }

  long tot = cnt;
  ParallelDescriptor::ReduceLongSum(tot);
  std::ofstream osf;
  std::string outfile="outfile.dat"; pp.query("outfile",outfile);
  for (int i=0; i<ParallelDescriptor::NProcs(); ++i) {
    if (i==ParallelDescriptor::MyProc()) {
      std::ios_base::openmode mode = (i==0 ? std::ios_base::trunc : std::ios_base::app);
      osf.open(outfile.c_str(),mode);
      std::cout << "Proc " << i << " writing " << cnt << " lines..." << std::endl;
      for (int j=0; j<cnt; ++j) {
        for (int n=0; n<nc+BL_SPACEDIM; ++n) {
          osf << vals[n][j] << " ";
        }
        osf << '\n';
      }
      osf.close();
    }
    ParallelDescriptor::Barrier();
  }
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Wrote a total of " << tot << " lines." << std::endl;
  }
  BoxLib::Finalize();
  return 0;
}
