TOP             = ${HOME}/src/CCSE
BOXLIB_HOME     = ${TOP}/BoxLib
AMRVIS_HOME     = ${BOXLIB_HOME}/Src/Extern/amrdata
COMBUSTION_HOME = ${TOP}/Combustion
CHEMISTRY_DIR   = ${COMBUSTION_HOME}/Chemistry

PRECISION      = DOUBLE
DEBUG	       = FALSE
PROFILE        = FALSE
DIM    	       = 2
DIM    	       = 3
COMP           = g++
FCOMP          = gfortran
USE_MPI        = FALSE
USE_MPI        = TRUE
NEEDS_CHEM     = TRUE
EBASE          = replaceXwithC
EBASE          = pmfTest
USE_SDC        = TRUE

CHEMISTRY_MODEL=CHEMH
CHEMISTRY_MODEL=PROPANE
CHEMISTRY_MODEL=INERT30
CHEMISTRY_MODEL=CH4-2STEP
CHEMISTRY_MODEL=DRM19
CHEMISTRY_MODEL=LUDME
CHEMISTRY_MODEL=LIDRY
CHEMISTRY_MODEL=WANGDODECANE
CHEMISTRY_MODEL=GRI30NON
CHEMISTRY_MODEL=GRI12
CHEMISTRY_MODEL=GLAR
CHEMISTRY_MODEL=GRI30
CHEMISTRY_MODEL=ALZETA

#CXXFLAGS += -fno-inline -ggdb 
CFLAGS +="-std=c99"

include ${BOXLIB_HOME}/Tools/C_mk/Make.defs 

fincludes=${includes}

# Chemistry
ifeq (${NEEDS_CHEM}, TRUE)
  Bdirs += ${CHEMISTRY_DIR}/src
endif

# BoxLib 
Bdirs   += ${BOXLIB_HOME}/Src/C_BaseLib
Bdirs   += ${BOXLIB_HOME}/Src/C_BoundaryLib

DEFINES += -DBL_PARALLEL_IO -DBL_NOLINEVALUES
Bdirs   += ${AMRVIS_HOME}

Bpack	+= $(foreach dir, $(Bdirs), $(dir)/Make.package)
Blocs	+= $(foreach dir, $(Bdirs), $(dir))

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs) ${BOXLIB_HOME}/Src/C_AMRLib
VPATH_LOCATIONS   += $(Blocs) ${BOXLIB_HOME}/Src/C_AMRLib

CEXE_sources += ${EBASE}.cpp
ifeq ($(NEEDS_FORT), TRUE)
  FEXE_sources += ${EBASE}_F.F
endif
FEXE_sources += FILCC_$(DIM)D.F
CEXE_sources += AppendToPlotFile.cpp
CEXE_headers += AppendToPlotFile.H

ifeq (${USE_SDC}, TRUE)
  XTRADEFS += -DLMC_SDC
endif

ifeq (${NEEDS_CHEM}, TRUE)
  ifeq (${CHEMISTRY_MODEL}, DRM19)
    cEXE_sources += drm19.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/gri
  endif
  ifeq (${CHEMISTRY_MODEL}, CHEMH)
    cEXE_sources += chem-H.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/chem-H
  endif
  ifeq (${CHEMISTRY_MODEL}, GRI30)
    cEXE_sources += grimech30.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/gri
  endif
  ifeq (${CHEMISTRY_MODEL}, GRI12)
    cEXE_sources += grimech12.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/gri
  endif
  ifeq (${CHEMISTRY_MODEL}, LIDRY)
    cEXE_sources += LiDryer.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/LiDryer
  endif
  ifeq (${CHEMISTRY_MODEL}, LIDRYMOD)
    cEXE_sources += LiDryerMOD.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/LiDryer
  endif
  ifeq (${CHEMISTRY_MODEL}, LUDME)
    cEXE_sources += LuDME.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/Lu
  endif
  ifeq (${CHEMISTRY_MODEL}, GLARSKEL)
    cEXE_sources += glarSkel.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/glar
  endif
  ifeq (${CHEMISTRY_MODEL}, GLAR)
    cEXE_sources += glar.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/glar
  endif
  ifeq (${CHEMISTRY_MODEL}, GRI30NON)
    cEXE_sources += grimech30-noArN.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/gri
  endif
  ifeq (${CHEMISTRY_MODEL}, WANGDODECANE)
    cEXE_sources += dodecane_wang.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/dodecane_wang
  endif
  ifeq (${CHEMISTRY_MODEL}, ALZETA)
    cEXE_sources += alzeta.c
    MODEL_DIR = ${CHEMISTRY_DIR}/data/Alzeta
  endif


  VPATH_LOCATIONS += ${MODEL_DIR}:${MODEL_DIR}/PMFs

  # Include mech name in exec
  chemOptionsSuffix = $(addprefix $(optionsSuffix), .$(CHEMISTRY_MODEL))
  executable        = $(addsuffix $(chemOptionsSuffix).ex, $(EBASE))

endif

vpath %.c   $(VPATH_LOCATIONS)
vpath %.cpp $(VPATH_LOCATIONS)
vpath %.h   $(VPATH_LOCATIONS)
vpath %.H   $(VPATH_LOCATIONS)
vpath %.F   $(VPATH_LOCATIONS)
vpath %.f   $(VPATH_LOCATIONS)
vpath %.f90 $(VPATH_LOCATIONS)

all: $(executable)

$(executable):

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)
