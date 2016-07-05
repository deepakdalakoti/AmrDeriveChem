# Microsoft Developer Studio Project File - Name="replaceXwithY3d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=replaceXwithY3d - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "replaceXwithY3d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "replaceXwithY3d.mak" CFG="replaceXwithY3d - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "replaceXwithY3d - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "replaceXwithY3d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /include:"Release/" /nologo /warn:nofileopt
# ADD F90 /assume:noaccuracy_sensitive /compile_only /debug:none /iface:nomixed_str_len_arg /iface:cref /include:"Release/" /math_library:fast /nologo /threads /tune:k7 /warn:nofileopt /unroll:4
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GR /GX /O2  /I "." /I "..\Parallel\BoxLib" /I "..\Parallel\LMC" /I "..\Parallel\pAmrvis" /I "..\Parallel\pAmrDerive" /I "C:\WMPI\include" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D "BL_SPACEDIM=3"  /D "BL_PRVERSION=5" /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_USE_CHEM" /D "BL_NOLINEVALUES" /D for="if(0);else for" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dformt.lib /nologo /subsystem:console /machine:I386 /include:"__matherr"

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /include:"Debug/" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /dbglibs /debug:full /iface:nomixed_str_len_arg /iface:cref /include:"Debug/" /libs:static /nologo /traceback /warn:argument_checking /optimize:0 /threads /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /MTd /W3 /Gm /GR /GX /ZI /Od  /I "." /I "..\Parallel\BoxLib" /I "..\Parallel\LMC" /I "..\Parallel\pAmrvis" /I "..\Parallel\pAmrDerive" /I "C:\WMPI\include" /D "_CONSOLE" /D "_MBCS" /D "_DEBUG" /D "WIN32" /D "BL_USE_DOUBLE" /D "BL_ARCH_IEEE" /D "BL_USE_NEW_HFILES" /D "BL_SPACEDIM=3"  /D "BL_PRVERSION=5" /D "BL_PARALLEL_IO" /D "BL_FORT_USE_UPPERCASE" /D "BL_LANG_CC" /D "BL_USE_CHEM" /D "BL_NOLINEVALUES" /D for="if(0);else for" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib dformt.lib /nologo /subsystem:console /debug /machine:I386 /include:"__matherr" /pdbtype:sept

!ENDIF 

# Begin Target

# Name "replaceXwithY3d - Win32 Release"
# Name "replaceXwithY3d - Win32 Debug"
# Begin Group "C++ Sources"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxLib.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\ParmParse.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Utility.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\UseCount.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\DistributionMapping.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\ParallelDescriptor.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\VisMF.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Arena.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BArena.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\CArena.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLProfiler.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLThread.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLWorkQueue.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FabConv.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FPC.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Box.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\IntVect.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\IndexType.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Orientation.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxList.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxArray.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxDomain.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FArrayBox.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BaseFab.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\MultiFab.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FabArray.cpp
# End Source File
# Begin Source File

SOURCE=replaceXwithY.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\chem-H.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\ChemDriver.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\DataServices.cpp
# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\AmrData.cpp
# End Source File
# End Group
# Begin Group "C++ Headers"

# PROP Default_Filter "H"
# Begin Source File

SOURCE=..\Parallel\BoxLib\SPECIALIZE_F.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\ChemDriver_F.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\cdwrk.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxLib.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLVERSION.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\ParmParse.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Utility.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLassert.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\ArrayLim.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\REAL.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\CONSTANTS.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\SPACE.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\SPACE_F.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\UseCount.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\DistributionMapping.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\ParallelDescriptor.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\VisMF.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Arena.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BArena.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\CArena.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Profiler.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLFort.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Thread.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\WorkQueue.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FabConv.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FPC.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Box.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\IntVect.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\IndexType.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Orientation.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxList.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxArray.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BoxDomain.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\FArrayBox.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\Looping.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\MultiFab.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\ChemDriver.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\DataServices.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\AmrData.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\XYPlotDataList.H
# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\AmrvisConstants.H
# End Source File
# End Group
# Begin Group "Fortran"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLutil_F.f
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLParmParse_F.f
# End Source File
# Begin Source File

SOURCE=..\Parallel\BoxLib\BLBoxLib_F.f
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\vode.f
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\EGSlib.f
# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\EGini.f
# End Source File
# Begin Group "TEMPS-Debug"

# PROP Default_Filter "FOR"
# Begin Source File

SOURCE=Debug\Fort\SPECIALIZE_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Debug\Fort\ChemDriver_F.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Debug\Fort\ChemDriver_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Debug\Fort\DERIVE_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Debug\Fort\FABUTIL_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

!ENDIF 
# End Source File
# End Group
# Begin Group "TEMPS-Release"

# PROP Default_Filter "FOR"
# Begin Source File

SOURCE=Release\Fort\SPECIALIZE_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Release\Fort\ChemDriver_F.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Release\Fort\ChemDriver_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Release\Fort\DERIVE_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 
# End Source File
# Begin Source File

SOURCE=Release\Fort\FABUTIL_3D.FOR

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 
# End Source File
# End Group
# Begin Source File

SOURCE=..\Parallel\BoxLib\SPECIALIZE_3D.F

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\BoxLib\SPECIALIZE_3D.F
InputName=SPECIALIZE_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\BoxLib\SPECIALIZE_3D.F
InputName=SPECIALIZE_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF

# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\ChemDriver_F.F

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\LMC\ChemDriver_F.F
InputName=ChemDriver_F

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\LMC\ChemDriver_F.F
InputName=ChemDriver_F

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF

# End Source File
# Begin Source File

SOURCE=..\Parallel\LMC\ChemDriver_3D.F

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\LMC\ChemDriver_3D.F
InputName=ChemDriver_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\LMC\ChemDriver_3D.F
InputName=ChemDriver_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF

# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrDerive\DERIVE_3D.F

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\pAmrDerive\DERIVE_3D.F
InputName=DERIVE_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\pAmrDerive\DERIVE_3D.F
InputName=DERIVE_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF

# End Source File
# Begin Source File

SOURCE=..\Parallel\pAmrvis\FABUTIL_3D.F

!IF  "$(CFG)" == "replaceXwithY3d - Win32 Release"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\pAmrvis\FABUTIL_3D.F
InputName=FABUTIL_3D

"Release\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Release\Fort\$(InputName).FOR

# End Custom Build

!ELSEIF  "$(CFG)" == "replaceXwithY3d - Win32 Debug"

# PROP Ignore_Default_Tool 1
# Begin Custom Build
InputPath=..\Parallel\pAmrvis\FABUTIL_3D.F
InputName=FABUTIL_3D

"Debug\Fort\$(InputName).FOR" : $(SOURCE) "$(INTDIR)" "$(OUTDIR)"
	fpp /m /ansi /nologo  /S. /S..\Parallel\BoxLib /S..\Parallel\LMC /S..\Parallel\pAmrvis /S..\Parallel\pAmrDerive /DBL_LANG_FORT /DBL_NOLINEVALUES /DBL_SPACEDIM=3 /DBL_USE_DOUBLE /DBL_NO_FORT_FLUSH /DBL_USE_CHEM  /DBL_PRVERSION=5 $(InputPath) | perl ..\Parallel\scripts\strip72 -c > Debug\Fort\$(InputName).FOR

# End Custom Build

!ENDIF

# End Source File
# End Group
# End Target
# End Project
