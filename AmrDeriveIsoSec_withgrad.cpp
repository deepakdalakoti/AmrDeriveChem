
#include <winstd.H>
#ifdef WIN32
#pragma warning(disable:4503)
#endif

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <list>
#include <map>
#include <set>
#include <iomanip>
#include<vector>
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::set;
using std::pair;
using std::map;
#ifndef WIN32
#include <unistd.h>
#endif

#include "ParmParse.H"
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Utility.H"
#include "FArrayBox.H"
#include "Geometry.H"
#include "AmrDeriveIsoSec_withgrad_F.H"
#include "ChemDriver.H"
#include "ChemDriver_F.H"
#include "Array.H"
#include "coeff.h"
static Real isoVal_DEF = 1090.;
static string isoCompName_DEF = "temp";
static Real epsilon_DEF = 1.e-12;

extern "C"
{
    void FORT_TRICUBIC_GET_COEFF_STACKED(Real a[64], const Real x[64])
    {
        int i,j;
        for (i=0;i<64;i++) {
            a[i]=(Real)(0.0);
            for (j=0;j<64;j++) {
                a[i]+=A[i][j]*x[j];
            }
        }
    }

    void FORT_TRICUBIC_GET_COEFF(Real       a[64],
                                 const Real f[8],
                                 const Real dfdx[8],
                                 const Real dfdy[8],
                                 const Real dfdz[8],
                                 const Real d2fdxdy[8],
                                 const Real d2fdxdz[8],
                                 const Real d2fdydz[8],
                                 const Real d3fdxdydz[8])
    {
        int i;
        Real x[64];
        for (i=0;i<8;i++) {
            x[0+i]=f[i];
            x[8+i]=dfdx[i];
            x[16+i]=dfdy[i];
            x[24+i]=dfdz[i];
            x[32+i]=d2fdxdy[i];
            x[40+i]=d2fdxdz[i];
            x[48+i]=d2fdydz[i];
            x[56+i]=d3fdxdydz[i];
        }
        FORT_TRICUBIC_GET_COEFF_STACKED(a,x);
    }

    int ijk2n(int i, int j, int k) {
        return(i+4*j+16*k);
    }

    void point2xyz(int p, int *x, int *y, int *z) {
        switch (p) {
        case 0: *x=0; *y=0; *z=0; break;
        case 1: *x=1; *y=0; *z=0; break;
        case 2: *x=0; *y=1; *z=0; break;
        case 3: *x=1; *y=1; *z=0; break;
        case 4: *x=0; *y=0; *z=1; break;
        case 5: *x=1; *y=0; *z=1; break;
        case 6: *x=0; *y=1; *z=1; break;
        case 7: *x=1; *y=1; *z=1; break;
        default:*x=0; *y=0; *z=0;
        }
    }
    void FORT_TRICUBIC_EVAL(const Real  a[64],
                            const Real* xp,
                            const Real* yp,
                            const Real* zp,
                            Real*       eval)
    {
        int i,j,k;
        Real ret=0.;
        Real x=*xp, y=*yp, z=*zp;

        Real xpow[4] = {1, x, x*x, x*x*x};
        Real ypow[4] = {1, y, y*y, y*y*y};
        Real zpow[4] = {1, z, z*z, z*z*z};

        /* TRICUBIC EVAL
           This is the short version of tricubic_eval. It is used to compute
           the value of the function at a given point (x,y,z). To compute
           partial derivatives of f, use the full version with the extra args.
        */
        for (i=0;i<4;i++) {
            for (j=0;j<4;j++) {
                for (k=0;k<4;k++) {
		  ret+=a[ijk2n(i,j,k)]*xpow[i]*ypow[j]*zpow[k];
		  //ret+=a[ijk2n(i,j,k)]*pow(x,i)*pow(y,j)*pow(z,k);
                }
            }
        }
        *eval = ret;
    }
}

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " inputs infile=<s> isoCompName=<s> isoVal=<v> [options] \n\tOptions:\n";
    std::cerr << "\t     infile=<s> where <s> is a pltfile\n";
    std::cerr << "\t     isoCompName=<s> where <s> is the quantity being contoured\n";
    std::cerr << "\t     isoVal=<v> where <v> is an isopleth value\n";
    std::cerr << "\t     Choosing quantities to interp to surface: \n";
    std::cerr << "\t       comps=int comp list [overrides sComp/nComp]\n";
    std::cerr << "\t       sComp=start comp[DEF->0]\n";
    std::cerr << "\t       nComp=number of comps[DEF->all]\n";
    std::cerr << "\t     finestLevel=<n> finest level to use in pltfile[DEF->all]\n";
    std::cerr << "\t     writeSurf=<1,0> output surface in binary MEF format [DEF->1]\n";
    std::cerr << "\t     outfile=<s> name of tecplot output file [DEF->gen'd]\n";
    std::cerr << "\t  NOTE: The remaining options are as yet untested\n";
    std::cerr << "\t     box=int list of LL+UR of subbox of source data (lev-0 coords) [DEF->all]>\n";
    std::cerr << "\t     bounds=real list of LL+UR of subbox of final surface [DEF->all]>\n";
    std::cerr << "\t     nGrowPer=<#> number of coarse-grid cells to tack onto periodic bndry [DEF->0]>\n";
    std::cerr << "\t     finestLevel=<#> finest level to use [DEF->pltfile finest]>\n";
    exit(1);
}

struct Edge
{
    Edge(const IntVect& lhs, const IntVect& rhs) {
        if (IntVect::Compare()(lhs,rhs)) {
            IV_l = lhs;
            IV_r = rhs;
        } else {
            IV_l = rhs;
            IV_r = lhs;
        }
    }

    bool operator== (const Edge& rhs) const
        {
            return IV_l==rhs.IV_l && IV_r==rhs.IV_r;
        }
    bool operator< (const Edge& rhs) const
        {
        if (IV_l==rhs.IV_l) {
            return IntVect::Compare()(IV_r,rhs.IV_r);
        } else {
            return IntVect::Compare()(IV_l,rhs.IV_l);
        }
    }

    IntVect IV_l,IV_r;
};

class IVLess
{
public:
    bool operator () (const IntVect& lhs,
                      const IntVect& rhs) const
        {
            return lhs.lexLT(rhs);
        }
};

typedef Array<Real> Point;

typedef std::map<Edge, Point> PMap;

typedef PMap::iterator PMapIt;

struct Segment
{
    Segment() : p(2), mLength(-1) {}
    Real Length ();
    const PMapIt& operator[] (int n) const { return p[n]; }
    PMapIt& operator[] (int n) { return p[n]; }
    void flip ();
    static int xComp,yComp;
    Array<PMapIt> p;
    
private:
    void my_length();
    Real mLength;
};

int Segment::xComp = 0;
int Segment::yComp = 1;

Real
Segment::Length()
{
    if (mLength<0)
        my_length();
    return mLength;
}

void
Segment::my_length()
{
#ifndef NDEBUG
    BL_ASSERT(xComp>=0 && yComp>=0);
    for (int i=0; i<p.size(); ++i)
    {
        BL_ASSERT(xComp<(*p[i]).second.size());
        BL_ASSERT(yComp<(*p[i]).second.size());
    }
#endif
    const Point& p0 = (*p[0]).second;
    const Point& p1 = (*p[1]).second;
    mLength = std::sqrt(((p1[xComp] - p0[xComp])*(p1[xComp] - p0[xComp]))
                        +((p1[yComp] - p0[yComp])*(p1[yComp] - p0[yComp])));
}    

void
Segment::flip()
{
    PMapIt ptmp = p[0];
    p[0] = p[1];
    p[1] = ptmp;
}

std::ostream&
operator<< (std::ostream& os, const Segment& seg)
{
    const Point& p0 = seg.p[0]->second;
    const Point& p1 = seg.p[1]->second;

    os << '[' << p0[0];
    for (int i = 1; i < BL_SPACEDIM; i++)
    {
        os << ", " << p0[i];
    }
    os << "] ... ";

    os << '[' << p1[0];
    for (int i = 1; i < BL_SPACEDIM; i++)
    {
        os << ", " << p1[i];
    }
    os << "]";

    return os;
}

typedef list<Segment> SegList;

struct PMapItCompare
{
    bool operator() (const PMapIt& lhs, const PMapIt& rhs) const {
        const Point& l = lhs->second;
        const Point& r = rhs->second;
        const void* vl=&l;
        const void* vr=&r;
        return vl<vr;
    }
};

struct Triangle
{
    Triangle() : p(3), mArea(-1) {}
    Real Area();
    const PMapIt& operator[] (int n) const { return p[n]; }
    PMapIt& operator[] (int n) { return p[n]; }
    static int xComp,yComp,zComp;
    Array<PMapIt> p;
    
private:
    void my_area();
    Real mArea;
};

int Triangle::xComp = 0;
int Triangle::yComp = 1;
int Triangle::zComp = 2;

Real
Triangle::Area()
{
    if (mArea<0)
        my_area();
    return mArea;
}

void
Triangle::my_area()
{
#ifndef NDEBUG
    for (int i=0; i<p.size(); ++i)
    {
        BL_ASSERT(xComp>=0 && xComp<(*p[i]).second.size());
        BL_ASSERT(yComp>=0 && yComp<(*p[i]).second.size());
        BL_ASSERT(zComp>=0 && zComp<(*p[i]).second.size());
    }
#endif

    const Point& p0 = (*p[0]).second;
    const Point& p1 = (*p[1]).second;
    const Point& p2 = (*p[2]).second;
    mArea = 0.5*sqrt(
        pow(  ( p1[yComp] - p0[yComp])*(p2[zComp]-p0[zComp]) 
              -(p1[zComp] - p0[zComp])*(p2[yComp]-p0[yComp]), 2)
        
        + pow(( p1[zComp] - p0[zComp])*(p2[xComp]-p0[xComp]) 
              -(p1[xComp] - p0[xComp])*(p2[zComp]-p0[zComp]), 2)
        
        + pow(( p1[xComp] - p0[xComp])*(p2[yComp]-p0[yComp]) 
              -(p1[yComp] - p0[yComp])*(p2[xComp]-p0[xComp]), 2));
}    

/*
   Linearly interpolate the position where an isosurface cuts
   an edge between two vertices, each with their own scalar value
*/
Point
VI_doIt(Real isoVal, int isoComp, const Point& p1, const Point& p2)
{
    BL_ASSERT(p1.size()>isoComp && p2.size()>isoComp);

    const Real valp1 = p1[isoComp];
    const Real valp2 = p2[isoComp];
    
    if (std::abs(isoVal-valp1) < epsilon_DEF)
        return p1;
    if (std::abs(isoVal-valp2) < epsilon_DEF)
        return p2;
    if (std::abs(valp1-valp2) < epsilon_DEF)
        return p1;
    
    Point res(p1.size());
    const Real mu = (isoVal - valp1) / (valp2 - valp1);
    for (int j=0; j<res.size(); ++j)
        res[j] = p1[j] + mu * (p2[j] - p1[j]);
    
    return res;
}

PMapIt VertexInterp(Real isoVal, int isoComp,
                    const IntVect& p1,const Point& p1d,
                    const IntVect& p2,const Point& p2d,
                    PMap& vertCache)
{
    PMapIt fwd,rev;
    Edge edge(p1,p2);
    fwd = vertCache.find(edge);
    if (fwd == vertCache.end())
    {
        rev = vertCache.find(Edge(p2,p1));
        
        if (rev == vertCache.end())
        {
            std::pair<Edge,Point> ent(edge,VI_doIt(isoVal,isoComp,p1d,p2d));
            std::pair<PMapIt,bool> it = vertCache.insert(ent);
            BL_ASSERT(it.second);
            return it.first;
        }
        else
        {
            return rev;
        }
    }
    return fwd;
}

#if BL_SPACEDIM==2
/*
   Given a grid cell and an isoVal, calculate the line segments
   required to represent the contour through the cell.
   Return an array of (at most 2) line segments
*/

Array<Segment> Segmentise(const FArrayBox& pts,
                          const FArrayBox& mask,
                          PMap&            vertCache,
                          const IntVect&   baseIV,
                          Real             isoVal,
                          int              isoComp)
{
   Array<PMapIt> vertlist(4);
   Array<Segment> segments;

   const IntVect& p0 = baseIV;
   const IntVect  p1 = p0 + BoxLib::BASISV(0);
   const IntVect  p2 = p1 + BoxLib::BASISV(1);
   const IntVect  p3 = p0 + BoxLib::BASISV(1);
   //
   // Bail if any of the point are masked out.
   //
   if (mask(p0)<0 || mask(p1)<0 || mask(p2)<0 || mask(p3)<0)
       return segments;

   const int nComp = pts.nComp();

   Point p0d(nComp); pts.getVal(p0d.dataPtr(),p0);
   Point p1d(nComp); pts.getVal(p1d.dataPtr(),p1);
   Point p2d(nComp); pts.getVal(p2d.dataPtr(),p2);
   Point p3d(nComp); pts.getVal(p3d.dataPtr(),p3);

   int segCase = 0;

   if (p0d[isoComp] < isoVal) segCase |= 1;
   if (p1d[isoComp] < isoVal) segCase |= 2;
   if (p2d[isoComp] < isoVal) segCase |= 4;
   if (p3d[isoComp] < isoVal) segCase |= 8;

   if (segCase==0 || segCase==15)
       return segments;

   if (segCase==5 || segCase==10)
   {
       segments.resize(2);
   } else
   {
       segments.resize(1);
   }

   switch (segCase)
   {
   case 1:
   case 14:
       vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
       break;
   case 2:
   case 13:
       vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
       break;
   case 3:
   case 12:
       vertlist[0] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
       break;
   case 4:
   case 11:
       vertlist[0] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
       break;
   case 6:
   case 9:
       vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
       break;
   case 7:
   case 8:
       vertlist[0] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
       break;
   case 5:
   case 10:
       vertlist[0] = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
       vertlist[1] = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
       vertlist[2] = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
       vertlist[3] = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
       break;
   }

   segments[0][0] = vertlist[0];
   segments[0][1] = vertlist[1];
   
   if (segCase==5 || segCase==10)
   {
       segments[1][0] = vertlist[2];
       segments[1][1] = vertlist[3];
   }

   return segments;
}

#else

/*
   Given a grid cell and an isoVal, calculate the triangular
   facets required to represent the isosurface through the cell.
   Return an array of (at most 5) triangular facets
*/
Array<Triangle> Polygonise(const FArrayBox& pts,
                           const FArrayBox& mask,
                           PMap&            vertCache,
                           const IntVect&   baseIV,
                           Real             isoVal,
                           int              isoComp)
{
   int cubeindex;
   Array<PMapIt> vertlist(12);
   Array<Triangle> triangles; // result

   const IntVect& p0 = baseIV;
   const IntVect p1 = p0 + BoxLib::BASISV(0);
   const IntVect p2 = p1 + BoxLib::BASISV(1);
   const IntVect p3 = p0 + BoxLib::BASISV(1);
   const IntVect p4 = p0 + BoxLib::BASISV(2);
   const IntVect p5 = p4 + BoxLib::BASISV(0);
   const IntVect p6 = p5 + BoxLib::BASISV(1);
   const IntVect p7 = p4 + BoxLib::BASISV(1);

   // Bail if any of the point are masked out
   if (mask(p0)<0 || mask(p1)<0 || mask(p2)<0 || mask(p3)<0 ||
       mask(p4)<0 || mask(p5)<0 || mask(p6)<0 || mask(p7)<0 )
       return triangles;

   const int nComp = pts.nComp();
   Point p0d(nComp); pts.getVal(p0d.dataPtr(),p0);
   Point p1d(nComp); pts.getVal(p1d.dataPtr(),p1);
   Point p2d(nComp); pts.getVal(p2d.dataPtr(),p2);
   Point p3d(nComp); pts.getVal(p3d.dataPtr(),p3);
   Point p4d(nComp); pts.getVal(p4d.dataPtr(),p4);
   Point p5d(nComp); pts.getVal(p5d.dataPtr(),p5);
   Point p6d(nComp); pts.getVal(p6d.dataPtr(),p6);
   Point p7d(nComp); pts.getVal(p7d.dataPtr(),p7);

   int edgeTable[256]={
       0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
       0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
       0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
       0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
       0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
       0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
       0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
       0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
       0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
       0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
       0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
       0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
       0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
       0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
       0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
       0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
       0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
       0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
       0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
       0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
       0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
       0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
       0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
       0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
       0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
       0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
       0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
       0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
       0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
       0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
       0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
       0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

   int triTable[256][16] =
   {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};

   /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
   */
   cubeindex = 0;
   if (p0d[isoComp] < isoVal) cubeindex |= 1;
   if (p1d[isoComp] < isoVal) cubeindex |= 2;
   if (p2d[isoComp] < isoVal) cubeindex |= 4;
   if (p3d[isoComp] < isoVal) cubeindex |= 8;
   if (p4d[isoComp] < isoVal) cubeindex |= 16;
   if (p5d[isoComp] < isoVal) cubeindex |= 32;
   if (p6d[isoComp] < isoVal) cubeindex |= 64;
   if (p7d[isoComp] < isoVal) cubeindex |= 128;

   /* Cube is entirely in/out of the surface */
   if (edgeTable[cubeindex] == 0)
       return triangles;

   /* Find the vertices where the surface intersects the cube */
   if (edgeTable[cubeindex] & 1)
       vertlist[0]  = VertexInterp(isoVal,isoComp,p0,p0d,p1,p1d,vertCache);
   if (edgeTable[cubeindex] & 2)
       vertlist[1]  = VertexInterp(isoVal,isoComp,p1,p1d,p2,p2d,vertCache);
   if (edgeTable[cubeindex] & 4)
       vertlist[2]  = VertexInterp(isoVal,isoComp,p2,p2d,p3,p3d,vertCache);
   if (edgeTable[cubeindex] & 8)
       vertlist[3]  = VertexInterp(isoVal,isoComp,p3,p3d,p0,p0d,vertCache);
   if (edgeTable[cubeindex] & 16)
       vertlist[4]  = VertexInterp(isoVal,isoComp,p4,p4d,p5,p5d,vertCache);
   if (edgeTable[cubeindex] & 32)
       vertlist[5]  = VertexInterp(isoVal,isoComp,p5,p5d,p6,p6d,vertCache);
   if (edgeTable[cubeindex] & 64)
       vertlist[6]  = VertexInterp(isoVal,isoComp,p6,p6d,p7,p7d,vertCache);
   if (edgeTable[cubeindex] & 128)
       vertlist[7]  = VertexInterp(isoVal,isoComp,p7,p7d,p4,p4d,vertCache);
   if (edgeTable[cubeindex] & 256)
       vertlist[8]  = VertexInterp(isoVal,isoComp,p0,p0d,p4,p4d,vertCache);
   if (edgeTable[cubeindex] & 512)
       vertlist[9]  = VertexInterp(isoVal,isoComp,p1,p1d,p5,p5d,vertCache);
   if (edgeTable[cubeindex] & 1024)
       vertlist[10]  = VertexInterp(isoVal,isoComp,p2,p2d,p6,p6d,vertCache);
   if (edgeTable[cubeindex] & 2048)
       vertlist[11]  = VertexInterp(isoVal,isoComp,p3,p3d,p7,p7d,vertCache);

   /* Create the triangles */
   int nTriang = 0;
   for (int i=0;triTable[cubeindex][i]!=-1;i+=3)
       nTriang++;

   triangles.resize(nTriang);
   for (int j=0; j<nTriang; ++j)
   {
       int j3 = 3*j;
       triangles[j][0] = vertlist[triTable[cubeindex][j3  ]];
       triangles[j][1] = vertlist[triTable[cubeindex][j3+1]];
       triangles[j][2] = vertlist[triTable[cubeindex][j3+2]];
   }

   return triangles;
}
#endif

struct Node
{
    Node() { m_vec = 0; }
    Node(const vector<Real>& vec, long i) : m_idx(i)
        {
            m_size = vec.size();
            m_vec = new Real[m_size];
            for (int i=0; i<m_size; ++i)
                m_vec[i] = vec[i];
        }
    ~Node() { delete [] m_vec; }
    int m_size;
    int m_idx;
    Real* m_vec;
    Real operator[] (int n) const { return m_vec[n]; }
    static bool myLessVerbose; 
    inline bool operator< (const Node& rhs) const {
        Real sum = 0.0;
        for (int i=0; i<BL_SPACEDIM; ++i) {
            sum += (m_vec[i]-rhs[i])*(m_vec[i]-rhs[i]);
        }

        if (myLessVerbose)
        {
            cerr << "In LT" << endl;
            cerr << "  lhs: " << m_vec[0] << " " << m_vec[1] << " " << m_vec[2] << endl;
            cerr << "  rhs: " << rhs[0] << " " << rhs[1] << " " << rhs[2] << endl;
            if (std::sqrt(sum) < epsilon_DEF) 
            {
                cerr << "under eps tol" << endl;
            }
            else
            { 
                if (m_vec[0]==rhs[0]) {
                    if (m_vec[1]==rhs[1]) {
                        cerr << "          based on z" << endl;
                    } else {
                        cerr << "          based on y" << endl;
                    }
                } else {
                    cerr << "          based on x" << endl;
                }
            }
        }
        if (std::sqrt(sum) < epsilon_DEF) return false;

        if (m_vec[0]==rhs[0]) {
            if (m_vec[1]==rhs[1]) {
                return m_vec[2] < rhs[2];
            } else {
                return m_vec[1] < rhs[1];
            }
        } else {
            return m_vec[0] < rhs[0];
        }
    }
    // copy ctr
    Node(const Node& rhs)
        {
            m_size = rhs.m_size;
            m_idx = rhs.m_idx;
            m_vec = new Real[m_size];
            for (int i=0; i<m_size; ++i)
                m_vec[i] = rhs.m_vec[i];
        }
};
bool Node::myLessVerbose=false;

struct Element
{
    Element (const vector<int>& vec)
        {
        m_vec[0] = vec[0]; 
        m_vec[1] = vec[1]; 
#if BL_SPACEDIM==3
        m_vec[2] = vec[2];
#endif
        //
        // Rotate elements so that first one is the smallest.
        //
        int smallest = 0;
        for (int i=1; i<BL_SPACEDIM; ++i)
        {
            if (m_vec[smallest]>m_vec[i]) smallest=i;
        }
        std::rotate(&m_vec[0],&m_vec[smallest],&m_vec[0]+BL_SPACEDIM);
    }
    int size () const {return BL_SPACEDIM;}
    int m_vec[BL_SPACEDIM];
    int operator[] (int n) const
        {
        BL_ASSERT(n >= 0 && n < BL_SPACEDIM);
        return m_vec[n];
    }
    bool operator< (const Element& rhs) const
        {
        if (m_vec[0]==rhs[0])
        {
#if BL_SPACEDIM==3
            if (m_vec[1]==rhs[1])
            {
                return m_vec[2] < rhs[2];
            }
            else
#endif
            {
                return m_vec[1] < rhs[1];
            }
        } else
        {
            return m_vec[0] < rhs[0];
        }
    }
};

void
Collate(Array<Real>& NodeRaw,
        Array<int>&  EltRaw,
        int          nCompPerNode)
{
#if BL_USE_MPI
    const int nProcs = ParallelDescriptor::NProcs(); 

    if (nProcs < 2) return;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    BL_ASSERT(IOProc==0);

    Array<int> nmdataR(nProcs,0);
    Array<int> offsetR(nProcs,0);
    //
    // Tell root CPU how many Real data elements each CPU will be sending.
    //
    int countR = NodeRaw.size();
    MPI_Gather(&countR,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmdataR.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());

    Array<Real> NodeRawT;
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1; i < nProcs; i++)
            offsetR[i] = offsetR[i-1] + nmdataR[i-1];

        NodeRawT.resize(offsetR[nProcs-1] + nmdataR[nProcs-1]);
    }
    //
    // Gather all the Real data to IOProc into NodeRaw.
    //
    MPI_Gatherv(NodeRaw.dataPtr(),
                countR,
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                NodeRawT.dataPtr(),
                nmdataR.dataPtr(),
                offsetR.dataPtr(),
                ParallelDescriptor::Mpi_typemap<Real>::type(),
                IOProc,
                ParallelDescriptor::Communicator());

    // Now communicate the element info
    Array<int> nmdataI(nProcs,0);
    Array<int> offsetI(nProcs,0);
    int countI = EltRaw.size();
    MPI_Gather(&countI,
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               nmdataI.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<int>::type(),
               IOProc,
               ParallelDescriptor::Communicator());

    Array<int> EltRawT;
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1; i < nProcs; i++)
            offsetI[i] = offsetI[i-1] + nmdataI[i-1];

        EltRawT.resize(offsetI[nProcs-1] + nmdataI[nProcs-1]);
    }
    //
    // Gather all the data to IOProc into EltRaw
    //
    MPI_Gatherv(EltRaw.dataPtr(),
                countI,
                ParallelDescriptor::Mpi_typemap<int>::type(),
                EltRawT.dataPtr(),
                nmdataI.dataPtr(),
                offsetI.dataPtr(),
                ParallelDescriptor::Mpi_typemap<int>::type(),
                IOProc,
                ParallelDescriptor::Communicator());
    //
    // Shift nodeIDs in element definitions
    //
    if (ParallelDescriptor::IOProcessor())
    {
        for (int i = 1; i < nProcs; i++)
        {
            const int nodeOffset = offsetR[i]/nCompPerNode;
            for (int j = 0; j < nmdataI[i]; j++)                
                EltRawT[offsetI[i]+j] += nodeOffset;
        }
    }

    std::swap(EltRawT,EltRaw);
    std::swap(NodeRawT,NodeRaw);
#endif
}

class FABdata
{
public:
    FABdata(size_t i, int n)
        {fab.resize(Box(IntVect::TheZeroVector(),
                        IntVect(D_DECL(i-1,0,0))),n);}
    Real& operator[](size_t i) {return fab.dataPtr()[i];}
    FArrayBox fab;
    size_t boxSize;
};


void IsoBase( const Real*  pts,
                         const Real* theta,
                         const int nNodes, const Real dtheta,
                         const int nBins, Real* thetaAvg,const int nComp,Real* Base)
{

     
     Array<Real> Min(nBins,1.0);
     Array<int> Count(nBins,0);
  for (int i=0 ; i < nNodes ; i++) {

     int k = (int)((theta[i]-0)/dtheta);
      if(k>nBins-1) k = nBins-1;
     thetaAvg[k] = thetaAvg[k]+theta[i];

     Count[k] = Count[k] +1;

     if(pts[i*nComp+2]<Min[k]) {
       Min[k] = pts[i*nComp+2];
   for(int ip=0; ip<nComp;ip++)  Base[k*nComp+ip] = pts[i*nComp+ip];
       }
  }
 for (int i=0 ; i< nBins; i++) {
     if(Count[i]!=0)
     thetaAvg[i] = thetaAvg[i]/Count[i];
   
          }
}

////////////////////////////////////////////////////////////////////////////


typedef pair<int,int> Edge2D;
typedef map<Edge2D, pair<Point,int> > PMap2D;
typedef PMap2D::iterator PMapIt2D;
static PMap2D vertCache2D;



struct Segment2D
{
  Segment2D() : p(2), mLength(-1) {}
  Real Length ();
  const PMapIt2D& operator[] (int n) const { return p[n]; }
  PMapIt2D& operator[] (int n) { return p[n]; }
  int ID_l () const {return (*p[0]).second.second;}
  int ID_r () const {return (*p[1]).second.second;}
  void flip ();
  static int xComp,yComp,zComp;
  Array<PMapIt2D> p;

private:
  void my_length();
  Real mLength;
};

typedef list<Segment2D> SegList2D;
SegList2D::iterator FindMySeg2D(SegList2D& segs, int idx, vector<PMapIt2D>& vertVec)
{
    const PMapIt2D& vertToFind = vertVec[idx];
    for (SegList2D::iterator it=segs.begin(); it!=segs.end(); ++it)
    {
        if ( ((*it)[0] == vertToFind) || ((*it)[1] == vertToFind) )
            return it;
    }
    return segs.end();
}


Point VI_doIt2D(Real isoVal,int isoComp,const vector<Point>& pts,int p1,int p2)
{
    BL_ASSERT(p1<pts.size() && p2<pts.size());
    const Point& pt1 = pts[p1];
    const Point& pt2 = pts[p2];

    BL_ASSERT(isoComp!=pts.size() && isoComp<pt1.size() && pt1.size()==pt2.size());

    const Real valp1 = pt1[isoComp];
    const Real valp2 = pt2[isoComp];

    if (std::abs(isoVal-valp1) < epsilon_DEF)
        return pt1;
    if (std::abs(isoVal-valp2) < epsilon_DEF)
        return pt2;
    if (std::abs(valp1-valp2) < epsilon_DEF)
        return pt1;

    Point res(pt1.size());
//      Point res(pt1.size());
    const Real mu = (isoVal - valp1) / (valp2 - valp1);
    for (int j=0; j<res.size(); ++j)
        res[j] = pt1[j] + mu * (pt2[j] - pt1[j]);
//    res[pt1.size()] = std::pow(res[0]*res[0] +res[1]*res[1],0.5);
//    Real tmp = atan(std::abs(res[1]/(res[0]+1e-5)))*180/3.14;
//    res[pt1.size()+1] = tmp;
//    if (res[0]<0.0 && res[1] >0.0) res[pt1.size()+1] = 90+tmp;
//    if (res[0]<0.0 && res[1] <0.0) res[pt1.size()+1] = 180+tmp;
//    if (res[0]>0.0 && res[1] <0.0) res[pt1.size()+1] = 270+tmp;
    return res;
}




PMapIt2D VertexInterp2D(Real isoVal,int isoComp,const vector<Point>& pts,int p1,int p2)
{
    Edge2D ppair(p1,p2);
    PMapIt2D fwd,rev;
    fwd = vertCache2D.find(ppair);
    if (fwd == vertCache2D.end())
    {
        rev = vertCache2D.find(Edge2D(p2,p1));

        if (rev == vertCache2D.end())
        {
            return vertCache2D.insert(
                std::make_pair(ppair,
                std::make_pair(VI_doIt2D(isoVal,isoComp,pts,p1,p2),0))).first;
        }
        else
        {
            return rev;
        }
    }
    return fwd;
}


Array<Segment2D> Segmentise2D(const vector<Point>&  pts,
                          const vector<int>&    elt,
                          Real                  isoVal,
                          int                   isoComp)
{
   BL_ASSERT(elt.size()==3);
   Array<PMapIt2D> vertlist(3);

   const int p0 = elt[0];
   const int p1 = elt[1];
   const int p2 = elt[2];
   bool lo_0 = pts[p0][isoComp] < isoVal;
   bool lo_1 = pts[p1][isoComp] < isoVal;
   bool lo_2 = pts[p2][isoComp] < isoVal;

   int count = 0;
   if (lo_0 ^ lo_1) {
     vertlist[count++] = VertexInterp2D(isoVal,isoComp,pts,p0,p1);
   }
   if (lo_1 ^ lo_2) {
     vertlist[count++] = VertexInterp2D(isoVal,isoComp,pts,p1,p2);
   }
   if (lo_2 ^ lo_0) {
     vertlist[count++] = VertexInterp2D(isoVal,isoComp,pts,p2,p0);
   }

   Array<Segment2D> segments(0);
   if (count > 0) {
     BL_ASSERT(count == 2);
     segments.resize(1);
     segments[0][0] = vertlist[0];
     segments[0][1] = vertlist[1];
   }

   return segments;
}


int Segment2D::xComp = 0;
int Segment2D::yComp = 1;
int Segment2D::zComp = 2;

Real
Segment2D::Length()
{
    if (mLength<0)
        my_length();
    return mLength;
}

void
Segment2D::my_length()
{
  const Point& p0 = (*p[0]).second.first;
  const Point& p1 = (*p[1]).second.first;
  mLength = std::sqrt(((p1[xComp] - p0[xComp])*(p1[xComp] - p0[xComp]))
                      +((p1[yComp] - p0[yComp])*(p1[yComp] - p0[yComp]))
                      +((p1[zComp] - p0[zComp])*(p1[zComp] - p0[zComp])));
}

void
Segment2D::flip()
{
    PMapIt2D ptmp = p[0];
    p[0] = p[1];
    p[1] = ptmp;
}



bool operator==(const Edge2D& lhs, const Edge2D& rhs)
{
    return (lhs.first == rhs.first && lhs.second == rhs.second) ||
        (lhs.first == rhs.second && lhs.second == rhs.first);
}


bool operator==(const Segment2D& lhs, const Segment2D& rhs)
{
    return lhs.ID_l() == rhs.ID_l() && lhs.ID_r() == rhs.ID_r();
}

bool operator!=(const Segment2D& lhs, const Segment2D& rhs)
{
    return !operator==(lhs,rhs);
}

















/////////////////////////////////////////////////






















int
main (int   argc,
      char* argv[])
{
    BoxLib::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    int verbose=0; pp.query("verbose",verbose);
    if (verbose>1) AmrData::SetVerbose(true);
    int nPlotFiles = pp.countval("infile");
    Array<std::string> infile(nPlotFiles); 
    pp.getarr("infile",infile);

    Amrvis::FileType fileType(Amrvis::NEWPLT);
    Array<DataServices *> dataServicesPtrArray(nPlotFiles);                                         // DataServices array for each plot
    Array<AmrData *>      amrDataPtrArray(nPlotFiles);                                              // DataPtrArray for each plot
    Array<Real>           time(nPlotFiles);

  for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Loading " << infile[iPlot] << std::endl;

    dataServicesPtrArray[iPlot] = new DataServices(infile[iPlot], fileType);               // Populate DataServices array

    if( ! dataServicesPtrArray[iPlot]->AmrDataOk())                                               // Check AmrData ok
      DataServices::Dispatch(DataServices::ExitRequest, NULL);                                    // Exit if not

    amrDataPtrArray[iPlot] = &(dataServicesPtrArray[iPlot]->AmrDataRef());                        // Populate DataPtrArray

    time[iPlot] = amrDataPtrArray[iPlot]->Time();

    if (verbose) std::cout << "Time = " << time[iPlot] << std::endl;
  }
  int Proc2;
  
  int finestLevel = amrDataPtrArray[0]->FinestLevel()+1;
  
  pp.query("finestLevel",finestLevel);
  int Nlev = finestLevel+1;
   ChemDriver cd ; 
  int nSpec = cd.numSpecies();
  Array<Real> Datapt;
  Array<int>  ConPt;
  Array<Real > Out;
// std::cout << nSpec << std::endl;
   Array<Real> DatOut(nPlotFiles*8);
//  int nVars = pp.countval("varNames"); 
   int nVars = nSpec + 7;
//     int nVars =1 ;
  Array<std::string> varNames(nVars);
 int Nsp =0 ;
  int isoComp2D=0;
  pp.query("isoComp2D", isoComp2D);
  Real isoVal2D = 0 ;
  pp.query("isoVal2D", isoVal2D);
   std::string whichSpec="OH";
   pp.query("whichSpec",whichSpec);

// varNames[0] = "temp";
//pp.getarr("varNames", varNames);
  for (int i =0 ; i < nSpec; i++) {
       varNames[i] = "Y("+cd.speciesNames()[i]+")";
       if(cd.speciesNames()[i]==whichSpec) Nsp = i;
     }
   varNames[nSpec] = "density";
  varNames[nSpec+1] = "temp";
   varNames[nSpec+2] = "x_velocity";
   varNames[nSpec+3] = "y_velocity";
   varNames[nSpec+4] = "z_velocity";
   varNames[nSpec+5] = "Z";
   varNames[nSpec+6] = "I_R[rhoY(OH)]";
 if(ParallelDescriptor::IOProcessor()) std::cout << " Check 1" << std::endl;
  Array<int> destFill(nVars);
  for (int i=0; i< nVars; ++i)
   {   
       destFill[i] = i;
   }
 
 if(ParallelDescriptor::IOProcessor()) std::cout << " Check 2" << std::endl;

  


   Array<string> varO;
    // Add physical coordinates to bottom end of this list, set index to isoCompName
//   Array <Real > DatOut(nPlotFiles*8);
    int numC=0; 
 
    Array<Real> Speed2(4,-1);
   for (int iPlot=0; iPlot< nPlotFiles ; iPlot++) {
    Real isoVal = isoVal_DEF; pp.query("isoVal",isoVal);
    string isoCompName = isoCompName_DEF; pp.query("isoCompName",isoCompName);

    Array<int> comps;
    if (int nc = pp.countval("comps"))
    {
        comps.resize(nc);
        pp.getarr("comps",comps,0,nc);
    }
    else
    {
        int sComp = amrDataPtrArray[0]->StateNumber(isoCompName);
        pp.query("sComp",sComp);
        int nComp = 1;
        pp.query("nComp",nComp);
//        BL_ASSERT(sComp+nComp <= amrData.NComp());
        comps.resize(nComp);
        for (int i=0; i<nComp; ++i)
            comps[i] = sComp + i;
    }
 
     int isoComp = -1;
    const Array<std::string>& pltNames = amrDataPtrArray[iPlot]->PlotVarNames();
    Array<string> varnames(comps.size());
    varO.resize(varnames.size());
    for (int i=0; i<varnames.size(); ++i) {
        if (comps[i]>=pltNames.size())
            BoxLib::Abort("At least one of the components requested is not in pltfile");
        varnames[i] = pltNames[comps[i]];
        varO[i] = varnames[i];
        if(ParallelDescriptor::IOProcessor())
          std::cout << varnames[i] << std::endl;
        if (varnames[i]==isoCompName) isoComp = i;
    }
    comps.resize(comps.size()+BL_SPACEDIM);
    for (int i=comps.size()-1; i>=BL_SPACEDIM; --i) {
        comps[i] = comps[i-BL_SPACEDIM] + BL_SPACEDIM;
    }
    for (int i=0; i<BL_SPACEDIM; ++i)
        comps[i] = i;
    isoComp += BL_SPACEDIM;
    if (isoComp<BL_SPACEDIM)
        BoxLib::Abort("isoCompName not in list of variables to read in");

    Box subbox;
    if (int nx=pp.countval("box"))
    {
        Array<int> barr;
        pp.getarr("box",barr,0,nx);
        int d=BL_SPACEDIM;
        BL_ASSERT(barr.size()==2*d);
        subbox=Box(IntVect(D_DECL(barr[0],barr[1],barr[2])),
                   IntVect(D_DECL(barr[d],barr[d+1],barr[d+2])));

    }
    else
    {
        subbox = amrDataPtrArray[iPlot]->ProbDomain()[0];
    }



    Array<BoxArray> gridArray(Nlev);
    Array<Box> subboxArray(Nlev);

        int nGrowPer = 0; pp.query("nGrowPer",nGrowPer);
    int levNGP = nGrowPer;
    PArray<Geometry> geoms(Nlev);
    for (int lev=0; lev<Nlev; ++lev)
    {
        subboxArray[lev]
            = (lev==0 ? subbox
               : Box(subboxArray[lev-1]).refine(amrDataPtrArray[iPlot]->RefRatio()[lev-1]));

        if (nGrowPer>0)
        {
            geoms.set(lev,new Geometry(amrDataPtrArray[iPlot]->ProbDomain()[lev]));
            if (lev>0)
                levNGP *= amrDataPtrArray[iPlot]->RefRatio()[lev-1];
            for (int i=0; i<BL_SPACEDIM; ++i)
            {
                if (geoms[lev].isPeriodic(i))
                {
                    if (subboxArray[lev].smallEnd()[i] == amrDataPtrArray[iPlot]->ProbDomain()[lev].smallEnd()[i])
                        subboxArray[lev].growLo(i,levNGP);
                    if (subboxArray[lev].bigEnd()[i] == amrDataPtrArray[iPlot]->ProbDomain()[lev].bigEnd()[i])
                        subboxArray[lev].growHi(i,levNGP);
                }
            }
        }

        gridArray[lev] = BoxLib::intersect(amrDataPtrArray[iPlot]->boxArray(lev), subboxArray[lev]);

        if (nGrowPer>0)
        {

            for (int i=0; i<BL_SPACEDIM; ++i)
            {
                if (geoms[lev].isPeriodic(i))
                {
                    const Box& probDomain = amrDataPtrArray[iPlot]->ProbDomain()[lev];
                    Box gPD = subboxArray[lev];
                    BoxList baL =
                        BoxList(BoxArray(amrDataPtrArray[iPlot]->boxArray(lev)).shift(
                            i,-probDomain.length(i))).intersect(gPD);
                    BoxList baR =
                        BoxList(BoxArray(amrDataPtrArray[iPlot]->boxArray(lev)).shift(
                            i,+probDomain.length(i))).intersect(gPD);
                    BoxList newBL = gridArray[lev].boxList();
                    newBL.join(baL);
                    newBL.join(baR);
                    gridArray.set(lev,BoxArray(newBL));
                }
            }
        }

        if (!gridArray[lev].size())
        {
            Nlev = lev;
            gridArray.resize(Nlev);
            subboxArray.resize(Nlev);
        }
    }

    set<Node> nodeSet;
    set<Element> eltSet;
    int nodeCtr = 0;

    const int nGrow = 1;
    const int nComp=comps.size();
    numC = nComp;
    for (int lev=0; lev<Nlev; ++lev)
    {
        const BoxArray gGridArray =
            BoxLib::intersect(BoxArray(gridArray[lev]).grow(1),amrDataPtrArray[iPlot]->ProbDomain()[lev]);

        MultiFab state(gGridArray,nComp,0);

        // Turns out we can't trust the DxLevel here (ASCII I/O precision issues), so we make our own
        //const Array<Real>& dxf = amrData.DxLevel()[lev];
        Array<Real> dxf(BL_SPACEDIM);
        for (int i=0; i<BL_SPACEDIM; ++i)
            dxf[i] = amrDataPtrArray[iPlot]->ProbSize()[i] / amrDataPtrArray[iPlot]->ProbDomain()[lev].length(i);

        const Array<Real>& plo = amrDataPtrArray[iPlot]->ProbLo();

        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = state[mfi];
            const Box& box = fab.box();

            if (lev!=0)
            {
                const int ratio = amrDataPtrArray[iPlot]->RefRatio()[lev-1];
                FORT_SETCLOC(box.loVect(), box.hiVect(),
                             fab.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()),
                             dxf.dataPtr(), plo.dataPtr(), &ratio);
            }
        }

        {
            MultiFab vstate(gridArray[lev],BL_SPACEDIM,0);
            for (MFIter mfi(vstate); mfi.isValid(); ++mfi)
            {
                const Box& vbox = mfi.validbox();
                FArrayBox& fab = vstate[mfi];
                const Box& box = fab.box();
                FORT_SETLOC(vbox.loVect(), vbox.hiVect(),
                            fab.dataPtr(), ARLIM(box.loVect()), ARLIM(box.hiVect()),
                            dxf.dataPtr(), plo.dataPtr());
            }
            state.copy(vstate,0,0,BL_SPACEDIM);
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            cerr << "Filling interp data at level " << lev << ": ";

        MultiFab tmp(gGridArray,1,0);
        for (int i=0; i<varnames.size(); ++i)
        {
            amrDataPtrArray[iPlot]->FillVar(tmp,lev,varnames[i],0);
            for (MFIter mfi(state); mfi.isValid(); ++mfi)
                state[mfi].copy(tmp[mfi],0,BL_SPACEDIM+i,1);
            amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(varnames[i]));
            if (verbose && ParallelDescriptor::IOProcessor())
                cerr << varnames[i] << " ";
        }
        if (verbose && ParallelDescriptor::IOProcessor()) cerr << '\n';

        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const FArrayBox& sfab = state[mfi];

            // Build a mask using tmp...we know it is built on the box array in state
            FArrayBox& mask = tmp[mfi];
            const Box& myGbox = mask.box();

            mask.setVal(1.0);
            if (lev<finestLevel)
            {
                const int ratio = amrDataPtrArray[iPlot]->RefRatio()[lev];
                const BoxArray& fineBoxes = gridArray[lev+1];
                for (int i=0; i<fineBoxes.size(); ++i) {
                    const Box cgFineBox = Box(fineBoxes[i]).coarsen(ratio);
                    const Box isect = myGbox & cgFineBox;
                    if (isect.ok())
                        mask.setVal(-1.0,isect,0);
                }
            }

            // For looping over quads/bricks, make set of "base" points
            Box loopBox(myGbox);
            for (int i=0; i<BL_SPACEDIM; ++i) {
                loopBox.growHi(i,-1);
            }

            PMap vertCache;

#if BL_SPACEDIM==2
            SegList segments;
            for (IntVect iv=loopBox.smallEnd(); iv<=loopBox.bigEnd(); loopBox.next(iv))
            {
                Array<Segment> eltSegs = Segmentise(sfab,mask,vertCache,iv,isoVal,isoComp);
                for (int i = 0; i < eltSegs.size(); i++)
                {
                    // std::cout << "adding segment: " << eltSegs[i] << '\n';
                    segments.push_back(eltSegs[i]);
                }
            }
#else
            list<Triangle> triangles;
            for (IntVect iv=loopBox.smallEnd(); iv<=loopBox.bigEnd(); loopBox.next(iv))
            {
                Array<Triangle> eltTris = Polygonise(sfab,mask,vertCache,iv,isoVal,isoComp);
                for (int i = 0; i < eltTris.size(); i++)
                triangles.push_back(eltTris[i]);
            }
#endif

            std::map<PMapIt,std::set<Node>::iterator,PMapItCompare> PMI_N_map;

#if 0
            for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it)
            {
                Node n(it->second,-1);
                std::set<Node>::iterator nodeIt = nodeSet.find(n);
                if (nodeIt == nodeSet.end())
                {
                    n.m_idx = nodeCtr++;
                    nodeIt = nodeSet.insert(n).first;
                }
                PMI_N_map[it] = nodeIt;
            }
#else
            for (PMapIt it=vertCache.begin(); it!=vertCache.end(); ++it)
            {
                Node n(it->second,nodeSet.size());
                std::pair<std::set<Node>::iterator,bool> nsit = nodeSet.insert(n);
                if (nsit.second)
                {
                    PMI_N_map[it] = nsit.first;
                }
                else
                {
                    PMI_N_map[it] = nodeSet.find(n);
                }
            }
#endif
            std::vector<int> v(BL_SPACEDIM);

#if BL_SPACEDIM==2
            for (std::list<Segment>::const_iterator it=segments.begin(); it!=segments.end(); ++it)
            {
                for (int k=0; k<BL_SPACEDIM; ++k)
                {
                    v[k] = PMI_N_map[(*it)[k]]->m_idx;
                }
                if (v[0] != v[1]) eltSet.insert(Element(v));
            }
#else
            for (std::list<Triangle>::const_iterator it=triangles.begin(); it!=triangles.end(); ++it)
            {
                for (int k=0; k<BL_SPACEDIM; ++k)
                {
                    v[k] = PMI_N_map[(*it)[k]]->m_idx;
                }
                const bool degenerate = (v[0]==v[1] || v[1]==v[2] || v[0]==v[2]);
                if (!degenerate) eltSet.insert(Element(v));
            }
#endif
        }
    }

    // simulate a sort of the node set to node ordering to be consistent with element pointers
    std::vector<std::set<Node>::iterator> sortedNodes(nodeSet.size());
    for (std::set<Node>::iterator it=nodeSet.begin(); it!=nodeSet.end(); ++it)
        sortedNodes[it->m_idx] = it;

    const int nReal = comps.size()*nodeSet.size();
    Array<Real> nodeRaw(nReal);
    for (long i=0; i<sortedNodes.size(); ++i)
    {
        const Real* vec = sortedNodes[i]->m_vec;
        const int N = sortedNodes[i]->m_size;
        for (int j=0; j<N; ++j) {
            nodeRaw[i*N+j] = vec[j];
        }
    }

    int cnt = 0;
    Array<int> eltRaw(BL_SPACEDIM*eltSet.size());
    for (std::set<Element>::const_iterator it = eltSet.begin(); it != eltSet.end(); ++it)
    {
        const Element& e = *it;
        if (e[0]>=sortedNodes.size() || e[1]>=sortedNodes.size() || e[0]<0 || e[1]<0
#if BL_SPACEDIM==3
              || e[2]<0 || e[2]>=sortedNodes.size()
#endif
            )
        {
           cerr << "***** " << ParallelDescriptor::MyProc() << " *** bad element" << '\n';
           cerr << e[0] << ' ' << e[1]
#if BL_SPACEDIM==3
                << ' ' << e[2]
#endif
                << '\n';
           BoxLib::Abort("Bailing...");
        }
        const int N = e.size();
        for (int j=0; j<N; ++j) {
            eltRaw[cnt*N+j] = e[j];
        }
        ++cnt;
    }

    // All relevant data now in "raw" arrays
    nodeSet.clear();
    eltSet.clear();
    sortedNodes.clear();

    // Communicate node and element info from all procs to IOProc
    Collate(nodeRaw,eltRaw,nComp);
    //int nnode = nodeRaw.size()/nComp;
     
    Real MinY=1e10;
/*    int index =0 ;
    for (int i=0 ; i < nnode ; i++)
     
    {    if(nodeRaw[i*nComp] < 1e-5 & nodeRaw[i*nComp]>-1e-5 & nodeRaw[i*nComp+1]<0) {
         if (nodeRaw[i*nComp+2] < MinY)
          {
              MinY = nodeRaw[i*nComp+2];
              index = i;
         }
      }
     }
 
   for (int i=0 ;i < nComp ; i++ )
       {
           DatOut[iPlot*nComp +i ] = nodeRaw[index*nComp+i];
           if(ParallelDescriptor::IOProcessor())
                 std::cout << nodeRaw[index*nComp+i] << std::endl;
             }
 
   if(ParallelDescriptor::IOProcessor())
         std::cout << "index is " << index << " Min is " << MinY << std::endl;
*/  

    if (ParallelDescriptor::IOProcessor())
    {
        // Uniquify nodes, and make elements consistent

        std::vector<Real> newData(nComp);
        int nNodes = nodeRaw.size()/nComp;

        int nodeCtr = 0;
        std::vector<std::set<Node>::iterator> nodeVec(nNodes);

	// Populate nodeSet and nodeVec
	BL_ASSERT(nodeSet.size()==0);

	if (verbose)
	  {
	    cout << "Number of collated nodes: " << nNodes << endl;
	    cout << "  size of nodeSet iterator: " << sizeof(std::set<Node>::iterator) << endl;
	  }

        for (int j=0; j<nNodes; ++j)
        {
            for (int k=0; k<nComp; ++k)
                newData[k] = nodeRaw[j*nComp+k];

#if 0
            Node n(newData,-1);
            std::set<Node>::iterator nodeIt = nodeSet.find(n);
            if (nodeIt == nodeSet.end())
            {
                n.m_idx = nodeCtr++;
                nodeIt = nodeSet.insert(n).first;
            }
            nodeVec[j] = nodeIt;
#else
            Node n(newData,nodeSet.size());
            std::pair<std::set<Node>::iterator,bool> nsit = nodeSet.insert(n);
            if (nsit.second)
            {
                nodeVec[j] = nsit.first;
            }
            else
            {
                nodeVec[j] = nodeSet.find(n);
            }
#endif
        }
	// Free up some memory
	nodeRaw.clear();

        if (verbose)
            cout << "  final node merge complete, total number nodes now " << nodeSet.size() << endl;

        // simulate a sort of the node set to node ordering to be consistent with element pointers
        sortedNodes.resize(nodeSet.size());
        for (std::set<Node>::iterator it=nodeSet.begin(); it!=nodeSet.end(); ++it)
            sortedNodes[it->m_idx] = it;

        std::vector<int> eltData(BL_SPACEDIM);
        int nElts = eltRaw.size()/BL_SPACEDIM;

	if (verbose)
	  {
	    cout << "Number of collated elements: " << nElts << endl;
	    cout << "  size of Element: " << sizeof(Element) << endl;
	  }

	BL_ASSERT(eltSet.size()==0);
        for (int j=0; j<nElts; ++j)
        {
            for (int k=0; k<BL_SPACEDIM; ++k)
                eltData[k] = (nodeVec[eltRaw[j*BL_SPACEDIM+k]])->m_idx;

            eltSet.insert(eltData);
        }
        nElts = eltSet.size();

	// Clear up some memory
        nodeVec.clear();
	eltRaw.clear();
        
        if (verbose)
            cout << "  final element renumbering complete, total number elements now " << eltSet.size() << endl;

#if 0
	int my_chunk_size=1; pp.query("chunk_size",my_chunk_size);
	bool good=true;
	size_t tot=0;
	if (ParallelDescriptor::IOProcessor())
	  {
	    while(good)
	      {
		tot += my_chunk_size;
		cout << "attempting " << tot << endl;
		char* tester = new char[tot];
		delete tester;
	      }
	  }
#endif

        bool writeSurf = true; pp.query("writeSurf",writeSurf);
        if (writeSurf)
        {
            cout << "...write surface in mef format (mef = Marcs element format)" << endl;

            cout << "      (Nelts,Nnodes):(" << nElts << ", " << sortedNodes.size() << ")" << endl;

            // eltRaw presently contains elts prior to uniquefying, shrink and reload correctly
            const int nodesPerElt = BL_SPACEDIM;
            eltRaw.resize(nElts*nodesPerElt);
            int icnt = 0;
            for (std::set<Element>::const_iterator eit=eltSet.begin(); eit!=eltSet.end(); ++eit)
            {
                for (int i=0; i<nodesPerElt; ++i)
                {
                    eltRaw[icnt++] = (*eit)[i] + 1;
                }
            }
            // Get back some memory
            eltSet.clear();

            int nNodeSize = sortedNodes[0]->m_size;

            // If the surface is large, write data to disk/clear mem/read up into a single fab
            bool surface_is_large = false; pp.query("surface_is_large",surface_is_large);
            int chunk_size = 32768; pp.query("chunk_size",chunk_size);
            FABdata* tmpDataP;
            if (surface_is_large)
            {
                std::string tmpFile="isoTEMPFILE"; pp.query("tmpFile",tmpFile);
                std::ofstream ost;
                ost.open(tmpFile.c_str());

                size_t Npts = sortedNodes.size();
                Box bigBox(IntVect::TheZeroVector(),IntVect(D_DECL(Npts-1,0,0)));
                BoxArray little_boxes(bigBox); little_boxes.maxSize(chunk_size);
                FArrayBox little_fab;
                if (verbose)
                    cout << "  staging vertex data to disk in " << little_boxes.size() << " chunks..." << endl;

                for (int j=0; j<little_boxes.size(); ++j)
                {
                    const Box& little_box = little_boxes[j];
                    little_fab.resize(little_box,nNodeSize);
                    Real* fdat = little_fab.dataPtr();
                    size_t icnt = 0;
                    for (int i=little_box.smallEnd(0); i<=little_box.bigEnd(0); ++i)
                    {
                        const Real* vec = sortedNodes[i]->m_vec;
                        for (int k=0; k<nNodeSize; ++k)
                        {
                            fdat[icnt++] = vec[k];
                        }
                    }
                    little_fab.writeOn(ost);
                }
                little_fab.clear();

                ost.close();
                if (verbose)
                    cout << "  ... data staged." << endl;

                // Clear out some memory now
                sortedNodes.clear(); 
                nodeSet.clear();

                //if (verbose)
                 //   cout << "Total memory presently allocated on I/O proc in fab data: " << BoxLib::total_bytes_allocated_in_fabs << " bytes" << endl;

                // Now allocate final fabdata, and populate with file data
                if (verbose)
                    cout << "  allocating final fab data structure (" << Npts*nNodeSize << " data elements) ..." << endl;
#if 0
                for (int i=1; i<=Npts/chunk_size + 1; ++i)
                {
                    int test_size = std::min(size_t(i)*chunk_size, Npts);
                    FABdata* junk = new FABdata(test_size,nNodeSize);
                    cout << "...successfully built structure for " << test_size << " nodes" << endl;
                    delete junk;
                }
#endif
                tmpDataP = new FABdata(Npts,nNodeSize);
                if (verbose)
                    cout << "  .... final fab structure allocated" << endl;
                FABdata& tmpData = *tmpDataP;
                std::ifstream ist;
                ist.open(tmpFile.c_str());

                if (verbose)
                    cout << "  retrieving vertex data back from disk...";
                Real* fdat = tmpData.fab.dataPtr();
                for (int j=0; j<little_boxes.size(); ++j)
                {
                    little_fab.readFrom(ist);
                    const Box& little_box = little_fab.box();
                    int icntL=0;
                    int icnt=little_fab.box().smallEnd(0) * nNodeSize;
                    const Real* fdatL = little_fab.dataPtr();
                    for (int i=little_box.smallEnd(0); i<=little_box.bigEnd(0); ++i)
                    {
                        for (int k=0; k<nNodeSize; ++k)
                        {
                            fdat[icnt++] = fdatL[icntL++];
                        }
                    }
                    BL_ASSERT(icntL==little_boxes[j].numPts() * nNodeSize);
                }                    
                ist.close();
                if (verbose)
                    cout << "  ... data retrieved" << endl;
            }
            else
            {
                tmpDataP = new FABdata(sortedNodes.size(),nNodeSize);
               if(ParallelDescriptor::IOProcessor()) std::cout << "sortednodes size " << sortedNodes.size() << " elt " << eltRaw.size() << std::endl;
                FABdata& tmpData = *tmpDataP;
                Datapt.resize(sortedNodes.size()*nNodeSize);
                size_t icnt = 0;
                Real* fdat = tmpData.fab.dataPtr();
                for (size_t i=0; i<sortedNodes.size(); ++i)
                {
                    const Real* vec = sortedNodes[i]->m_vec;
                    for (int k=0; k<nNodeSize; ++k)
                    {
                        fdat[icnt++] = vec[k];
                        Datapt[i*nNodeSize+k] = vec[k];
                    }
                }
              ConPt.resize(eltRaw.size()); 
              for (int i=0 ; i< eltRaw.size(); i++) 
                ConPt[i] = eltRaw[i];
                   
        }

        
    

            }
  }


      int sz= Datapt.size();
      int sz2 = ConPt.size(); 
      ParallelDescriptor::Bcast(&sz,1,ParallelDescriptor::IOProcessorNumber());
      ParallelDescriptor::Bcast(&sz2,1,ParallelDescriptor::IOProcessorNumber());
      if (!ParallelDescriptor::IOProcessor()){
             Datapt.resize(sz);
             ConPt.resize(sz2);
          }
      ParallelDescriptor::Bcast(Datapt.dataPtr(),sz,ParallelDescriptor::IOProcessorNumber());
      ParallelDescriptor::Bcast(ConPt.dataPtr(),sz2,ParallelDescriptor::IOProcessorNumber());

     int nNodes = sortedNodes.size();
     ParallelDescriptor::Bcast(&nNodes,1,ParallelDescriptor::IOProcessorNumber());
   
    vector<Point> GridPts(nNodes);
    for (int i=0 ; i < nNodes; i++) { 
        GridPts[i].resize(nComp);
        for(int j=0 ; j< nComp ; j++) GridPts[i][j] = Datapt[i*nComp+j];
       }
    int nElts2 = sz2/3;
  
//   std::cout << " nElts " << nElts << std::endl; 
   vector<vector<int> > GridElts(nElts2);
//   ParallelDescriptor::Bcast(eltRaw.dataPtr(),nElts*3,ParallelDescriptor::IOProcessorNumber());
     for(int i=0 ; i<nElts2; i++ ) {
       GridElts[i].resize(3); 
        for( int j=0 ; j< 3 ; j++ ) {
             GridElts[i][j] = ConPt[i*3+j]-1;
             if(GridElts[i][j]>nNodes) std::cout << " Problemo " << nNodes << " " << GridElts[i][j] << std::endl;
         }
        }
//     std::cout << GridElts[0][0] << " " << GridElts[0][1] << " " << GridElts[0][2]  << " " << ParallelDescriptor::MyProc() <<  std::endl;
   Real Length=0;
   SegList2D segments; 
  for (int i=0; i<nElts2; ++i) {

    Array<Segment2D> eltSegs = Segmentise2D(GridPts,GridElts[i],isoVal2D,isoComp2D);
    if (eltSegs.size() > 0) {


    for (int j=0; j<eltSegs.size(); ++j)
    {
      Length += eltSegs[j].Length();
      segments.push_back(eltSegs[j]);
    }
        }
  }
  cout << "Found " << segments.size() << " segments "  << endl;

  // Number the isosurface vertices
  int cnt2=0;
  for (PMapIt2D it=vertCache2D.begin(); it!=vertCache2D.end(); ++it)
    (*it).second.second = cnt2++;

  // Build a reverse vector into the vertCache
  vector<PMapIt2D> vertVec(vertCache2D.size());
  cnt2=0;
  for (PMapIt2D it=vertCache2D.begin(); it!=vertCache2D.end(); ++it)
    vertVec[cnt2++] = it;

  // Find a segment with the specified vertex as one of its endpoints,
  // then assemble the list of segments to form the contour line.  If
  // we finish, and segments remain, start a new line.
  SegList2D segList = segments;
  list<SegList2D> cLines;

if (segList.size()>0)
  {
    int idx = segList.front().ID_l();
    int newIdx;
    cLines.push_back(SegList2D());

    while (segList.begin() != segList.end())
    {
      SegList2D::iterator segIt = FindMySeg2D(segList,idx,vertVec);
      if (segIt != segList.end())
      {
        int idx_l = (*segIt).ID_l();
        int idx_r = (*segIt).ID_r();
        if ( idx_l == idx )
        {
          newIdx = idx_r;
          cLines.back().push_back(*segIt);
        }
        else
        {
          newIdx = idx_l;
          Segment2D newSeg = Segment2D(*segIt); newSeg.flip();
          cLines.back().push_back(newSeg);
        }

        segList.erase(segIt);
      }
      else
      {
        cLines.push_back(SegList2D());
        newIdx = segList.front().ID_l();
      }

      idx = newIdx;
    }

    // Connect up the line fragments as much as possible
    bool changed;
   do
    {
      changed = false;
      for (std::list<SegList2D>::iterator it = cLines.begin(); it!=cLines.end(); ++it)
      {
        if (!(*it).empty())
        {
          const int idx_l = (*it).front().ID_l();
          const int idx_r = (*it).back().ID_r();
          for (std::list<SegList2D>::iterator it1 = cLines.begin(); it1!=cLines.end(); ++it1)
          {
            if (!(*it1).empty() && (*it).front()!=(*it1).front())
            {
              if (idx_r == (*it1).front().ID_l())
              {
                (*it).splice((*it).end(),*it1);
                changed = true;
              }
              else if (idx_r == (*it1).back().ID_r())
              {
                (*it1).reverse();
                for (SegList2D::iterator it2=(*it1).begin(); it2!=(*it1).end(); ++it2)
                  (*it2).flip();
                (*it).splice((*it).end(),*it1);
                changed = true;
              }
              else if (idx_l == (*it1).front().ID_l())
              {
                (*it1).reverse();
                for (SegList2D::iterator it2=(*it1).begin(); it2!=(*it1).end(); ++it2)
                  (*it2).flip();
                (*it).splice((*it).begin(),*it1);
                changed = true;
              }
            }
          }
        }
      }
    } while(changed);
  }
for (std::list<SegList2D>::iterator it = cLines.begin(); it!=cLines.end();)
  {
    if ((*it).empty())
      cLines.erase(it++);
    else
      it++;
  }
 int szT=0;

  int count =0 ;
   Array<int> usecLine(cLines.size(),0);
   for (std::list<SegList2D>::iterator it = cLines.begin(); it!=cLines.end(); ++it) {
       
        Real MinZ(1e10);
      if(it->size()>3500) {

      for(SegList2D::iterator it2=it->begin(); it2!=it-> end(); ++it2) 
         {
           const Point& p0= (*it2)[0]->second.first;
           const Point& p1 = (*it2)[1]->second.first;
           if(std::max(p0[2],p1[2]) < MinZ)
             MinZ=std::max(p0[2],p1[2]);
          }
   if(ParallelDescriptor::IOProcessor())  std::cout << " MinZ " << MinZ << std::endl;
        if(MinZ < 0.0020) {
       usecLine[count] =1;
      szT = szT+it->size();
       }
     }
      count = count+1;
  } 
  std::cout << "new szT " <<  szT << std::endl;
  Array<Real> Iso2D(szT*nComp);

  ParallelDescriptor::Barrier();
//  cerr << "  number of contours " << cLines.size() << endl;
   int count2 = 0;
   int sztemp=0;
   count=0;
  for (std::list<SegList2D>::iterator it = cLines.begin(); it!=cLines.end(); ++it) {
//      std::cout << "useCLine[count2] " << usecLine[count2] << std::endl;
     if (it->size() > 3500 && usecLine[count2]) {
//        if(usecLine[count2]) {     
      for (SegList2D::iterator it2=it->begin(); it2!=it->end(); ++it2) {
        const Point& p0 = (*it2)[0]->second.first;
       for (int in=0; in<nComp;in++)
         Iso2D[count*nComp+in] = p0[in];
  
         count = count+1;
        }
        sztemp=sztemp+it->size();
        }

         count2=count2+1;
   }
  
//  std::cout << sztemp << " SIZE OF ISO2D " <<  "segments.size() "<< segments.size() << std::endl;
   Array<Real> Plo = amrDataPtrArray[iPlot]->ProbLo();
   Array<Real> dx = amrDataPtrArray[iPlot]->DxLevel()[finestLevel];
   Array<Real > Speed(szT*15);
   int Proc2;
    Array<Box> boxes(szT);

  for (int ipt=0; ipt<szT; ipt++ ) {
        Array<int> pos(3);
         pos[0] = (Iso2D[ipt*nComp]-Plo[0])/dx[0]-0.5;
         pos[1] = (Iso2D[ipt*nComp+1]-Plo[1])/dx[1]-0.5;
         pos[2] = (Iso2D[ipt*nComp+2]-Plo[2])/dx[2]-0.5;

          boxes[ipt].setSmall(IntVect(D_DECL(pos[0],pos[1],pos[2])));
          boxes[ipt].setBig(IntVect(D_DECL(pos[0],pos[1],pos[2])));
          boxes[ipt].grow(5);
        if(boxes[ipt].smallEnd(2) < 0) boxes[ipt].setSmall(IntVect(D_DECL(pos[0],pos[1],0)));

      }   
         const  BoxArray ba(boxes.dataPtr(),szT);
     
        MultiFab  MF(ba,nVars+1,0);  //should make appropriate ghost cells and fill them but for this application not needed
        amrDataPtrArray[iPlot]->FillVar(MF,finestLevel, varNames, destFill);
       DistributionMapping dm = MF.DistributionMap();
       MultiFab Diff(ba,nSpec+1,0);
        const  Array<int> localMap = dm.ProcessorMap();
        for (int ib=0; ib<szT; ib++) {
         FArrayBox tmp(boxes[ib],1);
         tmp.setVal(ib);
         if(ParallelDescriptor::MyProc()==localMap[ib])
          MF[ib].copy(tmp,0,nVars,1);
            }
  

      Array<Real> Speed2(15);
     int Speedcnt=0;
   for(MFIter mfi(MF); mfi.isValid(); ++mfi ) {
     Box vBox = mfi.validbox();
     Box FluxBox = Box(vBox).grow(-2);  //4th order for both first and second derivative
     Box SpeedBox = Box(FluxBox).grow(-2);
     FArrayBox RhoD(vBox,nSpec+1); 
      FArrayBox Reac(vBox,nSpec) ;
      FArrayBox rhoY(vBox,nSpec);
      FArrayBox rhoH(vBox,1);
      FArrayBox cp(vBox,1);
      int index = MF[mfi].max(nVars);
      int do_temp =1;
      int do_VelVisc = 0;
      const FArrayBox& Y = MF[mfi];
      const FArrayBox& T = MF[mfi];
      int idRho = nSpec;
      int idTst = nSpec+1;
      rhoY.copy(Y,0,0,nSpec);
      for (int i=0 ; i< nSpec; i++)
          rhoY.mult(Y,idRho,i,1);
      cd.getHmixGivenTY(rhoH,T,Y,vBox,idTst,0,0);
       rhoH.mult(Y,idRho,0,1);
      Real Patm = 60.0;
      cd.reactionRateRhoY(Reac,rhoY,rhoH,T,Patm,vBox,0,0,idTst,0);  
     Array<Real> loc(3),posLo(3);
       loc[0] = Iso2D[index*nComp];
       loc[1] = Iso2D[index*nComp+1];
       loc[2] = Iso2D[index*nComp+2];
     FORT_MIXAVG_RHODIFF_TEMP(vBox.loVect(),vBox.hiVect(),
                                     RhoD.dataPtr(),ARLIM(RhoD.loVect()),ARLIM(RhoD.hiVect()),
                                     T.dataPtr(idTst),ARLIM(T.loVect()),ARLIM(T.hiVect()),
                                     rhoY.dataPtr(),ARLIM(rhoY.loVect()),ARLIM(rhoY.hiVect()),

                                     &Patm, &do_temp, &do_VelVisc);
       Diff[mfi].copy(RhoD,0,0,nSpec+1);
// compute diffusivity from thermal conductivity        
   cd.getCpmixGivenTY(cp,T,Y,vBox,idTst,0,0);
   RhoD.divide(cp,0,nSpec,1);
//   RhoD.divide(Y,idRho,nSpec,1); need rho alpha inside the fortran
   if(RhoD.contains_nan() || cp.contains_nan() || MF[mfi].contains_nan())
       std :: cout << "NAN" << std::endl;
    FORT_COMPDIFFTERMS(vBox.loVect(),vBox.hiVect(),FluxBox.loVect(), FluxBox.hiVect(), 
                               SpeedBox.loVect(), SpeedBox.hiVect(),
                              ARLIM(FluxBox.loVect()), ARLIM(FluxBox.hiVect()), ARLIM(SpeedBox.loVect()),
                              ARLIM(SpeedBox.hiVect()),

                               Y.dataPtr(),         ARLIM(Y.loVect()),      ARLIM(Y.hiVect()),
                               RhoD.dataPtr(),      ARLIM(RhoD.loVect()),   ARLIM(RhoD.hiVect()),
                               T.dataPtr(nSpec+6), ARLIM(Reac.loVect()), ARLIM(Reac.hiVect()),
                               T.dataPtr(idRho), ARLIM(T.loVect()), ARLIM(T.hiVect()), T.dataPtr(idRho+1),
                               Speed2.dataPtr(),dx.dataPtr(),&Nsp,&nSpec,Y.dataPtr(nSpec+2),loc.dataPtr(),Y.dataPtr(nSpec+5),RhoD.dataPtr(nSpec));
//  std::cout << Speed2[0] << " " << Speed2[1] << " " << " " << loc[0] <<" " <<  ParallelDescriptor::MyProc()<< " " << index<<  std::endl; 
    for (int is=0; is<15;is++) Speed[index*15+is]=Speed2[is]; 
    
  }

   ParallelDescriptor::ReduceRealSum(Speed.dataPtr(),szT*15); 
//  for (int ip=0 ;ip<120;ip++) ParallelDescriptor::Bcast(&Proc2,1,ip);
//  ParallelDescriptor::Bcast(Speed.dataPtr(), nPlotFiles*4*nBins, Proc2);
  
   Array<Real > theta (szT);
   Array <Real > R(szT); 
   for (int i=0 ; i <szT; i++) {
     R[i] = std::pow(Iso2D[i*nComp]*Iso2D[i*nComp]+Iso2D[i*nComp+1]*Iso2D[i*nComp+1],0.5);
     Real ang = std::atan(std::abs(Iso2D[i*nComp+1]/Iso2D[i*nComp]))*180.0/3.14;
     if (Iso2D[i*nComp]<0.0 & Iso2D[i*nComp+1]>0.0) ang = 180-ang;
     if (Iso2D[i*nComp]<0.0 & Iso2D[i*nComp+1]<0.0) ang = 180+ang;
     if (Iso2D[i*nComp]>0.0 & Iso2D[i*nComp+1]<0.0) ang = 360-ang;
     theta[i] = ang;
       }

   if(ParallelDescriptor::IOProcessor()) {
      std::string vars = "X Y";
      #if BL_SPACEDIM==3
            vars += " ZC";
      #endif
        for (int j=0; j<varO.size(); ++j)
                vars += " " + varO[j];
          vars+=" Sd_D Sd_R Sd Sz Se k Se_Sd Se_Sz Sd_Dn Sd_Dt Sz_n Sz_t UT2 Reac gradY theta R";
	  FILE *file;
	  char filename[512];
            std::string outfile =infile[iPlot]+"_base.dat";
              sprintf(filename,"%s", outfile.c_str());
              
	  file = fopen(filename,"w");
          std::cout << "      filename = " << filename << std::endl;
              fprintf(file,"%s", vars.c_str());
              fprintf(file,"\n");
	      for (int i=0; i<szT; i++) {
                  for(int j=0 ; j<numC ; j++ )
  		  fprintf(file,"%e ",Iso2D[i*numC+j]);
                  for (int k=0 ; k<15; k++) 
                  fprintf(file,"%e ",Speed[i*15+k] );
                  fprintf(file,"%e ", theta[i]);
                  fprintf(file,"%e ",R[i]);
                  fprintf(file,"\n");
              }
	  fclose(file);
   } 



      int outsize = numC+15+2; // 4 speed components , 2 r and theta will change when I include more speed components

   FArrayBox toWrite(Box(IntVect(D_DECL(0,0,0)),IntVect(D_DECL(szT-1,0,0))),outsize);
   toWrite.setVal(0.0);
   Real* writeData = toWrite.dataPtr();
   for (int i =0 ; i <szT ; i++ )
       {
          for(int j=0 ; j < numC ; j++) 
             writeData[i*(numC+15+2)+j] = Iso2D[i*numC+j];
           for (int k=0 ; k<15 ; k++) 
             writeData[i*(numC+15+2)+numC+k] = Speed[i*15+k];
             writeData[i*(numC+15+2)+numC+15] = theta[i];
 	     writeData[i*(numC+15+2)+numC+16] = R[i];
       }

   Array<int> elts(szT*2);
     cnt2=0;
    int cnt3=0;
   int start=0;
   count2=0;
   for(std::list<SegList2D>::iterator it = cLines.begin(); it!=cLines.end(); ++it){  
         start =cnt2;
     if(it->size()>3500 && usecLine[count2]) {
//          if(usecLine[count2]) {
     for(SegList2D::iterator it2 = (*it).begin();it2!=(*it).end();++it2) {
          
            
          PMapIt2D s1 = ((*it2).p[0]);
          PMapIt2D s2 = ((*it2).p[1]);
          Point p1 =   (*s1).second.first; 
//          Point p2 =  (*s2).second.first;
//          elts[cnt2] = (*s1).second.second+1;
          elts[cnt2] = cnt3+1;
          cnt2++;
//          elts[cnt2] = (*s2).second.second+1;
           elts[cnt2] = cnt3+2;
          cnt2++;
          cnt3++;
            
       
           
	}
        
        elts[cnt2-1]= elts[start] ;
      
 
        }
        count2=count2+1;
   }


     if(ParallelDescriptor::IOProcessor()) { 
          std::string vars = "X Y ";
#if BL_SPACEDIM==3
            vars += "ZC ";
#endif
            for (int j=0; j<varO.size(); ++j)
                vars += varO[j]+" " ;
           vars+=" Sd_D Sd_R Sd Sz Se k Se_Sd Se_Sz Sd_Dn Sd_Dt Sz_n Sz_t UT2 Reac gradY theta R";
            // Build a label and a filename
            char buf[72];
            sprintf(buf,"%d",time[iPlot]);
            string label(buf);

            sprintf(buf, "%g", isoVal2D);
            std::string outfile_DEF = infile[iPlot]+".mef";
            std::string outfile(outfile_DEF); pp.query("outfile",outfile);
            std::cout << outfile << std::endl;
            cout << "  Writing the file..." << endl;
            std::ofstream ofs;
            std::ostream* os =
                (outfile=="-" ? (std::ostream*)(&std::cout) : (std::ostream*)(&ofs) );
            if (outfile!="-")
                ofs.open(outfile.c_str(),std::ios::out|std::ios::trunc|std::ios::binary);
            (*os) << label << endl;
            (*os) << vars << endl;
            (*os) << szT << " " << "2" << endl;
            toWrite.writeOn(*os);
            (*os).write((char*)elts.dataPtr(),sizeof(int)*elts.size());
            if (outfile!="-")
                ofs.close();
            cout << "            ...done" << endl;
        
  }
 

/*  ofstream os("out.dat");
  os << "VARIABLES =";
  for (int i=0; i<newnames.size(); ++i) {
    os << " " << newnames[i];
  }
  os << '\n';
  for (std::list<SegList>::iterator it = cLines.begin(); it!=cLines.end(); ++it) {
//    if (it->size() < 500) continue;
    os << "ZONE F=FEBLOCK ET=LINESEG N=" << it->size()+1 << " E=" << it->size() << endl;



    for (int n=0; n<nComp+2; ++n) {
      int count = 0;
      for (std::list<Segment>::iterator it2=it->begin(); it2!=it->end(); ++it2) {
        const Point& p0 = (*it2)[0]->second.first;
        os << p0[n] << " ";
        ++count;
        if (count%100 == 0) {
          os << '\n';
        }
      }
      os << (it->back()[1]->second.first)[n] << endl;
    }
    int count = 0;
    for (SegList::iterator it2=it->begin(); it2!=it->end(); ++it2) {
      os << count+1 << " " << count+2 << endl;
//        os << (*((*it2).p[0])).second.second+1 <<" "<<   (*((*it2).p[1])).second.second+1<<endl;;
    }

  }
  os.close();

*/

 for(int il=0; il<varNames.size(); il++) 
   amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(varNames[il]));
 for (int il=0; il<varnames.size(); il++)
  amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(varnames[il]));

 







}

     BoxLib::Finalize();
    return 0;
}
