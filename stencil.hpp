/*****************************************************************************

YASK: Yet Another Stencil Kernel
Copyright (c) 2014-2016, Intel Corporation

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

*****************************************************************************/

// This file defines functions, types, and macros needed for the stencil
// kernel.

#ifndef STENCIL_HPP
#define STENCIL_HPP

// Auto-generated macros from foldBuilder.
// It's important that this be included before real_vec_t.hpp
// to properly set the vector lengths.
#include "stencil_macros.hpp"

// Define a folded vector of reals.
#include "realv.hpp"

#include <map>
#include <set>
#include <vector>

#ifdef WIN32
#define _mm_clevict(p,h) true
#define _Pragma(x)
#endif

#if defined(__GNUC__) && !defined(__ICC)
#define __assume(x) true
#define __declspec(x)
#endif

#if (defined(__GNUC__) && !defined(__ICC)) || defined(WIN32)
#define restrict
#define __assume_aligned(p,n)
#endif

// VTune and stub macros.
#ifdef USE_VTUNE
#include "sampling_MIC.h"
#define SEP_PAUSE  VTPauseSampling()
#define SEP_RESUME VTResumeSampling()
#else
#define SEP_PAUSE
#define SEP_RESUME
#endif

// MPI.
#if !defined(USE_MPI) && defined(MPI_VERSION)
#define USE_MPI
#endif
#ifdef USE_MPI
#include "mpi.h"
#else
#define MPI_PROC_NULL (-1)
#define MPI_Barrier(comm)
#define MPI_Comm int
#endif

// OpenMP and stub functions.
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() (1)
#define omp_get_max_threads() (1)
#define omp_get_thread_num()  (0)
#define omp_set_num_threads(n) (void(0))
#endif

// Enable hardware-thread work crew if requested.
#if (__INTEL_CREW) && (defined(__MIC__) || (defined(__linux) && defined(__x86_64)))
#define USE_CREW 1
extern "C" {
    extern void kmp_crew_create();
    extern void kmp_crew_destroy();
    extern int kmp_crew_get_max_size();
}
#define CREW_FOR_LOOP _Pragma("intel_crew parallel for")
#else
#define USE_CREW 0
#define kmp_crew_create()  ((void)0)
#define kmp_crew_destroy() ((void)0)
#define kmp_crew_get_max_size() (1)
#define CREW_FOR_LOOP
#endif

// macro for debug message.
#ifdef TRACE
#define TRACE_MSG(fmt,...) (printf("YASK trace: " fmt "\n",__VA_ARGS__), fflush(0))
#else
#define TRACE_MSG(fmt,...) ((void)0)
#endif

// Size of time dimension required in allocated memory.
// TODO: calculate this per-grid based on dependency tree and
// traversal order.
// TODO: separate required time-step slots vs. those for
// temp work areas.
#ifndef TIME_DIM_SIZE
#define TIME_DIM_SIZE (1)
#endif

// Cluster sizes in vectors.
// This are defaults to override those generated by foldBuilder
// in stencil_macros.hpp.
#ifndef CLEN_T
#define CLEN_T (1)
#endif
#ifndef CLEN_N
#define CLEN_N (1)
#endif
#ifndef CLEN_X
#define CLEN_X (1)
#endif
#ifndef CLEN_Y
#define CLEN_Y (1)
#endif
#ifndef CLEN_Z
#define CLEN_Z (1)
#endif

// Cluster sizes in points.
// This is the number of scalar results calculated by each
// call to the calc_vector function(s) generated by foldBuilder
// in stencil_code.hpp.
#define CPTS_T (CLEN_T * VLEN_T)
#define CPTS_N (CLEN_N * VLEN_N)
#define CPTS_X (CLEN_X * VLEN_X)
#define CPTS_Y (CLEN_Y * VLEN_Y)
#define CPTS_Z (CLEN_Z * VLEN_Z)
#define CPTS (CLEN_T * CPTS_N * CPTS_X * CPTS_Y * CPTS_Z)

// default sizes.
#ifndef DEF_PROBLEM_SIZE
#define DEF_PROBLEM_SIZE (1024)
#endif
#ifndef DEF_WAVEFRONT_REGION_SIZE
#define DEF_WAVEFRONT_REGION_SIZE (512)
#endif
#ifndef DEF_BLOCK_SIZE
#define DEF_BLOCK_SIZE (8)
#endif

// Memory-accessing code.
#include "mem_macros.hpp"
#include "realv_grids.hpp"

namespace yask {

#ifdef MODEL_CACHE
    extern Cache cache;
#endif

    // Default grid layouts.
    // 3D dims are 1=x, 2=y, 3=z.
    // 4D dims are 1=n/t, 2=x, 3=y, 4=z.
    // Last number in 'Layout' name has unit stride, e.g.,
    // Layout321 & Layout1432 have unit-stride in x.
    // Layout123 & Layout1234 have unit-stride in z.
#ifndef LAYOUT_4D
#define LAYOUT_4D Layout_1234
#endif
#ifndef LAYOUT_3D
#define LAYOUT_3D Layout_123
#endif
    typedef RealVecGrid_XYZ<LAYOUT_3D> Grid_XYZ;
    typedef RealVecGrid_NXYZ<LAYOUT_4D> Grid_NXYZ;
    typedef RealVecGrid_TXYZ<LAYOUT_4D> Grid_TXYZ;
    typedef RealVecGrid_TNXYZ<LAYOUT_4D> Grid_TNXYZ; // T and N reduced to 1st dim.
}

#endif
