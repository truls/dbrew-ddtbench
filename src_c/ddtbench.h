// Copyright (c) 2012 The Trustees of University of Illinois. All rights reserved.
//  Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include <mpi.h>
#include "config.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#ifdef HAVE_DBREW
#include <dbrew.h>
#endif

#ifndef _DDTBENCH_H_
#define _DDTBENCH_H_

#if MY_MPI_VERSION == 1
#define MPI_Get_address MPI_Address
#define MPI_Type_create_hvector MPI_Type_hvector
#define MPI_Type_create_struct MPI_Type_struct
#endif

void print_pid(int mypid);

// Epochs
typedef enum Epochs Epochs;
enum Epochs {
  DDTCreate = 1,
  Rewrite,
  Pack,
  Comm,
  Unpack,
  DDTFree,
  EpochMax,
};

#ifdef HAVE_DBREW
Rewriter* get_unpack_rewriter(bool verbose);
Rewriter* get_pack_rewriter(bool verbose);
void free_rewriter(Rewriter* r);


typedef int (*MPI_Pack_t) (void* inbuf,
                           int incount,
                           MPI_Datatype datatype,
                           void* outbuf,
                           int outcount,
                           int* postition,
                           MPI_Comm comm);

typedef int (*MPI_Unpack_t) (void* inbuf,
                             int insize,
                             int* position,
                             void* outbuf,
                             int outcount,
                             MPI_Datatype datatype,
                             MPI_Comm comm);

// Macros for setting up and freeing DBrew writer/rewriter instances
#define INIT_VERIFIER(verifier, rank)             \
    Verifier* v = 0;                              \
    if (rank == 0)                                \
      v = verifier_new()
#define GET_PACK_REWRITER(rewriter, rank, verbose)                      \
  Rewriter* rewriter = 0;                                               \
  if (rank == 0)                                                        \
    rewriter = get_pack_rewriter(verbose)
#define GET_UNPACK_REWRITER(rewriter, rank, verbose)                    \
  Rewriter* rewriter = 0;                                               \
  if (rank == 0)                                                        \
    rewriter = get_unpack_rewriter(verbose)
#define REWRITER_REWRITE_PACK(rewriter, rewritten, pos, rank, verbose,  \
                              inbuf, incount, datatype, outbuf,         \
                              outsize, position, comm)                  \
  MPI_Pack_t rewritten = 0;                                             \
  if (rank == 0)  {                                                     \
    pos = 0;                                                            \
    rewritten = (MPI_Pack_t) dbrew_rewrite(rewriter, inbuf, incount,    \
                                           datatype, outbuf, outsize,   \
                                           position, comm);             \
    assert(rewritten);                                                  \
  }
#define REWRITER_REWRITE_UNPACK(rewriter, rewritten, pos, rank,         \
                                verbose, inbuf, insize, position,       \
                                outbuf, outcount, datatype, comm)       \
  MPI_Unpack_t rewritten = 0;                                           \
  if (rank == 0) {                                                      \
    rewriter = get_unpack_rewriter(verbose);                            \
    pos = 0;                                                            \
    rewritten = (MPI_Unpack_t) dbrew_rewrite(rewriter, inbuf, insize,   \
                                             position, outbuf,          \
                                             outcount, datatype, comm); \
    assert(rewritten);                                                  \
  }
#define REWRITE_PACK(rewriter, rewritten, pos, rank, verbose, inbuf,    \
                     incount, datatype, outbuf, outsize,                \
                     position, comm)                                    \
  Rewriter* rewriter = 0;                                               \
  MPI_Pack_t rewritten = 0;                                             \
  if (rank == 0)  {                                                     \
    rewriter = get_pack_rewriter(verbose);                              \
    pos = 0;                                                            \
    rewritten = (MPI_Pack_t) dbrew_rewrite(rewriter, inbuf, incount,    \
                                           datatype, outbuf, outsize,   \
                                           position, comm);             \
    assert(rewritten);                                                  \
  }
#define REWRITE_UNPACK(rewriter, rewritten, pos, rank, verbose, inbuf,  \
                       insize, position, outbuf, outcount,              \
                       datatype, comm)                                  \
  Rewriter* rewriter = 0;                                               \
  MPI_Unpack_t rewritten = 0;                                           \
  if (rank == 0) {                                                      \
    rewriter = get_unpack_rewriter(verbose);                            \
    pos = 0;                                                            \
    rewritten = (MPI_Unpack_t) dbrew_rewrite(rewriter, inbuf, insize,   \
                                             position, outbuf,          \
                                             outcount, datatype, comm); \
    assert(rewritten);                                                  \
    }
#endif

// Verifier definitions
typedef struct Verifier Verifier;
void verifier_reset(Verifier* v);
Verifier* verifier_new();
void verifier_free(Verifier* v);
void verifier_capture(Verifier* v, void* data, size_t size);
bool verifier_verify(Verifier* v, void* data, size_t size);
void verifier_report(void* buf1, void* buf2, size_t size, int maxerrs);
void verifier_test(Verifier* verify);

#ifdef VERIFY_
#define PACK_MAYBE_ASSERT_VALID(verifier, rewritten, pos, inbuf,         \
                                incount, datatype, outbuf, outsize,     \
                                position, comm)                         \
  do {                                                                  \
    pos = 0;                                                            \
    /*printf("outbuf before: %p\n", outsize);*/                         \
    memset(outbuf, 0, outsize);                                         \
    MPI_Pack(inbuf, incount, datatype, outbuf,                          \
             outsize, position, comm);                                  \
    /*verifier_report(outbuf, outbuf, outsize, 100000);*/               \
    /*printf("outbuf after1: %p\n", outsize);*/                         \
    verifier_capture(verifier, outbuf, outsize);                        \
    memset(outbuf, 0, outsize);                                         \
    pos = 0;                                                            \
    rewritten(inbuf, incount, datatype, outbuf,                         \
              outsize, position, comm);                                 \
    /*verifier_report(outbuf, outbuf, outsize, 100000);*/               \
    /*printf("outbuf after2: %p\n", pos);*/                             \
    assert(verifier_verify(verifier, outbuf, outsize));                 \
    verifier_reset(verifier);                                           \
  } while(0)
#define UNPACK_MAYBE_ASSERT_VALID(verifier, rewritten, pos, inbuf,      \
                                  insize, position, outbuf,             \
                                  outcount, datatype, comm)             \
  do {                                                                  \
    pos = 0;                                                            \
    MPI_Unpack(inbuf, insize, position, outbuf,                         \
               outcount, datatype, comm);                               \
    verifier_capture(verifier, outbuf, pos);                            \
    /*memset(outbuf, 0, pos);*/                                         \
    pos = 0;                                                            \
    rewritten(inbuf, insize, position, outbuf,                          \
              outcount, datatype, comm);                                \
    assert(verifier_verify(verifier, outbuf, pos));                     \
    verifier_reset(verifier);                                           \
  } while(0)
#else
#define PACK_MAYBE_ASSERT_VALID(verifier, rewritten, pos, inbuf,        \
                                incount, datatype, outbuf, outsize,     \
                                position, comm)                         \
  do {                                                                  \
    pos = 0;                                                            \
    rewritten(inbuf, incount, datatype, outbuf,                         \
              outsize, position, comm);                                 \
  } while(0)
#define UNPACK_MAYBE_ASSERT_VALID(verifier, rewritten, pos, inbuf,      \
                                  insize, position, outbuf, outcount,   \
                                  datatype, comm)                       \
  do {                                                                  \
    pos = 0;                                                            \
    rewritten(inbuf, insize, position, outbuf,                          \
              outcount, datatype, comm);                                \
  } while(0)
#endif

#ifdef USE_ALIGNED_MALLOC
#define ddtmalloc(size) aligned_alloc(32, ((size_t) size + 0x1f) & (((size_t) (-1)) - 0x1f))
#else
#define ddtmalloc(size) malloc(size);
#endif


void timing_close_file();
void timing_init( char* ptestname, char* pmethod, int pbytes );
void timing_open_file( char* filename );
void timing_print( int last );
void timing_record( Epochs id );
void timing_set_max_tests(int value);
#if TEST_TYPE > 1
  void init_papi();
  void cleanup_papi();
#endif
void timing_hrt_init();

void timing_basic_ping_pong_nelements( int DIM1, int loop, char* testname, MPI_Comm local_communicator);
void timing_basic_alltoall_nelements( int DIM1, int procs, int loop, char* testname, MPI_Comm local_communicator);

void timing_fft2d_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_fft2d_manual( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_fft2d_mpi_pack_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_lammps_atomic_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_atomic_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_atomic_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_lammps_atomic_mpi_pack_ddt_dbrew( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_lammps_full_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_full_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_full_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);
void timing_lammps_full_mpi_pack_ddt_dbrew( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator);

void timing_milc_su3_zdown_ddt( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_milc_su3_zdown_manual( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_milc_su3_zdown_mpi_pack_ddt( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_milc_su3_zdown_mpi_pack_ddt_dbrew( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_nas_lu_x_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_x_manual( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_x_mpi_pack_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_x_mpi_pack_ddt_dbrew( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_manual( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_mpi_pack_ddt( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_lu_y_mpi_pack_ddt_dbrew( int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_nas_mg_x_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_x_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_x_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_x_mpi_pack_ddt_dbrew( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_y_mpi_pack_ddt_dbrew( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_nas_mg_z_mpi_pack_ddt_dbrew( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_specfem3D_cm_ddt( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_cm_manual( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_cm_mpi_pack_ddt( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int* list_cm, int* list_ic, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_manual( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int *correct_flag, int *ptypesize, char *testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_mpi_pack_ddt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3d_mt_mpi_pack_ddt_dbrew( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_manual( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_mpi_pack_ddt( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_specfem3D_oc_mpi_pack_ddt_dbrew( int DIM1, int icount, int* list, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void timing_wrf_manual ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_sa_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_sa_mpi_pack_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_vec_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );
void timing_wrf_vec_mpi_pack_ddt ( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js,
  int je, int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug, MPI_Comm local_communicator );

void utilities_fill_unique_array_1D_float( float* array, int DIM1, int base );
void utilities_fill_unique_array_2D_float( float* array, int DIM1, int DIM2, int base );
void utilities_fill_unique_array_3D_float( float* array, int DIM1, int DIM2, int DIM3, int base );
void utilities_fill_unique_array_4D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int base );
void utilities_fill_unique_array_5D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5, int base );
void utilities_fill_unique_array_1D_double( double* array, int DIM1, int base );
void utilities_fill_unique_array_2D_double( double* array, int DIM1, int DIM2, int base );
void utilities_fill_unique_array_3D_double( double* array, int DIM1, int DIM2, int DIM3, int base );
void utilities_random_array_shuffle( int* index_list, int list_dim, int global_dim );

void wrapper_timing_fft( int DIM1, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_lammps_atomic( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator);
void wrapper_timing_lammps_full( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_milc_su3_zdown( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_nas_lu( int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_nas_mg( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_specfem3D_cm( int DIM2_cm, int DIM2_ic, int icount_cm, int icount_ic, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_specfem3d_mt( int DIM1, int DIM2, int DIM3, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );
void wrapper_timing_specfem3D_oc( int DIM1, int icount, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator);
void wrapper_timing_wrf( int number_2D, int number_3D, int number_4D, int ims, int ime, int jms, int jme, int kms, int kme, int* limit_4D_arrays, int is, int ie, int js, int je,
  int ks, int ke, int param_first_scalar, int outer_loop, int inner_loop, MPI_File filehandle_correctness, MPI_File filehandle_debug, char* ptestname, MPI_Comm local_communicator );

#endif // _DDTBENCH_H_
