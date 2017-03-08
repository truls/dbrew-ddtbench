// Copyright (c) 2012 The Trustees of University of Illinois. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include "ddtbench.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>

#include <mpi.h>

#include "config.h"

#ifdef HAVE_DBREW
#include <dbrew.h>
#endif

typedef enum _CmpState {
  StateCompare,
  StateCapture,
} CmpState;

struct Verifier {
  void* captured;
  uintptr_t capturedPtr;
  size_t capturedSize;
  CmpState captureState;
};

void print_pid(int myrank) {
  printf("Rank %d is at pid %d\n", myrank, getpid());
}

//! gives list_dim unique numbers in index list back (from the range of
//! 1:global_dim

void utilities_random_array_shuffle( int* index_list, int list_dim, int global_dim ) {

  int i, temp, irandom;
  int* shuffle_array;

  shuffle_array = malloc( global_dim * sizeof(int) );

  for( i=0 ; i<global_dim ; i++ ) {
//! it is i+1 to be compliant with the fortran code
    shuffle_array[i] = i+1;
  }

  for( i=0 ; i<global_dim ; i++ ) {
    irandom = (rand() % global_dim);

    temp = shuffle_array[i];
    shuffle_array[i] = shuffle_array[irandom];
    shuffle_array[irandom] = temp;
  }

  memcpy( &index_list[0], &shuffle_array[0], list_dim*sizeof(int) );

  free( shuffle_array );

}

//! the following subroutines initialize the elements of an (1D/2D/3D/4D) array
//! with a unique number beginning from base
//! for each needed datatype there is a extra subroutine

void utilities_fill_unique_array_1D_float( float* array, int DIM1, int base ) {

  int i;

  for( i=1 ; i<DIM1 ; i++ ) {
    array[i] = base + i;
  }
}

void utilities_fill_unique_array_2D_float( float* array, int DIM1, int DIM2, int base ) {

  int i, j;

  for( j=0 ; j<DIM2 ; j++ ) {
    for( i=0 ; i<DIM1 ; i++ ) {
      array[i+j*DIM1] = base++;
    }
  }
}

void utilities_fill_unique_array_3D_float( float* array, int DIM1, int DIM2, int DIM3, int base ) {

  int i, j, k;

  for( k=0 ; k<DIM3 ; k++ ) {
    for( j=0 ; j<DIM2 ; j++ ) {
      for( i=0 ; i<DIM1 ; i++ ) {
        array[i+DIM1*(j+DIM2*k)] = base++;
      }
    }
  }
}

void utilities_fill_unique_array_4D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int base ) {

  int i, j, k, l;

  for( l=0 ; l<DIM4 ; l++ ) {
    for( k=0 ; k<DIM3 ; k++ ) {
      for( j=0 ; j<DIM2 ; j++ ) {
        for( i=0 ; i<DIM1 ; i++ ) {
          array[i+DIM1*(j+DIM2*(k+DIM3*l))] = base++;
        }
      }
    }
  }
}

void utilities_fill_unique_array_5D_float( float* array, int DIM1, int DIM2, int DIM3, int DIM4, int DIM5, int base ) {

  int i, j, k, l, m;

  for( m=0 ; m<DIM5 ; m++ ) {
    for( l=0 ; l<DIM4 ; l++ ) {
      for( k=0 ; k<DIM3 ; k++ ) {
        for( j=0 ; j<DIM2 ; j++ ) {
          for( i=0 ; i<DIM1 ; i++ ) {
            array[i+DIM1*(j+DIM2*(k+DIM3*(l+DIM4*m)))] = base++;
          }
        }
      }
    }
  }
}

void utilities_fill_unique_array_1D_double( double* array, int DIM1, int base ) {

  int i;

  for( i=0 ; i<DIM1 ; i++) {
    array[i] = base++;
  }
}

void utilities_fill_unique_array_2D_double( double* array, int DIM1, int DIM2, int base ) {

  int i, j;

  for( j=0 ; j<DIM2 ; j++ ) {
    for( i=0 ; i<DIM1 ; i++ ) {
      array[i+j*DIM1] = base++;
    }
  }
}

void utilities_fill_unique_array_3D_double( double* array, int DIM1, int DIM2, int DIM3, int base ) {

  int i, j, k;

  for( k=0 ; k<DIM3 ; k++ ) {
    for( j=0 ; j<DIM2 ; j++ ) {
      for( i=0 ; i<DIM1 ; i++ ) {
        array[i+DIM1*(j+DIM2*k)] = base++;
      }
    }
  }
}

/*uint64_t rewrite(Rewriter* r, ...)
{
  va_list args;
  va_start(args, r);
  dbrew_re
  }*/

#ifdef HAVE_DBREW

static
void setup_rewriter_function(Rewriter* r) {

  uintptr_t mpirtp = dbrew_util_symname_to_ptr(r, "MPIR_ToPointer");
  //assert(mpirtp_dir == mpirtp);
  //uintptr_t printfptr = (uintptr_t) &printf;
  uintptr_t printfptr = dbrew_util_symname_to_ptr(r, "printf@plt");
  //uintptr_t memcpyptr = dbrew_util_symname_to_ptr(r, "memcpy@plt");
  uintptr_t memcpyptr = (uintptr_t) &memcpy;
  //uintptr_t  mpierr = (uintptr_t) &MPIR_Error;
  //uintptr_t  mpierr_dir = (uintptr_t) &MPIR_Error;
  uintptr_t  mpierr = dbrew_util_symname_to_ptr(r, "MPIR_Error");
  //assert(mpierr == mpierr_dir);
  //uintptr_t mpierrsm_dir = (uintptr_t) &MPIR_Err_setmsg;
  //uintptr_t mpierrsm = (uintptr_t) &MPIR_Err_setmsg;
  uintptr_t mpierrsm = dbrew_util_symname_to_ptr(r, "MPIR_Err_setmsg");

  dbrew_config_function_parcount(r, mpirtp, 1);
  dbrew_config_function_par_setname(r, mpirtp, 0, "MPIR_ToPointer arg1");
  dbrew_config_function_setname(r, mpirtp, "MPIR_ToPointer");
  dbrew_config_function_setflags(r, mpirtp,
                                 FC_BypassEmu | FC_SetRetKnownViral | FC_RetValueHint);

  dbrew_config_function_parcount(r, mpierr, 5);
  dbrew_config_function_setflags(r, mpierr, FC_BypassEmu | FC_KeepCallInstr);

  dbrew_config_function_parcount(r, mpierrsm, 6);
  dbrew_config_function_setflags(r, mpierrsm, FC_BypassEmu | FC_KeepCallInstr);

  dbrew_config_function_parcount(r, printfptr, 6);
  dbrew_config_function_setflags(r, printfptr, FC_BypassEmu | FC_KeepCallInstr);

  dbrew_config_function_parcount(r, memcpyptr, 3);
  dbrew_config_function_setflags(r, memcpyptr, FC_BypassEmu | FC_KeepCallInstr | FC_SetReturnDynamic);

  dbrew_return_orig_on_fail(r, false);
}

Rewriter* get_pack_rewriter(bool verbose)
{

  Rewriter* r = dbrew_new();

  uintptr_t pack = (uintptr_t) &MPI_Pack;


  dbrew_set_function(r, (uint64_t)pack);
  if (verbose) {
    dbrew_verbose(r, 1, 1, 1);
    dbrew_optverbose(r, 1);
  }
  dbrew_colorful_output(r, true);
  dbrew_config_function_setname(r, (uint64_t) pack, "MPI_Pack");
  dbrew_config_function_parmap(r, (uint64_t) pack, 7,
                               SPInt(2) // incount
                               | SPInt(3) // datatype
                               | SRPInt(5) // outcount
                               | SRPInt(6) // position
                               | SPInt(7) // comm

                               );
  setup_rewriter_function(r);
  return r;
}

Rewriter* get_unpack_rewriter(bool verbose)
{
  Rewriter* r = dbrew_new();

  uintptr_t unpack = (uintptr_t) &MPI_Unpack;

  dbrew_set_function(r, (uint64_t)unpack);
  if (verbose) {
    dbrew_verbose(r, 1, 1, 1);
    dbrew_optverbose(r, 1);
  }
  dbrew_colorful_output(r, true);
  dbrew_config_function_setname(r, (uint64_t) unpack, "MPI_Unpack");
  dbrew_config_function_parmap(r, (uint64_t) unpack, 7,
                                 SPInt(2) // insize
                               | SRPInt(3) // position
                               | SRPInt(5) // outsize
                               | SPInt(6) // datatype;
                               | SPInt(7) // comm
                               );


  setup_rewriter_function(r);

  return r;
}

void free_rewriter(Rewriter* r) {
  dbrew_free(r);
}

#endif // HAVE_DBREW

void verifier_reset(Verifier* v) {
  if (v->captured)
    free(v->captured);
  v->captured = 0;
  v->capturedPtr = 0;
  v->capturedSize = 0;
  v->captureState = StateCapture;
}

Verifier* verifier_new() {
  Verifier* v = malloc(sizeof(Verifier));
  v->captured = 0;
  verifier_reset(v);
  return v;
}

void verifier_free(Verifier* v) {
  if (v && v->captured)
    free(v->captured);
  free(v);
}

void verifier_capture(Verifier* v, void* data, size_t size) {
  //printf("Capturing data with size %lu\n", size);
  assert(v->captureState == StateCapture);
  v->captured = calloc(size, 1);
  v->capturedPtr = (uintptr_t) data;
  v->capturedSize = size;
  memcpy(v->captured, data, size);
  v->captureState = StateCompare;
}

/*static
void print_vectors(double* v1, double* v2, size_t size) {
  printf("Value of vector 1 was: ");
  for (int i = 0; i < size; i++) {
    printf("%f ", v1[i]);
  }
  printf("\n\n\n");
  printf("Value of vector 1 was: ");
  for (int i = 0; i < size; i++) {
    printf("%f ", v1[i]);
  }
  printf("\n");
  }*/

void verifier_report(void* buf1, void* buf2, size_t size, int maxerrs) {
  int curerr = 0;
  double* fb1 = (double*) buf1;
  double* fb2 = (double*) buf2;
  //float* fb1 = (float*) buf1;
  //float* fb2 = (float*) buf2;
  bool fb1_allzero = true;
  bool fb2_allzero = true;

  unsigned char* bb1 = (unsigned char*) buf1;
  unsigned char* bb2 = (unsigned char*) buf2;

  for (unsigned i = 0; i < size; i ++) {
    if (bb1[i] != bb2[i]) {
      printf("First mismatch found at pos %d, %d != %d\n", i, bb1[i], bb2[i]);
    }
  }
  printf("Comparison successful\n");

  // Check if either buffer is all zero
  for (size_t i = 0; i < size; i++) {
    if (fb1[i] != 0.0) {
      printf("Buffer1 at location %lu is %f buf 2 is %f\n", i, fb1[i], fb2[i]);
      fb1_allzero = false;
      break;
    }
  }
  for (size_t i = 0; i < size; i++) {
    if (fb2[i] != 0.0) {
      printf("Buffer2 at location %lu is %f buf 1 is %f\n", i, fb2[i], fb1[i]);
      fb2_allzero = false;
      break;
    }
  }
  for (size_t i = 0; i < size; i++) {
    if (fb1[i] != fb2[i]) {
      printf("Value mismatch at pos %lu: %f != %f\n", i, fb1[i], fb2[i]);
      curerr++;
      if (curerr == maxerrs) {
        printf("Reached %d errors, stopping comparison\n", maxerrs);
        if (fb1_allzero)
          printf("Buf 1 is all zeros\n");
        if (fb2_allzero)
          printf("Buf 2 is all zeros\n");
        return;
      }
    } else {
      printf("Value match at pos %lu: %f == %f\n", i, fb1[i], fb2[i]);
    }
  }
}

bool verifier_verify(Verifier* v, void* data, size_t size) {
  //printf("Verifying data with size %lu\n", size);
  if((uintptr_t) data != v->capturedPtr) {
    printf("Pointer capture mismatch: %lu, %lu\n", (uintptr_t) data, v->capturedPtr);
    assert((uintptr_t) data == v->capturedPtr);
  }
  if (size != v->capturedSize) {
    printf("Size capturedSize mismatch: %lu %lu\n", size, v->capturedSize);
    assert(size == v->capturedSize);
  }
  assert(v->captureState == StateCompare);
  if (memcmp(data, v->captured, size) != 0) {
    //print_vectors(data, v->captured, size);
    verifier_report(data, v->captured, size, 100000000);
    return false;
  } else {
    return true;
  }
}

void verifier_test(Verifier* v) {
  int count = 600;
  int* a = malloc(count * sizeof(int));

  // Fill buffers
  for (int i = 0; i < count; i++) {
    a[i] = i;
  }
  verifier_capture(v, a, 600);
  memset(a, 0, count * sizeof(int));
  for (int i = 0; i < count; i++) {
    a[i] = i;
  }
  assert(verifier_verify(v, a, count));
  printf("Verification test successful\n");
  free(a);
  verifier_reset(v);
}
