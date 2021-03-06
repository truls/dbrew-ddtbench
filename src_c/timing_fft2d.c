// Copyright (c) 2012 The Trustees of University of Illinois. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "ddtbench.h"

void timing_fft2d_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_array;

  int myrank;
  int i, j , base, typesize, bytes;

//! variables for the datatype construction
  MPI_Datatype dtype_vector_t, dtype_resize_t, dtype_scatter_t, dtype_gather_t, dtype_complex_t;
  MPI_Aint lb, extent;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  matrix = ddtmalloc( DIM1 * DIM1/procs * 2 * sizeof(double) );
  recv_array = ddtmalloc( DIM1 * DIM1/procs * 2 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_ddt" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1/procs * DIM1 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    MPI_Type_contiguous( 2, MPI_DOUBLE, &dtype_complex_t );
    MPI_Type_vector( DIM1/procs, 1, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = 2 * sizeof(double);
#if MY_MPI_VERSION >= 2
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
#else
    MPI_Aint disp[3] = {lb, 0, extent};
    MPI_Datatype types[3] = {MPI_LB, dtype_vector_t, MPI_UB};
    int blocklens[3] = {1, 1, 1};
    MPI_Type_struct(3, blocklens, disp, types, &dtype_resize_t);
#endif
    MPI_Type_commit(&dtype_resize_t);

    MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &dtype_scatter_t );
    MPI_Type_commit( &dtype_scatter_t  );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_resize_t );

    MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = DIM1/procs * 2 * sizeof(double);
#if MY_MPI_VERSION >= 2
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
#else
    MPI_Aint disp2[3] = {lb, 0, extent};
    MPI_Datatype types2[3] = {MPI_LB, dtype_vector_t, MPI_UB};
    int blocklens2[3] = {1, 1, 1};
    MPI_Type_struct(3, blocklens2, disp2, types2, &dtype_gather_t);
    //MPI_Type_commit(&dtype_resize_t);
#endif
    MPI_Type_commit( &dtype_gather_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_complex_t );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      MPI_Alltoall( &matrix[0], 1, dtype_gather_t, &recv_array[0], 1, dtype_scatter_t, local_communicator );
      MPI_Alltoall( &recv_array[0], 1, dtype_gather_t, &matrix[0], 1, dtype_scatter_t, local_communicator );
      if ( myrank == 0 ) {
        timing_record(Comm);
      }

    }  //! inner loop

    MPI_Type_free( &dtype_gather_t );
    MPI_Type_free( &dtype_scatter_t );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }

  }  //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(recv_array);
  free(matrix);
}

void timing_fft2d_manual( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_buffer;
  double* buffer;

  int myrank;
  int i, j, k, l, base, typesize, bytes;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug;

  matrix = ddtmalloc( 2 * DIM1 * DIM1/procs * sizeof(double) );
  recv_buffer = ddtmalloc( 2 * DIM1 * DIM1/procs * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "manual" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 2 * DIM1/procs * DIM1 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = ddtmalloc( 2 * DIM1 * DIM1/procs * sizeof(double) );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! pack the data
      for( k=0 ; k < procs ; k++ ) {
        for( l=0 ; l<DIM1/procs ; l++ ) {
          memcpy( &buffer[2 * (k * DIM1/procs * DIM1/procs + l * DIM1/procs)], &matrix[2 * (k * DIM1/procs + l * DIM1)], 2 * DIM1/procs * sizeof(double) );
        }
      }

      if ( myrank == 0 ) {
        timing_record(Pack);
      }

      MPI_Alltoall( &buffer[0], 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, &recv_buffer[0], 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, local_communicator );

//! unpack the data
      for( k=0 ; k<DIM1/procs ; k++ ) {
        for( l=0 ; l<DIM1 ; l++ ) {
          memcpy( &matrix[2*(l+k*DIM1)], &recv_buffer[2*(k + l * DIM1/procs)], 2 * sizeof(double));
        }
      }

//! pack the data
      for( k=0 ; k < procs ; k++ ) {
        for( l=0 ; l<DIM1/procs ; l++ ) {
          memcpy( &buffer[2 * (k * DIM1/procs * DIM1/procs + l * DIM1/procs)], &matrix[2 * (k * DIM1/procs + l * DIM1)], 2 * DIM1/procs * sizeof(double) );
        }
      }

      MPI_Alltoall( &buffer[0], 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, recv_buffer, 2*DIM1/procs*DIM1/procs, MPI_DOUBLE, local_communicator );

      if ( myrank == 0 ) {
        timing_record(Comm);
      }

//! unpack the data
      for( k=0 ; k<DIM1/procs ; k++ ) {
        for( l=0 ; l<DIM1 ; l++ ) {
          memcpy( &matrix[2*(l+k*DIM1)], &recv_buffer[2*(k + l * DIM1/procs)], 2 * sizeof(double));
        }
      }

      if ( myrank == 0 ) {
        timing_record(Unpack);
      }

    } //! inner loop

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(recv_buffer);
  free(matrix);
}

void timing_fft2d_mpi_pack_ddt( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_buffer;
  double* buffer;

  int myrank;
  int i, j, base, typesize, bytes, pos;

//! variables for the datatype construction
  MPI_Datatype dtype_vector_t, dtype_resize_t, dtype_scatter_t, dtype_gather_t, dtype_complex_t;
  MPI_Aint lb, extent;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  matrix = ddtmalloc( DIM1 * DIM1/procs * 2 * sizeof(double) );
  recv_buffer = ddtmalloc( DIM1 * DIM1/procs * 2 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1/procs * DIM1 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = ddtmalloc( 2* DIM1 * DIM1/procs * sizeof(double) );
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 2 * DIM1/procs * DIM1/procs * typesize;

    MPI_Type_contiguous( 2, MPI_DOUBLE, &dtype_complex_t );
    MPI_Type_vector( DIM1/procs, 1, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = 2 * sizeof(double);

    //MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
#if MY_MPI_VERSION >= 2
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
#else
    MPI_Aint disp[3] = {lb, 0, extent};
    MPI_Datatype types[3] = {MPI_LB, dtype_vector_t, MPI_UB};
    int blocklens[3] = {1, 1, 1};
    MPI_Type_struct(3, blocklens, disp, types, &dtype_resize_t);
#endif
    MPI_Type_commit(&dtype_resize_t);

    MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &dtype_scatter_t );
    MPI_Type_commit( &dtype_scatter_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_resize_t );

    MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = DIM1/procs * 2 * sizeof(double);
    //MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
#if MY_MPI_VERSION >= 2
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
#else
    MPI_Aint disp2[3] = {lb, 0, extent};
    MPI_Datatype types2[3] = {MPI_LB, dtype_vector_t, MPI_UB};
    int blocklens2[3] = {1, 1, 1};
    MPI_Type_struct(3, blocklens2, disp2, types2, &dtype_gather_t);
    //MPI_Type_commit(&dtype_resize_t);
#endif
    MPI_Type_commit( &dtype_gather_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_complex_t );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++ ) {
//! pack the data
      pos = 0;
      MPI_Pack( &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );
      if ( myrank == 0 ) {
        timing_record(Pack);
      }

      MPI_Alltoall( &buffer[0], bytes, MPI_PACKED, &recv_buffer[0], bytes, MPI_PACKED, local_communicator );

//! unpack the data
      pos = 0;
      MPI_Unpack( &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );

//! pack the data
      pos = 0;
      MPI_Pack( &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );

      MPI_Alltoall( &buffer[0], bytes, MPI_PACKED, &recv_buffer[0], bytes, MPI_PACKED, local_communicator );

      if ( myrank == 0 ) {
        timing_record(Comm);
      }
//! unpack the data
      pos = 0;
      MPI_Unpack( &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );
      if ( myrank == 0 ) {
        timing_record(Unpack);
      }

    } //! inner loop

    MPI_Type_free( &dtype_gather_t );
    MPI_Type_free( &dtype_scatter_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(recv_buffer);
  free(matrix);
}

void timing_fft2d_mpi_pack_ddt_dbrew( int DIM1, int procs, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ) {

  double* matrix;
  double* recv_buffer;
  double* buffer;

  int myrank;
  int i, j, base, typesize, bytes, pos;

//! variables for the datatype construction
  MPI_Datatype dtype_vector_t, dtype_resize_t, dtype_scatter_t, dtype_gather_t, dtype_complex_t;
  MPI_Aint lb, extent;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  matrix = ddtmalloc( DIM1 * DIM1/procs * 2 * sizeof(double) );
  recv_buffer = ddtmalloc( DIM1 * DIM1/procs * 2 * sizeof(double) );

  MPI_Comm_rank( local_communicator, &myrank );

//! ================= initialize the arrays =================

  base = myrank * DIM1 * DIM1/procs * 2 + 1;
  utilities_fill_unique_array_2D_double( &matrix[0], 2*DIM1, DIM1/procs, base );

  //int foo = 0;
  //print_pid(myrank);
  //if (myrank == 0) {
  //  while (foo == 0) {}
  //}


  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt_dbrew" );

    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = DIM1/procs * DIM1 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {

    buffer = ddtmalloc( 2* DIM1 * DIM1/procs * sizeof(double) );
    MPI_Type_size( MPI_DOUBLE, &typesize );
    bytes = 2 * DIM1/procs * DIM1/procs * typesize;

    MPI_Type_contiguous( 2, MPI_DOUBLE, &dtype_complex_t );
    MPI_Type_vector( DIM1/procs, 1, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = 2 * sizeof(double);

    //MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
#if MY_MPI_VERSION >= 2
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_resize_t );
#else
    MPI_Aint disp[3] = {lb, 0, extent};
    MPI_Datatype types[3] = {MPI_LB, dtype_vector_t, MPI_UB};
    int blocklens[3] = {1, 1, 1};
    MPI_Type_struct(3, blocklens, disp, types, &dtype_resize_t);
#endif
    MPI_Type_commit(&dtype_resize_t);

    MPI_Type_contiguous( DIM1/procs, dtype_resize_t, &dtype_scatter_t );
    MPI_Type_commit( &dtype_scatter_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_resize_t );

    MPI_Type_vector( DIM1/procs, DIM1/procs, DIM1, dtype_complex_t, &dtype_vector_t );
    lb = 0;
    extent = DIM1/procs * 2 * sizeof(double);
    //MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
#if MY_MPI_VERSION >= 2
    MPI_Type_create_resized( dtype_vector_t, lb, extent, &dtype_gather_t );
#else
    MPI_Aint disp2[3] = {lb, 0, extent};
    MPI_Datatype types2[3] = {MPI_LB, dtype_vector_t, MPI_UB};
    int blocklens2[3] = {1, 1, 1};
    MPI_Type_struct(3, blocklens2, disp2, types2, &dtype_gather_t);
    //MPI_Type_commit(&dtype_resize_t);
#endif
    MPI_Type_commit( &dtype_gather_t );

    MPI_Type_free( &dtype_vector_t );
    MPI_Type_free( &dtype_complex_t );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }
    INIT_VERIFIER(v, myrank);
    //printf("foo bar\n");
    //int foo = 0;
    //print_pid(myrank);
    //if (myrank == 0) {
    //  while (foo == 0) {}
    //}
    REWRITE_PACK(pr, rp, pos, myrank, false, &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );
    //exit(1);
    //return;
    REWRITE_UNPACK(ur, ru, pos, myrank, false, &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );

    timing_record(Rewrite);

    for( j=0 ; j<inner_loop ; j++ ) {
//! pack the data
      pos = 0;
      if (myrank == 0) {
        PACK_MAYBE_ASSERT_VALID(v, rp, pos, &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(Pack);
      } else {
        MPI_Pack( &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );
      }

      MPI_Alltoall( &buffer[0], bytes, MPI_PACKED, &recv_buffer[0], bytes, MPI_PACKED, local_communicator );

//! unpack the data
      pos = 0;
      MPI_Unpack( &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );

//! pack the data
      pos = 0;
      MPI_Pack( &matrix[0], 1, dtype_gather_t, &buffer[0], bytes, &pos, local_communicator );

      MPI_Alltoall( &buffer[0], bytes, MPI_PACKED, &recv_buffer[0], bytes, MPI_PACKED, local_communicator );

      if ( myrank == 0 ) {
        timing_record(Comm);
      }
//! unpack the data
      pos = 0;
      if (myrank == 0) {
        UNPACK_MAYBE_ASSERT_VALID(v, ru, pos, &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );
        timing_record(Unpack);
      } else {
        MPI_Unpack( &recv_buffer[0], bytes, &pos, &matrix[0], 1, dtype_scatter_t, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_gather_t );
    MPI_Type_free( &dtype_scatter_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }
    verifier_free(v);
    dbrew_free(ur);
    dbrew_free(pr);

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(recv_buffer);
  free(matrix);
}
