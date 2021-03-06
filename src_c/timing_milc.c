// Copyright (c) 2012 The Trustees of University of Illinois. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
#include "ddtbench.h"

#define itag 0

static inline int idx5D(int x, int y, int z, int t, int u, int DIM1, int DIM2, int DIM3, int DIM4) {
  return x+DIM1*(y+DIM2*(z+DIM3*(t+DIM4*u)));
}

void timing_milc_su3_zdown_ddt( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ) {

  float* array;

  int myrank;
  int i, j, base, bytes;

  MPI_Datatype dtype_su3_vector_t, dtype_temp_t, dtype_su3_zdown_t;
  int typesize;
  MPI_Aint stride;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = ddtmalloc( 2 * 3 * DIM2 * DIM3 * DIM4 * DIM5 * sizeof(float));

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1;
  utilities_fill_unique_array_5D_float( &array[0], 6, DIM2, DIM3, DIM4, DIM5, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2 * typesize;

    timing_init( testname, method, bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
//! modelling the zdown direction
    MPI_Type_contiguous( 6, MPI_FLOAT, &dtype_su3_vector_t );

    MPI_Type_vector( DIM5, DIM2*DIM3/2, DIM2*DIM3*DIM4/2, dtype_su3_vector_t, &dtype_temp_t );

    MPI_Type_size( dtype_su3_vector_t, &typesize );
    stride = typesize * DIM2 * DIM3 * DIM4 * DIM5/2;
    MPI_Type_create_hvector( 2, 1, stride, dtype_temp_t, &dtype_su3_zdown_t );
    MPI_Type_commit( &dtype_su3_zdown_t );

    MPI_Type_free( &dtype_su3_vector_t );
    MPI_Type_free( &dtype_temp_t );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        MPI_Send( &array[0], 1, dtype_su3_zdown_t, 1, itag, local_communicator );
        MPI_Recv( &array[0], 1, dtype_su3_zdown_t, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(Comm);
      } else {
        MPI_Recv( &array[0], 1, dtype_su3_zdown_t, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        MPI_Send( &array[0], 1, dtype_su3_zdown_t, 0, itag, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_su3_zdown_t );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free( array );
}

void timing_milc_su3_zdown_manual( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ) {

  float* array;
  float* buffer;

  int myrank;
  int i, j, k, l, m, n, o;
  int base, bytes, pos, dsize;

  int typesize;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug;

  array = ddtmalloc( 2 * 3 * DIM2 * DIM3 * DIM4 * DIM5 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1;
  utilities_fill_unique_array_5D_float( &array[0], 6, DIM2, DIM3, DIM4, DIM5, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "manual" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    dsize = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2;
    buffer = ddtmalloc(dsize * sizeof(float) );

//! modelling the zdown direction
    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++) {

      if ( myrank == 0 ) {
        pos = 0;
        for( k=0 ; k<DIM5 ; k++ ) {
          for( l=0 ; l<DIM4 ; l+=DIM4/2 ) {
            for( m=0 ; m<DIM3/2 ; m++ ) {
              for( n=0; n<DIM2 ; n++ ) {
                for( o=0 ; o<6 ; o++ ) {
                  buffer[pos++] = array[idx5D(o,n,m,l,k,6,DIM2,DIM3,DIM4)];
                }
              }
            }
          }
        }
        timing_record(Pack);
        MPI_Send( &buffer[0], dsize, MPI_FLOAT, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], dsize, MPI_FLOAT, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(Comm);
        pos = 0;
        for( k=0 ; k<DIM5 ; k++ ) {
          for( l=0 ; l<DIM4 ; l+=DIM4/2 ) {
            for( m=0 ; m<DIM3/2 ; m++ ) {
              for( n=0; n<DIM2 ; n++ ) {
                for( o=0 ; o<6 ; o++ ) {
                  array[idx5D(o,n,m,l,k,6,DIM2,DIM3,DIM4)] = buffer[pos++];
                }
              }
            }
          }
        }
        timing_record(Unpack);
      } else {
        MPI_Recv( &buffer[0], dsize, MPI_FLOAT, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        for( k=0 ; k<DIM5 ; k++ ) {
          for( l=0 ; l<DIM4 ; l+=DIM4/2 ) {
            for( m=0 ; m<DIM3/2 ; m++ ) {
              for( n=0; n<DIM2 ; n++ ) {
                for( o=0 ; o<6 ; o++ ) {
                  array[idx5D(o,n,m,l,k,6,DIM2,DIM3,DIM4)] = buffer[pos++];
                }
              }
            }
          }
        }
        pos = 0;
        for( k=0 ; k<DIM5 ; k++ ) {
          for( l=0 ; l<DIM4 ; l+=DIM4/2 ) {
            for( m=0 ; m<DIM3/2 ; m++ ) {
              for( n=0; n<DIM2 ; n++ ) {
                for( o=0 ; o<6 ; o++ ) {
                  buffer[pos++] = array[idx5D(o,n,m,l,k,6,DIM2,DIM3,DIM4)];
                }
              }
            }
          }
        }
        MPI_Send( &buffer[0], dsize, MPI_FLOAT, 0, itag, local_communicator );
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

  free(array);
}

void timing_milc_su3_zdown_mpi_pack_ddt( int DIM2, int DIM3, int DIM4, int DIM5, int outer_loop, int inner_loop, int* correct_flag, int* ptypesize, char* testname, MPI_File filehandle_debug __attribute__((unused)), MPI_Comm local_communicator ){

  float* array;
  float* buffer;

  int myrank;
  int i, j;
  int base, bytes, pos, dsize;

  MPI_Datatype dtype_su3_vector_t, dtype_temp_t, dtype_su3_zdown_t;
  int typesize;
  MPI_Aint stride;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = ddtmalloc( 2 * 3 * DIM2 * DIM3 * DIM4 * DIM5 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1;
  utilities_fill_unique_array_5D_float( &array[0], 6, DIM2, DIM3, DIM4, DIM5, base );

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    dsize = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2;
    buffer = ddtmalloc( dsize * sizeof(float) );
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = dsize * typesize;

//! modelling the zdown direction
    MPI_Type_contiguous( 6, MPI_FLOAT, &dtype_su3_vector_t );

    MPI_Type_vector( DIM5, DIM2*DIM3/2, DIM2*DIM3*DIM4/2, dtype_su3_vector_t, &dtype_temp_t );

    MPI_Type_size( dtype_su3_vector_t, &typesize );
    stride = typesize*DIM2*DIM3*DIM4*DIM5/2;
    MPI_Type_create_hvector( 2, 1, stride, dtype_temp_t, &dtype_su3_zdown_t );
    MPI_Type_commit( &dtype_su3_zdown_t );

    MPI_Type_free( &dtype_su3_vector_t );
    MPI_Type_free( &dtype_temp_t );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {
        pos = 0;
        MPI_Pack( &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator );
        timing_record(Pack);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(Comm);
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator );
        timing_record(Unpack);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator );
        pos = 0;
        MPI_Pack( &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( buffer, pos, MPI_PACKED, 0, itag, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_su3_zdown_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);

}

void timing_milc_su3_zdown_mpi_pack_ddt_dbrew( int DIM2, int DIM3, int DIM4,
                                               int DIM5, int outer_loop,
                                               int inner_loop, int* correct_flag,
                                               int* ptypesize, char* testname,
                                               MPI_File filehandle_debug __attribute__((unused)),
                                               MPI_Comm local_communicator ){

  float* array;
  float* buffer;

  int myrank;
  int i, j;
  int base, bytes, pos, dsize;

  MPI_Datatype dtype_su3_vector_t, dtype_temp_t, dtype_su3_zdown_t;
  int typesize;
  MPI_Aint stride;

  char method[50];

//! just some statements to prevent compiler warnings of unused variables
//! those parameter are included for future features
  *correct_flag = 0;
  *ptypesize = 0;
//  typesize = filehandle_debug

  array = ddtmalloc( 2 * 3 * DIM2 * DIM3 * DIM4 * DIM5 * sizeof(float) );

  MPI_Comm_rank( local_communicator, &myrank );

  base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1;
  utilities_fill_unique_array_5D_float( &array[0], 6, DIM2, DIM3, DIM4, DIM5, base );
  //for (i = 0; i < 100; i++) {
  //  printf("%f ", array[i]);
  //}
  //printf("\n");

  if ( myrank == 0 ) {
    snprintf( &method[0], 50, "mpi_pack_ddt_dbrew" );

    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2 * typesize;

    timing_init( testname, &method[0], bytes );
  }

  for( i=0 ; i<outer_loop ; i++ ) {
    dsize = 2 * DIM5 * DIM2*DIM3/2 * 3 * 2;
    buffer = ddtmalloc( dsize * sizeof(float) );
    MPI_Type_size( MPI_FLOAT, &typesize );
    bytes = dsize * typesize;

//! modelling the zdown direction
    MPI_Type_contiguous( 6, MPI_FLOAT, &dtype_su3_vector_t );

    MPI_Type_vector( DIM5, DIM2*DIM3/2, DIM2*DIM3*DIM4/2, dtype_su3_vector_t, &dtype_temp_t );

    MPI_Type_size( dtype_su3_vector_t, &typesize );
    stride = typesize*DIM2*DIM3*DIM4*DIM5/2;
    MPI_Type_create_hvector( 2, 1, stride, dtype_temp_t, &dtype_su3_zdown_t );
    MPI_Type_commit( &dtype_su3_zdown_t );

    /// Rewrite pack/unpack functions using DBrew
    INIT_VERIFIER(v, myrank);
    REWRITE_PACK(pr, rp, pos, myrank, false, &array[0], 1, dtype_su3_zdown_t,
                 &buffer[0], bytes, &pos, local_communicator);
    /*Rewriter* pr = 0;
    MPI_Pack_t rp = 0;
    Verifier* v = 0;
    if (myrank == 0) {
      v = verifier_new();
      pr = get_pack_rewriter(false);
      pos = 0;
      rp = (MPI_Pack_t) dbrew_rewrite(pr, &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator);
      }*/

    /*Rewriter* ur = 0;
    MPI_Unpack_t ru = 0;
    if (myrank == 0) {
      ur = get_unpack_rewriter(false);
      pos = 0;
      ru = (MPI_Unpack_t) dbrew_rewrite(ur, &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator);
      }*/
    REWRITE_UNPACK(ur, ru, pos, myrank, false, &buffer[0], bytes, &pos,
                   &array[0], 1, dtype_su3_zdown_t, local_communicator);

    MPI_Type_free( &dtype_su3_vector_t );
    MPI_Type_free( &dtype_temp_t );

    if ( myrank == 0 ) {
      timing_record(DDTCreate);
    }

    for( j=0 ; j<inner_loop ; j++ ) {

      if ( myrank == 0 ) {

        /*pos = 0;
        MPI_Pack( &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator );
        verifier_capture(v, &buffer[0], bytes);
        memset(&buffer[0], 0, bytes);
        pos = 0;
        rp( &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator );
        assert(verifier_verify(v, &buffer[0], bytes));
        verifier_reset(v);*/
        PACK_MAYBE_ASSERT_VALID(v, rp, pos, &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator);

        timing_record(Pack);
        MPI_Send( &buffer[0], pos, MPI_PACKED, 1, itag, local_communicator );
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 1, itag, local_communicator, MPI_STATUS_IGNORE );
        timing_record(Comm);

        /*pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator );
        verifier_capture(v, &array[0], pos);
        pos = 0;
        ru( &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator );
        assert(verifier_verify(v, &array[0], pos));
        verifier_reset(v);*/
        UNPACK_MAYBE_ASSERT_VALID(v, ru, pos, &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator);
        timing_record(Unpack);
      } else {
        MPI_Recv( &buffer[0], bytes, MPI_PACKED, 0, itag, local_communicator, MPI_STATUS_IGNORE );
        pos = 0;
        MPI_Unpack( &buffer[0], bytes, &pos, &array[0], 1, dtype_su3_zdown_t, local_communicator );
        pos = 0;
        MPI_Pack( &array[0], 1, dtype_su3_zdown_t, &buffer[0], bytes, &pos, local_communicator );
        MPI_Send( buffer, pos, MPI_PACKED, 0, itag, local_communicator );
      }

    } //! inner loop

    MPI_Type_free( &dtype_su3_zdown_t );

    free( buffer );

    if ( myrank == 0 ) {
      timing_record(DDTFree);
    }

    verifier_free(v);
    free_rewriter(pr);
    free_rewriter(ur);

  } //! outer loop

  if ( myrank == 0 ) {
    timing_print( 1 );
  }

  free(array);

}
