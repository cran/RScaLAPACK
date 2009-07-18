* ====================================================================
*           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
*              Oak Ridge National Laboratory, Oak Ridge TN.
*        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
*                            Srikanth Yoginath.
*     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
*     Contact: Guruprasad Kora; (865) 576-6210; koragh@ornl.gov
*                 Computer Science and Mathematics Division
*             Oak Ridge National Laboratory, Oak Ridge TN 37831 
*                   (C) 2004 All Rights Reserved
*
*                              NOTICE
*
* Permission to use, copy, modify, and distribute this software and
* its documentation for any purpose and without fee is hereby granted
* provided that the above copyright notice appear in all copies and
* that both the copyright notice and this permission notice appear in
* supporting documentation.
*
* Neither the Oak Ridge National Laboratory nor the Authors make any
* representations about the suitability of this software for any
* purpose.  This software is provided ``as is'' without express or
* implied warranty.
*
* RScaLAPACK (http://www.aspect-sdm.org/Parallel-R) was funded
* as part of the Scientific Data Management Center
* (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
* Scientific Discovery through Advanced Computing (DOE SciDAC) program
* (http://www.scidac.org ). 
* ======================================================================
*      Based on:
*
*      -- ScaLAPACK example code --
*      University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*      and University of California, Berkeley.
*
*      Written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*
* ==============================================================================

      SUBROUTINE CRDistData(A,DESCA,WORK )

*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*     Purpose:
*     =======
*     CRDISTDATA function is used by the spawned processes to receive the 
*     corresponding data chunk of the input data, the distribution of the
*     input data is performed by the Parallel Agent (PADistData)
*
*  =============================================================================
*
*     .. Parameters ..
      INTEGER            NIN
      PARAMETER          ( NIN = 11 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_,KK
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IB, ICTXT, ICURCOL, ICURROW, II, J, JB,
     $                   JJ, K, LDA, M, MYCOL, MYROW, N, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IWORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, DGERV2D, DGESD2D,
     $                   IGEBS2D, IGEBR2D
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

       M = DESCA( M_ )
       N = DESCA( N_ )
*
*
      KK = 1

      II = 1
      JJ = 1
      ICURROW = DESCA( RSRC_ )
      ICURCOL = DESCA( CSRC_ )
      LDA = DESCA( LLD_ )
*
*     Loop over column blocks
*
      DO 50 J = 1, N, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), N-J+1 )
         DO 40 H = 0, JB-1
*
*           Loop over block of rows
*
            DO 30 I = 1, M, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), M-I+1 )
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
                     CALL CRRecvVectorFromPA (IB,1, A(II+(JJ+H-1)*LDA),
     $               DESCA(NB_))
                  END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
   30       CONTINUE
*
            II = 1
            ICURROW = DESCA( RSRC_ )
   40    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   50 CONTINUE
      CALL BLACS_BARRIER( ICTXT, 'A')
*
*
      RETURN
*
*     End of CRDistData 
*
      END
