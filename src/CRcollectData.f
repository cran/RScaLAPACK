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
*     Based on:
*
*     -- ScaLAPACK example code --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     Written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*
* ==============================================================================

      SUBROUTINE CRCollectData(M, N, A, IA, JA, DESCA, WORK )

*     .. Scalar Arguments ..
      INTEGER            IA, ICWRIT, IRWRIT, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * ), WORK( * )
*     ..
*
*     Purpose
*      =======
*
*     Using CRCollectData all the processes of the Process Grid send the
*     Parallel Agent, the result obtained by performing a requested function
*     on their chunk of data 
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NOUT
      PARAMETER          ( NOUT = 13 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IACOL, IAROW, IB, ICTXT, ICURCOL,
     $                   ICURROW, II, IIA, IN, J, JB, JJ, JJA, JN, K,
     $                   LDA, MYCOL, MYROW, NPCOL, NPROW, KK
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_BARRIER, BLACS_GRIDINFO, INFOG2L,
     $                   DGERV2D, DGESD2D
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
      KK = 1

      ICTXT = DESCA( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

      CALL INFOG2L( IA, JA, DESCA, NPROW, NPCOL, MYROW, MYCOL,
     $              IIA, JJA, IAROW, IACOL )
      ICURROW = IAROW
      ICURCOL = IACOL
      II = IIA
      JJ = JJA
      LDA = DESCA( LLD_ )
*
*     Handle the first block of column separately
*
      JN = MIN( ICEIL( JA, DESCA( NB_ ) ) * DESCA( NB_ ), JA+N-1 )
      JB = JN-JA+1
      DO 60 H = 0, JB-1
         IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
         IB = IN-IA+1
            IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN
               CALL CRSendVectorToPA(IB, 1, A(II+(JJ+H-1)*LDA ), 
     $              DESCA(MB_))
            END IF

         IF( MYROW.EQ.ICURROW )
     $      II = II + IB
         ICURROW = MOD( ICURROW+1, NPROW )
*         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Loop over remaining block of rows
*
         DO 50 I = IN+1, IA+M-1, DESCA( MB_ )
            IB = MIN( DESCA( MB_ ), IA+M-I )
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN

                  CALL CRSendVectorToPA(IB, 1, A(II+(JJ+H-1)*LDA ), 
     $              DESCA(MB_))
               END IF

            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
*            CALL BLACS_BARRIER( ICTXT, 'All' )
   50    CONTINUE
*
        II = IIA
        ICURROW = IAROW
   60 CONTINUE
*
      IF( MYCOL.EQ.ICURCOL )
     $   JJ = JJ + JB
      ICURCOL = MOD( ICURCOL+1, NPCOL )
*      CALL BLACS_BARRIER( ICTXT, 'All' )
*
*     Loop over remaining column blocks
*
      DO 130 J = JN+1, JA+N-1, DESCA( NB_ )
         JB = MIN(  DESCA( NB_ ), JA+N-J )
         DO 120 H = 0, JB-1
            IN = MIN( ICEIL( IA, DESCA( MB_ ) ) * DESCA( MB_ ), IA+M-1 )
            IB = IN-IA+1
               IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN

                  CALL CRSendVectorToPA(IB, 1, A(II+(JJ+H-1)*LDA ), 
     $              DESCA(MB_))
               END IF
            IF( MYROW.EQ.ICURROW )
     $         II = II + IB
            ICURROW = MOD( ICURROW+1, NPROW )
*            CALL BLACS_BARRIER( ICTXT, 'All' )
*
*           Loop over remaining block of rows
*
            DO 110 I = IN+1, IA+M-1, DESCA( MB_ )
               IB = MIN( DESCA( MB_ ), IA+M-I )
                  IF( MYROW.EQ.ICURROW .AND. MYCOL.EQ.ICURCOL ) THEN

                     CALL CRSendVectorToPA(IB, 1, A(II+(JJ+H-1)*LDA ), 
     $                   DESCA(MB_))
                  END IF
               IF( MYROW.EQ.ICURROW )
     $            II = II + IB
               ICURROW = MOD( ICURROW+1, NPROW )
*               CALL BLACS_BARRIER( ICTXT, 'All' )
  110       CONTINUE
*
            II = IIA
            ICURROW = IAROW
  120    CONTINUE
*
         IF( MYCOL.EQ.ICURCOL )
     $      JJ = JJ + JB
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*         CALL BLACS_BARRIER( ICTXT, 'All' )
*
  130 CONTINUE
*
      RETURN
*
*     End of CRCollectData 
*
      END
