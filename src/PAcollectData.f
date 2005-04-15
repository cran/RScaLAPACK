* ====================================================================
*           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
*              Oak Ridge National Laboratory, Oak Ridge TN.
*        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
*                            Srikanth Yoginath.
*     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
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
      SUBROUTINE PACollectData (RETRESULT,M,N,PGINFO, WORKSIZ)
*
*     .. Scalar Arguments ..
      INTEGER           M
      INTEGER           N 
      INTEGER           WORKSIZ 
*     ..
*     .. Array Arguments ..

      DOUBLE PRECISION   RETRESULT(*)
      INTEGER            PGINFO(*)
*     ..
*
*     Purpose
*     =======
*
*     The Parallel Agent receives the local pieces sent by different processes 
*     of process grid, into a double precision array.    
*
*  ==============================================================================
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
     $                   LDA, MYCOL, MYROW, NPCOL, NPROW, KK, RANK

*     .. Local Arrays
      DOUBLE PRECISION  WORK ( WORKSIZ )
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
      MB = PGINFO(5)
      NB =MB
      NPROW = PGINFO(7)
      NPCOL = PGINFO(8)    
      KK = 1

      ICURROW = 0
      ICURCOL = 0

      IA = 1 
      JA = 1


*
*     Handle the first block of column separately
*
      JN = MIN( ICEIL( JA, NB )  * NB, JA+N-1 )
      JB = JN-JA+1
      DO 60 H = 0, JB-1
         IN = MIN( ICEIL( IA,  MB ) * MB, IA+M-1 )
         IB = IN-IA+1

                RANK = ICURROW * NPCOL + ICURCOL  
                CALL PARecvVectorFromCR (IB, 1, WORK, MB, RANK)
               DO 20 K = 1, IB
                   RETRESULT(KK) = WORK(K)
                   KK = KK + 1
   20          CONTINUE


         ICURROW = MOD( ICURROW+1, NPROW )

*         CALL BLACS_BARRIER( ICTXT, 'All' )
*
*        Loop over remaining block of rows
*
         DO 50 I = IN+1, IA+M-1, MB
            IB = MIN(  MB, IA+M-I )


                RANK = ICURROW * NPCOL + ICURCOL  
                CALL PARecvVectorFromCR (IB, 1, WORK, MB, RANK)
                  DO 40 K = 1, IB
                      RETRESULT(KK) = WORK ( K )
                      KK = KK + 1
   40             CONTINUE

            ICURROW = MOD( ICURROW+1, NPROW )
*            CALL BLACS_BARRIER( ICTXT, 'All' )
   50    CONTINUE
*
        ICURROW = 0 
   60 CONTINUE
*

      ICURCOL = MOD( ICURCOL+1, NPCOL )
*      CALL BLACS_BARRIER( ICTXT, 'All' )
*
*     Loop over remaining column blocks
*
      DO 130 J = JN+1, JA+N-1, NB
         JB = MIN(  NB , JA+N-J )
         DO 120 H = 0, JB-1
            IN = MIN( ICEIL( IA,  MB ) *  MB , IA+M-1 )
            IB = IN-IA+1

                  RANK = ICURROW * NPCOL + ICURCOL  
                  CALL PARecvVectorFromCR (IB, 1, WORK, MB, RANK)

                  DO 80 K = 1, IB

                     RETRESULT (KK) = WORK (K)
                     KK = KK + 1
   80             CONTINUE


            ICURROW = MOD( ICURROW+1, NPROW )
*            CALL BLACS_BARRIER( ICTXT, 'All' )
*
*           Loop over remaining block of rows
*
            DO 110 I = IN+1, IA+M-1,  MB 
               IB = MIN( MB, IA+M-I )

                  RANK = ICURROW * NPCOL + ICURCOL  
                  CALL PARecvVectorFromCR (IB, 1, WORK, MB, RANK)

                     DO 100 K = 1, IB
                         RETRESULT (KK) = WORK (K)
                         KK =KK + 1 
  100                CONTINUE

               ICURROW = MOD( ICURROW+1, NPROW )
*               CALL BLACS_BARRIER( ICTXT, 'All' )
  110       CONTINUE
*

            ICURROW = 0
  120    CONTINUE
*
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*         CALL BLACS_BARRIER( ICTXT, 'All' )
*
  130 CONTINUE

*
*
 9999 FORMAT( D30.18 )
*
      RETURN
*
*     End of PACollectData 
*
      END
