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
*       Based on:
*
*       -- ScaLAPACK example code --
*       University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*       and University of California, Berkeley.
*
*       Written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*
* ==============================================================================

      SUBROUTINE PADistData(GlOBALMATRIX, PGINFO, WORKSIZ, MATNUM)
      DOUBLE PRECISION GLOBALMATRIX (*)
      INTEGER WORKSIZ
      INTEGER PGINFO (8)
      INTEGER MATNUM
*
*     Purpose
*     =======
*
*     PADistData distributes the input data matrix among different processes 
*     of the process grid as required by SCALAPACK functions
*
*  ================================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_,KK
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
*     ..
*     .. Local Scalars ..
      INTEGER            H, I, IB, ICURCOL, ICURROW, II, J, JB,
     $                   JJ, K, M, N, NPCOL, NPROW, S2RANK
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION            WORK( WORKSIZ )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, INFOG2L, DGERV2D, DGESD2D,
     $                   IGEBS2D, IGEBR2D
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX,MIN
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*

      IF (MATNUM.EQ.1) THEN
          M = PGINFO ( 1 )
          N = PGINFO ( 2 )
      ELSE
          M = PGINFO (3)
          N = PGINFO (4)
      ENDIF

       MB = PGINFO (5)
       NB = MB
       NPROW = PGINFO(7)
       NPCOL = PGINFO(8)


      KK = 1

      II = 1
      JJ = 1
      ICURROW = 0
      ICURCOL = 0
*
*     Loop over column blocks
*
*      IPA = 1
*      WORKSIZ = NB 
      
      DO 50 J = 1, N,  NB
         JB = MIN(  NB, N-J+1 )
         DO 40 H = 0, JB-1
*
*           Loop over block of rows
*
            DO 30 I = 1, M,  MB
               IB = MIN( MB , M-I+1 )
                     DO 20 K = 1, IB
                        WORK( K ) = GLOBALMATRIX ( KK )
*                         PRINT *, 'in pa_distdata -- ', GLOBALMATRIX(KK)
                        KK= KK + 1
   20                CONTINUE

                      S2RANK = ICURROW * NPCOL + ICURCOL
*                         PRINT *, 'Rank: ', S2RANK, ' data : ', WORK(1)
                      CALL PASendVectorToCR (IB,1,WORK,MB,S2RANK)
               ICURROW = MOD( ICURROW+1, NPROW )
   30       CONTINUE
*
            II = 1
            ICURROW = 0 
   40    CONTINUE
*
         ICURCOL = MOD( ICURCOL+1, NPCOL )
*
   50 CONTINUE
*
*
      RETURN
*
*     End of PADistData
*
      END
