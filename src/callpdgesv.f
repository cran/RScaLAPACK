* =======================================================================
*           R-ScaLAPACK version 0.2:  ScaLAPACK interface to R
*              Oak Ridge National Laboratory, Oak Ridge TN.
*      Authors: David Bauer, Nagiza. F. Samatova, Srikanth Yoginath
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
* R-ScaLAPACK (http://www.aspect-sdm.org/R-ScaLAPACK) was funded
* as part of the Scientific Data Management Center
* (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
* Scientific Discovery through Advanced Computing (DOE SciDAC) program
* (http://www.scidac.org ). 
* ========================================================================
*     Based on:
*
*     -- ScaLAPACK example code --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*
*     Written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
*
* ==============================================================================

      SUBROUTINE CALLPDGESV(PGINFO, MEM, MEMSIZ )
      
      INTEGER PGINFO(9)
      DOUBLE PRECISION MEM(*)
      INTEGER MEMSIZ

*     Purpose:
*     ========
*     This subroutine solves a linear system by calling the ScaLAPACK    
*     routine PDGESV. The solution is returned back.    
*
* ===============================================================================
*
*     .. Parameters ..
      INTEGER            DBLESZ, INTGSZ
      PARAMETER          ( DBLESZ = 8, INTGSZ = 4 )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IAM, ICTXT, INFO, IPA, IPACPY, IPB, IPPIV, IPX,
     $                   IPW, LIPIV, MYCOL,MYROW,N, NB_BLK, NOUT, NPCOL,
     $                   NPROCS, NPROW, NP, NQ, NQRHS, NRHS, WORKSIZ,
     $                   M_MAT1, N_MAT1, M_MAT2, N_MAT2,NOUTMAT
      DOUBLE PRECISION   ANORM, BNORM, EPS, XNORM, RESID
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCX( DLEN_ ),
     $                   OUTDIM(3)
*     ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   DESCINIT, IGSUM2D, PDGESV,CRCollectData,
     $                   CRDistData
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, NUMROC 
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
*     PGINFO(1) - No. of rows in matrix A
*     PGINFO(2) - No. of cols in matrix A 
*     PGINFO(3) - No. of rows in matrix B 
*     PGINFO(4) - No. of cols in matrix B 
*     PGINFO(5) - Row Block size of matrix A
*     PGINFO(6) - Col Block size of matrix B 
*     PGINFO(7) - No. of process rows in the process grid - Row Block Size
*     PGINFO(8) - No. of process cols in the process grid - Col Block Size
*     PGINFO(9) - Function ID 
*
      CALL BLACS_PINFO( IAM, NPROCS )
        M_MAT1   = PGINFO(1)
        N_MAT1   = PGINFO(2)
        M_MAT2   = PGINFO(3)
        N_MAT2   = PGINFO(4)
        MB_BLK    = PGINFO(5) 
        NB_BLK    = PGINFO(6)
        NPROW = PGINFO(7)
        NPCOL = PGINFO(8)        

      NOUTMAT  = 1

        NOUT = 6
        N = PGINFO(1)
*
*     Define process grid
*
      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Go to bottom of process grid loop if this case doesn't use my
*     process
*
      IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $   GO TO 20
*
      NP    = NUMROC( M_MAT1, MB_BLK, MYROW, 0, NPROW )

      NQ    = NUMROC( N_MAT1, NB_BLK, MYCOL, 0, NPCOL )

      NQRHS = NUMROC( N_MAT2, NB_BLK, MYCOL, 0, NPCOL )
*
*     Initialize the array descriptor for the matrix A and B
*

      CALL DESCINIT(DESCA, M_MAT1, N_MAT1, MB_BLK, NB_BLK, 0, 0,
     $        ICTXT, MAX(1, NP ),INFO )


      CALL DESCINIT(DESCB,M_MAT2,N_MAT2, MB_BLK, NB_BLK, 0, 0, ICTXT,
     $          MAX( 1, NP ),INFO )


      CALL DESCINIT( DESCX,M_MAT2, N_MAT2, MB_BLK, NB_BLK, 0, 0, ICTXT,
     $           MAX( 1, NP ), INFO )
*
*     Assign pointers into MEM for SCALAPACK arrays, A is
*     allocated starting at position MEM( 1 )
*
      IPA = 1
*      IPACPY = IPA + DESCA( LLD_ )*NQ
      IPB = IPA + DESCA( LLD_ )*NQ
      IPX = IPB + DESCB( LLD_ )*NQRHS
      IPPIV = IPX + DESCB( LLD_ )*NQRHS
      LIPIV = ICEIL( INTGSZ*( NP+NB_BLK ), DBLESZ )
      IPW = IPPIV + MAX( NP, LIPIV )
*
      WORKSIZ = NB_BLK
*
*     Check for adequate memory for problem size
*
      INFO = 0
      IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
*         IF( IAM.EQ.0 )
*     $      WRITE( NOUT, FMT = 9998 ) 'test', ( IPW+WORKSIZ )*DBLESZ
*         WRITE( NOUT, FMT = 9998 ) IAM, 'test', ( IPW+WORKSIZ )*DBLESZ,
*     $        2*(IPACPY - IPA), 2*(IPX - IPB), IPW - IPPIV, NB_BLK
    
         INFO = 1
      END IF
*
*     Check all processes for an error
*
      CALL IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )
      IF( INFO.GT.0 ) THEN
         IF( IAM.EQ.0 )
     $      WRITE( NOUT, FMT = 9999 ) 'MEMORY'
         GO TO 10
      END IF
*
*     Read from GLOBALMATRIX  and distribute matrices A and B
*

      CALL CRDistData( MEM( IPA ), DESCA, MEM( IPW ) )


      CALL CRDistData( MEM( IPB ), DESCB, MEM( IPW ) )

*
*    Call PDGESV function
*
      CALL PDGESV( N, N_MAT2, MEM( IPA ), 1, 1, DESCA, MEM( IPPIV ),
     $             MEM( IPB ), 1, 1, DESCB, INFO )
*
*     PERFORM MPI_Send giving the number of output matrices sv outputs
*     PERFORM MPI_Send sending the type of collect to be made along
*     with the output matrix dimensions
*
      IF (IAM.EQ.0) THEN
          CALL CRSendIntToPA(NOUTMAT, 1, 202)
          OUTDIM(1) = 0
          OUTDIM(2) = N
          OUTDIM(3) = N_MAT2

          CALL CRSendIntToPA(OUTDIM, 3 , 300)
      ENDIF

      CALL CRCollectData(N, N_MAT2, MEM( IPB ), 1, 1, DESCB,
     $                MEM( IPW ) )
*
   10 CONTINUE
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
   20 CONTINUE
*
*      IF (IAM.NE.0) THEN
*            CALL BLACS_EXIT( 0 )
*      ENDIF
*
 9999 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
 9998 FORMAT( I4, ' Unable to perform ', A, ': need TOTMEM of at least',
     $        I11, ' (', I3, ', ', I3, ', ', I3, ', ', I3, ')')

*     End of CALLPDGESV

      RETURN

      END