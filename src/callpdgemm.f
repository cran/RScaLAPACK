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

      SUBROUTINE CALLPDGEMM(PGINFO, MEM, MEMSIZ )
      
      INTEGER PGINFO(9)
      DOUBLE PRECISION MEM(*)
      INTEGER MEMSIZ

*     Purpose:
*     ========
*     This subroutine solves a linear system by calling the ScaLAPACK    
*     routine PDGEMM. The solution is returned back.    
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

      DOUBLE PRECISION   ALPHA
      PARAMETER          ( ALPHA = 1.0D+0 )

      DOUBLE PRECISION   BETA
      PARAMETER          ( BETA = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IAM, ICTXT, INFO, IPA, IPB, 
     $                   IPW, MYCOL,MYROW,N, NB_BLK, NOUT, NPCOL,
     $                   NPROCS, NPROW, WORKSIZ,FAILFLAG,
     $                   M_MAT1, N_MAT1, M_MAT2, N_MAT2, NOUTMAT
      INTEGER            IPC, NP1, NQ1, NP2, NQ2, NP3, NQ3, 
     $                   M_MAT3, N_MAT3 
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), DESCB( DLEN_ ), DESCX( DLEN_ ),
     $                   OUTDIM(3)
      INTEGER            DESCC( DLEN_ )
*     ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   DESCINIT, IGSUM2D, PDGEMM, CRCollectData,
     $                   CRDistData
*     ..
*     .. External Functions ..
      INTEGER            ICEIL, NUMROC
      EXTERNAL           ICEIL, NUMROC
      DOUBLE PRECISION   MPI_WTIME, sTime, eTime
      EXTERNAL           MPI_WTIME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Executable Statements ..
*
*     Get starting information
*
*      sTime = MPI_WTIME()
*
*==============================================================
*     Assign process grid values

      M_MAT1   = PGINFO(1)
      N_MAT1   = PGINFO(2)
      M_MAT2   = PGINFO(3)
      N_MAT2   = PGINFO(4)
      MB_BLK    = PGINFO(5) 
      NB_BLK    = PGINFO(6)
      NPROW = PGINFO(7)
      NPCOL = PGINFO(8)        

      M_MAT3 = M_MAT1
      N_MAT3 = N_MAT2

      NOUTMAT  = 1
      FAILFLAG = 0

      NOUT = 6
      N = PGINFO(1)

*      PRINT *, 'Initialization successful'
*==============================================================
*     Set up processes and memory
*
      CALL BLACS_PINFO( IAM, NPROCS )

*     Define Process Grid

      CALL BLACS_GET( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*      PRINT *, 'BLACS Initialization successful'

*     Do not perform anything is not part of process grid
      IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL )
     $   GO TO 20
* =============================================================
*     Initialize the array descriptors

*     Leading dimension of A

      NP1   = NUMROC( M_MAT1, MB_BLK, MYROW, 0, NPROW )
      NQ1   = NUMROC( N_MAT1, NB_BLK, MYCOL, 0, NPCOL )

      NP2   = NUMROC( M_MAT2, MB_BLK, MYROW, 0, NPROW )
      NQ2   = NUMROC( N_MAT2, NB_BLK, MYCOL, 0, NPCOL )

      NP3 = NP1
      NQ3 = NQ2

      CALL DESCINIT(DESCA, M_MAT1, N_MAT1, MB_BLK, NB_BLK, 0, 0,
     $        ICTXT, MAX( 1, NP1 ), INFO )

      CALL DESCINIT(DESCB, M_MAT2, N_MAT2, MB_BLK, NB_BLK, 0, 0, 
     $        ICTXT, MAX( 1, NP2 ), INFO )

      CALL DESCINIT(DESCC, M_MAT3, N_MAT3, MB_BLK, NB_BLK, 0, 0, 
     $        ICTXT, MAX( 1, NP3 ), INFO )

*      PRINT *, 'DESC Initialization successful'
* ================================================================

*     Assign Pointers into MEM for scalapack arrays

      IPA = 1
      IPB = IPA + DESCA( LLD_ ) * NQ1
      IPC = IPB + DESCB( LLD_ ) * NQ2
      IPW = IPC + DESCC( LLD_ ) * NQ3

      WORKSIZ = NB_BLK

*      PRINT *, 'Pointer assigning Initialization successful'
* ================================================================

*     Get the optimum working size by doing a workspace query


      IF( IPW+WORKSIZ.GT.MEMSIZ ) THEN
        PRINT *, 'NOT ENOUGH MEMORY..', MEMSIZ, IPW+WORKSIZ
        FAILFLAG=1
      END IF

*      PRINT *, 'Check Fail Flag'

      CALL CRCheckFailFlag (FAILFLAG)

      IF ( IAM.EQ.0)
     $      CALL CRSendIntToPA(FAILFLAG, 1, 1202)

      IF ( FAILFLAG.EQ.1 )
     $     GO TO 20

*      PRINT *, 'Memory verification successful'
* =================================================================

*     Distribute the Input Matrix

      CALL CRDistData( MEM( IPA ), DESCA, MEM( IPW ) )
      CALL CRDistData( MEM( IPB ), DESCB, MEM( IPW ) )

*      PRINT *, 'Distribution  successful'
* =================================================================
*     Call PDGEMM function

      CALL PDGEMM( 'No transpose', 'No transpose', M_MAT1, N_MAT2, 
     $             N_MAT1, ALPHA, MEM( IPA ), 1, 1, DESCA, MEM( IPB ), 
     $             1, 1, DESCB, BETA, MEM( IPC ), 1, 1, DESCC )


*      PRINT *, 'Multiplication successful'
* ================================================================
*     Collect Result
*
      IF (IAM.EQ.0) THEN
          CALL CRSendIntToPA(NOUTMAT, 1, 202)
          OUTDIM(1) = 0
            OUTDIM(2) = M_MAT3
            OUTDIM(3) = N_MAT3

            CALL CRSendIntToPA(OUTDIM, 3 , 300)
       ENDIF

      CALL CRCollectData(M_MAT3, N_MAT3, MEM( IPC ), 1, 1, DESCC,
     $                MEM( IPW ) )

*
*      PRINT *, 'SOLVE result collection  successful'
* ==================================================================
*     Exit the Grid

   20 CONTINUE
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
*      eTime = MPI_WTIME()
*      WRITE(*,9997) MYROW, MYCOL, (eTime - sTime)
* 9997 FORMAT( I2, '/', I2, ' : Time taken by PARALLEL ROUTINE 
*     $ (pdgemm) = ', F12.8, ' Sec')

*     End of CALLPDGEMM

      RETURN

      END
