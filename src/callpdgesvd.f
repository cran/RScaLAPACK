* =======================================================================
*           R-ScaLAPACK version 0.3.x:  R interface to ScaLAPACK
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
      SUBROUTINE CALLPDGESVD (PGINFO, SINGVALS, MEM, MEMSIZ)
*
      INTEGER PGINFO(9)
      DOUBLE PRECISION SINGVALS ( * )
      DOUBLE PRECISION MEM ( * )
      INTEGER MEMSIZ

*     Purpose:
*     ========
*     This subroutine computes the SVD decomposition by calling the ScaLAPACK
*     routine PDGESVD. The computed result is returned back.
*
* ===============================================================================
*
*     .. Array Descriptor Parameters ..
     
      INTEGER   BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DT_,
     $            LLD_, MB_, NB_, N_, RSRC_

      PARAMETER     (BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $              CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $              RSRC_ = 7, CSRC_ = 8, LLD_ =  9 )
*
*     .. Miscellaneous parameters ..
      DOUBLE PRECISION   ONE, ZERO, MAXERROR
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.D0 )
*
*     .. I/O Parameters ..
      INTEGER        NIN, NSTDOUT
      PARAMETER     ( NIN =11, NSTDOUT = 6 )
*
*     .. Memory Parameters ..
      INTEGER       DBLESZ
      PARAMETER     ( DBLESZ = 8 )
*
*    .. Processor variables ..
      INTEGER        NP, NPROW, NPCOL, MB, NB
*
*    .. Local Scalars ..
      CHARACTER*80   FILENAME
      CHARACTER*1    U, V
      INTEGER        IAM, ICTXT, INFO, IPA, IPU, IPVT, IPW, IPS,
     $               IPTAU, TEMP,LOOP,
     $               NPROCS, MYCOL, MYROW, M, N, NKEEP, 
     $               NROWA, NROWS, NCOL, WORKSIZ, NOUT, SINDX,
     $               UINDX, MAXITER, RET_SUBS

      DOUBLE PRECISION   STARTTIME, ENDTIME
*
*    .. Local Arrays ..
      INTEGER       DESCA ( DLEN_ ), DESCU ( DLEN_ ), DESCVT ( DLEN_ ),
     $              DESCVTK ( DLEN_ ), OUTDIM(3) 
*
*    .. Local Scalars ..
      INTEGER       I, J, K, L, NOUTMAT
*    ..
*    .. External Subroutines ..
      EXTERNAL      BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, 
     $              BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $              DESCINIT, PDGESVD, CRDistData, CRCollectData

*    ..
*    .. External Functions ..
      INTEGER       NUMROC 
      EXTERNAL      NUMROC
      DOUBLE PRECISION   MPI_WTIME, dTime0, dTime1, dTime2, dTime3
      EXTERNAL           MPI_WTIME
*
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*    ..
      dTime0 = MPI_WTIME()

      RET_SUBS = 1
      NOUTMAT = 3
     
      M =  PGINFO (1) 
      N =  PGINFO (2)
      MB = PGINFO (5)
      NB = PGINFO (6)
      NPROW = PGINFO (7)
      NPCOL = PGINFO (8)

*	  Parameters 3 and 4 of the input (normally the dimensions of the second
*	      matrix) are used to store the NU and NV parameters.
      IF (PGINFO(3) == 0) THEN
          U = 'N'
      ELSE
          U = 'V'
      ENDIF
      IF (PGINFO(4) == 0) THEN
          V = 'N'
      ELSE
          V = 'V'
      ENDIF

*    ..
*    .. Executable Statements ..
*
*    .. Set up processes and memory
      CALL BLACS_PINFO ( IAM, NP )
*
*    .. Define Process Grid ..
*
      CALL BLACS_GET ( -1, 0, ICTXT )
      CALL BLACS_GRIDINIT ( ICTXT, 'Row-major', NPROW, NPCOL )
      CALL BLACS_GRIDINFO ( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     .. Do not perform anything if not part of the process grid ..
*
      IF ( MYROW.GE.NPROW .OR. MYCOl.GE.NPCOL)
     $    GO TO 20
*
      NROWA = NUMROC ( M, MB, MYROW, 0, NPROW )
      NROWS = NUMROC ( N, MB, MYROW, 0, NPROW )
      NCOL = NUMROC ( N, NB, MYCOL, 0, NPCOL )
*
*    .. Initialize the array descriptors for A, U and VT ..
*
      CALL DESCINIT ( DESCA, M, N, MB, NB, 0, 0, ICTXT, MAX (1,NROWA),
     $              INFO)
      CALL DESCINIT ( DESCU, M, N, MB, NB, 0, 0, ICTXT, MAX (1,NROWA),
     $              INFO)
      CALL DESCINIT ( DESCVT, N, N, MB, NB, 0, 0, ICTXT, MAX(1,NROWS),
     $              INFO)
*
*    .. Assign Pointers into MEM for SCALAPACK arrays, A is
*    .. allocated starting at position MEM (1)
*

      IPA = 1
      IPU = IPA + DESCA ( LLD_ ) * NCOL
      IF (PGINFO(3) == 0) THEN
      IPS = IPU
      ELSE
      IPS = IPU + DESCU ( LLD_ ) * NCOL
      ENDIF
      IPVT = IPS + N
      IF (PGINFO(4) == 0) THEN
          IPTAU = IPVT
      ELSE
          IPTAU = IPVT + DESCVT ( LLD_ ) * NCOL
      ENDIF
      IPW = IPTAU 
*
*    .. Get the optimum working size by doing a workspace query ..
*
      CALL PDGESVD (U, V, M, N, MEM (IPA), 1,1,DESCA, MEM(IPS),
     $              MEM(IPU), 1,1, DESCU, MEM(IPVT), 1,1, DESCVT,
     $              MEM(IPW), -1, INFO )

  
      IF ( INFO.EQ.0) THEN
          WORKSIZ = MEM (IPW)
          
      ELSE
        IF ( IAM.EQ.0)
     $          PRINT *, 'PDGESVD DRY RUN FAILED .. EXITING  = ', INFO
               GO TO 10
      ENDIF
*
      IF (WORKSIZ.LT.(NPROW*MB*N)) WORKSIZ = NPROW*MB*N
*
*    .. Check for adequate memory for problem size
*
      INFO = 0
     
      IF (IPW+WORKSIZ.GT.MEMSIZ) THEN
          IF(IAM.EQ.0)
     $          PRINT *, 'test', (IPW+WORKSIZ)*DBLESZ
          INFO = 1
      ELSE
          WORKSIZ = MEMSIZ - IPW + 1
      ENDIF
*
*    .. check all processes for an error
*
      CALL IGSUM2D ( ICTXT, 'All',' ',1,1,INFO,1,-1,0)
      IF ( INFO.GT.0) THEN
          IF ( IAM.EQ.0)
     $         PRINT *, 'MEMORY PROBLEM'
          GO TO 10 
      ENDIF
*
*    .. Read From file and distribute the matrix
*         

*     Distribute Data ...
      CALL CRDistData( MEM( IPA ), DESCA, MEM( IPW ) )

*
*    Call PDGESVD function
*
      dTime1 = MPI_WTIME()
      CALL PDGESVD (U, V, M, N, MEM (IPA), 1,1,DESCA, MEM(IPS),
     $              MEM(IPU), 1,1, DESCU, MEM(IPVT), 1,1, DESCVT,
     $              MEM(IPW), WORKSIZ, INFO )
      dTime2 = MPI_WTIME()

*
*    Collect Result
*     
      IF ( IAM.EQ.0) THEN

*         Send Parallel Agent 'number of output matrices'

          CALL CRSendIntToPA(NOUTMAT, 1, 202)

          LOOP = MIN (DESCU(N_), DESCU(M_))
          DO TEMP = 1, LOOP 
               SINGVALS(TEMP) =  MEM(IPS + TEMP-1)
          END DO
          
          OUTDIM(1) = 1 
          OUTDIM(2) = 1 
          OUTDIM(3) = LOOP 
          CALL CRSendIntToPA(OUTDIM, 3 , 300)
          CALL CRSendDoubleToPA(SINGVALS, LOOP, 400)
      ENDIF
*
      CALL BLACS_BARRIER (ICTXT,'ALL')

      IF ( M.LE.N ) THEN
*		  If there is a matrix U, then send it.
          IF ( PGINFO(3).NE.0 ) THEN
               OUTDIM(1) = 0 
               OUTDIM(2) = M  
               OUTDIM(3) = M
               IF(IAM.EQ.0) 
     $              CALL CRSendIntToPA(OUTDIM, 3 , 301 )

              CALL CRCollectData( M,M,MEM( IPU ), 1,1,
     $                DESCU, MEM( IPW ) )
                     
*		  Else, send the dimenions [0, 0] to signal a NULL matrix.
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
              CALL CRSendIntToPA(OUTDIM, 3, 301 )
          ENDIF

*		  If there is a matrix Vt, then send it.
          IF ( PGINFO(4).NE.0 ) THEN
              IF (M.EQ.N) THEN
                    OUTDIM(1) = 0
                    OUTDIM(2) = N 
                    OUTDIM(3) = N 
                    IF(IAM.EQ.0) 
     $                   CALL CRSendIntToPA(OUTDIM, 3 , 302 ) 
     
                  CALL CRCollectData( N,N,MEM( IPVT ),
     $                              1, 1, DESCVT,MEM( IPW ) )
              ELSE
                    OUTDIM(1) = 0
                    OUTDIM(2) = M 
                    OUTDIM(3) = N

                    IF(IAM.EQ.0) 
     $                 CALL CRSendIntToPA(OUTDIM, 3 , 302 )

                   CALL CRCollectData(M,N, MEM( IPVT ),
     $                               1, 1, DESCVT, MEM( IPW ) )
              ENDIF
*		  Else, send the dimensions [0, 0] to signal a NULL matrix.
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
              CALL CRSendIntToPA(OUTDIM, 3, 302 )
          ENDIF
     
      ELSE
          IF ( PGINFO(3).NE.0 ) THEN
               OUTDIM(1) = 0
               OUTDIM(2) = M
               OUTDIM(3) = N 
  
               IF(IAM.EQ.0)
     $              CALL CRSendIntToPA(OUTDIM, 3 , 301 )

               CALL CRCollectData( M,N, MEM( IPU ), 1,1,
     $                          DESCU, MEM( IPW ) )
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
              CALL CRSendIntToPA(OUTDIM, 3, 301 )
          ENDIF
          IF ( PGINFO(4).NE.0 ) THEN
               OUTDIM(1) = 0
               OUTDIM(2) = N 
               OUTDIM(3) = N

               IF (IAM.EQ.0)               
     $               CALL CRSendIntToPA(OUTDIM, 3 , 302 )

               CALL CRCollectData( N,N,MEM( IPVT ),1,1,
     $                           DESCVT, MEM( IPW ) ) 
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
              CALL CRSendIntToPA(OUTDIM, 3, 302 )
          ENDIF
      ENDIF

*
*    .. EXIT THE GRID
*
   10 CONTINUE

      CALL BLACS_GRIDEXIT( ICTXT )
*
   20 CONTINUE
*
*      IF (IAM.NE.0)
*         CALL BLACS_EXIT( 0 )
*
*
      dTime3 = MPI_WTIME()
*      print *, "Time for Fortran call =", dTime3 - dTime0, " (",
*     $        dTime1 - dTime0, ",", dTime2 - dTime1, ",",
*     $        dTime3 - dTime2, ")"
      RETURN
*
      END
 
