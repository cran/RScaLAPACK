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

      SUBROUTINE CALLPDPOTRF (PGINFO,  MEM, MEMSIZ)
      INTEGER  PGINFO(*)
      DOUBLE PRECISION MEM (*)
      INTEGER MEMSIZ

*     Purpose:
*     ========
*     This subroutine computes the Choleski factorization on a real symmetric
*     positive-definite square matrix by calling the ScaLAPACK routine PDPOTRF.
*     The computed result is returned back.
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
*    .. Processor variables ..
      INTEGER        NP, NPROW, NPCOL, MB, NB
*
*    .. Local Scalars ..
      CHARACTER*80   FILENAME
      INTEGER        IAM, ICTXT, INFO, IPA, IPU, IPVT, IPW, IPS,
     $               IPTAU, IPACPY, IPVTCPY,TEMP,LOOP,
     $               NPROCS, MYCOL, MYROW, M, N, NKEEP, NROWK,
     $               NROWA, NROWS, NCOL, WORKSIZ, NOUT, SINDX,
     $               UINDX, MAXITER, ROW, COL, MYMEMSIZ

*
*    .. Local Arrays ..
      INTEGER       DESCA ( DLEN_ ), OUTDIM (3)
*
*    .. Local Scalars ..
      INTEGER       I, J, K, L, NOUTMAT
*    ..
*    .. External Subroutines ..
      EXTERNAL      BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, 
     $              BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $              DESCINIT, PDPOTRF, CRDistData, CRCollectData

*    ..
*    .. External Functions ..
      INTEGER       NUMROC 
      EXTERNAL      NUMROC 
*
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*    ..
     
*     Number of output matrices
      NOUTMAT = 1

      ROW = 1
      COL = 1
      M = PGINFO(1) 
      N = PGINFO(2) 
      MB = PGINFO(5) 
      NB = PGINFO(6) 
      NPROW = PGINFO (7)
      NPCOL = PGINFO (8) 
       
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
*
*    .. Initialize the array descriptors for A, U and VT ..
*
      CALL DESCINIT ( DESCA, M, N, MB, NB, 0, 0, ICTXT, MAX (1,NROWA),
     $              INFO)
*
*    .. Assign Pointers into MEM for SCALAPACK arrays, A is
*    .. allocated starting at position MEM (1)
*
      IPA = 1
      IPW = IPA + DESCA ( LLD_ ) * NROWA 
      
      WORKSIZ =  MAX( NB , NP )
      MYMEMSIZ = IPW + WORKSIZ

      IF (MYMEMSIZ.GT.MEMSIZ) THEN
           PRINT *, 'NOT ENOUGH MEMORY .. EXITING'
           PRINT *,'MYMEMSIZ :', MYMEMSIZ,' MEMSIZ: ', MEMSIZ
           GOTO 10
      ENDIF
 
*
*     DISTRIBITE DATA ...
      CALL CRDistData( MEM( IPA ), DESCA, MEM( IPW ) )

*
*    Call PDPOTRF function
*
      CALL PDPOTRF ('U',M,MEM (IPA),1,1,DESCA,INFO)
*
*    COLLECT RESULT
      IF (IAM.EQ.0) THEN

*         Send number of output matrices to PA
          CALL CRSendIntToPA(NOUTMAT, 1, 202)

*         Send the dimensions of the output matrix to PA
          OUTDIM(1) = 0
          OUTDIM(2) = M 
          OUTDIM(3) = M 
          CALL CRSendIntToPA(OUTDIM, 3 , 300)
      ENDIF
*    
      CALL CRCollectData( M, M, MEM( IPA ), 1,1,DESCA,
     $                MEM( IPW ) )
*
*    .. EXIT THE GRID
*
   10 CONTINUE
*
      CALL BLACS_GRIDEXIT( ICTXT )
*
   20 CONTINUE
*
*      IF (IAM.NE.0)
*      CALL BLACS_EXIT( 0 )
*
*
      RETURN 
*
      END
