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

      SUBROUTINE CALLPDSYEVD ( PGINFO, RETVAL)
      INTEGER  PGINFO(*)
      DOUBLE PRECISION RETVAL (*)

*     Purpose:
*     ========
*     This subroutine computes the eigen values and eigen vectors for a given 
*     symmetric square matrix by calling the ScaLAPACK routine PDSYEVD. 
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
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.D0 )
*
*     .. I/O Parameters ..
      INTEGER        NIN, NSTDOUT
      PARAMETER     ( NIN =11, NSTDOUT = 6 )
*
*     .. Memory Parameters ..
      INTEGER       DBLESZ, INTGSZ, MEMSIZ, TOTMEM, MMAX, NMAX
      PARAMETER     ( DBLESZ = 8, INTGSZ = 4, TOTMEM = 128000000,
     $                MEMSIZ = TOTMEM / DBLESZ,
     $                MMAX = 7500, NMAX = 1102 )
*
*    .. Processor variables ..
      INTEGER        NP, NPROW, NPCOL, MB, NB, NB2
*
*    .. Local Scalars ..
      INTEGER        IAM, ICTXT, INFO, IPA, IPW,IPZ,IPWORK,IPIWORK,
     $               RET_SUBS,MYCOL, MYROW, M, N, NROWA, IZ, JZ, 
     $               IA, JA, LWORK, LIWORK, TRILWMIN, NCOLZ, TEMP,
     $               MYMEMSIZ, WORKSIZ, NOUTMAT
                        

*
*    .. Local Arrays ..
      INTEGER       DESCA ( DLEN_ ),  DESCZ (DLEN_), OUTDIM(3) 
      DOUBLE PRECISION   MEM ( MEMSIZ ) 
*
*    .. Local Scalars ..
*      INTEGER      IZ, JZ, IA, JA, NCOLZ, NPP, 
*     $              NQQ, TRILWMIN
*    ..
*    .. External Subroutines ..
      EXTERNAL      BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT, 
     $              BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $              DESCINIT, PDSYEVD, CRDistData, CRCollectData

*    ..
*    .. External Functions ..
      INTEGER       NUMROC, INDXL2G
      DOUBLE PRECISION   PDLANGE
      EXTERNAL      NUMROC, PDLANGE, INDXL2G
*
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*    ..

      
       RET_SUBS =1
       NOUTMAT =2

       M = PGINFO(1) 
       N = PGINFO(2) 
       MB = PGINFO (5) 
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
      NCOLA = NUMROC ( N, NB, MYCOL, 0, NPCOL )

      IZ = 1
      JZ = 1
      IA = 1 
      JA = 1 

      NCOLZ = NUMROC ( JZ + N - 1, NB, MYCOL, 0, NPCOL )
      TRILWMIN = 3*N + MAX( NB*( NROWA+1 ), 3*NB )
      LWORK = MAX( 1 + 6*N + 2*NROWA*NCOLA, TRILWMIN ) + 2*N
      LIWORK = 7*N + 8*NPCOL + 2

*
*    .. Initialize the array descriptors for A and Z ..
*
      CALL DESCINIT ( DESCA, M, N, MB, NB, 0, 0, ICTXT, MAX (1,NROWA),
     $              INFO)
      CALL DESCINIT ( DESCZ, M, N, MB, NB, 0, 0, ICTXT, MAX(1,NROWA),
     $              INFO)
*
*    .. Assign Pointers into MEM for SCALAPACK arrays, A is
*    .. allocated starting at position MEM (1)
*
      IPA = 1
      IPW = IPA + DESCA ( LLD_ ) * NROWA
      IPZ = IPW + N
      IPWORK = IPZ + DESCZ ( LLD_ ) * NCOLZ
      IPIWORK = IPWORK + LWORK
      IPEXWORK = IPIWORK + LIWORK

      WORKSIZ =  MAX( NB , NP )

      MYMEMSIZ = IPEXWORK + WORKSIZ

      IF (MYMEMSIZ.GT.MEMSIZ) THEN
           Print *, ' REQUIRED MEMORY COULD NOT BE ALLOCATED '
           GO TO 20
      ENDIF
      
*
*     Get the data from Parallel Agent  ...
      CALL CRDistData( MEM( IPA ), DESCA, MEM( IPEXWORK ) )
*
*    Call PDSYEVD function
*
      CALL PDSYEVD('V', 'U', N, MEM(IPA),1,1, DESCA, MEM(IPW), MEM(IPZ),
     $         IZ,JZ, DESCZ,MEM(IPWORK), LWORK, 
     $         MEM(IPIWORK), LIWORK, INFO)
*
*    Collect Result
*
      IF ( IAM.EQ.0) THEN

*         Send Parallel Agent 'number of output matrices'
          CALL CRSendIntToPA(NOUTMAT, 1, 202)

          DO TEMP = 1, N 
               RETVAL(TEMP) = MEM(IPW + TEMP-1)
          END DO

          OUTDIM(1) = 1
          OUTDIM(2) = 1
          OUTDIM(3) = N 
          CALL CRSendIntToPA(OUTDIM, 3 , 300)

          CALL CRSendDoubleToPA(RETVAL, N, 400)
          RET_SUBS = TEMP
      ENDIF

      CALL BLACS_BARRIER(ICTXT,'ALL')

     
      IF ( IAM.EQ.0) THEN
          OUTDIM(1) = 0
          OUTDIM(2) = M
          OUTDIM(3) = M
          CALL CRSendIntToPA(OUTDIM, 3 , 301 )
      ENDIF

      CALL CRCollectData( M, M, MEM( IPZ ), 1,1,DESCZ,
     $                MEM( IPEXWORK ) )

*
*    .. EXIT THE GRID
*
   10 CONTINUE

      CALL BLACS_GRIDEXIT( ICTXT )
*
   20 CONTINUE
*
*      CALL BLACS_EXIT( 0 )
*
      RETURN 
*
      END
