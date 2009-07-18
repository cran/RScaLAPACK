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
       
      SUBROUTINE CALLPDGEQRF (PGINFO, MEM, MEMSIZ)
*
      INTEGER PGINFO(9)
      DOUBLE PRECISION MEM (*)
      INTEGER MEMSIZ

* ==============================================================
*     Variable Declaration

*     .. Array Descriptor Parameters ..
      INTEGER BLOCK_CYCLIC_2D, DLEN_, DT_, CTXT_, M_,N_,
     $             MB_, NB_, RSRC_, CSRC_, LLD_

      PARAMETER   (BLOCK_CYCLIC_2D=1, DLEN_=9, DT_=1, CTXT_=2,
     $             M_=3, N_=4, MB_=5, NB_=6, RSRC_=7, CSRC_=8,
     $             LLD_=9)
      
*     .. Processor variables ..
      INTEGER M, N, MB, NB, NPROW, NPCOL, IAM, NP, ICTXT 
*     .. Local information ...
      INTEGER MYROW, MYCOL
*     .. Local Descriptors ...
      INTEGER DESCA (DLEN_),DESCTAU(DLEN_), LLD_TAU, 
     $          LLD_A, INFO, NUMR_A, NUMC_A, NUMR_TAU, NUMC_TAU
*     .. Memory Pointers ...
      INTEGER IPA, IPTAU, IPWORK, LWORK, LTAU
*     .. Output Params
      INTEGER I,J,K,L,FAILFLAG,NOUTMAT,OUTDIM(3), RDIM, ITER
      INTEGER ICALLER, HISROW, HISCOL
      DOUBLE PRECISION RANK

* =============================================================
*     Function Declarations

*     .. External Functions ..
      INTEGER NUMROC
      EXTERNAL NUMROC

*==============================================================
*     Assign process grid values

      M = PGINFO (1)
      N = PGINFO (2)
      MB = PGINFO (5)
      NB = PGINFO (6)
      NPROW = PGINFO (7)
      NPCOL = PGINFO (8)

      NOUTMAT = 3
      FAILFLAG = 0 

*      PRINT *, 'Initialization successful'
*==============================================================
*     Set up processes and memory

      CALL BLACS_PINFO ( IAM, NP )

*     Define Process Grid
      CALL BLACS_GET ( -1, 0, ICTXT ) 
      CALL BLACS_GRIDINIT ( ICTXT, 'Row-major', NPROW, NPCOL )
      CALL BLACS_GRIDINFO ( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*      PRINT *, 'BLACS Initialization successful'

*     Do not perform anything is not part of process grid
      IF ( MYROW.GE.NPROW .OR. MYCOl.GE.NPCOL)
     $    GO TO 20
* =============================================================
*     Initialize the array descriptors

*     Leading dimension of A
      NUMR_A =  NUMROC (M, MB, MYROW, 0, NPROW)    
      NUMC_A =  NUMROC (N, NB, MYCOL, 0, NPCOL) 

      NUMR_TAU =  NUMROC (1,MB,MYROW,0,NPROW)    
      NUMC_TAU =  NUMROC (MIN(M,N), NB, MYCOL, 0, NPCOL) 

*     .. get LLD_A .. leading dimension of A
      LLD_A = MAX (1, NUMR_A)
*     .. get LLD_TAU .. leading dimension of A
      LLD_TAU = MAX (1, NUMR_TAU)

      CALL DESCINIT ( DESCA, M, N, MB, NB, 0,0, ICTXT, LLD_A, INFO)
      CALL DESCINIT ( DESCTAU,1,MIN(M,N),MB,NB,0,0,ICTXT,LLD_TAU,INFO)

*      PRINT *, 'DESC Initialization successful'
* ================================================================
      IA = 1
      JA = 1

*      LWORK = NB * (MP0 + NQ0 + NB)
      LTAU = NUMROC(JA+MIN(M,N)-1, NB, MYCOL, DESCA(CSRC_),
     $                           NPCOL)
 
*     Assign Pointers into MEM for scalapack arrays

      IPA = 1
      IPTAU = IPA + DESCA(LLD_) *  NUMC_A 
      IPWORK = IPTAU + LTAU
      

*      PRINT *, 'Pointer assigning successful'
* ================================================================

*     Get the optimum working size by doing a workspace query
*
      CALL PDGEQRF( M, N, MEM( IPA ), 1, 1, DESCA, MEM(IPTAU),
     $             MEM( IPWORK ), -1 , INFO )


      IF ( INFO.EQ.0 ) THEN
            LWORK = MEM ( IPWORK )
*            PRINT *, 'LWORK  = ', LWORK
            
            IF (IPWORK+LWORK-1 .GT. MEMSIZ ) THEN
                  PRINT *, 'NOT ENOUGH MEMORY ..myRank',IAM, MEMSIZ 
     $             , IPWORK + LWORK
                   FAILFLAG = 1
*                  GO TO 20
            ENDIF 
      ELSE
            IF ( IAM.EQ.0)
     $             PRINT *, 'PDGEQRF DRY RUN FAILED .. 
     $                  = ', INFO
             FAILFLAG = 1
*            GO TO 20
      ENDIF

*      PRINT *, 'Check Fail Flag'

      CALL CRCheckFailFlag (FAILFLAG)

      IF ( IAM.EQ.0)
     $      CALL CRSendIntToPA( FAILFLAG, 1, 1202 )

      IF ( FAILFLAG.EQ.1 )
     $     GO TO 20



*      PRINT *, 'Dry run successful'
* =================================================================

*     Distribute the Input Matrix

      CALL CRDistData ( MEM(IPA), DESCA, MEM(IPWORK))

*      PRINT *, 'Distribution  successful'
* =================================================================
*     Call PDGEQRF function

      CALL PDGEQRF( M, N, MEM( IPA ), 1, 1, DESCA, MEM(IPTAU),
     $             MEM( IPWORK ), LWORK , INFO )

        
*      PRINT *, 'QR computation  successful'
* ================================================================
*     Collect Result
*
      IF (IAM.EQ.0) THEN
           CALL CRSendIntToPA(NOUTMAT, 1, 202)
      ENDIF
*
      IF (IAM.EQ.0) THEN
          OUTDIM(1) = 0
          OUTDIM(2) = M
          OUTDIM(3) = N
          CALL CRSendIntToPA(OUTDIM, 3 , 300 )
      ENDIF

      CALL CRCollectData(M, N, MEM( IPA ), 1, 1, DESCA,
     $                MEM( IPWORK ) )

      CALL BLACS_BARRIER (ICTXT, 'ALL')

*
*     SEND THE RANK
*

      IF (IAM.EQ.0) THEN

          OUTDIM(1) = 1
          OUTDIM(2) = 1
          OUTDIM(3) = 1
          RANK = MIN (M, N) 
          CALL CRSendIntToPA(OUTDIM, 3 , 301)
          CALL CRSendDoubleToPA(RANK, 1, 401)
      ENDIF

*
*     SEND TAU MATRIX      
*

      OUTDIM(1) = 0
      RDIM = MIN (M,N)
      OUTDIM(2) = 1 
      OUTDIM(3) = RDIM 

      CALL CRSendIntToPA(OUTDIM, 3 , 302)

      CALL CRCollectData(1,RDIM, MEM( IPTAU ), 1, 1, DESCTAU,
     $                MEM( IPWORK ) )
*      ENDIF

*      DO 15 I= 1, MIN(M,N)
*          PRINT *, 'TAU VALUES  .... RANK ',MEM(IPTAU+I-1),IAM   
*   15 CONTINUE

*      PRINT *, 'QR Decomposition result collection  successful'
* ==================================================================
*     Exit the Grid
   20 CONTINUE

      CALL BLACS_GRIDEXIT ( ICTXT)

*      PRINT *, 'Exiting BLACS grid  successful'

      RETURN
*
      END
