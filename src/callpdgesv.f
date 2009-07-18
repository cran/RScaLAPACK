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
       
      SUBROUTINE CALLPDGESV (PGINFO, MEM, MEMSIZ)
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
      INTEGER M_1, N_1, M_2, N_2, MB, NB, NPROW, NPCOL, IAM, 
     $                             NPROCS, ICTXT 
*     .. Local information ...
      INTEGER MYROW, MYCOL
*     .. Local Descriptors ...
      INTEGER DESCA (DLEN_), DESCB(DLEN_),
     $             LLD_A, LLD_B, INFO, NUMR_A, NUMC_A, NUMR_B,NUMC_B
*     .. Memory Pointers ...
      INTEGER IPA, IPB, IPPIV, IPWORK, WORKSIZ
*     .. Output Params
      INTEGER FAILFLAG,NOUTMAT,OUTDIM(3)

* =============================================================
*     Function Declarations

*     .. External Functions ..
      INTEGER NUMROC
      EXTERNAL NUMROC

*==============================================================
*     Assign process grid values

      M_1 = PGINFO (1)
      N_1 = PGINFO (2)
      M_2 = PGINFO (3)
      N_2 = PGINFO (4)
      MB = PGINFO (5)
      NB = PGINFO (6)
      NPROW = PGINFO (7)
      NPCOL = PGINFO (8)

      NOUTMAT = 1 
      FAILFLAG = 0

*      PRINT *, 'Initialization successful'
*==============================================================
*     Set up processes and memory

      CALL BLACS_PINFO ( IAM, NPROCS )

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
      NUMR_A =  NUMROC (M_1, MB, MYROW, 0, NPROW)    
      NUMC_A =  NUMROC (N_1, NB, MYCOL, 0, NPCOL) 
      

      NUMR_B =  NUMROC (M_2, MB, MYROW, 0, NPROW)    
      NUMC_B =  NUMROC (N_2, NB, MYCOL, 0, NPCOL) 

*     .. get LLD_A .. leading dimension of A
      LLD_A = MAX (1, NUMR_A)
*     .. get LLD_Z .. leading dimension of U 
      LLD_B = MAX (1, NUMR_B)

      CALL DESCINIT ( DESCA, M_1, N_1, MB, NB, 0,0, ICTXT, LLD_A, INFO)
      
      CALL DESCINIT ( DESCB, M_2, N_2, MB, NB, 0,0, ICTXT, LLD_B, INFO)
      
*      PRINT *, 'DESC Initialization successful'
* ================================================================

*     Assign Pointers into MEM for scalapack arrays

      IPA = 1
      IPPIV = IPA + DESCA(LLD_) * NUMC_A
      IPB = IPPIV + NUMR_A + MB
      IPWORK = IPB + DESCB(LLD_) * NUMC_B
      
      WORKSIZ = NB
      
*      PRINT *, 'Pointer assigning Initialization successful'
* ================================================================

*     Get the optimum working size by doing a workspace query

      IF (IPWORK+WORKSIZ-1 .GT. MEMSIZ ) THEN
               PRINT *, 'NOT ENOUGH MEMORY .. ', MEMSIZ 
     $              , IPWORK + WORKSIZ 
             FAILFLAG = 1
            GO TO 20
      ENDIF

*      PRINT *, 'Check Fail Flag'

      CALL CRCheckFailFlag (FAILFLAG)

      IF ( IAM.EQ.0)
     $      CALL CRSendIntToPA( FAILFLAG, 1, 1202 )

      IF ( FAILFLAG.EQ.1 )
     $     GO TO 20


*      PRINT *, 'Memory verification successful'
* =================================================================

*     Distribute the Input Matrix

      CALL CRDistData ( MEM(IPA), DESCA, MEM(IPWORK))

      CALL CRDistData ( MEM(IPB), DESCB, MEM(IPWORK))

*      PRINT *, 'Distribution  successful'
* =================================================================
*     Call PDGESV function

      CALL PDGESV( M_1, N_2, MEM( IPA ), 1, 1, DESCA, MEM( IPPIV ),
     $             MEM( IPB ), 1, 1, DESCB, INFO )

*      PRINT *, 'SOLVE computation  successful'
* ================================================================
*     Collect Result
*
      IF (IAM.EQ.0) THEN
          CALL CRSendIntToPA(NOUTMAT, 1, 202)
          OUTDIM(1) = 0
          OUTDIM(2) = M_1
          OUTDIM(3) = N_2

          CALL CRSendIntToPA(OUTDIM, 3 , 300)
      ENDIF

      CALL CRCollectData(M_1, N_2, MEM( IPB ), 1, 1, DESCB,
     $                MEM( IPWORK ) )
*
*      PRINT *, 'SOLVE result collection  successful'
* ==================================================================
*     Exit the Grid
   20 CONTINUE

      CALL BLACS_GRIDEXIT ( ICTXT)

*      PRINT *, 'Exiting BLACS grid  successful'

      RETURN
*
      END
