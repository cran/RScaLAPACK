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
       
      SUBROUTINE CALLPDPOTRF (PGINFO, MEM, MEMSIZ)
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
      INTEGER M, N, MB, NB, NPROW, NPCOL, IAM, NPROCS, ICTXT 
*     .. Local information ...
      INTEGER MYROW, MYCOL
*     .. Local Descriptors ...
      INTEGER DESCA (DLEN_), LLD_A, INFO, NUMR_A, NUMC_A
*     .. Memory Pointers ...
      INTEGER IPA, IPW, WORKSIZ
*     .. Output Params
      INTEGER FAILFLAG,NOUTMAT,OUTDIM(3)

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
      NUMR_A =  NUMROC (M, MB, MYROW, 0, NPROW)    
      NUMC_A =  NUMROC (N, NB, MYCOL, 0, NPCOL) 
      

*     .. get LLD_A .. leading dimension of A
      LLD_A = MAX (1, NUMR_A)


      CALL DESCINIT ( DESCA, M, N, MB, NB, 0,0, ICTXT, LLD_A, INFO)
      
      
*      PRINT *, 'DESC Initialization successful'
* ================================================================

*     Assign Pointers into MEM for scalapack arrays

      IPA = 1
      IPW = IPA + DESCA(LLD_) * NUMC_A

      WORKSIZ = MAX ( MB, NB)
      
*      PRINT *, 'Pointer assigning Initialization successful'
* ================================================================
*     Get the optimum working size by doing a workspace query

      IF ( IPW + WORKSIZ -1 .GT.MEMSIZ ) THEN
               PRINT *, 'NOT ENOUGH MEMORY .. Exiting ..', MEMSIZ 
     $             , IPW + WORKSIZ 
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

      CALL CRDistData ( MEM(IPA), DESCA, MEM(IPW))

*      PRINT *, 'Distribution  successful'
* =================================================================
*     Call PDPOTRF function

      CALL PDPOTRF ('U',M,MEM (IPA),1,1,DESCA,INFO)

*      PRINT *, 'CHOL computation  successful'
* ================================================================
*     Collect Result
*
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

*      PRINT *, 'CHOL result collection  successful'
* ==================================================================
*     Exit the Grid
   20 CONTINUE

      CALL BLACS_GRIDEXIT ( ICTXT)

*      PRINT *, 'Exiting BLACS grid  successful'

      RETURN
*
      END
