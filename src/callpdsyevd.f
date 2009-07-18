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
       
      SUBROUTINE CALLPDSYEVD (PGINFO, MEM, MEMSIZ)
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
      INTEGER DESCA (DLEN_), DESCZ(DLEN_),
     $             LLD_A, LLD_Z, INFO, NUMR_A, NUMC_A, NUMC_Z, LWORK,
     $             IZ, JZ, IA, JA, LIWORK, NP, NQ, TRILWMIN             
*     .. Memory Pointers ...
      INTEGER IPA, IPW, IPZ, IPWORK, IPIWORK, TEMP
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

      NOUTMAT = 2
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
      
      IZ = 1
      JZ = 1
      IA = 1
      JA = 1

      NP = NUMROC ( N, NB, MYROW, 0, NPROW )
      NQ = NUMROC ( N, NB, MYCOL, 0, NPCOL )
      TRILWMIN = 3 * N + MAX (NB * ( NP +1 ), 3 * NB)
      TEMP = MAX ( 1 + 6*N + 2*NP*NQ , TRILWMIN) + 2*N

      NUMC_Z = NUMROC ( JZ+N-1, NB, MYCOL, 0, NPCOL)
      LIWORK = 7*N + 8*NPCOL + 2


*     .. get LLD_A .. leading dimension of A
      LLD_A = MAX (1, NUMR_A)
*     .. get LLD_Z .. leading dimension of U 
      LLD_Z = MAX (1, NUMR_A)


      CALL DESCINIT ( DESCA, M, N, MB, NB, 0,0, ICTXT, LLD_A, INFO)
      
      CALL DESCINIT ( DESCZ, M, N, MB, NB, 0,0, ICTXT, LLD_Z, INFO)
      
*      PRINT *, 'DESC Initialization successful'
* ================================================================

*     Assign Pointers into MEM for scalapack arrays

      IPA = 1
      IPW = IPA + DESCA (LLD_) * NUMC_A
      IPZ = IPW + N
      IPIWORK = IPZ + DESCZ (LLD_) * NUMC_Z
      IPWORK = IPIWORK +  LIWORK

*      PRINT *, 'Rank=', IAM, ' IPA=', IPA, 
*     $     ' IPW=', IPW, ' IPZ=', IPZ, 
*     $     ' IPIWORK=', IPIWORK, ' IPWORK = ', IPWORK, 
*     $     ' LIWORK = ', LIWORK, '\n\n'
      

*      PRINT *, 'Pointer assigning successful'
* ================================================================

*     Get the optimum working size by doing a workspace query

      CALL PDSYEVD('V', 'U', N, MEM(IPA),1,1, DESCA, MEM(IPW), MEM(IPZ),
     $         IZ,JZ, DESCZ,MEM(IPWORK), -1, 
     $         MEM(IPIWORK), LIWORK, INFO)



      IF ( INFO.EQ.0 ) THEN
            LWORK = MEM ( IPWORK )

*            PRINT *, 'Rank=', IAM, ' LWORK=', LWORK, 
*     $      ' IPWORK=', IPWORK, 'MEMSIZ=', MEMSIZ,  '\n\n'

            IF ( LWORK.EQ.TEMP) THEN
*                 PRINT *, 'LWORK  EQ TEMP ', LWORK, TEMP
            ELSE
                 PRINT *, 'LWORK  Not EQ TEMP ', LWORK, TEMP
            ENDIF
            
            IF (IPWORK+LWORK-1 .GT. MEMSIZ ) THEN
                  PRINT *, 'RANK=', IAM, ' NOT ENOUGH MEMORY .. ',
     $             MEMSIZ, IPWORK + LWORK 
                  FAILFLAG = 1 
            ENDIF
      ELSE
            IF ( IAM.EQ.0)
     $             PRINT *, 'PDSYEVD DRY RUN FAILED .. 
     $                  = ', INFO
            FAILFLAG = 1     
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
*     Call PDSYEVD function

      CALL PDSYEVD('V', 'U', N, MEM(IPA),1,1, DESCA, MEM(IPW), MEM(IPZ),
     $         IZ,JZ, DESCZ,MEM(IPWORK), LWORK, 
     $         MEM(IPIWORK), LIWORK, INFO)
        
*      PRINT *, 'EIGEN computation  successful'
* ================================================================
*     Collect Result
*
      IF ( IAM.EQ.0) THEN

*         Send Parallel Agent 'number of output matrices'
          CALL CRSendIntToPA(NOUTMAT, 1, 202)

          OUTDIM(1) = 1
          OUTDIM(2) = 1
          OUTDIM(3) = N 
          CALL CRSendIntToPA(OUTDIM, 3 , 300)

          CALL CRSendDoubleToPA(MEM(IPW), N, 400)
      ENDIF

      CALL BLACS_BARRIER(ICTXT,'ALL')

     
      IF ( IAM.EQ.0) THEN
          OUTDIM(1) = 0
          OUTDIM(2) = M
          OUTDIM(3) = M
          CALL CRSendIntToPA(OUTDIM, 3 , 301 )
      ENDIF

      CALL CRCollectData( M, M, MEM( IPZ ), 1,1,DESCZ,
     $                MEM( IPWORK ) )


*      PRINT *, 'EIGEN result collection  successful'
* ==================================================================
*     Exit the Grid
   20 CONTINUE
      
      CALL BLACS_GRIDEXIT ( ICTXT)

*      PRINT *, 'Exiting BLACS grid  successful'

      RETURN
*
      END
