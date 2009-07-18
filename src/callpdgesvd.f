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
       
      SUBROUTINE CALLPDGESVD (PGINFO, MEM, MEMSIZ)
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
      INTEGER DESCA (DLEN_), DESCU (DLEN_), DESCVT (DLEN_),
     $             LLD_A, LLD_U, LLD_VT, INFO, NUMR_A, NUMC_A,
     $             TSIZE, SIZEQ, SIZEP
*     .. Memory Pointers ...
      INTEGER IPA, IPU, IPS, IPVT, IPW, WORKSIZ
*     .. Output Params
      INTEGER I,J,K,L,FAILFLAG,NOUTMAT,OUTDIM(3)

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

*     SVD related
      NU = PGINFO (3)
      NV = PGINFO (4)

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
      NUMC_A =  NUMROC (N, MB, MYCOL, 0, NPCOL)


*     .. get LLD_A .. leading dimension of A
      LLD_A = MAX (1, NUMR_A)
*     .. get LLD_U .. leading dimension of U 
      LLD_U = MAX (1, NUMR_A)
*     .. get LLD_VT ..leading dimension of VT
      LLD_VT = MAX (1, NUMROC (N, MB, MYROW, 0, NPROW))

      TSIZE = MIN (M,N)

      CALL DESCINIT ( DESCA, M, N, MB, NB, 0,0, ICTXT, LLD_A, INFO)
      
      CALL DESCINIT ( DESCU, M, TSIZE, MB, NB, 0,0, ICTXT, LLD_U, INFO)
      
      CALL DESCINIT ( DESCVT,TSIZE, N, MB, NB, 0,0, ICTXT, LLD_VT,INFO)


      SIZEQ = NUMROC (TSIZE, DESCU(NB_), MYCOL, DESCU(CSRC_), NPCOL)
      SIZEP = NUMROC (TSIZE, DESCVT(MB_), MYROW, DESCVT(RSRC_),NPROW)
         

*      PRINT *, 'DESC Initialization successful'
* ================================================================

*     Assign Pointers into MEM for scalapack arrays

      IPA = 1
      IPS = IPA + NUMR_A * NUMC_A
      IPU = IPS + N
      IPVT = IPU + NUMR_A * SIZEQ
      IPW = IPVT + SIZEP * NUMC_A

*      PRINT *, 'Pointer assigning Initialization successful'
* ================================================================

*     Get the optimum working size by doing a workspace query

      CALL PDGESVD ('V','V', M, N, MEM (IPA), 1, 1, DESCA, MEM (IPS),
     $                   MEM(IPU), 1,1, DESCU, MEM(IPVT), 1,1, DESCVT, 
     $                   MEM(IPW), -1, INFO)

      IF ( INFO.EQ.0 ) THEN
            WORKSIZ = MEM(IPW)
            IF (MEMSIZ.LT.IPW+WORKSIZ-1) THEN
                  PRINT*,'NOT ENOUGH MEMORY.. MEMSIZ.LT.IPW+WORKSIZ-1',
     $                             MEMSIZ, IPW+WORKSIZ-1
*                   GOTO 20
                  FAILFLAG = 1  
            ENDIF
      ELSE
            IF ( IAM.EQ.0)
     $             PRINT *, 'PDGESVD DRY RUN FAILED .. NOT ENOUGH MEMORY
     $                  = ', INFO
*            GO TO 20
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

      CALL CRDistData ( MEM(IPA), DESCA, MEM(IPW))

*      PRINT *, 'Distribution  successful'
* =================================================================
*     Call PDGESVD function

      IF ((NU.NE.0).AND.(NV.NE.0)) THEN

            CALL PDGESVD('V','V', M, N, MEM(IPA), 1, 1, DESCA, MEM(IPS),
     $            MEM(IPU), 1, 1, DESCU, MEM(IPVT), 1, 1, DESCVT,
     $            MEM(IPW), WORKSIZ, INFO)
      
      ELSE IF ((NU.EQ.0).AND.(NV.NE.0)) THEN
     
            CALL PDGESVD('N','V', M, N, MEM(IPA), 1, 1, DESCA, MEM(IPS),
     $            MEM(IPU), 1, 1, DESCU, MEM(IPVT), 1, 1, DESCVT,
     $            MEM(IPW), WORKSIZ, INFO)

        
      ELSE IF ((NU.NE.0).AND.(NV.EQ.0)) THEN
     
            CALL PDGESVD('V','N', M, N, MEM(IPA), 1, 1, DESCA, MEM(IPS),
     $            MEM(IPU), 1, 1, DESCU, MEM(IPVT), 1, 1, DESCVT,
     $            MEM(IPW), WORKSIZ, INFO)

      ELSE 
     
            CALL PDGESVD('N','N', M, N, MEM(IPA), 1, 1, DESCA, MEM(IPS),
     $            MEM(IPU), 1, 1, DESCU, MEM(IPVT), 1, 1, DESCVT,
     $            MEM(IPW), WORKSIZ, INFO)

      ENDIF

*      PRINT *, 'SVD computation  successful'
* ================================================================
*     Collect Result
*
      
      IF ( IAM.EQ.0) THEN
*     Send PA number of matrices
            CALL CRSendIntToPA (NOUTMAT, 1, 202) 
      
            OUTDIM(1) = 1
            OUTDIM(2) = 1
            OUTDIM(3) = DESCU(N_)

            CALL CRSendIntToPA (OUTDIM, 3, 300)
            CALL CRSendDoubleToPA (MEM(IPS), OUTDIM(3), 400)
      ENDIF

      CALL BLACS_BARRIER ( ICTXT,'ALL')
*      PRINT *, 'After sending singular vals ... Barrier '
*      PRINT *, 'NU ... NV ',NU,NV

      IF ( M.LE.N ) THEN

*         If there is a matrix U, then send it.
          IF ( NU.NE.0 ) THEN
               OUTDIM(1) = 0
               OUTDIM(2) = M
               OUTDIM(3) = M
               IF(IAM.EQ.0)
     $              CALL CRSendIntToPA(OUTDIM, 3 , 301 )

*               PRINT *, 'in nu.ne.0 .. when m.le.n '
              CALL CRCollectData( M,M,MEM( IPU ), 1,1,
     $                DESCU, MEM( IPW ) )

*         Else, send the dimenions [0, 0] to signal a NULL matrix.
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0

*               PRINT *, 'in nu.eq.0 .. when m.le.n '

              IF(IAM.EQ.0)
     $              CALL CRSendIntToPA(OUTDIM, 3 , 301 )

          ENDIF


*         If there is a matrix Vt, then send it.
          IF ( NV.NE.0 ) THEN
              IF (M.EQ.N) THEN
                    OUTDIM(1) = 0
                    OUTDIM(2) = N
                    OUTDIM(3) = N
                    IF(IAM.EQ.0)
     $                   CALL CRSendIntToPA(OUTDIM, 3 , 302 )

*                  PRINT *, 'in nv.ne.0 .. when m.eq.n '
                  CALL CRCollectData( N,N,MEM( IPVT ),
     $                              1, 1, DESCVT,MEM( IPW ) )
              ELSE
                    OUTDIM(1) = 0
                    OUTDIM(2) = M
                    OUTDIM(3) = N

                    IF(IAM.EQ.0)
     $                 CALL CRSendIntToPA(OUTDIM, 3 , 302 )

*                  PRINT *, 'in nv.ne.0 .. when m.lt.n '
                   CALL CRCollectData(M,N, MEM( IPVT ),
     $                               1, 1, DESCVT, MEM( IPW ) )
              ENDIF
*         Else, send the dimensions [0, 0] to signal a NULL matrix.
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
*              PRINT *, 'in nv.eq.0 .. when m.le.n '
              IF(IAM.EQ.0)
     $                 CALL CRSendIntToPA (OUTDIM, 3 , 302 )
          ENDIF


      ELSE

          IF ( NU.NE.0 ) THEN
               OUTDIM(1) = 0
               OUTDIM(2) = M
               OUTDIM(3) = N

               IF(IAM.EQ.0)
     $              CALL CRSendIntToPA(OUTDIM, 3 , 301 )

*               PRINT *, 'in nu.ne.0 .. when m.gt.n '
               CALL CRCollectData( M,N, MEM( IPU ), 1,1,
     $                          DESCU, MEM( IPW ) )
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
               IF(IAM.EQ.0)
     $              CALL CRSendIntToPA(OUTDIM, 3 , 301 )
*               PRINT *, 'in nu.eq.0 .. when m.gt.n '
          ENDIF

          IF ( NV.NE.0 ) THEN
               OUTDIM(1) = 0
               OUTDIM(2) = N
               OUTDIM(3) = N

               IF (IAM.EQ.0)
     $               CALL CRSendIntToPA(OUTDIM, 3 , 302 )
*               PRINT *, 'in nv.ne.0 .. when m.gt.n '

               CALL CRCollectData( N,N,MEM( IPVT ),1,1,
     $                           DESCVT, MEM( IPW ) )
          ELSE
              OUTDIM(1) = 1
              OUTDIM(2) = 0
              OUTDIM(3) = 0
*              PRINT *, 'in nv.eq.0 .. when m.gt.n '
              IF (IAM.EQ.0)
     $               CALL CRSendIntToPA(OUTDIM, 3 , 302 )
          ENDIF

      ENDIF

*      PRINT *, 'SVD result collection  successful'
* ==================================================================
*     Exit the Grid
   20 CONTINUE

      CALL BLACS_GRIDEXIT ( ICTXT)

*      PRINT *, 'Exiting BLACS grid  successful'

      RETURN
*
      END
