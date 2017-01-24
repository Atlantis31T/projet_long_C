*
*
      PROGRAM SAMPLE_PDSYEV_CALL
*
*
*  -- ScaLAPACK routine (version 1.2) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*
*     This routine contains a sample call to PDSYEV.
*     When compiled and run, it produces output which can be
*     pasted directly into matlab.
*
      IMPLICIT NONE
*     .. Parameters ..
      INTEGER            LWORK, MAXN
      DOUBLE PRECISION   ZERO
      PARAMETER          ( LWORK = 50000, MAXN = 400, ZERO = 0.0D+0 )
      INTEGER            LDA
      DOUBLE PRECISION   MONE
      INTEGER            MAXPROCS
      PARAMETER          ( LDA = MAXN, MONE = -1.0D+0, MAXPROCS = 512 )
*     ..
*     .. Local Scalars ..
      INTEGER            CONTEXT, I, IAM, INFO, MYCOL, MYROW, N, NB,
     $                   NPCOL, NPROCS, NPROW
      INTEGER            CONTEXT_SYS, CONTEXT0
      INTEGER            J, K
      DOUBLE PRECISION   VAL
*     ..
*     .. Local Arrays ..
      INTEGER            DESCA( 50 ), DESCZ( 50 )
      INTEGER            DESCA0( 50 )
      DOUBLE PRECISION   A( LDA, LDA ), W( MAXN ), WBIS( MAXN ),
     $                   WORK( LWORK ), Z( LDA, LDA )
      DOUBLE PRECISION   A0( LDA, LDA ), ABIS( LDA, LDA ),
     $                   ZBIS( LDA, LDA ), ZTER( LDA, LDA)
      DOUBLE PRECISION   WRS(MAXN), WIS(MAXN)
      DOUBLE PRECISION   VLS(LDA, LDA), VRS(LDA,LDA)
      CHARACTER          OP(2)
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,
     $                   BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,
     $                   BLACS_SETUP, DESCINIT, PDLAMODHILB, PDLAPRNT,
     $                   PDSYEV
*     ..
*     .. Executable Statements ..
*
*
*     Set up the problem
*
      N = 356
      NB = 1
      NPROW = 2
      NPCOL = 2
*
*
*     Initialize the BLACS
*
      CALL BLACS_PINFO( IAM, NPROCS )
      IF( ( NPROCS.LT.1 ) ) THEN
         CALL BLACS_SETUP( IAM, NPROW*NPCOL )
      END IF
*
*
*     Initialize a single BLACS context
*
      CALL BLACS_GET( 0, 0, CONTEXT_SYS )
      CONTEXT = CONTEXT_SYS
      CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
      WRITE(*,*) 'CONTEXT', IAM, MYROW, MYCOL
*

*     Bail out if this process is not a part of this context.
*
      IF( MYROW.EQ.-1 )
     $   GO TO 20
*
*
*     These are basic array descriptors
*
      CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
      CALL DESCINIT( DESCZ, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )

* Set up a process grid of size 1 with context ctxt_0
      CONTEXT0 = CONTEXT_SYS
      CALL BLACS_GRIDINIT( CONTEXT0, 'R', 1, 1)
      CALL BLACS_GRIDINFO( CONTEXT0, NPROW, NPCOL, MYROW, MYCOL )
      WRITE(*,*) 'CONTEXT0', IAM, MYROW, MYCOL
* initialize the descriptor desca_0 for the global nxn matrix
      DESCA0(2) = -1
      if (IAM.eq.0) THEN
        call DESCINIT( DESCA0, N, N, N, N, 0, 0, CONTEXT0, LDA, INFO )

        open(file="Asp3.txt", unit = 22)

        DO i=1,n
          DO j=1,n
            read(22,*) A0(i,j)
          ENDDO
        ENDDO
        close(22)

        !DO i=1,n
        !  WRITE(*,*) '====', i, A0(I,1:N)
        !ENDDO


      END IF

      CALL PDGEMR2D(N, N, A0, 1, 1, DESCA0, ABIS, 1, 1, DESCA, CONTEXT)

*
*     Ask PDSYEV to compute the entire eigendecomposition
*
      CALL PDSYEV( 'V', 'U', N, ABIS, 1, 1, DESCA, WBIS, ZBIS, 1, 1,
     $             DESCZ, WORK, LWORK, INFO )

      IF(IAM.EQ.0) THEN
        CALL DGEEV( 'N', 'V', N, A0, LDA, WRS, WIS, VLS, LDA, VRS,
     $         LDA, WORK, LWORK, INFO )

      END IF
      
*
*     Print out the eigenvectors
*
      OP(1) = 'Z'
      OP(2) = ' '
!      CALL PDLAPRNT( N, N, ZBIS, 1, 1, DESCZ, 0, 0, OP, 6, WORK )
      CALL PDGEMR2D(N, N, ZBIS, 1, 1, DESCA, 
     &              ZTER, 1, 1, DESCA0, CONTEXT)

      IF(IAM.EQ.0) THEN
        DO j = 1, n
          DO i = 1, n
!            write(*,*) 'ZTER(', i, ',', j, ') = ', ZTER(i,j)
          ENDDO
        ENDDO 
      END IF

      IF(IAM.EQ.0) THEN
        DO i = 1, n-1
          DO j = i+1, n
            IF (WBIS(i) < WBIS(j)) THEN
              val = WBIS(i)
              WBIS(i) = WBIS(j)
              WBIS(j) = val
              DO k = 1, n
               val = ZTER(k,i)
               ZTER(k,i) = ZTER(k,j)
               ZTER(k,j) = val
              ENDDO
            END IF
            IF (WRS(i) < WRS(j)) THEN
              val = WRS(i)
              WRS(i) = WRS(j)
              WRS(j) = val
              DO k = 1, n
                val = VRS(k,i)
                VRS(k,i) = VRS(k,j)
                VRS(k,j) = val
              ENDDO
            END IF
          ENDDO
        ENDDO

        write(*,*) 'Valeurs propres'
        DO i = 1, n
          write(*,*) 'W(', i, ') = ', WRS(i), WBIS(i)
        ENDDO
        DO j = 1, n
*          write(23,*) 'W(', j, ') = ', WRS(j), WBIS(j)
          DO i = 1, n
*            write(23,*) 'VP(', i, ',', j, ') = ', VRS(i,j), ZTER(i,j),
*     &            VRS(i,j)/ZTER(i,j)
            write(23,*) VRS(i,j), ZTER(i,j)
          ENDDO
        ENDDO

      END IF
*
*     Print out matlab code which will check the residual
*
      CALL BLACS_GRIDEXIT( CONTEXT )
*
   20 CONTINUE
*
      CALL BLACS_EXIT( 0 )
*
*
*     Uncomment this line on SUN systems to avoid the useless print out
*
*      CALL IEEE_FLAGS( 'clear', 'exception', 'underflow', '')
*
*
 9999 FORMAT( 'W=diag([', 4D16.12, ']);' )
*
      STOP
      END
