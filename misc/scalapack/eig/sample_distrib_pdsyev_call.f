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
      PARAMETER          ( LWORK = 264, MAXN = 100, ZERO = 0.0D+0 )
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

      INTEGER   PERM( MAXN ), INVPERM( MAXN )
      DOUBLE PRECISION VAL1( MAXN ), VAL2( MAXN )
      LOGICAL   TREATED( MAXN)
      INTEGER   ival
      LOGICAL   completed
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
      N = 10
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

        open(file="A10.txt", unit = 22)

        DO i=1,n
          DO j=1,n
            read(22,*) A0(i,j)
          ENDDO
        ENDDO
        close(22)

        DO i=1,n
          WRITE(*,*) '====', i, A0(I,1:N)
        ENDDO


      END IF

      CALL PDGEMR2D(N, N, A0, 1, 1, DESCA0, ABIS, 1, 1, DESCA, CONTEXT)

*
*     Build a matrix that you can create with
*     a one line matlab command:  hilb(n) + diag([1:-1/n:1/n])
*
      CALL PDLAMODHILB( N, A, 1, 1, DESCA, INFO )
      IF(IAM.EQ.3) THEN
        DO i = 1, n/2
          DO j = 1, n/2
            write(*,*) '***', i, j, A(i,j), ABIS(i,j)
          ENDDO
        ENDDO 
      END IF

*
*     Ask PDSYEV to compute the entire eigendecomposition
*
      CALL PDSYEV( 'V', 'U', N, A, 1, 1, DESCA, W, Z, 1, 1,
     $             DESCZ, WORK, LWORK, INFO )
      CALL PDSYEV( 'V', 'U', N, ABIS, 1, 1, DESCA, WBIS, ZBIS, 1, 1,
     $             DESCZ, WORK, LWORK, INFO )

      IF(IAM.EQ.0) THEN
        CALL DGEEV( 'N', 'V', N, A0, LDA, WRS, WIS, VLS, LDA, VRS,
     $         LDA, WORK, LWORK, INFO )

        write(*,*) 'Valeurs propres'
        DO i = 1, n
          write(*,*) 'W(', i, ') = ', W(i), WBIS(i), WRS(i)
        ENDDO
      END IF
      
*
*     Print out the eigenvectors
*
      OP(1) = 'Z'
      OP(2) = ' '
*      CALL PDLAPRNT( N, N, Z, 1, 1, DESCZ, 0, 0, OP, 6, WORK )
*      CALL PDLAPRNT( N, N, ZBIS, 1, 1, DESCZ, 0, 0, OP, 6, WORK )
      CALL PDGEMR2D(N, N, ZBIS, 1, 1, DESCA, 
     &              ZTER, 1, 1, DESCA0, CONTEXT)

*      IF(IAM.EQ.0) THEN
*        DO j = 1, n
*          DO i = 1, n
*            write(*,*) 'ZTER(', i, ',', j, ') = ', ZTER(i,j)
*          ENDDO
*        ENDDO 
*      END IF

      IF(IAM.EQ.0) THEN
        DO i = 1,n
          PERM(i) = i
        ENDDO
        DO i = 1, n-1
          DO j = i+1, n
            IF (W(i) < W(j)) THEN

              ival = PERM(i)
              PERM(i) = PERM(j)
              PERM(j) = ival

              val = W(i)
              W(i) = W(j)
              W(j) = val
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

        DO i = 1, n
          write(*,*) PERM(i)
          INVPERM(PERM(i)) = i
        ENDDO

        i = 1
        VAL1(:) = ZTER(:, PERM(i))
        TREATED(:) = .FALSE.
        completed = .FALSE.

 100    if(.not.completed) then

 200      if(.not.TREATED(PERM(i))) then
            TREATED(PERM(i)) = .TRUE.
            VAL2(:) = ZTER(:, i)
            ZTER(:, i) = VAL1(:)
            VAL1(:) = VAL2(:)
            i = INVPERM(i)
            goto 200
          end if

          j = 1
 300      if((j <= n) .and. TREATED(j)) then
            j = j +1
            goto 300
          end if

          if(j > n) then
            completed = .TRUE.
          else
            i = j
            VAL1(:) = ZTER(:, PERM(i))
          end if

          goto 100
        end if
        
        DO j = 1, n
          DO i = 1, n
            write(*,*) 'VP(', i, ',', j, ') = ', ZTER(i,j), VRS(i,j),
     &        abs(abs(ZTER(i,j)) - abs(VRS(i,j)))/abs(VRS(i,j))
          ENDDO
        ENDDO


      END IF
*
*     Print out matlab code which will check the residual
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         PRINT *, ' N =', N
         PRINT *, ' A =  hilb(N) + diag([1:-1/N:1/N])'
         DO 10 I = 1, N
            PRINT *, ' W(', I, ')=', W( I ), ';'
   10    CONTINUE
         PRINT *, ' backerror = A - Z * diag(W) * Z'' '
         PRINT *, ' resid = A * Z - Z * diag(W)'
         PRINT *, ' ortho = Z'' * Z - eye(N)'
         PRINT *, ' norm(backerror)'
         PRINT *, ' norm(resid)'
         PRINT *, ' norm(ortho)'
      END IF
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
*
      SUBROUTINE PDLAMODHILB( N, A, IA, JA, DESCA, INFO )
*
*  -- ScaLAPACK routine (version 1.2) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 10, 1996
*
*
*
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_, CTXT_, M_, N_,
     $                   MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,
     $                   CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                   RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PDELSET
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
*
*     The matlab code for a real matrix is:
*         hilb(n) + diag( [ 1:-1/n:1/n ] )
*     The matlab code for a complex matrix is:
*         hilb(N) + toeplitz( [ 1 (1:(N-1))*i ] )
*
*       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DT_*LLD_*MB_*M_*NB_*N_*
     $    RSRC_.LT.0 )RETURN
*
      INFO = 0
*
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
*
*
      IF( IA.NE.1 ) THEN
         INFO = -3
      ELSE IF( JA.NE.1 ) THEN
         INFO = -4
      END IF
*
      DO 20 J = 1, N
         DO 10 I = 1, N
            IF( I.EQ.J ) THEN
               CALL PDELSET( A, I, J, DESCA,
     $                       ( DBLE( N-I+1 ) ) / DBLE( N )+ONE /
     $                       ( DBLE( I+J )-ONE ) )
            ELSE
               CALL PDELSET( A, I, J, DESCA, ONE / ( DBLE( I+J )-ONE ) )
            END IF
   10    CONTINUE
   20 CONTINUE
*
*
      RETURN
*
*     End of PDLAMODHLIB
*
      END
