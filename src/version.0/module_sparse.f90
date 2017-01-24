module module_sparse
  use module_structure
  use module_solve
  use module_embed
contains

  !*****************************************
  !calcul des clusters
  subroutine sp_calculclusters(numproc,nblimit,nbideal,dataw,sigma)

    implicit integer(i,j,q)
    include 'mpif.h'
    type(type_data) :: dataw
    integer :: numproc,nbproc
    double precision :: sigma
    double precision,dimension(:,:),pointer :: cluster_center
    integer :: n, k, nbcluster
    double precision,dimension(:),pointer :: ratiomax,cluster_energy,&
         ratiomin,ratiomoy,ratiorii,ratiorij
    integer,dimension(:),pointer ::cluster,cluster_population,nbinfo
    integer :: nblimit,nbideal
    double precision :: norme,ratio,ratio1,ratio2,seuilrij
    character*30 :: num,files

! sparsification deb
    !double precision :: value
    double precision :: t1, t2, t_cons_a, t_cons_vp
    integer :: nnz, nnz2, nb
    double precision :: facteur
    integer :: l
    double precision :: treshold
    double precision,dimension(:),pointer :: AS
    integer, dimension(:),pointer :: IAS, JAS
    double precision, dimension(:),pointer :: D
    double precision,dimension(:,:),pointer :: Z
    double precision, dimension(:), pointer :: W

! sparsification fin

    !creation de la matrice
#if aff
    print *,numproc,'valeur du sigma',sigma
#endif
    n=dataw%nb

! sparsification deb
    nnz = 0
    ! valeur de treshold arbitraire -> paramètre du sp ou calcul interne
    !                                  (voir avec S.)
    ! TODO : mettre la valeur du facteur dans le fichier param
    facteur = 3.0
    treshold = facteur*sigma
    !write(*,*) '************** SIGMA', sigma

    t1 = MPI_WTIME();
    do i=1,n-1  ! borne ?
       do j=i+1,n ! borne ?

          norme=0.0

          do k=1,dataw%dim
             norme=norme+(dataw%point(i)%coord(k)-dataw%point(j)%coord(k))**2
          end do

          if(sqrt(norme) <= treshold) then
            nnz = nnz + 1;
          end if

       end do
    end do

    t2 = MPI_WTIME();
    t_cons_a = t2 - t1
    print *, numproc, 'surcoût A', t_cons_a

    t1 = MPI_WTIME();
    nnz2 = nnz*2

    allocate(AS(nnz2))
    allocate(IAS(nnz2))
    allocate(JAS(nnz2))
    l = 1
    do i=1,n-1
       do j=i+1,n
          norme=0.0
          do k=1,dataw%dim
             norme=norme+(dataw%point(i)%coord(k)-dataw%point(j)%coord(k))**2
          end do
          value=exp(-norme/sigma)
          ! on garde si value <= treshold
          ! (si on veut tout garder, commenter ligne if, end if)
          if(sqrt(norme) <= treshold) then
            AS(l) = value
            IAS(l) = i
            JAS(l) = j
            l = l+1
            !------
            AS(l) = value
            IAS(l) = j
            JAS(l) = i
            l = l+1
            !------
            ! autre rangement
            !AS(l+nnz) = value
            !IAS(l+nnz) = j
            !JAS(l+nnz) = i
            !l = l+1
          end if
       end do
    end do
    write(*,*) '========== facteur, n*n nnz2 = ', facteur, n*n, nnz2

  allocate(D(n)); D(:)=0.0
  do l=1, nnz2
    D(IAS(l)) = D(IAS(l)) + AS(l)
  end do

  do l=1, nnz2
    AS(l)=AS(l)/D(IAS(l))
  end do

  deallocate(D)

!  print *, "******************* A sparse ********************"
!  do l=1,nnz2
!    write(13,*)  IAS(l), JAS(l), AS(l)
!  end do

    ! nb et nblimit même valeur ?
    nb = 2*nblimit

    t1 = MPI_WTIME()
    call solve_arpack(AS, IAS, JAS, n, nnz2, nb, W, Z)
    print *, "---------- W -------------"
    do i=1,nb
       print *,'valeurs propres arpack brutes',i, W(i)
    end do
    

    !print *,numproc,'reordonne les vp...'  ! QUESTION nécessaire avec arpack ?
    do i=1,nb-1
       do j=i+1,nb
          if (W(i)<W(j)) then
             value=W(i); W(i)=W(j); W(j)=value
             do k=1,n
                value=Z(k,i); Z(k,i)=Z(k,j); Z(k,j)=value
             end do
          end if
       end do
    end do
    do i=1,nb
       print *,'valeurs propres arpack réordonnées',i, W(i)
    end do

    !Test spectral embedding avec different nbcluster   
    !***********************
    ! Spectral embedding
    !print *,numproc,'Spectral Embedding'

    if ((nbideal==0).and.(n>2)) then
       !** recherche du meilleur decoupage
       allocate(ratiomax(nblimit)); ratiomax(:)=0
       allocate(ratiomin(nblimit)); ratiomin(:)=0
       allocate(ratiomoy(nblimit)); ratiomoy(:)=0
       allocate(ratiorii(nblimit)); ratiorii(:)=0
       allocate(ratiorij(nblimit)); ratiorij(:)=0

       allocate(nbinfo(nblimit)); nbinfo(:)=0

       do nbcluster = 2 ,min(n,nblimit)
          !print *,numproc,'teste avec nbcluster=',nbcluster
          !call flush(6)

          allocate(cluster(n))
          cluster(:)=0.0
          allocate(cluster_center(nbcluster,nbcluster))
          cluster_center(:,:)=0.0
          allocate(cluster_population(nbcluster))
          cluster_population(:)=0.0
          allocate(cluster_energy(nbcluster))
          cluster_energy(:)=0.0

          call sp_spectral_embedding(nbcluster, n, Z, nnz2, AS, IAS, JAS, &
               ratiomax(nbcluster),cluster,cluster_center,cluster_population, &
               cluster_energy,nbinfo(nbcluster),numproc,ratiomoy(nbcluster), &
               ratiorij(nbcluster),ratiorii(nbcluster))

          deallocate(cluster);
          deallocate(cluster_center);
          deallocate(cluster_energy)
          deallocate(cluster_population);
       end do


#if aff
print *, 'ratio de frobenius'
#endif
       !*******************************
       ! Ratio de norme de frobenius
       !print *,numproc,'Ratio max par cluster',ratiomax(2:nblimit)
       !print *,numproc,'Ratio min par cluster',ratiomin(2:nblimit)
       !print *,numproc,'Ratio Rii par cluster',ratiorii(2:nblimit)
       !print *,numproc,'Ratio Rij par cluster',ratiorij(2:nblimit)
       ratio=ratiomax(nblimit)
       dataw%nbclusters=nblimit
       ratio1=0.0
       ratio2=1e+10

       do i=2,nblimit
          !if  ((nbinfo(i)==i).and.(ratio(i)<ratio)) then
          !if  ((nbinfo(i)==i).and.(ratiomoy(i)<ratio)) then
          if ((numproc==0).and.(nbproc>1)) then 
             seuilrij=1e-1
          else
             seuilrij=1e-4
          end if

          if ((ratiorii(i)>=0.95*ratio1).and.(ratiorij(i)-ratio2<=seuilrij)) then  
             !if (ratiomoy(i)-ratio1<=1e-4) then
             !2eme critère
             !(ratiorij(i)/ratiorii(i)<=1e-4)
             dataw%nbclusters=i
             ! ratio=ratiomax(i)
             ratio1=ratiorii(i)
             ratio2=ratiorij(i)
             !ratio1=ratiomoy(i)
          end if
       end do
      ! print *,numproc,'nb de clusters final :',dataw%nbclusters

    elseif ((nbideal==1).and.(n>nbideal)) then
       !** test avec un cluster impose
       allocate(nbinfo(nbideal))
       nbinfo(:) = 0
       allocate(ratiomin(1))
       ratiomin(:) = 0.0
       dataw%nbclusters = nbideal
    else
       !** cas d'un domaine avec moins de points que nbideal ou 1 seul point
       allocate(nbinfo(n))
       nbinfo(:)=0
       allocate(ratiomin(1))
       ratiomin(:)=0.0
       dataw%nbclusters=n
       allocate(ratiomax(n))
       ratiomax(:)=0
       allocate(ratiomoy(n))
       ratiomoy(:)=0
       allocate(ratiomin(n))
       ratiomin(:)=0
       allocate(ratiorii(n))
       ratiorii(:)=0
       allocate(ratiorij(n))
       ratiorij(:)=0
    endif
    ! cas où nbcluster==1
    if (dataw%nbclusters==2) then
       print *, 'difference ratio',ratiorij(2)/ratiorii(2)
       if (ratiomax(2)>=0.6) then 
          dataw%nbclusters=1
       else 
          dataw%nbclusters=2
       end if
    end if
#if aff
    print *,numproc,'cluster final obtenu : ',dataw%nbclusters
#endif

    !** calcul du clustering final
    if (dataw%nbclusters>1) then

       call sp_spectral_embedding(dataw%nbclusters, n, Z, nnz2, AS, IAS, JAS,ratio,cluster,&
            cluster_center,cluster_population,cluster_energy,&
            nbinfo(dataw%nbclusters),numproc,ratiomin(1),ratiorij(1),&
            ratiorii(1))

       do i=1,dataw%nb
          dataw%point(i)%cluster=cluster(i)
       enddo

       deallocate(cluster)
       deallocate(cluster_population)
       deallocate(ratiomax)
       deallocate(cluster_energy)
       deallocate(ratiomin)
       deallocate(ratiomoy)
       deallocate(ratiorii)
       deallocate(ratiorij)
       deallocate(cluster_center)

    else 
       !cluster_population(1)=dataw%nb
#if aff
       print *, numproc, 'ok'
#endif
       do i=1,dataw%nb
          dataw%point(i)%cluster=1
       end do
#if aff
       print *,numproc,'cluster'
#endif
    end if

    !affichage
    !do j=1,dataw%nbclusters
    !   print *,numproc,'centres des clusters',cluster_center(:,j)
    !end do
    !maj du cluster
  

    !deallocations
    deallocate(AS)
    deallocate(IAS)
    deallocate(JAS)
    deallocate(W)
    deallocate(Z)

    return
  end subroutine sp_calculclusters

    subroutine sp_spectral_embedding(nbcluster,n,Z, nnz, AS, IAS, JAS, ratio,cluster,&
       cluster_center,cluster_population,cluster_energy,nbinfo,numproc,&
       ratiomoy,ratiorij,ratiorii)

    !*****************************************
    ! spectral embedding
    !
    ! nbcluster = nbre de cluster
    ! dataw : points
    ! Z : matrice des vecteurs propres
    ! M : nbre de vp trouvéesx
    ! ratio : max des ration de frob sur matrice aff réordonnancée suivant
    ! les clusters
    ! cluster : appartenance des clusters
    ! cluster_center : centre des nbclusters clusters
    ! cluster_population : nbre de points par cluster
    ! cluster_energy : somme des énergies par cluster
    !

    implicit none
    !type(type_data) :: dataw
    double precision,dimension(:,:),pointer:: Z,cluster_center
    integer ::nbcluster,n,nbinfo,numproc
    double precision ::ratio,test,ratiomin,ratiorii,ratiorij, ratiomoy
    double precision,dimension(:),pointer :: cluster_energy,Z3
    integer,dimension(:),pointer ::cluster,cluster_population
    !integer,dimension(:),pointer::ordaffperclus
    double precision, dimension(:,:),pointer :: Frob
    double precision,dimension(:,:),pointer::Z1,Z2
    integer :: it_max,it_num,i,j,k
    integer,dimension(:,:),pointer :: clustercorresp
    integer :: ki,kj,ni,nj,ok,nbmax

    double precision,dimension(:),pointer:: AS
    integer,dimension(:),pointer:: IAS, JAS
    integer :: nnz

    integer :: l
    integer :: num1, num2

    allocate(cluster(n));
    allocate(cluster_center(nbcluster,nbcluster));
    allocate(cluster_population(nbcluster));
    allocate(cluster_energy(nbcluster));
    allocate(Z1(n,nbcluster));  allocate(Z2(nbcluster,n));
    allocate(Z3(n));Z3(:)=0.0

    print *, '************ sp_spectral_embedding *************'
    do i=1,n
       do j=1,nbcluster
          Z1(i,j)=Z(i,j)
          !if (i==1) print *,numproc,'matrice vp',j,W(j)
          Z3(i)=Z3(i)+Z1(i,j)**2
       end do
    end do

    do i=1,n
       test=0.0
       do j=1,nbcluster
          Z2(j,i)=Z1(i,j)/(sqrt(Z3(i)))
          test=test+Z2(j,i)**2
       end do
    end do

    print *, numproc,'methode kmeans'

    it_max=n*n !1000.0

    call kmeans_01 ( nbcluster, n, nbcluster, it_max, it_num, Z2,&
         cluster, cluster_center, cluster_population, cluster_energy, &
         numproc)

    !do i=1,nbcluster
    !   print *, 'Z2(',i,',:) ', Z2(i,:)
    !end do

    !print *,numproc,'fin de kmeans. nb d iterations effectuees : ',it_num

    !print *,numproc,'Nombre points par cluster', cluster_population
    ! print *,'vecteur cluster',cluster(1:5) 

    !*****************************
    ! Mesure de qualité
    !print *,'Indexation'

    nbmax=0
    do i=1,nbcluster
       nbmax=max(nbmax,cluster_population(i))
    end do
    print *, cluster_population
    allocate(clustercorresp(nbcluster,nbmax)); clustercorresp(:,:)=0
    do i=1,n
       j=cluster(i)
       ok=0;k=1
       do while(ok==0)
          if (clustercorresp(j,k)==0) then
             ok=1
          else
             k=k+1
          end if
       end do
       clustercorresp(j,k)=i
    end do


! sparsification début
    allocate(Frob(nbcluster,nbcluster)); Frob(:,:)=0.0
    do i=1, nnz
      num1 = cluster(IAS(i))
      num2 = cluster(JAS(i))
      Frob(num1, num2) = Frob(num1, num2) + AS(i)**2
    end do
! sparsification fin


! sparsification debut
    ratio=0.0; ratiomin=1.D+16;ratiorii=0.0;ratiorij=0.0
    ratiomoy = 0.0
    nbinfo=nbcluster
    do i=1,nbcluster
       if ((cluster_population(i)/=0).and.(Frob(i,i)/=0)) then
          do j=1,nbcluster
             if (i/=j) then
                ratio=ratio+Frob(i,j)/Frob(i,i)
                !ratio=max(ratio,Frob(i,j)/Frob(i,i))
                ratiomoy=ratiomoy+Frob(i,j)/Frob(i,i)
                ratiorij=ratiorij+Frob(i,j)
                ratiorii=ratiorii+Frob(i,i)
                ratiomin=min(ratiomin,Frob(i,j)/Frob(i,i))
             endif
          end do
       else
          nbinfo=nbinfo-1
       end if
       ratiorij=ratiorij*2/(nbcluster*(nbcluster-1))
       ratiomoy=ratiomoy*2/(nbcluster*(nbcluster-1))
       ratiorii=ratiorii!/nbcluster
    end do

    print *, "============= ratio ================", ratiomoy, ratiorij

    deallocate(Frob)
! sparsification fin

#if aff
    print *,numproc,'nbinfo=', nbinfo,' nbcluster=',nbcluster
#endif

    return 
  end subroutine sp_spectral_embedding

  subroutine sp_matvec(A, IA, JA, X, Y, n, nnz)
  implicit none

  double precision, intent(in), dimension(nnz) :: A
  integer, intent(in), dimension(nnz) :: IA, JA
  double precision, intent(in), dimension(n) :: X
  double precision, intent(out), dimension(n) :: Y
  integer, intent(in) :: n, nnz

  integer :: l

  Y(:) = dfloat(0)

  do l = 1, nnz
    Y(IA(l)) = Y(IA(l)) + A(l)*X(JA(l))
  enddo

  return

  end subroutine sp_matvec

  subroutine solve_arpack(A, IA, JA, ndim, nnz, nblimit, W, Z)

  double precision, intent(in), dimension(:) :: A
  integer, intent(in), dimension(:) :: IA, JA
  integer, intent(in) :: ndim, nnz, nblimit

  double precision, intent(out), pointer :: W(:)
  double precision, intent(out), pointer :: Z(:, :)
  
  integer :: maxn, maxnev, maxncv, ldv, i, nbite
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
  integer           iparam(11), ipntr(14)
  logical, dimension(:), allocatable :: select

!     d valeurs propres
!     v vecteurs propres

  Double precision, dimension(:), allocatable :: ax, resid, workd, workev, workl
  Double precision, dimension(:,:), allocatable :: d, v
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
  character         bmat*1, which*2
  integer           ido, n, nx, nev, ncv, lworkl, info, ierr, &
                    j, ishfts, maxitr, mode1, nconv
  Double precision  tol, sigmar, sigmai
  logical           first, rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
  Double precision  zero
  parameter         (zero = 0.0D+0)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  dlapy2, dnrm2
      external          dlapy2, dnrm2, daxpy 
!
!     %--------------------%
!     | Intrinsic function |
!     %--------------------%
!
      intrinsic         abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mnaupd = 1.                    |
!     %-------------------------------------------------%
!
      include 'debug.h'

! lien entre les tailles
      maxn = ndim
      ldv = maxn
      maxnev = nblimit
      maxncv = 2*maxnev + 1

! allocation memoire
      allocate(select(maxn))
      allocate(ax(maxn), resid(maxn), workd(3*maxn), &
               workev(3*maxncv), workl(3*maxncv*maxncv+6*maxncv))
      allocate(d(maxncv, 3), v(ldv, maxncv))

      ndigit = -3
      logfil = 6
      mnaitr = 0
      mnapps = 0
      mnaupd = 1
      mnaup2 = 0
      mneigh = 0
      mneupd = 0
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      n     = ndim

!     %-----------------------------------------------%
!     |                                               |
!     | Specifications for ARPACK usage are set       |
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |
!     |       computed.                               |
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization.                          |
!     |                                               |
!     |    3) This is a standard problem.             |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude.                      |
!     |         (indicated by which = 'LM')           |
!     |       See documentation in DNAUPD for the     |
!     |       other options SM, LR, SR, LI, SI.       |
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 2 <= NCV <= MAXNCV             |
!     |                                               |
!     %-----------------------------------------------%
!
      nev   = maxnev
      ncv   = maxncv
      bmat  = 'I'
      which = 'LR'
!
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DNAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from DNAUPD. (see usage below)                 |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DNAUPD.                                |
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     |
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID).  |
!     |                                                     |
!     | The work array WORKL is used in DNAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+6*ncv 
      tol    = 1.D-6
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting IPARAM(1) = 1).             |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode1 = 1
!
      iparam(1) = ishfts
!
      iparam(3) = maxitr
!
      iparam(7) = mode1
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
      nbite = 1
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dnaupd ( ido, bmat, n, which, nev, tol, resid, ncv, &
                       v, ldv, iparam, ipntr, workd, workl, lworkl, & 
                       info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- Op*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
            !call av (nx, workd(ipntr(1)), workd(ipntr(2)))
            call sp_matvec(A, IA, JA, workd(ipntr(1)), workd(ipntr(2)), &
                           ndim, nnz)

            !if(nbite .eq. 1) then
            !  do i = 0,19
            !    print *,'matvec workd ',i,  workd(ipntr(1)+i), workd(ipntr(2)+i)
            !  end do
            !end if
            nbite = nbite + 1
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DNAUPD again. |
!           %-----------------------------------------%
!
            go to 10
!
         endif
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation of _naupd'
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        |                                           |
!        | The routine DNEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1,)                                   |
!        |                                           |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv, &
              sigmar, sigmai, workev, bmat, n, which, nev, tol, & 
              resid, ncv, v, ldv, iparam, ipntr, workd, workl, &
              lworkl, ierr )
!
!        %------------------------------------------------%
!        | The real parts of the eigenvalues are returned |
!        | in the first column of the two dimensional     |
!        | array D, and the IMAGINARY part are returned   |
!        | in the second column of D.  The corresponding  |
!        | eigenvectors are returned in the first         |
!        | NCONV (= IPARAM(5)) columns of the two         |
!        | dimensional array V if requested.  Otherwise,  |
!        | an orthogonal basis for the invariant subspace |
!        | corresponding to the eigenvalues in D is       |
!        | returned in V.                                 |
!        %------------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DNEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '
!
         else
!
            first = .true.
            nconv =  iparam(5)
            do 20 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (IPARAM(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               if (d(j,2) .eq. zero)  then
!
!                 %--------------------%
!                 | Ritz value is real |
!                 %--------------------%
!
                  !call av(nx, v(1,j), ax)
                  call sp_matvec(A, IA, JA, v(1,j), ax, ndim, nnz)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  d(j,3) = d(j,3) / abs(d(j,1))
!
               else if (first) then
!
!                 %------------------------%
!                 | Ritz value is complex. |
!                 | Residual of one Ritz   |
!                 | value of the conjugate |
!                 | pair is computed.      |
!                 %------------------------%
!
                  !call av(nx, v(1,j), ax)
                  call sp_matvec(A, IA, JA, v(1,j), ax, ndim, nnz)
                  call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                  call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
                  d(j,3) = dnrm2(n, ax, 1)
                  !call av(nx, v(1,j+1), ax)
                  call sp_matvec(A, IA, JA, v(1,j+1), ax, ndim, nnz)
                  call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                  call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
                  d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                  d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                  d(j+1,3) = d(j,3)
                  first = .false.
               else
                  first = .true.
               end if
!
 20         continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call dmout(6, nconv, 3, d, maxncv, -6, &
                 'Ritz values (Real, Imag) and residual residuals')
         end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit &
                        Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
!
         print *, ' '
         print *, ' _NSIMP '
         print *, ' ====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv 
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program dnsimp. |
!     %---------------------------%
!
 9000 continue

      allocate(W(nblimit))
      allocate(Z(n, nblimit))

      W(1:nblimit) = d(1:nblimit, 1)

      do i = 1, nblimit
        Z(:,i) = v(:,i)
      end do

      deallocate(select)
      deallocate(ax, resid, workd, workev, workl)
      deallocate(d, v)

  end subroutine solve_arpack

  end module module_sparse
