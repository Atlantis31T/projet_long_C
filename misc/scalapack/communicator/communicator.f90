program hello
  implicit none
   include 'mpif.h'
   integer rank, size, ierror, tag, status(MPI_STATUS_SIZE)
   integer master_com, sous_com
   integer sous_rang, master_rang
   integer somme
   integer world_group, master_group

   integer p, i
   integer, allocatable :: masters(:)

   ! nombre de processus MPI par sous-communicateur
   p = 3
   
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

   ! répartition des processus par sous-communicateur
   ! rank/p : n° du sous-communicateur auquel est assigné le processus rank
   ! mod(rank/p) : rang dans le sous-communicateur
   call MPI_COMM_SPLIT(MPI_COMM_WORLD, rank/p, mod(rank, p), sous_com, ierror)
   ! récupération du rang dans le sous-communicateur (= mod(rank,p) ?)
   call MPI_COMM_RANK(sous_com, sous_rang, ierror)

   !print *, rank, mod(rank,p), sous_rang

   allocate(masters(0:size/p-1))

   masters(0) = 0
   do i = 1, size/p-1
     masters(i) = masters(i-1) + p
   enddo

   print *, masters

   ! groupe associé au communicateur MPI_COMM_WORLD
   call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierror)

   ! création du groupe des masters grâce au tableau masters
   !  size/p : nombre de processus dans le nouveau groupe
   !  masters : tableau des rangs dans le groupe des processus à sélectionner
   !  master_group : identifiant du nouveau groupe
   CALL MPI_GROUP_INCL(world_group, size/p, masters, master_group, ierror)

   ! création du communicateur associé au master_group
   CALL MPI_COMM_CREATE(MPI_COMM_WORLD, master_group, master_com, ierror)
   if(sous_rang == 0) then
     call MPI_COMM_RANK(master_com, master_rang, ierror)

     print *, master_rang
   end if

   ! opération entre processus du même sous-communicateur
   ! (somme des rangs des processus du sous-communicateur)
   somme = 0
   call mpi_allreduce(rank, somme, 1, mpi_integer, MPI_SUM, sous_com, ierror)
   print *, 'com mumps',  rank, somme

   ! opération au sein du groupe des masters
   ! (somme des rangs des processus masters)
   somme = 0
   if(sous_rang == 0) then
     call mpi_allreduce(rank, somme, 1, mpi_integer, MPI_SUM, master_com, ierror)
   end if
   ! résultat 0 pour les processus slaves
   print *, 'com master', rank, somme

   call MPI_FINALIZE(ierror)
   end
