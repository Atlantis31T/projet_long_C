all: communicator blax

communicator: communicator.f90
	${MPIF90} -o communicator communicator.f90 

blax: blax.f90 test_pdgesv.f
	${MPIF90} -o blax blax.f90 test_pdgesv.f ${INCLUDE} ${LSCALAPACK} ${LBLAS} -L${MKL} -lmkl_blacs_lp64

clean:
	rm -f *.o blax communicator

# appel de communicator ou blax avec un nombre de proc multiple de p
# mpirun -np # blax
