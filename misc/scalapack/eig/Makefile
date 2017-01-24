all: pdsyev spec_pdsyev

spec_pdsyev: spectral_pdsyev_call.f
	${MPIF90} -o spec_pdsyev spectral_pdsyev_call.f ${INCLUDE} ${LSCALAPACK} ${LBLAS} -L${MKL} -lmkl_blacs_lp64

pdsyev: sample_distrib_pdsyev_call.f
	${MPIF90} -o pdsyev sample_distrib_pdsyev_call.f ${INCLUDE} ${LSCALAPACK} ${LBLAS} -L${MKL} -lmkl_blacs_lp64
clean:
	rm -f *.o pdsyev spec_pdsyev

# appel de pdsyev avec un nombre de proc égal à NPROW*NPCOL
# mpirun -np # pdsyev
