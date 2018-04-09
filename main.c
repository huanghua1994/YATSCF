#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "CInt.h"

static void print_usage(char *exe_name)
{
	printf("Usage: %s <basis> <xyz> <niter>\n", exe_name);
}

int main(int argc, char **argv)
{
	if (argc < 4)
	{
		print_usage(argv[0]);
		return 255;
	}
	
	BasisSet_t basis;
	int natoms, nshells, nfunctions, niters;
	
	// Create and load basis set from input
	CInt_createBasisSet(&basis);
	CInt_loadBasisSet(basis, argv[1], argv[2]);
	
	natoms     = CInt_getNumAtoms(basis);
	nshells    = CInt_getNumShells(basis);
	nfunctions = CInt_getNumFuncs(basis);
	niters     = atoi(argv[3]);
	
	printf("Job information:\n");
	printf("    molecule  file    = %s\n", argv[1]);
	printf("    basis set file    = %s\n", argv[2]);
	printf("    charge            = %d\n", CInt_getTotalCharge(basis));
	printf("    # atoms           = %d\n", natoms);
	printf("    # shells          = %d\n", nshells);
	printf("    # basis functions = %d\n", nfunctions);
}