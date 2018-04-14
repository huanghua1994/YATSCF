/*
 * Copyright (c) 2013-2018 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * The GNU Lesser General Public License is included in this distribution
 * in the file COPYING.
 */

#ifndef __CINT_BASISSET_H__
#define __CINT_BASISSET_H__

#include "cint_config.h"

struct BasisSet
{
    // molecular information from xyz file
    int natoms;            // number of atoms in molecule
    int *eid;              // atomic numbers
    double *xn;            // x coords
    double *yn;            // y coords
    double *zn;            // z coords
    double *charge;        // double precision version of atomic numbers
    int nelectrons;        // sum of atomic numbers in molecule (not really num electrons)
    double **guess;        // initial guesses for each element in basis set (should not be in this section)
    int Q;                 // net charge read from xyz file (not related to nelectrons)
    double ene_nuc;        // nuclear energy (computed)

    // basis set information from gbs file
    int bs_nelements;      // max number of elements supported in basis set
    int bs_natoms;         // number of elements in basis set
    int basistype;         // Cartesian or spherical
    int *bs_eid;           // atomic numbers of elements in basis set
    int *bs_eptr;          // map atomic number to entry in basis set (array of len bs_nelements)
    int *bs_atom_start;    // start of element data in arrays of length nshells (array of length natoms+1)
    int bs_nshells;        // number of shells in the basis set (not the molecule)
    int bs_totnexp;        // total number of primitive functions in basis set
    int *bs_nexp;          // number of primitive functions for shell
    double **bs_exp;       // bs_exp[i] = orbital exponents for shell i
    double **bs_cc;        // bs_cc[i]  = contraction coefs for shell i
    double **bs_norm;      // bs_norm[i] = normalization constants for shell i
    int *bs_momentum;      // bs_momentum[i] = angular momentum for shell i
    
    // shell information for each shell in the given molecule
    uint32_t nshells;      // number of shells in given molecule
    uint32_t nfunctions;   // number of basis functions for molecule
    uint32_t *f_start_id;  // offset for first basis function for each shell
    uint32_t *f_end_id;    // offset for last basis function for each shell
    uint32_t *s_start_id;  // start of shell info for each atom
    uint32_t *nexp;        // number of primitives for each shell
    double **exp;          // exponents for each shell in molecule
    double *minexp;        // ?
    double **cc;           // contraction coefficients for each shell in molecule
    double **norm;         // ?
    uint32_t *momentum;    // angular momentum for each shell in molecule
    double *xyz0;          // centers for each shell in molecule, stored as linear array

    uint32_t maxdim;       // max number of functions among all shells in molecule
    uint32_t max_momentum;
    uint32_t max_nexp;
    uint32_t max_nexp_id;
    
    char str_buf[512];
};

typedef struct BasisSet *BasisSet_t;

#ifdef __cplusplus
extern "C" {
#endif

CIntStatus_t CInt_createBasisSet(BasisSet_t *basis);

CIntStatus_t CInt_destroyBasisSet(BasisSet_t basis);

CIntStatus_t CInt_loadBasisSet(BasisSet_t basis, char *bsfile, char *xyzfile );

int CInt_getNumAtoms(BasisSet_t basis);

int CInt_getNumShells(BasisSet_t basis);

int CInt_getNumFuncs(BasisSet_t basis);

int CInt_getNumOccOrb(BasisSet_t basis);

int CInt_getFuncStartInd(BasisSet_t basis, int shellid);

int CInt_getFuncEndInd(BasisSet_t basis, int shellid);

int CInt_getShellDim(BasisSet_t basis, int shellid);

int CInt_getMaxShellDim(BasisSet_t basis);

int CInt_getAtomStartInd(BasisSet_t basis, int atomid);

int CInt_getAtomEndInd(BasisSet_t basis, int atomid);

int CInt_getTotalCharge(BasisSet_t basis);

int CInt_getNneutral(BasisSet_t basis);

int CInt_getMaxMomentum(BasisSet_t basis);

int CInt_getMaxPrimid(BasisSet_t basis);

int CInt_getMaxnumExp(BasisSet_t basis);

double CInt_getNucEnergy(BasisSet_t basis);

void CInt_getInitialGuess(BasisSet_t basis, int atomid, double **guess, int *spos, int *epos);

void CInt_getShellxyz(BasisSet_t basis, int shellid, double *x, double *y, double *z);

// The following 4 functions are called by CInt_loadBasisSet, be careful and make sure you 
// understand what you are doing if you want to call any of them from external program
CIntStatus_t CInt_import_molecule(char *file, BasisSet_t basis);
CIntStatus_t CInt_import_basis(char *file, BasisSet_t basis);
CIntStatus_t CInt_parse_molecule(BasisSet_t basis);
CIntStatus_t CInt_import_guess(char *file, BasisSet_t basis);

#ifdef __cplusplus
}
#endif

#endif
