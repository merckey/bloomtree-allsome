#ifndef bitvector_subsets_H		// (prevent multiple inclusion)
#define bitvector_subsets_H

#include "utilities.h"

// establish ownership of global variables

#ifdef bitvector_subsets_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef bitvector_subsets_owner
int bitvector_subsets_dbgBitvectors = false;
#else
global int bitvector_subsets_dbgBitvectors;
#endif


//----------
//
// types--
//
//----------

// bitvector subset specifiers

typedef struct subset
	{
	u32		firstBit;			// *bit* count to first bit in the subset
								// .. (origin zero; must be a multiple of 8)
	u32		numBits;			// number of *bits* in the subset
								// .. (must be a multiple of 8)
	} subset;

typedef struct subsetspecs
	{
	u32		len;				// number of entries in subset[]
	u32		totalBits;			// number of *bits* in the subsets, combined
								// .. (must be a multiple of 8)
	subset	subset[1];			// variable-length array of subsets
	} subsetspecs;

// bit counts

typedef struct bitcounts
	{
	u32		xorCount;
	u32		andCount;
	u32		orCount;
	} bitcounts;

//----------
//
// prototypes for functions in this module--
//
//----------

subsetspecs* new_subset_array  (u32 numSubsets);
void         parse_subset      (char* s, subset* ss);
void         read_bitvector    (char* pathPrefix, char* filename, char* suffix,
                                subsetspecs* subsetSpecs, u8* bv);
void         compute_distances (subsetspecs* subsetSpecs, u8* bv1, u8* bv2,
                                bitcounts* countsVector);
void         compute_densities (subsetspecs* subsetSpecs, u8* bv,
                                bitcounts* countsVector);

#undef global
#endif // bitvector_subsets_H
