// bitvector_subsets.c-- support of contiguous subsets of bitvectors.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include "utilities.h"

#define  bitvector_subsets_owner// (make this the owner of its globals)
#include "bitvector_subsets.h"	// interface to this module

//----------
//
// new_subset_array--
//	Allocate an array to hold a specified number of subsets.
//
//----------
//
// Arguments:
//	u32		numSubsets:		The number of subsets needed.
//
// Returns:
//	A pointer to a subsetspecs record, allocated from the heap;  the caller is
//	responsible for deallocating this memory;  failure to allocate results in
//	program termination.
//
//----------

subsetspecs* new_subset_array
   (u32				_numSubsets)
	{
	subsetspecs*	specs = NULL;
	u32				numSubsets;
	u32				bytesNeeded;

	if (_numSubsets < 1) numSubsets = 1;
	                else numSubsets = _numSubsets;

	bytesNeeded = sizeof(subsetspecs) + (numSubsets-1) * sizeof(subset);

	specs = malloc (bytesNeeded);
	if (specs == NULL) goto cant_allocate;
	specs->len       = _numSubsets;
	specs->totalBits = 0;
	return specs;

	// failures

cant_allocate:
	fprintf (stderr, "failed to allocate %u bytes for %u subsets\n",
					 bytesNeeded, _numSubsets);
	exit (EXIT_FAILURE);
	return NULL; // (never reaches here)
	}

//----------
//
// parse_subset--
//	Convert a string into a bitvector subset.
//
//----------
//
// Arguments:
//	char*	s:	The string to parse.  This should be in one of the forms
//				<lo>..<hi> or <numbits>@<lo>, where <lo>, <hi> and <numbits>
//				are unsigned integers.
//	subset*	ss:	Place to store the result
//
// Returns:
//	nothing;  parsing failures result in program termination.
//
//----------

void parse_subset
   (char*	s,
	subset*	ss)
	{
	char	loString[20], numBitsString[20];
	char*	dots, *atSign;
	u32		loLen, numBitsLen, lo, hi, numBits;

	dots = strstr (s,"..");
	if (dots != NULL)
		{
		if (dots == s) goto no_lo;
		loLen = dots - s;
		if (loLen == strlen(s)-2) goto no_hi;

		if (loLen > sizeof(loString)-1) goto lo_too_long;
		strncpy (loString,s,loLen);
		loString[loLen] = 0;

		lo = string_to_u32(loString);
		hi = string_to_u32(dots+2);
		numBits = hi - lo;

		if (hi == lo) goto empty_subset;
		if (hi <  lo) goto null_subset;
		if (lo % 8 != 0) goto lo_not_mod_8;
		if (hi % 8 != 0) goto hi_not_mod_8;
		goto success;
		}

	atSign = strchr (s,'@');
	if (dots != NULL)
		{
		if (atSign == s) goto no_numbits;
		numBitsLen = atSign - s;
		if (numBitsLen == strlen(s)-1) goto no_lo;

		if (numBitsLen > sizeof(numBitsString)-1) goto numbits_too_long;
		strncpy (numBitsString,s,numBitsLen);
		numBitsString[numBitsLen] = 0;

		numBits = string_to_u32(numBitsString);
		lo = string_to_u32(atSign+1);
		hi = lo + numBits;

		if (numBits == 0) goto numbits_zero;
		if (numBits % 8 != 0) goto numbits_not_mod_8;
		if (lo % 8 != 0) goto lo_not_mod_8;
		if (hi < lo) goto overflow;
		goto success;
		}

	goto cant_parse;

success:
	ss->firstBit = lo;
	ss->numBits  = numBits;
	return;

	// failures

cant_parse:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

no_lo:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (no lo value)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

no_hi:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (no hi value)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

no_numbits:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (no numbits value)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

lo_too_long:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (lo value too long)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

numbits_too_long:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (numbits value too long)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

empty_subset:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (lo=hi)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

null_subset:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (hi<lo)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

numbits_zero:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (numbits=0)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

lo_not_mod_8:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (lo isn't a multiple of 8)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

hi_not_mod_8:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (hi isn't a multiple of 8)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

numbits_not_mod_8:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (numbits isn't a multiple of 8)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)

overflow:
	fprintf (stderr, "failed to parse \"%s\" as a subset specifier (lo+numbits>=2^32)\n", s);
	exit (EXIT_FAILURE);
	return; // (never reaches here)
	}

//----------
//
// read_bitvector--
//	Read a bitvector from a file.  The bitvector consists of the concatenation
//	of a series of subsets of the bitvector in the file.
//
// Bit vectors in files are stored with 4 leading bytes (little endian)
// indicating the number of bits in the file, and the next 4 bytes containing
// some unknown information.  E.g.
//	00000000  00 94 35 77 00 00 00 00  01 00 00 00 00 10 00 00
//	00000010  00 00 00 00 00 00 00 00  00 00 00 04 00 10 00 00
//	 ...
//
//----------
//
// Arguments:
//	char*			pathPrefix:		a prefix to prepend to each filename;  if
//									.. this doesn't end with a "/", one is
//									.. added;  this may be NULL or empty, in
//									.. which case no path is added
//	char*			filename:		name of the binary file which contains the
//									.. bitvector
//	char*			suffix:			a suffix to append to each filename;  this
//									.. may be NULL or empty, in which case no
//									.. suffix is added
//	subsetspecs*	subsetSpecs:	specifiers for the subsets to read
//	u8*				bv:				place to write the bitvector to;  this must
//									.. have at least as many bytes as are
//									.. required to satisfy the subsetSpecs
//
// Returns:
//	(nothing)
//
//----------

#define bitvectorFileHeaderBytes 8

void read_bitvector
   (char*			pathPrefix,
	char*			filename,
	char*			suffix,
	subsetspecs*	subsetSpecs,
	u8*				bv)
	{
	char			_combinedName[5001];
	char*			combinedName;
	u32				pathLen, suffLen, combinedLen;
	int				addPathSlash;
	char*			scan;
	FILE*			f;
	subset*			ss;
	u32				subsetIx;
	u32				subsetByteOffset, bytesToRead, bytesRead;
	int				err;

	// build the filename

	if (pathPrefix != NULL)
		{ if (*pathPrefix == 0) pathPrefix = NULL; }

	if (suffix != NULL)
		{ if (*suffix == 0) suffix = NULL; }

	if ((pathPrefix == NULL) && (suffix == NULL))
		combinedName = filename;
	else
		{
		pathLen = 0;
		addPathSlash = false;
		if (pathPrefix != NULL)
			{
			pathLen = strlen(pathPrefix);
			if (pathPrefix[pathLen-1] != '/')
				{ pathLen++;  addPathSlash = true; }
			}

		suffLen = 0;
		if (suffix != NULL)
			suffLen = strlen(suffix);

		combinedLen = pathLen + strlen(filename) + suffLen;
		if (combinedLen + 1 > sizeof(_combinedName)) goto name_too_big;

		scan = _combinedName;
		if (pathPrefix != NULL)
			{
			strcpy (scan, pathPrefix);
			if (addPathSlash) scan[pathLen-1] = '/';
			scan += pathLen;
			}

		strcpy (scan, filename);
		scan += strlen(filename);

		if (suffix != NULL)
			strcpy (scan, suffix);

		combinedName = _combinedName;
		}

	if (bitvector_subsets_dbgBitvectors)
		fprintf (stderr, "reading from \"%s\"\n", combinedName);

	// read the bitvector subsets from the file

	f = fopen (combinedName, "rb");
	if (f == NULL) goto cant_open_file;

	subsetByteOffset = 0;
	for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
		{
		ss = &subsetSpecs->subset[subsetIx];
		bytesToRead = ss->numBits / 8;
		err = fseek (f, bitvectorFileHeaderBytes + ss->firstBit/8, SEEK_SET);
		if (err != 0) goto seek_failed;
		bytesRead = fread (bv+subsetByteOffset, 1, bytesToRead, f);
		if (bytesRead != bytesToRead) goto read_failed;

		if (bitvector_subsets_dbgBitvectors)
			{
			u8* scan;
			u32 scanIx;

			fprintf (stderr, "subset %u (%u bits starting at bit %u)\n",
			                 1+subsetIx, ss->numBits, ss->firstBit);
			scan = bv + subsetByteOffset;
			for (scanIx=0 ; scanIx<bytesToRead ; scanIx++)
				{
				if      (scanIx      == 0) ;
				else if (scanIx % 16 == 0) fprintf (stderr, "\n");
				else if (scanIx %  4 == 0) fprintf (stderr, "  ");
				else                       fprintf (stderr, " ");
				fprintf (stderr, "%02X", scan[scanIx]);
				}
			fprintf (stderr, "\n");
			}

		subsetByteOffset += bytesToRead;
		}

	fclose (f);

	return;

	// failures

name_too_big:
	fprintf (stderr, "combined name for \"%s\" and \"%s\" exceeds buffer size\n",
	                 pathPrefix, filename);
	exit(EXIT_FAILURE);
	return; // (never reaches here)

cant_open_file:
	fprintf (stderr, "failed to open \"%s\"\n", combinedName);
	exit(EXIT_FAILURE);
	return; // (never reaches here)

seek_failed:
	fprintf (stderr, "failed to seek to %u in \"%s\"\n"
	                 "fseek() returned %d\n",
	                 bitvectorFileHeaderBytes + subsetByteOffset, combinedName,
	                 err);
	exit(EXIT_FAILURE);
	return; // (never reaches here)

read_failed:
	fprintf (stderr, "failed to read %u from \"%s\"\n"
	                 "fread() returned %u\n",
	                 bytesToRead, combinedName, bytesRead);
	exit(EXIT_FAILURE);
	return; // (never reaches here)
	}

//----------
//
// fill_pop_count--
//	Fill the pop-count lookup table, if it hasn't already been filled.
//
//----------

static int popCountFilled = false;
static u32 popCount8[0x100];

static void fill_pop_count(void);
static void fill_pop_count(void)
	{
	u32 ix;

	if (popCountFilled) return;

	popCount8[0] = 0;
	for (ix=1 ; ix<0x100 ; ix++)
		popCount8[ix] = popCount8[ix>>1] + (ix&1);

	popCountFilled = true;
	}

//----------
//
// compute_distances--
//	Compute the distances between two bitvectors.  A separate distance record
//	is computed for each subset.
//
//----------
//
// Arguments:
//	subsetspecs*	subsetSpecs:	specifiers for the subsets to compare
//	u8*				bv1, bv2:		bitvectors to compare
//	bitcounts*		countsVector:	place to return the distances for each
//									subset; this must have at least as many
//									entries as subsetSpecs->len
//
// Returns:
//	(nothing)
//
//----------

void compute_distances
   (subsetspecs*	subsetSpecs,
	u8*				bv1,
	u8*				bv2,
	bitcounts*		countsVector)
	{
	u8*				scan1, *scan2;
	u32				xorCount, andCount, orCount;
	u32				subsetIx, subsetBytes;
	subset*			ss;

	if (!popCountFilled) fill_pop_count();

	for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
		{
		ss = &subsetSpecs->subset[subsetIx];
		scan1 = bv1 + (ss->firstBit / 8);
		scan2 = bv2 + (ss->firstBit / 8);
		subsetBytes = ss->numBits / 8;

		xorCount = andCount = orCount = 0;
		while (subsetBytes-- > 0)
			{
			xorCount += popCount8[*scan1 ^ *scan2];
			andCount += popCount8[*scan1 & *scan2];
			orCount  += popCount8[*scan1 | *scan2];
			scan1++;  scan2++;
			}

		countsVector[subsetIx].xorCount = xorCount;
		countsVector[subsetIx].andCount = andCount;
		countsVector[subsetIx].orCount  = orCount;
		}
	}

//----------
//
// compute_densities--
//	Compute the densities of a bitvector.  A separate density record
//	is computed for each subset.
//
//----------
//
// Arguments:
//	subsetspecs*	subsetSpecs:	specifiers for the subsets to compare
//	u8*				bv:				bitvector to analyze
//	bitcounts*		countsVector:	place to return the densities for each
//									subset; this must have at least as many
//									entries as subsetSpecs->len
//
// Returns:
//	(nothing)
//
//----------

void compute_densities
   (subsetspecs*	subsetSpecs,
	u8*				bv,
	bitcounts*		countsVector)
	{
	u8*				scan;
	u32				count;
	u32				subsetIx, subsetBytes;
	subset*			ss;

	if (!popCountFilled) fill_pop_count();

	for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
		{
		ss = &subsetSpecs->subset[subsetIx];
		scan = bv + (ss->firstBit / 8);
		subsetBytes = ss->numBits / 8;

		count = 0;
		while (subsetBytes-- > 0)
			{ count += popCount8[*scan];  scan++; }

		countsVector[subsetIx].xorCount = count;
		countsVector[subsetIx].andCount = count;
		countsVector[subsetIx].orCount  = count;
		}
	}
