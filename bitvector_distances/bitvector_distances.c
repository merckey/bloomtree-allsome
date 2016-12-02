// bitvector_distances.c-- compute the (Hamming) distance between bitvectors.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "utilities.h"
#include "bitvector_subsets.h"

char* programName = "bitvector_distances";

//----------
//
// global data and types--
//
//----------

// string list (array)

typedef struct stringlist
	{
	u32		len;				// number of entries in subset[]
	char*	s[1];				// variable-length array of strings;  this is
								// .. carved from the same heap block
	} stringlist;

// command line options

char*			fileListName = NULL;
char*			filePath     = NULL;
char*			fileSuffix   = NULL;
subsetspecs*	subsetSpecs  = NULL;

int             reportDensity       = false;
int             writeHeader         = true;
int             showLoadProgress    = false;
int             showComputeProgress = false;

int				dbgDumpSubsets = false;
int				dbgStringList  = false;
int				dbgFilenames   = false;
int				dbgBitvectors  = false;
int				dbgIntDistance = false;
int				dbgShowAndOr   = false;

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);

// private functions

void         parse_options     (int argc, char** argv);
stringlist*  read_filenames    (char* listFilename);
void         print_distances   (FILE* f, subsetspecs* subsetSpecs,
                                bitcounts* countsVector);
void         print_densities   (FILE* f, subsetspecs* subsetSpecs,
                                bitcounts* countsVector);

//----------
//
// main program--
//
//----------

int main
   (int			argc,
	char**		argv)
	{
	stringlist*	bvFilenames  = NULL;
	u8**		bitVectors   = NULL;
	bitcounts*	countsVector = NULL;
	u32			subsetBytesNeeded;
	u32			subsetIx;
	subset*		ss;
	u32			bvIx, bv1Ix, bv2Ix;

	parse_options (argc, argv);

	bvFilenames = read_filenames (fileListName);

	// allocate bits for bit vector subsets

	bitVectors = malloc (bvFilenames->len * sizeof(u8*));
	if (bitVectors == NULL) goto cant_allocate_pointers;
	for (bvIx=0 ; bvIx<bvFilenames->len ; bvIx++)
		bitVectors[bvIx] = NULL;

	subsetBytesNeeded = subsetSpecs->totalBits / 8;
	if (dbgDumpSubsets)
		fprintf (stderr, "subsetBytesNeeded = %u\n", subsetBytesNeeded);
	for (bvIx=0 ; bvIx<bvFilenames->len ; bvIx++)
		{
		bitVectors[bvIx] = malloc (subsetBytesNeeded);
		if (bitVectors[bvIx] == NULL) goto cant_allocate_bitvector;
		if (dbgBitvectors)
			fprintf (stderr, "bitvector[%u] = %p\n", bvIx, bitVectors[bvIx]);
		}

	// read bit vector subsets

	for (bvIx=0 ; bvIx<bvFilenames->len ; bvIx++)
		{
		if (showLoadProgress)
			fprintf (stderr, "loading bitvector %s (%u of %u)\n",
			                 bvFilenames->s[bvIx], 1+bvIx, bvFilenames->len);
		read_bitvector (filePath, bvFilenames->s[bvIx], fileSuffix,
		                subsetSpecs, bitVectors[bvIx]);
		}

	// compute and print densities

	countsVector = malloc (subsetSpecs->len * sizeof(bitcounts));
	if (countsVector == NULL) goto cant_allocate_counts_vector;

	if (reportDensity)
		{
		if (writeHeader)
			{
			printf ("#vector");
			for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
				{
				ss = &subsetSpecs->subset[subsetIx];
				printf ("\t%u@%u", ss->numBits, ss->firstBit);
				}
			printf ("\n");
			}

		for (bv1Ix=0 ; bv1Ix<bvFilenames->len ; bv1Ix++)
			{
			u8*   bv1   = bitVectors[bv1Ix];
			char* name1 = bvFilenames->s[bv1Ix];

			if (showComputeProgress)
				fprintf (stderr, "computing all bitvectors vs. %s (%u of %u)\n",
								 name1, 1+bv1Ix, bvFilenames->len);

			compute_densities (subsetSpecs, bv1, countsVector);
			printf ("%s", name1);

			print_densities (stdout, subsetSpecs, countsVector);
			printf ("\n");
			}
		}

	// -OR- compute and print distances

	else
		{
		if (writeHeader)
			{
			printf ("#vector1\tvector2");
			for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
				{
				ss = &subsetSpecs->subset[subsetIx];
				if (dbgShowAndOr)
					printf ("\t%u@%u\t(and)\t(or)", ss->numBits, ss->firstBit);
				else
					printf ("\t%u@%u", ss->numBits, ss->firstBit);
				}
			printf ("\n");
			}

		for (bv1Ix=0 ; bv1Ix<bvFilenames->len-1 ; bv1Ix++)
			{
			u8*   bv1   = bitVectors[bv1Ix];
			char* name1 = bvFilenames->s[bv1Ix];

			if (showComputeProgress)
				fprintf (stderr, "computing all bitvectors vs. %s (%u of %u)\n",
								 name1, 1+bv1Ix, bvFilenames->len);

			for (bv2Ix=bv1Ix+1 ; bv2Ix<bvFilenames->len ; bv2Ix++)
				{
				u8*   bv2   = bitVectors[bv2Ix];
				char* name2 = bvFilenames->s[bv2Ix];

				compute_distances (subsetSpecs, bv1, bv2, countsVector);
				printf ("%s\t%s", name1, name2);
				print_distances (stdout, subsetSpecs, countsVector);
				printf ("\n");
				}
			}
		}

	// relinquish allocated memory

	if (fileListName != NULL) { free (fileListName);  fileListName = NULL; }
	if (filePath     != NULL) { free (filePath);      filePath     = NULL; }
	if (fileSuffix   != NULL) { free (fileSuffix);    fileSuffix   = NULL; }
	if (subsetSpecs  != NULL) { free (subsetSpecs);   subsetSpecs  = NULL; }

	if (bitVectors != NULL)
		{
		for (bvIx=0 ; bvIx<bvFilenames->len ; bvIx++)
			{ if (bitVectors[bvIx] != NULL) free (bitVectors[bvIx]); }
		free (bitVectors);
		bitVectors = NULL;
		}

	if (bvFilenames  == NULL) { free (bvFilenames);  bvFilenames  = NULL; }
	if (countsVector == NULL) { free (countsVector); countsVector = NULL; }

	// proclaim success

	return EXIT_SUCCESS;

	// failures

cant_allocate_pointers:
	fprintf (stderr, "failed to allocate pointers for %u bit vectors\n",
					 bvFilenames->len);
	exit (EXIT_FAILURE);
	return EXIT_FAILURE; // (never reaches here)

cant_allocate_bitvector:
	fprintf (stderr, "failed to allocate %u bytes for bitvector %u\n",
					 subsetBytesNeeded, bvIx);
	exit (EXIT_FAILURE);
	return EXIT_FAILURE; // (never reaches here)

cant_allocate_counts_vector:
	fprintf (stderr, "failed to allocate %u bytes for counts vector\n",
					 (u32) (subsetSpecs->len * sizeof(bitcounts)));
	exit (EXIT_FAILURE);
	return EXIT_FAILURE; // (never reaches here)
	}

//----------
//
// option parsing--
//
//----------

static void  chastise  (const char* format, ...);
static void  usage     (char* message);


static void chastise (const char* format, ...)
	{
	va_list	args;

	va_start (args, format);
	if (format != NULL)
		vfprintf (stderr, format, args);
	va_end (args);

	usage (NULL);
	}

static void usage (char* message)
	{
	if (message != NULL) fprintf (stderr, "%s\n", message);

	fprintf (stderr, "%s-- all-vs-all bitvector (Hamming) distance\n",
	                  programName);
	fprintf (stderr, "usage: | %s <filenames_file> [options]\n", programName);
	fprintf (stderr, "\n");
	//                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	fprintf (stderr, "  <filenames_file>    (required) file containing the filenames of the bit\n");
	fprintf (stderr, "                      vectors, one per line\n");
	fprintf (stderr, "  <lo>..<hi>          (cumulative) bit subset to compute distance over;  <lo>\n");
	fprintf (stderr, "                      and <hi> are counted in bits, are origin zero half open,\n");
	fprintf (stderr, "                      and must be multiples of 8\n");
	fprintf (stderr, "                      (at least one subset is required)\n");
	fprintf (stderr, "  <numbits>@<lo>      same as <lo>..<lo+numbits>\n");
	fprintf (stderr, "  --path=<string>     path to be prepended to every bitvector filename\n");
	fprintf (stderr, "                      this does NOT apply to <filenames_file>, only to the\n");
	fprintf (stderr, "                      files named inside it\n");
	fprintf (stderr, "  --suffix=<string>   suffix to be appended to every bitvector filename\n");
	fprintf (stderr, "                      this does NOT apply to <filenames_file>, only to the\n");
	fprintf (stderr, "                      files named inside it\n");
	fprintf (stderr, "  --density           report density instead of pairwise distances\n");
	fprintf (stderr, "  --noheader          don't include a header in the output\n");
	fprintf (stderr, "  --progress:load     report files as we load them\n");
	fprintf (stderr, "  --progress:compute  report bit vector names as we compute them\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Output is a tab-delimited table, like this:\n");
	fprintf (stderr, "  #vector1   vector2   500000@0 100000@500000 ...\n");
	fprintf (stderr, "  SRR1186053 SRR805782 0.014428 0.014090      ...\n");
	fprintf (stderr, "  SRR1186053 SRR837459 0.024678 0.024760      ...\n");
	fprintf (stderr, "  SRR1186053 SRR837458 0.024362 0.024400      ...\n");
	fprintf (stderr, "   ...\n");
	fprintf (stderr, "Column headers (other than vector1 and vector2) are of the form <bits>@<start>,\n");
	fprintf (stderr, "where <bits> is the size of the subset used, and <start> is the first bit in\n");
	fprintf (stderr, "subset (origin zero).\n");

	exit (EXIT_FAILURE);
	}


void parse_options (int _argc, char** _argv)
	{
	int		argc;
	char**	argv;
	int		argIx;
	char*	arg, *argVal;
	int		numSubsets, subsetNum;

	// skip program name

	//programName = _argv[0];
	argv = _argv+1;  argc = _argc - 1;

	if (argc <= 0)
		chastise (NULL);

	//////////
	// scan arguments
	//
	// first pass, counting subsets
	//////////

	numSubsets = 0;
	for (argIx=0 ; argIx<argc ; argIx++)
		{
		arg = argv[argIx];
		if ((strstr(arg,"..") != NULL) || (strchr(arg,'@') != NULL))
			numSubsets++;
		}

	// make sure we got at least one subset

	if (numSubsets == 0)
		chastise ("you have to give me at least one subset to compute distances over\n");

	// allocate subset array

	subsetSpecs = new_subset_array(numSubsets);

	//////////
	// scan arguments
	//
	// second pass, normal parsing
	//////////

	subsetNum = 0;
	for (argIx=0 ; argIx<argc ; argIx++)
		{
		arg = argv[argIx];
		if (arg[0] == 0) continue;
		argVal = strchr(arg,'=');
		if (argVal != NULL) argVal++;

		// --path=<string>

		if ((strcmp_prefix (arg, "--path=") == 0))
			{ filePath = copy_string (argVal);  continue; }

		// --suffix=<string>

		if ((strcmp_prefix (arg, "--suffix=") == 0))
			{ fileSuffix = copy_string (argVal);  continue; }

		// --density

		if (strcmp (arg, "--density") == 0)
			{ reportDensity = true;  continue; }

		// --noheader

		if (strcmp (arg, "--noheader") == 0)
			{ writeHeader = false;  continue; }

		// --progress:load

		if ((strcmp (arg, "--progress:load") == 0)
		 || (strcmp (arg, "--progress=load") == 0))
			{ showLoadProgress = true;  continue; }

		// --progress:compute

		if ((strcmp (arg, "--progress:compute") == 0)
		 || (strcmp (arg, "--progress=compute") == 0))
			{ showComputeProgress = true;  continue; }

		// --debug=<whatever>

		if (strcmp (arg, "--debug=subsets") == 0)
			{ dbgDumpSubsets = true;  continue; }

		if (strcmp (arg, "--debug=stringlist") == 0)
			{ dbgStringList = true;  continue; }

		if (strcmp (arg, "--debug=filenames") == 0)
			{ dbgFilenames = true;  continue; }

		if (strcmp (arg, "--debug=bitvectors") == 0)
			{
			dbgBitvectors = true;
			bitvector_subsets_dbgBitvectors = true;
			continue;
			}

		if (strcmp (arg, "--debug=distance:int") == 0)
			{ dbgIntDistance = true;  continue; }

		if ((strcmp (arg, "--debug=show:andor") == 0)
		 || (strcmp (arg, "--debug=andor") == 0))
			{ dbgShowAndOr = true;  continue; }

		// unrecognized --option

		if ((strcmp_prefix (arg, "--") == 0))
			{
			chastise ("unrecognized option: \"%s\"\n", arg);
			continue; // (never reaches here)
			}

		// <lo>..<hi> or <numbits>@<lo>

		if ((strstr(arg,"..") != NULL) || (strchr(arg,'@') != NULL))
			{
			subset* ss = &subsetSpecs->subset[subsetNum];
			parse_subset (arg, ss);
			subsetSpecs->totalBits += ss->numBits;
			subsetNum += 1;
			continue;
			}

		// <filenames_file>

		if (fileListName == NULL)
			{
			fileListName = copy_string (arg);
			continue;
			}

		// extra filename?

		chastise ("unrecognized option: \"%s\"\n", arg);
		// (never reaches here)
		}

	//////////
	// sanity checks
	//////////

	// make sure we got at least one subset

	if (fileListName == NULL)
		chastise ("you have to give me the name of a file that lists the bitvector filenames\n");

	if (dbgDumpSubsets)
		{
		u32 subsetIx;
		subset* ss;

		for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
			{
			ss = &subsetSpecs->subset[subsetIx];
			fprintf (stderr, "%u bits starting at bit %u\n", ss->numBits, ss->firstBit);
			}
		}

	return;
	}

//----------
//
// read_filenames--
//	Read a list of filenames and allocate an array to hold them.
//
//----------
//
// Arguments:
//	char*	listFilename:	name of a file which contains a list of filenames.
//
// Returns:
//	A pointer to a stringlist record, allocated from the heap;  the caller is
//	responsible for deallocating this memory;  failure to allocate results in
//	program termination.
//
//----------

stringlist* read_filenames
   (char*		listFilename)
	{
	char		line[5001];
	stringlist*	sl = NULL;
	FILE*		f;
	u32			lineNum;
	int			missingEol;
	char*		waffle;
	u32			len, nameLen;
	u32			numNames, nameBytes, nameIx;
	u32			bytesNeeded;
	char*		nameCharScan;

	// first pass through the file, to count items and bytes

	f = fopen (listFilename, "rt");
	if (f == NULL) goto cant_open_file;

	numNames = nameBytes = 0;

	lineNum = 0;
	missingEol = false;
	while (true)
		{
		// get the next line, if we need one;  we also check for lines getting
		// split by fgets (the final line in the file might not have a newline,
		// but no internal lines can be that way)

		if (fgets (line, sizeof(line), f) == NULL) break;
		lineNum++;

		if (missingEol)
			goto split_line;

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (len == 0) continue;
		missingEol = (line[len-1] != '\n');
		if (line[len-1] == '\n') line[--len] = 0;
		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;
		trim_string (line);
		if (line[0] == 0) continue;

		// process the line

		nameLen = strlen(line);
		numNames++;
		nameBytes += nameLen + 1;
		}

	if (dbgStringList)
		fprintf (stderr, "%u names  %u bytes\n", numNames, nameBytes);

	// allocate

	bytesNeeded = sizeof(stringlist) \
	            + (numNames-1) * sizeof(char*)
	            + nameBytes;

	if (dbgStringList)
		{
		fprintf (stderr, "sizeof(stringlist) = %u\n", (u32) sizeof(stringlist));
		fprintf (stderr, "sizeof(char*) = %u\n", (u32) sizeof(char*));
		fprintf (stderr, "bytesNeeded = %u\n", bytesNeeded);
		}

	sl = malloc (bytesNeeded);
	if (sl == NULL) goto cant_allocate;
	sl->len = numNames;

	// second pass through the file, copying items

	rewind (f);

	nameCharScan = (char*) &sl->s[numNames];
	if (dbgStringList)
		fprintf (stderr, "nameCharScan = %u\n", (u32) (nameCharScan - ((char*) sl)));
	for (nameIx=0; nameIx<numNames ; nameIx++)
		{
		fgets (line, sizeof(line), f);

		// trim blanks, end of line, and comments, and ignore blank lines

		len = strlen(line);
		if (len == 0) continue;
		if (line[len-1] == '\n') line[--len] = 0;
		waffle = strchr (line, '#');
		if (waffle != NULL) *waffle = 0;
		trim_string (line);
		if (line[0] == 0) continue;

		// process the line
		// $$$ also incorporate path

		if (dbgStringList)
			{
			fprintf (stderr, "[%u]\n", nameIx);
			fprintf (stderr, "nameCharScan = %u\n", (u32) (nameCharScan - ((char*) sl)));
			}

		sl->s[nameIx] = nameCharScan;
		nameLen = strlen(line);
		strcpy (nameCharScan, line);
		nameCharScan += nameLen + 1;
		}

	fclose (f);

	if (dbgFilenames)
		{
		u32 nameIx;

		for (nameIx=0 ; nameIx<sl->len ; nameIx++)
			fprintf (stderr, "[%u] \"%s\"\n", nameIx, sl->s[nameIx]);
		}

	return sl;

	// failures

cant_open_file:
	fprintf (stderr, "failed to open \"%s\"\n", listFilename);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

split_line:
	fprintf (stderr, "line is too long (%s: line %u)\n", listFilename, lineNum-1);
	exit(EXIT_FAILURE);
	return NULL; // (never reaches here)

cant_allocate:
	fprintf (stderr, "failed to allocate %u bytes for %u filenames\n",
					 bytesNeeded, numNames);
	exit (EXIT_FAILURE);
	return NULL; // (never reaches here)
	}

//----------
//
// print_distances--
//	Print the distances between two bitvectors.
//
//----------
//
// Arguments:
//	FILE*			f:				a file to print to
//	subsetspecs*	subsetSpecs:	specifiers for the subsets being compared
//	bitcounts*		countsVector:	list of the distances for each subset
//
// Returns:
//	(nothing)
//
//----------

void print_distances
   (FILE*			f,
	subsetspecs*	subsetSpecs,
	bitcounts*		countsVector)
	{
	u32				xorCount, andCount, orCount;
	u32				subsetIx;
	subset*			ss;

	for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
		{
		ss = &subsetSpecs->subset[subsetIx];

		xorCount = countsVector[subsetIx].xorCount;
		andCount = countsVector[subsetIx].andCount;
		orCount  = countsVector[subsetIx].orCount;

		if (dbgShowAndOr)
			fprintf (f, "\t%u/%u\t%u\t%u", xorCount, ss->numBits, andCount, orCount);
		else if (dbgIntDistance)
			fprintf (f, "\t%u/%u", xorCount, ss->numBits);
		else
			fprintf (f, "\t%.6f", ((float) xorCount) / ss->numBits);
		}
	}

//----------
//
// print_densities--
//	Print the densities of a bitvector.
//
//----------
//
// Arguments:
//	FILE*			f:				a file to print to
//	subsetspecs*	subsetSpecs:	specifiers for the subsets being analyzed
//	bitcounts*		countsVector:	list of the density for each subset
//
// Returns:
//	(nothing)
//
//----------

void print_densities
   (FILE*			f,
	subsetspecs*	subsetSpecs,
	bitcounts*		countsVector)
	{
	u32				count;
	u32				subsetIx;
	subset*			ss;

	for (subsetIx=0 ; subsetIx<subsetSpecs->len ; subsetIx++)
		{
		ss = &subsetSpecs->subset[subsetIx];

		count = countsVector[subsetIx].xorCount;
		if (dbgIntDistance)
			fprintf (f, "\t%u/%u", count, ss->numBits);
		else
			fprintf (f, "\t%.6f", ((float) count) / ss->numBits);
		}
	}

