// utilities.c-- miscellaneous utility functions.

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <inttypes.h>
#include <math.h>
#include <float.h>

#define  utilities_owner		// (make this the owner of its globals)
#include "utilities.h"			// interface to this module

//----------
//
// copy_string, copy_prefix--
//	Create (in the heap) a copy of a string or a prefix of a string.
//
//----------
//
// Arguments:
//	const char*	s:	The string to copy.
//	int			n:	(copy_prefix only) the number of characters to copy.
//
// Returns:
//	A pointer to new string;  failures result in fatality.
//
//----------

char* copy_string
   (const char*	s)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc (strlen(s) + 1);
	if (ss == NULL) goto cant_allocate;
	return strcpy (/*to*/ ss, /*from*/ s);

cant_allocate:
	fprintf (stderr, "failed to allocate copy of string \"%s\"\n", s);
	exit (EXIT_FAILURE);
	return NULL; // (never reaches here)
	}


char* copy_prefix
   (const char*	s,
	u32			n)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc (n + 1);
	if (ss == NULL) goto cant_allocate;
	memcpy (/*to*/ ss, /*from*/ s, /*how much*/ n);
	ss[n] = 0;
	return ss;

cant_allocate:
	fprintf (stderr, "failed to allocate %d bytes for copy of string \"%s\"\n", n, s);
	exit (EXIT_FAILURE);
	return NULL; // (never reaches here)
	}

//----------
//
// strcmp_prefix--
//	Determine if a string contains another as a prefix.
//
//----------
//
// Arguments:
//	const char*	str1:	The string.
//	const char*	str2:	The prefix string.
//
// Returns:
//	The same as strcmp(prefix1,str2) would, where prefix1 is str1 truncated
//	to be no longer than str2.
//
//----------

int strcmp_prefix
   (const char*	str1,
	const char*	str2)
	{
	return strncmp (str1, str2, strlen(str2));
	}

//----------
//
// trim_string--
//	Remove blanks (and end-of-line) from both ends of a string.
//
//----------
//
// Arguments:
//	char*	s:	The string.
//
// Returns:
//	The string (the same as s).  Leading blanks are removed by copying
//	characters forward.  Trailing blanks are removed by depositing a
//	terminating zero.
//
//----------

char* trim_string
   (char*	s)
	{
	char*	ss, *dd, *lastInk;

	// skip to first non-blank

	ss = s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;

	if 	(*ss == 0) // (string has nothing but blanks)
		{ *s = 0;  return s; }

	// copy the rest of the string (except the terminating zero)

	dd = lastInk = s;
	while (*ss != 0)
		{
		*(dd++) = *(ss++);

		if ((*ss != 0) && (*ss != ' ') && (*ss != '\t') && (*ss != '\n'))
			lastInk = dd;
		}

	// poke a terminating zero just past the last non-blank

	lastInk[1] = 0;

	return s;
	}

//----------
//
// string_to_u32--
//	Parse a string for the integer value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in program
//	termination.
//
//----------

u32 string_to_u32
   (const char*	s)
	{
	char*		ss;
	u32			v;
	char		extra;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// convert to number

	if (*ss == '-') goto not_an_integer;
	if (sscanf (ss, "%u%c", &v, &extra) != 1) goto not_an_integer;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	fprintf (stderr, "an empty string is not an unsigned integer\n");
	exit (EXIT_FAILURE);

not_an_integer:
	fprintf (stderr, "\"%s\" is not an unsigned integer\n", s);
	exit (EXIT_FAILURE);

	return 0;
	}

//----------
//
// string_to_unitized_u32--
//	Parse a string for the integer value it contains, allowing K, M, and G
//	suffixes.
//
//----------
//
// Arguments:
//	const char*	s:		The string to parse.
//	int	byThousands:	true  => K means one thousand
//						false => K means 1,024.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in program
//	termination.
//
//----------

u32 string_to_unitized_u32
   (const char*	s,
	int			byThousands)
	{
	char		ss[20];
	int			len = strlen (s);
	char*		parseMe;
	u32			v;
	float		vf;
	char		extra;
	u32			mult;
	int			isFloat;

	// convert to number

	mult = 1;

	if (len >= (int) sizeof (ss))
		parseMe = (char*) s;
	else
		{
		parseMe = ss;
		strcpy (ss, s);

		if (len > 0)
			{
			switch (ss[len-1])
				{
				case 'K': case 'k':
					mult = (byThousands)? 1000 : 1024;
					break;
				case 'M': case 'm':
					mult = (byThousands)? 1000000 : 1024L * 1024L;
					break;
				case 'G': case 'g':
					mult = (byThousands)? 1000000000 : 1024L * 1024L * 1024L;
					break;
				}

			if (mult != 1)
				ss[len-1] = 0;
			}
		}

	isFloat = false;
	if (sscanf (parseMe, "%u%c", &v, &extra) != 1)
		{
		if (sscanf (parseMe, "%f%c", &vf, &extra) != 1) goto bad;
		isFloat = true;
		}

	if (isFloat)
		{
		if  (vf < 0) goto bad;
		if ((vf > 0) && (vf*mult+.5 > u32Max)) goto overflow;
		v = vf*mult + .5;
		}
	else if (mult != 1)
		{
		if ((v > 0) && ( v > u32Max / mult)) goto overflow;
		v *= mult;
		}

	return v;

bad:
	fprintf (stderr, "\"%s\" is not an unsigned integer\n", s);
	exit (EXIT_FAILURE);

overflow:
	fprintf (stderr, "\"%s\" is out of range for an unsigned integer\n", s);
	exit (EXIT_FAILURE);
	}

