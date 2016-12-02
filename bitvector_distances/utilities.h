#ifndef utilities_H				// (prevent multiple inclusion)
#define utilities_H

// sizing types

#include <inttypes.h>
typedef int8_t    s8;
typedef uint8_t   u8;
typedef int32_t   s32;
typedef uint32_t  u32;
typedef int64_t   s64;
typedef uint64_t  u64;

#define u32Max ((u32) -1)

// macro to convince gnu c compiler not to complain about unusued function
// arguments

#ifdef __GNUC__
#define arg_dont_complain(arg) arg __attribute__ ((unused))
#else
#define arg_dont_complain(arg) arg
#endif // __GNUC__

// establish ownership of global variables

#ifdef utilities_owner
#define global
#else
#define global extern
#endif

//----------
//
// prototypes for functions in this module--
//
//----------

char* copy_string            (const char* s);
char* copy_prefix            (const char* s, u32 n);
int   strcmp_prefix          (const char* str1, const char* str2);
char* trim_string            (char* s);
u32   string_to_u32          (const char* s);
u32   string_to_unitized_u32 (const char* s, int byThousands);

//----------
//
// miscellany
//
//----------

#define round_up_16(b)  ((((u64) (b))+15)&(~15))
#define round_up_1K(b)  ((((u64) (b))+1023)&(~1023))

#undef global
#endif // utilities_H
