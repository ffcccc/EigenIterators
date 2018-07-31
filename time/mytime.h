#if !defined(timeutil_H__INCLUDED_)
#define timeutil_H__INCLUDED_
//	For an MFC app, use VC_EXTRALEAN, otherwise use WIN32_LEAN_AND_MEAN. Define them before including the header files
//	#define WIN32_LEAN_AND_MEAN
#include <iostream>
#include <stdio.h>

inline void doPause() {
	std::cout << "press any key to continue\n";
	getchar();
}

#ifdef WIN32
	#define VC_EXTRALEAN
	#define hrtime   unsigned long long
	#include <windows.h>

	inline hrtime gethrtime() {
	  return GetTickCount();
	}

#else
	#include <sys/time.h>
	#define hrtime   long int
	inline hrtime gethrtime() {
		//get the current number of microseconds since january 1st 1970
		timeval ts;
		gettimeofday(&ts,0);
		return (hrtime)(ts.tv_sec * 1000 + (ts.tv_usec / 1000));
	}

	/**
	 * C++ version 0.4 char* style "itoa":
	 * Written by Lukás Chmela
	 * Released under GPLv3.
	 */
	inline char* itoa(int value, char* result, int base) {
		// check that the base if valid
		if (base < 2 || base > 36) { *result = '\0'; return result; }

		char* ptr = result, *ptr1 = result, tmp_char;
		int tmp_value;

		do {
			tmp_value = value;
			value /= base;
			*ptr++ = "zyxwvutsrqponmlkjihgfedcba9876543210123456789abcdefghijklmnopqrstuvwxyz" [35 + (tmp_value - value * base)];
		} while ( value );

		// Apply negative sign
		if (tmp_value < 0) *ptr++ = '-';
		*ptr-- = '\0';
		while(ptr1 < ptr) {
			tmp_char = *ptr;
			*ptr--= *ptr1;
			*ptr1++ = tmp_char;
		}
		return result;
	}

#endif




#endif
