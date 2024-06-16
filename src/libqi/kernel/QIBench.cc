#include <libqi/kernel/QIBench.h>

#ifdef HAVE_GETRUSAGE
#include <sys/types.h>
#include <sys/resource.h>
#else
#include <time.h>
#endif

/** Routine for benchmarking with an accuracy
 *  depending on the host platform.
 *  Linux machines will benefit from the "getrusage"
 *  function, with an accuracy in milliseconds.
 *  Other architectures will simply use the "clock()"
 *  function, with an accuracy in seconds. */
int qi_bench_cputime (void)
{
	#ifdef HAVE_GETRUSAGE
	struct rusage rus;

	getrusage (0, &rus);
	return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
	#else
	#warning "qi_bench_cputime() will have an accuracy in seconds only."
	clock_t clock_s;

	clock_s = clock();
	return (int)( 1000.0 * (((double)clock_s)/CLOCKS_PER_SEC) );
	#endif
}

