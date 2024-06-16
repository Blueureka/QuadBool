#ifndef _qi_bench_h_
#define _qi_bench_h_


extern short QI_CPU_TIME_MS;

/** Routine for benchmarking with an accuracy
 *  depending on the host platform.
 *  Linux machines will benefit from the "getrusage"
 *  function, with an accuracy in milliseconds.
 *  Other architectures will simply use the "clock()"
 *  function, with an accuracy in seconds. */
int qi_bench_cputime (void);

#endif
