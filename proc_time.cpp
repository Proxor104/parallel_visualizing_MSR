#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include "headers.h"


double get_cpu_time()
	{
		struct rusage buf;
		getrusage(RUSAGE_THREAD, &buf);
		return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec/1.e6;
	}


double get_full_time() 
	{
		  struct timeval buf;
		  gettimeofday(&buf, 0);
		  return buf.tv_sec + buf.tv_usec / 1.e6;
	}
