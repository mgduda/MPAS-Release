#ifdef __linux__
#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#endif

void set_compute_affinity(int rank)
{
#ifdef __linux__
	cpu_set_t mask;

	CPU_ZERO(&mask);
/*
	if (rank % 35 < 17) {
		CPU_SET(rank % 34, &mask);
		fprintf(stderr, "Assigning rank %i to CPU %i\n", rank, rank % 34);
	} else {
		CPU_SET((rank % 34) + 1, &mask);
		fprintf(stderr, "Assigning rank %i to CPU %i\n", rank, (rank % 34) + 1);
	}
*/
	CPU_SET(rank % 35, &mask);
	fprintf(stderr, "Assigning rank %i to CPU %i\n", rank, rank % 35);
	sched_setaffinity(0, sizeof(cpu_set_t), &mask);	
#endif
}
