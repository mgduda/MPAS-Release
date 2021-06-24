#ifdef __linux__
#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#endif
#include "mpi.h"

void set_compute_affinity(int rank)
{
#ifdef __linux__
	cpu_set_t mask;
	int i;
	char filename[32];
	FILE *f;

	CPU_ZERO(&mask);
//define CPUS_PER_NODE 34
//define CPUS_PER_SOCKET 17
//define SOCKETS_PER_NODE 2

#define CPUS_PER_NODE 35
#define CPUS_PER_SOCKET 35
#define SOCKETS_PER_NODE 1

	if (rank % CPUS_PER_NODE < CPUS_PER_SOCKET) {
		CPU_SET(rank % CPUS_PER_NODE, &mask);
	} else {
		CPU_SET(rank % CPUS_PER_NODE + 1, &mask);
	}
	sched_setaffinity(0, sizeof(cpu_set_t), &mask);	

#endif
}
