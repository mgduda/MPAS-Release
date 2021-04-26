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

#if 0
	snprintf(filename, 32, "task_mask.%4.4i", rank);
	f = fopen(filename, "w");

        sched_getaffinity(0, sizeof(cpu_set_t), &mask);

MPI_Barrier(MPI_COMM_WORLD);

	fprintf(f, "  from OS:");
	for (i = 0; i < 176; i++) {
		if (i % 4 == 0) {
			fprintf(f, " ");
		}
		if (CPU_ISSET(i, &mask)) {
			fprintf(f, "1");
		} else {
			fprintf(f, "0");
		}
	}
	fprintf(f, "\n");
	fprintf(f, "          ");
	for (i = 0; i < 44; i++) {
		fprintf(f, " %2i  ", i);
	}
	fprintf(f, "\n\n");

#if 0
	CPU_ZERO(&mask);
	if (rank % 38 < 19) {
		CPU_SET(4*(rank % 38), &mask);
	} else {
		CPU_SET(4*((rank % 38) + 3), &mask);
	}

	fprintf(f, "requested:");
	for (i = 0; i < 176; i++) {
		if (i % 4 == 0) {
			fprintf(f, " ");
		}
		if (CPU_ISSET(i, &mask)) {
			fprintf(f, "1");
		} else {
			fprintf(f, "0");
		}
	}
	fprintf(f, "\n");
	fprintf(f, "          ");
	for (i = 0; i < 44; i++) {
		fprintf(f, " %2i  ", i);
	}
	fprintf(f, "\n\n");
	
	sched_setaffinity(0, sizeof(cpu_set_t), &mask);	

        sched_getaffinity(0, sizeof(cpu_set_t), &mask);

	fprintf(f, "      got:");
	for (i = 0; i < 176; i++) {
		if (i % 4 == 0) {
			fprintf(f, " ");
		}
		if (CPU_ISSET(i, &mask)) {
			fprintf(f, "1");
		} else {
			fprintf(f, "0");
		}
	}
	fprintf(f, "\n");
	fprintf(f, "          ");
	for (i = 0; i < 44; i++) {
		fprintf(f, " %2i  ", i);
	}
	fprintf(f, "\n");
#endif

	fclose(f);
#endif
#endif
}
