// Copyright (C) 2010, Steffen Knollmann
// Released under the terms of the GNU General Public License version 3.


#include "xmem.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>


size_t  global_allocated_bytes     = 0;
size_t  global_max_allocated_bytes = 0;
int64_t global_malloc_vs_free      = 0;


extern void *
xmalloc(size_t size)
{
	void *dummy;

#ifdef XMEM_TRACK_MEM
	dummy = malloc(size + 8);
#else
	dummy = malloc(size);
#endif
	if (dummy == NULL) {
		fprintf(stderr, "Could not allocate ");
		if (size < 1024) {
			fprintf(stderr, "%ib\n", (int)size);
		} else if (size < 1048576) {
			fprintf(stderr, "%fkb\n", size / 1024.);
		} else {
			fprintf(stderr, "%fMb\n", size / 1048576.);
		}
#ifdef XMEM_TRACK_MEM
		xmem_info(stderr);
#endif
		fprintf(stderr, "Exiting... :-(\n");
		exit(EXIT_FAILURE);
	}

#ifdef XMEM_TRACK_MEM
	global_allocated_bytes += size;
	global_malloc_vs_free++;
	if (global_allocated_bytes > global_max_allocated_bytes)
		global_max_allocated_bytes = global_allocated_bytes;
	*((uint64_t *)dummy) = (uint64_t)size;
	dummy                = (void *)(((uint64_t *)dummy) + 1);
#endif

	return dummy;
} /* xmalloc */

extern void
xfree(void *ptr)
{
#ifdef XMEM_TRACK_MEM
	if (ptr == NULL)
		return;

	if (global_malloc_vs_free <= 0) {
		fprintf(stderr,
		        "Calling free too often.\n");
		xmem_info(stderr);
		exit(EXIT_FAILURE);
	}
	global_allocated_bytes -= *(((uint64_t *)ptr) - 1);
	global_malloc_vs_free--;
	free(((uint64_t *)ptr) - 1);
#else
	free(ptr);
#endif

	return;
} /* xfree */

extern void *
xrealloc(void *ptr, size_t size)
{
	void     *dummy;
#ifdef XMEM_TRACK_MEM
	uint64_t old_size;
#endif

	if (ptr == NULL) {
		dummy = xmalloc(size);
		return dummy;
	}

#ifdef XMEM_TRACK_MEM
	dummy = realloc(((uint64_t *)ptr) - 1, size + 8);
#else
	dummy = realloc(ptr, size);
#endif
	if (dummy == NULL) {
		fprintf(stderr, "Could not re-allocate ");
		if (size < 1024) {
			fprintf(stderr, "%ib\n", (int)size);
		} else if (size < 1048576) {
			fprintf(stderr, "%fkb\n", size / 1024.);
		} else {
			fprintf(stderr, "%fMb\n", size / 1048576.);
		}
#ifdef XMEM_TRACK_MEM
		xmem_info(stderr);
#endif
		fprintf(stderr, "Exiting... :-(\n");
		exit(EXIT_FAILURE);
	}

#ifdef XMEM_TRACK_MEM
	old_size                = *((uint64_t *)dummy);
	global_allocated_bytes -= old_size;
	global_allocated_bytes += size;
	if (global_allocated_bytes > global_max_allocated_bytes)
		global_max_allocated_bytes = global_allocated_bytes;
	*((uint64_t *)dummy) = (uint64_t)size;
	dummy                = (void *)(((uint64_t *)dummy) + 1);
#endif

	return dummy;
} /* xrealloc */

#ifdef XMEM_TRACK_MEM
void
xmem_info(FILE *f)
{
	fprintf(f, "Currently holding: ");
	if (global_allocated_bytes < 1024) {
		fprintf(f, "%ib\n", (int)global_allocated_bytes);
	} else if (global_allocated_bytes < 1048576) {
		fprintf(f, "%fkb\n", global_allocated_bytes / 1024.);
	} else {
		fprintf(f, "%fMb\n", global_allocated_bytes / 1048576.);
	}
	fprintf(f, "Peak usage: ");
	if (global_max_allocated_bytes < 1024) {
		fprintf(f, "%ib\n", (int)global_max_allocated_bytes);
	} else if (global_max_allocated_bytes < 1048576) {
		fprintf(f, "%fkb\n", global_max_allocated_bytes / 1024.);
	} else {
		fprintf(f, "%fMb\n", global_max_allocated_bytes / 1048576.);
	}
	fprintf(f, "Malloc vs. free balance: %" PRIi64 "\n",
	        global_malloc_vs_free);

	return;
} /* xmem_info */

#endif
