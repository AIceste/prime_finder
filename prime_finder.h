#ifndef PRIME_FINDER_H
#define PRIME_FINDER_H

#include <gmp.h>
#include <memory>

struct pf_interval {
	mpz_t lower_bound, upper_bound;	
	size_t prime_count;
	mpz_t *primes;
};

struct pf_instance {
	size_t thread_count;
	size_t block_size;
	size_t interval_count;
	struct pf_interval *intervals;
};

void pf_process_interval(
	pf_instance const *const inst, mpz_t iterator, size_t interval_index
);

#endif
