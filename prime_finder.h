#ifndef PRIME_FINDER_H
#define PRIME_FINDER_H

#include <gmp.h>
#include <string.h>
#include <memory>
#include <atomic>

struct pf_interval {
	std::atomic_flag is_processed = ATOMIC_FLAG_INIT;
	mpz_t lower_bound, upper_bound;	
	size_t prime_count;
	mpz_t *primes;

	// For our very specific need, we need an efficient copy
	// operator that does not care about keeping is_processed
	// stater. Same goes for primes.
	pf_interval &operator=(pf_interval const &other) {
		memcpy(
			&lower_bound, &other.lower_bound,
			offsetof(pf_interval, primes) - offsetof(pf_interval, lower_bound)
		);
		return *this;
	}
	// Same goes for copy constructor
	pf_interval(pf_interval const &other) : is_processed(), primes() {
		*this = other;
	}
	// The idea is to use it as a standard C struct .-.
	pf_interval() : is_processed() {
		memset(
			&lower_bound, 0, 	
			offsetof(pf_interval, primes) - offsetof(pf_interval, lower_bound)
		);
	}
};

struct pf_instance {
	size_t thread_count;
	size_t block_size;
	size_t interval_count;
	struct pf_interval *intervals;
};

struct pf_param {
	size_t thread_id;
	struct pf_instance *instance;
};

void *pf_thread(void*);

#endif
