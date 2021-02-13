#include "prime_finder.h"

// Notice that current is modified
static size_t find_primes(mpz_t current, mpz_t const &up_to, mpz_t *const buf) {
	size_t count = 0;
	unsigned int const reps = 40u;
	while (mpz_cmp(current, up_to) <= 0) {
		if (mpz_probab_prime_p(current, reps))
			mpz_init_set(buf[count++], current);
		mpz_add_ui(current, current, 1lu);
	}
	return count;
}

// Return 0 unless a processed interval was met.
void pf_process_interval(
	pf_instance const *const inst, mpz_t iterator, size_t interval_index
) {
	pf_interval &interval = inst->intervals[interval_index];
	mpz_set(iterator, interval.lower_bound);
	interval.prime_count = find_primes(
		iterator, interval.upper_bound, interval.primes
	);
}

