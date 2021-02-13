#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <assert.h>
#include <math.h>

#include <omp.h>

#include "prime_finder.h"
#include "Chrono.hpp"

constexpr int base = 10;
constexpr size_t block_size = 1;

static void pf_instance_destroy(pf_instance &instance) {
	for (size_t i = 0; i < instance.interval_count; ++i) {
		auto &interval = instance.intervals[i];
		mpz_clear(interval.lower_bound);
		mpz_clear(interval.upper_bound);
		for (size_t j = 0; j < interval.prime_count; ++j)
			mpz_clear(interval.primes[j]);
	}
	free(instance.intervals[0].primes);
}

// Technically for our use case it would be faster to write and
// destroy at once, but much cleaner if separated.
static void pf_instance_print(pf_instance const &instance) {
	// Try to output as best as possible, too bad if there's an error.
	// At this point we don't save anything by detecting them.
	for (size_t i = 0; i < instance.interval_count; ++i) {
		auto const &interval = instance.intervals[i];
		for (size_t j = 0; j < interval.prime_count; ++j) {
			// Might want to make sure using both stdout and std::cout
			// is compatible with all systems.*
			mpz_out_str(stdout, base, interval.primes[j]);
			std::cout << std::endl;
		}
	}
	fflush(stdout);
}

static int pf_instance_read(
	pf_instance &instance, std::vector<pf_interval> &buf,
	int const argc, char const *const *const argv
) {
	long thread_count;
	FILE *input_file;
	pf_interval tmp;

	if (argc != 3)
		goto usage;
	thread_count = strtol(argv[1], nullptr, base);
	if (thread_count <= 0)
		goto usage;
	
	if (!(input_file = fopen(argv[2], "r")))
		goto file;
	
	mpz_init(tmp.lower_bound);
	if (!mpz_inp_str(tmp.lower_bound, input_file, base)) {
		mpz_clear(tmp.lower_bound);
		goto file;
	}
	while (true) {
		mpz_init(tmp.upper_bound);
		if (!mpz_inp_str(tmp.upper_bound, input_file, base)) {
			mpz_clear(tmp.upper_bound);
			mpz_clear(tmp.lower_bound);
			goto file;
		}
		
		// Will reallocate as seen fit
		buf.push_back(tmp);
	
		// Done later to reuse buf cleanup code
		if (mpz_cmp(tmp.lower_bound, tmp.upper_bound) > 0)
			goto file;

		mpz_init(tmp.lower_bound);
		if (!mpz_inp_str(tmp.lower_bound, input_file, base)) {
			mpz_clear(tmp.lower_bound);
			break;
		}	
	}

	// Don't care about errors at this point.
	fclose(input_file);

	instance.thread_count = thread_count;
	instance.block_size = block_size;
	instance.interval_count = buf.size();
	instance.intervals = &buf[0];

	return 0;

usage:
	std::cerr << "Usage: first parameter is number of thread, second is input file." << std::endl;
	goto error;
file:
	std::cerr << "There was an issue handling input file." << std::endl
		      << "Make sure it exists and consist lines describing integer ranges." << std::endl;
error:
	// In case we allocated stuff, nah.
	for (auto i = buf.begin(); i != buf.end(); ++i) {
		mpz_clear(i->lower_bound);
		mpz_clear(i->upper_bound);
	}
	buf.clear();

	return 1;
}

static void print_stuff(pf_interval *const intervals, size_t const count) {
	std::cout << "INTERVALS" << std::endl;
	for (size_t i = 0; i < count; ++i) {
		std::cout << "[";
		mpz_out_str(stdout, base, intervals[i].lower_bound);
		std::cout << ",";
		mpz_out_str(stdout, base, intervals[i].upper_bound);
		std::cout << "]" << std::endl;
	}
}

// Return middle-ish of the split array
static size_t split(pf_interval *const intervals, size_t const count) {
	pf_interval tmp;
	size_t const split = count - 1;
	size_t left = 0;
	for (size_t j = 0; j < count; ++j) {
		if (mpz_cmp(
				intervals[j].lower_bound,
				intervals[split].lower_bound
			) < 0
		) {
			tmp = intervals[left];
			intervals[left++] = intervals[j];
			intervals[j] = tmp;
		}
	}
	tmp = intervals[left];
	intervals[left] = intervals[split];
	intervals[split] = tmp;
	return left;
}

// Basically quicksort
static void sort(pf_interval *const intervals, size_t const count) {
	if (count <= 1)
		return;
	size_t const center = split(intervals, count);
	if (center)
		sort(intervals, center);
	sort(intervals + center + 1, count - center - 1);
}

// Works on sorted array of intervals
static void trim(pf_interval *const intervals, size_t const count) {
	// Used to keep a reference to last useful range, allows
	// handling case with e.g. [2, 5], [5, 5].
	mpz_t *last = &intervals[0].upper_bound;
	for (size_t i = 1; i < count; ++i) {
		if (mpz_cmp(*last, intervals[i].lower_bound) >= 0) {
			mpz_add_ui(intervals[i].lower_bound, *last, 1);
			if (mpz_cmp(
					intervals[i].lower_bound, intervals[i].upper_bound
				) <= 0
			)
				last = &intervals[i].upper_bound;
		}
	}
}

// Takes an instance, a beginning and a count.
// Return number of interval sorted and trimmed.
static size_t pf_instance_sort_and_trim(
	pf_instance &instance, size_t const start, size_t count
) {
	assert(start < instance.interval_count);
	
	pf_interval *const intervals = instance.intervals + start;
	if (instance.interval_count - start < count)
		count = instance.interval_count - start;

	// Of course, we could optimise that by merging all algorithmes
	// and having less loops.. but not.
	sort(intervals, count);
	trim(intervals, count);
	print_stuff(intervals, count);
	return count;
}


// Non-zero return on error 
static int pf_instance_preallocate(pf_instance &instance) {
	// The algo assumse that intervals are small enough to
	// be stored in memory, which really should hold true.

	pf_interval *const intervals = instance.intervals;
	mpz_t upper_count;
	mpz_init(upper_count);
	size_t total_count = 0;

	// First calculate an overestimated count.
	for (size_t i = 0; i < instance.interval_count; ++i) {
		int const cmp = mpz_cmp(
			intervals[i].lower_bound, intervals[i].upper_bound
		);
		if (cmp > 0) {
			intervals[i].prime_count = 0;
			continue;
		}
		if (cmp == 0) {
			intervals[i].prime_count = 1;
			++total_count;
			continue;
		}

		// Idea based of the prime number theorem*.
		// *Number of prime from 2 to N ~= N/ln(N)
		// There is a small error value, but since primes
		// are denser on low numbers, the method is sure to 
		// work with only the +1 in the formula.
		// We could get enhance precision by playing around with
		// gmp logarithms, but the latency would get increased so
		// much it's worth spending the memory.
		mpz_set(upper_count, intervals[i].upper_bound);
		mpz_sub(upper_count, upper_count, intervals[i].lower_bound);
		double const range_size = mpz_get_d(upper_count) + 1;		
		intervals[i].prime_count = (
			(range_size == 2.)
			? 1 
			: static_cast<size_t>(
				range_size / log(range_size) + 1
			)
		);

		assert(total_count + intervals[i].prime_count > total_count);
		total_count += intervals[i].prime_count;
	}
	mpz_clear(upper_count);

	// Allocate buffer (or try to)
	intervals[0].primes = (mpz_t*)malloc(sizeof(mpz_t) * total_count);
	if (!intervals[0].primes) {
		for (size_t i = 0; i < instance.interval_count; ++i)
			intervals[i].prime_count = 0;
		return 1;
	}

	// Then assign each interval to its buffer. 
	for (size_t i = 1; i < instance.interval_count; ++i) {
		pf_interval *const prev = &intervals[i -1];
		intervals[i].primes = prev->primes + prev->prime_count;
		prev->prime_count = 0;
	}	

	return 0;
}


int main(int argc, char **argv) {
	std::vector<pf_interval> buf;
	pf_instance instance = {};
	Chrono clock;
	double duration;

	if (pf_instance_read(instance, buf, argc, argv))
		exit(EXIT_FAILURE);

	// Start timer here
	duration = clock.get();

	// Prepare instance for process by threads. 
	// This is the part that creates lag and could be parallelized 
	// given enough time.
	pf_instance_sort_and_trim(instance, 0, instance.interval_count);
	pf_instance_preallocate(instance);

	print_stuff(instance.intervals, instance.interval_count);
	exit(0);

	mpz_t iterator;
	size_t interval;

	#pragma omp parallel private(iterator, interval) shared(instance) num_threads(instance.thread_count)
	{
		mpz_init(iterator);
		
		# pragma omp for schedule(guided, instance.block_size)
		for (interval = 0; interval < instance.interval_count; ++interval) {
			pf_process_interval(&instance, iterator, interval);
		}
	}

	// Stop timer here
	duration = clock.get() - duration;

	pf_instance_print(instance);
	std::cerr << duration << std::endl;

	pf_instance_destroy(instance);
	exit(EXIT_SUCCESS);
}

