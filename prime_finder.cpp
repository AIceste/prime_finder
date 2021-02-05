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
static int pf_thread_process_block(
	pf_instance const *const inst, mpz_t iterator, size_t interval_index
) {
	size_t const block_end = (
		 (interval_index + inst->block_size < inst->interval_count)
 		? interval_index + inst->block_size
		: inst->interval_count 
	);
	while (interval_index < block_end) {
		pf_interval &interval = inst->intervals[interval_index];
		// This here is the single race condition of the program, we could
		// try adding safety here if deemed necessary, tho per-block safety
		// would make more sense then. (1)
		if (interval.is_processed.test_and_set())
			return 1;
		mpz_set(iterator, interval.lower_bound);
		interval.prime_count = find_primes(
			iterator, interval.upper_bound, interval.primes
		);
		++interval_index;
	}
	return 0;
}

void *pf_thread(void *param) {
	pf_instance const *const inst = reinterpret_cast<pf_param*>(param)->instance;
	size_t const tid = reinterpret_cast<pf_param*>(param)->thread_id;

	// A lock-less approach is used to minimalise synchronisation overhead.
	// Each thread is assigned blocks according to their tid, from 0, meaning
	// that work distribution is entirely pre-defined by the instance and thus
	// does not need synchronisation. To prevent disparity between processed
	// intervals or CPU time given to each thread from slowing the process too
	// much by forcin it to wait after a single, overloaded thread, a work-steal
	// approach is used, as threads done with their assigned blocks try to find
	// an unfinished thread and process their block in opposite order. 

	// Process completion is detected using interval's is_processed 
	// field, which is a c++ STL atomic_flag, defined to be lockless. 

	// Various choices, such as consistently stealing work from the same thread,
	// aim at reducing the risk of collision between threads currently stealing
	// work.

	// Blocks are used for better integration to future instance initialisation
	// parallelisation.

	// For now, blocks are never subdivised to smaller units (interval), but we
	// could consider splitting the final work to further minimise the time 
	// spent waiting after the last thread.


	// Iterator for find_primes, declared here to minimise memory allocation
	// by keeping it alive accross all iterations.
	mpz_t it;
	mpz_init(it);

	// The first loop iterates through all blocks assigned to the thread or until
	// an already processed interval is met (or end of buffer, reached).
	size_t const block_interval = inst->block_size * inst->thread_count;
	size_t interval_index = inst->block_size * tid;
	while (interval_index < inst->interval_count) {
		if (pf_thread_process_block(inst, it, interval_index))
			break;
		interval_index += block_interval;
	}

	// The second loop handles the "work stealing" part of the algorithm.
	// We start by finding the beginning of the last block and which thread it
	// belongs to.
	size_t last_block = 
		inst->interval_count - inst->interval_count % inst->block_size;
	if (last_block == inst->interval_count)
		last_block -= inst->block_size;
	size_t const final_block = last_block - block_interval;

	// Want to keep going until only processed blocks remain, in which case
	// the algorithm is done and the thread can terminate.
	while (1) {
		// Search for an unfinished thread.
		interval_index = last_block;
		while (inst->intervals[interval_index].is_processed.test_and_set()
			&& (interval_index -= inst->block_size) > final_block);
		if (interval_index == final_block)
			break;

		// If an unfinished thread is found, backward-process all of its blocks.
		inst->intervals[interval_index].is_processed.clear(); // For code reuse.
		while (interval_index > inst->block_size) {
			if (pf_thread_process_block(inst, it, interval_index))
				break;
			interval_index -= block_interval;
		}
	} 

	mpz_clear(it);
	return 0;
}

