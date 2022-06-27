#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <pthread.h>

#include <htslib/sam.h>
#include <khash.h>

#include "functions.h"

// This program takes a bam file and a list of k-mers as input, and outputs all
// reads containing these kmers or their reverse complement in SAM format
// The output SAM file has no header in the current implementation
// The k-mer list must be a text file with one k-mer per line
int main(int argc, char* argv[]) {

	// Parsing the command-line arguments
	params args = {NULL, NULL, "-", 1, 100000, 300};
	if(parse_args(argc, argv, &args) != 0) return 1;

	// Processing the command-line arguments
	FILE *kmer_list = fopen(args.kmer_list, "r");
	char **forward_kmers, **reverse_kmers, **all_kmers;

	// Declaring simple variables
	int n_read, records_per_thread, executing_threads, n_kmers = 0, khash_return;
	long long int n_processed = 0;
	char **seq = (char**) malloc(args.bufsize * sizeof(char*));

	// Declaring an array containing the data processed by each thread
	thread_data *thread_chunks;
	// A struct containing the data for the decompression thread
	buffer_data bufdata;
	// Another array containing structs with data for the copying threads
	copy_data *copy_chunks;

	// Declaring the array of threads and a thread for decompressing the file
	pthread_t *threads = NULL, bufthread;

	// Declaring and initializing the hash table that will store the hashed values of the kmers, and its end iterator hash_end
	khash_t(kmer_hash) *kmer_table;
	kmer_table = kh_init(kmer_hash);
	khint_t hash_end;

	// Declaring and initializing htslib-related variables
	samFile *input = sam_open(args.input_file, "r");
	assert(input != NULL);

	samFile *output = hts_open(args.output_file, "w");
	assert(output != NULL);

	sam_hdr_t *header = sam_hdr_init();
	header = sam_hdr_read(input);
	assert(header != NULL);

        // Initializing buffers to store several alignment records
        bam1_t **bambuf = (bam1_t**) malloc(args.bufsize * sizeof(bam1_t*));
        bam1_t **tmpbuf = (bam1_t**) malloc(args.bufsize * sizeof(bam1_t*));

        for(int i = 0; i < args.bufsize; i++) {
                bambuf[i] = bam_init1();
                tmpbuf[i] = bam_init1();
                seq[i] = (char*) malloc((args.seq_alloc + 1) * sizeof(char));
        }

	// Allocating memory for the k-mer arrays and reading the k-mers from file
	forward_kmers = (char**) malloc(MAX_KMERS * sizeof(char*));
	reverse_kmers = (char**) malloc(MAX_KMERS * sizeof(char*));

	for(int i = 0; i < MAX_KMERS; i++) {
		forward_kmers[i] = (char*) malloc((KMER_LENGTH + 2) * sizeof(char));
		reverse_kmers[i] = (char*) malloc((KMER_LENGTH + 2) * sizeof(char));
	}

	while(fgets(forward_kmers[n_kmers], KMER_LENGTH + 2, kmer_list) != NULL) {
		// Removing the newline character from the string	
		forward_kmers[n_kmers][strcspn(forward_kmers[n_kmers], "\n")] = '\0';
		n_kmers++;
	}

	// Creating a reverse-complemented version of the k-mers
	for(int i = 0; i < n_kmers; i++) {
		if(revcomp(forward_kmers[i], reverse_kmers[i]) != 0) exit(1);
	}

	// Creating a common array with both forward and reverse k-mers together
	all_kmers = (char**) malloc(n_kmers * 2 * sizeof(char*));

	for(int i = 0; i < n_kmers; i++) {
		all_kmers[i * 2] = forward_kmers[i];
		all_kmers[i * 2 + 1] = reverse_kmers[i];
	}

	// Adding all k-mers to the hash table
	for(int i = 0; i < (2 * n_kmers); i++) kh_put(kmer_hash, kmer_table, all_kmers[i], &khash_return);

	// Getting the end position of the hash table
	hash_end = kh_end(kmer_table);

	// Creating an array of thread_data structs with as many elements as there are threads
	// The copy chunks which are meant for copying the data from tmpbuf to bambuf wiht multithreading are also processed here
	thread_chunks = (thread_data*) malloc(args.n_threads * sizeof(thread_data));
	copy_chunks = (copy_data*) malloc(args.n_threads * sizeof(copy_data));
	records_per_thread = args.bufsize / args.n_threads;

	// Initializing the chunk elements
	for(int i = 0; i < args.n_threads; i++) {
		thread_chunks[i].kmer_table = kmer_table;
		thread_chunks[i].hash_end = hash_end;
		thread_chunks[i].output = output;
		thread_chunks[i].header = header;
		thread_chunks[i].bambuf = bambuf + i * records_per_thread;
		thread_chunks[i].seq = seq + i * records_per_thread;
		thread_chunks[i].n_records = records_per_thread; // this value might change after filling the buffer

		copy_chunks[i].tmpbuf = tmpbuf + i * records_per_thread;
		copy_chunks[i].bambuf = bambuf + i * records_per_thread;
		copy_chunks[i].seq = seq + i * records_per_thread;
		copy_chunks[i].n_records = records_per_thread; // we won't care about adjusting this at the end of file for now
		copy_chunks[i].max_length = args.seq_alloc;
	}

	// Creating the array of analysis threads
	threads = (pthread_t*) malloc(args.n_threads * sizeof(pthread_t));

	// Initializing the lock
	if(pthread_mutex_init(&lock, NULL) != 0) {
		fprintf(stderr, "mutex init has failed\n");
		exit(1);
	}

	// Filling the temporary buffer by preparing the data and running the associated thread
	bufdata.bamfile = input;
	bufdata.header = header;
	bufdata.tmpbuf = tmpbuf;
	bufdata.max_read = args.bufsize;

	pthread_create(&bufthread, NULL, fillbuf, &bufdata);
	pthread_join(bufthread, NULL);

	n_read = bufdata.n_read;

	// Loop over the alignments in the bam file
        while(n_read > 0) {

		// Copying data from the temporary buffer to the one used for analysis
		// With multithreading
		for(int i = 0; i < args.n_threads; i++) pthread_create(&threads[i], NULL, copy_buf, &copy_chunks[i]);
		for(int i = 0; i < args.n_threads; i++) pthread_join(threads[i], NULL);

		// We no longer need the data in tmpbuf so we can launch the decompression thread
		pthread_create(&bufthread, NULL, fillbuf, &bufdata);

		// Generate the character string representation of the sequences in the bam file
		for(int i = 0; i < args.n_threads; i++) pthread_create(&threads[i], NULL, prepare_seq, &copy_chunks[i]);
		for(int i = 0; i < args.n_threads; i++) pthread_join(threads[i], NULL);

		// If we filled the buffer then the number of executing threads is args.n_threads
		if(n_read == args.bufsize) {
			executing_threads = args.n_threads;
		} else {
			// Otherwise there will be a certain number of threads that process
			// the usual number of records, and the last thread will process
			// the remaining records
			if(n_read % records_per_thread == 0) {
				executing_threads = n_read / records_per_thread;
			} else {
				executing_threads = n_read / records_per_thread + 1;
			}

			fprintf(stderr, "Buffer did not fill completely.\n"
					"Remaining records processed by %d threads processing "
					"%d records and one thread processing %d records.\n",
					executing_threads - 1,
					records_per_thread, 
					n_read - (executing_threads - 1) * records_per_thread);

			thread_chunks[executing_threads - 1].n_records = n_read - (executing_threads - 1) * records_per_thread;
		}

		// Launching the analysis threads
		for(int i = 0; i < executing_threads; i++) pthread_create(&threads[i], NULL, match_kmers, &thread_chunks[i]);
		for(int i = 0; i < executing_threads; i++) pthread_join(threads[i], NULL);

		n_processed += n_read;
		if(n_processed % 100000 == 0) fprintf(stderr, "%lld reads processed using %d threads\n", n_processed, executing_threads);

		// We wait for the buffer to be filled before starting the loop again
		pthread_join(bufthread, NULL);
		n_read = bufdata.n_read;
	}

	// Freeing the file handles and memory allocated
	sam_close(input);
	sam_close(output);
	fclose(kmer_list);

	sam_hdr_destroy(header);

	for(int i = 0; i < args.bufsize; i++) {
		bam_destroy1(bambuf[i]);
		bam_destroy1(tmpbuf[i]);
		free(seq[i]);
	}

	for(int i = 0; i < n_kmers; i++) {
		free(forward_kmers[i]);
		free(reverse_kmers[i]);
	}

	free(forward_kmers);
	free(reverse_kmers);
	free(all_kmers);
	free(thread_chunks);

	kh_destroy(kmer_hash, kmer_table);

	return 0;
}

