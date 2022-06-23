#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <pthread.h>

#include "sam.h"
#include "khash.h"

// The amount of memory allocated for the read sequence
#define READ_ALLOC 300

// The length of the k-mers
#define KMER_LENGTH 31

// The maximum number of k-mers that the program can process
#define MAX_KMERS 100000

// The maximum number of records that can be recorded in the buffer at one time
#define BUFSIZE 1000000

// The number of threads that the program can use
#define N_THREADS 5

// Initializing a hash table type named "kmer_hash" using klib's khash functions
KHASH_SET_INIT_STR(kmer_hash);

// This program takes a bam file and a list of k-mers as input, and outputs all
// reads containing these kmers or their reverse complement in SAM format
// The output SAM file has no header in the current implementation
// The k-mer list must be a text file with one k-mer per line

// Usage extract_kmers <in.bam> <kmer_list.txt> > out.sam
//
// Creating a global table used to translate 4-bit integers into chars
char nt_table[16] = {'*', 'A', 'C', '*',
                    'G', '*', '*', '*',
                    'T', '*', '*', '*',
                    '*', '*', '*', 'N'};

// A struct that holds the data required for a single thread to execute
typedef struct thread_data {
	khash_t(kmer_hash) *kmer_table; // Pointer to the hash table containing the k-mers
	khint_t hash_end; // the end of the hash table ; will not change after the table has been computed
	samFile *output; // Pointer to the sam output file
	sam_hdr_t *header; // Pointer to the header
	bam1_t **bambuf; // Pointer to the buffered records being processed
	int n_records; // the number of records that will be processed by this thread
	char **seq; // Pointer to the sequence being queried
} thread_data;

// A struct that holds the data for the decompression thread (the one that reads the bam file into memory)
typedef struct buffer_data {
	samFile *bamfile; // Pointer to the bam file that is being read from
	sam_hdr_t *header; // The header of the bam file being read
	bam1_t **tmpbuf; // The temporary buffer from which data will be copied to the main buffer
	int max_read; // The maximum number of records to read from the file
	int n_read; // The number of records that were read
} buffer_data;

// A struct that holds the data for copying the records from the temporary buffer to the main buffer
typedef struct copy_data {
	bam1_t **tmpbuf; // A pointer to the temporary buffer
	bam1_t **bambuf; // A pointer to the main buffer that is being processed
	int n_records; // The number of records that this thread must process
} copy_data;

// A global lock used to prevent threads from writing to the output file simultaneously
pthread_mutex_t lock;

// A function that takes the read sequence as encoded by 4-bit nucleotides
// in the BAM format and converts it to a string
// bamseq: a pointer to the start of the sequence in the BAM file
// seq: a pointer to an array of characters large enough to hold the read sequence
// seqlength: the length of the read sequence (in nucleotides)
void bamseq_to_char(uint8_t *bamseq, char *seq, int seqlength);

// A function that generates the reverse-complemented version of a k-mer
// forward: a pointer to an array of characters holding the forward k-mer
// reverse: a pointer to an array of characters that will hold the reverse k-mer (must be sufficiently large)
int revcomp(char *forward, char *reverse);

// A function that allows matching read sequence and k-mers with multithreading
void *match_kmers(void *input);

// A function that fills a temporary buffer with bam records
void *fillbuf(void* input);

// A function that copies alignments from the temporary buffer to the main one
void *copy_buf(void *input);

int main(int argc, char* argv[]) {

	// Checking the input
	if(argc != 3) {
		fprintf(stderr, "Usage: %s <in.bam> <kmer_list.txt> > out.sam\n", argv[0]);
		return 1;
	}

	// Processing the command-line arguments
	char *input_file = argv[1];
	FILE *kmer_list = fopen(argv[2], "r");
	char **forward_kmers, **reverse_kmers, **all_kmers;

	// Declaring simple variables
	int n_read, records_per_thread, executing_threads, n_kmers = 0, khash_return;
	long long int n_processed = 0;
	int32_t seqlength;
	uint8_t *bamseq = NULL;
	char **seq = (char**) malloc(BUFSIZE * sizeof(char*));

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
	samFile *input = sam_open(input_file, "r");
	assert(input != NULL);

	samFile *output = hts_open("-", "w");
	assert(output != NULL);

	sam_hdr_t *header = sam_hdr_init();
	header = sam_hdr_read(input);
	assert(header != NULL);

        // Initializing buffers to store several alignment records
        bam1_t **bambuf = (bam1_t**) malloc(BUFSIZE * sizeof(bam1_t*));
        bam1_t **tmpbuf = (bam1_t**) malloc(BUFSIZE * sizeof(bam1_t*));

        for(int i = 0; i < BUFSIZE; i++) {
                bambuf[i] = bam_init1();
                tmpbuf[i] = bam_init1();
                seq[i] = (char*) malloc(READ_ALLOC * sizeof(char));
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
		fprintf(stderr, "Read kmer %s\n", forward_kmers[n_kmers]);
		n_kmers++;
	}

	// Creating a reverse-complemented version of the k-mers
	for(int i = 0; i < n_kmers; i++) {
		if(revcomp(forward_kmers[i], reverse_kmers[i]) != 0) return 2;
		fprintf(stderr, "Generated reverse-complemented sequence %s for: %s\n", reverse_kmers[i], forward_kmers[i]);
	}

	// Creating a common array with both forward and reverse k-mers together
	all_kmers = (char**) malloc(n_kmers * 2 * sizeof(char*));

	for(int i = 0; i < n_kmers; i++) {
		all_kmers[i * 2] = forward_kmers[i];
		all_kmers[i * 2 + 1] = reverse_kmers[i];
	}

	// Adding all k-mers to the hash table
	for(int i = 0; i < (2 * n_kmers); i++) {
		kh_put(kmer_hash, kmer_table, all_kmers[i], &khash_return);
	}

	// Getting the end position of the hash table
	hash_end = kh_end(kmer_table);

	/* DEBUG to test the hash table
	// The string to check is the third argument on the command line
	if(kh_get(kmer_hash, kmer_table, argv[3]) == kh_end(kmer_table)) {
		fprintf(stderr, "k-mer sequence %s is not in k-mer hash table\n", argv[3]);
	} else {
		fprintf(stderr, "k-mer sequence %s is in k-mer hash table\n", argv[3]);
	}

	return 0;
	 DEBUG */

	/* DEBUG
	for(int i = 0; i < (n_kmers * 2); i++) {
		fprintf(stderr, "all_kmers[%d]: %s\n", i, all_kmers[i]);
	}
	DEBUG */

	// Setting the data for multithreading
	assert(BUFSIZE % N_THREADS == 0);

	// Creating an array of thread_data structs with as many elements as there are threads
	// The copy chunks which are meant for copying the data from tmpbuf to bambuf wiht multithreading are also processed here
	thread_chunks = (thread_data*) malloc(N_THREADS * sizeof(thread_data));
	copy_chunks = (copy_data*) malloc(N_THREADS * sizeof(copy_data));
	records_per_thread = BUFSIZE / N_THREADS;

	// Initializing the chunk elements
	for(int i = 0; i < N_THREADS; i++) {
		thread_chunks[i].kmer_table = kmer_table;
		thread_chunks[i].hash_end = hash_end;
		thread_chunks[i].output = output;
		thread_chunks[i].header = header;
		thread_chunks[i].bambuf = bambuf + i * records_per_thread;
		thread_chunks[i].seq = seq + i * records_per_thread;
		thread_chunks[i].n_records = records_per_thread; // this value might change after filling the buffer

		copy_chunks[i].tmpbuf = tmpbuf + i * records_per_thread;
		copy_chunks[i].bambuf = bambuf + i * records_per_thread;
		copy_chunks[i].n_records = records_per_thread; // we won't care about adjusting this at the end of file for now
	}

	/* DEBUG
	for(int i = 0; i < n_chunks; i++) {
		fprintf(stderr, "Chunk #%d:\n", i);
		for(int j = 0; j < thread_chunks[i].n_kmers; j++) {
			fprintf(stderr, "%s\n", thread_chunks[i].kmers[j]);
		}
	}

	return 0;
	DEBUG */

	// Creating the array of threads
	threads = (pthread_t*) malloc(N_THREADS * sizeof(pthread_t));

	// Initializing the lock
	if(pthread_mutex_init(&lock, NULL) != 0) {
		fprintf(stderr, "mutex init has failed\n");
		return 1;
	}

	// Filling the temporary buffer by preparing the data and running the associated thread
	bufdata.bamfile = input;
	bufdata.header = header;
	bufdata.tmpbuf = tmpbuf;
	bufdata.max_read = BUFSIZE; // The maximum number of records to read from the file

	pthread_create(&bufthread, NULL, fillbuf, &bufdata);
	pthread_join(bufthread, NULL);

	n_read = bufdata.n_read;

	// Loop over the alignments in the bam file
        while(n_read > 0) {

		// Copying data from the temporary buffer to the one used for analysis
		// With multithreading
		for(int i = 0; i < N_THREADS; i++) {
			pthread_create(&threads[i], NULL, copy_buf, &copy_chunks[i]);
		}

		// Joining the copying threads before going any further
		for(int i = 0; i < N_THREADS; i++) {
			pthread_join(threads[i], NULL);
		}

		// We no longer need the data in tmpbuf so we can launch the decompression thread
		pthread_create(&bufthread, NULL, fillbuf, &bufdata);

		// We need to prepare the seq array so that the threads can use it
		for(int i = 0; i < n_read; i++) {
			// Get a pointer to the sequence and the length of the sequence
			bamseq = bam_get_seq(bambuf[i]);
			seqlength = bambuf[i]->core.l_qseq;

			// Filling the seq memory segment with the character sequence
			bamseq_to_char(bamseq, seq[i], seqlength);
		}

		// If we filled the buffer then the number of executing threads is N_THREADS
		if(n_read == BUFSIZE) {
			executing_threads = N_THREADS;
		} else {
			// Otherwise there will be a certain number of threads that process
			// the usual number of records, and the last thread will process
			// the remaining records
			if(n_read % records_per_thread == 0) {
				executing_threads = n_read / records_per_thread;
			} else {
				executing_threads = n_read / records_per_thread + 1;
			}

			fprintf(stderr, "Buffer did not fill completely.");
			fprintf(stderr, "Rest of file will be processed by %d threads processing %d records each and one thread processing %d records.\n",
					executing_threads - 1, records_per_thread, (n_read - (executing_threads - 1) * records_per_thread));
			thread_chunks[executing_threads - 1].n_records = (n_read - (executing_threads - 1) * records_per_thread);
		}

		// Launching the threads
		for(int i = 0; i < executing_threads; i++) {
			pthread_create(&threads[i], NULL, match_kmers, &thread_chunks[i]);
		}

		// Joining the threads
		for(int i = 0; i < executing_threads; i++) {
			pthread_join(threads[i], NULL);
		}

		n_processed += n_read;
		if(n_processed % 1000000 == 0) fprintf(stderr, "%lld reads processed using %d threads\n", n_processed, executing_threads);

		// We wait for the buffer to be filled before starting the loop again
		pthread_join(bufthread, NULL);
		n_read = bufdata.n_read;
	}

	// Freeing the file handles and memory allocated
	sam_close(input);
	sam_close(output);
	fclose(kmer_list);

	sam_hdr_destroy(header);

	for(int i = 0; i < BUFSIZE; i++) {
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

// bamseq_to_char
void bamseq_to_char(uint8_t *bamseq, char *seq, int seqlength) {

	int i;

	for(i = 0; i < seqlength; i++) {
		seq[i] = nt_table[bam_seqi(bamseq, i)];
	}

	// Appending a NULL character so that this is understood as a string
	seq[i] = '\0';

	return;
}

// revcomp
int revcomp(char *forward, char *reverse) {

	for(int i = 0; i < KMER_LENGTH; i++) {
		switch (forward[i]) {
			case 'A':
				reverse[KMER_LENGTH - 1 - i] = 'T';
				break;

			case 'T':
				reverse[KMER_LENGTH - 1 - i] = 'A';
				break;

			case 'G':
				reverse[KMER_LENGTH - 1 - i] = 'C';
				break;

			case 'C':
				reverse[KMER_LENGTH - 1 - i] = 'G';
				break;

			default:
				fprintf(stderr, "Error! Unknown nucleotide found in k-mer sequence %s\n", forward);
				return 1;
		}


	}

	// Ending the array with a NULL character to indicate the end of the string
	reverse[KMER_LENGTH] = '\0';

	return 0;
}

// A function that checks whether the k-mers in a given read are found in the k-mer hash table
// This function is set for multithreading
// CAUTION: this function modifies the seq char** by setting the last character to a \0 character until the whole
// read was processed. This is OK in the current implementation because each record is accessed only by one thread
// However, this should be taken into account if the implementation were to change
void *match_kmers(void *input) {
	thread_data* chunk = (thread_data*) input;
	int write_op, seqlength;
	char *j, *seqstart;

	// Iterating over all the records
	for(int i = 0; i < chunk->n_records; i++) {

		// Getting the length of the read and a pointer to its start
		seqlength = strlen(chunk->seq[i]);
		seqstart = chunk->seq[i];

		// Iterating over all the k-mers in the read
		// j is the last character of the read
		for(j = seqstart + seqlength - 1; j > seqstart + KMER_LENGTH - 2; j--) {

			// Checking if the k-mer until the end of the read is found in the k-mer hash table
			if(kh_get(kmer_hash, chunk->kmer_table, j - KMER_LENGTH + 1) != chunk->hash_end) {
				pthread_mutex_lock(&lock);
				write_op = sam_write1(chunk->output, chunk->header, chunk->bambuf[i]);
				pthread_mutex_unlock(&lock);
			}

			// Setting the last character in the string to \0 so we can process the next k-mer
			*j = '\0';
		}
	}
}

// A function that fills a buffer with bam records and counts the number of records that were read
// This function is set such that it can be run as an independent thread that will decompress
// the input while the records that have already being read are being processed
void *fillbuf(void* input) {
	buffer_data *data = (buffer_data*) input;
        int read_op;

	// Resetting the number of records read
        data->n_read = 0;

	// Reading the records
        while(data->n_read < data->max_read && (read_op = sam_read1(data->bamfile, data->header, data->tmpbuf[data->n_read])) > 0) {
                data->n_read++;
        }
}

// A function that copies alignment records from the temporary buffer to the main one
void *copy_buf(void *input) {
	copy_data *data = (copy_data*) input;
	bam1_t *result;
	for(int i = 0; i < data->n_records; i++) {
		result = bam_copy1(data->bambuf[i], data->tmpbuf[i]);
	}
}

