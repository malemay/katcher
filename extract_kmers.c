#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <pthread.h>

#include "sam.h"

// The amount of memory allocated for the read sequence
#define READ_ALLOC 300

// The length of the k-mers
#define KMER_LENGTH 31

// The maximum number of k-mers that the program can process
#define MAX_KMERS 5000

// The maximum number of records that can be recorded in the buffer at one time
#define N_RECORDS 1000000

// The number of threads that the program can use
#define N_THREADS 5

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
// kmers: a pointer to the first k-mer associated with this thread
// n_kmers: the number of k-mers to process on this thread
// output: a pointer to the output SAM file where the results are written
// header: a pointer to the header of the input BAM file
// bambuffer: a pointer to the buffer of alignment records
// seq: a pointer to an array of character strings containing the nucleotide sequence of the reads in the buffer
typedef struct kmer_data {
	// Constant members to be set before traversing the file
	char **kmers;
	int n_kmers;
	samFile *output;
	sam_hdr_t *header;
	bam1_t **bambuffer;
	char **seq;

	// Member that must be updated each time the buffer is read
	int n_records;
} kmer_data;

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

// A function that allows matching sequence and k-mers with multithreading
void *match_kmers(void *input);

// A function that fills a buffer with bam records and returns the number of records read
int fillbuf(bam1_t **buffer, samFile *bamfile, sam_hdr_t *header, int max_read);

int main(int argc, char* argv[]) {

	// Checking the input
	if(argc != 3) {
		fprintf(stderr, "Usage: %s <in.bam> <kmer_list.txt> > out.sam\n", argv[0]);
		return 1;
	}

	// Processing the command-line arguments
	char *input_file = argv[1];
	FILE *kmer_list = fopen(argv[2], "r");
	char **forward_kmers, **reverse_kmers;

	// Declaring simple variables
	int write_op, read_op, n_chunks, kmers_per_chunk, n_kmers = 0;
	long long int n_processed = 0;
	int32_t seqlength;
	uint8_t *bamseq = NULL;
	char **seq = (char**) malloc(N_RECORDS * sizeof(char*));

	// Declaring an array containing the data processed by each thread
	kmer_data *thread_chunks;

	// Declaring and initializing htslib-related variables
	samFile *input = sam_open(input_file, "r");
	assert(input != NULL);

	samFile *output = hts_open("-", "w");
	assert(output != NULL);

	sam_hdr_t *header = sam_hdr_init();
	header = sam_hdr_read(input);
	assert(header != NULL);

        // Initializing a buffer to store several alignment records
        bam1_t **bambuffer = (bam1_t**) malloc(N_RECORDS * sizeof(bam1_t*));

        for(int i = 0; i < N_RECORDS; i++) {
                bambuffer[i] = bam_init1();
        }
	
	// Also initializing a buffer containing string sequences of the reads
        for(int i = 0; i < N_RECORDS; i++) {
                seq[i] = (char*) malloc(READ_ALLOC * sizeof(char));
        }
	

	// Declaring the array of threads
	pthread_t *threads = NULL;

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

	// Setting the data for multithreading
	assert((n_kmers * 2) % N_THREADS == 0);

	// Creating a common array with both forward and reverse k-mers together
	char **all_kmers = (char**) malloc(n_kmers * 2 * sizeof(char*));

	for(int i = 0; i < n_kmers; i++) {
		all_kmers[i * 2] = forward_kmers[i];
		all_kmers[i * 2 + 1] = reverse_kmers[i];
	}

	/* DEBUG
	for(int i = 0; i < (n_kmers * 2); i++) {
		fprintf(stderr, "all_kmers[%d]: %s\n", i, all_kmers[i]);
	}
	DEBUG */

	// Creating an array of kmer_data structs with as many elements as there are threads
	kmers_per_chunk = (n_kmers * 2) / N_THREADS;
	n_chunks = (n_kmers * 2) / kmers_per_chunk;
	thread_chunks = (kmer_data*) malloc(n_chunks * sizeof(kmer_data));

	// Initializing the chunk elements
	for(int i = 0; i < n_chunks; i++) {
		thread_chunks[i].kmers = all_kmers + i * kmers_per_chunk;
		thread_chunks[i].n_kmers = kmers_per_chunk;
		thread_chunks[i].output = output;
		thread_chunks[i].header = header;
		thread_chunks[i].bambuffer = bambuffer;
		thread_chunks[i].seq = seq;
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

	// Loop over the alignments in the bam file
        while((read_op = fillbuf(bambuffer, input, header, N_RECORDS)) > 0) {

		// We need to prepare the seq array so that the threads can use it
		for(int i = 0; i < read_op; i++) {
			// Get a pointer to the sequence and the length of the sequence
			bamseq = bam_get_seq(bambuffer[i]);
			seqlength = bambuffer[i]->core.l_qseq;

			// Filling the seq memory segment with the character sequence
			bamseq_to_char(bamseq, seq[i], seqlength);
		}

		// Setting the number of records for all chunks
		for(int i = 0; i < N_THREADS; i++) {
			thread_chunks[i].n_records = read_op;	
		}

		// Launching the threads
		for(int i = 0; i < N_THREADS; i++) {
			pthread_create(&threads[i], NULL, match_kmers, &thread_chunks[i]);
		}

		// Joining the threads
		for(int i = 0; i < N_THREADS; i++) {
			pthread_join(threads[i], NULL);
		}

		n_processed += read_op;
		if(n_processed % 100000 == 0) fprintf(stderr, "%lld reads processed\n", n_processed);
	}

	// Making sure that we reached the EOF
	assert(read_op == -1);
	
	// Freeing the file handles and memory allocated
	sam_close(input);
	sam_close(output);
	fclose(kmer_list);

	sam_hdr_destroy(header);

	for(int i = 0; i < N_RECORDS; i++) {
		bam_destroy1(bambuffer[i]);
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

// A function that checks the match between a sequence and a set of k-mers using multi-threading
void *match_kmers(void *input) {
	kmer_data* chunk = (kmer_data*) input;
	int write_op;

	// Iterating over all the records
	for(int i = 0; i < chunk->n_records; i++) {

		// Iterating over all the k-mers
		for(int j = 0; j < chunk->n_kmers; j++) {

			// Output to file if the k-mer is found within the sequence
			if(strstr(chunk->seq[i], chunk->kmers[j]) != NULL) {
				pthread_mutex_lock(&lock);
				write_op = sam_write1(chunk->output, chunk->header, chunk->bambuffer[i]);
				pthread_mutex_unlock(&lock);
			}
		}
	}
}

// A function that fills a buffer with bam records and returns the number of records read
int fillbuf(bam1_t **buffer, samFile *bamfile, sam_hdr_t *header, int max_read) {
        int read_op;
        int n_read = 0;

        while(n_read < max_read && ((read_op = sam_read1(bamfile, header, buffer[n_read])) > 0)) {
                n_read++;
        }

        if(n_read > 0) {
                return n_read;
        } else {
                return -1;
        }
}

