#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "sam.h"

// The amount of memory allocated for the read sequence
#define READ_ALLOC 500

// The length of the k-mers
#define KMER_LENGTH 31

// The maximum number of k-mers that the program can process
#define MAX_KMERS 1000

// Creating a global table used to translate 4-bit integers into chars
char nt_table[16] = {'*', 'A', 'C', '*',
                    'G', '*', '*', '*',
                    'T', '*', '*', '*',
                    '*', '*', '*', 'N'};

// This program takes a bam file and a list of k-mers as input, and outputs all
// reads containing these kmers or their reverse complement in SAM format
// The output SAM file has no header in the current implementation
// The k-mer list must be a text file with one k-mer per line

// Usage extract_kmers <in.bam> <kmer_list.txt> > out.sam

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
	int write_op, read_op, n_kmers = 0;
	int32_t seqlength;
	uint8_t *bamseq = NULL;
	char *seq = (char*) malloc(READ_ALLOC * sizeof(char));

	// Declaring and initializing htslib-related variables
	samFile *input = sam_open(input_file, "r");
	assert(input != NULL);

	samFile *output = hts_open("-", "w");
	assert(output != NULL);

	sam_hdr_t *header = sam_hdr_init();
	header = sam_hdr_read(input);
	assert(header != NULL);

	bam1_t *alignment = bam_init1();

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

	// Loop over the alignments in the bam file
	while((read_op = sam_read1(input, header, alignment)) > 0) {
		// Get a pointer to the sequence and the length of the sequence
		bamseq = bam_get_seq(alignment);
		seqlength = alignment->core.l_qseq;

		// Filling the seq memory segment with the character sequence
		bamseq_to_char(bamseq, seq, seqlength);

		// Iterating over all the k-mers
		for(int i = 0; i < n_kmers; i++) {

			// Output to file if the k-mer is found within the sequence
			if(strstr(seq, forward_kmers[i]) != NULL) {
				write_op = sam_write1(output, header, alignment);
			}

			if(strstr(seq, reverse_kmers[i]) != NULL) {
				write_op = sam_write1(output, header, alignment);
			}
		}

	}

	// Making sure that we reached the EOF
	assert(read_op == -1);
	
	// Freeing the file handles and memory allocated
	sam_close(input);
	sam_close(output);
	fclose(kmer_list);

	sam_hdr_destroy(header);
	bam_destroy1(alignment);

	free(seq);

	for(int i = 0; i < n_kmers; i++) {
		free(forward_kmers[i]);
		free(reverse_kmers[i]);
	}

	free(forward_kmers);
	free(reverse_kmers);

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

