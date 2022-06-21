#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "sam.h"

// The amount of memory allocated for the read sequence
#define READ_ALLOC 500
// The length of the kmers being queried
#define KMER_LENGTH 31

// Creating a global table used to translate 4-bit integers into chars
char n_table[16] = {'*', 'A', 'C', '*',
                    'G', '*', '*', '*',
                    'T', '*', '*', '*',
                    '*', '*', '*', 'N'};

// This program takes a bam file and a k-mer as input, and outputs all
// reads containing this kmer or its reverse complement in SAM format

// Usage extract_kmers <in.bam> <kmer_sequence> > out.sam

void intseq_to_char(uint8_t *intseq, char *seq, int seqlength);
void revcomp(char *forward, char *reverse);

int main(int argc, char* argv[]) {
	// Checking and processing the command-line arguments
	char *input_bam = argv[1];

	char *kmer_forward = kmer_forward = argv[2];

	char *kmer_reverse = (char*) malloc((KMER_LENGTH + 1) * sizeof(char));
	revcomp(kmer_forward, kmer_reverse);

	// Reading the input from the file on the command line
	samFile *bamfile = sam_open(input_bam, "r");

	// The output is written to stdout
	samFile *output = hts_open("-", "w");

	// Initializing the sam header
	sam_hdr_t *header = sam_hdr_init();

	// Initializing the alignment record (will be recycled)
	bam1_t *alignment = bam_init1();

	// The length of the sequence
	int32_t seqlength;

	// A pointer to the sequence itself
	uint8_t *intseq = NULL;

	// A pointer to a character version of the sequence
	char *seq = (char*) malloc(READ_ALLOC * sizeof(char));

	// Variables for holding the result of I/O operations
	int write_op;
	int read_op;

	// Reading the header and ouputting it as is
	header = sam_hdr_read(bamfile);
	assert(header != NULL);

	//write_op = sam_hdr_write(output, header);
	//assert(write_op == 0);

	// Loop over the alignments in the bam file
	while((read_op = sam_read1(bamfile, header, alignment)) > 0) {
		// Get a pointer to the sequence and the length of the sequence
		intseq = bam_get_seq(alignment);
		seqlength = alignment->core.l_qseq;

		// Filling the seq memory segment with the character sequence
		intseq_to_char(intseq, seq, seqlength);

		// Output to file if the k-mer is found within the sequence
		if(strstr(seq, kmer_forward) != NULL) {
			write_op = sam_write1(output, header, alignment);
		}

		if(strstr(seq, kmer_reverse) != NULL) {
			write_op = sam_write1(output, header, alignment);
		}

	}

	assert(read_op == -1);
	
	// Freeing the memory allocated
	sam_hdr_destroy(header);
	bam_destroy1(alignment);
	sam_close(bamfile);
	sam_close(output);
	free(kmer_reverse);
	free(seq);

	return 0;
}

// A function that takes a pointer to a sequence in a bam record
// (encoded with 4 bits) and uses another pointer (seq) to fill
// a seqlength-long sequence represented as characters
void intseq_to_char(uint8_t *intseq, char *seq, int seqlength) {

	int i;

	for(i = 0; i < seqlength; i++) {
		seq[i] = n_table[bam_seqi(intseq, i)];
	}

	// Appending a NULL character so that this is understood as a string
	seq[i] = '\0';

	return;
}

// A function that performs the reverse complement of a k-mer sequence
void revcomp(char *forward, char *reverse) {

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
				printf("Error! Unknown nucleotide found in k-mer sequence");
				break;
		}
	}

	// Ending the array with a NULL character to indicate the end of the string
	reverse[KMER_LENGTH] = '\0';
}

