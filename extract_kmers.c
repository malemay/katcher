#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "sam.h"

// The amount of memory allocated for the read sequence
#define READ_ALLOC 500

// The length of the k-mers
#define KMER_LENGTH 31

// Creating a global table used to translate 4-bit integers into chars
char nt_table[16] = {'*', 'A', 'C', '*',
                    'G', '*', '*', '*',
                    'T', '*', '*', '*',
                    '*', '*', '*', 'N'};

// This program takes a bam file and a k-mer as input, and outputs all
// reads containing this kmer or its reverse complement in SAM format
// The output SAM file has no header in the current implementation

// Usage extract_kmers <in.bam> <kmer_sequence> > out.sam

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
		fprintf(stderr, "Usage: %s <in.bam> <kmer_sequence> > out.sam\n", argv[0]);
		return 1;
	}

	// Processing the command-line arguments
	char *input_file = argv[1];
	char *kmer_forward = argv[2], *kmer_reverse = NULL;

	// Declaring simple variables
	int write_op, read_op;
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

	// Creating a reverse-complemented version of the k-mer
	kmer_reverse = (char*) malloc((KMER_LENGTH + 1) * sizeof(char));
	if(revcomp(kmer_forward, kmer_reverse) != 0) return 2;

	// Loop over the alignments in the bam file
	while((read_op = sam_read1(input, header, alignment)) > 0) {
		// Get a pointer to the sequence and the length of the sequence
		bamseq = bam_get_seq(alignment);
		seqlength = alignment->core.l_qseq;

		// Filling the seq memory segment with the character sequence
		bamseq_to_char(bamseq, seq, seqlength);

		// Output to file if the k-mer is found within the sequence
		if(strstr(seq, kmer_forward) != NULL) {
			write_op = sam_write1(output, header, alignment);
		}

		if(strstr(seq, kmer_reverse) != NULL) {
			write_op = sam_write1(output, header, alignment);
		}

	}

	// Making sure that we reached the EOF
	assert(read_op == -1);
	
	// Freeing the file handles and memory allocated
	sam_close(input);
	sam_close(output);

	sam_hdr_destroy(header);
	bam_destroy1(alignment);

	free(kmer_reverse);
	free(seq);

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
				fprintf(stderr, "Error! Unknown nucleotide found in k-mer sequence");
				return 1;
		}

		return 0;
	}

	// Ending the array with a NULL character to indicate the end of the string
	reverse[KMER_LENGTH] = '\0';
}

