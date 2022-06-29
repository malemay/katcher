#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <htslib/khash.h>
#include <htslib/sam.h>

// Declaring a hash table for storing the read IDs that are to be kept
KHASH_SET_INIT_STR(qtable);

// This program takes as input a bam file with a subset of reads
// of interest and another bam file with the complete set of reads
// from which the subset was obtained. The objective is to retrieve
// all the reads in the whole set that have the QNAME of any read
// in the subset. The idea is to eventually retrieve both reads in
// signficant pairs in order to use them for assembly.

void create_hash(samFile *input, sam_hdr_t *header, khash_t(qtable) *hash_table);

// Usage: ./extract_qname <subset.bam> <complete.bam> > <output.bam>
int main(int argc, char* argv[]) {

	// Declaring htslib-related variables
	samFile *subset_bam = hts_open(argv[1], "r"), *complete_bam = hts_open(argv[2], "r"), *output_bam = hts_open("-", "wb");
	sam_hdr_t *subset_header = sam_hdr_init(), *complete_header = sam_hdr_init();
	bam1_t *alignment = bam_init1();

	// Declaring khash-related variables
	khash_t(qtable) *hash_table = kh_init(qtable);
	khiter_t k;

	// Reading the headers
	subset_header = sam_hdr_read(subset_bam);
	complete_header = sam_hdr_read(complete_bam);

	// Populating the hash table from the subset bam
	create_hash(subset_bam, subset_header, hash_table);
	k = kh_end(hash_table);

	// Writing the header to the output bam
	if(sam_hdr_write(output_bam, subset_header) != 0) {
		fprintf(stderr, "ERROR: Unable to write header to output file. Aborting\n");
		exit(1);
	}

	// Iterating over the complete bam file
	while(sam_read1(complete_bam, complete_header, alignment) >= 0) {
		// Checking whether this QNAME is found in the subset, and writing if so
		if(kh_get(qtable, hash_table, bam_get_qname(alignment)) != k) {
			if(sam_write1(output_bam, subset_header, alignment) < 0) {
				fprintf(stderr, "ERROR: Unable to write to output file. Aborting\n");
				exit(1);
			}
		}
	}

	// Closing the files and freeing memory
	sam_close(subset_bam);
	sam_close(complete_bam);
	sam_close(output_bam);

	for (k = kh_begin(hash_table); k != kh_end(hash_table); k++)
		if (kh_exist(hash_table, k)) free((char*)kh_key(hash_table, k));

	kh_destroy(qtable, hash_table);
}

void create_hash(samFile *input, sam_hdr_t *header, khash_t(qtable) *hash_table) {
	// Declaring variables
	int khash_return;
	bam1_t *alignment = bam_init1();

	// Iterating over the alignments
	while(sam_read1(input, header, alignment) >= 0) {
		kh_put(qtable, hash_table, strdup(bam_get_qname(alignment)), &khash_return);
	}
}

