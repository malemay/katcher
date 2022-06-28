#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// The size (in bytes) of the buffer into which each line is read
#define BUFSIZE 50000

// The number of fields to initialize the fields vector with
#define N_FIELDS_INIT 100

// This program takes a PAV table and a sample name as input,
// and outputs the k-mers found in that sample to stdout
int main(int argc, char* argv[]) {
	// Processing the command-line arguments
	FILE *input = fopen(argv[1], "r");
	char *sample = argv[2];

	// Creating variables
	char buffer[BUFSIZE], **fields, *token;
	int sample_field, n_col = 0, n_match = 0, n_kmers = 0, row = 1, field_size = N_FIELDS_INIT, n_fields;

	// Initializing an array of strings for the fields
	fields = (char**) malloc(field_size * sizeof(char*));

	// Reading the first line of the file
	if(fgets(buffer, BUFSIZE, input) == NULL) {
		fprintf(stderr, "Unable to read input file %s. Aborting\n", argv[1]);
		exit(1);
	}

	// Removing the newline character from the input string
	buffer[strcspn(buffer, "\n")] = '\0';

	// Looping over the fields in the first column to find the sample that matches the input
	token = strtok(buffer, "\t");

	while(token != NULL) {

		// Adding the field to the fields array
		if(n_col + 1 >= field_size) {
			field_size *= 2;
			fields = (char**) realloc(fields, field_size * sizeof(char*));
		}

		if(strcmp(token, sample) == 0) {
			n_match++;
			sample_field = n_col;
		}

		token = strtok(NULL, "\t");

		n_col++;
	}

	// Checking that one and only one column matched
	if(n_match == 0) {
		fprintf(stderr, "ERROR: sample name %s matches no column name in %s.\n", sample, argv[1]);
		exit(1);
	}

	if(n_match > 1) {
		fprintf(stderr, "ERROR: sample name %s matches more than one column name in %s.\n", sample, argv[1]);
		exit(1);
	}

	// Debug printing
	fprintf(stderr, "Sample %s matched column %d in %s\n", sample, sample_field, argv[1]);

	// Looping over the rows containing the k-mers and the presence/absence data
	while(fgets(buffer, BUFSIZE, input) != NULL) {

		// Resetting the counter for the number of fields read in the row
		n_fields = 0;

		// Removing the newline character from the input string
		buffer[strcspn(buffer, "\n")] = '\0';

		// Looping over the fields in the first column to find the sample that matches the input
		token = strtok(buffer, "\t");

		while(token != NULL) {

			// Checking that we are not exceeding the number of columns read in the first row
			if(n_fields == n_col) {
				fprintf(stderr, "ERROR: Row number %d contains more fields than expected (%d).\n", row, n_col);
			}

			// Checking that the value in that field is allowed ("0" or "1")
			if(n_fields> 0 && strcmp(token, "0") != 0 && strcmp(token, "1") != 0) {
				fprintf(stderr, "Field number %d at row number %d contains disallowed string %s.\n", n_fields, row, token);
				exit(1);
			}

			// Adding that field to the string array
			fields[n_fields] = strdup(token);

			token = strtok(NULL, "\t");

			n_fields++;
		}

		// Checking that we stopped at the right number of fields
		if(n_fields != n_col) {
			fprintf(stderr, "ERROR: Row number %d contains %d fields whereas the expected number is %d.\n", row, n_fields, n_col);
		}

		// Now that the sanity checks have been done we can output the k-mer if it is found in the sample of interest
		if(strcmp(fields[sample_field], "1") == 0) {
			puts(fields[0]);
			n_kmers++;
		}

			row++;
	}

	// Debug printing
	fprintf(stderr, "%d out of %d k-mers were found for sample %s in %s. Analysis completed successfully.\n", n_kmers, row, sample, argv[1]);

	// Closing files and freeing memory
	fclose(input);

	for(int i = 0; i < n_fields; i++) {
		free(fields[i]);
	}

	free(fields);

	return 0;
}

