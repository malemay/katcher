CC = gcc
IDIR = include
ODIR = obj
SRC = src
CFLAGS = -O2 -I include -I ~/.local/include -L ~/.local/lib -L ~/.local/lib64
LIBS = -lhts -lpthread
DEPS = $(IDIR)/functions.h
OBJ = $(ODIR)/functions.o $(ODIR)/main.o

all: bin/katcher bin/list_kmers bin/add_pvalues bin/extract_qname

bin/katcher: $(OBJ)
	mkdir -p bin
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^

bin/add_pvalues: src/add_pvalues.c src/functions.c include/functions.h
	mkdir -p bin
	$(CC) $(CFLAGS) $(LIBS) -o bin/add_pvalues src/add_pvalues.c src/functions.c

bin/list_kmers: src/list_kmers.c
	mkdir -p bin
	$(CC) $(CFLAGS) $(LIBS) src/list_kmers.c -o bin/list_kmers

bin/extract_qname: src/extract_qname.c
	mkdir -p bin
	$(CC) $(CFLAGS) $(LIBS) -o bin/extract_qname src/extract_qname.c

$(ODIR)/%.o: $(SRC)/%.c $(DEPS)
	mkdir -p obj
	$(CC) $(CFLAGS) $(LIBS) -c -o $@ $<

test: bin/katcher
	bin/list_kmers test/significant_pav_table.txt SRR1533395 > test/input_kmers.txt
	bin/katcher -i test/input_test.bam -k test/input_kmers.txt -o test/reads_with_kmers.bam -t6 -b 600000 -m 300 
	bin/add_pvalues test/reads_with_kmers.bam test/kmer_pvalues.txt > test/reads_with_pvalues.bam
	bin/extract_qname test/reads_with_pvalues.bam test/input_test.bam > test/significant_pairs.bam

