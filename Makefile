CC=gcc
IDIR = include
ODIR = obj
SRC = src
CFLAGS=-I ~/.local/include -L ~/.local/lib64 -I include
LIBS=-lhts -lpthread
DEPS = $(IDIR)/functions.h
OBJ = $(ODIR)/functions.o $(ODIR)/main.o

all: bin/katcher bin/list_kmers

bin/katcher: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

bin/list_kmers: src/list_kmers.c
	gcc src/list_kmers.c -o bin/list_kmers

$(ODIR)/%.o: $(SRC)/%.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

test: bin/katcher
	bin/katcher -i ~/sv_gwas/usda_lines/merged_bams/SRR1533395/SRR1533395_merged.bam -k SRR1533395_all_kmers.txt -o SRR1533395_reads.bam -t6 -b 1200000 -m 100 

