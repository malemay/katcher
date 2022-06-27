CC=gcc
IDIR = include
ODIR = obj
SRC = src
CFLAGS=-I ~/.local/include -L ~/.local/lib64 -I include
LIBS=-lhts -lpthread
DEPS = $(IDIR)/functions.h
OBJ = $(ODIR)/functions.o $(ODIR)/main.o

katcher: $(OBJ)
		$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.o: $(SRC)/%.c $(DEPS)
		$(CC) -c -o $@ $< $(CFLAGS)

test: katcher
	./katcher -i ~/sv_gwas/usda_lines/merged_bams/SRR1533395/SRR1533395_merged.bam -k SRR1533395_all_kmers.txt -o SRR1533395_reads.sam -t6 -b 1200000 -m 100 

