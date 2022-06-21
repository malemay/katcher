GCC=/usr/bin/gcc

extract_kmers: extract_kmers.c
	$(GCC) -O3 -I/home/malem420/.local/include/htslib -L/home/malem420/.local/lib64 -lhts -o extract_kmers extract_kmers.c
