# ArraySplitter C - Makefile

CC = gcc
CFLAGS = -Wall -Wextra -O3 -pthread
LDFLAGS = -lz -lpthread

# Source files for arraysplitter_c
SRCS = main.c fasta.c sequence.c fstree.c decompose.c
OBJS = $(SRCS:.c=.o)

# Source files for arraysplitter_classify_c
CLASSIFY_SRCS = classify_main.c classify_kmer.c classify_cluster.c
CLASSIFY_OBJS = $(CLASSIFY_SRCS:.c=.o)

# Source files for arraysplitter_minhash_c (cascade clustering)
MINHASH_SRCS = minhash_main.c minhash_sig.c minhash_lsh.c minhash_cluster.c
MINHASH_OBJS = $(MINHASH_SRCS:.c=.o)

# Target binaries
TARGET = arraysplitter_c
CLASSIFY_TARGET = arraysplitter_classify_c
MINHASH_TARGET = arraysplitter_minhash_c

# Default target - build all
all: $(TARGET) $(CLASSIFY_TARGET) $(MINHASH_TARGET)

# Link arraysplitter_c
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Link arraysplitter_classify_c
$(CLASSIFY_TARGET): $(CLASSIFY_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Link arraysplitter_minhash_c
$(MINHASH_TARGET): $(MINHASH_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile arraysplitter sources
%.o: %.c arraysplitter.h
	$(CC) $(CFLAGS) -c $< -o $@

# Compile classify sources (depend on classify.h)
classify_%.o: classify_%.c classify.h
	$(CC) $(CFLAGS) -c $< -o $@

# Compile minhash sources (depend on minhash.h)
minhash_%.o: minhash_%.c minhash.h
	$(CC) $(CFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(CLASSIFY_OBJS) $(MINHASH_OBJS) $(TARGET) $(CLASSIFY_TARGET) $(MINHASH_TARGET)

# Install to project bin directory
install: all
	mkdir -p ../../bin
	cp $(TARGET) $(CLASSIFY_TARGET) $(MINHASH_TARGET) ../../bin/

# Install to user bin (requires sudo for /usr/local/bin)
install-user: all
	cp $(TARGET) $(CLASSIFY_TARGET) $(MINHASH_TARGET) ~/.local/bin/ 2>/dev/null || \
	cp $(TARGET) $(CLASSIFY_TARGET) $(MINHASH_TARGET) ~/bin/ 2>/dev/null || \
	echo "Install manually to your PATH"

# Debug build
debug: CFLAGS = -Wall -Wextra -g -O0 -pthread -DDEBUG
debug: clean all

# Static analysis
check:
	cppcheck --enable=all --std=c11 $(SRCS) $(CLASSIFY_SRCS) $(MINHASH_SRCS)

# Format code
format:
	clang-format -i $(SRCS) $(CLASSIFY_SRCS) $(MINHASH_SRCS) arraysplitter.h classify.h minhash.h

.PHONY: all clean install install-user debug check format
