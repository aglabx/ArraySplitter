/*
 * FASTA file reading functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "arraysplitter.h"

#define INITIAL_SEQ_CAPACITY 10000
#define INITIAL_ENTRIES_CAPACITY 1000

/* Check if file is gzipped */
static bool is_gzipped(const char *filename) {
    size_t len = strlen(filename);
    return (len > 3 && strcmp(filename + len - 3, ".gz") == 0);
}

/* Read line from gzipped or regular file */
typedef struct {
    gzFile gz;
    FILE *fp;
    bool is_gz;
} FileHandle;

static FileHandle* file_open(const char *filename) {
    FileHandle *fh = malloc(sizeof(FileHandle));
    if (!fh) return NULL;

    fh->is_gz = is_gzipped(filename);

    if (fh->is_gz) {
        fh->gz = gzopen(filename, "r");
        fh->fp = NULL;
        if (!fh->gz) {
            free(fh);
            return NULL;
        }
    } else {
        fh->fp = fopen(filename, "r");
        fh->gz = NULL;
        if (!fh->fp) {
            free(fh);
            return NULL;
        }
    }

    return fh;
}

static void file_close(FileHandle *fh) {
    if (!fh) return;
    if (fh->is_gz && fh->gz) {
        gzclose(fh->gz);
    } else if (fh->fp) {
        fclose(fh->fp);
    }
    free(fh);
}

static char* file_gets(char *buf, int size, FileHandle *fh) {
    if (!fh) return NULL;
    if (fh->is_gz) {
        return gzgets(fh->gz, buf, size);
    } else {
        return fgets(buf, size, fh->fp);
    }
}

/* Read all entries from FASTA file */
FastaArray* fasta_read_file(const char *filename) {
    FileHandle *fh = file_open(filename);
    if (!fh) {
        fprintf(stderr, "Error: Cannot open file %s\n", filename);
        return NULL;
    }

    FastaArray *fa = malloc(sizeof(FastaArray));
    if (!fa) {
        file_close(fh);
        return NULL;
    }

    fa->capacity = INITIAL_ENTRIES_CAPACITY;
    fa->count = 0;
    fa->entries = malloc(fa->capacity * sizeof(FastaEntry));
    if (!fa->entries) {
        free(fa);
        file_close(fh);
        return NULL;
    }

    char line[MAX_HEADER_LENGTH];
    char *current_header = NULL;
    char *current_seq = NULL;
    size_t seq_len = 0;
    size_t seq_capacity = 0;

    while (file_gets(line, sizeof(line), fh)) {
        /* Remove trailing newline/carriage return */
        size_t len = strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) {
            line[--len] = '\0';
        }

        if (line[0] == '>') {
            /* Save previous entry if exists */
            if (current_header && current_seq) {
                /* Expand array if needed */
                if (fa->count >= fa->capacity) {
                    fa->capacity *= 2;
                    fa->entries = realloc(fa->entries, fa->capacity * sizeof(FastaEntry));
                }

                fa->entries[fa->count].header = current_header;
                fa->entries[fa->count].sequence = current_seq;
                fa->entries[fa->count].length = seq_len;
                fa->count++;
            } else {
                /* Free incomplete entries */
                free(current_header);
                free(current_seq);
            }

            /* Start new entry */
            current_header = strdup(line + 1); /* Skip '>' */
            seq_capacity = INITIAL_SEQ_CAPACITY;
            current_seq = malloc(seq_capacity);
            current_seq[0] = '\0';
            seq_len = 0;

        } else if (current_header) {
            /* Append to current sequence */
            size_t line_len = len;

            /* Expand sequence buffer if needed */
            while (seq_len + line_len + 1 > seq_capacity) {
                seq_capacity *= 2;
                current_seq = realloc(current_seq, seq_capacity);
            }

            /* Convert to uppercase and append */
            for (size_t i = 0; i < line_len; i++) {
                current_seq[seq_len++] = toupper((unsigned char)line[i]);
            }
            current_seq[seq_len] = '\0';
        }
    }

    /* Save last entry */
    if (current_header && current_seq && seq_len > 0) {
        if (fa->count >= fa->capacity) {
            fa->capacity *= 2;
            fa->entries = realloc(fa->entries, fa->capacity * sizeof(FastaEntry));
        }

        fa->entries[fa->count].header = current_header;
        fa->entries[fa->count].sequence = current_seq;
        fa->entries[fa->count].length = seq_len;
        fa->count++;
    } else {
        free(current_header);
        free(current_seq);
    }

    file_close(fh);
    return fa;
}

/* Free FASTA array */
void fasta_free(FastaArray *fa) {
    if (!fa) return;

    for (size_t i = 0; i < fa->count; i++) {
        free(fa->entries[i].header);
        free(fa->entries[i].sequence);
    }
    free(fa->entries);
    free(fa);
}
