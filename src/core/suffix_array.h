#ifndef DTRA_SUFFIX_ARRAY_DATABASE
#define DTRA_SUFFIX_ARRAY_DATABASE

#include "gene.h"
#include <divsufsort.h>

#include <vector>

class SuffixArray {
public:
    SuffixArray(const char* text, size_t textLen);

    /**
    * @param pattern to search for
    * @param patternLen
    * @param numSolutions number of found patterns in the text
    * @returns the pointer to the index of the first occurrence; the second occurrence is on *(pointer+1) and so on
    */
    const int* search(const char* pattern, int patternLen, int* numSolutions);

private:
    const char* text;
    size_t textLen;
    std::vector<saidx_t> suffix_array_;
};

#endif

