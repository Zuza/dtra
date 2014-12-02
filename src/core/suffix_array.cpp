#include "suffix_array.h"

SuffixArray::SuffixArray(const char* text, size_t textLen) : text(text), textLen(textLen) {
    suffix_array_.clear();
    suffix_array_.resize(textLen);
    int err = divsufsort((const sauchar_t *)text, &suffix_array_[0], textLen);
    assert(err != -1); // no errors please
}

const int *SuffixArray::search(const char *pattern, int patternLen, int *numSolutions) {
    int firstIndex;
    *numSolutions = sa_search(
            (const sauchar_t *)text, textLen,
            (const sauchar_t *)pattern, patternLen,
            &suffix_array_[0], textLen, &firstIndex);
    assert(*numSolutions != -1);
    return &suffix_array_[firstIndex];
}
