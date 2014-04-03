#ifndef INNOCENTIVE_MAPPING
#define INNOCENTIVE_MAPPING

#include "database.h"
#include "read.h"
#include "gene.h"
#include <algorithm>
#include <memory>

#include "FmIndexWavelet/DnaIndex.hpp"

void performMapping(std::vector<std::shared_ptr<Gene> >& genes,
                    std::shared_ptr<DnaIndex> idx, 
                    const int seedLen,
                    std::shared_ptr<Read> read,
                    bool fillEditDistance);

#endif
