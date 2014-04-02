#ifndef INNOCENTIVE_MAPPING
#define INNOCENTIVE_MAPPING

#include "database.h"
#include "read.h"
#include "index.h"
#include "gene.h"
#include <algorithm>
#include <memory>

void performMapping(std::vector<std::shared_ptr<Gene> >& genes,
                    std::shared_ptr<Index> idx, 
                    const int seedLen,
                    std::shared_ptr<Read> read,
                    bool fillEditDistance);

#endif
