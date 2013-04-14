#ifndef INNOCENTIVE_MAPPING
#define INNOCENTIVE_MAPPING

#include "database.h"
#include "read.h"
#include "index.h"
#include <algorithm>
#include <memory>

void performMapping(std::shared_ptr<Index> idx, std::shared_ptr<Read> read);

#endif
