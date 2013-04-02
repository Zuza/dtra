#ifndef INNOCENTIVE_MAPPING
#define INNOCENTIVE_MAPPING

#include "database.h"
#include "read.h"
#include <memory>

struct MappingResult {
  int a;
};

void performMapping(MappingResult* mappingResult,
		    const Database& db, std::shared_ptr<Read> read);

#endif
