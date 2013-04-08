#ifndef INNOCENTIVE_MAPPING
#define INNOCENTIVE_MAPPING

#include "database.h"
#include "read.h"
#include <memory>

struct MappingResult {
  int databaseId;
  int pos;
  bool isReverseComplemented;
  double score;
};

void performMapping(MappingResult* mappingResult,
		    const Database& db, std::shared_ptr<Read> read);

#endif
