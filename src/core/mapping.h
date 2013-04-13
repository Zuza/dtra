#ifndef INNOCENTIVE_MAPPING
#define INNOCENTIVE_MAPPING

#include "database.h"
#include "read.h"
#include <memory>

struct MappingResult {
  int databaseId;
  int start, end;
  bool isReverseComplemented;
  double score;

  bool operator < (const MappingResult& other) {
    return score < other.score;
  }
};

void performMapping(MappingResult* mappingResult,
		    const Database& db, std::shared_ptr<Read> read);

#endif
