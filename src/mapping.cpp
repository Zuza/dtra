// #include "mapping.h"

// using namespace std;

// void performMapping(MappingResult* mappingResult,
// 		    const Database& db, shared_ptr<Read> read) {
//   unsigned hsh = 0;
//   unsigned andMask = (1<<(2*db.getSeedLen()))-1;

//   // probaj normalni
//   // TODO: reverse complement

//   for (int i = 0; i < read->data.size(); ++i) {
//     int next = baseToInt(getnome->data()[i]);
//     hsh = hsh*4+next;
//     hsh &= andMask;

//     if (i+1 >= db.getSeedLen()) {
//       // nekakvo brojanje
//     }
//   }
// }
