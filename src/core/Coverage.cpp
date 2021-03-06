#include "Coverage.h"
#include <cassert>
#include <algorithm>
#include <map>
#include <set>
using namespace std;

namespace {

  enum EventType { LEFT, RIGHT };
  const int kInfinity = 1000000000;

  struct Event {
    int x;
    int value;
    EventType type;
    int id;

    Event(int x, int value, EventType type, int id) :
      x(x), value(value), type(type), id(id) {
      // printf("pusham x=%d, value=%d, type=%s, id=%d\n", x, value,
      // 	     type == LEFT ? "left": "right", id);
    }

    bool operator < (const Event& a) const {
      if (x != a.x) return x < a.x;
      if (type != a.type) return type < a.type;
      if (value != a.value) return value < a.value;
      return id < a.id;
    }
  };

  // cvorovi dijele prostor po indeksima, pamtim 2 kljucne velicine:
  // a) najbolje rjesenje za intervale kojima sam obradio desni kraj, 
  //    tkzv. 'zatvorene' ili 'closed' intervale
  // b) najbolje rjesenje za intervale kojima sam obradio lijevi, ali jos
  //    nisam desni kraj, tkzv. 'otvorene' (opened) intervale
  class IntervalTree {
  public:
    IntervalTree(int n) {
      for (offset_ = 1; offset_ < n; offset_ *= 2);
      ivec_ = vector<Node> (offset_*2);
    }
    
    struct Result {
      int openedMax;
      int openedMaxLastId;

      int closedMax;
      int closedMaxLastId;
      
      Result() {}
      Result(int openedMax, int openedMaxLastId,
	     int closedMax, int closedMaxLastId) :
	openedMax(openedMax), openedMaxLastId(openedMaxLastId),
	closedMax(closedMax), closedMaxLastId(closedMaxLastId) {}
    };

    Result mergeResults(const Result& l, const Result& r) {
      Result combined;
      
      combined.openedMax = 
	(l.openedMax > r.openedMax ? l.openedMax : r.openedMax);
      combined.openedMaxLastId =
	(l.openedMax > r.openedMax ? l.openedMaxLastId : r.openedMaxLastId);

      combined.closedMax =
	(l.closedMax > r.closedMax ? l.closedMax : r.closedMax);
      combined.closedMaxLastId = 
	(l.closedMax > r.closedMax ? l.closedMaxLastId : r.closedMaxLastId);
      return combined;
    }
    
    void insertOpened(int value, int id, int x, int best) {
      NodeItem ni(best-x, x, value, id);
      insertOpened(0, 0, offset_-1, ni);
      idToNodeItem_[id] = ni;
    }

    Result get(int value) {
      return get(0, 0, offset_-1, value);
    }
    
    void convertOpenedToClosed(int id, int closingX) {
      const NodeItem& toRemove = idToNodeItem_[id];
      NodeItem toInsert = toRemove;
      toInsert.best += closingX+1;
      toInsert.x = closingX;
      // printf("CONVERTING id=%d, best: %d -> %d, value: %d->%d\n", 
      // 	     id, toRemove.best, toInsert.best,
      // 	     toRemove.value, toInsert.value);
      convertOpenedToClosed(0, 0, offset_-1, toRemove, toInsert);
    }
    
  private:

    struct NodeItem {
      int best;
      int x;
      int value;
      int id;

      NodeItem(int best = 0, int x = 0, int value = 0, int id = 0) : 
	best(best), x(x), value(value), id(id) {}
      
      bool operator < (const NodeItem& a) const {
	if (best != a.best) { return best < a.best; }
	if (x != a.x) { return x < a.x; }
	if (value != a.value) { return value < a.value; }
	return id < a.id;
      }
    };
    
    struct Node {
      multiset<NodeItem> opened; // best = najbolje do sad - lijevi kraj
      multiset<NodeItem> closed; // best = najbolje do sad
      
      int openedMax;
      int openedMaxLastId;

      int closedMax;
      int closedMaxLastId;

      Node() {
	openedMax = -kInfinity;
	openedMaxLastId = -1;
	closedMax = 0;
	closedMaxLastId = -1;
      }

      Result getResult() {
	return Result(openedMax, openedMaxLastId,
		      closedMax, closedMaxLastId);
      }
    };

    Result get(int node, int lo, int hi, int value) {
      if (value >= hi) return ivec_[node].getResult();
      if (lo > value) return Result(-kInfinity, -1, 0, -1);
      
      int mid = (lo+hi)/2;
      Result l = get(node*2+1, lo, mid, value);
      Result r = get(node*2+2, mid+1, hi, value);
      return mergeResults(l, r);
    }

    void updateNode(int node) { // non-leaf node
	// update opened
	if (ivec_[node*2+1].openedMax > ivec_[node*2+2].openedMax) {
	  ivec_[node].openedMax = ivec_[node*2+1].openedMax;
	  ivec_[node].openedMaxLastId = ivec_[node*2+1].openedMaxLastId;
	} else {
	  ivec_[node].openedMax = ivec_[node*2+2].openedMax;
	  ivec_[node].openedMaxLastId = ivec_[node*2+2].openedMaxLastId;
	}

	// update closed
	// printf("node=%d, lcm: %d, rcm: %d\n",
	//        node,
	//        ivec_[node*2+1].closedMax,
	//        ivec_[node*2+2].closedMax);
	if (ivec_[node*2+1].closedMax > ivec_[node*2+2].closedMax) {
	  ivec_[node].closedMax = ivec_[node*2+1].closedMax;
	  ivec_[node].closedMaxLastId = ivec_[node*2+1].closedMaxLastId;
	  // printf("1: %d %d\n", 
	  // 	 ivec_[node].closedMax, ivec_[node].closedMaxLastId);
	} else {
	  ivec_[node].closedMax = ivec_[node*2+2].closedMax;
	  ivec_[node].closedMaxLastId = ivec_[node*2+2].closedMaxLastId;
	  // printf("2: %d %d\n", 
	  // 	 ivec_[node].closedMax, ivec_[node].closedMaxLastId);
	}
    }

    void insertOpened(int node, int lo, int hi, const NodeItem& ni) {
      if (lo > ni.value) { return; }
      if (hi < ni.value) { return; }

      if (lo == hi) {
	ivec_[node].opened.insert(ni);
	ivec_[node].openedMax = ivec_[node].opened.rbegin()->best;
	ivec_[node].openedMaxLastId = ivec_[node].opened.rbegin()->id;

	if (!ivec_[node].closed.empty()) {
	  ivec_[node].closedMax = ivec_[node].closed.rbegin()->best;
	  ivec_[node].closedMaxLastId = ivec_[node].closed.rbegin()->id;
	} else {
	  ivec_[node].closedMax = 0;
	  ivec_[node].closedMaxLastId = -1;
	}
      } else {
	int mid = (lo+hi)/2;
	insertOpened(node*2+1, lo, mid, ni);
	insertOpened(node*2+2, mid+1, hi, ni);
	updateNode(node);
      }
    }

    void convertOpenedToClosed(int node, int lo, int hi,
			       const NodeItem& toRemove,
			       const NodeItem& toInsert) {
      int value = toRemove.value;
      if (lo > value) { return; }
      if (hi < value) { return; }

      if (lo == hi) {
	auto iter = ivec_[node].opened.find(toRemove);
#ifdef DEBUG
	assert(iter != ivec_[node].opened.end());
#endif
	ivec_[node].opened.erase(iter);
	ivec_[node].closed.insert(toInsert);

	if (!ivec_[node].opened.empty()) {
	  ivec_[node].openedMax = ivec_[node].opened.rbegin()->best;
	  ivec_[node].openedMaxLastId = ivec_[node].opened.rbegin()->id;
	} else {
	  ivec_[node].openedMax = -kInfinity;
	  ivec_[node].openedMaxLastId = -1;
	}
	ivec_[node].closedMax = ivec_[node].closed.rbegin()->best;
	ivec_[node].closedMaxLastId = ivec_[node].closed.rbegin()->id;
	// printf("updateao [%d, %d] na opened=%d,%d closed=%d,%d\n", 
	//        lo, hi,
	//        ivec_[node].openedMax, ivec_[node].openedMaxLastId,
	//        ivec_[node].closedMax, ivec_[node].closedMaxLastId);
      } else {
	int mid = (lo+hi)/2;
	convertOpenedToClosed(node*2+1, lo, mid, toRemove, toInsert);
	convertOpenedToClosed(node*2+2, mid+1, hi, toRemove, toInsert);
	updateNode(node);
	// printf("updateao [%d, %d] na opened=%d,%d closed=%d,%d\n", 
	//        lo, hi,
	//        ivec_[node].openedMax, ivec_[node].openedMaxLastId,
	//        ivec_[node].closedMax, ivec_[node].closedMaxLastId);
      }
    }

  private:
    map<int, NodeItem> idToNodeItem_;
    vector<Node> ivec_;
    int offset_;
  };

};

struct CoverSlowCmp {
  bool operator () (const Interval& a, const Interval& b) const {
    return a.left < b.left;
  }
};

void reconstruct(vector<Interval>* resultCoverage,
		 const vector<int>& dpRecon,
		 const int endIndex,
		 const vector<Interval>& intervals) {
  int leftBorder = -1;
  int rightBorder = -1;
  int minVal = -1;
    
  for (int curr = endIndex; curr != -1; ) {
    //printf("%d %d\n", curr, dpRecon[curr]);

    bool pushback = false;
    if (leftBorder == -1) {
      minVal = intervals[curr].value;
      leftBorder = intervals[curr].left;
      rightBorder = intervals[curr].right;
      curr = dpRecon[curr];
    } else if (intervals[curr].right >= leftBorder) {
      minVal = intervals[curr].value;
      leftBorder = intervals[curr].left;
      curr = dpRecon[curr];
    } else {
      pushback = true;
    }
    //printf("%d %d %d %d %d\n", curr, leftBorder, rightBorder, minVal, pushback);

    if (curr == -1 || pushback) {
      resultCoverage->push_back(Interval(leftBorder,rightBorder,minVal));
      leftBorder = -1;
      rightBorder = -1;
      minVal = -1;
    }
  }

  sort(resultCoverage->begin(), resultCoverage->end(), CoverSlowCmp());
}

void coverFast(vector<Interval>* resultCoverage,
	       int* result,
	       const vector<Interval>& intervals) {
  vector<Event> events;
  int maxValue = 0;
  int n = intervals.size();

  for (size_t i = 0; i < n; ++i) {
    events.push_back(Event(intervals[i].left, intervals[i].value, LEFT, i));
    events.push_back(Event(intervals[i].right, intervals[i].value, RIGHT, i));
    maxValue = max(maxValue, intervals[i].value);
  }

  sort(events.begin(), events.end());

  // OPASKA: ukoliko ce biti potrebno, ova se funkcija
  // moze dodatno ubrzati tako da se 'value' u intervalima
  // na ulazu sazme -> tad je potrebno i sazeti prije sweepanja
  // i vratiti prave vrijednosti kod rekonstrukcije
  IntervalTree itree(maxValue+1);
  *result = 0;
  vector<int> dpRecon(n, -1);
  int endIndex = -1;

  for (size_t i = 0; i < events.size(); ++i) {
    int value = events[i].value;
    int id = events[i].id;
    int x = events[i].x;

    if (events[i].type == LEFT) {
      IntervalTree::Result r = itree.get(value-1);
      // printf("LEFT: r.openedMax=%d, r.closedMax=%d\n",
      //     r.openedMax, r.closedMax);

      int best = r.openedMax + x;      
      if (best > r.closedMax) {
	dpRecon[id] = r.openedMaxLastId;
      } else {
	best = r.closedMax;
	dpRecon[id] = r.closedMaxLastId;
      }
      itree.insertOpened(value, id, x, best);
    } else { // events[i].type == RIGHT
      itree.convertOpenedToClosed(id, x);
      IntervalTree::Result r = itree.get(value);
      //printf("RIGHT: bestClosed=%d\n", r.closedMax);
      
      if (r.closedMax > *result) {
	*result = r.closedMax;
	endIndex = r.closedMaxLastId;
      }
    }
  }

  if (resultCoverage != NULL) {
    reconstruct(resultCoverage, dpRecon, endIndex, intervals);
  }
}

void coverSlow(vector<Interval>* resultCoverage,
	       int* result,
	       vector<Interval> intervals) {
  sort(intervals.begin(), intervals.end(), CoverSlowCmp());

  int n = intervals.size();
  vector<int> dp(n);
  vector<int> dpRecon(n, -1);

  *result = 0;
  int endIndex = -1;

  for (int i = 0; i < n; ++i) {
    const Interval& curr = intervals[i];
    dp[i] = curr.right - curr.left + 1;

    for (int j = i-1; j >= 0; --j) {
      const Interval& prev = intervals[j];

      if (prev.value < curr.value) {
	int extend = 0;
	if (prev.right < curr.left) {
	  extend = dp[j] + curr.right - curr.left + 1;
	} else if (curr.right >= prev.right) {
	  extend = dp[j] + curr.right - prev.right;
	}
	if (extend > dp[i]) {
	  dp[i] = extend;
	  dpRecon[i] = j;
	}
      }
    }

    if (dp[i] > *result) {
      *result = dp[i];
      endIndex = i;
    }
  }

  // slijedi rekonstrukcija
  // (*result sadrzi u ovom trenutku najbolje rjesenje)
  if (resultCoverage != NULL) {
    reconstruct(resultCoverage, dpRecon, endIndex, intervals);
  }
}

void cover(vector<Interval>* resultCoverage,
	   int* result,
	   const vector<Interval>& intervals) {
  // za razlicite brojeve intervala n razliciti algoritmi rade bolje:
  // coverSlow radi u O(n^2)
  // coverFast radi u O(n lg n), ali ima veliku konstantu koja blijedi
  //           tek kada n dovoljno naraste
  if (intervals.size() < 700) { // empirijski
    coverSlow(resultCoverage, result, intervals);
  } else {
    coverFast(resultCoverage, result, intervals);
  }
}
