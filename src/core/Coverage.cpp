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
      int closedMax;
      
      Result(int openedMax = 0, int closedMax = 0) :
	openedMax(openedMax), closedMax(closedMax) {}
    };
    
    void insertOpened(int value, int id, int x, int best) {
      NodeItem ni(best-x, x, value);
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
      convertOpenedToClosed(0, 0, offset_-1, toRemove, toInsert);
    }
    
  private:

    struct NodeItem {
      int best;
      int x;
      int value;

      NodeItem(int best = 0, int x = 0, int value = 0) : 
	best(best), x(x), value(value) {}
      
      bool operator < (const NodeItem& a) const {
	if (best != a.best) { return best < a.best; }
	if (x != a.x) { return x < a.x; }
	return value < a.value;
      }
    };
    
    struct Node {
      multiset<NodeItem> opened; // best = najbolje do sad - lijevi kraj
      multiset<NodeItem> closed; // best = najbolje do sad
      
      int openedMax;
      int closedMax;

      Node() {
	openedMax = -kInfinity;
	closedMax = 0;
      }

      Result getResult() {
	return Result(openedMax, closedMax);
      }
    };

    Result get(int node, int lo, int hi, int value) {
      if (value >= hi) return ivec_[node].getResult();
      if (lo > value) return Result(-kInfinity, 0);
      
      int mid = (lo+hi)/2;
      Result l = get(node*2+1, lo, mid, value);
      Result r = get(node*2+2, mid+1, hi, value);
      return Result(max(l.openedMax, r.openedMax),
		    max(l.closedMax, r.closedMax));
    }

    void insertOpened(int node, int lo, int hi, const NodeItem& ni) {
      if (lo > ni.value) { return; }
      if (hi < ni.value) { return; }

      if (lo == hi) {
	ivec_[node].opened.insert(ni);
	ivec_[node].openedMax = ivec_[node].opened.rbegin()->best;

	if (!ivec_[node].closed.empty()) {
	  ivec_[node].closedMax = ivec_[node].closed.rbegin()->best;
	} else {
	  ivec_[node].closedMax = 0;
	}
      } else {
	int mid = (lo+hi)/2;
	insertOpened(node*2+1, lo, mid, ni);
	insertOpened(node*2+2, mid+1, hi, ni);

	ivec_[node].openedMax = max(ivec_[node*2+1].openedMax,
				    ivec_[node*2+2].openedMax);
	ivec_[node].closedMax = max(ivec_[node*2+1].closedMax,
				    ivec_[node*2+2].closedMax);
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
	} else {
	  ivec_[node].openedMax = -kInfinity;
	}
	ivec_[node].closedMax = ivec_[node].closed.rbegin()->best;
      } else {
	int mid = (lo+hi)/2;
	convertOpenedToClosed(node*2+1, lo, mid, toRemove, toInsert);
	convertOpenedToClosed(node*2+2, mid+1, hi, toRemove, toInsert);

	ivec_[node].openedMax = max(ivec_[node*2+1].openedMax,
				    ivec_[node*2+2].openedMax);
	ivec_[node].closedMax = max(ivec_[node*2+1].closedMax,
				    ivec_[node*2+2].closedMax);
      }
    }

  private:
    map<int, NodeItem> idToNodeItem_;
    vector<Node> ivec_;
    int offset_;
  };

};

void cover(vector<Interval>* resultCoverage,
	   int* result,
	   const vector<Interval>& intervals) {
  vector<Event> events;
  int maxValue = 0;

  for (size_t i = 0; i < intervals.size(); ++i) {
    events.push_back(Event(intervals[i].left, intervals[i].value, LEFT, i));
    events.push_back(Event(intervals[i].right, intervals[i].value, RIGHT, i));
    maxValue = max(maxValue, intervals[i].value);
  }

  sort(events.begin(), events.end());
  IntervalTree itree(maxValue+1);
  *result = 0;

  for (size_t i = 0; i < events.size(); ++i) {
    int value = events[i].value;
    int id = events[i].id;
    int x = events[i].x;

    if (events[i].type == LEFT) {
      IntervalTree::Result r = itree.get(value-1);
      //printf("LEFT: r.openedMax=%d, r.closedMax=%d\n",
      //     r.openedMax, r.closedMax);
      int best = r.openedMax + x;
      best = max(best, r.closedMax);
      itree.insertOpened(value, id, x, best);
    } else { // events[i].type == RIGHT
      itree.convertOpenedToClosed(id, x);
      int bestClosed = itree.get(value).closedMax;
      //printf("RIGHT: bestClosed=%d\n", bestClosed);
      *result = max(*result, bestClosed);
    }
  }
}

struct CoverSlowCmp {
  bool operator () (const Interval& a, const Interval& b) const {
    return a.left < b.left;
  }
};

void coverSlow(vector<Interval>* resultCoverage,
	       int* result,
	       const vector<Interval>& intervals) {
  vector<Interval> sortedIntervals = intervals;
  sort(sortedIntervals.begin(), sortedIntervals.end(), CoverSlowCmp());

  int n = sortedIntervals.size();
  vector<int> dp(n);
  *result = 0;

  for (int i = 0; i < n; ++i) {
    const Interval& curr = sortedIntervals[i];
    dp[i] = curr.right - curr.left + 1;

    for (int j = i-1; j >= 0; --j) {
      const Interval& prev = sortedIntervals[j];

      if (prev.value < curr.value) {
	int extend = 0;
	if (prev.right < curr.left) {
	  extend = dp[j] + curr.right - curr.left + 1;
	} else if (curr.right >= prev.right) {
	  extend = dp[j] + curr.right - prev.right;
	}
	dp[i] = max(dp[i], extend);
      }
    }
    *result = max(*result, dp[i]);
  }
}
