CFLAGS = -O2 -Wall -std=c++0x
LDFLAGS = -ldivsufsort

OBJS = RankedBitmap.o WaveletTree.o FmIndex.o DnaIndex.o util.o
MAINS = main main_simple main_construct main_read_index main_wavelet_test main_dnaindex

all: $(MAINS)

%.o: %.cpp
	g++ -c $< -o $@ $(CFLAGS) $(LDFLAGS)

main: $(OBJS) main.cpp
	g++ $+ -o $@ $(CFLAGS) $(LDFLAGS)

main_simple: $(OBJS) main_simple.cpp
	g++ $+ -o $@ $(CFLAGS) $(LDFLAGS)

main_construct: $(OBJS) main_construct.cpp
	g++ $+ -o $@ $(CFLAGS) $(LDFLAGS)

main_read_index: $(OBJS) main_read_index.cpp
	g++ $+ -o $@ $(CFLAGS) $(LDFLAGS)

main_wavelet_test: $(OBJS) main_wavelet_test.cpp
	g++ $+ -o $@ $(CFLAGS) $(LDFLAGS)

main_dnaindex: $(OBJS) main_dnaindex.cpp
	g++ $+ -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *.out
	rm -f $(MAINS)
