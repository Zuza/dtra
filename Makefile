export INCLUDES := -I$(CURDIR)/src/
export MAIN_DIR := $(CURDIR)/src/main/
export TEST_DIR := $(CURDIR)/src/test/
export CORE_DIR := $(CURDIR)/src/core/
export SSW_DIR := $(CURDIR)/src/ssw/
export FMINDEX_DIR := $(CURDIR)/src/FmIndexWavelet/

CORE_H_FILES := $(wildcard $(CORE_DIR)/*.h)
CORE_CPP_FILES := $(wildcard $(CORE_DIR)/*.cpp)
CORE_OBJ_FILES := $(addprefix obj/core/,$(notdir $(CORE_CPP_FILES:.cpp=.o)))
SSW_H_FILES := $(wildcard $(SSW_DIR)/*.h)
SSW_CPP_FILES := $(wildcard $(SSW_DIR)/*.cpp)
SSW_OBJ_FILES := $(addprefix obj/ssw/,$(notdir $(SSW_CPP_FILES:.cpp=.o)))
MAIN_H_FILES := $(wildcard $(MAIN_DIR)/*.h)
MAIN_CPP_FILES := $(wildcard $(MAIN_DIR)/*.cpp)
MAIN_OBJ_FILES := $(addprefix obj/main/,$(notdir $(MAIN_CPP_FILES:.cpp=.o)))
TEST_H_FILES := $(wildcard $(TEST_DIR)/*.h)
TEST_CPP_FILES := $(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJ_FILES := $(addprefix obj/test/,$(notdir $(TEST_CPP_FILES:.cpp=.o)))
FMINDEX_H_FILES := $(wildcard $(FMINDEX_DIR)/*.hpp)
FMINDEX_CPP_FILES := $(filter-out $(wildcard $(FMINDEX_DIR)/main_*.cpp),$(wildcard $(FMINDEX_DIR)/*.cpp))  # remove main files
FMINDEX_OBJ_FILES := $(addprefix obj/FmIndexWavelet/,$(notdir $(FMINDEX_CPP_FILES:.cpp=.o)))

CC := mpiCC.openmpi
LD_FLAGS := -pthread -lgflags -ldivsufsort
CC_FLAGS := -fopenmp -O2 --std=c++0x -Wno-unused-result -D_FILE_OFFSET_BITS=64 $(INCLUDES)

all: client reducer lisa simulator test

forceall: clean all

client: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES) $(FMINDEX_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/FmIndexWavelet/*.o obj/main/client.o $(LD_FLAGS)

reducer: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/main/DataReducer.o $(LD_FLAGS)

lisa: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/main/lisa.o $(LD_FLAGS)

simulator: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/main/rand_lcs_distr.o $(LD_FLAGS)

test: test_coverage

test_coverage: $(TEST_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES)
	mkdir -p obj/test
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/test/testCoverage.o $(LD_FLAGS)

$(CORE_OBJ_FILES): $(CORE_CPP_FILES) $(CORE_H_FILES)
	mkdir -p obj/core
	$(CC) $(CC_FLAGS) -c -o $@ $(CORE_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(SSW_OBJ_FILES): $(SSW_CPP_FILES) $(SSW_H_FILES)
	mkdir -p obj/ssw
	$(CC) $(CC_FLAGS) -c -o $@ $(SSW_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(MAIN_OBJ_FILES): $(MAIN_CPP_FILES) $(MAIN_H_FILES)
	mkdir -p obj/main
	$(CC) $(CC_FLAGS) -c -o $@ $(MAIN_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(TEST_OBJ_FILES): $(TEST_CPP_FILES) $(TEST_H_FILES)
	mkdir -p obj/test
	$(CC) $(CC_FLAGS) -c -o $@ $(TEST_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(FMINDEX_OBJ_FILES): $(FMINDEX_CPP_FILES) $(FMINDEX_H_FILES)
	mkdir -p obj/FmIndexWavelet
	$(CC) $(CC_FLAGS) -c -o $@ $(FMINDEX_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

clean:
	rm -rf obj
	rm -f bin/client
	rm -f bin/reducer
