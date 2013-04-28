export INCLUDES := -I$(CURDIR)/src/
export MAIN_DIR := $(CURDIR)/src/main/
export CORE_DIR := $(CURDIR)/src/core/
export SSW_DIR := $(CURDIR)/src/ssw/

CORE_H_FILES := $(wildcard $(CORE_DIR)/*.h)
CORE_CPP_FILES := $(wildcard $(CORE_DIR)/*.cpp)
CORE_OBJ_FILES := $(addprefix obj/core/,$(notdir $(CORE_CPP_FILES:.cpp=.o)))

SSW_H_FILES := $(wildcard $(SSW_DIR)/*.h)
SSW_CPP_FILES := $(wildcard $(SSW_DIR)/*.cpp)
SSW_OBJ_FILES := $(addprefix obj/ssw/,$(notdir $(SSW_CPP_FILES:.cpp=.o)))

MAIN_H_FILES := $(wildcard $(MAIN_DIR)/*.h)
MAIN_CPP_FILES := $(wildcard $(MAIN_DIR)/*.cpp)
MAIN_OBJ_FILES := $(addprefix obj/main/,$(notdir $(MAIN_CPP_FILES:.cpp=.o)))

CC := mpic++ # to je omotac oko g++-a koji brine o MPI pathovima
LD_FLAGS := -pthread -lgflags
CC_FLAGS := -DDEBUG -fopenmp -O2 --std=c++0x -Wno-unused-result -D_FILE_OFFSET_BITS=64 $(INCLUDES)
#PAZI, IMA DEBUG!!!

all: client reducer

forceall: clean all

client: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/main/client.o $(LD_FLAGS)

reducer: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) $(SSW_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o bin/$@ obj/core/*.o obj/ssw/*.o obj/main/DataReducer.o $(LD_FLAGS)

$(CORE_OBJ_FILES): $(CORE_CPP_FILES) $(CORE_H_FILES)
	mkdir -p obj/core
	$(CC) $(CC_FLAGS) -c -o $@ $(CORE_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(MAIN_OBJ_FILES): $(MAIN_CPP_FILES) $(MAIN_H_FILES)
	mkdir -p obj/main
	$(CC) $(CC_FLAGS) -c -o $@ $(MAIN_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(SSW_OBJ_FILES): $(SSW_CPP_FILES) $(SSW_H_FILES)
	mkdir -p obj/ssw
	$(CC) $(CC_FLAGS) -c -o $@ $(SSW_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

clean:
	rm -rf obj
	rm -f bin/client
	rm -f bin/reducer
