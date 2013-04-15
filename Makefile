export INCLUDES := -I$(CURDIR)/src/
export MAIN_DIR := $(CURDIR)/src/main/
export CORE_DIR := $(CURDIR)/src/core/

CORE_CPP_FILES := $(wildcard $(CORE_DIR)/*.cpp)
CORE_OBJ_FILES := $(addprefix obj/core/,$(notdir $(CORE_CPP_FILES:.cpp=.o)))
MAIN_CPP_FILES := $(wildcard $(MAIN_DIR)/*.cpp)
MAIN_OBJ_FILES := $(addprefix obj/main/,$(notdir $(MAIN_CPP_FILES:.cpp=.o)))

LD_FLAGS := -pthread
CC_FLAGS := -O2 --std=c++0x -Wno-unused-result $(INCLUDES) -g -ggdb

all: client reducer

client: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES)
	g++ $(LD_FLAGS) -o $@ obj/core/*.o obj/main/client.o

reducer: $(MAIN_OBJ_FILES) $(CORE_OBJ_FILES) 
	g++ $(LD_FLAGS) -o $@ obj/core/*.o obj/main/DataReducer.o

$(CORE_OBJ_FILES): $(CORE_CPP_FILES)
	mkdir -p obj/core
	g++ $(CC_FLAGS) -c -o $@ $(CORE_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(MAIN_OBJ_FILES): $(MAIN_CPP_FILES)
	mkdir -p obj/main
	g++ $(CC_FLAGS) -c -o $@ $(MAIN_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

clean:
	rm -r obj
	rm client
	rm reducer
