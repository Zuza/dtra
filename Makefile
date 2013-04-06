CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
LD_FLAGS := -pthread
CC_FLAGS := -O2 --std=c++0x -Wno-unused-result -ggdb -g

client: $(OBJ_FILES)
	g++ $(LD_FLAGS) -o $@ $^

obj/%.o: src/%.cpp src/%.h
	mkdir -p obj
	g++ $(CC_FLAGS) -c -o $@ $<
