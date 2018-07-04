CXX=mpic++
CXXFLAGS+=-g -Wall -std=c++11
LIBRARY=-lpthread -L../io -lio

INC_DIR=-I./include -I./src -I../io/include

TEST_DIR=./Test
TEST_OBJ_DIR=${TEST_DIR}/obj

OBJ_DIR=./obj
SRC_DIR=./src
BIN_DIR=./bin

TARGET=$(basename ${TARGET_SRCS})
TARGET_SRCS=$(wildcard ${TEST_DIR}/*.cpp)
TARGET_OBJS=$(patsubst %.cpp, %.o, ${TARGET_SRCS})

SRCS=$(foreach dir, $(SRC_DIR), $(wildcard ${dir}/*.cpp))
OBJFILES=$(patsubst %.cpp, %.o, $(SRCS))
OBJS=$(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRCS)))

all: $(OBJFILES) $(TARGET_OBJS) Rotor-OOP-V1.0 

$(OBJFILES):%.o:%.cpp
	@echo "Compiling $^ ==> $@..."
	$(CXX) -c $^ -o $(OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

$(TARGET_OBJS):%.o:%.cpp
	@echo "Compiling $^ ==> $@..."
	$(CXX) -c $^ -o $(TEST_OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

Rotor-OOP-V1.0: $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/Rotor-OOP-V1.0.o $^ -o ${BIN_DIR}/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(TEST_OBJ_DIR):
	mkdir -p $(TEST_OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)
.PHONY : clean
clean :
	-rm $(OBJ_DIR)/* $(BIN_DIR)/* $(TEST_OBJ_DIR)/*

