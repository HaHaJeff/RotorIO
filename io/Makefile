CXX=mpic++
CXXFLAGS+=-g -w -std=c++11
LIBRARY=-lpthread ./third_party/lib/libhdf5.a ./third_party/lib/libhdf5_hl.a -lz -ldl -lm

INC_DIR=-I./include -I./include/IOStrategy -I./third_party/include

TEST_DIR=./Test
TEST_OBJ_DIR=${TEST_DIR}/obj

OBJ_DIR=./obj
SRC_DIR=./src ./src/IOStrategy
BIN_DIR=./bin

TARGET=$(basename ${TARGET_SRCS})
TARGET_SRCS=$(wildcard ${TEST_DIR}/*.cc)
TARGET_OBJS=$(patsubst %.cc, %.o, ${TARGET_SRCS})

SRCS=$(foreach dir, $(SRC_DIR), $(wildcard ${dir}/*.cc))
OBJFILES=$(patsubst %.cc, %.o, $(SRCS))
OBJS=$(patsubst %.cc, $(OBJ_DIR)/%.o, $(notdir $(SRCS)))

LIB=libio.a

all: init $(OBJFILES) $(TARGET_OBJS) $(LIB)
#TestAtomic TestThread TestMPIIO TestPOSIXIO TestHDF5IO

$(OBJFILES):%.o:%.cc
	@echo "Compiling $^ ==> $@..."
	$(CXX) -c $^ -o $(OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

$(TARGET_OBJS):%.o:%.cc
	@echo "Compiling $^ ==> $@..."
	$(CXX) -c $^ -o $(TEST_OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

$(LIB) : $(OBJS) ./third_party/lib/libhdf5.a ./third_party/lib/libhdf5_hl.a
	rm -f $@
	ar cr $@ $(^)

TestAtomic : $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestAtomic.o $^ -o ${BIN_DIR}/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

TestThread : $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestThread.o $^ -o ${BIN_DIR}/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS)

TestMPIIO : $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestMPIIO.o $^ -o ${BIN_DIR}/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS)

TestPOSIXIO : $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestPOSIXIO.o $^ -o ${BIN_DIR}/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS)

TestHDF5IO : $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestHDF5IO.o $^ -o ${BIN_DIR}/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS)

init:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(TEST_OBJ_DIR)
	mkdir -p $(BIN_DIR)

.PHONY : clean
clean :
	-rm -r $(OBJ_DIR) $(BIN_DIR) $(TEST_OBJ_DIR) $(LIB)

