CXX=mpicxx
CXXFLAGS+=-g -w -std=c++11 -fopenmp 

LIBRARY=-lpthread -L../io -lio -L../tinyxml -ltinyxml ../io/third_party/lib/libhdf5.a ../io/third_party/lib/libhdf5_hl.a -lz -ldl -lm

INC_DIR=-I./include -I../io/include -I../tinyxml/include -I../io/third_party/include

TEST_DIR=./Test
TEST_OBJ_DIR=${TEST_DIR}/obj

OBJ_DIR=./obj
SRC_DIR=./src
BIN_DIR=./bin

TARGET=$(basename ${TARGET_SRCS})
TARGET_SRCS=$(wildcard ${TEST_DIR}/*.cc)
TARGET_OBJS=$(patsubst %.cc, %.o, ${TARGET_SRCS})

SRCS=$(foreach dir, $(SRC_DIR), $(wildcard ${dir}/*.cc))
OBJFILES=$(patsubst %.cc, %.o, $(SRCS))
OBJS=$(patsubst %.cc, $(OBJ_DIR)/%.o, $(notdir $(SRCS)))

LIB=libck.a

all: init $(OBJFILES) $(LIB) $(TARGET_OBJS) TestCheckpoint TestCKConfigure TestRefcount 

$(OBJFILES): %.o:%.cc
	@echo "Compilinh $^ ==> $@..."
	$(CXX) -c $^ -o $(OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS)

$(TARGET_OBJS):%.o:%.cc
	@echo "Compiling $^ ==> $@..."
	$(CXX) -c $^ -o $(TEST_OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS)

$(LIB) : $(OBJS)
	rm -f $@
	ar cr $@ $(^)

TestCKConfigure : $(OBJS) 
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestCKConfigure.o $^ -o $(BIN_DIR)/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

TestRefcount: $(OBJS) 
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestRefcount.o $^ -o $(BIN_DIR)/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

TestCheckpoint : $(OBJS)
	@echo "Compiling $^ ==> $@..."
	$(CXX) $(TEST_OBJ_DIR)/TestCheckpoint.o $^ -o $(BIN_DIR)/$@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

init:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(TEST_OBJ_DIR)
	mkdir -p $(BIN_DIR)

.PHONY : clean
clean :
	-rm -r $(OBJ_DIR) $(BIN_DIR) $(TEST_OBJ_DIR) $(LIB)
