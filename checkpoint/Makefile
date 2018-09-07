CXX=g++
CXXFLAGS+=-g -w -std=c++11 -fopenmp 

LIBRARY=-lpthread -L../io -lio -L../tinyxml -ltinyxml ../io/third_party/lib/libhdf5.a ../io/third_party/lib/libhdf5_hl.a -lz -ldl -lm

INC_DIR=-I./include -I./src -I../io/include -I../tinyxml/include -I../io/third_party/include

all: TestCKConfigure

TestCKConfigure : ckconfigure.cc TestCKConfigure.cc
	@echo "Compiling $^ ==> $@..."
	$(CXX) $^ -o $@ $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 
