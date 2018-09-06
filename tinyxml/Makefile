CXX=mpic++
CXXFLAGS+=-g -Wall -std=c++11
LIBRARY=

INC_DIR=-I./include

SUBFIX=.cpp
OBJ_DIR=./obj
SRC_DIR=./src

SRCS=$(foreach dir, $(SRC_DIR), $(wildcard ${dir}/*.cpp))
OBJFILES=$(patsubst %.cpp, %.o, $(SRCS))
OBJS=$(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRCS)))

LIB=libtinyxml.a

all: init $(OBJFILES) $(LIB)

$(OBJFILES):%.o:%.cpp
	@echo "Compiling $^ ==> $@..."
	$(CXX) -c $^ -o $(OBJ_DIR)/$(notdir $@) $(LIBRARY) $(INC_DIR) $(CXXFLAGS) 

$(LIB) : $(OBJS) 
	rm -f $@
	ar cr $@ $(^)

init:
	mkdir -p $(OBJ_DIR)

.PHONY : clean
clean :
	-rm -rf $(OBJ_DIR) 

