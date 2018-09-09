#ifndef CONSTANT_H
#define CONSTANT_H

#include <vector>
#include "memory.h"

// TConstant argument must call Malloc function
class Constant {
public:
    using TConstant=int*;
public:
    Constant(const TConstant& constant) : constant_(constant) {}
    Constant();
    ~Constant();
private:
    //v_cycle iter iteration time_step
    TConstant constant_;
};

#endif
