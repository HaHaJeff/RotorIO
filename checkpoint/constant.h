#ifndef CONSTANT_H
#define CONSTANT_H

#include <vector>
#include "memory.h"
#include "refcount.h"

using TConstant=int*;

// TConstant argument must call Malloc function
class Constant {
public:
    Constant(const TConstant& constant);
    ~Constant();
    const TConstant& GetConstant() const;
private:
    //v_cycle iter iteration time_step
    TConstant constant_;
};

class RCConstant {
public:
    RCConstant(const TConstant& constant);
    const TConstant& GetConstant() const;
private:
    RCPtr<Constant> value_;
};

#endif
