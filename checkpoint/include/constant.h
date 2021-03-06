#ifndef CONSTANT_H
#define CONSTANT_H

#include <vector>
#include "memory.h"
#include "refcount.h"

using TConstant=int*;

// TConstant argument must call Malloc function
class Constant {
public:
    Constant(const TConstant& constant, int size);
    ~Constant();
    const TConstant& GetConstant() const;
    TConstant& GetConstant();
    size_t GetSize() const;
private:
    struct InnerData {
        //v_cycle iter iteration time_step
        TConstant constant;
        size_t size;
        InnerData(const TConstant& constant, int size);
        ~InnerData();
    };
    InnerData constant_;
};

// sequence of deconstruct
// RCPtr --> Constant -> RCObject
class RCConstant {
public:
    RCConstant(const TConstant& constant, int size);
    const TConstant& GetConstant() const;
    TConstant& GetConstant();
    size_t GetSize() const;
    Constant& operator*();
private:
    RCPtr<Constant> value_;
};

#endif
