#ifndef CONSTANT_H
#define CONSTANT_H

#include <vector>
#include "memory.h"

class Constant {
public:
    using TConstant=std::vector<int>;
public:
    Constant(const TConstant& constant) : constant_(constant) {}
    Constant() {}
    void ToFile(const char* filename);
    void ToMemory(const char* filename);
private:
    //v_cycle iter iteration time_step
    TConstant constant_;
};

#endif
