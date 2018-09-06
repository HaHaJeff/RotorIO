#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include "IOStrategy/HDF5IO.h"
#include "memory.h"

class Field {
public:
    using TField=double****;
public:
    Field(const TField& field) : field_(field) {}
    Field() {}
    void ToFile(const char* filename);
    void ToMemory(const char* filename);
private:
    //q11 q12 q13 q14 q15 q16
    TField field_;
};

#endif