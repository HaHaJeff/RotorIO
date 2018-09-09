#ifndef FIELD_H
#define FIELD_H

#include <vector>
#include "IOStrategy/HDF5IO.h"
#include "memory.h"

class Field {
public:
    using TField=double****;
public:
    Field(const TField& field);
    Field(const Field& field);
    Field();
    Field& operator=(const Field& field);
    ~Field();
    
private:
    //q11 q12 q13 q14 q15 q16
    TField field_;
};

#endif