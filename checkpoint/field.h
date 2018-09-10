#ifndef FIELD_H
#define FIELD_H

#include "IOStrategy/HDF5IO.h"
#include "memory.h"
#include "refcount.h"

using TField=double****;

class Field {
public:
    // field need generate from heap
    explicit Field(const TField& field);
    ~Field();
    const TField& GetField() const;

 private:
    //q11 q12 q13 q14 q15 q16
    TField field_;
};

class RCField {
public:
    RCField(const TField& field);
    const TField& GetField() const;
private:
    RCPtr<Field> value_;
};

#endif
