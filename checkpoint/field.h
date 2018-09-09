#ifndef FIELD_H
#define FIELD_H

#include "IOStrategy/HDF5IO.h"
#include "memory.h"
#include "refcount.h"

class Field {
public:
    using TField=double****;
public:
// field need generate from heap
    Field(const TField& field);
private:
    //q11 q12 q13 q14 q15 q16
    struct FieldValue : public RCObject {
        TField data;
        FieldValue(const TField& value);
        FieldValue(const FieldValue& rhs);
        ~FieldValue();
    };
    RCPtr<FieldValue> field_;
};

#endif
