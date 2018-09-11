#ifndef FIELD_H
#define FIELD_H

#include "IOStrategy/HDF5IO.h"
#include "memory.h"
#include "refcount.h"

#include <vector>

using TField=double*;

class Field {
public:
    // field need generate from heap
    explicit Field(const TField& field, const std::vector<int>& info);
    ~Field();
    const TField& GetField() const;
    const std::vector<int>& GetInfo() const;

 private:
    struct InnerData {
        //q11 q12 q13 q14 q15 q16
        TField field;
        //erver dim infomation and block_id
        std::vector<int> info;
        InnerData(const TField& t, const std::vector<int>& info); 
        ~InnerData();
    };
    InnerData field_;
};

class RCField {
public:
    RCField(const TField& field, const std::vector<int>& info);
    const TField& GetField() const;
    const std::vector<int>& GetInfo() const;
private:
    RCPtr<Field> value_;
};

#endif
