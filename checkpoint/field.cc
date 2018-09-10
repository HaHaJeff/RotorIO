#include "field.h"
#include <iostream>
Field::Field(const TField& field) : field_(field){
}

Field::~Field() {
    if (field_ != nullptr) {
        Free(field_);
        std::cout<< "Free" << std::endl;
    }
    field_ = nullptr;
}

// need Copy?
// solution: make copy task to caller
const TField& Field::GetField() const {
    return field_;
}

RCField::RCField(const TField& field) : value_(new Field(field)) {
}

const TField& RCField::GetField() const {
    return value_->GetField();
}
