#include "field.h"
#include <iostream>

Field::Field(const TField& field, const std::vector<int>& info) : field_(field, info){
}

Field::~Field() {
}

// need Copy?
// solution: make copy task to caller
const TField& Field::GetField() const {
    return field_.field;
}

const std::vector<int>& Field::GetInfo() const {
    return field_.info;
}

Field::InnerData::InnerData(const TField& t, const std::vector<int>& i) : field(t), info(i){
}

Field::InnerData::~InnerData() {
    if (field != nullptr) {
        Free(field);
        std::cout << "Free(field)" << std::endl;
    }
    field = nullptr;
}

RCField::RCField(const TField& field, const std::vector<int>& info) : value_(new Field(field, info)) {
}

const TField& RCField::GetField() const {
    return value_->GetField();
}

const std::vector<int>& RCField::GetInfo() const {
    return value_->GetInfo();
}
