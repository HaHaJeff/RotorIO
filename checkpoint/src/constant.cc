#include "constant.h"
#include <iostream>

Constant::Constant(const TConstant& constant, int size) : constant_(constant, size){}

Constant::~Constant() {
}

const TConstant& Constant::GetConstant() const {
    return constant_.constant;
}

TConstant& Constant::GetConstant() {
    return constant_.constant;
}

size_t Constant::GetSize() const {
    return constant_.size;
}

Constant::InnerData::InnerData(const TConstant& constant, int size) : constant(constant), size(size) {
}

Constant::InnerData::~InnerData() {
    if (constant != nullptr) {
       // Free(constant);
        std::cout << "Free(constant)" << std::endl;
    }
    constant = nullptr;
}

RCConstant::RCConstant(const TConstant& constant, int size) : value_(new Constant(constant, size)) {
}

const TConstant& RCConstant::GetConstant() const {
    return value_->GetConstant();
}

TConstant& RCConstant::GetConstant(){
    return value_->GetConstant();
}

Constant& RCConstant::operator*() {
    return *value_;
}

size_t RCConstant::GetSize() const {
    return value_->GetSize();
}
