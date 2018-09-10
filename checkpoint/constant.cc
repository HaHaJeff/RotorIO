#include "constant.h"

Constant::Constant(const TConstant& constant) : constant_(constant){}

Constant::~Constant() {
    if (constant_ != nullptr) {
        Free(constant_);
    }
    constant_ = nullptr;
}

const TConstant& Constant::GetConstant() const {
    return constant_;
}

RCConstant::RCConstant(const TConstant& constant) : value_(new Constant(constant)) {
}

const TConstant& RCConstant::GetConstant() const {
    return value_->GetConstant();
}