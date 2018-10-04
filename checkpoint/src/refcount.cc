#include "refcount.h"
#include <iostream>

void RCObject::AddReference() {
    ++refcount_;
}

void RCObject::RemoveReference() {
    if (--refcount_ == 0) {
    //    std::cout << "delete this" << std::endl;
        delete this;
    }
}

void RCObject::MarkUnshareable() {
    shareable_ = false;
}

bool RCObject::IsShareable() const {
    return shareable_;
}

bool RCObject::IsShared() const {
    return refcount_ > 0;
}

RCObject::RCObject() : refcount_(0), shareable_(true) {
}

RCObject::RCObject(const RCObject& rhs) : refcount_(0), shareable_(true) {
}

RCObject& RCObject::operator=(const RCObject& rhs) {
    return *this;
}

RCObject::~RCObject() {
 // std::cout << "~RCObjecet" << std::endl;
}

