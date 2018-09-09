#ifndef REFCOUNT_IMPL_H
#define REFCOUNT_IMPL_H

#include "refcount.h"

void RCObject::AddReference() {
    ++refcount_;
}

void RCObject::RemoveReference() {
    if (--refcount_ == 0) {
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

RCObject::RCObject(const RCObjecet& rhs) : refcount_(0), shareable_(true) {
}

RCObject& RCObject::operator(const RCObject& rhs) {
    return *this;
}

RCObject::~RCObject() {
}


template<class T>
void RCPtr::Init_() {
    if (pointee_ == NULL) return;
    if (pointee_->IsShareable() == false) {
        pointee_ = new T(*pointee_);
    }

    pointee_->AddReference();
}

template<class T>
RCPtr<T>::RCPtr(T* realptr) : pointee_(realptr){
    Init_();
}

template<class T>
RCPtr<T>::~RCPtr() {
    if (pointee_ != NULL) {
        pointee_->RemoveReference();
    }
}

template<class T>
RCPtr<T>& RCPtr<T>::operator=(const RCPtr& rhs) {
    if (pointee_ != rhs.pointee_) {
        if (pointee_ != NULL) {
            pointee_ = rhs.pointee_;
            Init_();
        }
    }
    return *this;
}

template<class T>
T* RCPtr<T>::operator->() const {
    return pointee_;
}

template<class T>
T& RCPtr<T>::operator*() const {
    return *pointee_;
}

#endif