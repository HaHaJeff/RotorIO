#ifndef REFCOUNT_H
#define REFCOUNT_H

class RCObject {
public:
    void AddReference();
    void RemoveReference();

    void MarkUnshareable();
    bool IsShareable() const;
    bool IsShared() const;

protected:
    RCObject();
    RCObject(const RCObject& rhs);
    RCObject& operator=(const RCObject& rhs);
    //Means RCObject should be a base class;
    virtual ~RCObject() = 0;

private:
    int refcount_;
    bool shareable_;
};

// make operatoring refcount more convenient
// T must inherit RCObject
template<class T>
class RCPtr {
public:
    RCPtr(T* realptr = nullptr);
    RCPtr(const RCPtr& rhs);
    ~RCPtr();

    RCPtr& operator=(const RCPtr& rhs);
    T* operator->() const;
    T& operator*() const;

private:
    T* pointee_;
    void Init_();
};

#include "refcount_impl.h"

#endif
