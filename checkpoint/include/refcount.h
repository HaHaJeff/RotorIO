#ifndef REFCOUNT_H
#define REFCOUNT_H
#include <iostream>

class RCObject {
public:
    void AddReference();
    void RemoveReference();

    void MarkUnshareable();
    bool IsShareable() const;
    bool IsShared() const;


// forbid object generate from stack
// object must generate from heap, so make construct or deconstruct be proctect
// can not make them be private, because it cause RCObject base class can not be
// inherited.
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
    const T* operator->() const;
    T* operator->();
    const T& operator*() const;
    T& operator*();

private:
    // Add an intermediate layer
    struct CountHolder : public RCObject {
        ~CountHolder() {
           // std::cout << "~CountHolder" << std::endl;
            delete pointee;
        }
        T* pointee;
    };
    CountHolder *counter_;
    void MakeCopy_();
    void Init_();
};

#include "refcount_impl.h"

#endif
