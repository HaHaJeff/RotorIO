template<class T>
void RCPtr<T>::Init_() {
    if (pointee_ == nullptr) return;
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
    if (pointee_ != nullptr) {
        pointee_->RemoveReference();
    }
}

template<class T>
RCPtr<T>& RCPtr<T>::operator=(const RCPtr& rhs) {
    if (pointee_ != rhs.pointee_) {
        if (pointee_ != nullptr) {
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
