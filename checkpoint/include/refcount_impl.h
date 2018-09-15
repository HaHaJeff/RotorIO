template<class T>
void RCPtr<T>::Init_() {
    if (counter_->IsShareable() == false) {
        T* oldvalue = counter_->pointee;
        counter_->pointee = new T(*oldvalue);
    }

    counter_->AddReference();
}

template<class T>
RCPtr<T>::RCPtr(T* realptr) : counter_(new CountHolder){
    counter_->pointee = realptr;
    Init_();
}

template<class T>
RCPtr<T>::RCPtr(const RCPtr& rhs) : counter_(rhs.counter_){
    Init_();
}

template<class T>
RCPtr<T>::~RCPtr() {
    counter_->RemoveReference();
}

template<class T>
RCPtr<T>& RCPtr<T>::operator=(const RCPtr& rhs) {
    if (counter_!= rhs.counter_) {
        if (counter_!= nullptr) {
            counter_->RemoveReference();
            counter_ = rhs.counter_;
            Init_();
        }
    }
    return *this;
}

template<class T>
T* RCPtr<T>::operator->(){
    return counter_->pointee;
}

template<class T>
const T* RCPtr<T>::operator->() const {
    return counter_->pointee;
}

template<class T>
T& RCPtr<T>::operator*(){
    return *(counter_->pointee);
}

template<class T>
const T& RCPtr<T>::operator*() const{
    return *(counter_->pointee);
}
