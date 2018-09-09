#include "field.h"

Field::Field(const TField& field) : field_(field) {}

Field::Field(const Field& field) {}

Field::Filed() : field_(NULL) {}

Field&::operator(const Field& field) {
    if (this->field_ != field.field_) {
        
    }
}

Field::~Field() {
    if (NULL != field_) {
        Free(field_);
    }
    field_ = NULL;
}