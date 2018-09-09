#include "field.h"
#include <iostream>
Field::Field(const TField& field) : field_(new FieldValue(field)){
}

Field::FieldValue::FieldValue(const TField& value) : data(value) {
}

Field::FieldValue::FieldValue(const FieldValue& rhs) {
    data = rhs.data;
}

Field::FieldValue::~FieldValue() {
    Free(data);
    std::cout << "~FieldValue" << std::endl;
}
