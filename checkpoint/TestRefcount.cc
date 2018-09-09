#include "field.h"

int main()
{
  Field::TField value;
  Malloc(value, 1, 2, 3, 4);
  Field field(value);
}
