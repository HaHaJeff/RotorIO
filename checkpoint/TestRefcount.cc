#include "field.h"
#include "constant.h"

int main()
{
  TField value;
  Malloc(value, 1, 2, 3, 4);
  RCField field(value);

  TConstant value_1;
  Malloc(value_1, 1);
  RCConstant constant(value_1);

}
