#include "field.h"
#include "constant.h"

int main()
{
  TField value;
  Malloc(value, 1);
  std::vector<int> info = {1, 2, 3, 4};
  RCField field(value, info);

  TConstant value_1;
  Malloc(value_1, 1);
  RCConstant constant(value_1, 1);

}
