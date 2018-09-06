#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "constant.h"
#include "field.h"

class Checkpoint {
public:
  Checkpoint(const Constant& constant, const Field& field) : constant_(constant), field_(field) {}
  Checkpoint() {}
  void RestoreConstant(const char* filename);
  void RestoreField(const char* filename);
  void SaveConstant(const char* filename);
  void SaveField(const char* filename);

private:
  Constant constant_;
  Field    field_;
};

#endif
