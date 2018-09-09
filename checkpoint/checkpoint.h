#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "constant.h"
#include "field.h"

// TODO: 把数据直接交给checkpoint对象，由checkpoint完成数据的恢复与存储
// 采用double***对三维数据进行存储，那么是不是需要考虑数据的拷贝，是否可以考虑引入引用计数技术
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
