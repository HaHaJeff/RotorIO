#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "constant.h"
#include "field.h"
#include "refcount.h"
#include "IOStrategy/POSIXIO.h"

// TODO: 把数据直接交给checkpoint对象，由checkpoint完成数据的恢复与存储
// 采用double***对三维数据进行存储，那么是不是需要考虑数据的拷贝，是否可以考虑引入引用计数技术
class Checkpoint {
public:
  //Constant: int* size, Field: double*, vector<int>(x,y,z,k,block_id)
  Checkpoint(Constant& constant, Field& field);
  void RestoreConstant(Strategy& io);
  void RestoreField(Strategy& io);
  void SaveConstant(Strategy& io);
  void SaveField(Strategy& io);
  const Constant& GetConstant();
  const Field& GetField();

private:
  RCConstant constant_;
  RCField    field_;
};

void SetCheckpoint(Constant& constant, Field& field);
void Restart(Constant& constant, Field& field);

#endif
