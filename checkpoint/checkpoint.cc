#include "checkpoint.h"
#include "IOStrategy/Data.h"


Checkpoint::Checkpoint(Constant& constant, Field& field) : constant_(constant.GetConstant(), constant.GetSize()),field_(field.GetField(), field.GetInfo()) {}
 
const Constant& Checkpoint::GetConstant() {
    return *constant_;
}

const Field& Checkpoint::GetField() {
    return *field_;
}

// before write or read, need set file view
void Checkpoint::RestoreConstant(Strategy& io) {

}

void Checkpoint::RestoreField(Strategy& io) {

}

void Checkpoint::SaveConstant(Strategy& io) {

}

void Checkpoint::SaveField(Strategy& io) {
    std::vector<int> info = field_.GetInfo();
    int num = info[0];
    int x = info[1];
    int y = info[2];
    int z = info[3];
    int block_id = info[4];
    int count = num * x * y * z;
    double* field = field_.GetField();
    Data_3D data(field, num*x, y, z, "1", block_id);
    io.Lseek(block_id*count);
    io.Write(data, true);
}