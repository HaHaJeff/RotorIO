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
    std::vector<int> info = field_.GetInfo();
    int num = info[0];
    int x = info[1];
    int y = info[2];
    int z = info[3];
    int block_id = info[4];
    int count = num * x * y * z;
    double* field = field_.GetField();
    Data_3D data(field, num*x*y*z, block_id, "");
    io.Lseek(block_id*count*sizeof(double));
    io.Read(data);
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
    Data_3D data(static_cast<void*>(field), num*x* y* z, block_id, "");
    io.Lseek(block_id*count*sizeof(double));
    io.Write(data, false);
}

void SetCheckpoint(Constant& constant, Field& field) {
    Checkpoint ck(constant, field);
    POSIXIO io;
    Strategy* strategy = io.GetIOStrategy(static_cast<TYPE>(0));
    char ckfld[12];
    char ckcst[15];
    int block_id = field.GetInfo()[4];

    // save field information
    snprintf(ckfld, sizeof(ckfld), "ck_field_%02d", block_id);
    strategy->Open(ckfld);
    ck.SaveField(*strategy);
    strategy->Close(); 

    // save constant information
    snprintf(ckcst, sizeof(ckcst), "ck_constant_%02d", block_id);
    strategy->Open(ckcst);
    ck.SaveConstant(*strategy);
    strategy->Close();

    delete strategy;
}

void Restart(Constant& constant, Field& field) {

}