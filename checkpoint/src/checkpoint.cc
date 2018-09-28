#include "checkpoint.h"
#include "IOStrategy/Data.h"


Checkpoint::Checkpoint(Constant& constant, Field& field, int group) : constant_(constant.GetConstant(), constant.GetSize()),field_(field.GetField(), field.GetInfo()), group_(group) {}

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

    int block_id = field_.GetInfo()[4];
    int count = constant_.GetSize();
    double* cst = new double[count];
    TConstant constant = constant_.GetConstant();

    for (int i = 0; i < count; i++) {
        cst[i] = static_cast<double>(constant[i]);
    }

    Data_3D data(cst, count, block_id, "");
    io.Lseek(block_id*count*sizeof(double));
    io.Write(data, true); 
    delete[] cst;
    cst = nullptr;
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
    Data_3D data(field, num*x* y* z, block_id, "");
    io.Lseek(block_id*count*sizeof(double));
    io.Write(data, true);
}

void SetCheckpoint(Constant& constant, Field& field) {
    Checkpoint ck(constant, field);
    POSIXIO io;
    Strategy* strategy = io.GetIOStrategy(static_cast<TYPE>(0));
    char ckfld[12];
    char ckcst[15];
    int block_id = field.GetInfo()[4];

    char tmpfld[15];
    char tmpcst[18];

    snprintf(tmpfld, sizeof(tmpfld), "tmp_field_%02d", block_id);
    snprintf(tmpcst, sizeof(tmpcst), "tmp_constant_%02d", block_id);

    // save field information
    snprintf(ckfld, sizeof(ckfld), "ck_field_%02d", block_id);
    strategy->Open(tmpfld);
    ck.SaveField(*strategy);
    strategy->Close();

    // save constant information
    snprintf(ckcst, sizeof(ckcst), "ck_constant_%02d", block_id);
    strategy->Open(tmpcst);
    ck.SaveConstant(*strategy);
    strategy->Close();

    delete strategy;

    rename(tmpfld, ckfld);
    rename(tmpcst, ckcst);
}

void Restart(Constant& constant, Field& field) {

}
