#include "ckconfigure.h"

CKConfigure::CKConfigure(const char* filename) : config_(new TiXmlDocument()) {
    config_->LoadFile(filename);
}

CKConfigure::~CKConfigure() {
    if (config_ != NULL) delete config_; 

    config_ = NULL;
}


bool CKConfigure::Init() {

}

void CKConfigure::SetFrequent(int frequent) {

}

int CKConfigure::GetFrequent() {
    return 0;
}

void CKConfigure::SetGroupSize(int num) {

}

int CKConfigure::GetGroupSize() {
    return 0;
}

void CKConfigure::SetBlockSize() {

}

int CKConfigure::GetBlockSize() {
    return 0;
}