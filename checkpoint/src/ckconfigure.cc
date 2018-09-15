#include "ckconfigure.h"
#include <cstdlib>
#include <stdio.h>
#include <iostream>

CKConfigure::CKConfigure(const char* filename) : config_(new TiXmlDocument()) {
  if (!config_->LoadFile(filename)) {
    std::cerr << config_->ErrorDesc() << std::endl;
    exit(-1);
  }
  if ((root_ = config_->RootElement()) == NULL) {
    std::cerr << "Failed to load file: No root element." << std::endl;
    exit(-1);
  }
  properties_ = root_->FirstChildElement();
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
    return atoi(properties_->Attribute("frequent"));
}

void CKConfigure::SetGroupSize(int num) {
}

int CKConfigure::GetGroupSize() {
    return atoi(properties_->Attribute("group_size"));
}

void CKConfigure::SetBlockSize(int dims[3]) {
}

int* CKConfigure::GetBlockSize() {
    int x = atoi(properties_->Attribute("block_x"));
    int y = atoi(properties_->Attribute("block_y"));
    int z = atoi(properties_->Attribute("block_z"));
    int* ret = (int*)malloc(3*sizeof(int));
    ret[0] = x;
    ret[1] = y;
    ret[2] = z;
    return ret;
}
