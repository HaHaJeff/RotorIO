#ifndef CKCONFIGURE_H
#define CKCONFIGURE_H

#include "tinyxml.h"

class CKConfigure {
public:
    // tinyxml loadfile
    CKConfigure(const char* filename);
    
    // generate empty xml file(filename is "CKConfigure.xml")
    // frequent, group_size, root_dir
    bool Init();

    void SetFrequent(int frequent);
    int GetFrequent();

    void SetGroupSize(int num);
    int GetGroupSize();

    void SetRootDir(const char* dirname);
    const char* GetRootDir();

    void SetBlockSize(int dims[3]);
    int* GetBlockSize();

    ~CKConfigure();

private:
    TiXmlDocument* config_;
};

#endif