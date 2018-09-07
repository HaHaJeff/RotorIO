#include "ckconfigure.h"
#include <iostream>

int main()
{
    CKConfigure config("checkpoint.xml");
    std::cout << config.GetFrequent() << std::endl;
    std::cout << config.GetGroupSize() << std::endl;
    int* num = config.GetBlockSize();
    std::cout << num[0] << "\t" << num[1] << "\t" << num[2] << std::endl;
    free(num);
}

