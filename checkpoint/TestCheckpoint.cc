#include "checkpoint.h"
#include "memory.h"
#include <iostream>
int main()
{
    double ****ptr;
    int a = 1;
    Malloc(ptr, 3,2,2,3);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                for (int z = 0; z < 3; z++) {
                    ptr[i][j][k][z] = a++;
                    std::cout << ptr[i][j][k][z] << std::endl;
                }
            }
        }
    }
    std::vector<int> info = {3,2,2,3,0};
    Field field(&ptr[0][0][0][0],info);
    int *iPtr;
    Malloc(iPtr, 3);
    Constant constant(iPtr, 3);
    Checkpoint ck(constant, field);


    POSIXIO io;
    Strategy* strategy = io.GetIOStrategy(static_cast<TYPE>(0));
    strategy->Open("test");
    ck.SaveField(*strategy);
}
