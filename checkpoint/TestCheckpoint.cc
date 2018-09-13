#include "checkpoint.h"
#include "memory.h"

#include <iostream>
#include <string.h>
#include <thread>

void Func(Constant& constant, Field& field) {
    Checkpoint ck(constant, field);
    POSIXIO io;
    Strategy* strategy = io.GetIOStrategy(static_cast<TYPE>(0));
    strategy->Open("test");
    ck.RestoreField(*strategy);

    const Field& f= ck.GetField();
    double *p = f.GetField();
    std::vector<int> info = f.GetInfo();

    std::cout << std::this_thread::get_id() << std::endl;
    for (int i = 0; i < info[0]*info[1]*info[2]*info[3]; ++i) {
      std::cout << p[i] << std::endl;
    }
    std::cout << "block_id: " << info[4] << std::endl;
}

int main()
{
    double ****ptr;
    int a = 1;
    Malloc(ptr, 3,2,2,3);

    double ****ptr1;
    Malloc(ptr1, 3, 2, 2, 3);
    memcpy(&ptr1[0][0][0][0], &ptr[0][0][0][0], 3*2*2*3*sizeof(double));

    std::vector<int> info = {3,2,2,3,0};
    Field field(&ptr[0][0][0][0],info);

    std::vector<int> info1 = {3,2,2,3,1};
    Field field1(&ptr1[0][0][0][0],info1);

    int *iPtr, *iPtr1;
    Malloc(iPtr, 3);
    Malloc(iPtr1, 3);
    Constant constant(iPtr, 3);
    Constant constant1(iPtr1, 3);

    // if Func do not add std::ref
    std::thread t1(Func, std::ref(constant), std::ref(field));
    std::thread t2(Func, std::ref(constant1), std::ref(field1));

    t1.join();
    t2.join();
}
