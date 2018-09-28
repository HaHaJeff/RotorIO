#include "checkpoint.h" #include "memory.h"

#include <iostream>
#include <string.h>
#include <thread>

void Func(Constant& constant, Field& field) {
  /*
    Checkpoint ck(constant, field, 1);
    POSIXIO io;
    Strategy* strategy = io.GetIOStrategy(static_cast<TYPE>(0));
    strategy->Open("test");
    */

    //ck.RestoreField(*strategy);
    //ck.SaveField(*strategy);
    SetCheckpoint(constant, field, 1);

    /*
    const Field& f= ck.GetField();
    double *p = f.GetField();
    std::vector<int> info = f.GetInfo();

    std::cout << std::this_thread::get_id() << std::endl;
    for (int i = 0; i < info[0]*info[1]*info[2]*info[3]; ++i) {
      std::cout << p[i] << std::endl;
    }
    std::cout << "block_id: " << info[4] << std::endl;

    const Constant& c = ck.GetConstant();
    int *p1 = c.GetConstant();

    for (int i = 0; i < c.GetSize(); i++) {
      std::cout << p1[i] << std::endl;
    }
    */
}

int main()
{
    double ****ptr;
    int a = 1;
    Malloc(ptr, 3,2,2,3);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
          for (int z = 0; z < 3; ++z) {
            ptr[i][j][k][z] = a++;
          }
        }
      }

    }

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

    for (int i = 0; i < 3; i++) {
      iPtr[i] = i;
      iPtr1[i] = i;
    }
    Constant constant(iPtr, 3);
    Constant constant1(iPtr1, 3);

    // if Func do not add std::ref
    std::thread t1(Func, std::ref(constant), std::ref(field));
//    std::thread t2(Func, std::ref(constant1), std::ref(field1));

    t1.join();
//    t2.join();
}
