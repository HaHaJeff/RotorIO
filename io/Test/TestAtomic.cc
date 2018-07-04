#include "AtomicPointer.h"
#include <iostream>
#include <thread>

using namespace std;

void fun1(void* arg) {
    AtomicPointer* ap = static_cast<AtomicPointer*>(arg);

}

void fun2(void* arg) {

}

int main() {
	int id = 0;
	AtomicPointer ap(&id);

	int* ptr = static_cast<int*>(ap.NoBarrier_Load());
	std::cout << *ptr << std::endl;
}
