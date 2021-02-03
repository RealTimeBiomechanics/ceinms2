#include <string>
#include <iostream>
namespace ceinms {
template<typename T>
bool runTest(T testFun, std::string name) {
    bool failed = testFun();
    std::cout << name << ": ";
    if (failed)
        std::cout << "FAILED\n";
    else
        std::cout << "PASSED\n";
    return failed;
}
}// namespace ceinms
