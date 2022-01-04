#include "ceinms2/Mileusnic2006MuscleSpindle.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>

int test1() {
    ceinms::Mileusnic2006MuscleSpindle spindle;
    return 1;
}

int main() {
    bool failed = false;
    failed |= ceinms::runTest(&test1, "Correctness test");
    return failed;
}