#include "ceinms2/ElectromechanicalDelay.h"

namespace ceinms {

void connectSocket(const Excitation &parent, ElectromechanicalDelay &child) {
    child.setInput(parent);
}
}