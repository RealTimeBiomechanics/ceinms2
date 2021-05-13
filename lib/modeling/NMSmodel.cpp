#include "ceinms2/NMSmodel.h"

namespace ceinms {
std::ostream &operator<<(std::ostream &os, const Socket &sock) {
    os << sock.getName();
    if (sock.hasSlot()) { os << "." << sock.getSlot(); }
    return os;
}
}