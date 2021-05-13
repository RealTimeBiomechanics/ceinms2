#include "ceinms2/LinearActuator.h"

namespace ceinms {

	void connectSocket(const Activation &parent, LinearActuator &child) {
		child.setInput(parent);
	}

}