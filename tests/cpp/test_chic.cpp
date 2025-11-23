#include "chic.hpp"
#include <cassert>

int main() {
    chic::Example ex;
    assert(ex.add(2, 3) == 5);
    return 0;
}