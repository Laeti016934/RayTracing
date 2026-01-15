#ifndef HIT_H
#define HIT_H

#include <limits>
#include "PrimitiveRef.h"
#include "Vec3.h"

struct Hit {
    bool hit = false;
    float t = std::numeric_limits<float>::infinity();
    Vec3 position;
    Vec3 normal;

    PrimitiveRef primitive;
};

#endif