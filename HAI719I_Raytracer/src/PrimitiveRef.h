#ifndef PRIMITIVEREF_H
#define PRIMITIVEREF_H

enum class PrimitiveType {
    Triangle,
    Sphere,
    Mesh
};

struct PrimitiveRef {
    PrimitiveType type;
    int objectIndex;     // index dans scene.getMeshes / getSpheres / getSquares
    int primitiveIndex;  // index du triangle (si Mesh), -1 sinon

    bool operator==(const PrimitiveRef& other) const {
        return type == other.type
            && objectIndex == other.objectIndex
            && primitiveIndex == other.primitiveIndex;
    }
};

#endif
