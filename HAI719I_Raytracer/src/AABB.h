#ifndef AABB_H
#define AABB_H

#include <limits>
#include <cmath>
#include <algorithm>
#include "Ray.h"
#include "Vec3.h"

struct AABB{
    Vec3 maxCoor;
    Vec3 minCoor;

    AABB(){
        maxCoor = Vec3(
            -std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity()
        );

        minCoor = Vec3(
            std::numeric_limits<float>::infinity(),
            std::numeric_limits<float>::infinity(),
            std::numeric_limits<float>::infinity()
        );
    }

    void expandByPoint(const Vec3& p);

    bool intersectRayAABB(const Ray& ray, float& t_entry, float& t_exit)
    {
        // Intervalle global du rayon
        t_entry = -std::numeric_limits<float>::infinity();
        t_exit  =  std::numeric_limits<float>::infinity();

        // Pour chaque axe : X = 0, Y = 1, Z = 2
        for (int axis = 0; axis < 3; ++axis) {

            float o = ray.origin()[axis];
            float d = ray.direction()[axis];
            float minA = minCoor[axis];
            float maxA = maxCoor[axis];

            // Rayon parallèle aux plans de la slab
            if (std::abs(d) < 1e-8f) {
                // Origine hors de la slab → pas d'intersection
                if (o < minA || o > maxA) {
                    return false;
                }
                // Sinon : axe ignoré → aucune contrainte
                continue;
            }

            // Rayon non parallèle : calcul des t
            float t1 = (minA - o) / d;
            float t2 = (maxA - o) / d;

            // Ordonner t1 et t2
            if (t1 > t2) std::swap(t1, t2);

            // Intersection avec l’intervalle global
            t_entry = std::max(t_entry, t1);
            t_exit  = std::min(t_exit,  t2);

            // Rejet anticipé
            if (t_entry > t_exit) {
                return false;
            }
        }

        // Intersection valide seulement si elle est devant le rayon
        return t_exit >= 0.0f;
    }
};

inline void AABB::expandByPoint(const Vec3& p) {
    minCoor[0] = std::min(minCoor[0], p[0]);
    minCoor[1] = std::min(minCoor[1], p[1]);
    minCoor[2] = std::min(minCoor[2], p[2]);

    maxCoor[0] = std::max(maxCoor[0], p[0]);
    maxCoor[1] = std::max(maxCoor[1], p[1]);
    maxCoor[2] = std::max(maxCoor[2], p[2]);
}
    
#endif