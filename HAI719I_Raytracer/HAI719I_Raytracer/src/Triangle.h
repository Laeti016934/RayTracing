#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "Vec3.h"
#include "Ray.h"
#include "Plane.h"

struct RayTriangleIntersection{
    bool intersectionExists;
    float t;
    float w0,w1,w2;
    unsigned int tIndex;
    Vec3 intersection;
    Vec3 normal;
};

class Triangle {
private:
    Vec3 m_c[3] , m_normal;
    float area;
public:
    Triangle() {}
    Triangle( Vec3 const & c0 , Vec3 const & c1 , Vec3 const & c2 ) {
        m_c[0] = c0;
        m_c[1] = c1;
        m_c[2] = c2;
        updateAreaAndNormal();
    }
    void updateAreaAndNormal() {
        Vec3 nNotNormalized = Vec3::cross( m_c[1] - m_c[0] , m_c[2] - m_c[0] );
        float norm = nNotNormalized.length();
        m_normal = nNotNormalized / norm;
        area = norm / 2.f;
    }
    void setC0( Vec3 const & c0 ) { m_c[0] = c0; } // remember to update the area and normal afterwards!
    void setC1( Vec3 const & c1 ) { m_c[1] = c1; } // remember to update the area and normal afterwards!
    void setC2( Vec3 const & c2 ) { m_c[2] = c2; } // remember to update the area and normal afterwards!
    
    Vec3 const & normal() const { return m_normal; }
    
    Vec3 projectOnSupportPlane( Vec3 const & p ) const {
        Vec3 result;

        Line line(m_c[0], m_normal);  // Crée une ligne passant par un sommet et parallèle à la normale du triangle
        result = line.project(p);     // Projette le point p sur cette ligne

        return result;
    }

    float squareDistanceToSupportPlane( Vec3 const & p ) const {
        float result;
        
        // On projette le point p sur le plan du triangle
        Vec3 projection = projectOnSupportPlane(p);
        result = (p - projection).length(); // Retourne la distance carrée

        return result;
    }

    float distanceToSupportPlane( Vec3 const & p ) const { return sqrt( squareDistanceToSupportPlane(p) ); }
    
    bool isParallelTo( Line const & L ) const {
        bool result;
        
        // On vérifie si la normale au triangle est parallèle à la direction de la ligne.
        float dotProduct = Vec3::dot(m_normal, L.direction());
        
        result = std::abs(dotProduct) < 1e-6f; // Si le produit scalaire est proche de zéro, ils sont parallèles.

        return result;
    }

    Vec3 getIntersectionPointWithSupportPlane( Line const & L ) const {
        // you should check first that the line is not parallel to the plane!
        Vec3 result;
        
        // On vérifie d'abord si la ligne est parallèle au plan du triangle
        if (isParallelTo(L)) {
            // Si elle est parallèle, on retourne un vecteur nulle ou un autre indicateur
            std::cerr << "La ligne est parallèle au plan." << std::endl;
            return Vec3(); // Retourne un point nul
        }
        
        // Calcul du facteur t pour l'intersection de la ligne avec le plan
        float t = Vec3::dot(m_c[0] - L.origin(), m_normal) / Vec3::dot(L.direction(), m_normal);
        
        // Retourne le point d'intersection
        result = L.origin() + t * L.direction();

        return result;
    }

    void computeBarycentricCoordinates( Vec3 const & p , float & u0 , float & u1 , float & u2 ) const {
        // Vecteurs de l'élément du triangle
        Vec3 v0 = m_c[1] - m_c[0];
        Vec3 v1 = m_c[2] - m_c[0];
        Vec3 v2 = p - m_c[0];

        // Calcul du déterminant de la matrice pour l'inverser et obtenir les barycentriques
        float d00 = Vec3::dot(v0, v0);
        float d01 = Vec3::dot(v0, v1);
        float d11 = Vec3::dot(v1, v1);
        float d20 = Vec3::dot(v2, v0);
        float d21 = Vec3::dot(v2, v1);
        float denom = d00 * d11 - d01 * d01;

        // Calcul des coordonnées barycentriques
        u1 = (d11 * d20 - d01 * d21) / denom;
        u2 = (d00 * d21 - d01 * d20) / denom;
        u0 = 1.0f - u1 - u2;
    }

    RayTriangleIntersection getIntersection( Ray const & ray ) const {
        RayTriangleIntersection result;
        result.intersectionExists = false;
        
        // 1) check that the ray is not parallel to the triangle:
        float denominator = Vec3::dot(m_normal, ray.direction());
        if (fabs(denominator) < 1e-6f) {
            return result; // Pas d'intersection, rayon parallèle au triangle
        }

        // 2) check that the triangle is "in front of" the ray:
        Vec3 edgeToRayOrigin = m_c[0] - ray.origin();
        float t = Vec3::dot(m_normal, edgeToRayOrigin) / denominator;

        if (t < 1e-6f) { // Ajouter une petite tolérance
            return result; // Intersection derrière ou trop proche
        }

        // 3) check that the intersection point is inside the triangle:
        // CONVENTION: compute u,v such that p = w0*c0 + w1*c1 + w2*c2, check that 0 <= w0,w1,w2 <= 1
        float u0, u1, u2;
        Vec3 intersectionPoint = ray.origin() + t * ray.direction();
        computeBarycentricCoordinates(intersectionPoint, u0, u1, u2);

        // Vérifier que le point est à l'intérieur du triangle
        if (u0 < 0 || u1 < 0 || u2 < 0) {
            return result; // Intersection en dehors du triangle
        }

        // 4) Finally, if all conditions were met, then there is an intersection! :
        result.intersectionExists = true;
        result.t = t;
        result.w0 = u0;
        result.w1 = u1;
        result.w2 = u2;
        result.intersection = intersectionPoint;
        result.normal = m_normal; // Normalisée par construction

        return result;
    }
};
#endif
