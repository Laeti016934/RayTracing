#ifndef Sphere_H
#define Sphere_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySphereIntersection{
    bool intersectionExists;
    float t;
    float theta,phi;
    Vec3 intersection;
    Vec3 secondintersection;
    Vec3 normal;
};

static
Vec3 SphericalCoordinatesToEuclidean( Vec3 ThetaPhiR ) {
    return ThetaPhiR[2] * Vec3( cos(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[0]) * cos(ThetaPhiR[1]) , sin(ThetaPhiR[1]) );
}
static
Vec3 SphericalCoordinatesToEuclidean( float theta , float phi ) {
    return Vec3( cos(theta) * cos(phi) , sin(theta) * cos(phi) , sin(phi) );
}

static
Vec3 EuclideanCoordinatesToSpherical( Vec3 xyz ) {
    float R = xyz.length();
    float phi = asin( xyz[2] / R );
    float theta = atan2( xyz[1] , xyz[0] );
    return Vec3( theta , phi , R );
}



class Sphere : public Mesh {
public:
    Vec3 m_center;
    float m_radius;

    Sphere() : Mesh() {}
    Sphere(Vec3 c , float r) : Mesh() , m_center(c) , m_radius(r) {}

    void build_arrays(){
        unsigned int nTheta = 20 , nPhi = 20;
        positions_array.resize(3 * nTheta * nPhi );
        normalsArray.resize(3 * nTheta * nPhi );
        uvs_array.resize(2 * nTheta * nPhi );
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta ; ++thetaIt ) {
            float u = (float)(thetaIt) / (float)(nTheta-1);
            float theta = u * 2 * M_PI;
            for( unsigned int phiIt = 0 ; phiIt < nPhi ; ++phiIt ) {
                unsigned int vertexIndex = thetaIt + phiIt * nTheta;
                float v = (float)(phiIt) / (float)(nPhi-1);
                float phi = - M_PI/2.0 + v * M_PI;
                Vec3 xyz = SphericalCoordinatesToEuclidean( theta , phi );
                positions_array[ 3 * vertexIndex + 0 ] = m_center[0] + m_radius * xyz[0];
                positions_array[ 3 * vertexIndex + 1 ] = m_center[1] + m_radius * xyz[1];
                positions_array[ 3 * vertexIndex + 2 ] = m_center[2] + m_radius * xyz[2];
                normalsArray[ 3 * vertexIndex + 0 ] = xyz[0];
                normalsArray[ 3 * vertexIndex + 1 ] = xyz[1];
                normalsArray[ 3 * vertexIndex + 2 ] = xyz[2];
                uvs_array[ 2 * vertexIndex + 0 ] = u;
                uvs_array[ 2 * vertexIndex + 1 ] = v;
            }
        }
        triangles_array.clear();
        for( unsigned int thetaIt = 0 ; thetaIt < nTheta - 1 ; ++thetaIt ) {
            for( unsigned int phiIt = 0 ; phiIt < nPhi - 1 ; ++phiIt ) {
                unsigned int vertexuv = thetaIt + phiIt * nTheta;
                unsigned int vertexUv = thetaIt + 1 + phiIt * nTheta;
                unsigned int vertexuV = thetaIt + (phiIt+1) * nTheta;
                unsigned int vertexUV = thetaIt + 1 + (phiIt+1) * nTheta;
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuv );
                triangles_array.push_back( vertexUV );
                triangles_array.push_back( vertexuV );
            }
        }
    }


    RaySphereIntersection intersect(const Ray &ray) const {
        RaySphereIntersection intersection;
        //TODO calcul l'intersection rayon sphere
        
        //Equation de la sphère : ||x+c||^2 - r^2 = 0
        //Equation rayon-sphère : t^2 d.d + 2t d.(o-c) + ||o-c||^2 - r^2 = 0
        
        //Equation Quadratique du second degré en t :
        //Delta = b^2 - 4ac
        //Avec a = d.d
        //     b = 2.d.(o-c)
        //     c = ||o-c||^2 - r^2

        Vec3 rayDirection = ray.direction();
        Vec3 rayOrigin = ray.origin();
        Vec3 diff_Origin_Center = rayOrigin - m_center;

        float a = Vec3::dot(rayDirection, rayDirection);
        float b = Vec3::dot(2 * rayDirection, (diff_Origin_Center));
        float c = (diff_Origin_Center.norm() * diff_Origin_Center.norm()) - (m_radius * m_radius);

        float delta = b*b - 4*a*c;

        if(delta > 0){
            float sqrtDelta = sqrt(delta);
            float a2 = 2*a;
            float t1 = (-b - sqrtDelta) / (a2); //t1 = (-b - racineCarre(b^2 - 4ac)) / 2a
            float t2 = (-b + sqrtDelta) / (a2); //t2 = (-b + racineCarre(b^2 - 4ac)) / 2a

            if (t1 > 0 && t2 > 0) {
                //Le plus petit t positif
                intersection.t = std::min(t1, t2);
                intersection.intersectionExists = true;

                //Calcule de l'intersection la plus proche
                intersection.intersection = rayOrigin + intersection.t * rayDirection;
                intersection.normal = (intersection.intersection - m_center);
                intersection.normal.normalize();

                //Calcule du deuxième point d'intersection
                intersection.secondintersection = rayOrigin + std::max(t1, t2) * rayDirection;

            } else if (t1 > 0) {
                intersection.t = t1;
                intersection.intersectionExists = true;
                intersection.intersection = rayOrigin + t1 * rayDirection;
                intersection.normal = (intersection.intersection - m_center);
                intersection.normal.normalize();

            } else if (t2 > 0) {
                intersection.t = t2;
                intersection.intersectionExists = true;
                intersection.intersection = rayOrigin + t2 * rayDirection;
                intersection.normal = (intersection.intersection - m_center);
                intersection.normal.normalize();

            } else {
                // Si les deux sont négatifs, pas d'intersection "visible"
                intersection.intersectionExists = false;
            }
        }
        else if(delta == 0){
            //t = -b / 2a
            intersection.t = -b / (2 * a);
            intersection.intersectionExists = (intersection.t > 0);
            if (intersection.intersectionExists) {
                intersection.intersection = rayOrigin + intersection.t * rayDirection;
                intersection.normal = (intersection.intersection - m_center);
                intersection.normal.normalize();
            }
        }
        else{
            //Pas d'intersection
            intersection.intersectionExists = false; 
        }

        return intersection;
    }
};
#endif
