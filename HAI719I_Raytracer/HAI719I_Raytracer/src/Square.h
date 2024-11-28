#ifndef SQUARE_H
#define SQUARE_H
#include "Vec3.h"
#include <vector>
#include "Mesh.h"
#include <cmath>

struct RaySquareIntersection{
    bool intersectionExists;
    float t;
    float u,v;
    Vec3 intersection;
    Vec3 normal;
};


class Square : public Mesh {
public:
    Vec3 m_normal;
    Vec3 m_bottom_left;
    Vec3 m_right_vector;
    Vec3 m_up_vector;

    Square() : Mesh() {}
    Square(Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
           float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) : Mesh() {
        setQuad(bottomLeft, rightVector, upVector, width, height, uMin, uMax, vMin, vMax);
    }

    void setQuad( Vec3 const & bottomLeft , Vec3 const & rightVector , Vec3 const & upVector , float width=1. , float height=1. ,
                  float uMin = 0.f , float uMax = 1.f , float vMin = 0.f , float vMax = 1.f) {
        m_right_vector = rightVector;
        m_up_vector = upVector;
        m_normal = Vec3::cross(rightVector , upVector);
        m_bottom_left = bottomLeft;

        m_normal.normalize();
        m_right_vector.normalize();
        m_up_vector.normalize();

        m_right_vector = m_right_vector*width;
        m_up_vector = m_up_vector*height;

        vertices.clear();
        vertices.resize(4);
        vertices[0].position = bottomLeft;                                      vertices[0].u = uMin; vertices[0].v = vMin;
        vertices[1].position = bottomLeft + m_right_vector;                     vertices[1].u = uMax; vertices[1].v = vMin;
        vertices[2].position = bottomLeft + m_right_vector + m_up_vector;       vertices[2].u = uMax; vertices[2].v = vMax;
        vertices[3].position = bottomLeft + m_up_vector;                        vertices[3].u = uMin; vertices[3].v = vMax;
        vertices[0].normal = vertices[1].normal = vertices[2].normal = vertices[3].normal = m_normal;
        triangles.clear();
        triangles.resize(2);
        triangles[0][0] = 0;
        triangles[0][1] = 1;
        triangles[0][2] = 2;
        triangles[1][0] = 0;
        triangles[1][1] = 2;
        triangles[1][2] = 3;


    }

    RaySquareIntersection intersect(const Ray &ray) const {
        //TODO calculer l'intersection rayon quad

        //Equation du plan
        // x.n - D = 0
        // x : point sur le plan
        // n : vecteur normal au plan
        // D : distance du plan par rapport à l'origine dans la direction n

        //Distance entre entre l'origine du rayon et le point d'intersection avec le plan du carré le long de la direction du rayon (t)
        // t = (D - o.n) / d.n
        // o : origine du rayon
        // d : direction unitaire du rayon        

        RaySquareIntersection intersection;

        //Largeur et hauteur du quad
        Vec3 vec_width = (this->vertices[1].position - this->vertices[0].position);
        Vec3 vec_height = (this->vertices[3].position - this->vertices[0].position);

        //Normale du quad (perpendiculaire au plan)
        Vec3 normal = Vec3::cross(vec_width , vec_height);

        //Distance t entre l'origine du rayon et le point d'intersection avec le plan du carré
        float t = Vec3::dot(this->vertices[0].position - ray.origin(), normal) / Vec3::dot(ray.direction(), normal);

        //Intersection derrière la caméra
        if(t < 0){
            intersection.intersectionExists = false;
        }
        //Intersection devant ou à l'origine de la caméra
        if(t >= 0){
            Vec3 inter_pt = ray.origin() + t * ray.direction();
            Vec3 bottomLeftInterPt = inter_pt - this->vertices[0].position;
            
            float width = vec_width.length();
            float height = vec_height.length();

            float proj_width = Vec3::dot(bottomLeftInterPt, vec_width)/width;
            float proj_height = Vec3::dot(bottomLeftInterPt, vec_height)/height;

            if( 0 < proj_width && proj_width < width && 0 < proj_height && proj_height < height){
                intersection.intersectionExists = true;
                intersection.intersection = inter_pt;
                intersection.t = t;
                intersection.normal = normal;
                intersection.normal.normalize();
            }else{
                intersection.intersectionExists = false;
            }
        }
        return intersection;
    }

};
#endif // SQUARE_H
