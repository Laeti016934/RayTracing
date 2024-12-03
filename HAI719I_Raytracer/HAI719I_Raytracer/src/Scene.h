#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include <iostream>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>


using namespace std;

const int MESHE = 0;
const int SPHERE = 1;
const int SQUARE = 2;

const Vec3 globalAmbientLight(0.2, 0.2, 0.2);

enum LightType {
    LightType_Spherical,
    LightType_Quad
};


struct Light {
    Vec3 material;
    bool isInCamSpace;
    LightType type;

    Vec3 pos;
    float radius;

    Mesh quad;

    float powerCorrection;

    Light() : powerCorrection(1.0) {}

};

struct RaySceneIntersection{
    bool intersectionExists;
    unsigned int typeOfIntersectedObject;
    unsigned int objectIndex;
    float t;
    RayTriangleIntersection rayMeshIntersection;
    RaySphereIntersection raySphereIntersection;
    RaySquareIntersection raySquareIntersection;
    RaySceneIntersection() : intersectionExists(false) , t(FLT_MAX) {}
};



class Scene {
    std::vector< Mesh > meshes;
    std::vector< Sphere > spheres;
    std::vector< Square > squares;
    std::vector< Light > lights;

public:


    Scene() {
    }

    void draw() {
        // iterer sur l'ensemble des objets, et faire leur rendu :
        for( unsigned int It = 0 ; It < meshes.size() ; ++It ) {
            Mesh const & mesh = meshes[It];
            mesh.draw();
        }
        for( unsigned int It = 0 ; It < spheres.size() ; ++It ) {
            Sphere const & sphere = spheres[It];
            sphere.draw();
        }
        for( unsigned int It = 0 ; It < squares.size() ; ++It ) {
            Square const & square = squares[It];
            square.draw();
        }
    }




    RaySceneIntersection computeIntersection(Ray const & ray) {
        RaySceneIntersection result;
        //TODO calculer les intersections avec les objets de la scene et garder la plus proche

        for (int i = 0; i < meshes.size(); i++) {
            RayTriangleIntersection RSI = meshes[i].intersect(ray);
            if (RSI.intersectionExists && RSI.t < result.t) {
                result.rayMeshIntersection = RSI;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 0;
                result.t = RSI.t;
            }
        }

        for(int i = 0; i < spheres.size(); i++){
            RaySphereIntersection RSI = spheres[i].intersect(ray);
            if(RSI.intersectionExists == true && RSI.t < result.t ){
                result.raySphereIntersection = RSI;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 1;
                result.t = RSI.t;
            }
        }

        for (int i = 0; i < squares.size(); i++) {
            RaySquareIntersection RSI = squares[i].intersect(ray);
            if (RSI.intersectionExists && RSI.t < result.t) {
                result.raySquareIntersection = RSI;
                result.intersectionExists = true;
                result.objectIndex = i;
                result.typeOfIntersectedObject = 2;
                result.t = RSI.t;
            }
        }

        return result;
    }

    Vec3 calculateAmbient(float Ka) {
        return Ka * globalAmbientLight;
    }

    Vec3 phong(RaySceneIntersection intersectionObjet) {
        Vec3 color(0., 0., 0.);
        Vec3 ambient(0., 0., 0.);
        Vec3 ambientColor(0., 0., 0.);
        Vec3 sommeSpecular(0., 0., 0.);
        Vec3 sommeDifraction(0., 0., 0.);

        // Réflexion Ambiante
        float Ka = 0.35; // Coefficient de réflexion ambiante

        ambientColor = calculateAmbient(Ka);

        if (!intersectionObjet.intersectionExists) {
            return Vec3(0., 0., 0.); // Pas d'intersection
        }

        //P : Point d'intersection entre le rayon et l'objet
        //(L.N) : angle entre la source de lumière et la normale
        //(R.V) : angle entre les directions de réflexion et de la vue
        Vec3 P, L, N, V, R;

        //Parcours des sources de lumière
        for (int i_Light = 0; i_Light < lights.size(); i_Light++) {
            Light &light = lights[i_Light];

            // Calcul des contributions pour chaque type d'objet
            if (intersectionObjet.typeOfIntersectedObject == SPHERE) {
                RaySphereIntersection result_tmp = intersectionObjet.raySphereIntersection;
                Sphere &sph = spheres[intersectionObjet.objectIndex];

                P = result_tmp.intersection;
                N = result_tmp.normal;
                N.normalize();

            } else if (intersectionObjet.typeOfIntersectedObject == SQUARE) {
                RaySquareIntersection result_tmp = intersectionObjet.raySquareIntersection;
                Square &squ = squares[intersectionObjet.objectIndex];

                P = result_tmp.intersection;
                N = result_tmp.normal;
                N.normalize();

            } else {
                continue; // Autres types d'objets
            }

            L = (light.pos - P); // Direction vers la lumière
            L.normalize();

            V = P * -1;           // Direction vers la caméra
            V.normalize();

            R = (2. * Vec3::dot(N, L) * N - L); // Vecteur réfléchi
            R.normalize();

            // Vérifie si le point P est dans l'ombre de la lumière
            Ray shadowRay(P + N * 0.001, L); // Rayon légèrement décalé pour éviter l'autointersection
            RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);

            if (shadowIntersection.intersectionExists && 
                shadowIntersection.t < (light.pos - P).length()) {
                // Si un objet bloque la lumière, ignorer les contributions diffuse et spéculaire
                continue;
            }

            for (int i = 0; i < 3; i++) {
                float Isd = light.material[i]; // Intensité de la lumière diffuse
                float Kd, Ks; // Coefficients de réflexion diffuse et spéculaire
                float shininess;

                if (intersectionObjet.typeOfIntersectedObject == SPHERE) {
                    Sphere &sph = spheres[intersectionObjet.objectIndex];
                    Kd = sph.material.diffuse_material[i];
                    Ks = sph.material.specular_material[i];
                    shininess = sph.material.shininess;

                } else if (intersectionObjet.typeOfIntersectedObject == SQUARE) {
                    Square &squ = squares[intersectionObjet.objectIndex];
                    Kd = squ.material.diffuse_material[i];
                    Ks = squ.material.specular_material[i];
                    shininess = squ.material.shininess;
                }

                sommeDifraction[i] += Isd * Kd * std::max(0.f, Vec3::dot(L, N));
                sommeSpecular[i] += Isd * Ks * pow(std::max(0.f, Vec3::dot(R, V)), shininess);
            }
        }

        for (int i = 0; i < 3; i++) {
            color[i] = ambientColor[i] + sommeDifraction[i] + sommeSpecular[i];
        }

        return color;
    }


    Vec3 rayTraceRecursive( Ray ray , int NRemainingBounces ) {
        Vec3 color;
    
        //Calcule de l'intersection avec la scène
        RaySceneIntersection raySceneIntersection = computeIntersection(ray);
        
        if(!raySceneIntersection.intersectionExists){
                color = Vec3(0.,0.,0.);
                return color;
        }

        //PHASE 1
        /*
        if(raySceneIntersection.intersectionExists){
            if(raySceneIntersection.typeOfIntersectedObject == 0){
                color = meshes[raySceneIntersection.objectIndex].material.diffuse_material;

            }else if(raySceneIntersection.typeOfIntersectedObject == 1){
                color = spheres[raySceneIntersection.objectIndex].material.diffuse_material;

            }else if(raySceneIntersection.typeOfIntersectedObject == 2){
                color = squares[raySceneIntersection.objectIndex].material.diffuse_material;

            }
        }
        */
        
        //PHASE 2
        if (raySceneIntersection.intersectionExists)
        {
            color = phong(raySceneIntersection);
        }
        
        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
        //TODO appeler la fonction recursive
        int initialBounceCount = 10;

        Vec3 color = rayTraceRecursive(rayStart, initialBounceCount);

        return color;
    }

    void setup_single_sphere() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1. , 0. , 0.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
        {
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(0. , 0. , 1.);
            s.m_radius = 1.f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3( 0.2,0.2,0.2 );
            s.material.shininess = 20;
        }
    }

    void setup_single_square() {
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        {
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3(-5,5,5);
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }

        {
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1.5, -1., 0.), Vec3(0.5, 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.,1.,0. );
            s.material.specular_material = Vec3( 0.8,0.8,0.8 );
            s.material.shininess = 20;
        }
    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        { //Light
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( 0.0, 1.5, 0.0 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1,1,1);
            light.isInCamSpace = false;
        }

        { //Light 2
            lights.resize( lights.size() + 1 );
            Light & light = lights[lights.size() - 1];
            light.pos = Vec3( -1.9, -1.9, 1.9 );
            light.radius = 2.5f;
            light.powerCorrection = 2.f;
            light.type = LightType_Spherical;
            light.material = Vec3(1.,0.,0.);
            light.isInCamSpace = false;
        }

        { //Back Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.5,0. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
        }

        { //Left Wall

            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.scale(Vec3(2., 2., 1.));
            s.translate(Vec3(0., 0., -2.));
            s.rotate_y(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.,0.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
        }

        { //Right Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,1.0,0.0 );
            s.material.specular_material = Vec3( 0.0,1.0,0.0 );
            s.material.shininess = 16;
        }

        { //Floor
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(-90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.5,0.5,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }

        { //Ceiling
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_x(90);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 0.0,0.,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }
        /*
        { //Front Wall
            squares.resize( squares.size() + 1 );
            Square & s = squares[squares.size() - 1];
            s.setQuad(Vec3(-1., -1., 0.), Vec3(1., 0, 0.), Vec3(0., 1, 0.), 2., 2.);
            s.translate(Vec3(0., 0., -2.));
            s.scale(Vec3(2., 2., 1.));
            s.rotate_y(180);
            s.build_arrays();
            s.material.diffuse_material = Vec3( 1.0,1.0,1.0 );
            s.material.specular_material = Vec3( 1.0,1.0,1.0 );
            s.material.shininess = 16;
        }*/


        { //GLASS Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,0. );
            s.material.specular_material = Vec3( 1.,0.,0. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }


        { //MIRRORED Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,0.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.;
            s.material.index_medium = 0.;
        }

    }

};



#endif
