#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include "Mesh.h"
#include "Sphere.h"
#include "Square.h"


#include <GL/glut.h>


using namespace std;

const int MESH = 0;
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

    void initialize_AABB(const Scene& scene){
        const std::vector<Mesh>& meshes = scene.getMeshes();
        const std::vector<Square>& squares = scene.getSquares();
        const std::vector<Sphere>& spheres = scene.getSpheres();

        for(int i = 0; i < meshes.size(); i++ ){
            for(int y = 0; y < meshes[i].vertices.size(); y++){

                // X Axis
                if(meshes[i].vertices[y].position[0] < minCoor[0]){
                    minCoor[0] = meshes[i].vertices[y].position[0];
                }
                if(meshes[i].vertices[y].position[0] > maxCoor[0]){
                    maxCoor[0] = meshes[i].vertices[y].position[0];
                }

                // Y Axis
                if(meshes[i].vertices[y].position[1] < minCoor[1]){
                    minCoor[1] = meshes[i].vertices[y].position[1];
                }
                if(meshes[i].vertices[y].position[1] > maxCoor[1]){
                    maxCoor[1] = meshes[i].vertices[y].position[1];
                }

                // Z Axis
                if(meshes[i].vertices[y].position[2] < minCoor[2]){
                    minCoor[2] = meshes[i].vertices[y].position[2];
                }
                if(meshes[i].vertices[y].position[2] > maxCoor[2]){
                    maxCoor[2] = meshes[i].vertices[y].position[2];
                }
            }
        }

        for(int i = 0; i < squares.size(); i++){
            for (int y = 0; y < squares[i].vertices.size(); y++){

                // Axis X
                if(squares[i].vertices[y].position[0] < minCoor[0]){
                    minCoor[0] = squares[i].vertices[y].position[0];
                }
                if(squares[i].vertices[y].position[0] > maxCoor[0]){
                    maxCoor[0] = squares[i].vertices[y].position[0];
                }

                // Axis Y
                if(squares[i].vertices[y].position[1] < minCoor[1]){
                    minCoor[1] = squares[i].vertices[y].position[1];
                }
                if(squares[i].vertices[y].position[1] > maxCoor[1]){
                    maxCoor[1] = squares[i].vertices[y].position[1];
                }

                // Axis Z
                if(squares[i].vertices[y].position[2] < minCoor[2]){
                    minCoor[2] = squares[i].vertices[y].position[2];
                }
                if(squares[i].vertices[y].position[2] > maxCoor[2]){
                    maxCoor[2] = squares[i].vertices[y].position[2];
                }
            }
        }

        for(int i = 0; i < spheres.size(); i++){
            Vec3 minSph;
            minSph[0] = spheres[i].m_center[0] - spheres[i].m_radius;
            minSph[1] = spheres[i].m_center[1] - spheres[i].m_radius;
            minSph[2] = spheres[i].m_center[2] - spheres[i].m_radius;

            Vec3 maxSph;
            maxSph[0] = spheres[i].m_center[0] + spheres[i].m_radius;
            maxSph[1] = spheres[i].m_center[1] + spheres[i].m_radius;
            maxSph[2] = spheres[i].m_center[2] + spheres[i].m_radius;

            // Axis X
            if(minSph[0] < minCoor[0]){
                minCoor[0] = minSph[0];
            }
            if(maxSph[0] > maxCoor[0]){
                maxCoor[0] = maxSph[0];
            }

            // Axis Y
            if(minSph[1] < minCoor[1]){
                minCoor[1] = minSph[1];
            }
            if(maxSph[1] > maxCoor[1]){
                maxCoor[1] = maxSph[1];
            }

            // Axis Z
            if(minSph[2] < minCoor[2]){
                minCoor[2] = minSph[2];
            }
            if(maxSph[2] > maxCoor[2]){
                maxCoor[2] = maxSph[2];
            }
        }
    }

    bool intersectRayAABB(const Ray& ray, const AABB& aabb,
                        float& t_entry, float& t_exit)
    {
        // Intervalle global du rayon
        t_entry = -std::numeric_limits<float>::infinity();
        t_exit  =  std::numeric_limits<float>::infinity();

        // Pour chaque axe : X = 0, Y = 1, Z = 2
        for (int axis = 0; axis < 3; ++axis) {

            float o = ray.origin()[axis];
            float d = ray.direction()[axis];
            float minA = aabb.minCoor[axis];
            float maxA = aabb.maxCoor[axis];

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

enum class PrimitiveType {
    Triangle,
    Sphere,
    Mesh
};

struct PrimitiveRef {
    PrimitiveType type;

    // Index vers les données réelles
    int objectIndex;     // index du mesh / sphere dans la scène
    int primitiveIndex;  // index du triangle (inutile pour une sphère)
};

struct Hit {
    bool hit = false;
    float t = std::numeric_limits<float>::infinity();
    Vec3 position;
    Vec3 normal;

    PrimitiveRef primitive;
};

struct Node
{
    // Spatial data : AABB
    AABB boundingBox;

    // If internal node
    int splittingAxis;    // Splitting axis (X, Y or Z)
    float cuttingPos;       // Cutting position
    Node* leftChild;        // Left child 
    Node* rightChild;       // Right child

    // If leaf node
    std::vector<PrimitiveRef> primitives;    // List of primitives 

    bool leaf; // Node State 

    Node()
        : splittingAxis(-1),
          cuttingPos(0.f),
          leftChild(nullptr),
          rightChild(nullptr),
          leaf(false)
    {}
};


struct KDTree
{
    Node* root;                 // Root node
    
    int maxDepth;               // Maximum depth
    int maxPrimitivesPerLeaf;   // Max number of primitives per leaf
    float epsilon;              // Epsilon of minimal size

    KDTree();

    // Construction globale
    void buildKDTree(const Scene& scene){
        // Récupération des primitives
        std::vector<PrimitiveRef> allPrimitives;

        // Spheres
        const auto& spheres = scene.getSpheres();
        for (int i = 0; i < spheres.size(); i++) {
            PrimitiveRef prim;
            prim.type = PrimitiveType::Sphere;
            prim.objectIndex = i;
            prim.primitiveIndex = -1; // inutile pour sphère
            allPrimitives.push_back(prim);
        }

        // Squares
        const auto& squares = scene.getSquares();
        for (int i = 0; i < squares.size(); i++) {
            PrimitiveRef prim;
            prim.type = PrimitiveType::Triangle;
            prim.objectIndex = i;
            prim.primitiveIndex = -1;
            allPrimitives.push_back(prim);
        }

        // Mesh
        const auto& meshes = scene.getMeshes();
        for (int i = 0; i < meshes.size(); i++) {
            PrimitiveRef prim;
            prim.type = PrimitiveType::Mesh;
            prim.objectIndex = i;
            prim.primitiveIndex = -1;
            allPrimitives.push_back(prim);
        }

        // AABB globale de la scene
        AABB globalAABB;
        globalAABB.initialize_AABB(scene);

        // Construction récursive du KDTree
        root = buildNode(globalAABB, allPrimitives, scene, 0);
    }

    bool intersectPrimitive(
    const PrimitiveRef& prim,
    const Scene& scene,
    const Ray& ray,
    Hit& hit
    ) const {
        bool hitRes = false;
        switch (prim.type) {
            case PrimitiveType::Sphere: {
                const Sphere& s = scene.getSpheres()[prim.objectIndex];
                RaySphereIntersection res = s.intersect(ray);
                hitRes = res.intersectionExists;
                if (hitRes) {
                    hit.hit = true;
                    hit.t = res.t;
                    hit.position = res.intersection;
                    hit.normal = res.normal;
                    hit.primitive = prim;
                }
                break;
            }

            case PrimitiveType::Triangle: {
                const Square& sq = scene.getSquares()[prim.objectIndex];
                RaySquareIntersection res = sq.intersect(ray);
                hitRes = res.intersectionExists;
                if (hitRes) {
                    hit.hit = true;
                    hit.t = res.t;
                    hit.position = res.intersection;
                    hit.normal = res.normal;
                    hit.primitive = prim;
                }
                break;
            }

            case PrimitiveType::Mesh: {
                const Mesh& m = scene.getMeshes()[prim.objectIndex];
                RayTriangleIntersection res = m.intersect(ray);
                hitRes = res.intersectionExists;
                if (hitRes) {
                    hit.hit = true;
                    hit.t = res.t;
                    hit.position = res.intersection;
                    hit.normal = res.normal;
                    hit.primitive = prim;
                }
                break;
            }
        }
        return hitRes;
    }



private :
    // Construction récursive
    Node* buildNode(AABB &aabb, std::vector<PrimitiveRef> &prim, const Scene& scene, int depth){
        Node* n = new Node();
        n->boundingBox = aabb;

        //primitives ≤ seuil ?
        if(prim.size() <= maxPrimitivesPerLeaf){
            n->leaf = true;
            n->primitives = prim;
            return n;
        }
        //profondeur max atteinte ?
        if(depth >= maxDepth){
            n->leaf = true;
            n->primitives = prim;
            return n;
        }

        float sizeX = aabb.maxCoor[0] - aabb.minCoor[0];
        float sizeY = aabb.maxCoor[1] - aabb.minCoor[1];
        float sizeZ = aabb.maxCoor[2] - aabb.minCoor[2];
        float aabbSize = max(sizeX, sizeY, sizeZ);

        if(aabbSize < epsilon){
            n->leaf = true;
            n->primitives = prim;
            return n;
        }

        // Sinon : nœud interne
        n->leaf = false;

        // Axe de coupe
        if(sizeX >= sizeY && sizeX >= sizeZ){
            n->splittingAxis = 0;
        }
        else if(sizeY >= sizeX && sizeY >= sizeZ){
            n->splittingAxis = 1;
        }
        else{
            n->splittingAxis = 2;
        }

        // Calcul de la position de coupe (Couper au milieu de la boîte englobante sur l’axe choisi)
        switch (n->splittingAxis){
            case 0:             // X Axis
                n->cuttingPos = (aabb.minCoor[0] + aabb.maxCoor[0]) / 2;
                break;
            
            case 1:             // Y Axis
                n->cuttingPos = (aabb.minCoor[1] + aabb.maxCoor[1]) / 2;
                break;

            case 2:             // Z Axis
                n->cuttingPos = (aabb.minCoor[2] + aabb.maxCoor[2]) / 2;
                break;
            
            default:
                break;
        }

        // Construction des bounding boxes des enfants
        AABB aabbLeftChild = aabb;
        aabbLeftChild.maxCoor[n->splittingAxis] = n->cuttingPos;

        AABB aabbRightChild = aabb;
        aabbRightChild.minCoor[n->splittingAxis] = n->cuttingPos;

        std::vector<PrimitiveRef> leftPrimitives;
        std::vector<PrimitiveRef> rightPrimitives;
        std::vector<Mesh> const& meshes    = scene.getMeshes();
        std::vector<Square> const& squares = scene.getSquares();
        std::vector<Sphere> const& spheres = scene.getSpheres();

        for(int i = 0; i < prim.size(); i++){
            /*
            calculer pMin (sa borne minimale sur l’axe de coupe), pMax (sa borne maximale sur l’axe de coupe)
            tester les trois cas
            ajouter la primitive à la/les bonne(s) liste(s)
            */

            // Calculs de pMin et pMax 
            float pMin, pMax;
            int axis = n->splittingAxis;

            switch (prim[i].type) {
                case PrimitiveType::Sphere: {
                    /*
                    pMin = center[axis] - radius
                    pMax = center[axis] + radius
                    */
                    const Sphere& s = spheres[prim[i].objectIndex];
                    pMin = s.m_center[axis] - s.m_radius;
                    pMax = s.m_center[axis] + s.m_radius;
                    break;
                }
                case PrimitiveType::Triangle: {
                    /*
                    prendre les 3 sommets
                    pMin = min(v0, v1, v2)[axis]
                    pMax = max(v0, v1, v2)[axis]
                    */
                    const Square& t = squares[prim[i].objectIndex];
                    pMin = std::min({ t.vertices[0].position[axis],
                                    t.vertices[1].position[axis],
                                    t.vertices[2].position[axis]});

                    pMax = std::max({ t.vertices[0].position[axis],
                                    t.vertices[1].position[axis],
                                    t.vertices[2].position[axis]});
                    break;
                }
                default: // Mesh
                    /*
                    utiliser son AABB pré-calculée
                    prendre min/max sur l’axe
                    */
                    const Mesh& m = meshes[prim[i].objectIndex];
                    
                    // Initialiser pMin et pMax
                    pMin =  std::numeric_limits<float>::infinity();
                    pMax = -std::numeric_limits<float>::infinity();

                    // Parcourir tous les vertices pour trouver min/max sur l'axe choisi
                    for (int v = 0; v < m.vertices.size(); v++) {
                        float val = m.vertices[v].position[axis];
                        if (val < pMin) pMin = val;
                        if (val > pMax) pMax = val;
                    }

                    break;
            }

            // Répartition des primitives
            if (pMax <= n->cuttingPos) {
                leftPrimitives.push_back(prim[i]);
            } else if (pMin >= n->cuttingPos) {
                rightPrimitives.push_back(prim[i]);
            } else {
                // la primitive traverse la coupe → ajouter aux deux
                leftPrimitives.push_back(prim[i]);
                rightPrimitives.push_back(prim[i]);
            }
        }

        // Gestion des cas dégénérés
        bool degenerate = false;

        // Cas 1 : aucune séparation réelle
        if (leftPrimitives.size() == prim.size() || rightPrimitives.size() == prim.size()) {
            degenerate = true;
        }

        // Cas 2 : une des listes est vide
        if (leftPrimitives.empty() || rightPrimitives.empty()) {
            degenerate = true;
        }

        // Cas 3 : duplication totale
        if (leftPrimitives == prim && rightPrimitives == prim) {
            degenerate = true;
        }

        // Si partition dégénérée → créer une feuille
        if (degenerate) {
            n->leaf = true;
            n->primitives = prim;
            return n;
        }

        n->leftChild  = buildNode(aabbLeftChild, leftPrimitives, scene, depth + 1);
        n->rightChild = buildNode(aabbRightChild, rightPrimitives, scene, depth + 1);

        return n;
    }

    bool intersectNode(Node* node, const Scene& scene, const Ray& ray, float t_entry, float t_exit, Hit& hit) const {
        int axis = node->splittingAxis;
        // Feuille
        if (node->leaf) {
            bool hitSomething = false;
            for (const PrimitiveRef& p : node->primitives) {
                Hit tempHit;
                if (intersectPrimitive(p, scene, ray, tempHit)) {
                    hitSomething = true;
                    if (!hit.hit || tempHit.t < hit.t) {
                        hit = tempHit;
                    }
                }
            }
            return hitSomething;
        }
        
        //Noeud interne

        //Calcul de t_split (la valeur de t pour laquelle le rayon coupe ce plan)
        //t_split sert à savoir si, quand et dans quel ordre le rayon traverse les deux enfants du KD-tree.
        //Est-ce que le rayon passe du premier enfant au second avant ou après avoir quitté le nœud ?
        float t_split = (node->cuttingPos - ray.origin()[axis]) / ray.direction()[axis];
        /*
        Quand tu testes l’intersection du rayon avec l’AABB du nœud, tu obtiens :
            t_entry → entrée dans la boîte
            t_exit → sortie de la boîte

        Donc le rayon n’existe dans ce nœud que pour :
            t∈[t_entry,t_exit]
        */
        /* Étape 2 : définir :
            near child
            far child

            if (ray.direction()[axis] >= 0) {
                near = leftChild;
                far  = rightChild;
            } else {
                near = rightChild;
                far  = leftChild;
            }

            Parce que :
                leftChild = coordonnées plus petites
                rightChild = coordonnées plus grandes

            Si tu avances vers les valeurs croissantes (D > 0) :
                tu touches d’abord la zone des petites valeurs

            Si tu avances vers les valeurs décroissantes (D < 0) :
                tu touches d’abord la zone des grandes valeurs

        */

        Node* nearChild;
        Node* farChild;

        if (ray.direction()[axis] >= 0) {
            nearChild = node->leftChild;
            farChild  = node->rightChild;
        } else {
            nearChild = node->rightChild;
            farChild  = node->leftChild;
        }

        // Étape 3 : tester les cas

        // Cas 1 : Le rayon est déjà du bon côté dès qu’il entre dans le nœud (t_split <= t_entry)
        if(t_split <= t_entry){
            return intersectNode(farChild, scene, ray, t_entry, t_exit, hit);
        }
        // Cas 2 : Le rayon quitte le nœud avant de traverser le plan de coupe (t_split >= t_exit)
        else if(t_split >= t_exit){
            return intersectNode(nearChild, scene, ray, t_entry, t_exit, hit);
        }
        // Cas 3 : Le rayon traverse le plan de séparation (t_entry < t_split < t_exit)
        else{
            // Tester d’abord le near child
            bool hitNear = intersectNode(nearChild, scene, ray, t_entry, t_split, hit);
            bool hitFar = false;
            // Si on a touché quelque chose dans le near child et que c’est avant t_split
            if (!hitNear || hit.t > t_split) {
                hitFar = intersectNode(farChild, scene, ray, t_split, t_exit, hit);
            }
            
            return hitNear || hitFar;
        }
    }
};



void initialize_quad_light(Light &light, float width, float height) {
    if (light.type != LightType_Quad) {
        std::cerr << "La lumière doit être de type Quad pour utiliser cette fonction." << std::endl;
        return;
    }

    float halfWidth = width * 0.5f;
    float halfHeight = height * 0.5f;

    light.quad.vertices.clear();
    light.quad.vertices.push_back(MeshVertex(light.pos + Vec3(-halfWidth, -halfHeight, 0), Vec3(0, 0, 1))); // Bas gauche
    light.quad.vertices.push_back(MeshVertex(light.pos + Vec3(halfWidth, -halfHeight, 0), Vec3(0, 0, 1)));  // Bas droit
    light.quad.vertices.push_back(MeshVertex(light.pos + Vec3(halfWidth, halfHeight, 0), Vec3(0, 0, 1)));   // Haut droit
    light.quad.vertices.push_back(MeshVertex(light.pos + Vec3(-halfWidth, halfHeight, 0), Vec3(0, 0, 1)));  // Haut gauche

    light.quad.triangles.clear();
    light.quad.triangles.push_back(MeshTriangle(0, 1, 2)); // Triangle 1
    light.quad.triangles.push_back(MeshTriangle(0, 2, 3)); // Triangle 2

    light.quad.build_arrays();
}


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

    const std::vector<Mesh>& getMeshes() const;

    const std::vector<Sphere>& getSpheres() const;

    const std::vector<Square>& getSquares() const;

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


//-------------------------Eclairage Phong et ombres ------------------------------------


    Vec3 samplePointOnLight(const Light &light) {
        if (light.type != LightType_Quad) {
            std::cerr << "Échantillonnage impossible : la lumière n'est pas de type Quad." << std::endl;
            return light.pos; // Retourne simplement la position pour les autres types de lumière.
        }

        // Générer des coordonnées barycentriques aléatoires
        float u = static_cast<float>(rand()) / RAND_MAX;
        float v = static_cast<float>(rand()) / RAND_MAX;

        // Assurer que u + v <= 1 pour le premier triangle
        if (u + v > 1.0f) {
            u = 1.0f - u;
            v = 1.0f - v;
        }

        const Mesh &quad = light.quad;

        // Calculer la position interpolée dans le triangle (0, 1, 2)
        Vec3 p0 = quad.vertices[quad.triangles[0][0]].position;
        Vec3 p1 = quad.vertices[quad.triangles[0][1]].position;
        Vec3 p2 = quad.vertices[quad.triangles[0][2]].position;

        return (1 - u - v) * p0 + u * p1 + v * p2;
    }

    Vec3 calculateAmbient(float Ka) {
        return Ka * globalAmbientLight;
    }


    void getMaterialProperties(int objectType, int objectIndex, int channel, float &Kd, float &Ks, float &shininess) {
        if (objectType == SPHERE) {
            Sphere &sph = spheres[objectIndex];
            Kd = sph.material.diffuse_material[channel];
            Ks = sph.material.specular_material[channel];
            shininess = sph.material.shininess;
        } else if (objectType == SQUARE) {
            Square &squ = squares[objectIndex];
            Kd = squ.material.diffuse_material[channel];
            Ks = squ.material.specular_material[channel];
            shininess = squ.material.shininess;
        } else if (objectType == MESH) {
            Mesh &mesh = meshes[objectIndex];
            Material material = mesh.material;
            Kd = material.diffuse_material[channel];
            Ks = material.specular_material[channel];
            shininess = material.shininess;
        }
    }


    void computeLightContribution(Vec3 &summedColor, Light &light, Vec3 &P, Vec3 &N, Vec3 &V, Vec3 &L, int objectType, int objectIndex) {
        Vec3 R = (2. * Vec3::dot(N, L) * N - L);
        R.normalize();

        for (int i = 0; i < 3; i++) {
            float Kd, Ks, shininess;
            getMaterialProperties(objectType, objectIndex, i, Kd, Ks, shininess);

            float Isd = light.material[i];
            summedColor[i] += Isd * Kd * std::max(0.f, Vec3::dot(L, N));
            summedColor[i] += Isd * Ks * pow(std::max(0.f, Vec3::dot(R, V)), shininess);
        }
    }


    Vec3 phong(RaySceneIntersection intersectionObjet) {
        Vec3 color(0., 0., 0.);
        Vec3 ambientColor(0., 0., 0.);
        Vec3 sommeSpecular(0., 0., 0.);
        Vec3 sommeDiffuse(0., 0., 0.);

        // Réflexion Ambiante
        float Ka = 0.35; // Coefficient de réflexion ambiante
        ambientColor = calculateAmbient(Ka);

        if (!intersectionObjet.intersectionExists) {
            return Vec3(0., 0., 0.); // Pas d'intersection
        }

        Vec3 P, L, N, V;

        for (int i_Light = 0; i_Light < lights.size(); i_Light++) {
            Light &light = lights[i_Light];

            // Récupération des informations en fonction du type d'objet intersecté
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

            } else if (intersectionObjet.typeOfIntersectedObject == MESH) {
                RayTriangleIntersection result_tmp = intersectionObjet.rayMeshIntersection;
                Mesh &mesh = meshes[intersectionObjet.objectIndex];
                P = result_tmp.intersection;
                N = result_tmp.normal;
                N.normalize();

            } else {
                continue; // Type d'objet non géré
            }

            if (light.type == LightType_Quad) {
                // Lumière étendue (Quad)
                Vec3 summedColor(0.0f, 0.0f, 0.0f);
                int numSamples = 16; // Échantillonnage Monte Carlo

                for (int i = 0; i < numSamples; ++i) {
                    Vec3 sampledLightPosition = samplePointOnLight(light);
                    L = (sampledLightPosition - P);
                    L.normalize();

                    V = P * -1;
                    V.normalize();

                    // Vérifier l'ombre pour ce point échantillonné
                    Ray shadowRay(P + N * 0.001, L);
                    RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);

                    if (shadowIntersection.intersectionExists &&
                        shadowIntersection.t < (sampledLightPosition - P).length()) {
                        continue; // Rayon bloqué par un objet
                    }

                    computeLightContribution(summedColor, light, P, N, V, L, 
                                            intersectionObjet.typeOfIntersectedObject, intersectionObjet.objectIndex);
                }

                summedColor /= numSamples; // Moyenne des échantillons

                // Ajouter la couleur calculée
                sommeDiffuse += summedColor;
                sommeSpecular += summedColor;

            } else {
                // Lumière ponctuelle
                L = (light.pos - P);
                L.normalize();

                V = P * -1;
                V.normalize();

                // Vérifier si le point est dans l'ombre
                Ray shadowRay(P + N * 0.001, L);
                RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);

                if (shadowIntersection.intersectionExists && 
                    shadowIntersection.t < (light.pos - P).length()) {
                    continue; // Rayon bloqué par un objet
                }

                computeLightContribution(sommeDiffuse, light, P, N, V, L, 
                                        intersectionObjet.typeOfIntersectedObject, intersectionObjet.objectIndex);
            }
        }

        // Finalisation de la couleur
        color = ambientColor + sommeDiffuse + sommeSpecular;

        return color;
    }




/* // Ancien phong
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

        // P : Point d'intersection entre le rayon et l'objet
        //(L.N) : angle entre la source de lumière et la normale
        //(R.V) : angle entre les directions de réflexion et de la vue
        Vec3 P, L, N, V, R;

        // Parcours des sources de lumière
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

            // Pour les lumières de type Quad
            if (light.type == LightType_Quad) {
                Vec3 summedColor(0.0f, 0.0f, 0.0f);
                int numSamples = 16; // Nombre d'échantillons pour la lumière étendue

                // Moyenne des échantillons pour simuler des ombres douces
                for (int i = 0; i < numSamples; ++i) {
                    // Échantillonnage d'un point sur la lumière
                    Vec3 sampledLightPosition = samplePointOnLight(light);

                    // Calcul de la direction de la lumière
                    L = (sampledLightPosition - P);
                    L.normalize();

                    // Calcul de la direction vers la caméra
                    V = P * -1;
                    V.normalize();

                    // Calcul de la direction réfléchie
                    R = (2. * Vec3::dot(N, L) * N - L);
                    R.normalize();

                    // Vérifier l'ombre pour ce point échantillonné
                    Ray shadowRay(P + N * 0.001, L); // Rayon légèrement décalé pour éviter l'autointersection
                    RaySceneIntersection shadowIntersection = computeIntersection(shadowRay);

                    if (shadowIntersection.intersectionExists &&
                        shadowIntersection.t < (sampledLightPosition - P).length()) {
                        // Si un objet bloque la lumière, ignorer cette contribution
                        continue;
                    }

                    // Ajouter la contribution diffuse et spéculaire pour ce point échantillonné
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

                        summedColor[i] += Isd * Kd * std::max(0.f, Vec3::dot(L, N));
                        summedColor[i] += Isd * Ks * pow(std::max(0.f, Vec3::dot(R, V)), shininess);
                    }
                }

                // Moyenne des contributions de la lumière étendue
                summedColor /= numSamples;

                // Ajouter la couleur calculée à la couleur finale
                for (int i = 0; i < 3; i++) {
                    sommeDifraction[i] += summedColor[i];
                    sommeSpecular[i] += summedColor[i];
                }

            } else {
                // Pour les lumières ponctuelles (Spherical)
                L = (light.pos - P);
                L.normalize();

                V = P * -1; // Direction vers la caméra
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
        }

        // Finalisation de la couleur
        for (int i = 0; i < 3; i++) {
            color[i] = ambientColor[i] + sommeDifraction[i] + sommeSpecular[i];
        }

        return color;
    }
*/

//----------------------------------------------------------------------------------------



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

        //PHASE 3 et 4
        if (NRemainingBounces > 0) {
            
            if (raySceneIntersection.typeOfIntersectedObject == SPHERE) {

                Sphere &sph = spheres[raySceneIntersection.objectIndex];

                Vec3 P = raySceneIntersection.raySphereIntersection.intersection; // Point d'intersection
                Vec3 N = raySceneIntersection.raySphereIntersection.normal;       // Normale à la surface
                Vec3 I = ray.direction();                                         // Direction du rayon incident

                // Si l'objet a un matériau miroir (réfléchissant)
                if (sph.material.type == Material_Mirror) {

                    Vec3 R = I - 2.0f * Vec3::dot(I, N) * N; // Direction du rayon réfléchi

                    Ray reflectedRay(P + N * 0.001f, R); // Déplace légèrement le point d'intersection pour éviter les auto-intersections

                    // Appelle récursivement la fonction rayTrace pour obtenir la couleur réfléchie
                    Vec3 reflectedColor = rayTraceRecursive(reflectedRay, NRemainingBounces - 1);

                    // Ajoute la contribution de la réflexion à la couleur finale
                    color[0] += reflectedColor[0] * sph.material.specular_material[0];
                    color[1] += reflectedColor[1] * sph.material.specular_material[1];
                    color[2] += reflectedColor[2] * sph.material.specular_material[2];

                    // Normalisation de la couleur pour éviter que les réflexions n'augmentent trop la luminosité
                    // Normalisation uniquement de la partie réflexion
                    color = Vec3::clamp(color, 0.0f, 1.0f);
                }

                // Si l'objet a un matériau en verre (transparent)
                if (sph.material.type == Material_Glass) {
                    Vec3 P = raySceneIntersection.raySphereIntersection.intersection; // Point d'intersection
                    Vec3 N = raySceneIntersection.raySphereIntersection.normal; // Normale
                    Vec3 I = ray.direction(); // Direction incidente

                    // Indices de réfraction
                    float n1 = 1.0f; // Air
                    float n2 = sph.material.index_medium; // Indice du verre

                    // Si le rayon est sortant, inverser les indices et la normale
                    if (Vec3::dot(I, N) > 0.0f) {
                        std::swap(n1, n2);
                        N = N * -1.0f; // Inverser manuellement la normale
                    }

                    // Calcul des coefficients de Fresnel
                    float cosI = fabs(Vec3::dot(I, N));
                    float sinT2 = (n1 / n2) * (n1 / n2) * (1.0f - cosI * cosI);
                    float fresnelReflectance;

                    if (sinT2 > 1.0f) {
                        // Réflexion totale interne
                        fresnelReflectance = 1.0f;
                    } else {
                        float cosT = sqrt(1.0f - sinT2);
                        float Rs = pow((n1 * cosI - n2 * cosT) / (n1 * cosI + n2 * cosT), 2);
                        float Rp = pow((n1 * cosT - n2 * cosI) / (n1 * cosT + n2 * cosI), 2);
                        fresnelReflectance = (Rs + Rp) / 2.0f;
                    }

                    // Rayon réfléchi
                    Vec3 R = I - 2.0f * Vec3::dot(I, N) * N; // Direction réfléchie
                    Ray reflectedRay(P + N * 0.001f, R); // Déplacement pour éviter les auto-intersections
                    Vec3 reflectedColor = rayTraceRecursive(reflectedRay, NRemainingBounces - 1);

                    // Rayon réfracté (si pas de réflexion totale interne)
                    Vec3 refractedColor(0.0f, 0.0f, 0.0f);
                    if (fresnelReflectance < 1.0f) {
                        Vec3 T = (I - N * cosI) * (n1 / n2) - N * sqrt(1.0f - sinT2); // Direction réfractée
                        Ray refractedRay(P - N * 0.001f, T); // Déplacement léger vers l'intérieur
                        refractedColor = rayTraceRecursive(refractedRay, NRemainingBounces - 1);
                    }

                    // Combiner les couleurs
                    color = fresnelReflectance * reflectedColor + (1.0f - fresnelReflectance) * refractedColor;
                    color = Vec3::clamp(color, 0.0f, 1.0f); // Éviter les dépassements
                }

            }
        }

        return color;
    }


    Vec3 rayTrace( Ray const & rayStart ) {
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

    void setup_cornell_box_meshes(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        // Lumière de type Quad
        {
            lights.resize(lights.size() + 1);
            Light &light = lights.back();
            light.pos = Vec3(0., 1.5, 0.);
            light.type = LightType_Quad;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;

            // Initialiser la géométrie du quad
            initialize_quad_light(light, 4.0f, 4.0f); // Taille 4x4
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

        // Initialisation du fichier OFF
        {
            meshes.resize(meshes.size() + 1);
            Mesh &mesh = meshes.back();

            // Chemin vers le fichier OFF
            std::string offFilePath = "img/star.off";

            // Chargement du fichier OFF
            mesh.loadOFF(offFilePath);

            // Recentre et met à l'échelle le maillage
            mesh.centerAndScaleToUnit();

            // Recalcule les normales des sommets
            mesh.recomputeNormals();

            // Hauteur minimale du maillage (axe Y)
            float minY = std::numeric_limits<float>::max();
            for (const auto &vertex : mesh.vertices) {
                if (vertex.position[1] < minY) {
                    minY = vertex.position[1];
                }
            }

            // Translation verticale pour aligner le maillage sur le sol
            Vec3 translation(0.0, -minY - 1.8, 0.0);
            for (auto &vertex : mesh.vertices) {
                vertex.position += translation;
            }
            
            std::cout << "Maillage chargé : " << offFilePath << std::endl;
        }

    }

    void setup_cornell_box(){
        meshes.clear();
        spheres.clear();
        squares.clear();
        lights.clear();

        // Lumière de type Quad
        {
            lights.resize(lights.size() + 1);
            Light &light = lights.back();
            light.pos = Vec3(0., 1.5, 0.);
            light.type = LightType_Quad;
            light.material = Vec3(1, 1, 1);
            light.isInCamSpace = false;

            // Initialise la géométrie du quad
            initialize_quad_light(light, 4.0f, 4.0f); // Taille 4x4
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
        
        { //MIRRORED Sphere

            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(1.0, -1.25, 0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Mirror;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3( 1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 1.0;
            s.material.index_medium = 1.4;
        }

        { //GLASS Sphere
            spheres.resize( spheres.size() + 1 );
            Sphere & s = spheres[spheres.size() - 1];
            s.m_center = Vec3(-1.0, -1.25, -0.5);
            s.m_radius = 0.75f;
            s.build_arrays();
            s.material.type = Material_Glass;
            s.material.diffuse_material = Vec3( 1.,1.,1. );
            s.material.specular_material = Vec3(  1.,1.,1. );
            s.material.shininess = 16;
            s.material.transparency = 0.9;
            s.material.index_medium = 1.5;
        }
        

    }

};



#endif
