#include "KDTree.h"
#include "Scene.h"
#include <algorithm>
#include <limits>
#include <cmath>

// ---------------- Node ----------------

Node::Node()
    : splittingAxis(-1),
      cuttingPos(0.f),
      leftChild(nullptr),
      rightChild(nullptr),
      isLeaf(false) {}

// ---------------- KDTree ----------------

KDTree::KDTree()
    : root(nullptr),
      maxDepth(15),
      maxPrimitivesPerLeaf(30),
      epsilon(0.01f)
{}

KDTree::~KDTree() {
    deleteNode(root);
    root = nullptr;
}

void KDTree::deleteNode(Node* node){
    if (!node) return;
    deleteNode(node->leftChild);
    deleteNode(node->rightChild);
    delete node;
}

void KDTree::setRoot(Node* _root){
    root = _root;
}

Node* KDTree::getRoot(){
    return root;
}

void KDTree::setMaxDepth(int _maxDepth){
    maxDepth = _maxDepth;
}

int KDTree::getMaxDepth(){
    return maxDepth;
}

void KDTree::setMaxPrimitivesPerLeaf(int _maxPrimitivesPerLeaf){
    maxPrimitivesPerLeaf = _maxPrimitivesPerLeaf;
}

int KDTree::getMaxPrimitivesPerLeaf(){
    return maxPrimitivesPerLeaf;
}

void KDTree::setEpsilon(float _epsilon){
    epsilon = _epsilon;
}

float KDTree::getEpsilon(){
    return epsilon;
}

/*
// buildNode : V1
Node* KDTree::buildNode(AABB &aabb, std::vector<PrimitiveRef> &prim, const Scene& scene, int depth){
    Node* n = new Node();
    n->boundingBox = aabb;

    //primitives ≤ seuil ?
    if(prim.size() <= maxPrimitivesPerLeaf){
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }
    //profondeur max atteinte ?
    if(depth >= maxDepth){
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    float sizeX = aabb.maxCoor[0] - aabb.minCoor[0];
    float sizeY = aabb.maxCoor[1] - aabb.minCoor[1];
    float sizeZ = aabb.maxCoor[2] - aabb.minCoor[2];
    float aabbSize = std::max({sizeX, sizeY, sizeZ});

    if(aabbSize < epsilon){
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    // Sinon : nœud interne
    n->isLeaf = false;

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
            
            //calculer pMin (sa borne minimale sur l’axe de coupe), pMax (sa borne maximale sur l’axe de coupe)
            //tester les trois cas
            //ajouter la primitive à la/les bonne(s) liste(s)
            

        // Calculs de pMin et pMax 
        float pMin, pMax;
        int axis = n->splittingAxis;

        switch (prim[i].type) {
            case PrimitiveType::Sphere: {
                    
                    //pMin = center[axis] - radius
                    //pMax = center[axis] + radius
                    
                const Sphere& s = spheres[prim[i].objectIndex];
                pMin = s.m_center[axis] - s.m_radius;
                pMax = s.m_center[axis] + s.m_radius;
                break;
            }
            case PrimitiveType::Triangle: {
                    
                    //prendre les 3 sommets
                    //pMin = min(v0, v1, v2)[axis]
                    //pMax = max(v0, v1, v2)[axis]
                    
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
                    
                    //utiliser son AABB pré-calculée
                    //prendre min/max sur l’axe
                    
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
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    n->leftChild  = buildNode(aabbLeftChild, leftPrimitives, scene, depth + 1);
    n->rightChild = buildNode(aabbRightChild, rightPrimitives, scene, depth + 1);

    return n;
}
*/

/*
// buildNode : V2
Node* KDTree::buildNode(AABB &aabb, std::vector<PrimitiveRef> &prim, const Scene& scene, int depth){
    Node* n = new Node();
    n->boundingBox = aabb;

    // Si trop peu de primitives ou profondeur max atteinte, on crée une feuille
    if (prim.size() <= maxPrimitivesPerLeaf || depth >= maxDepth) {
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    // Axe de coupe = axe le plus long
    float sizeX = aabb.maxCoor[0] - aabb.minCoor[0];
    float sizeY = aabb.maxCoor[1] - aabb.minCoor[1];
    float sizeZ = aabb.maxCoor[2] - aabb.minCoor[2];

    if (std::max({sizeX, sizeY, sizeZ}) < epsilon) { // feuille trop petite
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    int axis = 0;
    if (sizeY >= sizeX && sizeY >= sizeZ) axis = 1;
    else if (sizeZ >= sizeX && sizeZ >= sizeY) axis = 2;
    n->splittingAxis = axis;

    // ---------------- SAH-lite : coupe selon le centre de gravité des primitives ----------------
    float cutPos = 0.f;
    {
        // calculer le centre moyen des primitives sur l'axe
        for (const auto& p : prim) {
            switch (p.type) {
                case PrimitiveType::Sphere:
                    cutPos += scene.getSpheres()[p.objectIndex].m_center[axis];
                    break;
                case PrimitiveType::Triangle:
                    cutPos += (scene.getSquares()[p.objectIndex].vertices[0].position[axis] +
                               scene.getSquares()[p.objectIndex].vertices[1].position[axis] +
                               scene.getSquares()[p.objectIndex].vertices[2].position[axis]) / 3.0f;
                    break;
                case PrimitiveType::Mesh: {
                    const Mesh& m = scene.getMeshes()[p.objectIndex];
                    float sum = 0.f;
                    for (const auto& v : m.vertices) sum += v.position[axis];
                    cutPos += sum / m.vertices.size();
                    break;
                }
            }
        }
        cutPos /= prim.size(); // centre moyen
    }
    n->cuttingPos = cutPos;

    // ---------------- Répartition des primitives ----------------
    std::vector<PrimitiveRef> leftPrimitives, rightPrimitives;
    for (const auto& p : prim) {
        float pMin, pMax;

        switch (p.type) {
            case PrimitiveType::Sphere: {
                const Sphere& s = scene.getSpheres()[p.objectIndex];
                pMin = s.m_center[axis] - s.m_radius;
                pMax = s.m_center[axis] + s.m_radius;
                break;
            }
            case PrimitiveType::Triangle: {
                const Square& t = scene.getSquares()[p.objectIndex];
                pMin = std::min({ t.vertices[0].position[axis], t.vertices[1].position[axis], t.vertices[2].position[axis] });
                pMax = std::max({ t.vertices[0].position[axis], t.vertices[1].position[axis], t.vertices[2].position[axis] });
                break;
            }
            case PrimitiveType::Mesh: {
                const Mesh& m = scene.getMeshes()[p.objectIndex];
                pMin = std::numeric_limits<float>::infinity();
                pMax = -std::numeric_limits<float>::infinity();
                for (const auto& v : m.vertices) {
                    float val = v.position[axis];
                    if (val < pMin) pMin = val;
                    if (val > pMax) pMax = val;
                }
                break;
            }
        }

        // --------- assigner selon le centre, pas la duplication ---------
        float center = (pMin + pMax) / 2.0f;
        if (center <= cutPos) leftPrimitives.push_back(p);
        else rightPrimitives.push_back(p);
    }

    // Si partition dégénérée → feuille
    if (leftPrimitives.empty() || rightPrimitives.empty()) {
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    n->isLeaf = false;

    // Construction récursive
    AABB aabbLeft = aabb;
    aabbLeft.maxCoor[axis] = cutPos;
    AABB aabbRight = aabb;
    aabbRight.minCoor[axis] = cutPos;

    n->leftChild  = buildNode(aabbLeft, leftPrimitives, scene, depth + 1);
    n->rightChild = buildNode(aabbRight, rightPrimitives, scene, depth + 1);

    return n;
}
*/

// buildNode : V3
Node* KDTree::buildNode(AABB &aabb, std::vector<PrimitiveRef> &prim, const Scene& scene, int depth){
    Node* n = new Node();
    n->boundingBox = aabb;

    // Feuille si primitives peu nombreuses ou profondeur max atteinte
    if (prim.size() <= maxPrimitivesPerLeaf || depth >= maxDepth) {
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    // Axe le plus long
    float sizeX = aabb.maxCoor[0] - aabb.minCoor[0];
    float sizeY = aabb.maxCoor[1] - aabb.minCoor[1];
    float sizeZ = aabb.maxCoor[2] - aabb.minCoor[2];

    if (std::max({sizeX, sizeY, sizeZ}) < epsilon) { // nœud trop petit
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    int axis = 0;
    if (sizeY >= sizeX && sizeY >= sizeZ) axis = 1;
    else if (sizeZ >= sizeX && sizeZ >= sizeY) axis = 2;
    n->splittingAxis = axis;

    // ---------------- SAH simplifié : trouver centre moyen des AABB ----------------
    float cutPos = 0.f;
    for (const auto& p : prim) {
        AABB primAABB;
        switch (p.type) {
            case PrimitiveType::Sphere: {
                const Sphere& s = scene.getSpheres()[p.objectIndex];
                primAABB.minCoor = s.m_center - Vec3(s.m_radius, s.m_radius, s.m_radius);
                primAABB.maxCoor = s.m_center + Vec3(s.m_radius, s.m_radius, s.m_radius);
                break;
            }
            case PrimitiveType::Triangle: {
                const Square& t = scene.getSquares()[p.objectIndex];
                primAABB.minCoor = primAABB.maxCoor = t.vertices[0].position;
                for (int v=1; v<3; v++){
                    for (int d=0; d<3; d++){
                        if (t.vertices[v].position[d] < primAABB.minCoor[d]) primAABB.minCoor[d] = t.vertices[v].position[d];
                        if (t.vertices[v].position[d] > primAABB.maxCoor[d]) primAABB.maxCoor[d] = t.vertices[v].position[d];
                    }
                }
                break;
            }
            case PrimitiveType::Mesh: {
                primAABB = scene.getMeshes()[p.objectIndex].aabb; // AABB précalculée
                break;
            }
        }
        cutPos += (primAABB.minCoor[axis] + primAABB.maxCoor[axis]) * 0.5f; // centre de la primitive
    }
    cutPos /= prim.size();
    n->cuttingPos = cutPos;

    // ---------------- Partitionnement ----------------
    std::vector<PrimitiveRef> leftPrimitives, rightPrimitives;
    for (const auto& p : prim) {
        AABB primAABB;
        switch (p.type) {
            case PrimitiveType::Sphere: {
                const Sphere& s = scene.getSpheres()[p.objectIndex];
                primAABB.minCoor = s.m_center - Vec3(s.m_radius, s.m_radius, s.m_radius);
                primAABB.maxCoor = s.m_center + Vec3(s.m_radius, s.m_radius, s.m_radius);
                break;
            }
            case PrimitiveType::Triangle: {
                const Square& t = scene.getSquares()[p.objectIndex];
                primAABB.minCoor = primAABB.maxCoor = t.vertices[0].position;
                for (int v=1; v<3; v++){
                    for (int d=0; d<3; d++){
                        if (t.vertices[v].position[d] < primAABB.minCoor[d]) primAABB.minCoor[d] = t.vertices[v].position[d];
                        if (t.vertices[v].position[d] > primAABB.maxCoor[d]) primAABB.maxCoor[d] = t.vertices[v].position[d];
                    }
                }
                break;
            }
            case PrimitiveType::Mesh:
                primAABB = scene.getMeshes()[p.objectIndex].aabb;
                break;
        }

        // placer selon le centre de l'AABB
        float center = 0.5f * (primAABB.minCoor[axis] + primAABB.maxCoor[axis]);
        if (center <= cutPos) leftPrimitives.push_back(p);
        else rightPrimitives.push_back(p);
    }

    // Si partition dégénérée → feuille
    if (leftPrimitives.empty() || rightPrimitives.empty()) {
        n->isLeaf = true;
        n->primitives = prim;
        return n;
    }

    // Construction récursive des enfants
    n->isLeaf = false;
    AABB aabbLeft = aabb;
    aabbLeft.maxCoor[axis] = cutPos;
    AABB aabbRight = aabb;
    aabbRight.minCoor[axis] = cutPos;

    n->leftChild  = buildNode(aabbLeft, leftPrimitives, scene, depth+1);
    n->rightChild = buildNode(aabbRight, rightPrimitives, scene, depth+1);

    return n;
}



void KDTree::buildKDTree(const Scene& scene){
    deleteNode(root);
    root = nullptr;

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
    AABB globalAABB = scene.computeGlobalAABB();

    // Construction récursive du KDTree
    root = buildNode(globalAABB, allPrimitives, scene, 0);
}


bool KDTree::intersectPrimitive(const PrimitiveRef& prim, const Scene& scene, const Ray& ray, Hit& hit) const {
    
    switch (prim.type) {

        case PrimitiveType::Sphere: {
            const Sphere& s = scene.getSpheres()[prim.objectIndex];
            RaySphereIntersection res = s.intersect(ray);

            if (res.intersectionExists) {
                if (!hit.hit || res.t < hit.t) {
                    hit.hit = true;
                    hit.t = res.t;
                    hit.position = res.intersection;
                    hit.normal = res.normal;
                    hit.primitive = prim;
                }
                return true;
            }
            break;
        }

        case PrimitiveType::Triangle: {
            const Square& sq = scene.getSquares()[prim.objectIndex];
            RaySquareIntersection res = sq.intersect(ray);

            if (res.intersectionExists) {
                if (!hit.hit || res.t < hit.t) {
                    hit.hit = true;
                    hit.t = res.t;
                    hit.position = res.intersection;
                    hit.normal = res.normal;
                    hit.primitive = prim;
                }
                return true;
            }
            break;
        }

        case PrimitiveType::Mesh: {
            const Mesh& m = scene.getMeshes()[prim.objectIndex];
            RayTriangleIntersection res = m.intersect(ray);

            if (res.intersectionExists) {
                if (!hit.hit || res.t < hit.t) {
                    hit.hit = true;
                    hit.t = res.t;
                    hit.position = res.intersection;
                    hit.normal = res.normal;
                    hit.primitive = prim;
                }
                return true;
            }
            break;
        }
    }

    return false;
}
    

bool KDTree::intersectNode(Node* node, const Scene& scene, const Ray& ray, float t_entry, float t_exit, Hit& hit) const {
        
    int axis = node->splittingAxis;
    // Feuille
    if (node->isLeaf) {
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
    float t_split;
    if(fabs(ray.direction()[axis]) < epsilon){
        t_split = std::numeric_limits<float>::infinity();
    } else {
        t_split = (node->cuttingPos - ray.origin()[axis]) / ray.direction()[axis];
    }
    
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

bool KDTree::intersect(const Scene& scene, const Ray& ray, Hit& hit) const {
    float tEntry, tExit;
    if (!root) return false;

    if (!root->boundingBox.intersectRayAABB(ray, tEntry, tExit))
        return false;

    return intersectNode(root, scene, ray, tEntry, tExit, hit);
}
