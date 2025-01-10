#ifndef LINE_H
#define LINE_H

#include <cmath>
#include <iostream>
#include "Vec3.h"

class Line {
private:
    Vec3 m_origin , m_direction;
public:
    Line() {}
    Line( Vec3 const & o , Vec3 const & d )  {
        m_origin = o;
        m_direction = d; m_direction.normalize();
    }
    Vec3 & origin() { return m_origin; }
    Vec3 const & origin() const { return m_origin; }
    Vec3 & direction() { return m_direction; }
    Vec3 const & direction() const { return m_direction; }
    
    Vec3 project( Vec3 const & p ) const {
        Vec3 result;
        
        Vec3 OP = p - m_origin; // Vecteur entre l'origine de la ligne et le point
        float t = Vec3::dot(OP, m_direction); // Projection scalaire sur la direction
        
        result = m_origin + t * m_direction;   // Position du point projeté
        
        return result;
    }

    float squareDistance( Vec3 const & p ) const {
        float result;
        // Calcul du point projeté
        Vec3 projectedPoint = project(p);

        // Distance au carré entre le point d'origine et le point projeté
        result = (p - projectedPoint).length();
        
        return result;
    }

    float distance( Vec3 const & p ) const {
        return sqrt( squareDistance(p) );
    }
};

static inline std::ostream & operator << (std::ostream & s , Line const & l) {
    s << l.origin() << " " << l.direction();
    return s;
}

#endif
