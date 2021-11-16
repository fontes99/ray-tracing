#ifndef HIPERBOL_H
#define HIPERBOL_H

#include "hittable.h"
#include "vec3.h"

class hiperbol : public hittable {
    public:
        hiperbol() {}
        hiperbol(point3 cen, double r, double h, int o, shared_ptr<material> m) // parametrização
            : center(cen), radius(r), height(h), orientation(o), mat_ptr(m) {};

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;

    public: // parametrização
        point3 center;
        double radius;
        double height;
        int orientation; // 0 = x, 1 = y, 2 = z
        shared_ptr<material> mat_ptr;
};

bool hiperbol::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    
 /*

    x**2      +   y**2      -   z**2       =   R**2
   (x-Cx)**2  +  (y-Cy)**2  -  (z-Cz)**2   =   R**2
   
    P = A + t*d
    (x, y, z) = (Ax, Ay, Az) + t * (dx, dy, dz)
    (x, y, z) = (Ax + t*dx , Ay + t*dy , Az + t*dz)

    (Ax + t*dx - Cx)**2  +  (Ay + t*dy - Cy)**2  -  (Az + t*dz - Cz)**2   =   R**2

    e = A - C
    (ex, ey, ez) = Ax-Cx, Ay-Cy, Az-Cz, 

    (t*dx + ex)**2  +  (t*dy + ey)**2  -  (t*dz + ez)**2  -  R**2  = 0
    t**2 * dx**2 + 2*t*dx*ex + ex**2   +   t**2 * dy**2 + 2*t*dy*ey + ey**2   -   t**2 * dz**2 - 2*t*dz*ez - ez**2)  -  R**2  =  0

    (dx**2 + dy**2 - dz**2) * t**2   +   2*(dx*ex + dy*ey - dz*ez)*t  +  (ex**2 + ey**2 - ez**2) - R**2 = 0
 
 */

    //     e   =      A     -    C
    vec3 vec_e = r.origin() - center;

    auto a = r.direction().length_squared();
    auto half_b = dot(vec_e, r.direction());
    auto c = vec_e.length_squared() - radius*radius;

    auto delta = half_b*half_b - a*c; // Delta Bhaskara

    if (delta < 0) return false;
    auto sqrtd = sqrt(delta);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.p = r.at(rec.t);
    vec3 outward_normal = (rec.p - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;   

    return true;
}

#endif
