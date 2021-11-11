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
    
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius*radius;

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
