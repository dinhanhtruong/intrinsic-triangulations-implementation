#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <BVH/Object.h>
#include <util/tiny_obj_loader.h>

class Triangle : public Object
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Triangle();
    Triangle(Eigen::Vector3f v1, Eigen::Vector3f v2, Eigen::Vector3f v3,
             Eigen::Vector3f n1, Eigen::Vector3f n2, Eigen::Vector3f n3,
             int index);

    bool getIntersection(const Ray &ray, IntersectionInfo *intersection) const override;

    Eigen::Vector3f getNormal(const IntersectionInfo &I) const override;
    virtual Eigen::Vector3f getNormal(const Eigen::Vector3f &p) const;

    BBox getBBox() const override;

    Eigen::Vector3f getCentroid() const override;

    int getIndex() const;

    tinyobj::material_t getMaterial() const;
    tinyobj::MaterialType getMaterialType() const;
    void setMaterial(const tinyobj::material_t &material);

    Eigen::Vector3<Eigen::Vector3f> getVertices() const { return Eigen::Vector3<Eigen::Vector3f>(_v1, _v2, _v3); }
    Eigen::Vector3<Eigen::Vector3f> getNormals()  { return Eigen::Vector3<Eigen::Vector3f>(_n1, _n2, _n3); }

    Eigen::Vector3f emittedRadiance(const Eigen::Vector3f &rayOut) const;
    Eigen::Vector3f brdf(const Eigen::Vector3f &w_i, const Eigen::Vector3f &w_o, const Eigen::Vector3f &normal) const;
    float getTriangleArea();
    Eigen::Vector3f sampleTrianglePoint();

private:
    Eigen::Vector3f _v1, _v2, _v3;
    Eigen::Vector3f _n1, _n2, _n3;

    tinyobj::material_t m_material;
    tinyobj::MaterialType m_materialType;

    int m_index;

    BBox _bbox;

    Eigen::Vector3f _centroid;

};

#endif // TRIANGLE_H
