#pragma once

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>


#include "Eigen/StdVector"
#include "Eigen/Geometry"
#include "Eigen/Dense"

//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2d);
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3d);
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3i);

typedef struct ExHalfedge ExHalfedge;
typedef struct ExVertex ExVertex;
typedef struct ExEdge ExEdge;
typedef struct ExFace ExFace;
typedef struct InVertex InVertex;

typedef struct ExHalfedge {
    std::shared_ptr<ExVertex> v;
    std::shared_ptr<ExHalfedge> next;
    std::shared_ptr<ExHalfedge> twin;
    std::shared_ptr<ExEdge> edge;
    std::shared_ptr<ExFace> face;
    double angle;
} ExHalfedge;

typedef struct ExVertex {
    std::shared_ptr<ExHalfedge> halfedge;
    Eigen::Vector3d pos;
    std::shared_ptr<InVertex> inVertex;
    double bigTheta;
} ExVertex;

typedef struct ExEdge {
    std::shared_ptr<ExHalfedge> halfedge;
    double length;
} ExEdge;

typedef struct ExFace {
    std::shared_ptr<ExHalfedge> halfedge;
} ExFace;

class HEmesh
{
public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    void init(int numV, int numF);
    void validate();
    std::shared_ptr<ExVertex> makeVertex(std::shared_ptr<ExHalfedge> halfedge, Eigen::Vector3d pos, std::shared_ptr<InVertex> inVertex);
    std::shared_ptr<ExHalfedge> makeHalfedge(std::shared_ptr<ExVertex> v, std::shared_ptr<ExHalfedge> next, std::shared_ptr<ExHalfedge> twin, std::shared_ptr<ExEdge> edge, std::shared_ptr<ExFace> face);
    std::shared_ptr<ExEdge> makeEdge(std::shared_ptr<ExHalfedge> halfedge);
    std::shared_ptr<ExFace> makeFace(std::shared_ptr<ExHalfedge> halfedge);
    std::shared_ptr<ExFace> getExTriangle(int index);

private:
    void checkCircular(const std::shared_ptr<ExHalfedge> &halfedge);
    void checkTwin(const std::shared_ptr<ExHalfedge> &halfedge);
    void checkFaces();
    void checkVertices();
//    int getDegree(const std::shared_ptr<ExVertex> &v);
//    Eigen::Vector3d getNormal(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3);
//    double getArea(Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3);

    std::unordered_set<std::shared_ptr<ExVertex>> _verts;
    std::unordered_set<std::shared_ptr<ExEdge>> _edges;
    std::unordered_set<std::shared_ptr<ExHalfedge>> _halfedges;
    std::vector<std::shared_ptr<ExFace>> _faces;

    std::unordered_map<std::shared_ptr<ExFace>, int> _faceColors;
};
