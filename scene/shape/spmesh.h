#pragma once

#include "hemesh.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "Eigen/StdVector"
#include "Eigen/Geometry"
#include "Eigen/Dense"

typedef struct InHalfedge InHalfedge;
typedef struct InEdge InEdge;
typedef struct InFace InFace;

typedef struct InHalfedge {
    std::shared_ptr<InVertex> v;
    std::shared_ptr<InHalfedge> next;
    std::shared_ptr<InHalfedge> twin;
    std::shared_ptr<InEdge> edge;
    std::shared_ptr<InFace> face;
} InHalfedge;

typedef struct InVertex {
    std::shared_ptr<InHalfedge> halfedge;
    Eigen::Vector3f pos;
    std::shared_ptr<ExVertex> exVertex;
    std::shared_ptr<ExFace> exFace;
} InVertIn;

typedef struct InEdge {
    std::shared_ptr<InHalfedge> halfedge;
//    float error;
//    Eigen::Vector3f newPoint;
//    float length;
//    bool collapse;
} InEdge;

typedef struct InFace {
    std::shared_ptr<InHalfedge> halfedge;
} InFace;

class SPmesh
{
public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    void initFromVectors(const std::vector<Eigen::Vector3f> &vertices,
                         const std::vector<Eigen::Vector3i> &faces);
    void loadFromFile(const std::string &filePath);
    void loadHalfEdges();
    void validate();

private:
    void checkCircular(const std::shared_ptr<InHalfedge> &halfedge);
    void checkTwin(const std::shared_ptr<InHalfedge> &halfedge);
    void checkFaces();
    void checkVertices();
//    int getDegree(const std::shared_ptr<InVertex> &v);
//    Eigen::Vector3f getNormal(Eigen::Vector3f &v1, Eigen::Vector3f &v2, Eigen::Vector3f &v3);
//    float getArea(Eigen::Vector3f &v1, Eigen::Vector3f &v2, Eigen::Vector3f &v3);

    std::vector<Eigen::Vector3f> _vertices;
    std::vector<Eigen::Vector3i> _facesList;
    std::unordered_set<std::shared_ptr<InVertex>> _verts;
    std::unordered_set<std::shared_ptr<InEdge>> _edges;
    std::unordered_set<std::shared_ptr<InHalfedge>> _halfedges;
    std::unordered_set<std::shared_ptr<InFace>> _faces;

    std::unordered_set<std::shared_ptr<InVertex>> _newVerts;
    std::unordered_set<std::shared_ptr<InEdge>> _newEdges;
    std::unordered_set<std::shared_ptr<InHalfedge>> _newHalfedges;
    std::unordered_set<std::shared_ptr<InFace>> _newFaces;
    std::unordered_set<std::shared_ptr<InEdge>> _newMiddleEdges;
    HEmesh _exMesh;
};
