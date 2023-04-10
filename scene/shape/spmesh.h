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
    float angle; // angle relative to the reference direction of the vertex at the base of the half angle
} InHalfedge;

typedef struct InVertex {
    std::shared_ptr<InHalfedge> halfedge;
    Eigen::Vector3f pos;
    std::shared_ptr<ExVertex> exVertex; // (only set for original vertices of the input)
    std::shared_ptr<ExFace> exFace; // null (only set for new vertices without extrinsic counterpart)
    float angleSum; // (big theta in equation 1)
    Eigen::Vector3f barycentricPos; //(of this vertex in the extrinsic triangle it’s contained in)
} InVertIn;

typedef struct InEdge {
    std::shared_ptr<InHalfedge> halfedge;
    float length;
} InEdge;

typedef struct InFace {
    std::shared_ptr<InHalfedge> halfedge;
} InFace;



class SPmesh
{
public:
    // SPmesh();
    void initFromVectors(const std::vector<Eigen::Vector3f> &vertices,
                         const std::vector<Eigen::Vector3i> &faces);
    void loadHalfEdges();
    void validate();

private:
    // validator helpers
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
