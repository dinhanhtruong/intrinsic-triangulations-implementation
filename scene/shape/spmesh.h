#pragma once

#include "hemesh.h"
#include "triangle.h"

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
//    float length; // should be the same as InEdge->length
} InHalfedge;

typedef struct InVertex {
    std::shared_ptr<InHalfedge> halfedge;
    std::shared_ptr<ExVertex> exVertex; // (only set for original vertices of the input)
    std::shared_ptr<ExFace> exFace; // null (only set for new vertices without extrinsic counterpart)
    float bigTheta; // (big theta in equation 1)
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
    void initSignpost();
    void validate();

    // algo 11
    std::tuple<std::shared_ptr<InFace>, Eigen::Vector3f> pointQuery(std::shared_ptr<ExFace> xyz, Eigen::Vector3f& p);

    // visualization functions
    void assignColors();
    int getColor(const Triangle* tri, Eigen::Vector3f point);

private:
    /// validator helpers
    void checkCircular(const std::shared_ptr<InHalfedge> &halfedge);
    void checkTwin(const std::shared_ptr<InHalfedge> &halfedge);
    void checkFaces();
    void checkVertices();

    /// signpost algos & helpers
    // old code
//    int getDegree(const std::shared_ptr<InVertex> &v);
//    Eigen::Vector3f getNormal(Eigen::Vector3f &v1, Eigen::Vector3f &v2, Eigen::Vector3f &v3);
//    float getArea(Eigen::Vector3f &v1, Eigen::Vector3f &v2, Eigen::Vector3f &v3);

    // math/connectivity heleprs
    float getAngle(Eigen::Vector3f v1, Eigen::Vector3f v2);
    float getAngleFromEdgeLengths(float l_ij, float l_jk, float l_ki);
    float baseLength(float a, float b, float theta);
    float angleBetween(float a, float b);
    float argument(Eigen::Vector2f u, Eigen::Vector2f v);
    Eigen::Vector3f getBaryCoords(Eigen::Vector3f &p, Eigen::Vector3f &v1, Eigen::Vector3f &v2, Eigen::Vector3f &v3);
    Eigen::Vector3f getVPos(std::shared_ptr<InVertex> v);
    std::shared_ptr<InEdge> getEdge(std::shared_ptr<InVertex> v0, std::shared_ptr<InVertex> v1) const;
    std::shared_ptr<InHalfedge> getHalfEdgeWithSource(std::shared_ptr<InEdge> edge, std::shared_ptr<InVertex> sourceVertex) const;
    void eraseTriangle(std::shared_ptr<InFace> tri);
    std::shared_ptr<InFace> insertTriangle(std::shared_ptr<InVertex> v0, std::shared_ptr<InVertex> v1, std::shared_ptr<InVertex> v2);

    // algos (in order)
    void updateSignpost(std::shared_ptr<InHalfedge> h_ij);
    /// (note from Anh): my part needs traceFromVertex to return an additional output not explicitly specified in the pseudocode.
    /// The inputs are:
    ///     (1) starting point of the trace query at an intrinsic vertex v_i
    ///     (2) distance to trace
    ///     (3) direction to trace along as an angle in [0, 2*pi) relative to the reference direction at v_i
    /// The outputs in order should be:
    ///     (1) pointer to the extrinsic face containing the end point of the trace,
    ///     (2) barycentric coords of the end point
    ///     (3)** (this is the extra one my function needs):
    ///         direction vector u of the trace in the local 2D coordinate system of the destination EXTRINSIC triangle. First two
    std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3f, Eigen::Vector2f> traceFromIntrinsicVertex(std::shared_ptr<InVertex> v_i, float distance, float angle);
    std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3f, Eigen::Vector2f> traceVectorExtrinsic(std::shared_ptr<ExHalfedge> base, Eigen::Vector3f baryCoords, float distance, float angle);
    std::tuple<std::shared_ptr<InFace>, Eigen::Vector3f> traceFromExtrinsicVertex(std::shared_ptr<ExVertex> v_i, float distance, float angle);
    std::tuple<std::shared_ptr<InFace>, Eigen::Vector3f> traceVectorIntrinsic(std::shared_ptr<InHalfedge> base, Eigen::Vector3f baryCoords, float distance, float angle);
    void updateVertex(std::shared_ptr<InVertex> i);
    void flipEdge(std::shared_ptr<InEdge> ij);
    float distance(float l_12, float l_23, float l_31, const Eigen::Vector3f p, const Eigen::Vector3f q);
    void insertVertex(std::shared_ptr<InFace> face, Eigen::Vector3f& barycentricCoords);
    std::pair<float, float> vectorToPoint(float l_ij, float l_jk, float l_ki, const Eigen::Vector3f &i, const Eigen::Vector3f &j, const Eigen::Vector3f &p);
    void moveVertex(std::shared_ptr<InVertex> i, std::shared_ptr<InFace> iab, const Eigen::Vector3f &p);


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
    std::unordered_map<std::shared_ptr<InFace>, int> _faceColors;
};
