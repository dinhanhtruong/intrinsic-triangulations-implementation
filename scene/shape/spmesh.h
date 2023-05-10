#pragma once

#include "hemesh.h"
#include "triangle.h"

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <QImage>

#include "Eigen/StdVector"
#include "Eigen/Geometry"
#include "Eigen/Dense"

const std::vector<Eigen::Vector3d> colors = {
    Eigen::Vector3d(0.863, 0.314, 0.306),
    Eigen::Vector3d(1, 0.667, 0),
    Eigen::Vector3d(1, 0.957, 0.373),
    Eigen::Vector3d(0.557, 0.8, 0.106),
    Eigen::Vector3d(0, 0.659, 0.839),
    Eigen::Vector3d(0.792, 0.635, 1)
};

class Scene;
class PathTracer;

typedef struct InHalfedge InHalfedge;
typedef struct InEdge InEdge;
typedef struct InFace InFace;

typedef struct InHalfedge {
    std::shared_ptr<InVertex> v;
    std::shared_ptr<InHalfedge> next;
    std::shared_ptr<InHalfedge> twin;
    std::shared_ptr<InEdge> edge;
    std::shared_ptr<InFace> face;
    double angle; // angle relative to the reference direction of the vertex at the base of the half angle
//    double length; // should be the same as InEdge->length
} InHalfedge;

typedef struct InVertex {
    std::shared_ptr<InHalfedge> halfedge;
    std::shared_ptr<ExVertex> exVertex; // (only set for original vertices of the input)
    std::shared_ptr<ExFace> exFace; // null (only set for new vertices without extrinsic counterpart)
    double bigTheta; // (big theta in equation 1)
    Eigen::Vector3d barycentricPos; //(of this vertex in the extrinsic triangle it’s contained in)
} InVertIn;

typedef struct InEdge {
    std::shared_ptr<InHalfedge> halfedge;
    double length;
} InEdge;

typedef struct InFace {
    std::shared_ptr<InHalfedge> halfedge;
} InFace;

class SPmesh
{
public:
    void initFromVectors(const std::vector<Eigen::Vector3f> &vertices,
                         const std::vector<Eigen::Vector3i> &faces);
    void loadHalfEdges();
    void initSignpost();
    void validate();
    void validateSignpost();

    // algo 11
    std::tuple<std::shared_ptr<InFace>, Eigen::Vector3d> pointQuery(std::shared_ptr<ExFace> xyz, Eigen::Vector3d& p);

    // visualization functions
    void setRenderInfo(Scene* scene, PathTracer* tracer, QImage* image, QString outDir);
    void assignColors(std::unordered_set<std::shared_ptr<InFace>>& faces);
    void renderImage(int frame);
    Eigen::Vector3d getColor(const Triangle* tri, Eigen::Vector3d point, const Eigen::Vector3d &camPos);
    void renderFlipping();

private:
    // validator helpers
    void checkCircular(const std::shared_ptr<InHalfedge> &halfedge);
    void checkTwin(const std::shared_ptr<InHalfedge> &halfedge);
    void checkFaces();
    void checkVertices();
    void checkTriangleInequality(const std::shared_ptr<InFace> face);

    // math heleprs
    double getAngle(Eigen::Vector3d v1, Eigen::Vector3d v2);
    double getAngleFromEdgeLengths(double l_ij, double l_jk, double l_ki);
    double baseLength(double a, double b, double theta);
    double angleBetween(double a, double b);
    bool isEqual(double a, double b, double epsilon=0.00001);
    double argument(Eigen::Vector2d u, Eigen::Vector2d v);
    double getArea(std::shared_ptr<InFace> face);
    Eigen::Vector3d getBaryCoords(Eigen::Vector3d &p, Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3);
    double distanceToEdge(Eigen::Vector3d &p, Eigen::Vector3d &v1, Eigen::Vector3d &v2, double l_ij, double l_jk, double l_ki);

    // mesh helpers
    int getDegree(const std::shared_ptr<InVertex> &v);
    Eigen::Vector3d getVPos(std::shared_ptr<InVertex> v);
    std::shared_ptr<InEdge> getEdge(std::shared_ptr<InVertex> v0, std::shared_ptr<InVertex> v1) const;
    std::shared_ptr<InHalfedge> getHalfEdgeWithSource(std::shared_ptr<InEdge> edge, std::shared_ptr<InVertex> sourceVertex) const;
    void eraseTriangle(std::shared_ptr<InFace> tri);
    std::shared_ptr<InFace> insertTriangle(std::shared_ptr<InVertex> v0, std::shared_ptr<InVertex> v1, std::shared_ptr<InVertex> v2);
    std::shared_ptr<InVertex> insertCircumcenter(std::shared_ptr<InFace> face);

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
    ///         direction vector u of the trace in the local 2D coordinate system of the destination EXTRINSIC triangle.
    std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3d, Eigen::Vector2d> traceFromIntrinsicVertex(std::shared_ptr<InVertex> v_i, double distance, double angle);
    std::tuple<std::shared_ptr<InFace>, Eigen::Vector3d> traceFromExtrinsicVertex(std::shared_ptr<ExVertex> v_i, double distance, double angle);
    template <typename T>
    std::tuple<std::shared_ptr<T>, Eigen::Vector3d, Eigen::Vector2d> traceVector(auto base, Eigen::Vector3d baryCoords, double distance, double angle);
    void updateVertex(std::shared_ptr<InVertex> i);
    std::shared_ptr<InEdge> flipEdge(std::shared_ptr<InEdge> ij);
    double distance(double l_12, double l_23, double l_31, const Eigen::Vector3d p, const Eigen::Vector3d q);
    std::shared_ptr<InVertex> insertVertex(std::shared_ptr<InFace> face, Eigen::Vector3d barycentricCoords);
    std::pair<double, double> vectorToPoint(double l_ij, double l_jk, double l_ki, const Eigen::Vector3d &i, const Eigen::Vector3d &j, const Eigen::Vector3d &p);

    // triangulation
    bool edgeIsDelaunay(std::shared_ptr<InEdge> edge);
    bool shouldRefine(std::shared_ptr<InFace> tri, double minAngle);
    void flipToDelaunay(std::unordered_set<std::shared_ptr<InEdge>>& edgesToCheck, double minAngle, int maxFlips);
    void delaunayRefine(double minAngle, int maxInsertions);

    std::vector<Eigen::Vector3d> _vertices;
    std::vector<Eigen::Vector3i> _facesList;
    std::unordered_set<std::shared_ptr<InVertex>> _verts;
    std::unordered_set<std::shared_ptr<InEdge>> _edges;
    std::unordered_set<std::shared_ptr<InHalfedge>> _halfedges;
    std::unordered_set<std::shared_ptr<InFace>> _faces;
    // cache (partial) 1-ring neighborhood during erase/insertTriangles to determine edge existence without re-traversal
    // map from ordered vertex pairs (i,j) where i <= j to the InEdge connecting vertex i and j, if one exists.
    std::unordered_map<
        std::shared_ptr<InVertex>,
        std::unordered_map<std::shared_ptr<InVertex>, std::shared_ptr<InEdge>>
    > _vertPairToEdge;

    HEmesh _exMesh;

    std::unordered_map<std::shared_ptr<InFace>, int> _faceColors;
    Scene* _scene;
    PathTracer* _tracer;
    QImage* _image;
    QRgb* _data;
    QString _outDir;
    double _meanIntrinsicEdgeLength = 0;
    void computeMeanIntrinsicEdgeLength();
};

