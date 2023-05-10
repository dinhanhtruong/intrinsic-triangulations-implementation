#include "spmesh.h"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>

#include <QFileInfo>
#include <QString>
#include <assert.h>

#define TINYOBJLOADER_IMPLEMENTATION
#include "util/tiny_obj_loader.h"


using namespace Eigen;
using namespace std;

void SPmesh::initFromVectors(const vector<Vector3f> &vertices,
                           const vector<Vector3i> &faces)
{
    // Copy vertices and faces into internal vector
    for (Vector3f v: vertices) {
        _vertices.push_back(v.cast<double>());
    }
    _facesList    = faces;

    _verts = unordered_set<shared_ptr<InVertex>>();
    _edges = unordered_set<shared_ptr<InEdge>>();
    _faces = unordered_set<shared_ptr<InFace>>();
    _halfedges = unordered_set<shared_ptr<InHalfedge>>();

    _verts.reserve(vertices.size());
    _faces.reserve(faces.size());
    _edges.reserve(faces.size() * 1.5);
    _halfedges.reserve(faces.size() * 3);

    _exMesh.init(vertices.size(), faces.size());

    loadHalfEdges();
    _exMesh.validate();
    cout << "extrinsic mesh validated" << endl;

    initSignpost();
    validate();
    cout << "intrinsic mesh validated" << endl;
}

void SPmesh::loadHalfEdges() {
    cout << "vertices " << _vertices.size() << endl;
    unordered_map<int, shared_ptr<InVertex>> inVerts;
    for (int i = 0; i < _vertices.size(); i++) {
        shared_ptr<InVertex> v = make_shared<InVertex>(InVertex{nullptr, nullptr});
        shared_ptr<ExVertex> exv = _exMesh.makeVertex(nullptr, _vertices[i], v);
        v->exVertex = exv;
        inVerts.insert({i, v});
        _verts.insert(v);
    }
    unordered_map<shared_ptr<InVertex>, unordered_set<shared_ptr<InHalfedge>>> unmatched;
    unordered_map<shared_ptr<InHalfedge>, shared_ptr<ExHalfedge>> inToExHalfedge;
    for (int i = 0; i < _facesList.size(); i++) {
        Vector3i verts = _facesList[i];
        shared_ptr<InVertex> v1 = inVerts.at(verts[0]);
        shared_ptr<InVertex> v2 = inVerts.at(verts[1]);
        shared_ptr<InVertex> v3 = inVerts.at(verts[2]);
        shared_ptr<InHalfedge> h3, h2, h1;
        shared_ptr<ExHalfedge> exh3, exh2, exh1;

        // h3
        shared_ptr<InHalfedge> twin = nullptr;
        if (unmatched.find(v1) != unmatched.end()) {
            for (const shared_ptr<InHalfedge> &halfedge: unmatched.at(v1)) {
                if (halfedge->next != nullptr && halfedge->next->v == v3) {
                    twin = halfedge;
                }
            }
        }
        if (twin != nullptr) {
            h3 = make_shared<InHalfedge>(InHalfedge{v3, nullptr, twin, twin->edge, nullptr});
            shared_ptr<ExHalfedge> exTwin = inToExHalfedge.at(twin);
            exh3 = _exMesh.makeHalfedge(v3->exVertex, nullptr, exTwin, exTwin->edge, nullptr);
            twin->twin = h3;
            exTwin->twin = exh3;
            unmatched.at(v1).erase(twin);
        } else {
            h3 = make_shared<InHalfedge>(InHalfedge{v3, nullptr, nullptr, nullptr, nullptr});
            exh3 = _exMesh.makeHalfedge(v3->exVertex, nullptr, nullptr, nullptr, nullptr);
            shared_ptr<InEdge> edge = make_shared<InEdge>(InEdge{h3});
            shared_ptr<ExEdge> exEdge = _exMesh.makeEdge(exh3);
            h3->edge = edge;
            exh3->edge = exEdge;
            _edges.insert(edge);
            unmatched[v3].insert(h3);
        }
        v3->halfedge = h3;
        v3->exVertex->halfedge = exh3;
        shared_ptr<InFace> face = make_shared<InFace>(InFace{h3});
        shared_ptr<ExFace> exFace = _exMesh.makeFace(exh3);
        _faces.insert(face);

        // h2
        twin = nullptr;
        if (unmatched.find(v3) != unmatched.end()) {
            for (const shared_ptr<InHalfedge> &halfedge: unmatched.at(v3)) {
                if (halfedge->next != nullptr && halfedge->next->v == v2) {
                    twin = halfedge;
                }
            }
        }
        if (twin != nullptr) {
            h2 = make_shared<InHalfedge>(InHalfedge{v2, h3, twin, twin->edge, face});
            shared_ptr<ExHalfedge> exTwin = inToExHalfedge.at(twin);
            exh2 = _exMesh.makeHalfedge(v2->exVertex, exh3, exTwin, exTwin->edge, exFace);
            twin->twin = h2;
            exTwin->twin = exh2;
            unmatched.at(v3).erase(twin);
        } else {
            h2 = make_shared<InHalfedge>(InHalfedge{v2, h3, nullptr, nullptr, face});
            exh2 = _exMesh.makeHalfedge(v2->exVertex, exh3, nullptr, nullptr, exFace);
            shared_ptr<InEdge> edge = make_shared<InEdge>(InEdge{h2});
            shared_ptr<ExEdge> exEdge = _exMesh.makeEdge(exh2);
            h2->edge = edge;
            exh2->edge = exEdge;
            _edges.insert(edge);
            unmatched[v2].insert(h2);
        }
        v2->halfedge = h2;
        v2->exVertex->halfedge = exh2;

        // h1
        twin = nullptr;
        if (unmatched.find(v2) != unmatched.end()) {
            for (const shared_ptr<InHalfedge> &halfedge: unmatched.at(v2)) {
                if (halfedge->next != nullptr && halfedge->next->v == v1) {
                    twin = halfedge;
                }
            }
        }
        if (twin != nullptr) {
            h1 = make_shared<InHalfedge>(InHalfedge{v1, h2, twin, twin->edge, face});
            shared_ptr<ExHalfedge> exTwin = inToExHalfedge.at(twin);
            exh1 = _exMesh.makeHalfedge(v1->exVertex, exh2, exTwin, exTwin->edge, exFace);
            twin->twin = h1;
            exTwin->twin = exh1;
            unmatched.at(v2).erase(twin);
        } else {
            h1 = make_shared<InHalfedge>(InHalfedge{v1, h2, nullptr, nullptr, face});
            exh1 = _exMesh.makeHalfedge(v1->exVertex, exh2, nullptr, nullptr, exFace);
            shared_ptr<InEdge> edge = make_shared<InEdge>(InEdge{h1});
            shared_ptr<ExEdge> exEdge = _exMesh.makeEdge(exh1);
            h1->edge = edge;
            exh1->edge = exEdge;
            _edges.insert(edge);
            unmatched[v1].insert(h1);
        }
        v1->halfedge = h1;
        v1->exVertex->halfedge = exh1;
        // updates
        h3->next = h1;
        exh3->next = exh1;
        h3->face = face;
        exh3->face = exFace;
        _halfedges.insert(h1);
        _halfedges.insert(h2);
        _halfedges.insert(h3);

        inToExHalfedge.insert({h3, exh3});
        inToExHalfedge.insert({h2, exh2});
        inToExHalfedge.insert({h1, exh1});
    }
    cout << _edges.size() << endl;
}


// algo 4
void SPmesh::initSignpost() {
    for (shared_ptr<InEdge> e : _edges) {
        // get vertex positions of endpoints
        shared_ptr<InVertex> v0 = e->halfedge->v;
        shared_ptr<InVertex> v1 = e->halfedge->twin->v;
        e->length = (getVPos(v1) - getVPos(v0)).norm();
    }
    
    for (shared_ptr<InVertex> v: _verts) {
        double bigTheta_i = 0.0;
        shared_ptr<InHalfedge> curr = v->halfedge->next->next->twin;
        Vector3d rightEdge = getVPos(v->halfedge->twin->v) - getVPos(v);
        do {
            Vector3d leftEdge = getVPos(curr->twin->v) - getVPos(curr->v);
            double angle = getAngle(rightEdge, leftEdge);
            bigTheta_i += angle;
            rightEdge = leftEdge;
            curr = curr->next->next->twin;
        } while (curr != v->halfedge->next->next->twin);
        v->bigTheta = bigTheta_i;

        double phi_ij0 = 0.0;
        curr = v->halfedge;
        curr->angle = 0.0;
        while (curr->next->next->twin != v->halfedge) {
            updateSignpost(curr);
            curr = curr->next->next->twin;
        } 
    }

    for (shared_ptr<InVertex> v: _verts) {
        v->exVertex->bigTheta = v->bigTheta;
        shared_ptr<InHalfedge> curr = v->halfedge;
        shared_ptr<ExHalfedge> exCurr = v->exVertex->halfedge;
        do {
            exCurr->angle = curr->angle;
            exCurr->edge->length = curr->edge->length;
            curr = curr->next->next->twin;
            exCurr = exCurr->next->next->twin;
        } while (curr != v->halfedge);
    }

    cout << "initialized signpost values" << endl;
    validate();
}


// algo 1: updates the signpost angle of halfedge ik in triangle ijk
// ij must already have a valid angle (bc angle of ik depends on angle of ij). Edge lengths of ijk must also be valid
void SPmesh::updateSignpost(shared_ptr<InHalfedge> h_ij) {
    // given the 'right' halfedge ij of the face originiating from vertex i, get the left halfedge ik
    shared_ptr<InHalfedge> h_ik = h_ij->next->next->twin;
    // get edge lengths of triangle
    double l_ij = h_ij->edge->length;
    double l_jk = h_ij->next->edge->length;
    double l_ki = h_ij->next->next->edge->length;
    // update angle (phi) of halfedge ik
    double theta_i_jk = getAngleFromEdgeLengths(l_ij, l_jk, l_ki); // euclidean angle of opposite edge ik relative to ij (i.e. interior angle at vertex i between ij and ik)
    assert(theta_i_jk < M_PI); // else triangle angles sum to >= pi
    double phi_ik = h_ij->angle +  ((2*M_PI) * theta_i_jk)/h_ij->v->bigTheta; // == phi_ij + offset
    assert(h_ij->angle >= 0);
    assert(((2*M_PI) * theta_i_jk)/h_ij->v->bigTheta > 0);

    h_ik->angle = (phi_ik >= 2*M_PI) ? phi_ik - 2*M_PI : phi_ik; // constrain to [0, 2pi)
}

// angle should be the flat intrinsic angle (NOT phi)
std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3d, Eigen::Vector2d> SPmesh::traceFromIntrinsicVertex(std::shared_ptr<InVertex> v_i, double distance, double traceAngle) {
    shared_ptr<ExHalfedge> base;
    Vector3d baryCoords;
    double traceAngleRelativeTheta;
    if (v_i->exVertex == nullptr) {
        baryCoords = v_i->barycentricPos;
        base = v_i->exFace->halfedge;
        traceAngleRelativeTheta = traceAngle;
    } else {
        baryCoords = Vector3d(1.0, 0.0, 0.0);
        base = v_i->exVertex->halfedge;
        shared_ptr<ExHalfedge> right;
        shared_ptr<ExHalfedge> left;
        // in the general case, loop will terminate when the trace angle is between the phi angles of the left and right halfedges of the curr triangle
        // i.e. when right <= trace angle < left
        do {
            right = base; // right halfedge in perspective of v_i
            left = right->next->next->twin; // one step ccw

            // edge case: reference DIRECTION (0 degrees) is inside the current triangle => left halfedge's phi < right halfedge's phi due to wrap-around
            if (left->angle < right->angle) {
                // target/trace direction is contained in this triangle iff traceAngle <= left XOR > right. Cannot simultaneously satisfy (right < traceAngle <= left) in this case.
                if (traceAngle >= right->angle ||  traceAngle < left->angle) {
                    break;
                }
            }

            // general case: left phi > right phi in the triangle
            if (traceAngle < right->angle) {
                // move clockwise to adjacent face
                base = base->twin->next;
            } else if (traceAngle >= left->angle) {
                // move ccw
                base = base->next->next->twin;
            }
        } while (traceAngle < right->angle || traceAngle >=  left->angle);

        // unproject relative angle from [0, 2pi) -> [0, bigTheta)
        traceAngleRelativeTheta = angleBetween(traceAngle, base->angle) * (v_i->exVertex->bigTheta/(2*M_PI));
        assert(traceAngleRelativeTheta >= 0);
    }
    return traceVector<ExFace>(base, baryCoords, distance, traceAngleRelativeTheta);
}

// algo 2: see note in header file
// the angle argument is the normalized/projected phi angle in the range [0, 2pi] relative to the reference dir of v_i
std::tuple<std::shared_ptr<InFace>, Eigen::Vector3d> SPmesh::traceFromExtrinsicVertex(std::shared_ptr<ExVertex> v_i, double distance, double traceAngle) {
    shared_ptr<InHalfedge> base = v_i->inVertex->halfedge;

    // find the base vector whose source is at b_i and whose intrinsic triangle contains the dir to trace along (defined by the angle argument)
    // also find the flat/unprojected trace angle relative to base
    shared_ptr<InHalfedge> right;
    shared_ptr<InHalfedge> left;
    // in the general case, loop will terminate when the trace angle is between the phi angles of the left and right halfedges of the curr triangle
    // i.e. when right <= trace angle < left
    do {
        right = base; // right halfedge in perspective of v_i
        left = right->next->next->twin; // one step ccw

        // edge case: reference DIRECTION (0 degrees) is inside the current triangle => left halfedge's phi < right halfedge's phi due to wrap-around
        if (left->angle < right->angle) {
            // target/trace direction is contained in this triangle iff traceAngle <= left XOR > right. Cannot simultaneously satisfy (right < traceAngle <= left) in this case.
            if (traceAngle >= right->angle ||  traceAngle < left->angle) {
                break;
            }
        }

        // general case: left phi > right phi in the triangle
        if (traceAngle < right->angle) {
            // move clockwise to adjacent face
            base = base->twin->next;
        } else if (traceAngle >= left->angle) {
            // move ccw
            base = base->next->next->twin;
        }
    } while (traceAngle < right->angle || traceAngle >=  left->angle);

    // unproject relative angle from [0, 2pi) -> [0, bigTheta)
    double traceAngleRelativeTheta = angleBetween(traceAngle, base->angle) * (v_i->inVertex->bigTheta/(2*M_PI));
    assert(traceAngleRelativeTheta >= 0);
    auto [face, bary, dir] = traceVector<InFace>(base, Vector3d(1.0, 0.0, 0.0), distance, traceAngleRelativeTheta);
    return make_tuple(face, bary);
}

// angle should be the flat intrinsic angle (NOT phi) between base and the desired trace direction
template <typename T>
std::tuple<std::shared_ptr<T>, Eigen::Vector3d, Eigen::Vector2d> SPmesh::traceVector(auto base, Eigen::Vector3d baryCoords, double distance, double traceAngleRelativeTheta) {
    auto top = base->next->next;
    assert(base->edge->length > 0);
    Vector2d fi = Vector2d(0.0, 0.0);
    Vector2d fj = Vector2d(base->edge->length, 0.0);
    // compute interior flat/unprojected angle at vertex i (source of base): see fig 13 left of tutorial
    double theta_i = angleBetween(top->twin->angle, base->angle) * base->v->bigTheta/(2*M_PI);
    Vector2d fk = top->edge->length * Vector2d(cos(theta_i), sin(theta_i));
    Vector3d bary = baryCoords;

    if (traceAngleRelativeTheta > theta_i) {
        // this is only reached due to numerical imprecision: reflect the trace angle about base to guarantee that the trace path is in the triangle
        // or else you get negative bary coords
        traceAngleRelativeTheta -= abs(traceAngleRelativeTheta - theta_i) + 0.00001;
    }
    assert(traceAngleRelativeTheta <= theta_i);
    assert(theta_i >= 0 && theta_i < M_PI);
    // get the trace direction vector in 2D homogeneous local coords
    Vector3d dir;
    if (traceAngleRelativeTheta > M_PI_2 && traceAngleRelativeTheta < (3.0 * M_PI_2)) {
        // trace angle is in 2nd quadrant: reflect u about vertical axis of 2D local coordinate system so that trace angle lives in 1st quadrant for tan computation
        dir = Vector3d(-base->edge->length, tan(traceAngleRelativeTheta) * -base->edge->length, 0.0);
    } else {
        dir = Vector3d(base->edge->length, tan(traceAngleRelativeTheta) * base->edge->length, 0.0);
    }
    dir.normalize();

    // compute transformation from bary coords to local 2D (homogeneous) and its inverse
    Matrix3d A;
    A << fi(0), fj(0), fk(0),
         fi(1), fj(1), fk(1),
         1.0, 1.0, 1.0;
    bool invertible = false;
    Matrix3d inverse; // 2D hom to bary
    A.computeInverseWithCheck(inverse, invertible);
    assert(invertible);
    Vector3d baryDir = inverse * dir;

    // find ray-edge intersection with triangle edges within distance (if any)
    double t = distance;
    int minT = -1; // for triangle whose base is ij: minT=0 => closest/intersected edge is jk, minT=1 => ki closest, minT=2 => ij closest
    for (int i = 0; i < 3; i++) {
        if (baryDir(i) != 0.0) { // else direction is parallel to curr edge => no intersection
            double ti = -bary(i)/baryDir(i); // see pg 27 of tutorial
            if (ti > 0.0 && ti < t) {
                t = ti;
                minT = i;
            }
        }
    }

    Vector2d newDir = -Vector2d::Ones(); // for easy debugging
    Vector2d p = -Vector2d::Ones();
    // move to neighboring triangle and repeat trace until the target trace distance is reached
    int numEdgesIntersected =0;
    while (t < distance) {
        numEdgesIntersected++;
        Vector3d edgeIntersectBary = bary + t * baryDir;
        assert(isEqual(edgeIntersectBary[0], 0) || isEqual(edgeIntersectBary[1], 0) || isEqual(edgeIntersectBary[2], 0)); // one of the coords should be 0 since we're on an edge
        Vector2d intersectedEdge; // 2D local coords

        // parallel transport direction vector to next triangle: see pg 27 of tutorial
        if (minT == 0) {
            intersectedEdge = fk - fj; // 2D local coords of edge relative to current base (treat source of base as origin)
            base = base->next->twin; // base on next triangle (across the intersected edge)
        } else if (minT == 1) {
            intersectedEdge = fi - fk;
            base = base->next->next->twin;
        } else if (minT == 2) {
            intersectedEdge = fj - fi;
            base = base->twin;
        }
        intersectedEdge.normalize();
        // get new 2D local coords of the new triangle vertices relative to the new base
        top = base->next->next;
        fi = Vector2d(0.0, 0.0);
        fj = Vector2d(base->edge->length, 0.0);
        theta_i = angleBetween(top->twin->angle, base->angle) * base->v->bigTheta/(2*M_PI);
        fk = top->edge->length * Vector2d(cos(theta_i), sin(theta_i));

        Vector2d tijk = Vector2d(-intersectedEdge(1), intersectedEdge(0));
        Vector2d newEdge = fj - fi; // base edge (where intersection happened)
        newEdge.normalize();
        Vector2d tnew = Vector2d(-newEdge(1), newEdge(0));
        Vector2d dir2d = Vector2d(dir(0), dir(1));
        assert (dir(2) == 0.0);

        newDir = -((dir2d.dot(intersectedEdge) * newEdge) + (dir2d.dot(tijk) * tnew));
        newDir.normalize();
        if (isEqual(edgeIntersectBary(0), 0, 0.000001)) {
            p = edgeIntersectBary(1) * fj + edgeIntersectBary(2) * fi;
        } else if (isEqual(edgeIntersectBary(1), 0, 0.000001)) {
            p = edgeIntersectBary(0) * fi + edgeIntersectBary(2) * fj;
        } else if (isEqual(edgeIntersectBary(2), 0, 0.000001)) {
            p = edgeIntersectBary(0) * fj + edgeIntersectBary(1) * fi;
        }

        dir = Vector3d(newDir(0), newDir(1), 0.0);
        distance -= t;

        A << fi(0), fj(0), fk(0),
             fi(1), fj(1), fk(1),
             1.0, 1.0, 1.0;
        invertible = false;
        Matrix3d inverse;
        A.computeInverseWithCheck(inverse, invertible);
        assert(invertible);
        bary = inverse * Vector3d(p(0), p(1), 1.0);
        baryDir = inverse * dir;

        // intersect with triangle edges again
        t = distance;
        minT = -1;
        for (int i = 0; i < 3; i++) {
            if (baryDir(i) != 0.0) {
                double ti = -bary(i)/baryDir(i);
                if (ti > 0.0 && ti < t) {
                    t = ti;
                    minT = i;
                }
            }
        }
    }

    Vector3d newPointBary = bary + t * baryDir;
    Vector3d newDirBary = baryDir;
    // convert to barycentric coords wrt orientation of the face defined by its representative halfedge
    if (base != base->face->halfedge) {
        if (base->next == base->face->halfedge) {
            newPointBary = Vector3d(newPointBary(1), newPointBary(2), newPointBary(0));
            newDirBary = Vector3d(newDirBary(1), newDirBary(2), newDirBary(0));
        } else if (base->next->next == base->face->halfedge) {
            newPointBary = Vector3d(newPointBary(2), newPointBary(0), newPointBary(1));
            newDirBary = Vector3d(newDirBary(2), newDirBary(0), newDirBary(1));
        }
        base = base->face->halfedge;
        top = base->next->next;

        fi = Vector2d(0.0, 0.0);
        fj = Vector2d(base->edge->length, 0.0);
        // compute interior flat/unprojected angle at vertex i (source of base): see fig 13 left of tutorial
        theta_i = angleBetween(top->twin->angle, base->angle) * base->v->bigTheta/(2*M_PI);
        fk = top->edge->length * Vector2d(cos(theta_i), sin(theta_i));

        A << fi(0), fj(0), fk(0),
             fi(1), fj(1), fk(1),
             1.0, 1.0, 1.0;
    }

    Vector3d newDirLocal = A * newDirBary;

    // sanity check: barycentric coords must be positive and sum to 1
    assert(newPointBary[0] >= 0 && newPointBary[1] >= 0 && newPointBary[2] >= 0);
    double barySum = newPointBary[0] + newPointBary[1] + newPointBary[2];
    assert(isEqual(barySum, 1));
    double dirBarySum = newDirBary[0] + newDirBary[1] + newDirBary[2];
    assert(isEqual(dirBarySum, 0));

    return make_tuple(base->face, newPointBary, Vector2d(newDirLocal(0), newDirLocal(1)));
}

// algo 3: updates the signpost angles for every edge incident to the given vertex vi.
// specifically, updates both the angle phi_ij and phi_ji for every edge (i,j)
void SPmesh::updateVertex(shared_ptr<InVertex> v_i) {
//    cout << "updating vertex..." << endl;
    /// update incoming angles phi_ji (of halfedges that point to vi)
    shared_ptr<InHalfedge> currHalfEdge = v_i->halfedge;
    do {
        // update h_ji
        updateSignpost(currHalfEdge->next);
        currHalfEdge = currHalfEdge->twin->next;
    } while (currHalfEdge != v_i->halfedge);


    /// update outgoing angles phi_ij with respect to the fixed reference direction of the extrinsic triangle (set during initialization)
    // 1) update angle of canonical halfedge (h_ij0) by measuring the angle between the extrinsic face's reference direction and e_ij0
    shared_ptr<InHalfedge> h_ij0 = v_i->halfedge;
    shared_ptr<InHalfedge> h_j0i = h_ij0->twin;
    // get extrinsic face containing i by tracing from j0 to i along EXTRINSIC mesh
    auto [exTriangle, barycentricCoords, uTransformed] = traceFromIntrinsicVertex(
                h_j0i->v,
                h_j0i->edge->length,
                h_j0i->angle); // <- u vector in section 3.2.3 in the local polar coordinate system at v_j0

    // NOTE: reference direction e_xyz of extrinsic face is along edge xy where x=(1,0,0) and y=(0,1,0) in barycentric coords
    // In 2D local coords, it is the vector from (0,0) to (1,0)
    h_ij0->angle = argument(Vector2d(1,0), -uTransformed);
    // set barycentric coords of vi (located inside the extrinsic triangle)
    v_i->barycentricPos = barycentricCoords;
    v_i->exFace = exTriangle;

    // 2) update remaining outgoing angles phi_ij in ccw order
    currHalfEdge = v_i->halfedge;
    while (currHalfEdge->next->next->twin != v_i->halfedge) {
        updateSignpost(currHalfEdge);
        currHalfEdge = currHalfEdge->next->next->twin;
    }
}

// algo 5: takes an edge ij with opposite vertices k,l and flips it to be kl
// replaces triangles ijk, jil with klj,lki
// returns a nullptr if the edge is not flippable due to insufficient vertex degree, or if edge is non-delaunay
std::shared_ptr<InEdge> SPmesh::flipEdge(std::shared_ptr<InEdge> ij) {
    assert(ij);
    if (getDegree(ij->halfedge->v) < 3 || getDegree(ij->halfedge->twin->v) < 3) {
        cout << "cannot flip: endpoint has degree < 3" << endl;
        return nullptr;
    }

    std::shared_ptr<InHalfedge> lj = ij->halfedge->twin->next->next;
    std::shared_ptr<InHalfedge> ki = ij->halfedge->next->next;
    double l_ij = ij->length;
    double l_jk = ij->halfedge->next->edge->length;
    double l_ki = ki->edge->length;
    double l_il = ij->halfedge->twin->next->edge->length;
    double l_lj = lj->edge->length;

    // theta_i = total internal angle of vertex i in the diamond (i.e. angle kil)
    double theta = getAngleFromEdgeLengths(l_ij, l_jk, l_ki) + getAngleFromEdgeLengths(l_il, l_lj, l_ij);

    std::shared_ptr<InVertex> vi = ij->halfedge->v;
    std::shared_ptr<InVertex> vj = ij->halfedge->next->v;
    std::shared_ptr<InVertex> vk = ki->v;
    std::shared_ptr<InVertex> vl = lj->v;

    // update signpost mesh connectivity
    /// 1) remove the existing 2 intrinsic faces adjacent to edge
    _vertPairToEdge.clear();
    eraseTriangle(ij->halfedge->face);
    eraseTriangle(ij->halfedge->twin->face);
    /// 2) insert 2 new intrinsic faces adjacent to flipped edge
    shared_ptr<InFace> tri_ilk = insertTriangle(vi, vl, vk);
    shared_ptr<InFace> tri_jkl = insertTriangle(vj, vk, vl);

    // get new length of flipped edge
    double l_kl = baseLength(l_ki, l_il, theta);
    shared_ptr<InHalfedge> kl = tri_jkl->halfedge->next;
    kl->edge->length = l_kl;
    assert(kl->twin->edge->length == l_kl);

    // update signposts
    /// update angle of HE lk using lj
    assert(lj->next->next->twin == tri_ilk->halfedge->next);
    updateSignpost(lj);
    /// update angle of HE kl using ki
    assert(ki->next->next->twin == kl);
    updateSignpost(ki);

    validate();

    return kl->edge;
}

// algo 6: returns the distance between points p and q inside a triangle with the three given edge lengths.
// p and q are barycentric coords.
double SPmesh::distance(double l_12, double l_23, double l_31, const Vector3d p, const Vector3d q) {
    Vector3d u = q - p;
    double d = -(
        pow(l_12, 2)*u[0]*u[1] +
        pow(l_23, 2)*u[1]*u[2] +
        pow(l_31, 2)*u[2]*u[0]
    );
    return sqrt(d);
}

// algo 7: inserts a new intrinsic vertex in the given intrinsic face ijk at the position specified using barycentric coords
// barycentric coords must be positive and sum to 1
shared_ptr<InVertex> SPmesh::insertVertex(std::shared_ptr<InFace> face, Vector3d barycentricCoords) {
//    cout << "inserting vertex..." << endl;
    // get face's verts (ccw orientation)
    std::shared_ptr<InVertex> v_i = face->halfedge->v;
    std::shared_ptr<InVertex> v_j = face->halfedge->next->v;
    std::shared_ptr<InVertex> v_k = face->halfedge->next->next->v;
    // get face's edge lengths
    double l_ij = face->halfedge->edge->length;
    double l_jk = face->halfedge->next->edge->length;
    double l_ki = face->halfedge->next->next->edge->length;

    // make new vertex
    std::shared_ptr<InVertex> p = make_shared<InVertex>();
    p->exVertex = nullptr; // new intrinsic vertices do not correspond to any extrinsic vertex
    p->bigTheta = (2*M_PI); // p is inside an intrinsic triangle => no curvature at p => no need to project; same angle sum as if p were in the plane (2pi)
    _verts.insert(p);

    // update signpost mesh connectivity
    /// 1) remove the existing intrinsic face
    _vertPairToEdge.clear(); // clear 1-ring neighborhood cache
    eraseTriangle(face);
    /// 2) insert 3 new intrinsic faces around p and update half edge connectivity (will update signposts afterwards)
    shared_ptr<InFace> ijp = insertTriangle(v_i, v_j, p);
    p->halfedge = ijp->halfedge->next->next; // assign halfedge pi to p
    shared_ptr<InFace> jkp = insertTriangle(v_j, v_k, p);
    shared_ptr<InFace> kip = insertTriangle(v_k, v_i, p);

    // compute and set edge lengths of all new edges incident to new vertex p
    // barycentric coordinates of vertices are: vi = (1,0,0),   vj = (0,1,0),   vk = (0,0,1).
    std::shared_ptr<InEdge> e_ip = p->halfedge->edge;
    e_ip->length = distance(l_ij, l_jk, l_ki, barycentricCoords, Vector3d(1,0,0)); // distance along intrinsic (flattened) triangle from vi to p
    std::shared_ptr<InEdge> e_jp = p->halfedge->next->next->twin->edge;
    e_jp->length = distance(l_ij, l_jk, l_ki, barycentricCoords, Vector3d(0,1,0)); // vj to p
    std::shared_ptr<InEdge> e_kp = p->halfedge->next->next->twin->next->next->twin->edge;
    e_kp->length = distance(l_ij, l_jk, l_ki, barycentricCoords, Vector3d(0,0,1)); // vk to p

    checkTriangleInequality(ijp);
    checkTriangleInequality(jkp);
    checkTriangleInequality(kip);

    // update signpost angles around p and affected neighbors
    updateVertex(p);

    validate();

    return p;
}

// algo 9: given barycentric coords p for a point in triangle ijk returns polar coords r,phi of the vector from i to p
// phi is relative to edge ij
std::pair<double, double> SPmesh::vectorToPoint(double l_ij, double l_jk, double l_ki, const Eigen::Vector3d &i, const Eigen::Vector3d &j, const Eigen::Vector3d &p) {
    double r_pi = distance(l_ij, l_jk, l_ki, i, p); // distance from p to i
    double r_jp = distance(l_ij, l_jk, l_ki, j, p); // distance from j to p
    double phi_ip = getAngleFromEdgeLengths(l_ij, r_jp, r_pi); // angle from ij to ip
    assert(!isnan(phi_ip));
    return std::make_pair(r_pi, phi_ip);
}

// algo 11: takes the barycentric coords of point on extrinsic face and returns the barycentric coordinates of that position on its intrinsic face
std::tuple<std::shared_ptr<InFace>, Eigen::Vector3d> SPmesh::pointQuery(std::shared_ptr<ExFace> xyz, Eigen::Vector3d& p) {
    std::shared_ptr<ExHalfedge> xy = xyz->halfedge;

    std::shared_ptr<ExVertex> x = xy->v;
    std::shared_ptr<ExVertex> y = xy->next->v;
    std::shared_ptr<ExVertex> z = xy->next->next->v;

    double l_xy = (y->pos - x->pos).norm();
    double l_yz = (z->pos - y->pos).norm();
    double l_zx = (x->pos - z->pos).norm();

    // vector from x to p relative to xy
    std::pair<double, double> rPhi = vectorToPoint(l_xy, l_yz, l_zx, Vector3d(1, 0, 0), Vector3d(0, 1, 0), p);

    double angle = (rPhi.second / x->bigTheta) * (2*M_PI) + xy->angle;

    return traceFromExtrinsicVertex(x, rPhi.first, angle);
}


//////////////////////////////////////////////
//////////////// HELPERS /////////////////////
//////////////////////////////////////////////

// HELPER :get the extrinsic position of an ORIGINAL intrinsic vertex
Vector3d SPmesh::getVPos(shared_ptr<InVertex> v) {
    return v->exVertex->pos;
}

// HELPER: for a triangle with edge lengths l_ij, l_jk, l_ki, returns the interior angle at vertex i
double SPmesh::getAngleFromEdgeLengths(double l_ij, double l_jk, double l_ki) {
    // see fig. 9 of paper
    double cosTheta = (pow(l_ij, 2) + pow(l_ki, 2) - pow(l_jk, 2)) / (2*l_ij*l_ki);
    cosTheta = std::clamp(cosTheta, -1., 1.);
    assert(-1 <= cosTheta && cosTheta <=1);
    return acos( cosTheta );
}

// HELPER: returns the angle between two 3-vectors in [0, pi]
double SPmesh::getAngle(Vector3d u, Vector3d v) {
    return acos(u.normalized().dot(v.normalized()));
}


// HELPER: returns the third side length of a triangle with side lengths a,b meeting at angle theta (law of cos)
double SPmesh::baseLength(double a, double b, double theta) {
    assert(theta < M_PI);
    assert(a>0 && b>0);
    return sqrt(a*a + b*b - 2*a*b*cos(theta));
}

// HELPER: returns smallest unsigned angle between the points on the unit circle given by angles a,b
double SPmesh::angleBetween(double a, double b) {
    double diff = abs(a - b); // get positive diff
    if (diff > M_PI) diff = (2*M_PI) - diff; // smaller angle should be in [0, pi]
    return diff;
}


// HELPER: returns the ccw angle FROM vector u TO vector v in the range [0, 2*pi).
//Imagine fixing u as the x-axis in R^2 and going ccw to find v.
double SPmesh::argument(Vector2d u, Vector2d v) {
    Vector2d a = u.normalized();
    Vector2d b = v.normalized();
    // adapted from  https://stackoverflow.com/questions/40286650/how-to-get-the-anti-clockwise-angle-between-two-2d-vectors
    double angle = atan2(a[0]*b[1] - a[1]*b[0], a[0]*b[0] + a[1]*b[1]);
    if (angle < 0) {
        angle += (2*M_PI);
    }
    return angle;
}

// HELPER: calculates the area of a face
double SPmesh::getArea(std::shared_ptr<InFace> face) {
    double a = face->halfedge->edge->length;
    double b = face->halfedge->next->edge->length;
    double c = face->halfedge->next->next->edge->length;
    double s = (a + b + c) / 2.0;

    return sqrt(s * (s - a) * (s - b) * (s - c));
}

pair<shared_ptr<InVertex>, shared_ptr<InVertex>> getOrderedVertexPair(shared_ptr<InVertex> v0, shared_ptr<InVertex> v1) {
    return minmax(v0, v1);
}

// HELPER: removes the given intrinsic triangle.
// Also removes any incident InEdges that would no longer border two faces after the current removal
// All of the face's edges and their halfedges are preserved
void SPmesh::eraseTriangle(shared_ptr<InFace> tri) {
    // store adjacencies of vertices in the triangle
    shared_ptr<InHalfedge> startHalfEdge = tri->halfedge;
    shared_ptr<InHalfedge> currHalfEdge = startHalfEdge;
    do {
        // mark v_curr and v_twin as neighbors of each other (may be negated if edge between them gets deleted)
        // i.e. store shared edge incident to v_curr and v_twin
        shared_ptr<InEdge> e_curr = currHalfEdge->edge;
        shared_ptr<InVertex> v_curr = currHalfEdge->v;
        shared_ptr<InVertex> v_twin = currHalfEdge->twin->v;
        assert(e_curr && _edges.contains(e_curr));
        assert(v_curr && _verts.contains(v_curr));
        assert(v_twin && _verts.contains(v_twin));
        auto [v0, v1] = getOrderedVertexPair(v_curr, v_twin);
        _vertPairToEdge[v0][v1] = e_curr;
        currHalfEdge = currHalfEdge->next;
    } while(currHalfEdge != startHalfEdge);


    // removal: set face pointer of each halfedge to null
    startHalfEdge = tri->halfedge;
    currHalfEdge = startHalfEdge;
    do {
        currHalfEdge->face = nullptr;
        // if edge is no longer adjacent to any faces
        if (!currHalfEdge->twin->face) {
            // 1) re-assign representative halfedge of the endpoint vertices v0 and v1, if necessary
            // the reference direction does not change; all other halfedges' phi angles are unchanged.
            shared_ptr<InHalfedge> currTwin = currHalfEdge->twin;
            shared_ptr<InVertex> v0 = currHalfEdge->v;
            shared_ptr<InVertex> v1 = currTwin->v;

            if (v0->halfedge == currHalfEdge) {
                // pick next halfedge clockwise around the curr vertex
                v0->halfedge = currHalfEdge->twin->next;
            }
            if (v1->halfedge == currTwin) {
                // pick next halfedge clockwise around the curr vertex
                v1->halfedge = currTwin->twin->next;
            }

            // 2) the endpoints of the edge to remove will no longer be neighbors
            shared_ptr<InEdge> currEdge = currHalfEdge->edge;
            auto [v_a, v_b] = getOrderedVertexPair(v0, v1);
            _vertPairToEdge[v_a][v_b] = nullptr;

            // 3) remove curr InEdge entirely since it is no longer adjacent to any faces (e.g. if two adjacent faces are erased, shared edge is also erased on the 2nd call)
            _edges.erase(currHalfEdge->edge);
            _halfedges.erase(currHalfEdge);
            _halfedges.erase(currHalfEdge->twin);
        }

        currHalfEdge = currHalfEdge->next;
    } while(currHalfEdge != startHalfEdge);

    _faces.erase(tri);
}

// HELPER: retrieves the halfedge belonging to the given edge with the matching source vertex, else nullptr
shared_ptr<InHalfedge> SPmesh::getHalfEdgeWithSource(shared_ptr<InEdge> edge, shared_ptr<InVertex> sourceVertex) const {
    assert(sourceVertex);
    assert(_edges.contains(edge));
    shared_ptr<InHalfedge> he_0 = edge->halfedge;
    assert(he_0);
    shared_ptr<InHalfedge> he_1 = he_0->twin;
    return (he_0->v == sourceVertex) ? he_0 : he_1;
}

// HELPER: inserts an intrinsic triangle whose vertices are v0, v1, v2. Assumes v0,v1,v2 are provided in ccw order. halfedge (v0, v1) will be assigned as the representative of the face
// does not assume that there exist faces adjacent to the triangle (v0,v1,v2). Can be called in any order to insert multiple adjacent faces.
shared_ptr<InFace> SPmesh::insertTriangle(shared_ptr<InVertex> v0, shared_ptr<InVertex> v1, shared_ptr<InVertex> v2) {
    // make empty face
    shared_ptr<InFace> newFace = make_shared<InFace>(InFace{}); // pick arbitrary half edge representative
    _faces.insert(newFace);
    // construct new edges and their corresponding half-edges on this face if necessary.
    shared_ptr<InHalfedge> newHalfEdge;
    shared_ptr<InHalfedge> prevHalfEdge;
    auto vertices = vector{v0, v1, v2};
    shared_ptr<InHalfedge> he_0;
    shared_ptr<InHalfedge> he_1;
    shared_ptr<InHalfedge> he_2;
    for (int i = 0; i < 3; i++) {
        shared_ptr<InVertex> currVertex = vertices[i];
        shared_ptr<InVertex> nextVertex = vertices[(i+1) % 3];
        shared_ptr<InEdge> edge = getEdge(currVertex, nextVertex); // edge (i, i->next) containing halfedge to be constructed (can be null)
        if (edge) {
            // grab the halfedge contained in this edge whose source is at currVertex (if it exists)
            shared_ptr<InHalfedge> existingHalfEdge = getHalfEdgeWithSource(edge, currVertex);
            if (existingHalfEdge) {
                existingHalfEdge->face = newFace;
                newHalfEdge = existingHalfEdge; // use the existing halfedge as the 'new halfedge'
            } else {
                // edge and twin of the desired halfedge exist (but not the halfedge itself): make a new halfedge and glue to the existing twin (the single representative)
                // initialize angle to -1 to indicate invalid
                newHalfEdge = make_shared<InHalfedge>(InHalfedge{currVertex, nullptr, edge->halfedge, edge, nullptr, -1});
                _halfedges.insert(newHalfEdge);
                newHalfEdge->twin->twin = newHalfEdge;
                newHalfEdge->edge = edge;
                newHalfEdge->face = newFace;
            }
        } else {
            newHalfEdge = make_shared<InHalfedge>(InHalfedge{currVertex, nullptr, nullptr, nullptr, nullptr, -1});
            _halfedges.insert(newHalfEdge);
            // edge does not exist: make one and connect it to the new halfEdge
            edge = make_shared<InEdge>(InEdge{newHalfEdge, -1}); // initialize length to -1 to indicate invalid
            _edges.insert(edge);
            newHalfEdge->edge = edge;
            newHalfEdge->face = newFace;
            // mark the endpoints of the new edge as neighbors
            auto [va, vb] = getOrderedVertexPair(currVertex, nextVertex);
            _vertPairToEdge[va][vb] = edge;
        }

        // record current half edge for setting next pointers later. Set he_i = current halfedge (whose base is at v_i)
        switch(i) {
            case 0: he_0 = newHalfEdge; break;
            case 1: he_1 = newHalfEdge; break;
            case 2: he_2 = newHalfEdge; break;
        }
    }

    // set next pointers
    he_0->next = he_1;
    he_1->next = he_2;
    he_2->next = he_0;

    // set face representative as halfedge whose source is v0 (bary coords = (1,0,0))
    newFace->halfedge = he_0;
    return newFace;
}

// HELPER: returns the InEdge containing v0 and v1 as endpoints or nullptr if it doesn't exist
shared_ptr<InEdge> SPmesh::getEdge(shared_ptr<InVertex> v0, shared_ptr<InVertex> v1) const {
    auto [va, vb] = getOrderedVertexPair(v0, v1);
    // check that the ordered pair (va, vb) exists in the map
    if (!_vertPairToEdge.contains(va) || !_vertPairToEdge.at(va).contains(vb))
        return nullptr;

    return _vertPairToEdge.at(va).at(vb);
}

// HELPER: returns vertex degree
int SPmesh::getDegree(const shared_ptr<InVertex> &v) {
    int i = 0;
    shared_ptr<InHalfedge> curr = v->halfedge;
    do {
        curr = curr->twin->next;
        i += 1;
    } while (curr != v->halfedge);
    return i;
}


//////////////////////////////////////////////
////////////// TRIANGULATION /////////////////
//////////////////////////////////////////////

// HELPER: returns true if the opposite angles on the adjacent faces of the given edge have sum <= pi. See diagram of equation 4.
bool SPmesh::edgeIsDelaunay(shared_ptr<InEdge> edge) {
    shared_ptr<InHalfedge> halfedge = edge->halfedge;
    shared_ptr<InHalfedge> twin = halfedge->twin;
    // get side lengths of the face lji incident to the representative halfedge of the given edge
    double l_ji = edge->length;
    double l_il = halfedge->next->edge->length;
    double l_lj = halfedge->next->next->edge->length;
    // get side lengths of the other face kij
    double l_ij = l_ji;
    double l_jk = twin->next->edge->length;
    double l_ki = twin->next->next->edge->length;

    double theta_l = getAngleFromEdgeLengths(l_lj, l_ji, l_il);
    double theta_k = getAngleFromEdgeLengths(l_ki, l_ij, l_jk);

    return theta_l + theta_k <= M_PI;
}

// HELPER: returns true if the triangle does not satisfy the minimum angle bound
// course algorithm 4
bool SPmesh::shouldRefine(std::shared_ptr<InFace> tri, double minAngle) {
    double lengths[3] = {
        tri->halfedge->edge->length,
        tri->halfedge->next->edge->length,
        tri->halfedge->next->next->edge->length
    };

    shared_ptr<InVertex> vertices[3] = {
        tri->halfedge->v,
        tri->halfedge->next->v,
        tri->halfedge->next->next->v
    };

    for (int i = 0; i < 3; i++) {
        if (getAngleFromEdgeLengths(lengths[i], lengths[(i + 1) % 3], lengths[(i + 2) % 3]) < minAngle && getDegree(vertices[i]) > 1) {
            return true;
        }
    }

    return false;
}

shared_ptr<InVertex> SPmesh::insertCircumcenter(shared_ptr<InFace> face) {
    // calculate barycentric circumcenter using edge lengths
    double l_ij = face->halfedge->edge->length;
    double l_jk = face->halfedge->next->edge->length;
    double l_ki = face->halfedge->next->next->edge->length;

    double v_i = l_jk*l_jk * (l_ij*l_ij + l_ki*l_ki - l_jk*l_jk);
    double v_j = l_ki*l_ki * (l_jk*l_jk + l_ij*l_ij - l_ki*l_ki);
    double v_k = l_ij*l_ij * (l_ki*l_ki + l_jk*l_jk - l_ij*l_ij);

    Vector3d circumcenter = Vector3d(v_i, v_j, v_k) / (v_i + v_j + v_k);

    assert(isEqual(circumcenter[0] + circumcenter[1] + circumcenter[2], 1));

    // trace to find the face and coordinates where the circumcenter is
    shared_ptr<InHalfedge> base = face->halfedge;
    Vector3d barycenter = Vector3d(1.0/3, 1.0/3, 1.0/3);

    // convert to 2d local to get angle
    Vector2d fi = Vector2d(0.0, 0.0);
    Vector2d fj = Vector2d(l_ij, 0.0);
    double theta_i = angleBetween(base->next->next->twin->angle, base->angle) * base->v->bigTheta/(2*M_PI);
    Vector2d fk = l_ki * Vector2d(cos(theta_i), sin(theta_i));

    Matrix3d A;
    A << fi(0), fj(0), fk(0),
         fi(1), fj(1), fk(1),
         1.0, 1.0, 1.0;

    Vector3d barycenterLocal = A * barycenter;
    Vector3d circumcenterLocal = A * circumcenter;

    Vector2d traceDir = Vector2d(circumcenterLocal[0] - barycenterLocal[0], circumcenterLocal[1] - barycenterLocal[1]);

    double dist = traceDir.norm();
    double angle = argument(fj, traceDir);

    auto [endFace, circumBary, dir] = traceVector<InFace>(base, barycenter, dist, angle);

    // don't insert if we're on an edge
    if (isEqual(circumBary[0], 0, 0.0005) || isEqual(circumBary[1], 0, 0.0005) || isEqual(circumBary[2], 0, 0.0005)) return nullptr;

    return insertVertex(endFace, circumBary);
}

// course algorithm 1
void SPmesh::flipToDelaunay(unordered_set<shared_ptr<InEdge>>& edgesToCheck, double minAngle, int maxFlips) {
    // enqueue all edges
    queue<shared_ptr<InEdge>> edgeQ;
    for (shared_ptr<InEdge> edge: edgesToCheck) {
        edgeQ.push(edge);
    }

    // keep track of faces that need refinement
//    unordered_set<shared_ptr<InFace>> facesToCheck;

    // iterate until queue is empty or hit cap
    int i = 0;
    while (!edgeQ.empty() && i < maxFlips) {
        shared_ptr<InEdge> nextEdge = edgeQ.front();
        edgeQ.pop();
        edgesToCheck.erase(edgesToCheck.find(nextEdge));
        // confirm that edge still exists and is not delaunay
        if (_edges.contains(nextEdge) && !edgeIsDelaunay(nextEdge)) {
            // if not delaunay then flip
            shared_ptr<InEdge> flipped = flipEdge(nextEdge);
            i++;

            // enqueue adjacent edges that aren't already in the queue
            shared_ptr<InEdge> adjacentEdges[4] = {
                flipped->halfedge->next->edge,
                flipped->halfedge->next->next->edge,
                flipped->halfedge->twin->next->edge,
                flipped->halfedge->twin->next->next->edge
            };
            for (shared_ptr<InEdge> adjacent: adjacentEdges) {
                if (!edgesToCheck.contains(adjacent)) {
                    edgeQ.push(adjacent);
                    edgesToCheck.insert(adjacent);
                }
            }

//            // track adjacent faces that need refinement
//            shared_ptr<InFace> adjacentFace = flipped->halfedge->face;
//            if (!facesToCheck.contains(adjacentFace) && shouldRefine(adjacentFace, minAngle)) {
//                facesToCheck.insert(adjacentFace);
//            }
//            adjacentFace = flipped->halfedge->twin->face;
//            if (!facesToCheck.contains(adjacentFace) && shouldRefine(adjacentFace, minAngle)) {
//                facesToCheck.insert(adjacentFace);
//            }
        }
    }

    cout << "mesh is delaunay after " << i << "/" << maxFlips << " flips" << endl;

//    return facesToCheck;
}

// course algorithm 5
void SPmesh::delaunayRefine(double minAngle, int maxInsertions) {
    int flips = 0;
    int insertions = 0;

    // enqueue all edges
    deque<shared_ptr<InEdge>> edgeQ;
    unordered_set<shared_ptr<InEdge>> inEdgeQ = _edges;
    for (shared_ptr<InEdge> edge: _edges) {
        edgeQ.push_back(edge);
    }

    // enqueue all faces that need refinement
    /// since faces are only enqueued adjacent to an insert or flip the pointers should always be new so there's
    /// no point in keeping a set to check which faces are already in the queue
    priority_queue<pair<double, shared_ptr<InFace>>> faceQ;
    for (shared_ptr<InFace> face: _faces) {
        if (shouldRefine(face, minAngle)) {
            faceQ.push(make_pair(getArea(face), face));
        }
    }

    cout << "beginning delaunay refinement with " << faceQ.size() << "/" << _faces.size() << "faces needing refinement" << endl;

    // main loop that keeps inserting and flipping
    do  {

        // inner loop for flipping to delaunay
        /// this is redundant given the above flipToDelaunay but the geometry-central codebase also includes this redundancy
        /// it's also easier to just add directly to the face queue here rather than returning
        /// the faces from flipToDelaunay
        while (!edgeQ.empty()) {
            shared_ptr<InEdge> nextEdge = edgeQ.front();
            edgeQ.pop_front();
            inEdgeQ.erase(nextEdge);
            // if edge no longer exists or is delaunay don't need to flip
            if (!_edges.contains(nextEdge) || edgeIsDelaunay(nextEdge)) continue;

            // flip edge
            shared_ptr<InEdge> flipped = flipEdge(nextEdge);
            flips++;

            // enqueue adjacent faces that need refinement
            shared_ptr<InFace> adjacentFace = flipped->halfedge->face;
            if (shouldRefine(adjacentFace, minAngle)) {
                faceQ.push(make_pair(getArea(adjacentFace), adjacentFace));
            }
            adjacentFace = flipped->halfedge->twin->face;
            if (shouldRefine(adjacentFace, minAngle)) {
                faceQ.push(make_pair(getArea(adjacentFace), adjacentFace));
            }

            // enqueue adjacent edges that aren't already in the queue
            shared_ptr<InEdge> adjacentEdges[4] = {
                flipped->halfedge->next->edge,
                flipped->halfedge->next->next->edge,
                flipped->halfedge->twin->next->edge,
                flipped->halfedge->twin->next->next->edge
            };
            for (shared_ptr<InEdge> adjacent: adjacentEdges) {
                if (!inEdgeQ.contains(adjacent)) {
                    edgeQ.push_back(adjacent);
                    inEdgeQ.insert(adjacent);
                }
            }
        }

        // break if max insertions has been reached
        if (insertions >= maxInsertions) break;

        // insert a single vertex
        if (!faceQ.empty()) {
            // get biggest face
            shared_ptr<InFace> nextFace = faceQ.top().second;
            double area = faceQ.top().first;
            faceQ.pop();

            // check if face still exists, needs refinement, and hasn't changed its area
            /// I don't know if there should ever be a case where area changes but I'm checking it for now to be safe
            if (!_faces.contains(nextFace) || area != getArea(nextFace) || !shouldRefine(nextFace, minAngle)) continue;

            // insert circumcenter
            shared_ptr<InVertex> inserted = insertCircumcenter(nextFace);

            // make sure something was actually inserted (does nothing if trace lands on an edge)
            if (inserted) {
                insertions++;

                // check new faces next to inserted vertex
                shared_ptr<InHalfedge> curr = inserted->halfedge;
                do {
                    shared_ptr<InFace> adjFace = curr->face;

                    // enqeue if necessary
                    if (shouldRefine(adjFace, minAngle)) {
                        faceQ.push(make_pair(getArea(adjFace), adjFace));
                    }

                    // check adjacent face edges and enqueue if necessary
                    shared_ptr<InHalfedge> faceHE = curr;
                    do {
                        if (!inEdgeQ.contains(faceHE->edge)) {
                            edgeQ.push_back(faceHE->edge);
                            inEdgeQ.insert(faceHE->edge);
                        }

                        faceHE = faceHE->next;
                    } while (faceHE != curr);

                    curr = curr->next->next->twin;
                } while (curr != inserted->halfedge);
            }
        }

        if (insertions % 100 == 0) {
            cout << insertions << "/" << maxInsertions << " insertions" << endl;
        }

    } while (!edgeQ.empty() || !faceQ.empty());

    cout << "finished refining after " << flips << " flips and " << insertions << " insertions" << endl;

    int badFaces = 0;
    for (shared_ptr<InFace> face: _faces) {
        if (shouldRefine(face, minAngle)) badFaces++;
    }
    cout << badFaces << "/" << _faces.size() << " faces still violate min angle bound" << endl;
}


//////////////////////////////////////////////
/////////////// VALIDATION ///////////////////
//////////////////////////////////////////////

void SPmesh::checkCircular(const shared_ptr<InHalfedge> &halfedge) {
    assert(halfedge == halfedge->next->next->next);
    assert(halfedge->face == halfedge->next->face);
    assert(halfedge->face == halfedge->next->next->face);
}

void SPmesh::checkTwin(const shared_ptr<InHalfedge> &halfedge) {
    assert(halfedge == halfedge->twin->twin);
    assert(halfedge->v == halfedge->twin->next->v);
    assert(halfedge->edge == halfedge->twin->edge);
}

void SPmesh::checkFaces() {
    for (const shared_ptr<InFace> &face: _faces) {
        assert(_halfedges.contains(face->halfedge));
        assert(_halfedges.contains(face->halfedge->next));
        assert(_halfedges.contains(face->halfedge->next->next));
        assert(face->halfedge->face == face);
        assert(face->halfedge->next->face == face);
        assert(face->halfedge->next->next->face == face);
    }
    for (const shared_ptr<InEdge> &edge: _edges) {
        assert(edge->halfedge->edge == edge);
        assert(edge->halfedge->twin->edge == edge);
    }
    assert(_halfedges.size() == _faces.size() * 3);
    assert(_edges.size() * 2 == _faces.size() * 3);
}

void SPmesh::checkVertices() {
    for (const shared_ptr<InVertex> &v: _verts) {
        assert(getDegree(v) >= 2);
        shared_ptr<InHalfedge> curr = v->halfedge->twin->next;
        while (curr != v->halfedge) {
            assert(curr->v == v);
            curr = curr->twin->next;
        }
        assert(curr == v->halfedge);
        assert(v->bigTheta > 0);
    }
}

bool SPmesh::isEqual(double a, double b, double epsilon) {
    return abs(a-b) < epsilon;
}

void SPmesh::checkTriangleInequality(const shared_ptr<InFace> face) {
    double l_ij = face->halfedge->edge->length;
    double l_jk = face->halfedge->next->edge->length;
    double l_ki = face->halfedge->next->next->edge->length;
    assert(l_ij + l_jk >= l_ki);
    assert(l_jk + l_ki >= l_ij);
    assert(l_ki + l_ij >= l_jk);
}

void SPmesh::validateSignpost() {
    for (const shared_ptr<InHalfedge> &halfedge: _halfedges) {
        assert(halfedge->angle >= 0);
        assert(halfedge->angle < 2*M_PI);
        // strictly positive edge lengths
        assert(halfedge->edge->length > 0);
        // halfedge and its immediate radial neighbors should not have the same phi angle (or else there's a degenerate triangle).
            // can be violated by bad edge flips
        assert(!isEqual(halfedge->angle, halfedge->twin->next->angle));
        assert(!isEqual(halfedge->angle, halfedge->next->next->twin->angle));
        // the flat angle between a halfedge and its next radial neighbor should always be <= pi (or else there exists a triangle with interior angle sum > pi)
        double rightAngle = halfedge->angle;
        double leftAngle = halfedge->next->next->twin->angle;
        if (leftAngle < rightAngle)
            leftAngle += 2*M_PI;
        assert( angleBetween(leftAngle, rightAngle) * halfedge->v->bigTheta/(2*M_PI)  <= M_PI);
    }
    for (const shared_ptr<InFace> &face: _faces) {
        // edges satisfy triangle inequalities
        checkTriangleInequality(face);
    }
}

void SPmesh::validate() {
    for (const shared_ptr<InHalfedge> &halfedge: _halfedges) {
        checkCircular(halfedge);
        checkTwin(halfedge);
        assert(halfedge->face);
        assert(_edges.contains(halfedge->edge));
        assert(_verts.contains(halfedge->v));
        assert(_faces.contains(halfedge->face));
    }
    for (const shared_ptr<InEdge> &edge: _edges) {
        assert(edge->halfedge); // non-null representative halfedge
        assert(_halfedges.contains(edge->halfedge));
    }
    checkFaces();
    checkVertices();
    validateSignpost(); // signpost-specific assertions
}


//////////////////////////////////////////////
///////////// VISUALIZATION //////////////////
//////////////////////////////////////////////

void SPmesh::assignColors() {
    // adapted from https://www.geeksforgeeks.org/graph-coloring-set-2-greedy-algorithm/#
    _faceColors.reserve(_faces.size());
    int numColors = colors.size();
    int colorsUsed[numColors];
    for (int i = 0; i < numColors; i++) colorsUsed[i] = false;

    auto face = _faces.begin();
    _faceColors[*face] = rand() % numColors; // assign first color randomly

    for (face = _faces.begin()++; face != _faces.end(); face++) {
        // flag colors already used by adjacent faces
        std::shared_ptr<InFace> adj1 = (*face)->halfedge->twin->face;
        std::shared_ptr<InFace> adj2 = (*face)->halfedge->next->twin->face;
        std::shared_ptr<InFace> adj3 = (*face)->halfedge->next->next->twin->face;
        if (_faceColors.contains(adj1)) colorsUsed[_faceColors[adj1]] = true;
        if (_faceColors.contains(adj2)) colorsUsed[_faceColors[adj2]] = true;
        if (_faceColors.contains(adj3)) colorsUsed[_faceColors[adj3]] = true;

        // select random color, if its used, loop over array until an unused one is found
        // mod to wrap around to beginning
        int firstAvailableColor;
        for (firstAvailableColor = rand() % numColors; firstAvailableColor < numColors; firstAvailableColor = (firstAvailableColor + 1) % numColors) {
            if (!colorsUsed[firstAvailableColor]) break;
        }

        // assign color to face
        _faceColors[*face] = firstAvailableColor;

        // reset color flags
        for (int i = 0; i < numColors; i++) colorsUsed[i] = false;
    }
}

Eigen::Vector3d SPmesh::getBaryCoords(Eigen::Vector3d &p, Eigen::Vector3d &v1, Eigen::Vector3d &v2, Eigen::Vector3d &v3) {
    // courtesy of https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    Vector3d a = v2 - v1;
    Vector3d b = v3 - v1;
    Vector3d c = p - v1;
    double d00 = a.dot(a);
    double d01 = a.dot(b);
    double d11 = b.dot(b);
    double d20 = c.dot(a);
    double d21 = c.dot(b);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1. - v - w;
    return Vector3d(u, v, w);
}

// gets distance from p to the edge v1-v2 of triangle with side lengths l_ij, l_jk, l_ki
// p, v1, v2 are all barycentric coords
double SPmesh::distanceToEdge(Eigen::Vector3d &p, Eigen::Vector3d &v1, Eigen::Vector3d &v2, double l_ij, double l_jk, double l_ki) {
    Vector3d u = p - v1;
    Vector3d v = v2 - v1;
    Vector3d projected = v1 + v * (u.dot(v) / v.dot(v));
    assert(isEqual(p[0] + p[1] + p[2], 1));
    assert(isEqual(projected[0] + projected[1] + projected[2], 1));
    assert(projected[0] >= 0 && projected[1] >= 0 && projected[2] >= 0);
    assert(isEqual(projected[0], 0) || isEqual(projected[1], 0) || isEqual(projected[2], 0)); // projection is on an edge: one bary coord == 0
    return distance(l_ij, l_jk, l_ki, p.cwiseMax(0), projected);

//    Vector2d fi = Vector2d(0.0, 0.0);
//    Vector2d fj = Vector2d(l_ij, 0.0);
//    float theta_i = getAngleFromEdgeLengths(l_ij, l_jk, l_ki);
//    Vector2d fk = l_ki * Vector2d(cos(theta_i), sin(theta_i));
//    Matrix3d A;
//    A << fi(0), fj(0), fk(0),
//         fi(1), fj(1), fk(1),
//         1.0, 1.0, 1.0;

//    Vector3d p_loc = A * p;
//    Vector3d v1_loc = A * v1;
//    Vector3d v2_loc = A * v2;


//    Vector3d b_a = v2_loc - v1_loc;
//    Vector3d c_a = p_loc - v1_loc;

//    return b_a.cross(c_a).norm() / sqrt(b_a.head<2>().norm());
}

void SPmesh::computeMeanIntrinsicEdgeLength() {
    double totalLen = 0;
    for (shared_ptr<InEdge> e : _edges) {
        totalLen += e->length;
    }
    _meanIntrinsicEdgeLength = totalLen/_edges.size();
}

Vector3d SPmesh::getColor(const Triangle* tri, Eigen::Vector3d point, const Eigen::Vector3d &camPos) {
    shared_ptr<ExFace> exFace = _exMesh.getExTriangle(tri->getIndex());

    // calculate barycentric coords on exFace
    Vector3d v_i = exFace->halfedge->v->pos;
    Vector3d v_j = exFace->halfedge->next->v->pos;
    Vector3d v_k = exFace->halfedge->next->next->v->pos;

    Vector3d p = getBaryCoords(point, v_i, v_j, v_k);

    // trace to get corresponding intrinsic face
    tuple<shared_ptr<InFace>, Vector3d> intrinsic = pointQuery(exFace, p);

    Vector3d i(1, 0, 0);
    Vector3d j(0, 1, 0);
    Vector3d k(0, 0, 1);

    // paint extrinsic edges black
    double l_ij = exFace->halfedge->edge->length;
    double l_jk = exFace->halfedge->next->edge->length;
    double l_ki = exFace->halfedge->next->next->edge->length;

    double d_ij = distanceToEdge(p, i, j, l_ij, l_jk, l_ki);
    double d_jk = distanceToEdge(p, j, k, l_ij, l_jk, l_ki);
    double d_ki = distanceToEdge(p, k, i, l_ij, l_jk, l_ki);

    if (!_meanIntrinsicEdgeLength)
        computeMeanIntrinsicEdgeLength();
    double EXTRINSIC_OUTLINE = _meanIntrinsicEdgeLength/20.;
    double INTRINSIC_OUTLINE = EXTRINSIC_OUTLINE/3.;

    if (d_ij < EXTRINSIC_OUTLINE || d_jk < EXTRINSIC_OUTLINE || d_ki < EXTRINSIC_OUTLINE) {
        return Vector3d(0, 0, 0);
    }

    // paint intrinsic edges white
    Vector3d intrinsicBary = get<1>(intrinsic);
    shared_ptr<InFace> inFace = get<0>(intrinsic);

    l_ij = inFace->halfedge->edge->length;
    l_jk = inFace->halfedge->next->edge->length;
    l_ki = inFace->halfedge->next->next->edge->length;

    d_ij = distanceToEdge(intrinsicBary, i, j, l_ij, l_jk, l_ki);
    d_jk = distanceToEdge(intrinsicBary, j, k, l_ij, l_jk, l_ki);
    d_ki = distanceToEdge(intrinsicBary, k, i, l_ij, l_jk, l_ki);

    if (d_ij < INTRINSIC_OUTLINE || d_jk < INTRINSIC_OUTLINE || d_ki < INTRINSIC_OUTLINE) {
        return Vector3d(1, 1, 1);
    }

    // return color of that face
    return colors[_faceColors[inFace]];
}
