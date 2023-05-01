#include "spmesh.h"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>

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
   _vertices = vertices;
   _facesList    = faces;

    _verts = unordered_set<shared_ptr<InVertex>>();
    _edges = unordered_set<shared_ptr<InEdge>>();
    _faces = unordered_set<shared_ptr<InFace>>();
    _halfedges = unordered_set<shared_ptr<InHalfedge>>();

    _newVerts = unordered_set<shared_ptr<InVertex>>();
    _newEdges = unordered_set<shared_ptr<InEdge>>();
    _newFaces = unordered_set<shared_ptr<InFace>>();
    _newHalfedges = unordered_set<shared_ptr<InHalfedge>>();
    _newMiddleEdges = unordered_set<shared_ptr<InEdge>>();

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
        float bigTheta_i = 0.f;
        shared_ptr<InHalfedge> curr = v->halfedge->next->next->twin;
        Vector3f rightEdge = getVPos(v->halfedge->twin->v) - getVPos(v);
        do {
            Vector3f leftEdge = getVPos(curr->twin->v) - getVPos(curr->v);
            float angle = getAngle(rightEdge, leftEdge);
            bigTheta_i += angle;
            rightEdge = leftEdge;
            curr = curr->next->next->twin;
        } while (curr != v->halfedge->next->next->twin);
        v->bigTheta = bigTheta_i;

        float phi_ij0 = 0.f;
        curr = v->halfedge;
        curr->angle = 0.f;
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


    int numFlips = 8;
    for (int i = 0; i < numFlips; i++) {
        shared_ptr<InEdge> flippedEdge = nullptr;
        auto itr = _edges.begin();
        shared_ptr<InEdge> edge;
        while (flippedEdge == nullptr) {
            edge = *itr;
            flippedEdge = flipEdge(edge);
            itr++; // skip till find flippable edge
        }
        validate();
    }
//    shared_ptr<InEdge> edge = *_edges.begin();
//    shared_ptr<InEdge> flippedEdge = flipEdge(edge);

    validate();
}


// algo 1: updates the signpost angle of halfedge ik in triangle ijk
// ij must already have a valid angle (bc angle of ik depends on angle of ij). Edge lengths of ijk must also be valid
void SPmesh::updateSignpost(shared_ptr<InHalfedge> h_ij) {
    // given the 'right' halfedge ij of the face originiating from vertex i, get the left halfedge ik
    shared_ptr<InHalfedge> h_ik = h_ij->next->next->twin;
    // get edge lengths of triangle
    float l_ij = h_ij->edge->length;
    float l_jk = h_ij->next->edge->length;
    float l_ki = h_ij->next->next->edge->length;
    // update angle (phi) of halfedge ik
    float theta_i_jk = getAngleFromEdgeLengths(l_ij, l_jk, l_ki); // euclidean angle of opposite edge ik relative to ij (i.e. interior angle at vertex i between ij and ik)
    float phi_ik = h_ij->angle +  ((2*M_PI) * theta_i_jk)/h_ij->v->bigTheta; // == phi_ij + offset
    h_ik->angle = (phi_ik >= 2*M_PI) ? phi_ik - 2*M_PI : phi_ik; // constrain to [0, 2pi)
}

// angle should be the flat intrinsic angle (NOT phi)
std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3f, Eigen::Vector2f> SPmesh::traceFromIntrinsicVertex(std::shared_ptr<InVertex> v_i, float distance, float angle) {
    // TODO: convert to new base-finding algo in traceFromExtrinsic

    shared_ptr<ExHalfedge> base;
    Vector3f baryCoords;
    if (v_i->exVertex == nullptr) {
        baryCoords = v_i->barycentricPos;
        base = v_i->exFace->halfedge;
    } else {
        baryCoords = Vector3f(1.f, 0.f, 0.f);
        base = v_i->exVertex->halfedge;
        while (base->next->next->twin != v_i->exVertex->halfedge && base->next->next->twin->angle < angle) {
            base = base->next->next->twin;
        }
    }
    float a = fmod(angle - base->angle, (2*M_PI)) * v_i->exVertex->bigTheta/(2*M_PI);
    return traceVectorExtrinsic(base, baryCoords, distance, a);
}

// the angle argument is the normalized/projected phi angle in the range [0, 2pi]
std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3f, Eigen::Vector2f> SPmesh::traceVectorExtrinsic(std::shared_ptr<ExHalfedge> base, Eigen::Vector3f baryCoords, float distance, float angle) {
    shared_ptr<ExHalfedge> top = base->next->next;
    Vector2f fi = Vector2f(0.f, 0.f);
    Vector2f fj = Vector2f(0.f, base->edge->length);
    float theta = (top->twin->angle - base->angle) * base->v->bigTheta/(2*M_PI);
    Vector2f fk = top->edge->length * Vector2f(cos(theta), sin(theta));
    Vector3f bary = baryCoords;
    Vector3f dir;
    if (angle - M_PI_2 < 0.f && angle - M_PI_2 > M_PI) {
        dir = Vector3f(base->edge->length, atan(angle) * base->edge->length, 0.f);
    } else {
        dir = Vector3f(-base->edge->length, atan(angle) * -base->edge->length, 0.f);
    }
    dir.normalize();
    Matrix3f A;
    A << fi(0), fj(0), fk(0),
         fi(1), fj(1), fk(1),
         1.f, 1.f, 1.f;
    bool invertible = false;
    Matrix3f inverse;
    A.computeInverseWithCheck(inverse, invertible);
    assert(invertible);
    Vector3f baryDir = inverse * dir;

    float t = distance;
    int minT = 0;
    for (int i = 0; i < 3; i++) {
        if (baryDir(i) != 0.f) {
            float ti = -bary(i)/baryDir(i);
            if (ti > 0.f && ti < t) {
                t = ti;
                minT = i;
            }
        }
    }
//    cout << "t:" << t << " d:" << distance << endl;
//    cout << minT << endl;

    while (t < distance) {
        Vector3f edgeIntersectBary = bary + t * baryDir;
        Vector2f edge;
        // I'm not totally sure if this logic here is right. Based on the minT I think these are the edges we are
        // intersecting with but I'm actually not sure. an easy way to debug is to see which minT is being chose
        // technically in this case minT==1 should never be chosen since we will never be intersecting back with
        // the edge we're on. if you find this not to be true switch the edge and base assignments around based
        // on which minT isn't being chosen
        if (minT == 0) {
            edge = fk - fj;
            base = base->next->twin;
        } else if (minT == 1) {
            edge = fi - fk;
            base = base->next->next->twin;
        } else if (minT == 2) {
            edge = fj - fi;
            base = base->twin;
        }
        edge.normalize();
        top = base->next->next;
        Vector2f newfi = Vector2f(0.f, 0.f);
        Vector2f newfj = Vector2f(base->edge->length, 0.f);
        theta = (top->twin->angle - base->angle) * base->v->bigTheta/(2*M_PI);
        Vector2f newfk = top->edge->length * Vector2f(cos(theta), sin(theta));

        Vector2f tijk = Vector2f(-edge(1), edge(0));
        Vector2f newEdge = newfj - newfi;
        newEdge.normalize();
        Vector2f tnew = Vector2f(-newEdge(1), newEdge(0));
        Vector2f dir2d = Vector2f(dir(0), dir(1));
        assert (dir(2) == 0.f);

        Vector2f newDir = -((dir2d.dot(edge) * newEdge) + (dir2d.dot(tijk) * tnew));
        newDir.normalize();
        Vector2f p;
        if (edgeIntersectBary(0) == 0.f) {
            p = edgeIntersectBary(1) * newfj + edgeIntersectBary(2) * newfi;
        } else if (edgeIntersectBary(1) == 0.f) {
            p = edgeIntersectBary(0) * newfi + edgeIntersectBary(2) * newfj;
        } else if (edgeIntersectBary(2) == 0.f) {
            p = edgeIntersectBary(0) * newfj + edgeIntersectBary(1) * newfi;
        }

        fi = newfi;
        fj = newfj;
        fk = newfk;
        dir = Vector3f(newDir(0), newDir(1), 0.f);
        distance -= t;

        A << fi(0), fj(0), fk(0),
             fi(1), fj(1), fk(1),
             1.f, 1.f, 1.f;
        invertible = false;
        Matrix3f inverse;
        A.computeInverseWithCheck(inverse, invertible);
        assert(invertible);
        bary = inverse * Vector3f(p(0), p(1), 1.f);
        baryDir = inverse * dir;

        t = distance;
        minT = 0;
        for (int i = 0; i < 3; i++) {
            if (baryDir(i) != 0.f) {
                float ti = -bary(i)/baryDir(i);
                if (ti > 0.f && ti < t) {
                    t = ti;
                    minT = i;
                }
            }
        }
    }

    Vector3f newPointBary = bary + t * baryDir;
    if (base->next == base->face->halfedge) {
        newPointBary = Vector3f(newPointBary(1), newPointBary(2), newPointBary(0));
    } else if (base->next->next == base->face->halfedge) {
        newPointBary = Vector3f(newPointBary(2), newPointBary(0), newPointBary(1));
    }

    return make_tuple(base->face, newPointBary, Vector2f(baryDir(0), baryDir(1)));
}

// algo 2: see note in header file
// the angle argument is the normalized/projected phi angle in the range [0, 2pi] relative to the reference dir of v_i
std::tuple<std::shared_ptr<InFace>, Eigen::Vector3f> SPmesh::traceFromExtrinsicVertex(std::shared_ptr<ExVertex> v_i, float distance, float traceAngle) {
    shared_ptr<InHalfedge> base = v_i->inVertex->halfedge;

    // find the base vector whose source is at b_i and whose intrinsic triangle contains the dir to trace along (defined by the angle argument)
    // also find the flat/unprojected trace angle relative to base
    shared_ptr<InHalfedge> right;
    shared_ptr<InHalfedge> left;
    // in the general case, loop will terminate when the trace angle is between the phi angles of the left and right halfedges of the curr triangle
    do {
        right = base; // right halfedge in perspective of v_i
        left = right->next->next->twin; // one step ccw

        // edge case: reference DIRECTION (0 degrees) is inside the current triangle => left halfedge's phi < right halfedge's phi due to wrap-around
        if (left->angle < right->angle) {
            // target/trace direction is contained in this triangle iff traceAngle <= left XOR > right. Cannot simultaneously satisfy (right < traceAngle <= left) in this case.
            if (traceAngle > right->angle ||  traceAngle <= left->angle) {
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
    float traceAngleRelativeTheta = angleBetween(traceAngle, base->angle) * (v_i->inVertex->bigTheta/(2*M_PI));
    return traceVectorIntrinsic(base, Vector3f(1.f, 0.f, 0.f), distance, traceAngleRelativeTheta);
}


// angle should be the flat intrinsic angle (NOT phi)
std::tuple<std::shared_ptr<InFace>, Eigen::Vector3f> SPmesh::traceVectorIntrinsic(std::shared_ptr<InHalfedge> base, Eigen::Vector3f baryCoords, float distance, float angle) {
    shared_ptr<InHalfedge> top = base->next->next;
    Vector2f fi = Vector2f(0.f, 0.f);
    Vector2f fj = Vector2f(base->edge->length, 0.f);
    float theta = (top->twin->angle - base->angle) * base->v->bigTheta/(2*M_PI);
    Vector2f fk = top->edge->length * Vector2f(cos(theta), sin(theta));
    Vector3f bary = baryCoords;
    Vector3f dir = Vector3f(base->edge->length, tan(angle) * base->edge->length, 0.f);
    dir.normalize();
    Matrix3f A;
    A << fi(0), fj(0), fk(0),
         fi(1), fj(1), fk(1),
         1.f, 1.f, 1.f;
    bool invertible = false;
    Matrix3f inverse;
    A.computeInverseWithCheck(inverse, invertible);
    assert(invertible);
    Vector3f baryDir = inverse * dir;

    float t = distance;
    int minT = 0;
    for (int i = 0; i < 3; i++) {
        if (baryDir(i) != 0.f) {
            float ti = -bary(i)/baryDir(i);
            if (ti > 0.f && ti < t) {
                t = ti;
                minT = i;
            }
        }
    }
//    cout << "t:" << t << " d:" << distance << endl;
//    cout << minT << endl;

    while (t < distance) {
        Vector3f edgeIntersectBary = bary + t * baryDir;
        Vector2f edge;
        // I'm not totally sure if this logic here is right. Based on the minT I think these are the edges we are
        // intersecting with but I'm actually not sure. an easy way to debug is to see which minT is being chose
        // technically in this case minT==1 should never be chosen since we will never be intersecting back with
        // the edge we're on. if you find this not to be true switch the edge and base assignments around based
        // on which minT isn't being chosen
        if (minT == 0) {
            edge = fk - fj;
            base = base->next->twin;
        } else if (minT == 1) {
            edge = fi - fk;
            base = base->next->next->twin;
        } else if (minT == 2) {
            edge = fj - fi;
            base = base->twin;
        }
        edge.normalize();
        top = base->next->next;
        Vector2f newfi = Vector2f(0.f, 0.f);
        Vector2f newfj = Vector2f(base->edge->length, 0.f);
        theta = (top->twin->angle - base->angle) * base->v->bigTheta/(2*M_PI);
        Vector2f newfk = top->edge->length * Vector2f(cos(theta), sin(theta));

        Vector2f tijk = Vector2f(-edge(1), edge(0));
        Vector2f newEdge = newfj - newfi;
        newEdge.normalize();
        Vector2f tnew = Vector2f(-newEdge(1), newEdge(0));
        Vector2f dir2d = Vector2f(dir(0), dir(1));
        assert (dir(2) == 0.f);

        Vector2f newDir = -((dir2d.dot(edge) * newEdge) + (dir2d.dot(tijk) * tnew));
        newDir.normalize();
        Vector2f p;
        if (edgeIntersectBary(0) == 0.f) {
            p = edgeIntersectBary(1) * newfj + edgeIntersectBary(2) * newfi;
        } else if (edgeIntersectBary(1) == 0.f) {
            p = edgeIntersectBary(0) * newfi + edgeIntersectBary(2) * newfj;
        } else if (edgeIntersectBary(2) == 0.f) {
            p = edgeIntersectBary(0) * newfj + edgeIntersectBary(1) * newfi;
        }

        fi = newfi;
        fj = newfj;
        fk = newfk;
        dir = Vector3f(newDir(0), newDir(1), 0.f);
        distance -= t;

        A << fi(0), fj(0), fk(0),
             fi(1), fj(1), fk(1),
             1.f, 1.f, 1.f;
        invertible = false;
        Matrix3f inverse;
        A.computeInverseWithCheck(inverse, invertible);
        assert(invertible);
        bary = inverse * Vector3f(p(0), p(1), 1.f);
        baryDir = inverse * dir;

        t = distance;
        minT = 0;
        for (int i = 0; i < 3; i++) {
            if (baryDir(i) != 0.f) {
                float ti = -bary(i)/baryDir(i);
                if (ti > 0.f && ti < t) {
                    t = ti;
                    minT = i;
                }
            }
        }
//        cout << "t:" << t << " d:" << distance << endl;
    }

    Vector3f newPointBary = bary + t * baryDir;
    if (base->next == base->face->halfedge) {
        newPointBary = Vector3f(newPointBary(1), newPointBary(2), newPointBary(0));
    } else if (base->next->next == base->face->halfedge) {
        newPointBary = Vector3f(newPointBary(2), newPointBary(0), newPointBary(1));
    }

    return make_tuple(base->face, newPointBary);
}

// algo 3: updates the signpost angles for every edge incident to the given vertex vi.
// specifically, updates both the angle phi_ij and phi_ji for every edge (i,j)
void SPmesh::updateVertex(shared_ptr<InVertex> v_i) {
    cout << "updating vertex..." << endl;
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
    h_ij0->angle = argument(Vector2f(1,0), -uTransformed);
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
    if (getDegree(ij->halfedge->v) < 2 || getDegree(ij->halfedge->twin->v) < 2) {
        cout << "cannot flip: endpoint has degree < 2" << endl;
        return nullptr;
    }
    if (!edgeIsDelaunay(ij)) {
        cout << "cannot flip: edge is not delaunay" << endl;
        return nullptr;
    }
    std::shared_ptr<InHalfedge> lj = ij->halfedge->twin->next->next;
    std::shared_ptr<InHalfedge> ki = ij->halfedge->next->next;
    float l_ij = ij->length;
    float l_jk = ij->halfedge->next->edge->length;
    float l_ki = ki->edge->length;
    float l_il = ij->halfedge->twin->next->edge->length;
    float l_lj = lj->edge->length;

    float theta = getAngleFromEdgeLengths(l_ij, l_jk, l_ki) + getAngleFromEdgeLengths(l_il, l_lj, l_ij);

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
    float l_kl = baseLength(l_ki, l_il, theta);
    shared_ptr<InHalfedge> kl = tri_jkl->halfedge->next;
    kl->edge->length = l_kl;

    // update signposts
    /// update angle of HE lk using lj
    updateSignpost(lj);
    /// update angle of HE kl using ki
    updateSignpost(ki);

    cout << "flipped edge" << endl;
    return kl->edge;
}

// algo 6: returns the distance between points p and q inside a triangle with the three given edge lengths.
// p and q are barycentric coords.
float SPmesh::distance(float l_12, float l_23, float l_31, const Vector3f p, const Vector3f q) {
    Vector3f u = q - p;
    float d = -(
        pow(l_12, 2)*u[0]*u[1] +
        pow(l_23, 2)*u[1]*u[2] +
        pow(l_31, 2)*u[2]*u[0]
    );
    return sqrt(d);
}

// algo 7: inserts a new intrinsic vertex in the given intrinsic face ijk at the position specified using barycentric coords
// barycentric coords must be positive and sum to 1
void SPmesh::insertVertex(std::shared_ptr<InFace> face, Vector3f& barycentricCoords) {
    cout << "inserting vertex..." << endl;
    // get face's verts (ccw orientation)
    std::shared_ptr<InVertex> v_i = face->halfedge->v;
    std::shared_ptr<InVertex> v_j = face->halfedge->next->v;
    std::shared_ptr<InVertex> v_k = face->halfedge->next->next->v;
    // get face's edge lengths
    float l_ij = v_i->halfedge->edge->length;
    float l_jk = v_j->halfedge->edge->length;
    float l_ki = v_k->halfedge->edge->length;

    // make new vertex
    std::shared_ptr<InVertex> p = make_shared<InVertex>();
    p->exVertex = nullptr; // new intrinsic vertices do not correspond to any extrinsic vertex
    p->barycentricPos = barycentricCoords;
    p->bigTheta = (2*M_PI); // p is inside an intrinsic triangle => no curvature at p => no need to project; same angle sum as if p were in the plane (2pi)
    _verts.insert(p);

    // update signpost mesh connectivity
    /// 1) remove the existing intrinsic face
    _vertPairToEdge.clear(); // clear 1-ring neighborhood cache
    eraseTriangle(face);
    /// 2) insert 3 new intrinsic faces around p and update half edge connectivity (will update signposts afterwards)
    insertTriangle(v_i, v_j, p);
    insertTriangle(v_j, v_k, p);
    insertTriangle(v_k, v_i, p);
    validate();

    // compute and set edge lengths of all new edges incident to new vertex p
    // barycentric coordinates of vertices are: vi = (1,0,0),   vj = (0,1,0),   vk = (0,0,1).
    std::shared_ptr<InEdge> e_ip = p->halfedge->edge;
    e_ip->length = distance(l_ij, l_jk, l_ki, p->barycentricPos, Vector3f(1,0,0)); // distance along intrinsic (flattened) triangle from vi to p
    std::shared_ptr<InEdge> e_jp = p->halfedge->edge;
    e_jp->length = distance(l_ij, l_jk, l_ki, p->barycentricPos, Vector3f(0,1,0)); // vj to p
    std::shared_ptr<InEdge> e_kp = p->halfedge->edge;
    e_kp->length = distance(l_ij, l_jk, l_ki, p->barycentricPos, Vector3f(0,0,1)); // vk to p

    // update signpost angles around p and affected neighbors
    updateVertex(p);
}

// algo 9: given barycentric coords p for a point in triangle ijk returns polar coords r,phi of the vector from i to p
// phi is relative to edge ij
std::pair<float, float> SPmesh::vectorToPoint(float l_ij, float l_jk, float l_ki, const Eigen::Vector3f &i, const Eigen::Vector3f &j, const Eigen::Vector3f &p) {
    float r_pi = distance(l_ij, l_jk, l_ki, i, p); // distance from p to i
    float r_jp = distance(l_ij, l_jk, l_ki, j, p); // distance from j to p
    float phi_ip = getAngleFromEdgeLengths(l_ij, r_jp, r_pi); // angle from ij to ip
    return std::make_pair(r_pi, phi_ip);
}

// algo 10: moves inserted vertex i to p given by nonnegative barycentric coords inside triangle iab
void SPmesh::moveVertex(std::shared_ptr<InVertex> i, std::shared_ptr<InFace> iab, const Eigen::Vector3f &p) {
    std::shared_ptr<InHalfedge> ia = iab->halfedge;
    float l_ia = ia->edge->length;
    float l_ab = ia->next->edge->length;
    float l_bi = ia->next->next->edge->length;

    std::shared_ptr<InVertex> a = ia->next->v; /// TODO: (mehek) has no barycentric pos for initial vertices but can't hardcode coords

    std::pair<float, float> rPhi = vectorToPoint(l_ia, l_ab, l_bi, i->barycentricPos, Vector3f(0, 1, 0), p); // vector from i to p relative to ia
    float r = rPhi.first;
    float phi = rPhi.second; // angle of ip relative to ia
    // since i is inserted, all adjacent faces are coplanar and bigTheta_i=2*pi which means:
    // angle refHE->ip = angle ia->ip + angle refHE->ia
    phi = phi + ia->angle;

    // update lengths for adjacent neighbors
    std::shared_ptr<InHalfedge> ij = i->halfedge;
    do {
        float alpha = angleBetween(phi, ij->angle);
        float l_ij = ij->edge->length;
        /// note from Mehek: this looks like what baseLength does but for some reason they don't use it?
        float l_pj = sqrt(r * r + l_ij * l_ij - 2 * r * l_ij * cos(alpha)); // length from j to p (new pos for i)
        ij->edge->length = l_pj;
        ij = ij->next->next->twin;
    } while (ij != i->halfedge);

    // update signpost angles for adjacent halfedges
    updateVertex(i);
}

// algo 11: takes the barycentric coords of point on extrinsic face and returns the barycentric coordinates of that position on its intrinsic face
std::tuple<std::shared_ptr<InFace>, Eigen::Vector3f> SPmesh::pointQuery(std::shared_ptr<ExFace> xyz, Eigen::Vector3f& p) {
    std::shared_ptr<ExHalfedge> xy = xyz->halfedge;

    std::shared_ptr<ExVertex> x = xy->v;
    std::shared_ptr<ExVertex> y = xy->next->v;
    std::shared_ptr<ExVertex> z = xy->next->next->v;

    float l_xy = (y->pos - x->pos).norm();
    float l_yz = (z->pos - y->pos).norm();
    float l_zx = (x->pos - z->pos).norm();

    // vector from x to p relative to xy
    std::pair<float, float> rPhi = vectorToPoint(l_xy, l_yz, l_zx, Vector3f(1, 0, 0), Vector3f(0, 1, 0), p);

    float angle = (rPhi.second / x->bigTheta) * (2*M_PI) + xy->angle;

    return traceFromExtrinsicVertex(x, rPhi.first, angle);
}


//////////////////////////////////////////////
//////////////// HELPERS /////////////////////
//////////////////////////////////////////////

// HELPER :get the extrinsic position of an ORIGINAL intrinsic vertex
Vector3f SPmesh::getVPos(shared_ptr<InVertex> v) {
    return v->exVertex->pos;
}

// HELPER: for a triangle with edge lengths l_ij, l_jk, l_ki, returns the interior angle at vertex i
float SPmesh::getAngleFromEdgeLengths(float l_ij, float l_jk, float l_ki) {
    // see fig. 9 of paper
    return acos( (pow(l_ij, 2) + pow(l_ki, 2) - pow(l_jk, 2)) / (2*l_ij*l_ki) );
}

// HELPER: returns the angle between two 3-vectors in [0, pi]
float SPmesh::getAngle(Vector3f u, Vector3f v) {
    return acos(u.normalized().dot(v.normalized()));
}


// HELPER: returns the third side length of a triangle with side lengths a,b meeting at angle theta (law of cos)
float SPmesh::baseLength(float a, float b, float theta) {
    return sqrt(a * a + b * b - 2 * a * b * cos(theta));
}

// HELPER: returns smallest unsigned angle between the points on the unit circle given by angles a,b
float SPmesh::angleBetween(float a, float b) {
    float diff = abs(a - b); // get positive diff
    if (diff > M_PI) diff = (2*M_PI) - diff; // smaller angle should be in [0, pi]
    return diff;
}


// HELPER: returns the ccw angle FROM vector u TO vector v in the range [0, 2*pi).
//Imagine fixing u as the x-axis in R^2 and going ccw to find v.
float SPmesh::argument(Vector2f u, Vector2f v) {
    Vector2f a = u.normalized();
    Vector2f b = v.normalized();
    // adapted from  https://stackoverflow.com/questions/40286650/how-to-get-the-anti-clockwise-angle-between-two-2d-vectors
    float angle = atan2(a[0]*b[1] - a[1]*b[0], a[0]*b[0] + a[1]*b[1]);
    if (angle < 0) {
        angle += (2*M_PI);
    }
    return angle;
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

            // TODO: insert into _halfedges or _newHalfedges????

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

int SPmesh::getDegree(const shared_ptr<InVertex> &v) {
    int i = 0;
    shared_ptr<InHalfedge> curr = v->halfedge;
    do {
        curr = curr->twin->next;
        i += 1;
    } while (curr != v->halfedge);
    return i;
}

// returns true if the opposite angles on the adjacent faces of the given edge have sum <= pi. See diagram of equation 4.
bool SPmesh::edgeIsDelaunay(shared_ptr<InEdge> edge) {
    shared_ptr<InHalfedge> halfedge = edge->halfedge;
    shared_ptr<InHalfedge> twin = halfedge->twin;
    // get side lengths of the face lji incident to the representative halfedge of the given edge
    float l_ji = edge->length;
    float l_il = halfedge->next->edge->length;
    float l_lj = halfedge->next->next->edge->length;
    // get side lengths of the other face kij
    float l_ij = l_ji;
    float l_jk = twin->next->edge->length;
    float l_ki = twin->next->next->edge->length;

    float theta_l = getAngleFromEdgeLengths(l_lj, l_ji, l_il);
    float theta_k = getAngleFromEdgeLengths(l_ki, l_ij, l_jk);

    return theta_l + theta_k <= M_PI;
}

//Vector3f SPmesh::getNormal(Vector3f &v1, Vector3f &v2, Vector3f &v3) {
//    Vector3f ab = v3 - v1;
//    Vector3f ac = v2 - v1;
//    Vector3f normal = ac.cross(ab);
//    normal.normalize();
//    return normal;
//}


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
        assert(v->exVertex->inVertex == v);
        assert(v->bigTheta > 0);
    }
}

bool isEqual(float a, float b, float epsilon=0.001) {
    return abs(a-b) < epsilon;
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
    }
    for (const shared_ptr<InFace> &face: _faces) {
        // edges satisfy triangle inequalities
        float l_ij = face->halfedge->edge->length;
        float l_jk = face->halfedge->next->edge->length;
        float l_ki = face->halfedge->next->next->edge->length;
        assert(l_ij + l_jk >= l_ki);
        assert(l_jk + l_ki >= l_ij);
        assert(l_ki + l_ij >= l_jk);


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
    cout<<"passed validation"<<endl;
}


//////////////////////////////////////////////
///////////// VISUALIZATION //////////////////
//////////////////////////////////////////////

void SPmesh::assignColors() {
    // adapted from https://www.geeksforgeeks.org/graph-coloring-set-2-greedy-algorithm/#
    _faceColors.reserve(_faces.size());
    int colorsUsed[4] = {false, false, false, false};

    auto face = _faces.begin();
    _faceColors[*face] = 0; // assign first color

    for (face = _faces.begin()++; face != _faces.end(); face++) {
        // flag colors already used by adjacent faces
        std::shared_ptr<InFace> adj1 = (*face)->halfedge->twin->face;
        std::shared_ptr<InFace> adj2 = (*face)->halfedge->next->twin->face;
        std::shared_ptr<InFace> adj3 = (*face)->halfedge->next->next->twin->face;
        if (_faceColors.contains(adj1)) colorsUsed[_faceColors[adj1]] = true;
        if (_faceColors.contains(adj2)) colorsUsed[_faceColors[adj2]] = true;
        if (_faceColors.contains(adj3)) colorsUsed[_faceColors[adj3]] = true;

        // find lowest unused color
        int firstAvailableColor;
        for (firstAvailableColor = 0; firstAvailableColor < 4; firstAvailableColor++) {
            if (!colorsUsed[firstAvailableColor]) break;
        }

        // assign color to face
        _faceColors[*face] = firstAvailableColor;

        // reset color flags
        for (int i = 0; i < 4; i++) colorsUsed[i] = false;
    }
}

Eigen::Vector3f SPmesh::getBaryCoords(Eigen::Vector3f &p, Eigen::Vector3f &v1, Eigen::Vector3f &v2, Eigen::Vector3f &v3) {
    // courtesy of https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
    Vector3f a = v2 - v1;
    Vector3f b = v3 - v1;
    Vector3f c = p - v1;
    float d00 = a.dot(a);
    float d01 = a.dot(b);
    float d11 = b.dot(b);
    float d20 = c.dot(a);
    float d21 = c.dot(b);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return Vector3f(u, v, w);
}

int SPmesh::getColor(const Triangle* tri, Eigen::Vector3f point) {
    shared_ptr<ExFace> exFace = _exMesh.getExTriangle(tri->getIndex());

    // calculate barycentric coords on exFace
    Vector3f v1 = exFace->halfedge->v->pos;
    Vector3f v2 = exFace->halfedge->next->v->pos;
    Vector3f v3 = exFace->halfedge->next->next->v->pos;

    Vector3f p = getBaryCoords(point, v1, v2, v3);

    // trace to get corresponding intrinsic face
    tuple<shared_ptr<InFace>, Vector3f> intrinsic = pointQuery(exFace, p);

    // paint intrinsic edges
    Vector3f intrinsicBary = get<1>(intrinsic);
    if (intrinsicBary[0] < 0.025 || intrinsicBary[1] < 0.025 || intrinsicBary[2] < 0.025) {
        return -3;
    }

    // ignore color and draw outline for points close to exFace edges
    if (p[0] < 0.05 || p[1] < 0.05) {
        return -2;
    }
    else if (p[2] < 0.05) {
        return -1;
    }

    // return color of that face
    return _faceColors[get<0>(intrinsic)];
}
