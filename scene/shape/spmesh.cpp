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
    initSignpost();
}

Vector3f SPmesh::getVPos(shared_ptr<InVertex> v) {
    return v->exVertex->pos;
}

// updates the signpost angle of halfedge ik in triangle ijk
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
    h_ik->angle = h_ij->angle +  (2*M_PI * theta_i_jk)/h_ij->v->bigTheta;
}

// for a triangle with edge lengths l_ij, l_jk, l_ki, returns the interior angle at vertex i
float SPmesh::getAngleFromEdgeLengths(float l_ij, float l_jk, float l_ki) {
    // see fig. 9 of paper
    return acos( (pow(l_ij, 2) + pow(l_ki, 2) - pow(l_jk, 2)) / (2*l_ij*l_ki) );
}

// returns the angle between two 3-vectors in [0, pi]
float SPmesh::getAngle(Vector3f u, Vector3f v) {
    return acos(u.normalized().dot(v.normalized()));
}

// returns the ccw angle FROM vector u TO vector v in the range [0, 2*pi). Imagine fixing u as the x-axis in R^2 and going ccw to find v.
float SPmesh::argument(Vector2f u, Vector3f v) {
    // adapted from  https://stackoverflow.com/questions/40286650/how-to-get-the-anti-clockwise-angle-between-two-2d-vectors
    float angle = atan2(u[0]*v[1] - u[1]*v[0], u[0]*v[0] + u[1]*v[1]);
    if (angle < 0) {
        angle += 2*M_PI;
    }
    return angle;
}

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

    cout << "initialized signpost values" << endl;
}



void SPmesh::eraseTriangle(shared_ptr<InFace> tri) {
    // remove all half edges in face and set any references to them to null
    shared_ptr<InHalfedge> startHalfEdge = tri->halfedge;
    shared_ptr<InHalfedge> currHalfEdge = startHalfEdge;
    do {
        currHalfEdge->twin->twin = nullptr;
        // set the twin as the representative half edge of InEdge
        currHalfEdge->edge->halfedge = currHalfEdge->twin;
        // TODO remove inEdge entirely if it no longer contains half edges (e.g. if two adjacent faces are erased, middle edge is also erased)


        _halfedges.erase(currHalfEdge);
        currHalfEdge = currHalfEdge->next;
    } while(currHalfEdge != startHalfEdge);

    _faces.erase(tri);
}

// HELPER: returns the InEdge containing v0 and v1 as endpoints or nullptr if it doesn't exist
shared_ptr<InEdge> SPmesh::getEdge(shared_ptr<InVertex> v0, shared_ptr<InVertex> v1) const {
    // traverse 1-ring neighborhood of v0 and check incident vertices
    shared_ptr<InHalfedge> currHalfEdge = v0->halfedge;
    do {
        // return current half edge if it belongs to edge (v0, v1)
        if (currHalfEdge->next->v == v1)
            return currHalfEdge->edge;
        // advance clockwise
        currHalfEdge = currHalfEdge->twin->next;
    } while (currHalfEdge != v0->halfedge);
    return nullptr;
}

// HELPER: inserts an intrinsic triangle whose vertices are v0, v1, v2. Assumes v0,v1,v2 are provided in ccw order
// does not assume that there exist faces adjacent to the triangle (v0,v1,v2)
shared_ptr<InFace> SPmesh::insertTriangle(shared_ptr<InVertex> v0, shared_ptr<InVertex> v1, shared_ptr<InVertex> v2) {
    // make 3 half edges. link half edges to their existing twins (if any)
    shared_ptr<InHalfedge> newHalfEdge;
    auto vertices = vector{v0, v1, v2};
    for (int i = 0; i < 3; i++) {
        shared_ptr<InVertex> currVertex = vertices[i];
        shared_ptr<InEdge> edge = getEdge(vertices[i], vertices[(i+1) % 3]); // edge (i, i->next) containing halfedge to be constructed (can be null)
        shared_ptr<InHalfedge> twin = (edge == nullptr) ? nullptr : edge->halfedge; // twin of halfedge to be constructed (can be null)
        newHalfEdge = make_shared<InHalfedge>(InHalfedge{currVertex, nullptr, twin, edge, nullptr});
        if (edge == nullptr) {
            // edge does not exist: make one and connect it to the new halfEdge
            edge = make_shared<InEdge>(InEdge{newHalfEdge});
            newHalfEdge->edge = edge;
        } else {
            // link new half edge to its existing twin
            twin->twin = newHalfEdge;
        }
        if (currVertex->halfedge == nullptr)
            currVertex->halfedge = newHalfEdge; // set representative (don't touch existing representatives)
    }


    // set next pointers
    v0->halfedge->next = v1->halfedge;
    v1->halfedge->next = v2->halfedge;
    v2->halfedge->next = v0->halfedge;

    // make face
    shared_ptr<InFace> newFace = make_shared<InFace>(InFace{newHalfEdge}); // pick arbitrary half edge representative
    _faces.insert(newFace);

    return newFace;
}

float SPmesh::distance(float l_12, float l_23, float l_31, const Vector3f p, const Vector3f q) {
    Vector3f u = q - p;
    float d = -(
        pow(l_12, 2)*u[0]*u[1] +
        pow(l_23, 2)*u[1]*u[2] +
        pow(l_31, 2)*u[2]*u[0]
    );
    return sqrt(d);
}

// updates the signpost angles for every edge incident to the given vertex i.
// specifically, updates both the angle phi_ij and phi_ji for every edge (i,j)
//void SPmesh::updateVertex(shared_ptr<InVertex> i) {
//    // update incoming angles phi_ji
//    shared_ptr<InHalfedge> currHalfEdge = i->halfedge;
//    do {
//        updateSignpost();
//        currHalfEdge = currHalfEdge->twin->next;
//    } while (currHalfEdge != i->halfedge);

//    auto [exTriangle, barycentricCoords] = traceFromVertex();
//}

// inserts a new intrinsic vertex in the given intrinsic face ijk at the position specified using barycentric coords
// barycentric coords must be positive and sum to 1
void SPmesh::insertVertex(std::shared_ptr<InFace> face, Vector3f& barycentricCoords) {
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
    _verts.insert(p);

    // update signpost mesh connectivity
    /// 1) remove the existing intrinsic face (and its half-edges, but not its verts)
    eraseTriangle(face);
    /// 2) insert 3 new intrinsic faces around p and update half edge connectivity (will update signpost afterwards)
    std::shared_ptr<InFace> tri_ijp = insertTriangle(v_i, v_j, p);
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


    // update signpost angles
    updateVertex();
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

    cout << "loaded" << endl;
    validate();
    cout << "other" << endl;
    _exMesh.validate();
    cout << "validated" << endl;
}

//int SPmesh::getDegree(const shared_ptr<InVertex> &v) {
//    int i = 1;
//    shared_ptr<InHalfedge> curr = v->halfedge->twin->next;
//    while (curr != v->halfedge) {
//        curr = curr->twin->next;
//        i += 1;
//    }
//    return i;
//}

//Vector3f SPmesh::getNormal(Vector3f &v1, Vector3f &v2, Vector3f &v3) {
//    Vector3f ab = v3 - v1;
//    Vector3f ac = v2 - v1;
//    Vector3f normal = ac.cross(ab);
//    normal.normalize();
//    return normal;
//}

//float SPmesh::getArea(Vector3f &v1, Vector3f &v2, Vector3f &v3) {
//    Vector3f AC = v3 - v1;
//    Vector3f AB = v2 - v1;
//    return AB.cross(AC).norm()/2.f;
//}

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
        shared_ptr<InHalfedge> curr = v->halfedge->twin->next;
        while (curr != v->halfedge) {
            assert(curr->v == v);
            curr = curr->twin->next;
        }
        assert(curr == v->halfedge);
    }
}

void SPmesh::validate() {
    for (const shared_ptr<InHalfedge> &halfedge: _halfedges) {
        checkCircular(halfedge);
        checkTwin(halfedge);
    }
    checkFaces();
    checkVertices();
}
