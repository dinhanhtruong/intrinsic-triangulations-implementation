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
// ij must already have a valid angle (bc angle of ik depends on angle of ij)
void SPmesh::updateSignpost(shared_ptr<InHalfedge> ij) {
    // given the 'right' edge of the face originiating from vertex i, get the left edge
    Vector3f ijEdge = getVPos(ij->twin->v) - getVPos(ij->v);
    shared_ptr<InHalfedge> ik = ij->next->next->twin;
    Vector3f ikEdge = getVPos(ik->v) - getVPos(ik->v);
    float theta_jk = getAngle(ijEdge, ikEdge); // (relative) angle of opposite edge jk (i.e. angle between ij and ik)
    ik->angle = ij->angle +  (2*M_PI * theta_jk)/ij->v->bigTheta; // set angle of jk relative to reference of this vertex
}


float SPmesh::getAngle(Vector3f v1, Vector3f v2) {
    return acos(v1.normalized().dot(v2.normalized()));
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
