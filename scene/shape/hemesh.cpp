#include "hemesh.h"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include <assert.h>

using namespace Eigen;
using namespace std;

void HEmesh::init(int numV, int numF)
{
    _verts = unordered_set<shared_ptr<ExVertex>>();
    _edges = unordered_set<shared_ptr<ExEdge>>();
    _faces = unordered_set<shared_ptr<ExFace>>();
    _halfedges = unordered_set<shared_ptr<ExHalfedge>>();

    _verts.reserve(numV);
    _faces.reserve(numF);
    _edges.reserve(numF * 1.5);
    _halfedges.reserve(numF * 3);
}

shared_ptr<ExVertex> HEmesh::makeVertex(shared_ptr<ExHalfedge> halfedge, Vector3f pos, shared_ptr<InVertex> inVertex) {
    shared_ptr<ExVertex> v = make_shared<ExVertex>(ExVertex{halfedge, pos, inVertex});
    _verts.insert(v);
    return v;
}

shared_ptr<ExHalfedge> HEmesh::makeHalfedge(shared_ptr<ExVertex> v, shared_ptr<ExHalfedge> next, shared_ptr<ExHalfedge> twin, shared_ptr<ExEdge> edge, shared_ptr<ExFace> face) {
    shared_ptr<ExHalfedge> halfedge = make_shared<ExHalfedge>(ExHalfedge{v, next, twin, edge, face});
    _halfedges.insert(halfedge);
    return halfedge;
}
shared_ptr<ExEdge> HEmesh::makeEdge(shared_ptr<ExHalfedge> halfedge) {
    shared_ptr<ExEdge> edge = make_shared<ExEdge>(ExEdge{halfedge});
    _edges.insert(edge);
    return edge;
}
shared_ptr<ExFace> HEmesh::makeFace(shared_ptr<ExHalfedge> halfedge) {
    shared_ptr<ExFace> face = make_shared<ExFace>(ExFace{halfedge});
    _faces.insert(face);
    return face;
}

//int HEmesh::getDegree(const shared_ptr<ExVertex> &v) {
//    int i = 1;
//    shared_ptr<ExHalfedge> curr = v->halfedge->twin->next;
//    while (curr != v->halfedge) {
//        curr = curr->twin->next;
//        i += 1;
//    }
//    return i;
//}

//Vector3f HEmesh::getNormal(Vector3f &v1, Vector3f &v2, Vector3f &v3) {
//    Vector3f ab = v3 - v1;
//    Vector3f ac = v2 - v1;
//    Vector3f normal = ac.cross(ab);
//    normal.normalize();
//    return normal;
//}

//float HEmesh::getArea(Vector3f &v1, Vector3f &v2, Vector3f &v3) {
//    Vector3f AC = v3 - v1;
//    Vector3f AB = v2 - v1;
//    return AB.cross(AC).norm()/2.f;
//}

void HEmesh::checkCircular(const shared_ptr<ExHalfedge> &halfedge) {
    assert(halfedge == halfedge->next->next->next);
    assert(halfedge->face == halfedge->next->face);
    assert(halfedge->face == halfedge->next->next->face);
}

void HEmesh::checkTwin(const shared_ptr<ExHalfedge> &halfedge) {
    assert(halfedge == halfedge->twin->twin);
    assert(halfedge->v == halfedge->twin->next->v);
    assert(halfedge->edge == halfedge->twin->edge);
}

void HEmesh::checkFaces() {
    for (const shared_ptr<ExFace> &face: _faces) {
        assert(face->halfedge->face == face);
        assert(face->halfedge->next->face == face);
        assert(face->halfedge->next->next->face == face);
    }
    for (const shared_ptr<ExEdge> &edge: _edges) {
        assert(edge->halfedge->edge == edge);
        assert(edge->halfedge->twin->edge == edge);
    }
    assert(_halfedges.size() == _faces.size() * 3);
    assert(_edges.size() * 2 == _faces.size() * 3);
}

void HEmesh::checkVertices() {
    for (const shared_ptr<ExVertex> &v: _verts) {
        shared_ptr<ExHalfedge> curr = v->halfedge->twin->next;
        while (curr != v->halfedge) {
            assert(curr->v == v);
            curr = curr->twin->next;
        }
        assert(curr == v->halfedge);
    }
}

void HEmesh::validate() {
    for (const shared_ptr<ExHalfedge> &halfedge: _halfedges) {
        checkCircular(halfedge);
        checkTwin(halfedge);
    }
    checkFaces();
    checkVertices();
}

void HEmesh::assignColors() {
    // adapted from https://www.geeksforgeeks.org/graph-coloring-set-2-greedy-algorithm/#
    _faceColors.reserve(_faces.size());
    int colorsUsed[4] = {false, false, false, false};

    auto face = _faces.begin();
    _faceColors[*face] = 0; // assign first color

    for (face = _faces.begin()++; face != _faces.end(); face++) {
        // flag colors already used by adjacent faces
        std::shared_ptr<ExFace> adj1 = (*face)->halfedge->twin->face;
        std::shared_ptr<ExFace> adj2 = (*face)->halfedge->next->twin->face;
        std::shared_ptr<ExFace> adj3 = (*face)->halfedge->next->next->twin->face;
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
