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

    cout << "initialized signpost values" << endl;
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
    h_ik->angle = h_ij->angle +  (2*M_PI * theta_i_jk)/h_ij->v->bigTheta;
}

// algo 2: see note in header file
std::tuple<std::shared_ptr<ExFace>, Eigen::Vector3f, float> SPmesh::traceFromVertex(std::shared_ptr<InVertex> v_i, float distance, float angle) {
    return make_tuple(nullptr, Vector3f(0,0,0), 0);
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
    auto [exTriangle, barycentricCoords, uTransformed] = traceFromVertex(
                h_j0i->v,
                h_j0i->edge->length,
                h_j0i->angle); // <- u vector in section 3.2.3 in the local coordinate system at v_j0
    // get -u in coordinate system of the extrinsic triangle (180 rotation of u from the trace query)
    float e_ij0_dir = (uTransformed < M_PI) ? (uTransformed + M_PI) : (uTransformed - M_PI); // in [0, 2*pi)
    h_ij0->angle = e_ij0_dir;
    // set barycentric coords of vi in the extrinsic triangle
    v_i->barycentricPos = barycentricCoords;

    // 2) update remaining outgoing angles phi_ij in ccw order
    currHalfEdge = v_i->halfedge;
    while (currHalfEdge->next->next->twin != v_i->halfedge) {
        updateSignpost(currHalfEdge);
        currHalfEdge = currHalfEdge->next->next->twin;
    }
}

// algo 5: takes an edge ij with opposite vertices k,l and flips it to be kl
// replaces triangles ijk, jil with klj,lki
void SPmesh::flipEdge(std::shared_ptr<InEdge> ij) {
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
    eraseTriangle(ij->halfedge->face);
    eraseTriangle(ij->halfedge->twin->face);
    /// 2) insert 2 new intrinsic faces adjacent to flipped edge
    insertTriangle(vi, vl, vk);
    insertTriangle(vj, vk, vl);

    // get new length of flipped edge
    float l_kl = baseLength(l_ki, l_il, theta);

    // update signposts
    /// update angle of HE lk using lj
    updateSignpost(lj);
    /// update angle of HE kl using ki
    updateSignpost(ki);
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
    _verts.insert(p);

    // update signpost mesh connectivity
    /// 1) remove the existing intrinsic face (and its half-edges, but not its verts)
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
    // need to get vertex a but don't know what iab points to
    // loop until we get the HE starting at i (which is HE ia)
    std::shared_ptr<InHalfedge> ia = iab->halfedge;
    while (ia->v != i) {
        ia = ia->next;
    }
    std::shared_ptr<InVertex> a = ia->next->v;

    float l_ia = ia->edge->length;
    float l_ab = ia->next->edge->length;
    float l_bi = ia->next->next->edge->length;

    std::pair<float, float> rPhi = vectorToPoint(l_ia, l_ab, l_bi, i->barycentricPos, a->barycentricPos, p); // vector from i to p
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

// HELPER: returns the third side length of a triangle with side lengths a,b meeting at angle theta (law of cos)
float SPmesh::baseLength(float a, float b, float theta) {
    return sqrt(a * a + b * b + 2 * a * b * cos(theta));
}

/// note from Mehek: need to clarify what this is supposed to do. I implemented it based on the way it was described
/// but I don't think that matches what 3.3.3 is using it for
// HELPER: returns smallest unsigned angle between the points on the unit circle given by angles a,b
float SPmesh::angleBetween(float a, float b) {
    float diff = abs(a - b); // get positive diff
    if (diff > M_PI) diff = 2.f * M_PI - diff; // smaller angle should be in [0, pi]
    return diff;
}

// HELPER: returns the angle between two 3-vectors in [0, pi]
float SPmesh::getAngle(Vector3f u, Vector3f v) {
    return acos(u.normalized().dot(v.normalized()));
}

// HELPER: returns the ccw angle FROM vector u TO vector v in the range [0, 2*pi). Imagine fixing u as the x-axis in R^2 and going ccw to find v.
float SPmesh::argument(Vector2f u, Vector3f v) {
    // adapted from  https://stackoverflow.com/questions/40286650/how-to-get-the-anti-clockwise-angle-between-two-2d-vectors
    float angle = atan2(u[0]*v[1] - u[1]*v[0], u[0]*v[0] + u[1]*v[1]);
    if (angle < 0) {
        angle += 2*M_PI;
    }
    return angle;
}

// HELPER: removes the given intrinsic triangle.
// Also removes any incident InEdges that would no longer border two faces after the current removal
// Face's halfedges are preserved (unless the edge was removed, then the halfedges would also be removed)
void SPmesh::eraseTriangle(shared_ptr<InFace> tri) {
    shared_ptr<InHalfedge> startHalfEdge = tri->halfedge;
    shared_ptr<InHalfedge> currHalfEdge = startHalfEdge;
    do {
        currHalfEdge->face = nullptr;
        if (!currHalfEdge->twin->face) {
            // remove curr InEdge entirely since it is no longer adjacent to any faces (e.g. if two adjacent faces are erased, shared edge is also erased on the 2nd call)
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
    assert(_edges.contains(edge));
    assert(edge->halfedge);
    shared_ptr<InHalfedge> he_0 = edge->halfedge;
    shared_ptr<InHalfedge> he_1 = he_0 ? he_0->twin : nullptr;
    return (he_0->v == sourceVertex) ? he_0 : he_1;
}

// HELPER: inserts an intrinsic triangle whose vertices are v0, v1, v2. Assumes v0,v1,v2 are provided in ccw order.
// does not assume that there exist faces adjacent to the triangle (v0,v1,v2). Can be called in any order to insert multiple adjacent faces.
shared_ptr<InFace> SPmesh::insertTriangle(shared_ptr<InVertex> v0, shared_ptr<InVertex> v1, shared_ptr<InVertex> v2) {
    // make empty face
    shared_ptr<InFace> newFace = make_shared<InFace>(InFace{}); // pick arbitrary half edge representative
    _faces.insert(newFace);
    // construct new edges and their corresponding half-edges on this face if necessary. Link half edges of existing edges to their existing twins
    shared_ptr<InHalfedge> newHalfEdge;
    shared_ptr<InHalfedge> newFaceRepresentativeHalfEdge = nullptr;
    auto vertices = vector{v0, v1, v2};
    shared_ptr<InHalfedge> he_0;
    shared_ptr<InHalfedge> he_1;
    shared_ptr<InHalfedge> he_2;
    for (int i = 0; i < 3; i++) {
        shared_ptr<InVertex> currVertex = vertices[i];
        shared_ptr<InEdge> edge = getEdge(vertices[i], vertices[(i+1) % 3]); // edge (i, i->next) containing halfedge to be constructed (can be null)
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
                newHalfEdge->twin->twin = newHalfEdge;
                newHalfEdge->edge = edge;
                newHalfEdge->face = newFace;
            }
        } else {
            newHalfEdge = make_shared<InHalfedge>(InHalfedge{currVertex, nullptr, nullptr, nullptr, nullptr, -1});
            // edge does not exist: make one and connect it to the new halfEdge
            edge = make_shared<InEdge>(InEdge{newHalfEdge, -1}); // initialize length to -1 to indicate invalid
            newHalfEdge->edge = edge;
            newHalfEdge->face = newFace;
        }

        // pick an arbitrary representative halfedge for the new face
        if (!newFaceRepresentativeHalfEdge)
            newFaceRepresentativeHalfEdge = newHalfEdge;

        // record current half edge for setting next pointers later. Set he_i = current halfedge (whose base is at v_i)
        switch(i) {
            case 0: he_0 = newHalfEdge; break;
            case 1: he_1 = newHalfEdge; break;
            case 2: he_2 = newHalfEdge; break;
        }

/* old
        newHalfEdge = make_shared<InHalfedge>(InHalfedge{currVertex, nullptr, twin, edge, nullptr, -1}); // initialize angle to -1 to indicate invalid
        if (edge == nullptr) {
            // edge does not exist: make one and connect it to the new halfEdge
            edge = make_shared<InEdge>(InEdge{newHalfEdge, -1}); // initialize length to -1 to indicate invalid
            newHalfEdge->edge = edge;
        } else {
            // link new half edge to its existing twin
            twin->twin = newHalfEdge;
        }
        if (currVertex->halfedge == nullptr)
            currVertex->halfedge = newHalfEdge; // set representative (don't touch existing representatives)
          */
    }


    // set next pointers
    he_0->next = he_1;
    he_1->next = he_2;
    he_2->next = he_0;

    // set face representative
    newFace->halfedge = newFaceRepresentativeHalfEdge;
    return newFace;
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
        assert(halfedge->angle >= 0);
        assert(halfedge->angle < 2*M_PI);
        //assert(halfedge->edge->length > 0); // TODO: this is failing

        checkCircular(halfedge);
        checkTwin(halfedge);
    }
    checkFaces();
    checkVertices();
}
