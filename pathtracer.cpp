#include "pathtracer.h"

#include <iostream>

#include <Eigen/Dense>

#include <util/CS123Common.h>

using namespace Eigen;
//std::default_random_engine generator;
//std::uniform_real_distribution<float> uniformDist(0.0,1.0);

PathTracer::PathTracer(int width, int height)
    : m_width(width), m_height(height)
{
}

void PathTracer::traceScene(QRgb *imageData, Scene& scene) {
    std::vector<Vector3f> intensityValues(m_width * m_height);
    Matrix4f invViewMat = (scene.getCamera().getScaleMatrix() * scene.getCamera().getViewMatrix()).inverse();  // scale???
    int numSamplesPerPixel = USE_STRATIFIED_SUBPIXEL_SAMPLING ? pow(floor(sqrt(NUM_SAMPLES_PER_PIXEL)), 2) // get closest smaller perfect square
                                                              : NUM_SAMPLES_PER_PIXEL;
    #pragma omp parallel for
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            // trace N ray samples through each pixel and accumulate the radiance
            Vector3f currPixelRadiance = Vector3f::Zero();
            for(int i=0; i < numSamplesPerPixel; i++) {
                currPixelRadiance += tracePixel(x, y, scene, invViewMat, i, numSamplesPerPixel);
            }
            // store the averaged radiance
            int offset = x + (y * m_width);
            intensityValues[offset] = currPixelRadiance / numSamplesPerPixel;
        }
    }

    toneMap(imageData, intensityValues);
}

Ray PathTracer::sampleRandomRayThruPixel(Vector3f camOrigin, int imgRow, int imgCol) {
    // uniformly sample random point on the current pixel of the img plane by translating from the top left corner
    float randXOffset = generateRandFloat01(); // in [0,1]
    float randYOffset = -generateRandFloat01(); // in [-1,0]
    // direction is offset randomly from top left corner of pixel
    float imgPlaneDepth = DEPTH_OF_FIELD <= 0 ? 1 : DEPTH_OF_FIELD;
    Vector3f imagePlaneLoc(
                (2.f * (imgCol + randXOffset) / m_width) - 1,
                1 - (2.f * (imgRow + randYOffset) / m_height),
                -1 // must be on image plane (z=-1)
            );
    // rescale image plane based on its depth along look
    imagePlaneLoc *= imgPlaneDepth;
    Vector3f d = imagePlaneLoc - camOrigin;
    // construct camera-space ray
    return Ray(camOrigin, d.normalized());
}

Ray PathTracer::stratifiedSubpixelSampling(Vector3f camOrigin, int imgRow, int imgCol, int numCellsPerSide, int pixelCellIdx) {
    // basic idea: split pixel into numCellsPerSide*numCellsPerSide cells and sample one ray from each
    int cellRow = pixelCellIdx / numCellsPerSide;
    int cellCol = pixelCellIdx % numCellsPerSide;
    // assume pixels have width 1 => pixel cells have width 1/numCellsPerSide
    float cellWidth = 1.f/numCellsPerSide;
    float cellXOffset = cellWidth * cellCol;
    float cellYOffset = -cellWidth * cellRow;
    // randomly offset ray from the top left of the cell
    float randXOffset = generateRandFloat01() * cellWidth; // in [0,cellWidth]
    float randYOffset = -generateRandFloat01() * cellWidth; // in [-cellWidth,0]
    float imgPlaneDepth = DEPTH_OF_FIELD <= 0 ? 1 : DEPTH_OF_FIELD;
    Vector3f imagePlaneLoc(
                (2.f * (imgCol + cellXOffset + randXOffset)/m_width ) - 1,
                 1 - (2.f * (imgRow + cellYOffset + randYOffset)/m_height),
                -1 // must be on image plane (z=-1)
            );
    // rescale image plane based on its depth along look
    imagePlaneLoc *= imgPlaneDepth;
    Vector3f d = imagePlaneLoc - camOrigin;
    // construct camera-space ray
    return Ray(camOrigin, d.normalized());
}

Vector3f PathTracer::tracePixel(int x, int y, Scene& scene, const Matrix4f &invViewMatrix, int currentSampleIdx, int numSamples) {
    Vector3f origin(0, 0, 0);
    if (DEPTH_OF_FIELD > 0) {
        // scatter the eye location over a circular planar lens normal to look
        Vector2f circlePoint = sampleCirclePoint(0.35);
        origin = Vector3f(circlePoint.x(), circlePoint.y(), 0);
    }

    Ray r(origin, Vector3f::Zero());
    // get camera space ray from the integer image plane coords
    if (USE_STRATIFIED_SUBPIXEL_SAMPLING) {
        // shoot one ray through each sub-pixel 'cell'
        r = stratifiedSubpixelSampling(origin, y, x, sqrt(numSamples), currentSampleIdx);
    } else {
        r = sampleRandomRayThruPixel(origin, y, x);
    }

    // transform ray to world space
    r = r.transform(invViewMatrix);
    return traceRay(r, scene, true); // always count direct lighting toward camera
}

/**
 * @brief PathTracer::traceRay
 * @param r world space ray
 * @param scene to be rendered
 * @return incoming radiance Li(r_o, r_d) toward the origin of the given ray r from its first intersection point
 */
Vector3f PathTracer::traceRay(const Ray& r, Scene& scene, bool countEmitted) {
    IntersectionInfo insct;
    Ray rayOut(r);

    // find first intersection (if any) of ray w/ scene geometry
    if(scene.getIntersection(rayOut, &insct)) {
        const Triangle *tri = static_cast<const Triangle *>(insct.data); // cast intersected obj to tri since handling tri meshes
        SPmesh* spmesh = scene.getSPMesh();
        Vector3f color = spmesh->getColor(tri, insct.hit.cast<double>()).cast<float>();

        float kd = 0.5f;
        float ka = 0.3f;

        Vector3f normal = tri->getNormal(insct.hit); // world space normal at insct point
        Vector3f directionalLight(-1.f, -1.f, -1.f);

        return ka * color + kd * normal.dot(-directionalLight) * color;
    } else {
        return Vector3f(0, 0, 0);
    }
}

void PathTracer::toneMap(QRgb *imageData, std::vector<Vector3f> &intensityValues) {
    // find the maximum luminance value in the entire image
    float L_max = 0;
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);
            Vector3f radiance = intensityValues[offset].cwiseMax(0);
            float curr_L = RGBradianceToLuminance(radiance);
            L_max = std::max(L_max, curr_L);
        }
    }
    // write to the output image
    for(int y = 0; y < m_height; ++y) {
        for(int x = 0; x < m_width; ++x) {
            int offset = x + (y * m_width);

            // clip to [0,infty]
            Vector3f radiance = intensityValues[offset].cwiseMax(0);
            if (radiance.isZero()) {
                imageData[offset] = qRgb(0, 0, 0); // avoid divide by 0 issues later
                continue;
            }

            // convert to luminance using convex combination of color channels
            float L_world = RGBradianceToLuminance(radiance);
            // apply (extended) reinhard tone map to luminance
            float numerator = L_world * (1.f + (L_world / (L_max * L_max)));
            float L_display = numerator / (1.f + L_world);
            // scale original radiance by luminance ratios
            Vector3f tonemappedRadiance = (L_display / L_world) * radiance;

            // gamma correction
            tonemappedRadiance = tonemappedRadiance.array().pow(1.f/2);
            // clip to [0,1]
            tonemappedRadiance = tonemappedRadiance.cwiseMin(1);

            // scale to [0,255] for displaying
            tonemappedRadiance *= 255.f;

            // cast to int RGB vals
            imageData[offset] = qRgb((int)tonemappedRadiance[0], (int)tonemappedRadiance[1], (int)tonemappedRadiance[2]);
        }
    }
}

