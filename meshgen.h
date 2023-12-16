#pragma once

#include "util.h"

#define ANSI_DECLARATORS
#define TRILIBRARY
// #define SINGLE
extern "C" {
#include "Triangle/triangle.c"
}


#include <vector>
#include <unordered_set>


void triangleWrapper(
    const std::vector<glm::vec2>& points,
    const std::vector<glm::ivec2>& segments,
    const std::vector<glm::vec2> &holes,
    std::vector<glm::vec2>& outputVerts,
    std::vector<glm::ivec3>& outputTrigs
) {

    triangulateio in, out;

    in.numberofpoints = (int)points.size();
    in.numberofpointattributes = 0;
    in.pointlist = (REAL*)malloc(points.size()*2*sizeof(REAL));
    in.pointmarkerlist = (int*)malloc(points.size()*sizeof(int));
    for (int i = 0; i < (int)points.size(); i++) {
        in.pointlist[2*i] = points[i].x;
        in.pointlist[2*i+1] = points[i].y;
        in.pointmarkerlist[i] = 1;
    }

    in.numberofsegments = (int)segments.size();
    in.segmentlist = (int*)malloc(2*segments.size()*sizeof(int));
    in.segmentmarkerlist = nullptr;
    for (int i = 0; i < (int)segments.size(); i++) {
        in.segmentlist[2*i] = i;
        in.segmentlist[2*i+1] = segments[i].y;
    }

    in.numberofholes = (int)holes.size();
    in.holelist = (REAL*)malloc(holes.size()*2*sizeof(REAL));
    for (int i = 0; i < (int)holes.size(); i++) {
        in.holelist[2*i] = holes[i].x;
        in.holelist[2*i+1] = holes[i].y;
    }

    in.numberofregions = 0;
    in.numberoftriangles = 0;
    in.numberofcorners = 0;
    in.numberoftriangleattributes = 0;
    in.trianglelist = nullptr;
    in.triangleattributelist = nullptr;

    out.pointlist = nullptr;
    out.trianglelist = nullptr;
    out.segmentlist = nullptr;

    triangulate("pzBqa.001Q", &in, &out, nullptr);

    outputVerts.clear();
    outputTrigs.clear();
    for (int i = 0; i < out.numberofpoints; i++) {
        const REAL *p = &out.pointlist[2*i];
        outputVerts.push_back(vec2(p[0], p[1]));
    }
    for (int i = 0; i < out.numberoftriangles; ++i) {
        const int *p = &out.trianglelist[3*i];
        outputTrigs.push_back(vec3(p[0], p[1], p[2]));
    }

    // Clean up
    free(in.pointlist);
    free(in.pointmarkerlist);
    free(in.segmentlist);
    free(in.holelist);
    free(out.pointlist);
    free(out.segmentlist);
    free(out.trianglelist);
}


void generateMesh(
    std::vector<vec2> verts,
    const std::vector<std::vector<int>> &boundary,
    std::vector<glm::vec2>& outputVerts,
    std::vector<glm::ivec3>& outputTrigs
) {
    outputVerts.clear();
    outputTrigs.clear();
    if (boundary.empty()) return;

    float time0 = getTimePast();

    std::vector<vec2> points;
    std::vector<ivec2> segments;
    std::vector<vec2> holes;
    for (std::vector<int> b0 : boundary) {
        // remove too close verts
        int bn0 = (int)b0.size();
        if (!(bn0 >= 3)) continue;
        float totLength = 0.0f;
        for (int i = 0; i < bn0; i++)
            totLength += length(verts[b0[(i+1)%bn0]]-verts[b0[i]]);
        float th = 0.2f * totLength / (float)bn0;
        std::vector<int> b;
        for (int i : b0) {
            if (!b.empty() && length(verts[i]-verts[b.back()]) < th)
                continue;
            b.push_back(i);
        }
        while (!b.empty() && length(verts[b[0]]-verts[b.back()]) < th)
            b.pop_back();
        int bn = (int)b.size();
        if (!(bn >= 6)) continue;

        // smooth boundary
        std::vector<vec2> ps(b.size());
        for (int _ = 0; _ < 2; _++) {
            // float k = _ % 2 == 0 ? 0.5f : -0.5f;
            float k = 0.2f;
            for (int i = 0; i < bn; i++) {
                vec2 p1 = verts[b[(i+1)%bn]];
                vec2 p0 = verts[b[(i+bn-1)%bn]];
                vec2 pc = verts[b[i]];
                ps[i] = pc + k * (p0+p1-2.0f*pc);
            }
            for (int i = 0; i < bn; i++)
                verts[b[i]] = ps[i];
        }

        // check if it's a hole
        float area2 = 0.0;
        for (int i = 0; i < bn; i++)
            area2 += determinant(mat2(ps[i], ps[(i+1)%bn]));
        if (area2 < 0.0) {
            vec2 p = 0.5f*(ps[1]+ps[0]);
            vec2 d = 0.5f*(ps[1]-ps[0]);
            holes.push_back(p-0.01f*vec2(-d.y,d.x));
        }

        // add points and segments
        int i0 = (int)points.size();
        for (int i = 0; i < bn; i++) {
            points.push_back(ps[i]);
            segments.push_back(ivec2(i0+i, i0+(i+1)%bn));
        }
        // printf("%d %d %f\n", bn0, bn, 0.5*area2);
    }

    float time1 = getTimePast();

    triangleWrapper(points, segments, holes, outputVerts, outputTrigs);

    float time2 = getTimePast();
    printf("generateMesh: %.2g + %.2g = %.2g secs\n",
        time1-time0, time2-time1, time2-time0);
}


void splitBridgeEdges(
    std::vector<glm::vec2> &verts,
    std::vector<glm::ivec3> &trigs
) {
    float time0 = getTimePast();

    // find boundary
    std::unordered_set<ivec2> bedges;
    for (ivec3 t : trigs) {
        for (int _ = 0; _ < 3; _++) {
            ivec2 e(t[(_+1)%3], t[_]);
            if (bedges.find(e) != bedges.end())
                bedges.erase(e);
            else bedges.insert(vec2(e.y, e.x));
        }
    }
    std::vector<bool> isBoundary(verts.size(), false);
    for (ivec2 e : bedges)
        isBoundary[e.x] = isBoundary[e.y] = true;

    float time1 = getTimePast();

    // find and split cross edges
    std::unordered_map<ivec2, int> cedges;
    std::vector<ivec3> trigs1;
    const int LUT[8][12] = {
        { 3,4,5, -1 },
        { 3,0,5, 0,4,5, -1 },
        { 3,4,1, 3,1,5, -1 },
        { 3,0,5, 5,0,1, 1,0,4, -1 },
        { 3,4,2, 2,4,5, -1 },
        { 3,0,2, 2,0,5, 5,0,4, -1 },
        { 3,4,1, 3,1,2, 2,1,5, -1 },
        { 3,0,2, 1,0,4, 2,1,5, 2,0,1 }
    };
    for (ivec3 t : trigs) {
        int tverts[6];
        int scase = 0;
        for (int _ = 0; _ < 3; _++) {
            tverts[_+3] = t[_];
            ivec2 e(t[_], t[(_+1)%3]);
            if (isBoundary[e.x] && isBoundary[e.y] &&
                (bedges.find(e) == bedges.end() &&
                 bedges.find(vec2(e.y,e.x)) == bedges.end())
            ) {
                if (e.x > e.y) std::swap(e.x, e.y);
                if (cedges.find(e) == cedges.end()) {
                    cedges[e] = (int)verts.size();
                    verts.push_back(0.5f*(verts[e.x]+verts[e.y]));
                }
                tverts[_] = cedges[e];
            }
            else tverts[_] = -1;
            scase |= int(tverts[_] != -1) << _;
        }
        const int* lut = LUT[scase];
        for (int i = 0; i < 12 && lut[i] != -1; i += 3) {
            trigs1.push_back(ivec3(
                tverts[lut[i]], tverts[lut[i+1]], tverts[lut[i+2]]
            ));
        }
    }
    trigs = trigs1;

    float time2 = getTimePast();
    printf("splitBridgeEdges: %.2g + %.2g = %.2g secs\n",
        time1-time0, time2-time1, time2-time0);
}
