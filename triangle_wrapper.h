#pragma once

#include "util.h"

#define ANSI_DECLARATORS
#define TRILIBRARY
#define SINGLE
extern "C" {
#include "Triangle/triangle.c"
}


#include <vector>


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


void triangleGenerateMesh(
    const std::vector<vec2> &verts,
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
        float th = 0.4f * totLength / (float)bn0;
        std::vector<int> b;
        for (int i : b0) {
            if (!b.empty() && length(verts[i]-verts[b.back()]) < th)
                continue;
            b.push_back(i);
        }
        while (!b.empty() && length(verts[b[0]]-verts[b.back()]) < th)
            b.pop_back();
        int bn = (int)b.size();
        if (!(bn >= 3)) continue;

        // check if it's a hole
        float area2 = 0.0;
        for (int i = 0; i < bn; i++)
            area2 += determinant(mat2(verts[b[i]], verts[b[(i+1)%bn]]));
        if (area2 < 0.0) {
            vec2 p = 0.5f*(verts[b[1]]+verts[b[0]]);
            vec2 d = 0.5f*(verts[b[1]]-verts[b[0]]);
            holes.push_back(p-0.01f*vec2(-d.y,d.x));
        }

        // add points and segments
        int i0 = (int)points.size();
        for (int i = 0; i < bn; i++) {
            points.push_back(verts[b[i]]);
            segments.push_back(ivec2(i0+i, i0+(i+1)%bn));
        }
        // printf("%d %d %f\n", bn0, bn, 0.5*area2);
    }

    float time1 = getTimePast();

    triangleWrapper(points, segments, holes, outputVerts, outputTrigs);

    float time2 = getTimePast();
    printf("triangleGenerateMesh: %.2g + %.2g = %.2g secs\n",
        time1-time0, time2-time1, time2-time0);
}
