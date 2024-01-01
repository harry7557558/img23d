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


void resamplePolygon(
    const std::vector<vec2> &boundary,
    std::vector<vec2> &output
) {
    // find largest angle
    int bn0 = (int)boundary.size();
    float maxc = -1.0f;
    int i0 = 0;
    for (int i = 0; i < bn0; i++) {
        vec2 v = boundary[i];
        vec2 v0 = boundary[(i+bn0-1)%bn0];
        vec2 v1 = boundary[(i+1)%bn0];
        float c = dot(v0-v,v1-v)/(length(v0-v)*length(v1-v));
        if (c > maxc)
            maxc = c, i0 = i;
    }

    // cubic interpolation
    const float tau = 0.5f;
    mat4 catmull_rom(
        0.0f, -tau, 2.0f*tau, -tau,
        1.0f, 0.0f, tau-3.0f, 2.0f-tau,
        0.0f, tau, 3.0f-2.0f*tau, tau-2.0f,
        0.0f, 0.0f, -tau, tau
    );
    const float quad_samples[4] = { .06943184420297371238,.33000947820757186759,.66999052179242813240,.93056815579702628761 };
    const float quad_weights[4] = { .17392742256872692868,.32607257743127307131,.32607257743127307131,.17392742256872692868 };
    std::vector<mat2x4> weights;
    weights.reserve(bn0);
    std::vector<float> length_psa;
    length_psa.reserve(bn0+1);
    length_psa.push_back(0.0f);
    for (int i = 0; i < bn0; i++) {
        // weights
        vec2 v0 = boundary[(i0+i+bn0-1)%bn0];
        vec2 v1 = boundary[(i0+i)%bn0];
        vec2 v2 = boundary[(i0+i+1)%bn0];
        vec2 v3 = boundary[(i0+i+2)%bn0];
        mat4x2 v(v0, v1, v2, v3);
        mat2x4 w = catmull_rom * transpose(v);
        weights.push_back(w);
        // length
        float l = 0.0f;
        for (int _ = 0; _ < 4; _++) {
            float t = quad_samples[_];
            vec4 u(0.0f, 1.0f, 2.0f*t, 3.0f*t*t);
            float dl = length(u*w);
            l += quad_weights[_] * dl;
        }
        length_psa.push_back(length_psa.back()+l);
    }
    for (int i = 0; i < bn0; i++)
        assert(length_psa[i] < length_psa[i+1]);
    float clength = length_psa.back();
    // printf("%f / %d = %f\n", length, bn0, length/bn0);
    auto curve = [&](float s) -> vec2 {
        s = clamp(s, 0.0f, 0.999999f*clength);
        int i = std::upper_bound(length_psa.begin(), length_psa.end()-1, s) - length_psa.begin() - 1;
        float t = (s-length_psa[i])/(length_psa[i+1]-length_psa[i]);
        assert(t >= 0.0 && t <= 1.0);
        t = clamp(t, 0.0f, 1.0f);
        return vec4(1.0f, t, t*t, t*t*t) * weights[i];
    };

    // resample
    int bn = std::max(bn0/3, std::min(bn0, 6));
    float dl = clength / bn0;
    float tol = 0.1f * dl;
    float ltol = 0.01f;  // sqrt(4/sqrt(3)*0.001)=0.05 for equilateral
    auto segmentSdf = [](vec2 p, vec2 a, vec2 b) {
        vec2 ba = b - a, pa = p - a;
        float h = dot(pa, ba) / dot(ba, ba);
        vec2 dp = pa - clamp(h, 0.0f, 1.0f) * ba;
        return length(dp);
    };
    auto errorBetween = [&](float s1, float s2) {
        vec2 p1 = curve(s1);
        vec2 p2 = curve(s2);
        int ns = (int)((s2-s1)/dl+0.51f);
        float ds = (s2-s1) / ns;
        float err = 0.0;
        for (int i = 0; i < ns; i++) {
            float s = mix(s1, s2, (i+0.5f)/ns);
            vec2 p = curve(s);
            float sdf = segmentSdf(p, p1, p2);
            // err += sdf * ds;
            err = fmax(err, sdf);
        }
        // return err / (s2-s1);
        return err;
    };
    std::vector<vec2> stack;
    stack.push_back(vec2(0.0f, clength));
    stack.push_back(vec2(0.0f, 0.5f*clength));  // prevent instant termination
    output.clear();
    while (!stack.empty()) {
        vec2 s = stack.back();
        float err = errorBetween(s.x, s.y);
        float l = length(curve(s.y) - curve(s.x));
        if (err > 9.0f * tol || l > 3.0f * ltol) {
            stack.push_back(vec2(s.x, 0.5f*(s.x+s.y)));
            continue;
        }
        output.push_back(curve(s.x));
        if (err > 4.0f * tol || l > 2.0f * ltol) {
            output.push_back(curve(mix(s.x, s.y, 1.0f/3.0f)));
            output.push_back(curve(mix(s.x, s.y, 2.0f/3.0f)));
        }
        else if (err > tol) {
            output.push_back(curve(0.5f*(s.x+s.y)));
        }
        while (!stack.empty() && s.y == stack.back().y)
            stack.pop_back();
        if (!stack.empty())
            stack.push_back(vec2(s.y, stack.back().y));
    }
}


float polygonArea(const std::vector<vec2>& polygon) {
    int n = (int)polygon.size();
    float area2 = 0.0f;
    for (int i = 0; i < n; i++)
        area2 += determinant(mat2(polygon[i], polygon[(i+1)%n]));
    return 0.5f * area2;
}


std::vector<vec2> holeLocations(
    std::vector<std::vector<vec2>> boundary
) {
    int num_segments = 0;
    std::vector<std::vector<vec2>*> hole_polygons;
    for (auto polygon = boundary.begin(); polygon < boundary.end(); polygon++) {
        if (polygonArea(*polygon) < 0.0f)
            hole_polygons.push_back(&(*polygon));
        num_segments += (int)polygon->size();
    }
    printf("%d contours, %d holes, %d segments\n",
        (int)boundary.size(), (int)hole_polygons.size(), num_segments);
    bool use_acceleration =
        1e-3f * (float)hole_polygons.size() * (float)num_segments >
        ((float)hole_polygons.size()+0.01f*(float)num_segments) * log((float)num_segments+1.0f);
    use_acceleration = true;

    // get a list of all segments
    std::vector<vec4> segments;
    segments.reserve(num_segments);
    for (auto polygon = boundary.begin(); polygon < boundary.end(); polygon++) {
        int n = (int)polygon->size();
        for (int i = 0; i < n; ++i) {
            vec2 v1(polygon->at(i));
            vec2 v2(polygon->at((i+1)%n));
            if (v2.y < v1.y) std::swap(v1, v2);
            segments.push_back(vec4(v1, v2));
        }
    }

    // build acceleration structure
    std::vector<float> sorted;
    for (auto polygon = boundary.begin(); polygon < boundary.end(); polygon++)
        for (vec2 p : *polygon)
            sorted.push_back(p.y);
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    int sn = (int)sorted.size();

    std::vector<std::vector<int>> intervals(sn);
    for (int si = 0; si < num_segments; si++) {
        vec4 seg = segments[si];
        vec2 v1(seg.x, seg.y), v2(seg.z, seg.w);
        int i1 = std::lower_bound(sorted.begin(), sorted.end(), v1.y) - sorted.begin();
        int i2 = std::lower_bound(sorted.begin(), sorted.end(), v2.y) - sorted.begin();
        assert(i1 >= 0 && i1 < sn);
        assert(i2 >= 0 && i2 < sn);
        assert(sorted[i1] == v1.y);
        assert(sorted[i2] == v2.y);
        for (int i = std::min(i1,i2); i < std::max(i1,i2); i++) {
            intervals[i].push_back(si);
        }
    }

    // random points + test if outside boundary
    std::vector<vec2> holes;
    for (const std::vector<vec2>* polygon : hole_polygons) {
        int pn = (int)polygon->size();
        srand(0);
        for (int guess = 0; guess < 20; guess++) {
            // generate random point
            int i = (int)((float)rand()/(float)RAND_MAX * (float)pn);
            int j = (i+1)%pn;
            vec2 dp = polygon->at(j)-polygon->at(i);
            vec2 dn = vec2(-dp.y, dp.x);
            float u = (float)rand()/(float)RAND_MAX;
            float v = (float)rand()/(float)RAND_MAX;
            vec2 p = polygon->at(i)+(0.3f+0.4f*u)*dp - (0.0f+0.4f*v*v)*dn;

            // add hole if the point is outside the boundary
            int crossings = 0;
            #if 0
            for (vec4 seg : segments) {
                vec2 v1(seg.x, seg.y), v2(seg.z, seg.w);
                if ((v1.y > p.y) != (v2.y > p.y) &&
                    (p.x - v1.x) * (v2.y - v1.y) < (v2.x - v1.x) * (p.y - v1.y)) {
                    crossings++;
                }
            }
            #else
            int ii = std::upper_bound(sorted.begin(), sorted.end(), p.y) - sorted.begin() - 1;
            if (ii < 0 || ii >= sn)
                continue;
            for (int si : intervals[ii]) {
                vec4 seg = segments[si];
                vec2 v1(seg.x, seg.y), v2(seg.z, seg.w);
                assert(v1.y <= v2.y);
                if ((v1.y > p.y) != (v2.y > p.y) &&
                    (p.x - v1.x) * (v2.y - v1.y) < (v2.x - v1.x) * (p.y - v1.y)) {
                    crossings++;
                }
            }
            #endif
            if (crossings % 2 == 0) {
                holes.push_back(p);
                break;
            }
        }
    }
    return holes;
}


void generateMesh(
    std::vector<vec2> verts,
    const std::vector<std::vector<int>> &boundary,
    std::vector<glm::vec2>& outputVerts,
    std::vector<glm::ivec3>& outputTrigs,
    bool resample
) {
    outputVerts.clear();
    outputTrigs.clear();
    if (boundary.empty()) return;

    float time0 = getTimePast();

    std::vector<vec2> points;
    std::vector<ivec2> segments;
    std::vector<std::vector<vec2>> boundary_r;
    for (std::vector<int> b0 : boundary) {
        std::vector<vec2> b;
        b.reserve(b0.size());
        for (int i : b0)
            b.push_back(verts[i]);
        std::vector<vec2> ps = b;
        if (resample)
            resamplePolygon(b, ps);
        int bn = (int)ps.size();
        if (bn < 3) continue;
        boundary_r.push_back(ps);

        // add points and segments
        int i0 = (int)points.size();
        for (int i = 0; i < bn; i++) {
            points.push_back(ps[i]);
            segments.push_back(ivec2(i0+i, i0+(i+1)%bn));
        }
    }

    float time1 = getTimePast();

    std::vector<vec2> holes = holeLocations(boundary_r);

    float time2 = getTimePast();

    triangleWrapper(points, segments, holes, outputVerts, outputTrigs);

    float time3 = getTimePast();
    printf("generateMesh: %.2g + %.2g + %.2g = %.2g secs\n",
        time1-time0, time2-time1, time3-time2, time3-time0);
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
