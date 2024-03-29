// WASM, compile with Emscripten

#define SUPPRESS_ASSERT 1

#include <emscripten/emscripten.h>

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#error "Please compile this as a C++ file."
#endif


#include <cstdio>
#include <random>

#include "solver.h"
#include "render.h"

#include "meshgen.h"

#include "marching_squares.h"

#include "write_model.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"


namespace Img23d {

int channel = 0;
bool reverseAlpha = false;
uint8_t threshold = 127;
bool resampleBoundary = true;

std::string name;
uint8_t *data = nullptr;
uint8_t *dataWithoutAlpha = nullptr;
uint8_t *dataAlpha = nullptr;
int w = 0, h = 0;

bool showEdges = false;
bool smoothShading = true;
bool doubleSided = true;
bool showTexture = true;
float zsPos = 1.0f, zsNeg = 1.0f;
bool zClip = false;
float zcPos = 1.0f, zcNeg = 1.0f;

float getPixel(int x, int y) {
    x = clamp(x, 0, w-1);
    y = clamp(y, 0, h-1);
    int idx = y * w + x;
    return (float)data[4*idx+3] / 255.0f;
}

void computeAlpha() {
    if (dataAlpha)
        delete dataAlpha;
    dataAlpha = new uint8_t[w*h];
    uint8_t th = reverseAlpha ? ~threshold : threshold;
    bool hasNonempty = false, hasEmpty = false;
    for (int i = 0; i < w*h; i++) {
        uint8_t c[4] = {
            data[4*i], data[4*i+1], data[4*i+2],
            data[4*i+3]
        };
        dataAlpha[i] = 0;
        if (channel == 0)  // alpha
            dataAlpha[i] = c[3];
        if (channel == 1)  // luma
            dataAlpha[i] = (uint8_t)((
                299*(int)c[0] +
                587*(int)c[1] +
                114*(int)c[2]) / 1000);
        if (channel == 2)  // chroma
            dataAlpha[i] = std::max(std::max(c[0], c[1]), c[2])
                - std::min(std::min(c[0], c[1]), c[2]);
        //if (channel == 3)  // min/max
        //    dataAlpha[i] = (uint8_t)(
        //        ((int)std::min(std::min(c[0], c[1]), c[2]) * 256)
        //        / (int)std::max(std::max(c[0], c[1]), std::max(c[2], (uint8_t)1)));
        if (channel == 3)  // min
            dataAlpha[i] = std::min(std::min(c[0], c[1]), c[2]);
        if (channel == 4)  // max
            dataAlpha[i] = std::max(std::max(c[0], c[1]), c[2]);
        if (channel == 5)  // red
            dataAlpha[i] = c[0];
        if (channel == 6)  // green
            dataAlpha[i] = c[1];
        if (channel == 7)  // red
            dataAlpha[i] = c[2];
        if (channel != 0)
            dataAlpha[i] = (uint8_t)(reverseAlpha ?
                255 - (((255-(int)dataAlpha[i]) * (int)c[3]) >> 8) :
                ((int)dataAlpha[i] * (int)c[3]) >> 8);
        if (!reverseAlpha)
            dataAlpha[i] = ~dataAlpha[i];
        hasNonempty |= (dataAlpha[i] <= th);
        hasEmpty |= (dataAlpha[i] > th);
    }
    if (!hasNonempty) {
        emscripten_run_script("onError('warning: empty image')");
    }
    if (!hasEmpty) {
        if (Img23d::channel == 0)
            emscripten_run_script("onError('warning: image has no alpha')");
        else
            emscripten_run_script("onError('warning: completely filled image')");
    }
}

void computeImageWithoutAlpha() {
    if (dataWithoutAlpha)
        delete dataWithoutAlpha;
    dataWithoutAlpha = new uint8_t[4*w*h];
    memcpy(dataWithoutAlpha, data, 4*w*h);
    for (int i = 0; i < w*h; i++)
        dataWithoutAlpha[4*i+3] = data[4*i+3] > threshold ? 255 : 0;
    std::vector<ivec2> visited;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
            if (dataWithoutAlpha[4*(y*w+x)+3])
                visited.push_back(ivec2(x, y));
    while (!visited.empty()) {
        std::vector<ivec2> visited1;
        for (ivec2 ij : visited) {
            int i0 = ij.x, j0 = ij.y;
            for (int di = -1; di <= 1; di++)
                for (int dj = -1; dj <= 1; dj++) {
                    int i = i0+di, j = j0+dj;
                    if (i < 0 || j < 0 || i >= w || j >= h)
                        continue;
                    if (dataWithoutAlpha[4*(j*w+i)+3])
                        continue;
                    visited1.push_back(ivec2(i, j));
                    ((uint32_t*)dataWithoutAlpha)[j*w+i] =
                        ((uint32_t*)dataWithoutAlpha)[j0*w+i0];
                }
        }
        visited = visited1;
    }
}

DiscretizedModel<float, float> model;

}


DiscretizedModel<float, float> imageTo3D() {
    using namespace Img23d;

    float time0 = getTimePast();

    vec2 scale = vec2(w, h) / (float)fmax(w, h);
    vec2 bc = vec2(0), br = scale;
    renderModel.bound = br;

    std::vector<vec2> verts;
    std::vector<ivec3> trigs;

    Img23d::computeAlpha();
    uint8_t *alphas = new uint8_t[w*h];

#if 1
    const int R = 1;
    const int filter[2*R+1][2*R+1] = {
        { 1, 2, 1 }, { 2, 20, 2 }, { 1, 2, 1 }
    };
    // const int filter[2*R+1][2*R+1] = {
    //     { 1, 2, 1 }, { 2, 4, 2 }, { 1, 2, 1 }
    // };
    // const int filter[2*R+1][2*R+1] = {
    //     { 1,4,7,4,1 }, { 4,16,26,16,4 }, { 7,26,41,26,7 }, { 4,16,26,16,4 }, { 1,4,7,4,1 }
    // };
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int totw = 0, totv = 0;
            for (int dy = -R; dy <= R; dy++) {
                for (int dx = -R; dx <= R; dx++) {
                    if (y+dy >= 0 && y+dy < h && x+dx >= 0 && x+dx < w) {
                        totw += filter[dy+1][dx+1];
                        totv += int(dataAlpha[(y+dy)*w+(x+dx)]) * filter[dy+1][dx+1];
                    }
                }
            }
            alphas[((h-1-y)*w+x)] = uint8_t((totv+totw/2)/totw);
        }
    }
#endif

    float time1 = getTimePast();

    std::vector<std::vector<int>> boundary;
    uint8_t th = reverseAlpha ? ~threshold : threshold;
    marchingSquaresEdges(w, h, alphas, th, verts, boundary);
    for (int i = 0; i < (int)verts.size(); i++)
        verts[i] = (2.0f*verts[i]/vec2(w-1,h-1)-1.0f)*br;
    delete[] alphas;

    float time2 = getTimePast();

    generateMesh(verts, boundary, verts, trigs, resampleBoundary);
    splitBridgeEdges(verts, trigs);

    float time3 = getTimePast();

    DiscretizedModel<float, float> res = solveLaplacian(
        verts, std::vector<float>(verts.size(), 4.0f), trigs);
    for (int i = 0; i < res.N; i++)
        res.U[i] = 1.0f*sqrt(fmax(res.U[i], 0.0f));
    float maxu = 0.0; for (int i = 0; i < res.N; i++) maxu = fmax(maxu, res.U[i]);

    float time4 = getTimePast();
    printf("imageTo3D: %.2g + %.2g + %.2g + %.2g = %.2g secs\n",
        time1-time0, time2-time1, time3-time2, time4-time3, time4-time0);

    return res;
}



void prepareMesh(const DiscretizedModel<float, float> &model) {
    vec2 scale = vec2(Img23d::w, Img23d::h) / (float)fmax(Img23d::w, Img23d::h);

    int vn = (int)model.X.size();
    int fn = (int)model.E.size();

    // compute height threshold (based on RMS height)
    float totHeight = 0.0, totArea = 0.0;
    for (vec3 t : model.E) {
        float dA = abs(0.5*determinant(mat2(
            vec2(model.X[t[1]]-model.X[t[0]]),
            vec2(model.X[t[2]]-model.X[t[0]])
        )));
        float h[3];
        for (int _ = 0; _ < 3; _++)
            h[_] = model.U[t[_]] * model.U[t[_]];
        totArea += dA;
        totHeight += dA * (h[0]+h[1]+h[2])/3.0f;
    }
    float rmsh = sqrt(totHeight / totArea);
    float ztpos = rmsh * Img23d::zcPos;
    float ztneg = rmsh * Img23d::zcNeg;
    auto clipZ = [&](float z) {
        float a = z > 0.0f ? ztpos : ztneg;
        return z / sqrt(1.0+(z*z)/(a*a));
    };

    // determine which verts are repeated
    std::vector<bool> isBoundary(vn, false);
    std::vector<int> vmap;
    int vmapsum = 0;
    if (Img23d::doubleSided) {
        std::unordered_set<uint64_t> boundaryEdges;
        for (int fi = 0; fi < fn; fi++) {
            ivec3 t = model.E[fi];
            for (int _ = 0; _ < 3; _++) {
                ivec2 e(t[_], t[(_+1)%3]);
                if (e.x > e.y) std::swap(e.x, e.y);
                uint64_t i = *(uint64_t*)&e;
                if (boundaryEdges.find(i) == boundaryEdges.end())
                    boundaryEdges.insert(i);
                else boundaryEdges.erase(i);
            }
        }
        for (uint64_t e : boundaryEdges)
            isBoundary[(int)e] = isBoundary[(int)(e>>32)] = true;
        vmap.resize(vn);
        for (int i = 0; i < vn; i++) {
            vmap[i] = vmapsum;
            vmapsum += !isBoundary[i];
        }
    }

    int fn1 = Img23d::doubleSided ? 2*fn : fn;
    std::vector<int> vsmap;
    if (Img23d::smoothShading) {
        // vertices
        renderModel.vertices.resize(vn+vmapsum);
        renderModel.normals = std::vector<vec3>(vn+vmapsum, vec3(0));
        renderModel.texcoords.resize(vn+vmapsum);
        for (int i = 0; i < vn; i++) {
            renderModel.vertices[i] = vec3(model.X[i], model.U[i]);
            renderModel.texcoords[i] = 0.5f + 0.5f * vec2(1,-1) * model.X[i] / scale;
        }
        if (Img23d::doubleSided) {
            for (int i = 0; i < vn; i++) if (!isBoundary[i]) {
                renderModel.vertices[vn+vmap[i]] = renderModel.vertices[i] * vec3(1, 1, -Img23d::zsNeg);
                renderModel.texcoords[vn+vmap[i]] = renderModel.texcoords[i];
            }
        }
        for (int i = 0; i < vn; i++)
            renderModel.vertices[i].z *= Img23d::zsPos;
        if (Img23d::zClip) {
            for (int i = 0; i < int(renderModel.vertices.size()); i++)
                renderModel.vertices[i].z = clipZ(renderModel.vertices[i].z);
        }
        // triangles
        renderModel.indicesF.resize(fn1);
        for (int i = 0; i < fn; i++) {
            ivec3 t = model.E[i];
            renderModel.indicesF[i] = t;
            if (Img23d::doubleSided) {
                for (int _ = 0; _ < 3; _++)
                    if (!isBoundary[t[_]])
                        t[_] = vn + vmap[t[_]];
                std::swap(t[0], t[2]);
                renderModel.indicesF[fn+i] = t;
            }
        }
        // compute normal
        for (int i = 0; i < fn1; i++) {
            ivec3 t = renderModel.indicesF[i];
            vec3 n = cross(
                renderModel.vertices[t[1]]-renderModel.vertices[t[0]],
                renderModel.vertices[t[2]]-renderModel.vertices[t[0]]);
            n = normalize(n);
            for (int _ = 0; _ < 3; _++)
                renderModel.normals[t[_]] += n;
        }
        for (int i = 0; i < (int)renderModel.normals.size(); i++)
            renderModel.normals[i] = normalize(renderModel.normals[i]);
    }
    else {
        vsmap.resize(vn);
        // vertices
        renderModel.vertices.resize(fn1*3);
        renderModel.normals = std::vector<vec3>(fn1*3, vec3(0));
        renderModel.texcoords.resize(fn1*3);
        // triangles
        renderModel.indicesF.resize(fn1);
        for (int i = 0; i < fn1; i++) {
            ivec3 t = model.E[i%fn];
            vec3 vs[3];
            for (int _ = 0; _ < 3; _++) {
                int _i = Img23d::doubleSided && i >= fn ? 2-_ : _;
                vec3 v = vec3(model.X[t[_]], model.U[t[_]]);
                if (i < fn) vsmap[t[_]] = 3*i+_i;
                if (Img23d::doubleSided && i >= fn) v[2] *= -Img23d::zsNeg;
                else v[2] *= Img23d::zsPos;
                if (Img23d::zClip) v[2] = clipZ(v[2]);
                renderModel.vertices[3*i+_i] = v;
                renderModel.texcoords[3*i+_i] = 0.5f + 0.5f * vec2(1,-1) * vec2(v) / scale;
                vs[_i] = v;
            }
            vec3 n = normalize(cross(vs[1]-vs[0], vs[2]-vs[0]));
            for (int _ = 0; _ < 3; _++)
                renderModel.normals[3*i+_] = n;
            renderModel.indicesF[i] = ivec3(3*i, 3*i+1, 3*i+2);
        }
    }

    // edges
    renderModel.indicesE.clear();
    if (Img23d::showEdges) {
        std::vector<uint64_t> edges;
        for (int fi = 0; fi < fn; fi++) {
            ivec3 t = model.E[fi];
            for (int _ = 0; _ < 3; _++) {
                ivec2 e(t[_], t[(_+1)%3]);
                if (e.x > e.y) std::swap(e.x, e.y);
                edges.push_back(*(uint64_t*)&e);
            }
        }
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        for (uint64_t ei : edges) {
            ivec2 e = ivec2(int(ei), int(ei>>32));
            if (!Img23d::smoothShading)
                e[0] = vsmap[e[0]], e[1] = vsmap[e[1]];
            renderModel.indicesE.push_back(e);
            if (Img23d::doubleSided && Img23d::smoothShading) {
                for (int _ = 0; _ < 2; _++)
                    if (!isBoundary[e[_]])
                        e[_] = vn + vmap[e[_]];
                if (e[0] > vn || e[1] > vn)
                    renderModel.indicesE.push_back(e);
            }
            else if (Img23d::doubleSided) {
                for (int _ = 0; _ < 2; _++) {
                    if (e[_] % 3 == 0) e[_] += 2;
                    else if (e[_] % 3 == 2) e[_] -= 2;
                }
                renderModel.indicesE.push_back(e+3*fn);
            }
        }
    }

    // texture
    Img23d::computeImageWithoutAlpha();
    uint32_t white = 0xffe0e0e0;
    if (Img23d::showTexture)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
            Img23d::w, Img23d::h, 0, GL_RGBA, GL_UNSIGNED_BYTE,
            Img23d::dataWithoutAlpha);
    else
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
            1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE,
            &white);
}


void mainGUICallback() {
    // inside requestAnimationFrame
}

int main() {
    initWindow();
    emscripten_run_script("onReady()");

    mainGUI(mainGUICallback);
    return 0;
}


EXTERN EMSCRIPTEN_KEEPALIVE
void setChannel(int channel) {
    bool isUpdated = Img23d::channel == channel;
    Img23d::channel = channel;
    if (!isUpdated) {
        Img23d::model = imageTo3D();
        prepareMesh(Img23d::model);
        RenderParams::viewport->renderNeeded = true;
    }
}
EXTERN EMSCRIPTEN_KEEPALIVE
void setThreshold(uint8_t threshold) {
    bool isUpdated = Img23d::threshold == threshold;
    Img23d::threshold = threshold;
    if (!isUpdated) {
        Img23d::model = imageTo3D();
        prepareMesh(Img23d::model);
        RenderParams::viewport->renderNeeded = true;
    }
}
EXTERN EMSCRIPTEN_KEEPALIVE
void setAlphaReverse(bool reverse) {
    bool isUpdated = Img23d::reverseAlpha == reverse;
    Img23d::reverseAlpha = reverse;
    if (!isUpdated) {
        Img23d::model = imageTo3D();
        prepareMesh(Img23d::model);
        RenderParams::viewport->renderNeeded = true;
    }
}
EXTERN EMSCRIPTEN_KEEPALIVE
void setResampleBoundary(bool resample) {
    bool isUpdated = Img23d::resampleBoundary == resample;
    Img23d::resampleBoundary = resample;
    if (!isUpdated) {
        Img23d::model = imageTo3D();
        prepareMesh(Img23d::model);
        RenderParams::viewport->renderNeeded = true;
    }
}

EXTERN EMSCRIPTEN_KEEPALIVE
void setMeshEdge(bool edge) {
    bool isUpdated = Img23d::showEdges == edge;
    Img23d::showEdges = edge;
    if (!isUpdated) prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}
EXTERN EMSCRIPTEN_KEEPALIVE
void setMeshNormal(bool normal) {
    bool isUpdated = Img23d::smoothShading == normal;
    Img23d::smoothShading = normal;
    if (!isUpdated) prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}
EXTERN EMSCRIPTEN_KEEPALIVE
void setMeshDoubleSided(bool double_sided) {
    bool isUpdated = Img23d::doubleSided == double_sided;
    Img23d::doubleSided = double_sided;
    if (!isUpdated) prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}
EXTERN EMSCRIPTEN_KEEPALIVE
void setMeshTexture(bool texture) {
    bool isUpdated = Img23d::showTexture == texture;
    Img23d::showTexture = texture;
    if (!isUpdated) prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}

EXTERN EMSCRIPTEN_KEEPALIVE
void setZScale(float zspos, float zsneg) {
    zspos = zspos / (1.0f-zspos);
    zsneg = zsneg / (1.0f-zsneg);
    bool isUpdated = Img23d::zsPos == zspos && Img23d::zsNeg == zsneg;
    Img23d::zsPos = zspos;
    Img23d::zsNeg = zsneg;
    if (!isUpdated) prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}

EXTERN EMSCRIPTEN_KEEPALIVE
void setZClip(bool zc, float zcpos, float zcneg) {
    zcpos = zcpos / (1.0f-zcpos);
    zcneg = zcneg / (1.0f-zcneg);
    bool isUpdated = Img23d::zClip == zc &&
        Img23d::zcPos == zcpos && Img23d::zcNeg == zcneg;
    Img23d::zClip = zc;
    Img23d::zcPos = zcpos;
    Img23d::zcNeg = zcneg;
    if (!isUpdated) prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}

EXTERN EMSCRIPTEN_KEEPALIVE
void resizeWindow(int w, int h, float x, float y) {
    RenderParams::iResolution = glm::ivec2(w, h);
    RenderParams::viewport->center = glm::vec3(x, y, 0.0f);
    RenderParams::viewport->renderNeeded = true;
}


EXTERN EMSCRIPTEN_KEEPALIVE
void updateImage(const char* name, int w, int h, uint8_t *data) {
    printf("Image size: %d x %d\n", w, h);
    Img23d::name = name;
    Img23d::w = w;
    Img23d::h = h;
    if (Img23d::data)
        delete Img23d::data;
    Img23d::data = data;

    Img23d::model = imageTo3D();
    prepareMesh(Img23d::model);
    RenderParams::viewport->renderNeeded = true;
}


std::vector<uint8_t> fileBuffer;

EXTERN EMSCRIPTEN_KEEPALIVE
bool isModelEmpty() {
    return renderModel.vertices.empty()
        || renderModel.indicesF.empty();
}

EXTERN EMSCRIPTEN_KEEPALIVE
size_t getFileSize() {
    return fileBuffer.size();
}

EXTERN EMSCRIPTEN_KEEPALIVE
uint8_t* generateSTL() {
    fileBuffer = writeSTL_(
        renderModel.vertices,
        renderModel.indicesF
    );
    return fileBuffer.data();
}

EXTERN EMSCRIPTEN_KEEPALIVE
uint8_t* generatePLY() {
    if (Img23d::smoothShading) {
        fileBuffer = writePLY_(
            renderModel.vertices,
            renderModel.indicesF,
            renderModel.normals
        );
    }
    else {
        Img23d::smoothShading = true;
        prepareMesh(Img23d::model);
        fileBuffer = writePLY_(
            renderModel.vertices,
            renderModel.indicesF
        );
        Img23d::smoothShading = false;
        prepareMesh(Img23d::model);
    }
    return fileBuffer.data();
}

EXTERN EMSCRIPTEN_KEEPALIVE
uint8_t* generateOBJ() {
    if (Img23d::smoothShading) {
        fileBuffer = writeOBJ_(
            Img23d::name.data(),
            renderModel.vertices,
            renderModel.indicesF,
            renderModel.normals
        );
    }
    else {
        Img23d::smoothShading = true;
        prepareMesh(Img23d::model);
        fileBuffer = writeOBJ_(
            Img23d::name.data(),
            renderModel.vertices,
            renderModel.indicesF
        );
        Img23d::smoothShading = false;
        prepareMesh(Img23d::model);
    }
    return fileBuffer.data();
}

EXTERN EMSCRIPTEN_KEEPALIVE
uint8_t* generateGLB() {
    int nbytes;
    uint8_t *imgbytes_raw = stbi_write_png_to_mem(
        Img23d::dataWithoutAlpha, 4*Img23d::w, Img23d::w, Img23d::h, 4, &nbytes);
    std::vector<uint8_t> imgbytes(imgbytes_raw, imgbytes_raw+nbytes);
    free(imgbytes_raw);

    std::vector<vec2> emptyTexcoords;
    std::vector<uint8_t> emptyBytes;
    if (Img23d::smoothShading) {
        fileBuffer = writeGLB_(
            Img23d::name.data(),
            renderModel.vertices,
            renderModel.indicesF,
            renderModel.normals,
            Img23d::showTexture ? renderModel.texcoords : emptyTexcoords,
            Img23d::showTexture ? imgbytes : emptyBytes
        );
    }
    else {
        Img23d::smoothShading = true;
        prepareMesh(Img23d::model);
        fileBuffer = writeGLB_(
            Img23d::name.data(),
            renderModel.vertices,
            renderModel.indicesF,
            std::vector<vec3>(),
            Img23d::showTexture ? renderModel.texcoords : emptyTexcoords,
            Img23d::showTexture ? imgbytes : emptyBytes
        );
        Img23d::smoothShading = false;
        prepareMesh(Img23d::model);
    }
    return fileBuffer.data();
}
