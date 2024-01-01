# Image to 3D Model

Drag and drop a transparent-background image in browser to get a rounded 3D shape as if the opaque part of the image is "inflated" like inflating a balloon. You can also apply transparency masks by per-pixel attributes like luma, chroma, and individual color channels. You can export inflated models to popular 3D formats like STL, PLY, OBJ, and GLB.

![](screenshot.png)

# How it works

Img23d inflates images by solving Poisson's equation. Let the image be a Cartesian plane with $x$ and $y$ axes, it solves the equation $\nabla^2u+4=0$ with $u(x,y)=0$ at boundary of the region to inflate, and the 3D shape is defined by $z=\pm\sqrt{u}$, where $z$ is the axis perpendicular to the image plane. This method produces smooth surfaces and is consistent with translation, rotation, and uniform scaling.

The boundary is extracted using marching squares, and Poisson's equation is solved using finite element method on a mesh generated using [Triangle](https://www.cs.cmu.edu/~quake/triangle.html). Image texture is applied at the same $(x,y)$ coordinates.

# Tips

**Image quality.** Clean and losslessly compressed illustration work the best for good boundary extraction. Noise and lossy compression artifacts can greatly reduce speed and stability. If you are applying Img23d on a photo of a print, take the photo in a bright and even lighting and avoid shadows and highlights. Keep background clean and make sure details like thin lines are clear without noise or blur.

**Stability.** If you are see missing or extra parts or a completely flat 3D model, which can happen when the boundary has noise or fine details, try to select a slightly different threshold or check/uncheck the "simplify boundary" checkbox and see if the problem resolves.

**Background removal.** To remove image background, mask image by luma for dark background and reversed luma for light background. To mask out grayscale background and/or strokes but preserve colored subjects, mask image by chroma.

**Getting nice-looking results.** For realistic results, the content of the image should not contain overlapping objects, 3D perspective effects, or strong lighting and shadowing. For example, a side photo of several waxed fruits likely won't produce an as visually pleasing 3D model as most binary-color logos do.

**Issues with drag-drop.** Drag-drop between browser windows and some SVG images may not work in some browsers. If you are experiencing issues in these cases, either switch to another browser (tested no problem on Google Chrome), save the image as a local file and drag-drop from file browser, or convert the image into a compatible format.

**3D model format choice.** Among supported 3D model formats, only GLB preserves color, and STL does not support smooth shading. PLY and GLB files are relatively small in size. STL is commonly used in 3D printing, while GLB and OBJ have good compatibility with 3D animation software. If you don't have software to view exported 3D models locally, you can web search "3D model viewer", "GLB viewer", etc. to find online tools.

**Touchscreen support.** As far as I know, Img23d does not work on touchscreen devices. I lack motivation to add support since I don't have personal need for this feature. Feel free to submit an issue or a pull request if you want to request or contribute to this feature.

# Licensing

Img23d software is distributed under the GPLv3 license. The [old version](https://github.com/harry7557558/Graphics/tree/master/modeling/img23d) and the [older version](https://github.com/harry7557558/Graphics/tree/master/modeling/png2obj) that do not depend on Triangle are also available under the MIT license. The default [hermit crab image](https://harry7557558.github.io/img23d/hermit_crab.svg) is an illustration I created for a school project in 2019, which is licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). The creator of Img23d does not claim right of 3D models generated from images that I do not own.
