# Image to 3D Model

Drag and drop a transparent-background image in browser to get a rounded 3D shape as if the opaque part of the image is "inflated" like inflating a balloon. You can also apply transparency masks by per-pixel attributes like luma, chroma, and individual color channels. You can export inflated models to popular 3D formats like STL, OBJ, and GLB.

![](screenshot.png)

# How it works

This tool works by solving Poisson's equation. let the image be a Cartesian plane with $x$ and $y$ axes, the algorithm solves the equation $\nabla^2u+4=0$ with $u(x,y)=0$ at boundary of the region to inflate, and the 3D shape is defined by $z=\pm\sqrt{u}$, where $z$ is the axis perpendicular to the image plane.

The boundary is extracted using marching squares, and Poisson's equation is solved using finite element method on a mesh generated using [Triangle](https://www.cs.cmu.edu/~quake/triangle.html). Image texture is applied at the same $(x,y)$ coordinates.

# Tips

For good boundary extraction, PNG and SVG graphics art work the best. Noise and lossy compression artifacts can reduce speed and stability. If you are applying Img23d on a photo (of a print, sticker, shirt, etc.), take the photo in a bright and even lighting and avoid shadows and highlights. Keep the background clean and make sure details like thin lines are clearly visible without noise or blur.

If you are see missing/extra parts, which can happen when the boundaries have noise or fine details, try to uncheck the "simplify boundary" checkbox and see if the problem resolves.

To mask out white or black background, mask image by (reversed) luma. To mask out grayscale background and/or strokes but preserve colorful subjects, mask image by chroma.

For realistic results, the content of the image should not contain overlapping objects, 3D perspective effects, or strong lighting and shadowing. For example, a side photo of several waxed fruits likely won't produce an as visually pleasing 3D model as most binary-color logos do.

Drag-drop between browser windows and some SVG images may not work in some browsers. If you are experiencing issues in these cases, either switch to another browser (tested no problem on Google Chrome), save the image as a local file and drag-drop from file browser, or convert the image into a compatible format.

As far as I know, Img23d does not work on touchscreen devices. I lack motivation to add support since I don't have personal need for this feature. Feel free to submit an issue or pull request if you want to request or contribute to this feature.

# Licensing

Img23d software is distributed under the GPLv3 license. The [old version](https://github.com/harry7557558/Graphics/tree/master/modeling/img23d) and the [older version](https://github.com/harry7557558/Graphics/tree/master/modeling/png2obj) that do not depend on Triangle are also available under the MIT license. The default [hermit crab image](https://harry7557558.github.io/img23d/hermit_crab.svg) is an illustration I created for a school project in 2019, which is licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). The creator of Img23d does not claim right of 3D models generated from images that I do not own.
