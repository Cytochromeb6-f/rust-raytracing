# Polynomial raytracing with Rust

This is my first ever project in Rust. I have mainly used this repository as a workspace to experiment with and learn the language. This is why the code struture is a bit messy and the comments are sparse. But I think the end result turned out okay.

## How it works (simplified)

The main players are the camera and a set of surfaces. As is common for ray tracing algorithms, the light rays travel in reverse.

### The camera

The camera is consists of a focal point and a rectangular array of pixels. Each pixel corresponds to a pixel in the final image and has a specific physical location in the scene. The program starts by shooting "light rays" from the focal point through each pixel and onward into the scene. The latest version uses multithreading and can simulate several ray paths in parallell.

### Polynomial surfaces

Each surface object has three variable polynomial $p(x,y,z)$ that specifies its shape. A light ray hits a surface when the $p$ evaluated at its position switches from positive to negative. The new direction of the ray after a reflection is calculated from the gradient $\nabla p$ using [Clifford algebra](https://github.com/Cytochromeb6-f/space_alg). Some surfaces are labeled to be light sources. When a ray hits a light source, it will multiply absorptive losses from earlier reflections with the color of the light source and add the that product to its corresponding pixel.

## Sample images

![A heart using a sixth degree polynomial with 22 terms](img\heart.png)
![The smaller cubes are made using p(x,y,z) = x^8 + y^8 + z^8 - c](img\rubikscube200.png)

More images in the folder [img](https://github.com/Cytochromeb6-f/rust-raytracing/tree/master/img).
