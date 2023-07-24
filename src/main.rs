// use std::io;
use std::{fmt, f32::consts::{PI, TAU}, time::Instant};
use rand::random;

use png_encode_mini;

use space_alg::{Multivector, E_X, E_Y};

// use std::io::Cursor;
// use image::io::Reader as ImageReader;

// use std::collections::HashMap;
// use rand::Rng;


// Types, Integer and Float should have the same number of bits
type Float = f32;       // General float used in e.g. polynomial coefficients
type Integer = i32;     // General integer used for exponents in the monomials

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
struct Monomial {
    p_x: Integer,
    p_y: Integer,
    p_z: Integer
}
#[allow(dead_code)]
impl Monomial {
    fn new(p_x: Integer, p_y: Integer, p_z: Integer) -> Monomial {
        Monomial {p_x, p_y, p_z}
    }

    fn eval(&self, x: Float, y: Float, z: Float) -> Float {
        x.powi(self.p_x) * y.powi(self.p_y) * z.powi(self.p_z)
    }

    fn is_constant(&self) -> bool {
        (self.p_x==0) & (self.p_y==0) & (self.p_z==0)
    }

    fn derivative(&self, var: char) -> (Integer, Monomial) {
        match var {
            'x' => (self.p_x, Monomial::new(self.p_x - 1, self.p_y, self.p_z)),
            'y' => (self.p_y, Monomial::new(self.p_x, self.p_y - 1, self.p_z)),
            'z' => (self.p_z, Monomial::new(self.p_x, self.p_y, self.p_z - 1)),
            _ => (0, Monomial::new(0,0,0))
        }
    }

}
impl fmt::Display for Monomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut output: String = String::from("");

        for (var, expo) in [('x', self.p_x), ('y', self.p_y), ('z', self.p_z)] {
            if expo == 1 {
                output = format!("{}{}", output, var);
            } else if expo > 1 {
                output = format!("{}{}^{}", output, var, expo);
            } 
            
        }
        write!(f, "{output}")
    }
}





#[derive(PartialEq)]
struct Polynomial {
    terms: Vec<(Float, Monomial)>
}
#[allow(dead_code)]
impl Polynomial {
    fn new(terms: Vec<(Float, Monomial)>) -> Polynomial {
        let mut non_zero: Vec<(Float, Monomial)> = Vec::new();

        for (coeff, monom) in terms {
            if coeff == 0. {
                continue;
            }
            non_zero.push((coeff, monom));
        }
        
        Polynomial {terms: non_zero}
    }


    fn eval(&self, x: Float, y: Float, z: Float) -> Float {
        let mut sum: Float = 0.;
        for (coeff, monom) in &self.terms {
            sum += coeff*monom.eval(x,y,z)
        }
        sum
    }

    fn derivative(&self, var: char) -> Polynomial {
        let mut terms: Vec<(Float, Monomial)> = Vec::new();
        
        for (coeff1, monom) in &self.terms {
            let (coeff2, monom_prime) = monom.derivative(var);

            let coeff = coeff1*(coeff2 as Float);

            terms.push((coeff, monom_prime))
        }
        Polynomial::new(terms)
    }

}
impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        
        let mut output: String = String::from("");
        
        if self.terms.len() == 0 {
            output = 0.to_string();
        } else {
            // First term
            let (coeff, monom) = &self.terms[0];
            if monom.is_constant() {
                output = format!("{}{}{}", output, coeff, monom);
            } else if coeff == &1. {
                output = format!("{}{}", output, monom);
            } else if coeff == &-1. {
                output = format!("{}-{}", output, monom);
            } else {
                output = format!("{}{}{}", output, coeff, monom);
            }
            
            // Subsequent terms
            for (coeff, monom) in &self.terms[1..] {
                if coeff < &0. {
                    if coeff == &-1. {
                        output = format!("{} - {}", output, monom);
                    } else {
                        output = format!("{} - {}{}", output, -coeff, monom);
                    }
                } else if coeff > &0. {
                    if coeff == &1. {
                        output = format!("{} + {}", output, monom);
                    } else {
                        output = format!("{} + {}{}", output, coeff, monom);
                    }
                } 
            }
        }
        
        write!(f, "{output}")
    }

        
}

// #[derive(Copy)]
struct Solid {
    // A solid which has a surface given by the zero set of a polynomial. The interior is where the polynomial is negative.
    surface: Polynomial,
    gradient: [Polynomial; 3],
    origin: [Float; 3],
    bounding_box: [Float; 6]        // This is relative to the Solids own origin
}

#[allow(dead_code)]
impl Solid {
    fn new(surface: Polynomial, origin: [Float; 3], bounding_box: [Float; 6]) -> Solid {
        let gradient: [Polynomial; 3] = [surface.derivative('x'),
                                         surface.derivative('y'), 
                                         surface.derivative('z')];
        Solid {surface, gradient, origin, bounding_box}
    }

    fn new_wall(point: [Float; 3], normal: [Float; 3], bounding_box: [Float; 6]) -> Solid {
        // A half-space with a surface plane that is coincident with a given point and has a given normal direction

        let x1 = Monomial::new(1,0,0);
        let y1 = Monomial::new(0,1,0);
        let z1 = Monomial::new(0,0,1);

        let surface = Polynomial::new(vec![(normal[0], x1), (normal[1], y1), (normal[2], z1)]);

        Solid::new(surface, point, bounding_box)
    }

    fn new_sphere(center: [Float; 3], radius: Float) -> Solid {
        // A sphere with given a center and radius

        let bounding_box = [-radius, radius, -radius, radius, -radius, radius];

        let x2 = Monomial::new(2,0,0);
        let y2 = Monomial::new(0,2,0);
        let z2 = Monomial::new(0,0,2);
        let c = Monomial::new(0,0,0);

        let surface = Polynomial::new(vec![(1., x2), (1., y2), (1., z2), (-radius.powi(2), c)]);

        Solid::new(surface, center, bounding_box)
    }
    fn new_heart(origin: [Float;3]) -> Solid {

        let bounding_box = [-0.7, 0.7, -1.2, 1.2, -1., 1.3];

        let x6 = Monomial::new(6,0,0);
        let x4y2 = Monomial::new(4,2,0);
        let x4z2 = Monomial::new(4,0,2);
        let x4 = Monomial::new(4,0,0);
        let x2y4 = Monomial::new(2,4,0);
        let x2y2z2 = Monomial::new(2,2,2);
        let x2y2 = Monomial::new(2,2,0);
        let x2z4 = Monomial::new(2,0,4);
        let x2z3 = Monomial::new(2,0,3);
        let x2z2 = Monomial::new(2,0,2);
        let x2 = Monomial::new(2,0,0);
        let y6 = Monomial::new(0,6,0);
        let y4z2 = Monomial::new(0,4,2);
        let y4 = Monomial::new(0,4,0);
        let y2z4 = Monomial::new(0,2,4);
        let y2z3 = Monomial::new(0,2,3);
        let y2z2 = Monomial::new(0,2,2);
        let y2 = Monomial::new(0,2,0);
        let z6 = Monomial::new(0,0,6);
        let z4 = Monomial::new(0,0,4);
        let z2 = Monomial::new(0,0,2);
        let c = Monomial::new(0,0,0);

        let surface = Polynomial::new(vec![
           (1., y6), (6.75, x2y4), (3., y4z2), (-3., y4), (15.1875, x4y2), (13.5, x2y2z2), (-13.5, x2y2), 
            (3., y2z4), (-1., y2z3), (-6., y2z2), (3., y2), (11.3906, x6), (15.1875, x4z2), (-15.1875, x4), 
            (6.75, x2z4), (-0.1125, x2z3), (-13.5, x2z2), (6.75, x2), (1., z6), (-3., z4), (3., z2), (-1., c)
        ]);
        
        // let surface = Polynomial::new(vec![
        //     (1., x6), (6.75, x4y2), (3., x4z2), (-3., x4), (15.1875, x2y4), (13.5, x2y2z2), (-13.5, x2y2), 
        //     (3., x2z4), (-1., x2z3), (-6., x2z2), (3., x2), (11.3906, y6), (15.1875, y4z2), (-15.1875, y4), 
        //     (6.75, y2z4), (-0.1125, y2z3), (-13.5, y2z2), (6.75, y2), (1., z6), (-3., z4), (3., z2), (-1., c)
        // ]);

        Solid::new(surface, origin, bounding_box)
    }
    
    fn is_inside(&self, pos: Multivector) -> bool {

        let x = pos.comps()[1] - self.origin[0];
        let y = pos.comps()[2] - self.origin[1];
        let z = pos.comps()[3] - self.origin[2];

        match (x,y,z) {
            (x, _, _) if !(self.bounding_box[0]..self.bounding_box[1]).contains(&x) => false,
            (_, y, _) if !(self.bounding_box[2]..self.bounding_box[3]).contains(&y) => false,
            (_, _, z) if !(self.bounding_box[4]..self.bounding_box[5]).contains(&z) => false,
            _ => self.surface.eval(x,y,z) <= 0.
        }        
    }

    fn gradient_at(&self, pos: Multivector) -> Multivector {
        let x = pos.comps()[1] - self.origin[0];
        let y = pos.comps()[2] - self.origin[1];
        let z = pos.comps()[3] - self.origin[2];
        
        let g_x = self.gradient[0].eval(x,y,z);
        let g_y = self.gradient[1].eval(x,y,z);
        let g_z = self.gradient[2].eval(x,y,z);

        Multivector::new([0., g_x, g_y, g_z, 0., 0., 0., 0.])
    }

    fn refl(&self, pos: Multivector, vel: Multivector) -> Multivector {
        let normal = self.gradient_at(pos);

        -normal*vel*normal
    }

    fn refl_diffuse(&self, pos: Multivector) -> Multivector {    // Indepenent of incident velocity
        // Start with the velocity as the unit normal at pos
        let unit_normal = self.gradient_at(pos).unit_vector();

        // Rotate it down with a random angle given by Lambert's cosine law
        let vertical_plane = (Self::BLIND_SPOT^unit_normal).unit_bivector();
        let theta = random::<Float>().asin();
        let vel1 = unit_normal.scaled(theta.cos()) + (unit_normal<<vertical_plane).scaled(theta.sin());

        // Then rotate it by a random angle in the tangent plane of the surface
        let horizontal_plane = unit_normal.complement();
        let phi = random::<Float>()*TAU;    // Horizontal rotation angle from unifrom distribution
        vel1.scaled(phi.cos()) + (vel1<<horizontal_plane).scaled(phi.sin())

        // Will break if the normal is not linearly independent from BLIND_SPOT.
        // This is a probability 0 event for reals and still extremely unlikely when using floats.
    }
    const BLIND_SPOT: Multivector = Multivector::new_grade1([0.157631487, -0.697615504, 0.698916971]);

}
// fn get_custom_solid() -> Vec<char> {
//     let mut input_str: String = String::new();
//     io::stdin()
//         .read_line(&mut input_str)
//         .expect("F");
//     input_str.chars()
//         .collect()
// }



#[derive(Clone)]
struct Pixel {
    ray_pos: Vec<Multivector>,
    ray_vel: Vec<Multivector>,
    r: u32,                     // Index of this pixels red-value in the image
    g: u32,                     // etc.
    b: u32
}
impl Pixel {
    fn new(ray_pos: Vec<Multivector>, ray_vel: Vec<Multivector>, r:u32, g:u32, b:u32) -> Pixel {
        Pixel {ray_pos, ray_vel, r, g, b}
    }
}

#[allow(dead_code)]
struct Camera {
    focal_point: Multivector,
    im_w: u32,
    im_h: u32,
    pixels: Vec<Pixel>,     // Array which holds the pixel objects
    image: Vec<Float>       // The actual rgba-array

}
impl Camera {
    fn new(focal_point: [Float; 3], screen_center: [Float; 3], im_w: u32, im_h: u32, 
        pix_size: Float, pix_w: u8
    ) -> Camera {

        let now = Instant::now();          // Start of timed section

        let mut pixels: Vec<Pixel> = Vec::new();
        let mut image: Vec<Float> = Vec::new();

        let focal_point: Multivector = Multivector::new_grade1(focal_point);
        let screen_center: Multivector = Multivector::new_grade1(screen_center);
        let sub_pix_size = pix_size/(pix_w as Float);
        
        // Construct the screen geometry
        let camera_direction: Multivector = screen_center - focal_point;       // Vector in the direction of the camera
        let e_z = Multivector::new_grade1([0.,0.,1.]);       // Vector pointing straigt up 
        let vertical_plane = e_z^camera_direction;                // Bivector plane
        let e_w = vertical_plane.complement().unit_vector();      // Unit vector in the screens width-direction
        let horizontal_plane = e_w^camera_direction;              // Bivector plane
        let e_h = -horizontal_plane.complement().unit_vector();   // Unit vector in the screens height-direction

        
        // Populate the screen with pixels
        for i in 0..im_h {
            let y_pix = pix_size*((i as Float) - 0.5*(im_h-1) as Float);
            
            for j in 0..im_w {
                let x_pix = pix_size*((j as Float) - 0.5*(im_w-1) as Float);

                let r_index = 4*(im_w*i + j);   // Index of this pixels red-channel in the image
                image.append(&mut vec![0., 0., 0., 255.]);          // Black non-transparent starting image
                
                let mut ray_pos = Vec::new();
                let mut ray_vel = Vec::new();
                // Populate the pixels with start sites for rays
                for k in 0..pix_w {
                    let y = y_pix + sub_pix_size*((k as Float) - 0.5*(pix_w-1) as Float);
                    
                    for l in 0..pix_w {     // Pixels are squares so pix_h==pix_w
                        let x = x_pix + sub_pix_size*((l as Float) - 0.5*(pix_w-1) as Float);
                        
                        let position = screen_center + e_w.scaled(x) + e_h.scaled(y);

                        ray_pos.push(position);
                        ray_vel.push((position - focal_point).unit_vector())
                    }
                }

                pixels.push(Pixel::new(ray_pos, ray_vel, r_index, r_index+1, r_index+2));
            }
        }
        let elapsed_time = now.elapsed(); // End of timed section
        println!("In {} seconds: camera with {} pixels created.", elapsed_time.as_secs_f32(), im_w*im_h);
        
        Camera {focal_point, im_w, im_h, pixels, image}
    }

    fn photo(self, name: &str) -> () {
        let now = Instant::now();          // Start of timed section
        
        let mut f = std::fs::File::create(format!("{name}.png")).unwrap();

        let mut image_u8: Vec<u8> = Vec::new();
        for value in self.image {
            image_u8.push(value as u8)
        }
        let result = png_encode_mini::write_rgba_from_u8(&mut f, &image_u8, self.im_w, self.im_h);
        
        let elapsed_time = now.elapsed(); // End of timed section

        match result {
            Ok(_) => println!("In {} seconds: image written to \"{name}.png\".", elapsed_time.as_secs_f32()),
            Err(e) => println!("Error {:?} after {} seconds.", e, elapsed_time.as_secs())
        }

    }
}

struct Scene {
    solids: Vec<Solid>,
    mirrors: Vec<Solid>,
    lights: Vec<Solid>,
    camera: Camera
}
impl Scene {
    fn new(camera: Camera) -> Scene {
        Scene {solids: Vec::new(), mirrors: Vec::new(), lights: Vec::new(), camera}
    }

    fn add_solid(&mut self, solid: Solid) {
        self.solids.push(solid);
    }
    fn add_mirror(&mut self, mirror: Solid) {
        self.mirrors.push(mirror);
    }
    fn add_light(&mut self, light: Solid) {
        self.lights.push(light);
    }
    
    

    fn run(&mut self, time_step: Float, duration: Float, rays_per_sub_pixel: u8) {
        let now = Instant::now();          // Start of timed section
        // Create rays
        let mut rays: Vec<(&Pixel, Multivector, Multivector)> = Vec::new();
        for pix in &self.camera.pixels {
            // Pixel grid
            for (pos, vel) in std::iter::zip(&pix.ray_pos, &pix.ray_vel) {
                // Sub-pixel grid
                for _ in 0..rays_per_sub_pixel {
                    //
                    rays.push((pix, pos.to_owned(), vel.to_owned()));
                }
            }
            
        }
        let mut absorbed = Vec::new();

        let elapsed_time = now.elapsed(); // End of timed section
        println!("In {} seconds: {} rays produced.", elapsed_time.as_secs_f32(), rays.len());

        let now = Instant::now();          // Start of timed section
        // Simulate light transport
        let mut t = 0.;

        while t < duration  {
            
            'ray_evolution: for i in 0..rays.len() {
                // First test light hit
                for light in &self.lights {
                    if light.is_inside(rays[i].1) {
                        self.camera.image[rays[i].0.r as usize] += 3.;
                        self.camera.image[rays[i].0.g as usize] += 1.5;
                        self.camera.image[rays[i].0.b as usize] += 10.;
                        
                        // Mark for deletion
                        absorbed.push(i);

                        // Skip directly to next the ray
                        continue 'ray_evolution
                    }
                }

                
                'reflections: {
                    // Then test for diffuse surface hit
                    for solid in &self.solids {
                        if solid.is_inside(rays[i].1) {
                            rays[i].2 = solid.refl_diffuse(rays[i].1);
                            
                            // break 'reflections // Skip the other surfaces
                        }
                    }
                    // Then test for mirror reflection
                    for mirror in &self.mirrors {
                        if mirror.is_inside(rays[i].1) {
                            rays[i].2 = mirror.refl(rays[i].1, rays[i].2);
                            
                            break 'reflections // Skip the other mirrors
                        }
                    }
                }
                

                // Lastly move forward
                let vel = rays[i].2;
                rays[i].1 += vel.scaled(time_step);
            }

            // Deletes all the rays that hit a light
            absorbed.reverse();             //Iterate in reverse so hat swap_remove doesn't disturb itself
            for i in &absorbed[..] {
                // Removes in constant time: by replacing with the last element
                rays.swap_remove(*i);
            }
            absorbed.clear();
            

            // Increase time
            t += time_step
        }
        let elapsed_time = now.elapsed(); // End of timed section
        println!("In {} seconds: ray transport simulated." , elapsed_time.as_secs_f32());
    }
}

#[allow(dead_code)]
fn poly_testing() {
    let a = Monomial::new(0,0,0);
    let b = Monomial::new(1,0,0);
    let c = Monomial::new(1,1,1);
    let d = Monomial::new(1,0,3);
    let e = Monomial::new(1,2,0);

    let p = Polynomial::new(vec![(-1., a), (-2.3, b), (0., c), (-1., d), (1., e)]);

    // println!("{}",a.eval(2.,2.,2.));
    // println!("{}",b);
    println!("{}",p);
    println!("d/dw: {}",p.derivative('w'));
    println!("d/dx: {}",p.derivative('x'));
    println!("d/dy: {}",p.derivative('y'));
    println!("d/dz: {}",p.derivative('z'));
}

#[allow(dead_code)]
fn png_testing() {

    let mut f = std::fs::File::create("png_test.png").unwrap();

    // pixel order in image: [(bottom left...bottom right)...(top left... right)]
    let image_width: u32 = 3;
    let image_height: u32 = 2;
    /* let image = [
        // R     G     B     A
        0xff, 0x00, 0x00, 0xff,
        0x00, 0xff, 0x00, 0xff,
        0x00, 0x00, 0xff, 0xff,

        0x80, 0x00, 0x00, 0xff,
        0x00, 0x80, 0x00, 0xff,
        0x00, 0x00, 0x80, 0xff,
    ]; */
    let image = [
        // R     G     B     A
        0x20, 0x20, 0x20, 0xff,
        0x40, 0x40, 0x40, 0xff,
        0x60, 0x60, 0x60, 0xff,

        0x80, 0x80, 0x80, 0xff,
        0xa0, 0xa0, 0xa0, 0xff,
        0xc0, 0xc0, 0xc0, 0xff,
    ];

   

    match png_encode_mini::write_rgba_from_u8(&mut f, &image, image_width, image_height) {
        Ok(_) => println!("Written image!"),
        Err(e) => println!("Error {:?}", e),
    }
    
}

#[allow(dead_code)]
fn screen_testing() {
    let mut camera = Camera::new([0., 0., 0.], [1., 0., 0.],
        128, 128, 1e-2, 2
    );
    let mut blue = 0.;
    for pix in &camera.pixels {
        camera.image[pix.b as usize] = blue;
        blue += 0.01;
    }
    camera.photo("screen_test")

}

#[allow(dead_code)]
fn diffuse_testing() {
    let rep = 4_000_000 ;

    let mut normal = Vec::new();
    for _ in 0..rep {
        normal.push(Multivector::new_grade1([random(), random(), random()]));
    }

    
    let now = Instant::now();          // Start of timed section
    for i in 0..rep {
    
        // Construct an orthonormal basis for the tangent plane
        let mut e_1 = (E_Y^normal[i])*normal[i];
        if e_1.is_zero(){                                // Happens if normal==E_X
            e_1 = (E_X^normal[i])*normal[i];
        }
        e_1 = e_1.unit_vector();
        
        
        let unit_normal = normal[i].unit_vector();
        let e_2 = (e_1*unit_normal).complement();
        
        // Make a random vector from a circle distribution in the tangent plane.
        let phi = random::<Float>()*TAU;
        let rand_vect = e_1.scaled(phi.cos()) + e_2.scaled(phi.sin());

        // Make unit-bivector to rotate the normal along
        // let rot_plane = unit_normal^rand_vect;

        // Rotation angle from Lambert's cosine law
        let r = random::<Float>().asin();
        

        // Get the velocity by rotating the normal and scaling it to unit length to set speed to 1
        let _out = unit_normal.scaled(r.cos()) + (unit_normal<<(unit_normal^rand_vect)).scaled(r.sin());
    }
    let elapsed_time = now.elapsed(); // End of timed section
    println!("Method 1 took {} seconds." , elapsed_time.as_secs_f32());
    

    let now = Instant::now();          // Start of timed section
    for i in 0..rep {

        // Make a random vector from a spherical distribution. It is used to get the rotation direction
        let theta = random::<Float>()*PI;
        let phi = random::<Float>()*TAU;
        let rand_vect = Multivector::new_grade1([phi.cos()*theta.sin(), phi.sin()*theta.sin(), theta.cos()]);
        
        // Make unit-bivector to rotate the normal along
        let rot_plane = (normal[i]^rand_vect).unit_bivector();
        
        // Rotation angle from Lambert's cosine law
        let r = random::<Float>().asin();
        
        
        // Get the velocity by rotating the normal and scaling it to unit length to set speed to 1
        let _out = (normal[i].scaled(r.cos()) + normal[i]<<rot_plane).scaled(r.sin()).unit_vector();
    }
    let elapsed_time = now.elapsed(); // End of timed section
    println!("Method 2 took {} seconds." , elapsed_time.as_secs_f32());

    let blind_spot = Multivector::new_grade1([0.157631487, -0.697615504, 0.698916971]);
    let now = Instant::now();          // Start of timed section
    for i in 0..rep {
        
        let unit_normal =  normal[i].unit_vector();
        
        // Make a unit-vector with angle to the surface normal given by Lambert's cosine law
        let vertical_plane = (blind_spot^unit_normal).unit_bivector();
        let theta = random::<Float>().asin();
        let vel1 = unit_normal.scaled(theta.cos()) + (unit_normal<<vertical_plane).scaled(theta.sin());

        // Then rotate it by a random angle in the tangent plane of the surface
        let horizontal_plane = unit_normal.complement();
        let phi = random::<Float>()*TAU;    // Horizontal rotation angle from unifrom distribution
        
        let _out = vel1.scaled(phi.cos()) + (vel1<<horizontal_plane).scaled(phi.sin());
    }
    let elapsed_time = now.elapsed(); // End of timed section
    println!("Method 3 took {} seconds." , elapsed_time.as_secs_f32());
    

    let now = Instant::now();          // Start of timed section
    for i in 0..rep {
        
        let mut vel =  normal[i].unit_vector();
        
        // Make a unit-vector with angle to the surface normal given by Lambert's cosine law
        let vertical_plane = (blind_spot^vel).unit_bivector();
        let theta = random::<Float>().asin();
        vel = vel.scaled(theta.cos()) + (vel<<vertical_plane).scaled(theta.sin());

        // Then rotate it by a random angle in the tangent plane of the surface
        let horizontal_plane = vel.complement();
        let phi = random::<Float>()*TAU;    // Horizontal rotation angle from unifrom distribution
        
        let _out = vel.scaled(phi.cos()) + (vel<<horizontal_plane).scaled(phi.sin());
    }
    let elapsed_time = now.elapsed(); // End of timed section
    println!("Method 4 took {} seconds." , elapsed_time.as_secs_f32());
}

#[allow(dead_code)]
fn transport_testing() {
    let focal_point = [-0.9, 0., 0.];
    let screen_center = [0., 0., 0.];

    let pixel_size = 5e-3;
    let im_w = 256;                 // Number of rays created is im_w*im_h*pix_w^2*rays_per_sub_pixel rays
    let im_h =  128;
    let pix_w = 3;
    let rays_per_sub_pixel = 1;

    
    let duration = 20.;
    let time_step = 0.2;

    let camera = Camera::new(focal_point, screen_center, im_w, im_h, pixel_size, pix_w);
    let mut scene = Scene::new(camera);
    
    scene.add_light(Solid::new_sphere([3., 1.2, 0.1], 1.0));
    // scene.add_solid(Solid::new_wall([-1., 0., 0.], [1., 0., 0.]));


    // scene.add_solid(Solid::new_sphere([3., -1.2, 0.], 1.));
    // scene.add_solid(Solid::new_wall([0., 0., -1.], [0., 0., 1.]));
    // scene.add_solid(Solid::new_wall([0., 0., 1.5], [0., 0., -1.]));
    // scene.add_solid(Solid::new_wall([6., 0., 0.], [-1., 0., 0.]));
    
    
    scene.run(time_step, duration, rays_per_sub_pixel);

    scene.camera.photo("transport_test");
}

#[allow(dead_code)]
fn time_testing() {
    let focal_point = [-0.9, 0., 0.];
    let screen_center = [0., 0., 0.];

    let pixel_size = 5e-3;
    let im_w = 256;                 // Number of rays created is im_w*im_h*pix_w^2*rays_per_sub_pixel rays
    let im_h =  128;
    let pix_w = 2;
    let rays_per_sub_pixel = 1;

    
    let duration = 10.;
    let time_step = 0.5;

    let camera = Camera::new(focal_point, screen_center, im_w, im_h, pixel_size, pix_w);
    let mut scene = Scene::new(camera);
    
    scene.add_light(Solid::new_sphere([3., 0.,0.], 1.));

    for _i in 0..10 {
        // scene.add_solid(Solid::new_wall([-1.-(i as Float), 0.,0.],[1., 0., 0.]));
        // scene.add_solid(Solid::new_sphere([-3.-(i as Float), 0.,0.], 1.))
        // scene.add_solid(Solid::new_heart([-3.-(i as Float), 0.,0.]))
    }
    // scene.add_solid(Solid::new_wall([-1.-(2 as Float), 0.,0.],[1., 0., 0.]));
    // scene.add_solid(Solid::new_wall([-1.-(2 as Float), 0.,0.],[1., 0., 0.]));

    scene.run(time_step, duration, rays_per_sub_pixel);

    scene.camera.photo("time_test");
}

#[allow(dead_code)]
fn scene_testing() {
    let focal_point = [-0.9, 0., 0.];
    let screen_center = [0., 0., 0.];

    let pixel_size = 5e-3;
    let im_w = 256;                 // Number of rays created is im_w*im_h*pix_w^2*rays_per_sub_pixel rays
    let im_h =  128;
    let pix_w = 2;
    let rays_per_sub_pixel = 1;

    
    let duration = 40.;
    let time_step = 0.5;

    let camera = Camera::new(focal_point, screen_center, im_w, im_h, pixel_size, pix_w);
    let mut scene = Scene::new(camera);
    
    scene.add_solid(Solid::new_wall([0., 0., -1.], [0., 0., 1.], [-5.,5.,-5.,5.,-2.,0.]));
    
    scene.add_solid(Solid::new_heart([3., 0., -0.1]));
    scene.add_light(Solid::new_sphere([2., -1., 6.], 4.));
    // scene.add_light(Solid::new_sphere([0., -6., 6.], 4.));
    
    scene.run(time_step, duration, rays_per_sub_pixel);

    scene.camera.photo("scene_test");
}

#[allow(dead_code)]
fn vertical() {
    let focal_point = [-4., 0., 0.];
    let screen_center = [-3., 0., 0.];

    let pixel_size = 2e-3;
    let im_w = 720;                 // Number of rays created is im_w*im_h*pix_w^2*rays_per_sub_pixel rays
    let im_h =  1280;
    let pix_w = 5;
    let rays_per_sub_pixel = 1;

    
    let duration = 50.;
    let time_step = 0.05;

    let camera = Camera::new(focal_point, screen_center, im_w, im_h, pixel_size, pix_w);
    let mut scene = Scene::new(camera);
    
    scene.add_light(Solid::new_wall([0., 0., 4.], [0., 0., -1.], [-5.,5.,-5.,5.,0.,2.]));
    scene.add_solid(Solid::new_wall([0., 0., -4.], [0., 0., 1.], [-5.,5.,-5.,5.,-2.,0.]));
    scene.add_solid(Solid::new_wall([4., 0., 0.], [-1., 0., 0.], [0.,2.,-5.,5.,-5.,5.]));
    scene.add_solid(Solid::new_wall([-4., 0., 0.], [1., 0., 0.], [-2.,0.,-5.,5.,-5.,5.]));
    scene.add_solid(Solid::new_wall([0., -4., 0.], [0., 1., 0.], [-5.,5.,-2.,0.,-5.,5.]));
    scene.add_solid(Solid::new_wall([0., 4., 0.], [0., -1., 0.], [-5.,5.,0.,2.,-5.,5.]));
    
    scene.add_mirror(Solid::new_sphere([0., 0., 0.], 1.));
    
    scene.run(time_step, duration, rays_per_sub_pixel);

    scene.camera.photo("vertical");
}

fn main() {
    
    

    // poly_testing();
    // png_testing();
    // screen_testing();
    // diffuse_testing();
    // transport_testing();
    // time_testing();
    // scene_testing();
    vertical();

    
}

