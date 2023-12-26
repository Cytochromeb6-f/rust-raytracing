
use std::{fmt, env, thread, f32::consts::{PI, TAU}, time::Instant, io::Write};
use rand::random;

// use kdam::tqdm;
use png_encode_mini;
use serde::{Serialize, Deserialize};
use serde_json::Result;

use space_alg::{Multivector, E_X, E_Y};



// Types, Integer and Float should have the same number of bits
type Float = f32;       // General float used in e.g. polynomial coefficients
type Integer = i32;     // General integer used for exponents in the monomials

#[derive(PartialEq, Eq, Hash, Clone, Copy, Serialize, Deserialize)]
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





#[derive(PartialEq, Clone, Serialize, Deserialize)]
struct Polynomial {
    terms: Vec<(Float, Monomial)>
}
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

#[derive(Clone, Serialize, Deserialize)]
struct Solid {
    // A solid which has a surface given by the zero set of a polynomial. The interior is where the polynomial is negative.
    surface: Polynomial,
    gradient: [Polynomial; 3],
    origin: [Float; 3],
    b_box: [Float; 6],        // This is relative to the Solids own origin
    color: [u8; 3],
    gloss: Float,
}

#[allow(dead_code)]
impl Solid {
    fn new(surface: Polynomial, origin: [Float; 3], b_box: [Float; 6], color: [u8; 3], gloss: Float) -> Solid {
        // Gloss=0 always give diffuse reflection and gloss=1 always give specular reflection
        let gradient: [Polynomial; 3] = [surface.derivative('x'),
                                         surface.derivative('y'), 
                                         surface.derivative('z')];
        Solid {surface, gradient, origin, b_box, color, gloss}
    }

    fn new_wall(point: [Float; 3], normal: [Float; 3], b_box: [Float; 6], color: [u8; 3], gloss: Float) -> Solid {
        // A half-space with a surface plane that is coincident with a given point and has a given normal direction

        let x1 = Monomial::new(1,0,0);
        let y1 = Monomial::new(0,1,0);
        let z1 = Monomial::new(0,0,1);

        let surface = Polynomial::new(vec![(normal[0], x1), (normal[1], y1), (normal[2], z1)]);

        Solid::new(surface, point, b_box, color, gloss)
    }

    fn new_sphere(center: [Float; 3], radius: Float, color: [u8; 3], gloss: Float) -> Solid {
        // A sphere with given a center and radius

        let b_box = [-radius, radius, -radius, radius, -radius, radius];

        let x2 = Monomial::new(2,0,0);
        let y2 = Monomial::new(0,2,0);
        let z2 = Monomial::new(0,0,2);
        let c = Monomial::new(0,0,0);

        let surface = Polynomial::new(vec![(1., x2), (1., y2), (1., z2), (-radius.powi(2), c)]);

        Solid::new(surface, center, b_box, color, gloss)
    }

     fn new_cube(center: [Float; 3], side: Float, color: [u8; 3], gloss: Float) -> Solid {
        // A superellipsoid with given a center and radius

        let radius = side/2.;

        let b_box = [-radius, radius, -radius, radius, -radius, radius];

        let x2 = Monomial::new(8,0,0);
        let y2 = Monomial::new(0,8,0);
        let z2 = Monomial::new(0,0,8);
        let c = Monomial::new(0,0,0);

        let surface = Polynomial::new(vec![(1., x2), (1., y2), (1., z2), (-radius.powi(8), c)]);

        Solid::new(surface, center, b_box, color, gloss)
    }

    fn new_heart(origin: [Float;3], color: [u8; 3], gloss: Float) -> Solid {

        let b_box = [-0.7, 0.7, -1.2, 1.2, -1., 1.3];

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

        Solid::new(surface, origin, b_box, color, gloss)
    }
    
    fn is_inside(&self, pos: Multivector) -> bool {

        let x = pos.comps()[1] - self.origin[0];
        let y = pos.comps()[2] - self.origin[1];
        let z = pos.comps()[3] - self.origin[2];

        match (x,y,z) {
            (x, _, _) if !(self.b_box[0]..self.b_box[1]).contains(&x) => false,
            (_, y, _) if !(self.b_box[2]..self.b_box[3]).contains(&y) => false,
            (_, _, z) if !(self.b_box[4]..self.b_box[5]).contains(&z) => false,
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

        (-normal*vel*normal).unit_vector()
    }

    fn refl_diffuse(&self, pos: Multivector) -> Multivector {    // Indepenent of incident velocity
        // Start with the velocity in the normal direction
        let normal = self.gradient_at(pos);

        // Make a vector from a spherical distribution. It is used to get a random rotation direction
        let theta = random::<Float>()*PI;
        let phi = random::<Float>()*TAU;
        let rand_vect = Multivector::new_grade1([phi.cos()*theta.sin(), phi.sin()*theta.sin(), theta.cos()]);
        
        // Make unit-bivector to rotate the normal along
        let rot_plane = (normal^rand_vect).unit_bivector();
        
        // Rotation angle from Lambert's cosine law
        let r = random::<Float>().asin();
        
        // Get the reflected velocity by rotating the normal and scaling it to unit length to set speed to 1
        (normal.scaled(r.cos()) + (normal<<rot_plane).scaled(r.sin())).unit_vector()
    }

}


#[derive(Clone)]
struct RaySite {
    pos: Multivector,
    ray_vel: Multivector,
    r_idx: usize,                     // Indices for the r,g,b channels of the pixel that this site belongs to
    g_idx: usize,                     //
    b_idx: usize                      //
}
impl RaySite {
    fn new(pos: Multivector, ray_vel: Multivector, r_idx: usize, g_idx: usize, b_idx: usize) -> RaySite {
        RaySite {pos, ray_vel, r_idx, g_idx, b_idx}
    }
}


#[derive(Clone, Copy)]
struct Ray<'a> {
    site: &'a RaySite,
    pos: Multivector,
    vel: Multivector,
    color: [Float; 3],
    rays_per_pixel: Float,
    has_split: bool,
    is_diffuse: bool,
}
impl Ray<'_> {
    fn new(site: &RaySite, pos: Multivector, vel: Multivector, rays_per_pixel: u16) -> Ray {
        Ray {site, pos, vel, color: [1.; 3], rays_per_pixel: (rays_per_pixel as Float), has_split: false, is_diffuse: false}
    }
    fn interact(&mut self, solid: &Solid, local_speed: Float, diffuse_speed: Float) -> bool {
        // Returns true if splitting happens

        if self.has_split {

            // Color absorption
            self.color[0] *= (solid.color[0] as Float)/256.;
            self.color[1] *= (solid.color[1] as Float)/256.;
            self.color[2] *= (solid.color[2] as Float)/256.;
            
                if random::<Float>() < solid.gloss {
                    // Specular reflection
                    self.vel = solid.refl(self.pos, self.vel);   
                    self.pos += self.vel.scaled(local_speed);   
                } else {
                    // Diffuse reflection
                    self.vel = solid.refl_diffuse(self.pos);
                    self.pos += self.vel.scaled(diffuse_speed);
                    self.is_diffuse = true;
                }
        } else {
            // Split into many less bright rays at first hit
            self.color[0] /= self.rays_per_pixel;
            self.color[1] /= self.rays_per_pixel;
            self.color[2] /= self.rays_per_pixel;
            
            self.has_split = true;
            return true
        }
        false
    }

    fn final_color(self, light_source: &Solid, brightness: &Float) -> [Float; 3] {
        [
            brightness*self.color[0]*(light_source.color[0] as Float),
            brightness*self.color[1]*(light_source.color[1] as Float),
            brightness*self.color[2]*(light_source.color[2] as Float)
        ]
    }
}

#[allow(dead_code)]
#[derive(Clone)]
struct Camera {
    focal_point: Multivector,
    screen_center: Multivector,
    im_w: u32,
    im_h: u32,
    sites: Vec<RaySite>,     // Array which holds the ray generation sites
    image: Vec<Float>,        // The actual rgba-array
}
impl Camera {
    fn new(focal_point: Multivector, screen_center: Multivector, im_w: u32, im_h: u32, pix_size: Float, site_s: u16) -> Camera {

        let now = Instant::now();          // Start of timed section

        
        // Construct the screen geometry
        let camera_direction: Multivector = screen_center - focal_point;       // Vector in the direction of the camera
        let e_z = Multivector::new_grade1([0.,0.,1.]);       // Vector pointing straigt up 
        let vertical_plane = e_z^camera_direction;                // Bivector plane
        let e_w = vertical_plane.complement().unit_vector();      // Unit vector in the screens width-direction
        let horizontal_plane = e_w^camera_direction;              // Bivector plane
        let e_h = -horizontal_plane.complement().unit_vector();   // Unit vector in the screens height-direction
        
        
        let site_size = pix_size/(site_s as Float);

        // Populate the screen with ray sites
        let mut sites: Vec<RaySite> = Vec::new();
        for i in 0..im_h {
            let y_pix = pix_size*((i as Float) - 0.5*(im_h-1) as Float);
            
            for j in 0..im_w {
                let x_pix = pix_size*((j as Float) - 0.5*(im_w-1) as Float);
                
                let r_index = 4*(im_w*i + j) as usize;   // Index the red-channel of the pixel that this site belongs t
                

                // Make a square grid of sites for each pixel
                for k in 0..site_s {
                    let site_y = y_pix + site_size*((k as Float) - 0.5*(site_s-1) as Float);
                    
                    for l in 0..site_s {     // Pixels are squares so pix_h==pix_w
                        let site_x = x_pix + site_size*((l as Float) - 0.5*(site_s-1) as Float);
                        
                        let pos = screen_center + e_w.scaled(site_x) + e_h.scaled(site_y);
                        let ray_vel = (pos - focal_point).unit_vector();
                        
                        sites.push(RaySite::new(pos, ray_vel, r_index, r_index+1, r_index+2));
                    }
                }
                
            }
        }
        
        // Make rgba screen vector (black non-transparent starting image)
        let image: Vec<Float> = [0., 0., 0., 255.].repeat((im_w*im_h) as usize);
        
        let elapsed_time = now.elapsed(); // End of timed section
        println!(
            "In {} seconds: camera with {} pixels and {} ray sites created.",
             elapsed_time.as_secs_f32(), im_w*im_h, im_w*im_h*(site_s as u32).pow(2)
        );
        
        Camera {focal_point, screen_center, im_w, im_h, sites, image}
    }

    fn photo(self, file_name: &str) -> () {
        let now = Instant::now();          // Start of timed section
        
        let mut f = std::fs::File::create(format!("{file_name}.png")).unwrap();

        let mut image_u8: Vec<u8> = Vec::new();
        for value in self.image {
            image_u8.push(value as u8)
        }
        let result = png_encode_mini::write_rgba_from_u8(&mut f, &image_u8, self.im_w, self.im_h);
        
        let elapsed_time = now.elapsed(); // End of timed section

        match result {
            Ok(_) => println!("In {} seconds: image written to \"{file_name}.png\".", elapsed_time.as_secs_f32()),
            Err(e) => println!("Error {:?} after {} seconds.", e, elapsed_time.as_secs())
        }

    }
}

#[derive(Serialize, Deserialize, Clone)]
struct Scene {
    focal_point: Multivector,
    screen_center: Multivector,
    im_w: u32,
    im_h: u32,
    pix_size: Float,
    n_ray_clones: u16,                      // Number of clones created when each original ray first hits something
    site_s: u16,                            // Square root of the number of ray-sites per pixel
    lights: Vec<(Solid, Float)>,
    solids: Vec<Solid>,
    objects: Vec<([Float; 6], Vec<Solid>)>,
    speed_zones: Vec<(Float, [Float; 6])>,  // The earliest will be used if two zones overlap
}
impl Scene {
    fn new(
        focal_point: [Float; 3], screen_center: [Float; 3], im_w: u32, im_h: u32,
        pix_size: Float, n_ray_clones: u16, site_s: u16
    ) -> Scene {

        let focal_point: Multivector = Multivector::new_grade1(focal_point);
        let screen_center: Multivector = Multivector::new_grade1(screen_center);

        Scene {
            focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, site_s,
            lights: Vec::new(), solids: Vec::new(), objects: Vec::new(), speed_zones: Vec::new()
        }
    }

    #[allow(dead_code)]
    fn new_from_json(path: &str) -> Scene {
        
        let data = std::fs::read_to_string(path).unwrap();

        serde_json::from_str(&data).unwrap()
    }
    #[allow(dead_code)]
    fn save_to_json(&self, file_name: &str) -> Result<()> {
        let file = std::fs::File::create(format!("{file_name}.json")).unwrap();
        let mut writer = std::io::BufWriter::new(file);
        
        serde_json::to_writer_pretty(&mut writer, self)?;
        writer.flush().unwrap(); 

        Ok(())
    }

    fn add_speed_zone(&mut self, speed: Float, s_box: [Float; 6]) {
        self.speed_zones.push((speed, s_box));

        // Let u be the speed at the surface of a solid which gives the desired precision.
        // A zone with speed v may then be safely placed at a distance s=v from the surface.

        // The optimal value of v depends on the width of the zone L, and is given by 
        // the geometric mean:          v = sqrt(u*L)
        // If the zone is between two surfaces at distance L, both with safe speeds u, then
        // the optimal v is: v = sqrt((u*L)/2)

        // Additional speed zones may be nested with diminishing returns
        
    }

    fn add_light(&mut self, solid: Solid, brightness: Float) {
        self.lights.push((solid, brightness));
    }
    fn add_solid(&mut self, solid: Solid) {
        self.solids.push(solid);
    }
    fn add_object(&mut self, object: ([Float; 6], Vec<Solid>)) {
        self.objects.push(object);
    }

    fn speed_at(&self, pos: Multivector, fallback: Float) -> Float {

        for (speed, s_box) in &self.speed_zones {
            match pos.comps()[1..=3] {
                [x, _, _] if !(s_box[0]..s_box[1]).contains(&x) => continue,
                [_, y, _] if !(s_box[2]..s_box[3]).contains(&y) => continue,
                [_, _, z] if !(s_box[4]..s_box[5]).contains(&z) => continue,
                _ => return *speed
            }
        }
        fallback
    }
    

    // fn run(&mut self, n_steps: u32, precise_speed: Float, diffuse_speed: Float, file_name: &str) {

    //     // Create camera
    //     let mut camera = Camera::new(self.focal_point, self.screen_center, self.im_w, self.im_h, self.pix_size, self.site_s);

    //     let now = Instant::now();          // Start of timed section

    //     // Create rays
    //     let rays_per_pixel = self.n_ray_clones + 1;
    //     let mut rays: Vec<Ray> = Vec::new();
    //     for site in &camera.sites {
    //         // Pixel grid
    //         rays.push(Ray::new(site, site.pos, site.ray_vel, rays_per_pixel));
            
            
    //     }
    //     let ray_total = rays.len();
    //     let mut rays_to_remove = Vec::new();
    //     let mut rays_to_add = Vec::new();

    //     let elapsed_time = now.elapsed(); // End of timed section
    //     println!("In {} seconds: {} rays produced.", elapsed_time.as_secs_f32(), ray_total);

    //     let now = Instant::now();          // Start of timed section
    //     // Simulate light transport
    //     for _ in tqdm!(0..n_steps) {

    //         'ray_iterator: for i in 0..rays.len() {
                
    //             // Check what speed is appropriate at current position
    //             let local_speed = if rays[i].is_diffuse {
    //                 diffuse_speed
    //             } else {
    //                 self.speed_at(rays[i].pos, precise_speed)
    //             };


    //             // Get next position
    //             let next_pos = rays[i].pos + rays[i].vel.scaled(local_speed);


    //             // Test for light hit
    //             for (light_source, brightness) in &self.lights {
    //                 if light_source.is_inside(rays[i].pos) {

    //                     let [red, green, blue] = rays[i].final_color(light_source, brightness);

    //                     camera.image[rays[i].site.r_idx] += red;
    //                     camera.image[rays[i].site.g_idx] += green;
    //                     camera.image[rays[i].site.b_idx] += blue;
                        
    //                     // Mark for deletion
    //                     rays_to_remove.push(i);

    //                     // Skip directly to next the ray
    //                     continue 'ray_iterator
    //                 }
    //             }
                
    //             // Test for object hit
    //             for (b_box, solids) in &self.objects {
    //                 match next_pos.comps()[1..=3] {
    //                     [x, _, _] if !(b_box[0]..b_box[1]).contains(&x) => continue,
    //                     [_, y, _] if !(b_box[2]..b_box[3]).contains(&y) => continue,
    //                     [_, _, z] if !(b_box[4]..b_box[5]).contains(&z) => continue,
    //                     _ => ()
    //                 }
    //                 for solid in solids {
    //                     if solid.is_inside(next_pos) {
    //                         if rays[i].interact(solid, local_speed, diffuse_speed) {
    //                             for _ in 1..(rays[i].rays_per_pixel as usize) {
    //                                 let mut clone = rays[i].clone();
    //                                 clone.interact(solid, local_speed, diffuse_speed);
    //                                 rays_to_add.push(clone);
    //                             }
    //                         }
                                
                            
    //                         continue 'ray_iterator         // Skip directly to next the ray
    //                     }
    //                 }
    //             }
    //             // Test for free surface hit
    //             for solid in &self.solids {
    //                 if solid.is_inside(next_pos) {
    //                     if rays[i].interact(solid, local_speed, diffuse_speed) {
    //                         for _ in 1..(rays[i].rays_per_pixel as usize) {
    //                             let mut clone = rays[i].clone();
    //                             clone.interact(solid, local_speed, diffuse_speed);
    //                             rays_to_add.push(clone);
    //                         }
    //                     }
    //                     continue 'ray_iterator             // Skip directly to next the ray
    //                 }
    //             }
                
                
    //             // Move forward as planned if nothing is in the way
    //             rays[i].pos = next_pos;
    //         }

    //         // Deletes all the rays that hit a light
    //         rays_to_remove.reverse();         // Iterate in reverse so hat swap_remove doesn't disturb itself
    //         for i in &rays_to_remove[..] {
    //             // Removes in constant time: by replacing with the last element
    //             rays.swap_remove(*i);
    //         }
    //         rays_to_remove.clear();
    //         rays.append(&mut rays_to_add);
            
            
    //     }
    //     let elapsed_time = now.elapsed(); // End of timed section
    //     println!(
    //         "\nIn {} seconds: ray transport simulated, {} rays unabsorbed ({}%).", 
    //         elapsed_time.as_secs_f32(), rays.len(), (100.*(rays.len() as Float))/(ray_total as Float)
    //     );
    //     camera.photo(file_name)
    // }


    fn run_mt(&mut self, n_threads: usize, n_steps: u32, precise_speed: Float, diffuse_speed: Float, file_name: &str) {
        // Multithreaded version

        // Create camera
        let mut camera = Camera::new(self.focal_point, self.screen_center, self.im_w, self.im_h, self.pix_size, self.site_s);

        let now = Instant::now();          // Start of timed section

        
        
        // Partition the screen and spawn threads
        let mut pixel_shift = 0; 
        let mut handles = Vec::new();
        for k in 0..n_threads {
            let me = self.clone();
            let cc = camera.clone();
            
            let sites_per_pixel = me.site_s.pow(2) as usize;
            let next_shift = (cc.sites.len()/sites_per_pixel*k/n_threads..cc.sites.len()/sites_per_pixel*(k+1)/n_threads).len();

            handles.push(thread::spawn(move || {
                // Get correct slice of camera.image and camera.sites
                let mut thread_img = cc.image[cc.image.len()*k/n_threads..cc.image.len()*(k+1)/n_threads].to_vec();
                let mut thread_sites = cc.sites[cc.sites.len()*k/n_threads..cc.sites.len()*(k+1)/n_threads].to_vec();

                for pix in &mut thread_sites { // Ensures that the rays get linked to the correct pixel
                    pix.r_idx -= 4*pixel_shift;
                    pix.g_idx -= 4*pixel_shift;
                    pix.b_idx -= 4*pixel_shift;
                }

                // Create rays
                let rays_per_site = me.n_ray_clones + 1;
                let rays_per_pixel = rays_per_site*me.site_s.pow(2);
                let mut thread_rays = Vec::new();
                for site in &thread_sites {

                    thread_rays.push(Ray::new(site, site.pos, site.ray_vel, rays_per_pixel));
                }
                // Simulate light transport
                let mut i = 0;
                'ray_iterator: while i < thread_rays.len() {
                    
                    'time_iterator: for _f in 0..n_steps {
                        
                        // match thread_rays[i].pos.comps()[1..=3] {
                        //     [x, _, _] if !(-5.0..5.).contains(&x) => continue,
                        //     [_, y, _] if !(-5.0..5.).contains(&y) => continue,
                        //     [_, _, z] if !(-5.0..5.).contains(&z) => continue,
                        //     _ => ()
                        // }
                        // Check what speed is appropriate at current position
                        let local_speed = if thread_rays[i].is_diffuse {
                            diffuse_speed
                        } else {
                            me.speed_at(thread_rays[i].pos, precise_speed)
                        };
        
        
                        // Get next position
                        let next_pos = thread_rays[i].pos + thread_rays[i].vel.scaled(local_speed);
        
        
                        // Test for light hit
                        for (light_source, brightness) in &me.lights {
                            if light_source.is_inside(thread_rays[i].pos) {
        
                                let [red, green, blue] = thread_rays[i].final_color(light_source, brightness);
        
                                thread_img[thread_rays[i].site.r_idx] += red;
                                thread_img[thread_rays[i].site.g_idx] += green;
                                thread_img[thread_rays[i].site.b_idx] += blue;
                                
                                // Skip directly to next the ray
                                i += 1;
                                continue 'ray_iterator
                            }
                        }
                        
                        // Test for object hit
                        for (b_box, solids) in &me.objects {
                            match next_pos.comps()[1..=3] {
                                [x, _, _] if !(b_box[0]..b_box[1]).contains(&x) => continue,
                                [_, y, _] if !(b_box[2]..b_box[3]).contains(&y) => continue,
                                [_, _, z] if !(b_box[4]..b_box[5]).contains(&z) => continue,
                                _ => ()
                            }
                            for solid in solids {
                                if solid.is_inside(next_pos) {
                                    if thread_rays[i].interact(solid, local_speed, diffuse_speed) {
                                        for _ in 1..rays_per_site{
                                            let mut clone = thread_rays[i].clone();
                                            clone.interact(solid, local_speed, diffuse_speed);
                                            thread_rays.push(clone);
                                        }
                                    }
                                        
                                    
                                    continue 'time_iterator         // Skip to next timestep
                                }
                            }
                        }
                        // Test for free surface hit
                        for solid in &me.solids {
                            if solid.is_inside(next_pos) {
                                if thread_rays[i].interact(solid, local_speed, diffuse_speed) {
                                    for _ in 1..rays_per_site {
                                        let mut clone = thread_rays[i].clone();
                                        clone.interact(solid, local_speed, diffuse_speed);
                                        thread_rays.push(clone);
                                    }
                                }
                                continue 'time_iterator         // Skip to next timestep
                            }
                        }
                        
                        // Move forward as planned if nothing is in the way
                        thread_rays[i].pos = next_pos;
                        // if _f==n_steps-1 {
                        //     println!("TIMEOUT: {i}, {}", thread_rays[i].pos)
                        // }
                    }
                    i += 1;
                }
                thread_img
            }));
            pixel_shift += next_shift;
        }

        // Recombine the image from the threads and save it to the camera
        camera.image = handles.into_iter()
            .map(|h| h.join().unwrap())
            .collect::<Vec<Vec<f32>>>()
            .concat();
    

        let elapsed_time = now.elapsed(); // End of timed section
        println!("\nIn {} seconds: ray transport simulated.", elapsed_time.as_secs_f32());
        camera.photo(file_name)
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
    let focal_point = Multivector::new_grade1([0., 0., 0.]);
    let screen_center = Multivector::new_grade1([1., 0., 0.]);

    let mut camera = Camera::new(focal_point, screen_center, 128, 128, 1e-2, 1);

    let mut blue = 0.;
    for pix in &camera.sites {
        camera.image[pix.b_idx] = blue;
        blue += 0.01;
    }
    camera.photo("screen_test")

}

#[allow(dead_code)]
fn diffuse_testing() {
    let rep = 40_000_000 ;

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
        let _out = (normal[i].scaled(r.cos()) + (normal[i]<<rot_plane).scaled(r.sin())).unit_vector();
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

        // Then rotate it by a random angle in the tangent plane of the surface.
        let horizontal_plane = unit_normal.complement();
        let phi = random::<Float>()*TAU;    // Horizontal rotation angle from unifrom distribution
        
        let c = Multivector::new_grade0(phi.cos());
        let s = horizontal_plane.scaled(phi.sin());
        let _out =  (c-s)*vel1*(c+s);
    }
    let elapsed_time = now.elapsed(); // End of timed section
    println!("Method 3 took {} seconds." , elapsed_time.as_secs_f32());
    

    let now = Instant::now();          // Start of timed section
    for i in 0..rep {
        // Currently broken because vel is not perpendicular to horizontal_plane. Needs the general rotation expression
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

    let pix_size = 5e-3;
    let im_w = 256;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
    let im_h =  128;
    let n_ray_clones = 8;
    

    
    let n_steps = 100;
    let precise_speed = 0.2;

    let mut scene = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, 1);
    
    // scene.add_light(Solid::new_sphere([3., 1.2, 0.1], 1.0));
    // scene.add_solid(Solid::new_wall([-1., 0., 0.], [1., 0., 0.]));


    // scene.add_solid(Solid::new_sphere([3., -1.2, 0.], 1.));
    // scene.add_solid(Solid::new_wall([0., 0., -1.], [0., 0., 1.]));
    // scene.add_solid(Solid::new_wall([0., 0., 1.5], [0., 0., -1.]));
    // scene.add_solid(Solid::new_wall([6., 0., 0.], [-1., 0., 0.]));
    
    
    scene.run_mt(8, n_steps, precise_speed, precise_speed, "transport_test");

}

#[allow(dead_code)]
fn time_testing() {
    let focal_point = [-0.9, 0., 0.];
    let screen_center = [0., 0., 0.];

    let pix_size = 5e-3;
    let im_w = 256;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
    let im_h =  128;
    let n_ray_clones = 3;

    
    let n_steps = 20;
    let precise_speed = 0.5;

   
    let mut scene = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, 1);
    
    // scene.add_light(Solid::new_sphere([3., 0.,0.], 1.));

    for _i in 0..10 {
        // scene.add_solid(Solid::new_wall([-1.-(i as Float), 0.,0.],[1., 0., 0.]));
        // scene.add_solid(Solid::new_sphere([-3.-(i as Float), 0.,0.], 1.))
        // scene.add_solid(Solid::new_heart([-3.-(i as Float), 0.,0.]))
    }
    // scene.add_solid(Solid::new_wall([-1.-(2 as Float), 0.,0.],[1., 0., 0.]));
    // scene.add_solid(Solid::new_wall([-1.-(2 as Float), 0.,0.],[1., 0., 0.]));

    scene.run_mt(8, n_steps, precise_speed, precise_speed, "time_test");

}

#[allow(dead_code)]
fn scene_testing() {
    let focal_point = [-0.9, 0., 0.];
    let screen_center = [0., 0., 0.];

    let pix_size = 5e-3;
    let im_w = 256;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
    let im_h =  128;
    let n_ray_clones = 3;

    
    let n_steps = 400;
    let precise_speed = 0.1;

    
    let mut scene = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, 1);
    
    // scene.add_solid(Solid::new_wall([0., 0., -1.], [0., 0., 1.], [-5.,5.,-5.,5.,-2.,0.]));
    
    // scene.add_solid(Solid::new_heart([3., 0., -0.1]));
    // scene.add_light(Solid::new_sphere([2., -1., 6.], 4.));
    // scene.add_light(Solid::new_sphere([0., -6., 6.], 4.));
    
    scene.run_mt(8, n_steps, precise_speed, precise_speed, "scene_test");

}

#[allow(dead_code)]
fn color_testing() {
    let phi = PI*0.1;
    let theta = PI*0.32;
    let r1 = 4.3;
    let r2 = 3.2;
    let focal_point = [r1*phi.cos()*theta.sin(), r1*phi.sin()*theta.sin(), r1*theta.cos()];
    let screen_center = [r2*phi.cos()*theta.sin(), r2*phi.sin()*theta.sin(), r2*theta.cos()];

    let scale = 1;
    let pix_size = 2e-3*(scale as Float);
    let im_w = 1280/scale;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
    let im_h =  1280/scale;
    let n_ray_clones = 25;

    
    let n_steps = 2000;
    let precise_speed: Float = 0.01;

    let mut scene = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, 1);
    
    
    let free_speed: Float = ((4.-1.)*precise_speed/2.).sqrt();
    
    let red = [0xff, 0x80, 0x80];
    let green = [0x80, 0xff, 0x80];
    let blue = [0x40, 0x40, 0x90];
    let orange = [0xff, 0xa5, 0x00];


    scene.add_solid(Solid::new_heart([0., 0., 0.], red, 0.));
    scene.add_speed_zone(precise_speed, [
        -0.7-free_speed,0.7+free_speed, -1.2-free_speed,1.2+free_speed, -1.-free_speed,1.3+free_speed
    ]);
    
    scene.add_light(Solid::new_wall([0., 0., 4.], [0., 0., -1.], [-5.,5., -5.,5., 0.,2.], [0xff; 3], 0.), 1.);
    
    scene.add_solid(Solid::new_wall([0., 0., -4.], [0., 0., 1.], [-5.,5., -5.,5.,-2.,0.], [0xff; 3], 0.));
    scene.add_solid(Solid::new_wall([4., 0., 0.], [-1., 0., 0.], [0.,2., -5.,5., -5.,5.], orange, 0.));
    scene.add_solid(Solid::new_wall([-4., 0., 0.], [1., 0., 0.], [-2.,0., -5.,5., -5.,5.], [0xff; 3], 0.4));
    scene.add_solid(Solid::new_wall([0., -4., 0.], [0., 1., 0.], [-5.,5., -2.,0., -5.,5.], green, 0.));
    scene.add_solid(Solid::new_wall([0., 4., 0.], [0., -1., 0.], [-5.,5., 0.,2., -5.,5.], blue, 0.));
    scene.add_speed_zone(free_speed, [
        -4.+free_speed,4.-free_speed, -4.+free_speed,4.-free_speed, -4.+free_speed,4.-free_speed
    ]);
    
    
    scene.run_mt(8, n_steps, precise_speed, precise_speed, "color_test");

}

// #[allow(dead_code)]
// fn json_testing() {

//     let focal_point = [-4., 0., 0.];
//     let screen_center = [-3., 0., 0.];
//     let pix_size = 2e-3*8.;
//     let im_w = 720/8;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
//     let im_h =  1280/8;
//     let n_ray_clones = 35;
//     let n_steps = 2000;
//     let precise_speed = 0.01;
    
//     let mut scene0 = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, 1);
//     scene0.add_light(Solid::new_sphere([0., 0., 0.], 1.5, [0xff;3], 0.), 1.);
//     let json_str = serde_json::to_string(&scene0).unwrap();

//     let scene1: Scene = serde_json::from_str(json_str.as_str()).unwrap();
//     scene1.save_to_json("test").unwrap(); 

//     let mut scene2 = Scene::new_from_json("test.json");
//     scene2.run(n_steps, precise_speed, precise_speed, "json_test");
    
// }

#[allow(dead_code)]
fn vertical() {
    let focal_point = [-4., 0., 0.];
    let screen_center = [-3., 0., 0.];

    let pix_size = 2e-3*8.;
    let im_w = 720/8;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
    let im_h =  1280/8;
    let n_ray_clones = 25;
    

    
    let n_steps = 2000;
    let precise_speed = 0.01;

    let mut scene = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, 1);
    
    let radius: Float = 1.2;
    let free_speed: Float = ((4.-radius)*precise_speed/2.).sqrt();
    
    
    scene.add_solid(Solid::new_sphere([0., 0., 0.], radius, [0xff;3], 1.));
    scene.add_speed_zone(precise_speed, [
        -radius-free_speed,radius+free_speed, -radius-free_speed,radius+free_speed, -radius-free_speed,radius+free_speed
    ]);
    
    scene.add_speed_zone(free_speed, [
        -4.+free_speed,4.-free_speed, -4.+free_speed,4.-free_speed, -4.+free_speed,4.-free_speed
    ]);

    scene.add_light(Solid::new_wall([0., 0., 4.], [0., 0., -1.], [-5.,5., -5.,5., 0.,2.], [0xff; 3], 0.), 1.);
    
    scene.add_solid(Solid::new_wall([0., 0., -4.], [0., 0., 1.], [-5.,5., -5.,5.,-2.,0.], [0xff; 3], 0.));
    scene.add_solid(Solid::new_wall([4., 0., 0.], [-1., 0., 0.], [0.,2., -5.,5., -5.,5.], [0xff; 3], 0.));
    scene.add_solid(Solid::new_wall([-4., 0., 0.], [1., 0., 0.], [-2.,0., -5.,5., -5.,5.], [0xff; 3], 0.));
    scene.add_solid(Solid::new_wall([0., -4., 0.], [0., 1., 0.], [-5.,5., -2.,0., -5.,5.], [0xff; 3], 0.));
    scene.add_solid(Solid::new_wall([0., 4., 0.], [0., -1., 0.], [-5.,5., 0.,2., -5.,5.], [0xff; 3], 0.));

    
    // scene.run(n_steps, precise_speed, precise_speed, "vertical");
    scene.run_mt(10, n_steps, precise_speed, precise_speed, "vertical");

}

#[allow(dead_code)]
fn rubikscube() {

    // Screen geometry
    let focal_point = [-3.5, 1.9, 1.6];
    let screen_center = [-2.5, 1.2, 1.2];

    let scale = 8;
    let pix_size = 2e-3*(scale as Float);
    let im_w = 1280/scale;                 // Number of rays created is im_w*im_h*(1 + n_ray_clones) rays
    let im_h =  720/scale;
    let n_ray_clones = 25;
    let site_s = 1;
    

    // Simulation parameters
    let n_steps = 2000;
    let precise_speed = 0.01;
    let diffuse_speed = precise_speed;// 1.0;
    

    // Colors
    let blue = [0x01, 0x69, 0xa8];
    let green = [0x02, 0xbf, 0xa3];
    let yellow = [0xfd, 0xd0, 0x28];
    let white = [0xe8, 0xea, 0xe9];
    let orange = [0xfd, 0x79, 0x26];
    let red = [0xd3, 0x37, 0x38];
    
    let black = [0x27, 0x29, 0x28];
    let l_grey = [0xaa; 3];
    
    
    // Scene geometry 
    let mut scene = Scene::new(focal_point, screen_center, im_w, im_h, pix_size, n_ray_clones, site_s);
    
    let side: Float = 0.6;
    
    let mut cube_object: ([f32; 6], Vec<Solid>) = ([-side*1.5,side*1.5, -side*1.5,side*1.5, -side*1.5,side*1.5], Vec::new());
    for i in -1..=1 {
        for j in -1..=1 {
            for k in -1..=1 {
                if i != 0 || j != 0 || k != 0 {             // Empty in the middle
                    let x = (i as Float)*side; 
                    let y = (j as Float)*side; 
                    let z = (k as Float)*side;

                    // Stickers
                    if i==-1 {
                        cube_object.1.push(Solid::new_wall([x-side/2., y, z], [-1.,0.,0.], 
                            [0.,0.1, -side*0.4,side*0.4, -side*0.4,side*0.4], blue, 0.2));
                    }
                    if i==1 {
                        cube_object.1.push(Solid::new_wall([x+side/2., y, z], [1.,0.,0.], 
                            [-0.1,0., -side*0.4,side*0.4, -side*0.4,side*0.4], green, 0.2));
                    }
                    if j==-1 {
                        cube_object.1.push(Solid::new_wall([x, y-side/2., z], [0.,-1.,0.], 
                            [-side*0.4,side*0.4, 0.,0.1, -side*0.4,side*0.4], white, 0.2));
                    }
                    if j==1 {
                        cube_object.1.push(Solid::new_wall([x, y+side/2., z], [0.,1.,0.], 
                            [-side*0.4,side*0.4, -0.1,0., -side*0.4,side*0.4], yellow, 0.2));
                    }
                    if k==-1 {
                        cube_object.1.push(Solid::new_wall([x, y, z-side/2.], [0.,0.,-1.], 
                            [-side*0.4,side*0.4, -side*0.4,side*0.4, 0.,0.1], orange, 0.2));
                    }
                    if k==1 {
                        cube_object.1.push(Solid::new_wall([x, y, z+side/2.], [0.,0.,1.], 
                            [-side*0.4,side*0.4, -side*0.4,side*0.4, -0.1,0.], red, 0.2));
                    }
                    // Cube
                    cube_object.1.push(Solid::new_cube([x, y, z], side, black, 0.9));
                }
            }
        }
    }
    scene.add_object(cube_object);

    let free_speed: Float = ((4.-side*1.5)*precise_speed/2.).sqrt();
    let mid_speed: Float = (free_speed*precise_speed).sqrt();

    scene.add_speed_zone(precise_speed, [
        -side*1.5-mid_speed,side*1.5+mid_speed, -side*1.5-mid_speed,side*1.5+mid_speed, -side*1.5-mid_speed,side*1.5+mid_speed
    ]);
    scene.add_speed_zone(mid_speed, [
        -side*1.5-free_speed,side*1.5+free_speed, -side*1.5-free_speed,side*1.5+free_speed, -side*1.5-free_speed,side*1.5+free_speed
    ]);
    scene.add_speed_zone(free_speed, [
        -4.+free_speed,4.-free_speed, -4.+free_speed,4.-free_speed, -4.+free_speed,4.-free_speed
    ]);
    scene.add_speed_zone(mid_speed, [
        -4.+mid_speed,4.-mid_speed, -4.+mid_speed,4.-mid_speed, -4.+mid_speed,4.-mid_speed
    ]);
    
    scene.add_light(Solid::new_sphere([-4., -3., 4.], 2., [0x8f; 3], 1.), 5.);
    scene.add_light(Solid::new_wall([0., 0., 4.], [0., 0., -1.], [-5.,5., -5.,5., 0.,2.],
        l_grey, 0.), 2.);
    
    // Walls
    scene.add_solid(Solid::new_wall([0., 0., -4.], [0., 0., 1.], [-5.,5., -5.,5.,-2.,0.],
        l_grey, 0.));
    scene.add_solid(Solid::new_wall([4., 0., 0.], [-1., 0., 0.], [0.,2., -5.,5., -5.,5.],
        l_grey, 0.));
    scene.add_solid(Solid::new_wall([-4., 0., 0.], [1., 0., 0.], [-2.,0., -5.,5., -5.,5.],
        l_grey, 0.));
    scene.add_solid(Solid::new_wall([0., -4., 0.], [0., 1., 0.], [-5.,5., -2.,0., -5.,5.],
        l_grey, 0.));
    scene.add_solid(Solid::new_wall([0., 4., 0.], [0., -1., 0.], [-5.,5., 0.,2., -5.,5.],
        l_grey, 0.));
    
    // Mirror
    scene.add_solid(Solid::new_wall([3., -3., 0.], [-1., 1., 0.], [-1.,5., -5.,1., -5.,5.],
        [0xf8; 3], 0.99));

    

    // scene.run(n_steps, precise_speed, diffuse_speed, "rubikscube");
    scene.run_mt(8, n_steps, precise_speed, diffuse_speed, "rubikscube");
    
}




fn main() {
    env::set_var("RUST_BACKTRACE", "1");

    // poly_testing();
    // png_testing();
    // screen_testing();
    // diffuse_testing();
    // transport_testing();
    // time_testing();
    // scene_testing();
    color_testing();
    // json_testing();
    
    // vertical();
    // rubikscube();
    

}