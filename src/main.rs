use std::io;
use std::fmt;




use png_encode_mini;
use space_alg;

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
impl Monomial {
    fn new(p_x: Integer, p_y: Integer, p_z: Integer) -> Monomial {
        Monomial {p_x: p_x, p_y: p_y, p_z: p_z}
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

struct Surface {
    scal_field: Polynomial,
    gradient: [Polynomial; 3]
}
impl Surface {
    fn new(scal_field: Polynomial) -> Surface {
        let gradient: [Polynomial; 3] = [scal_field.derivative('x'),
                                         scal_field.derivative('y'), 
                                         scal_field.derivative('z')];
        Surface {scal_field: scal_field, gradient: gradient}
    }

    fn hit_test(&self, p: Multivector, tol: Float) -> bool {
        let x = p.blades[1];
        let y = p.blades[2];
        let z = p.blades[3];
        
        self.scal_field.eval(x,y,z).abs() < tol

    }

    fn refl(&self, p: Multivector, velocity: Multivector) -> Multivector {
        let x = p.blades[1];
        let y = p.blades[2];
        let z = p.blades[3];

        let normal = Multivector::new(
            [0., self.gradient[0].eval(x,y,z), self.gradient[1].eval(x,y,z), self.gradient[2].eval(x,y,z), 0., 0., 0., 0.]
        );

        -normal*velocity*normal
    }
}
// fn get_sc_field() -> Vec<char> {
//     let mut input_str: String = String::new();
//     io::stdin()
//         .read_line(&mut input_str)
//         .expect("F");
//     input_str.chars()
//         .collect()
// }

struct Pixel {
    position: Multivector,
    ray_velocity: Multivector,
    r: Float,
    g: Float,
    b: Float
}

struct Camera {
    focal_point: Multivector,
    pixel_array: Vec<Pixel>
}
impl Camera {
    fn new(focal_point: [Float; 3], screen_center: [Float; 3], pixel_w: u16, pixel_h: u16, 
        pixel_spacing: Float, roll_angle: Float) -> Camera {

        let mut pixel_array: Vec<Pixel> = Vec::new();
        let focal_point: Multivector = Multivector::new_grade1(focal_point);
        let screen_center: Multivector = Multivector::new_grade1(screen_center);
        
        //Construct the screen geometry


        // Populate the screen with pixels
        for i in 1..=pixel_h {
            for j in 1..=pixel_w {
                i-(pixel_w+1)/2
                
            }
        }

        Camera {focal_point: focal_point, pixel_array: pixel_array}
    }
}

struct Scene {
    surfaces: Vec<Surface>,
    lights: Vec<Surface>,
    camera: Camera
}

fn poly_testing() {
    let a = Monomial::new(0,0,0);
    let b = Monomial::new(1,0,0);
    let c = Monomial::new(1,1,1);
    let d = Monomial::new(1,0,3);
    let e = Monomial::new(1,2,0);

    let p = Polynomial::new(vec![(-1., b), (-2.3, b), (0., c), (-1., d), (1., e)]);

    // println!("{}",a.eval(2.,2.,2.));
    // println!("{}",b);
    println!("{}",p);
    println!("d/dw: {}",p.derivative('w'));
    println!("d/dx: {}",p.derivative('x'));
    println!("d/dy: {}",p.derivative('y'));
    println!("d/dz: {}",p.derivative('z'));
}

fn cliff_testing() {
    
    let a = Multivector::new([-1.,-2.,3.,-5.,7.,11.,13.,-4.]);
    let b = Multivector::new([-7.,1.,16.,-24.,93.,-1.2,2.2,31.]);
    let c = a*b;

    println!("({})*({}) = {}", a, b, c)

}

fn png_testing() {

    let mut f = std::fs::File::create("test.png").unwrap();

    // image from bottom to top 3x2
    let image_width: u32 = 3;
    let image_height: u32 = 2;
    let image: Vec<u8> = vec!(
        // R     G     B     A
        0xff, 0x00, 0x00, 0xff,
        0x00, 0xff, 0x00, 0xff,
        0x00, 0x00, 0xff, 0xff,

        0x80, 0x00, 0x00, 0xff,
        0x00, 0x80, 0x00, 0xff,
        0x00, 0x00, 0x80, 0xff,
    );

    match png_encode_mini::write_rgba_from_u8(&mut f, &image, image_width, image_height) {
        Ok(_) => println!("Written image!"),
        Err(e) => println!("Error {:?}", e),
    }
    
}

fn transport_testing() {

}

fn main() {

    // poly_testing();
    // cliff_testing();
    png_testing();
    transport_testing();

    

    
}

