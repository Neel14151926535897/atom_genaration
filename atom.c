#include "atom.h"
int q_atoi(const char *s) {
    int sign = 1;
    int result = 0;

    // Skip leading whitespace
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r' || *s == '\v' || *s == '\f') {
        s++;
    }

    // Handle optional sign
    if (*s == '-') {
        sign = -1;
        s++;
    } else if (*s == '+') {
        s++;
    }

    // Convert digits
    while (*s >= '0' && *s <= '9') {
        result = result * 10 + (*s - '0');
        s++;
    }

    return sign * result;
}
double wrap_angle(double x) {
    while (x > PI)  x -= TWO_PI;
    while (x < -PI) x += TWO_PI;
    return x;
}
float q_atof(const char *s) {
    float sign = 1.0f;
    float result = 0.0f;
    float divisor = 1.0f;

    // 1. Skip leading whitespace
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r' || *s == '\v' || *s == '\f') {
        s++;
    }

    // 2. Handle optional sign
    if (*s == '-') {
        sign = -1.0f;
        s++;
    } else if (*s == '+') {
        s++;
    }

    // 3. Convert integer part
    while (*s >= '0' && *s <= '9') {
        result = result * 10.0f + (*s - '0');
        s++;
    }

    // 4. Convert fractional part
    if (*s == '.') {
        s++;
        while (*s >= '0' && *s <= '9') {
            result = result * 10.0f + (*s - '0');
            divisor *= 10.0f;
            s++;
        }
    }

    return (sign * result) / divisor;
}
double q_sin(double x) {
    x = wrap_angle(x);
    
    // Using a 7-term Taylor expansion for stability
    double x2 = x * x;
    double term = x;
    double sum = x;
    
    // We hard-code the divisors to avoid large factorials
    term *= -x2 / 6.0;                   // x^3 / 3!
    sum += term;
    term *= -x2 / 20.0;                  // x^5 / 5!
    sum += term;
    term *= -x2 / 42.0;                  // x^7 / 7!
    sum += term;
    term *= -x2 / 72.0;                  // x^9 / 9!
    sum += term;
    
    return sum;
}
double q_cos(double x) {
    x = wrap_angle(x);
    
    double x2 = x * x;
    double term = 1.0;
    double sum = 1.0;
    
    term *= -x2 / 2.0;                   // x^2 / 2!
    sum += term;
    term *= -x2 / 12.0;                  // x^4 / 4!
    sum += term;
    term *= -x2 / 30.0;                  // x^6 / 6!
    sum += term;
    term *= -x2 / 56.0;                  // x^8 / 8!
    sum += term;
    
    return sum;
}
float q_fmodf(float x, float y) {
    // 1. Handle "Not a Number" (NaN) and Infinity
    if (y == 0.0f || x != x || y != y) {
        return (x * y) / (x * y); // Returns NaN
    }

    // 2. Get absolute values (we'll re-apply the sign later)
    float ux = (x < 0) ? -x : x;
    float uy = (y < 0) ? -y : y;

    if (ux < uy) return x; // If x is smaller than y, the remainder is just x

    // 3. The "Loop": Subtract y from x until x < y
    // Note: Real libraries use bit-shifting on the exponent for speed!
    while (ux >= uy) {
        ux -= uy;
    }

    // 4. Return with the original sign of x
    return (x < 0) ? -ux : ux;
}
void q_strcpy(char* dest, const char* src) {
    while ((*dest++ = *src++));
}
float q_sqrtf(float x) {
    float xhalf = 0.5f * x;
    int i = *(int*)&x;            // Get bits of float as integer
    i = 0x5f3759df - (i >> 1);    // The "Magic Number" bit-shift
    x = *(float*)&i;              // Convert bits back to float
    x = x * (1.5f - xhalf * x * x); // One iteration of Newton-Raphson
    return 1.0f / x;
}
float q_fabsf(float x) {
    return (x < 0.0f) ? -x : x;
}
double power_by_squaring(double x, long long y) {
    if (y < 0) {
        x = 1.0 / x;
        y = -y;
    }
    double res = 1.0;
    while (y > 0) {
        if (y & 1) res *= x;
        x *= x;
        y >>= 1;
    }
    return res;
}
double my_log2(double x) {
    // This is a very basic approximation. 
    // For true precision without stdlib, you'd need a Remez polynomial.
    union { double d; long long i; } u = {x};
    double ex = (double)((u.i >> 52) & 0x7FF) - 1023;
    return ex; // Returns the exponent part (very rough log2)
}
double my_exp2(double x) {
    // Basic approximation: 2^x
    union { long long i; double d; } u;
    u.i = (long long)((x + 1023) * (1LL << 52));
    return u.d;
}
double q_pow(double x, double y) {
    if (y == 0) return 1.0;
    if (x == 0) return 0.0;
    
    // Check if y is an integer (e.g., 2.0, 3.0)
    long long iy = (long long)y;
    if (y == (double)iy) {
        return power_by_squaring(x, iy);
    }
    // Floating point power: x^y = 2^(y * log2(x))
    return my_exp2(y * my_log2(x));
}
float q_powf(float x, float y) {
    // 1. Special cases (The "Boring" but necessary part)
    if (y == 0) return 1.0f;
    if (x == 0) return 0.0f;
    if (x < 0 && (int)y != y) return 0.0f/0.0f; // NaN: Root of negative

    // 2. The Identity: x^y = exp(y * log(x))
    // We use double precision internally to prevent the "needle" jitter
    return (float)my_exp2(y * my_log2(x));
}
double binding_energy(atom a) {
    if (a.A <= 1) return 0;
    // Coefficients in amu and MeV
    double av = 15.8;    // Volume
    double as = 18.3;    // Surface
    double ac = 0.714;   // Coulomb
    double aa = 23.2;    // Asymmetry
    double ap;
    if (a.A % 2 != 0) ap = 0;                  
    else if (q_fmodf(a.Z, 2.0) == 0) ap = 12.0;          
    else ap = -12.0;                         
    double term1 = av * a.A;
    double term2 = as * q_pow(a.A, 2.0/3.0);
    double term3 = ac * (a.Z * (a.Z - 1)) / q_pow(a.A, 1.0/3.0);
    double term4 = aa * q_pow(a.A - 2 * a.Z, 2) / a.A;
    double term5 = ap / q_sqrtf(a.A);

    return term1 - term2 - term3 - term4 + term5; 
}
double energy(double mass) {
    const double c = 299792458.0; 
    return mass * c * c; 
}
Color GetColorFromEV(float deltaE) {
    if (deltaE <= 0.0f) return BLANK;
    float nm_val = HC_EV_NM / deltaE; 
    float r = 0.0f, g = 0.0f, b = 0.0f;
    // Spectral Mapping
    if (nm_val >= 380 && nm_val < 440) {
        r = (440 - nm_val) / (440 - 380);
        b = 1.0f;
    } else if (nm_val >= 440 && nm_val < 490) {
        g = (nm_val - 440) / (490 - 440);
        b = 1.0f;
    } else if (nm_val >= 490 && nm_val < 510) {
        g = 1.0f;
        b = (510 - nm_val) / (510 - 490);
    } else if (nm_val >= 510 && nm_val < 580) {
        r = (nm_val - 510) / (580 - 510);
        g = 1.0f;
    } else if (nm_val >= 580 && nm_val < 645) {
        r = 1.0f;
        g = (645 - nm_val) / (645 - 580);
    } else if (nm_val >= 645 && nm_val <= 780) {
        r = 1.0f;
    } else {
        if (nm_val < 380) return (Color){ 120, 0, 255, 221 }; // UV
        return (Color){ 100, 0, 0, 221 };                    // IR
    }
    // Intensity Fall-off
    float factor = 0.0f;
    if (nm_val >= 380 && nm_val < 420) factor = 0.3f + 0.7f * (nm_val - 380) / (420 - 380);
    else if (nm_val >= 420 && nm_val < 701) factor = 1.0f;
    else if (nm_val >= 701 && nm_val <= 780) factor = 0.3f + 0.7f * (780 - nm_val) / (780 - 701);
    const float gamma = 0.4545f;
    return (Color){
        (unsigned char)(q_powf(r * factor, gamma) * 255),
        (unsigned char)(q_powf(g * factor, gamma) * 255),
        (unsigned char)(q_powf(b * factor, gamma) * 255),
        221 
    };
}
float to_nm(float eV) {
    return (eV <= 0.0f) ? 0.0f : HC_EV_NM / eV;
}
float to_eV(float nm) {
    return (nm <= 0.0f) ? 0.0f : HC_EV_NM / nm;
}
float frequency(float eV) {
    return (eV <= 0.0f) ? 0.0f : eV / H_PLANCK_EV; 
}
float orbital_Energy(int n, int Z) {
    if (n <= 0) return 0.0f;
    return -(((float)(Z*Z) * RYDBERG_ENERGY) / (float)(n*n)); 
}
float orbital_radius(int n, int Z, float scale) {
    if (n <= 0 || Z <= 0) return 0.0f;
    return (((float)(n * n) / (float)Z) * BOHR_RADIUS_ANGSTROMS)* scale;
}
bool is_stable(atom a) {
    if (a.A < 1) return 0;
    if (a.Z == 1) return 1;
    float ratio = (a.A - a.Z) / a.Z;
    if (a.Z <= 20) {
        return (ratio >= 1.0 && ratio <= 1.16);
    } 
    if (a.Z <= 82) {
        return (ratio >= 1.1 && ratio <= 1.6);
    } 
    return 0; 
}
atom define(int num,int isotope){
    atom result;
    double protonsamu = num*1.007276466879;
    double electronsamu = num*0.000548579909070;
     static const unsigned short stable[118] = {1, 4, 7, 9, 11, 12, 14, 16, 19, 20,23, 24, 27, 28, 31, 32, 35, 40, 39, 40,45, 48, 51, 52, 55, 56, 59, 59, 64, 65,70, 73, 75, 79, 80, 84, 85, 88, 89, 91,93, 96, 98, 101, 103, 106, 107, 112, 115, 119,121, 128, 127, 131, 133, 137, 139, 140, 141, 144,145, 150, 152, 157, 159, 163, 165, 167, 169, 173,175, 178, 181, 184, 186, 190, 192, 195, 197, 201,204, 208, 209, 209, 210, 222, 223, 226, 227, 232,231, 238, 237, 244, 243, 247, 247, 251, 252, 257,258, 259, 266, 267, 268, 269, 270, 269, 278, 281,282, 285, 286, 289, 290, 293, 294, 294};
    double neutronsamu;
    int A;
    if(isotope <= 0){
     A= stable[num-1];
     neutronsamu = (A-num)*1.00866491588;
    }else{
     A = isotope;
     neutronsamu = ((double)isotope-num)*1.00866491588;
    }
    double energy_amu = binding_energy((atom){"He",num, A}) / c_square;
    double atomic_mass = (protonsamu + electronsamu + neutronsamu) - energy_amu;
    static const char* symbols = "H  He Li Be B  C  N  O  F  Ne Na Mg Al Si P  S  Cl Ar K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og ";
    result.A = A;
    result.Z = (double)num; 
    const char* src = &symbols[(num - 1) * 3];
    result.symbol[0] = src[0];
    result.symbol[1] = (src[1] == ' ') ? '\0' : src[1]; 
    result.symbol[2] = '\0'; 
    return result;
}
float radii(int orbital, int atomic_number, float scale) {
    float bohr_units = orbital_radius(orbital, atomic_number, scale);
    return bohr_units * scale;
}
void wave(float nm, float scale, float startX, float startY, float length, float time, float angledeg) {
    float wavelength = nm * scale;
    Color color = GetColorFromEV(to_eV(nm));
    float anglerad = angledeg * DEG2RAD;
    
    float dirX = q_cos(anglerad);
    float dirY = q_sin(anglerad);
    float perpX = -dirY; 
    float perpY = dirX;

    Vector2 points[500]; 
    int count = (int)length;
    if (count > 500) count = 500; 

    for (int i = 0; i < count; i++) {
        float baseX = startX + (i * dirX);
        float baseY = startY + (i * dirY);

        float phase = (2.0f * PI / wavelength) * i - time;
        float wiggle = q_sin(phase) * 20.0f;

        points[i].x = baseX + (wiggle * perpX);
        points[i].y = baseY + (wiggle * perpY);
    }

    DrawLineStrip(points, count, color);
}


HitResult CheckElectronHit(atom particle, float current_time_angle, float startX, float startY, float zoom) {
    HitResult result = { -1.0f, 0 };
    int re = (int)particle.Z;
    int shell = 1;

    while (re > 0) {
        int max_in_shell = 2 * (int)q_pow((double)shell, 2.0);
        int shelle = (re > max_in_shell) ? max_in_shell : re;
        float radius = radii(shell + 3, (int)particle.Z, zoom) * 1.2f;
        float angle_step = TWO_PI / (float)shelle;

        for (int i = 0; i < shelle; i++) {
            float ea = current_time_angle + (i * angle_step);
            float ex = (width / 2.0f) + (float)q_cos((double)ea) * radius;
            float ey = (height / 2.0f) + (float)q_sin((double)ea) * radius;

            if (q_fabsf(ey - startY) < 8.0f && ex > startX && q_fabsf(ex - (startX + 50)) < 15.0f) {
                result.angle = ea * RAD2DEG;
                result.orbital = shell;
                return result;
            }
        }
        re -= shelle;
        shell++;
    }
    return result;
}

