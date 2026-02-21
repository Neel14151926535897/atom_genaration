/* copyright 2022-2026 */
#ifndef PARTICLE_H
#define PARTICLE_H
#include <raylib.h>

// --- Constants ---
#define c_square 931.494061             // MeV/amu
#define RYDBERG_ENERGY 13.605693122994f // eV
#define BOHR_RADIUS_ANGSTROMS 0.529177f // Angstroms
#define H_PLANCK_EV 4.135667696e-15f    // eV·s
#define HC_EV_NM 1239.84193f    // eV·nm (More precise than 1240)
#define TWO_PI 6.28318530717958647692
double wrap_angle(double x) {
    while (x > PI)  x -= TWO_PI;
    while (x < -PI) x += TWO_PI;
    return x;
}
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
/**
 * Calculates the Binding Energy of a nucleus in amu and MeV 
 * using the Semi-Empirical Mass Formula (SEMF).
 */
float q_sqrtf(float x) {
    float xhalf = 0.5f * x;
    int i = *(int*)&x;            // Get bits of float as integer
    i = 0x5f3759df - (i >> 1);    // The "Magic Number" bit-shift
    x = *(float*)&i;              // Convert bits back to float
    x = x * (1.5f - xhalf * x * x); // One iteration of Newton-Raphson
    return 1.0f / x;
}
//for reading without any librarys for our pe file

//for reading without any librarys for our pe file
float q_fabsf(float x) {
    return (x < 0.0f) ? -x : x;
}
// --- HELPER FUNCTIONS FOR NOSTDLIB ---

// Fast integer power: x^y using binary exponentiation
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

// Simple Log2 approximation (Needed for non-integer powers)
double my_log2(double x) {
    // This is a very basic approximation. 
    // For true precision without stdlib, you'd need a Remez polynomial.
    union { double d; long long i; } u = {x};
    double ex = (double)((u.i >> 52) & 0x7FF) - 1023;
    return ex; // Returns the exponent part (very rough log2)
}

// Simple Exp2 approximation
double my_exp2(double x) {
    // Basic approximation: 2^x
    union { long long i; double d; } u;
    u.i = (long long)((x + 1023) * (1LL << 52));
    return u.d;
}

// Your Main Pow Function
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
static inline double binding_energy(int Z, int A) {
    if (A <= 1) return 0;
    // Coefficients in amu and MeV
    double av = 15.8;    // Volume
    double as = 18.3;    // Surface
    double ac = 0.714;   // Coulomb
    double aa = 23.2;    // Asymmetry
    double ap;
    if (A % 2 != 0) ap = 0;                  
    else if (Z % 2 == 0) ap = 12.0;          
    else ap = -12.0;                         
    double term1 = av * A;
    double term2 = as * q_pow(A, 2.0/3.0);
    double term3 = ac * (Z * (Z - 1)) / q_pow(A, 1.0/3.0);
    double term4 = aa * q_pow(A - 2 * Z, 2) / A;
    double term5 = ap / q_sqrtf(A);

    return term1 - term2 - term3 - term4 + term5; 
}
/**
 * E=mc^2 (Energy in Joules)
 */
static inline double energy(double mass) {
    const double c = 299792458.0; 
    return mass * c * c; 
}

/**
 * Converts Energy (eV) to a Raylib Color.
 * Uses exact spectral mapping for the visible range.
 */
static inline Color GetColorFromEV(float deltaE) {
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
static inline float to_nm(float eV) {
    return (eV <= 0.0f) ? 0.0f : HC_EV_NM / eV;
}
static inline float to_eV(float nm) {
    return (nm <= 0.0f) ? 0.0f : HC_EV_NM / nm;
}
static inline float frequency(float eV) {
    return (eV <= 0.0f) ? 0.0f : eV / H_PLANCK_EV; 
}
static inline float orbital_Energy(int n, int Z) {
    if (n <= 0) return 0.0f;
    return -(((float)(Z*Z) * RYDBERG_ENERGY) / (float)(n*n)); 
}
static inline float orbital_radius(int n, int Z, float scale) {
    if (n <= 0 || Z <= 0) return 0.0f;
    return (((float)(n * n) / (float)Z) * BOHR_RADIUS_ANGSTROMS)* scale;
}
bool is_stable(int Z, float A) {
    if (A < 1) return 0;
    if (Z == 1) return 1;
    float ratio = (A - (float)Z) / (float)Z;
    if (Z <= 20) {
        return (ratio >= 1.0 && ratio <= 1.16);
    } 
    if (Z <= 82) {
        return (ratio >= 1.1 && ratio <= 1.6);
    } 
    return 0; 
}
#endif