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
typedef struct atom {
    char symbol[3];
    int A;
    double Z;
}atom;
typedef struct {
    Vector2 position;
    bool active;
} Emitter;
typedef struct {
    float angle;
    int orbital;
} HitResult;
extern const int width;
extern const int height;
// --- HELPER FUNCTIONS FOR NOSTDLIB ---
double wrap_angle(double x);
int q_atoi(const char *s);
float q_atof(const char *s);
double q_sin(double x);
double q_cos(double x);
float q_fmodf(float x, float y);
void q_strcpy(char* dest, const char* src);
float q_sqrtf(float x);
float q_fabsf(float x);
double power_by_squaring(double x, long long y);
double my_log2(double x);
double my_exp2(double x);
double q_pow(double x, double y);
float q_powf(float x, float y);
// --- ATOMIC PHYSICS FUNCTIONS ---
/**
 * Calculates the Binding Energy of a nucleus in amu and MeV 
 * using the Semi-Empirical Mass Formula (SEMF).
 */
double binding_energy(atom a);
/**
 * E=mc^2 (Energy in Joules)
 */
//for later use
double energy(double mass);
/**
 * Converts Energy (eV) to a Raylib Color.
 * Uses exact spectral mapping for the visible range.
 */
static inline Color GetColorFromEV(float deltaE);
float to_nm(float eV);
float to_eV(float nm);
float frequency(float eV);
float orbital_Energy(int n, int Z);
float orbital_radius(int n, int Z, float scale);
bool is_stable(atom a);
atom define(int num,int isotope);
void wave(float nm, float scale, float startX, float startY, float length, float time, float angledeg);
HitResult CheckElectronHit(atom particle, float current_time_angle, float startX, float startY, float zoom);
float radii(int orbital, int atomic_number, float scale);
#endif