const int width = 800;
const int height = 800;
typedef struct atom {
    char symbol[3];
    int A;
    double Z;
}atom;
#include "atom.h"
#include "raylib.h"
static atom define(int num,int isotope){
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
    double energy_amu = binding_energy(num, A) / c_square;
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
static float radii(int orbital, int atomic_number, float scale) {
    return 53.0f * orbital_radius(orbital, atomic_number, scale);
}
void electron(float angle, float radius, float size) {
    float ex = (width / 2.0f) + (float)q_cos((double)angle) * radius;
    float ey = (height / 2.0f) + (float)q_sin((double)angle) * radius;
    
    DrawCircle((int)ex, (int)ey, size, BLUE);
    DrawText("-", (int)ex - 2, (int)ey - 5, 10, WHITE);
}
void orbital(float radius) {
    DrawCircleLines(width / 2, height / 2, (int)radius, WHITE);
}

void proton(float radius, float px, float py) {
    DrawCircle((int)px, (int)py, radius, RED);
    DrawText("+", (int)px-3, (int)py-5, radius/10, WHITE);
}

void neutron(float radius, float px, float py) {
    DrawCircle((int)px, (int)py, radius,GetColor(0xFFAD00EE));
   DrawText("n", (int)px - 3, (int)py - 5, radius/10, WHITE);
}

void nucleus(atom atom, float radius) {
     Image icon = LoadImage("atom.png"); 
     SetWindowIcon(icon);
     UnloadImage(icon);
     int total_particles = (int)atom.A;
    if (total_particles <= 0) return;
    int sqrt_particles = (int)q_sqrtf((float)total_particles);
    if (sqrt_particles * sqrt_particles < total_particles) sqrt_particles++;
    float max_nucleus_width = radius; 
    float size = max_nucleus_width/sqrt_particles ;
    if (size > 8.0f) size = 8.0f; 
    float space = size*1.8f; 
    int k = 0; 
    int l = 0; 
    float fs = ((sqrt_particles - 1) * space) / 2.0f;

    for (int i = 0; i < sqrt_particles; i++) {
        for (int j = 0; j < sqrt_particles; j++) {
            if (k + l >= total_particles) break;

            float nx = (width / 2.0f) - fs + (j * space);
            float ny = (height / 2.0f) - fs + (i * space);

            if (k < (int)atom.Z) {
                proton(size, nx, ny);
                k++;
            } else {
                neutron(size, nx, ny);
                l++;
            }
        }
    }
}


void Atom(atom particle, float start_angle, float size)
{
    int re = (int)particle.Z;
    int shell = 1;
    float angle = start_angle;  

    while (re > 0)
    {
        int max_in_shell = 2 * (int)q_pow((double)shell, 2.0); 
        int shelle = (re > max_in_shell) ? max_in_shell : re;
        float radius = radii(shell+3, (int)particle.Z, size);
        if(shell == 1){
          nucleus(particle,radius);
        }
        orbital(radius*1.2f);
        if (shelle > 0) {
            float angle_step = TWO_PI / (float)shelle;    
            for (int i = 0; i < shelle; i++) {
                float ea = start_angle + (i * angle_step);
                electron(ea, radius*1.2f, 4);
            }
        }
        re -= shelle;
        shell++;
    }
    DrawText(TextFormat("\nSymbol:%s\nAtomic Number:%d\nMass:%d", particle.symbol,(int)particle.Z, particle.A),(width/2)-100, height - 100, 20, WHITE);
}
int main(int argc, char** argv) { 
    int element_number = 1;
        if (argc >= 2) {
        element_number = q_atoi(argv[1]);
    }
    atom he=define(element_number,-1);
    int n = 1;
    float speed = 1.0f;
    float angle = 0.0f;
    const float k = 1.0f;
    InitWindow(width, height, "atom");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        angle += speed * dt;
        BeginDrawing();
        ClearBackground(BLACK);
        Atom(he, angle,width/q_atof(argv[2]));
        DrawText(TextFormat("~ %f picometer zoom",q_atof(argv[2])/10), 10, 10, 20, WHITE);

        EndDrawing();
    }

    CloseWindow();
    return 0; 
} 


