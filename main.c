#include "raylib.h"
#include "atom.h"
const int width = 800;
const int height = 800;
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
Vector2 ray = { 0, 0 };
bool rayActive = false;
float emitAngle = -1;

int main(int argc, char** argv) { 
    int element_number = (argc >= 2) ? q_atoi(argv[1]) : 1;
    int isotope = (argc >=4) ? q_atoi(argv[2]) : -1;
    float final_scale = (argc == 3) ? q_atof(argv[2]):(argc == 4) ? q_atof(argv[3]): 50.0f; 
    atom he = define(element_number, -1);
    float angle = 0.0f;
    float speed = 1.0f;
    InitWindow(width, height, "Atom");
    SetTargetFPS(60);
    while (!WindowShouldClose()) {
        angle += speed * GetFrameTime();
        BeginDrawing();
            ClearBackground(BLACK);
            Atom(he, angle, final_scale/16.6);
            DrawText(TextFormat("Scale: %.0f px per Bohr Radius", final_scale), 10, 10, 20, RAYWHITE);
            DrawLine(10, 50, 10 + (int)final_scale, 50, SKYBLUE);
            DrawText("1 Bohr Unit Distance", 10, 60, 15, SKYBLUE);

        EndDrawing();
    }
    CloseWindow();
    return 0;
}


