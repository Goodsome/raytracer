#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>

using namespace std;

struct Vec {
    double x, y, z;

    explicit Vec(double x = 0, double y =0, double z = 0) : x(x), y(y), z(z) {}

    Vec operator+(const Vec &b) const {
        return Vec(x + b.x, y + b.y, z + b.z);
    }

    Vec operator-(const Vec &b) const {
        return Vec(x - b.x, y - b.y, z - b.z);
    }

    Vec operator*(const double b) const {
        return Vec(x * b, y * b, z * b);
    }

    friend Vec operator*(const Vec &a, const Vec &b) {
        return Vec(a.x * b.x, a.y * b.y, a.z * b.z);
    }

    Vec norm() const {
        return *this * (1.0 / sqrt(x * x + y * y + z * z));
    }

    double operator%(const Vec &b) const {
        return x * b.x + y * b.y + z * b.z;
    }

    Vec operator^(const Vec &b) const {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

struct Face {
	int a, b, c, d;

	explict Face(int a = 0, int b = 0, int c = 0, int d = 0) : a(a), b(b), c(c), d(d) {}
}

double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

void ReadObj() {
    Vec vertex[1000];
    int vertex_index = 0 , faces_index = 0;
    char line[100];
    char sep[] = " \n";

    FILE *read_obj = fopen("scene01.obj", "r");
    if(!read_obj)
        cout << "FILE NOT OPEN!" << endl;

    for (int i = 0; i < 10000; ++i) {
        if(!fgets(line, 100, read_obj)) {
            cout << i << endl;
            break;
        }
        if (line[0] == 'v' & line[1] == ' ') {
            strtok(line, sep);
            vertex[vertex_index].x = atof(strtok(nullptr, sep));
            vertex[vertex_index].y = atof(strtok(nullptr, sep));
            vertex[vertex_index].z = atof(strtok(nullptr, sep));
            vertex_index++;
        }

        if (line[0] == 'f') {
            strtok(line, sep);
        }
    }
    fclose(read_obj);
}


int main(int argc, char *argv[]) {

    const int width = 800;
    const int height = 600;
    const int samples = argc == 2 ? atoi(argv[1]) / 4 : 25;
    auto *color = new Vec[width * height];

    ReadObj();

    FILE *draw_pixels = fopen("scene01.ppm", "w");
    fprintf(draw_pixels, "P3\n%d %d \n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(draw_pixels, "%d %d %d\n", toInt(color[i].x), toInt(color[i].y), toInt(color[i].z));
    fclose(draw_pixels);
}
