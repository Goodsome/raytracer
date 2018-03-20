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

struct Ray {
    Vec origin, direction;

    Ray(const Vec &origin, const Vec &direction) : origin(origin), direction(direction) {}
};

enum Reflection_Type {
    DIFF,
    SPEC,
    REFR
};
struct Face {
	Vec a, b, c, d, n, color, emission;
    Reflection_Type reflection_type;

	explicit Face(const Vec &a, const Vec &b, const Vec &c, const Vec &d, const Vec &n,
                  const Vec &color, const Vec &emission, Reflection_Type reflection_type) :
            a(a), b(b), c(c), d(d), n(n),
            color(color), emission(emission), reflection_type(reflection_type) {}

    double intersect(const Ray &ray) const {

    }
};

double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}
static Vec vertex[1000];
static Face faces[1000];

void ReadObj() {
    Vec vertex_normal[1000];
    int vertex_index = 0 , faces_index = 0, normal_index = 0;
    char line[100];
    char sep[] = " \n";
    char *point;

    FILE *read_obj = fopen("scene01.obj", "r");
    if(!read_obj)
        cout << "FILE NOT OPEN!" << endl;

    for (int i = 0; i < 1000; ++i) {
        if(!fgets(line, 100, read_obj)) {
            cout << i << endl;
            break;
        }
        if (line[0] == 'v') {
            if (line[1] == ' '){
                strtok(line, sep);
                vertex[vertex_index].x = atof(strtok(nullptr, sep));
                vertex[vertex_index].y = atof(strtok(nullptr, sep));
                vertex[vertex_index].z = atof(strtok(nullptr, sep));
                vertex_index++;
            }
            else if (line[1] == 'n') {
                strtok(line, sep);
                vertex_normal[normal_index].x = atof(strtok(nullptr, sep));
                vertex_normal[normal_index].y = atof(strtok(nullptr, sep));
                vertex_normal[normal_index].z = atof(strtok(nullptr, sep));
            }
        }

        else if (line[0] == 'f') {
            strtok(line, sep);
            faces[faces_index].a = static_cast<int>(atof(strtok(nullptr, sep)) - 1);
            faces[faces_index].b = static_cast<int>(atof(strtok(nullptr, sep)) - 1);
            faces[faces_index].c = static_cast<int>(atof(strtok(nullptr, sep)) - 1);
            faces[faces_index].d = static_cast<int>(atof(strtok(nullptr, sep)) - 1);
            faces_index++;
        }
    }
    for (int i = 0; i < 11; i ++){
        cout << faces[i].a << " " << faces[i].b << " " << faces[i].c << " " << faces[i].d << endl;
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
