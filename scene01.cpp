#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>

using namespace std;

struct Vec {
    double x, y, z;

    static const Vec Zero;

    explicit Vec(double x = 0, double y = 0, double z = 0) : x(x), y(x), z(x) {}

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

const Vec Vec::Zero(0, 0, 0);

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

	explicit Face(const Vec &a = Vec::Zero, const Vec &b = Vec::Zero, const Vec &c = Vec::Zero, const Vec &d = Vec::Zero, const Vec &n = Vec::Zero,
                  const Vec &color = Vec::Zero, const Vec &emission = Vec::Zero, Reflection_Type reflection_type = DIFF) :
            a(a), b(b), c(c), d(d), n(n),
            color(color), emission(emission), reflection_type(reflection_type) {}

    double intersect(const Ray &ray); 
    
};

double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

static Vec vertex[1000];
static Face faces[11];

void ReadObj() {
    Vec vertex_normal[1000];
    int vertex_index = 0 , faces_index = 0, normal_index = 0;
    char line[100];
    char sep[] = " \n";
    char *point;
    double a, b, c, d;

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
            a = atof(strtok(nullptr, sep)) - 1;
            b = atof(strtok(nullptr, sep)) - 1;
            c = atof(strtok(nullptr, sep)) - 1;
            d = atof(strtok(nullptr, sep)) - 1;
            cout << a << endl;
            faces[faces_index] = Face(vertex[(int)a], vertex[(int)b], vertex[(int)c], vertex[(int)d], Vec(5, 5, 5), Vec(1, 1, 1), Vec(.5, .5, .5), DIFF);
            cout << faces[faces_index].a.y << endl;
            faces_index++;
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
