#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>

using namespace std;

struct Vec {
    double x, y, z;

    static const Vec Zero;

    explicit Vec(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

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

    Vec normalization() const {
        return *this * (1.0 / sqrt(x * x + y * y + z * z));
    }

    double dot(const Vec &b) const {
        return x * b.x + y * b.y + z * b.z;
    }

    Vec cross(const Vec &b) const {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }

    bool operator==(const Vec &vec) const {
        return (x == vec.x)&&(y == vec.y)&&(z == vec.z);
    }
};

const Vec Vec::Zero(0, 0, 0);

struct Ray {
    Vec origin, direction;

    Ray(const Vec &origin, const Vec &direction) : origin(origin), direction(direction) {}
};

enum Material : char{
    blinn1SG,
    blinn2SG,
    initialShadingGroup,
    lambert2SG,
    lambert3SG
};

bool InFace(const Vec &a, const Vec &b, const Vec &c, const Vec &p) {
    return  (b - a).cross(c - a).normalization() == (b - a).cross(p - a).normalization();
}

struct Face {
	Vec a, b, c, d, n;
    char *material;

	explicit Face(const Vec &a = Vec::Zero, const Vec &b = Vec::Zero,
                  const Vec &c = Vec::Zero, const Vec &d = Vec::Zero,
                  const Vec &n = Vec::Zero, char *material = (char*) initialShadingGroup) :
    a(a), b(b), c(c), d(d), n(n), material(material) {}

    double intersect(const Ray &ray) const {
        double t = (a.dot(n) - ray.origin.dot(n)) / (ray.direction.dot(n));
        Vec point = ray.origin + ray.direction * t;
        const double eps = 1e-4;

        if (t < eps)
            return 1e20;
        else
            return InFace(a, b, c, point) && InFace(b, c, a, point) && InFace(c, a, b, point)
                   && InFace(a, c, d, point) && InFace(c, d, a, point) && InFace(d, a, c, point) ? t : 1e20;
    }
    
};

static Face faces[11];

bool intersect(const Ray &ray, double &t, int &id) {
    double n = sizeof(faces) / sizeof(Face), d, inf = t = 1e20;
    for(auto i = int(n); i--;) {
        if((d = faces[i].intersect(ray)) && d < t){
            t = d;
            id = i;
        }
    }
    return t < inf;
}
double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

void ReadObj() {
    Vec vertex[100];
    Vec vertex_normal[100];
    int vertex_index = 0 , faces_index = 0, normal_index = 0;
    char line[100];
    char sep[] = " \n";
    char *material, mtl[20];
    double a, b, c, d;

    FILE *read_obj = fopen("scene01.obj", "r");
    if(!read_obj)
        cout << "FILE NOT OPEN!" << endl;


    while(fgets(line, 50, read_obj)) {
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
                normal_index++;
            }
        }
        if (line[0] == 'u') {
            strtok(line, sep);
            material = strtok(nullptr, sep);
            strcpy(mtl, material);
        }
        if (line[0] == 'f') {
            strtok(line, sep);
            a = atof(strtok(nullptr, sep)) - 1;
            b = atof(strtok(nullptr, sep)) - 1;
            c = atof(strtok(nullptr, sep)) - 1;
            d = atof(strtok(nullptr, sep)) - 1;
            faces[faces_index] = Face(vertex[(int)a], vertex[(int)b],
                                      vertex[(int)c], vertex[(int)d],
                                      vertex_normal[faces_index * 4], mtl);
            faces_index++;
        }
    }
    fclose(read_obj);
}

int main(int argc, char *argv[]) {

    const int width = 800;
    const int height = 600;
    const int samples = argc == 2 ? atoi(argv[1]) / 4 : 25;
    const Ray camera(Vec(0, 6, 10), Vec(0, -1, -9).normalization());
    const Vec cx(width * .5 / height);
    const Vec cy = (cx.cross(camera.direction)).normalization() * .5;
    auto *color = new Vec[width * height];

    ReadObj();

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            const int i = (height - y - 1) *  width + x;
            unsigned short Xi[3] = {0, 0, y * y * y};

            for (int sy = 0; sy < 2; sy++){
                Vec r = Vec::Zero;
                for (int sx = 0; sx < 2; sx++) {
                    for (int s = 0; s < samples; s++) {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / height - .5) + camera.direction;

                        //r = r + radiance(Ray(camera.origin + d * 140, d.normalization()), 0, Xi) * (1 / samples):;
                    }
                }
                color[i] = color[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
            }

        }
    }

    FILE *draw_pixels = fopen("scene01.ppm", "w");
    fprintf(draw_pixels, "P3\n%d %d \n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(draw_pixels, "%d %d %d\n", toInt(color[i].x), toInt(color[i].y), toInt(color[i].z));
    fclose(draw_pixels);
}
