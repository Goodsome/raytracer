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

    double norm() const {
        return sqrt(x * x + y * y + z * z);
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

    double max() const {
        return x > y && y > z ? x : y > z ? y : z;
    }

    void show() const {
        cout << x << " " << y << " " << z << endl;
    }
};

const Vec Vec::Zero(0, 0, 0);

struct Ray {
    Vec origin, direction;

    Ray(const Vec &origin, const Vec &direction) : origin(origin), direction(direction) {}
};

struct Mtl {
    Vec Kd, Ka, Tf, Ks;
    double Ni;
    static const Mtl Zero;

    Mtl (const Vec &Kd, const Vec &Ka, const Vec &Tf, double Ni, const Vec &Ks) :
            Kd(Kd), Ka(Ka), Tf(Tf), Ni(Ni), Ks(Ks) {}
};

const Mtl Mtl::Zero(Vec::Zero, Vec::Zero, Vec::Zero, 0, Vec::Zero);

static Mtl blinn1SG {
        Mtl (Vec(0.52, 0.52, 0.52),
             Vec(0.27, 0.27, 0.27),
             Vec(0.00, 0.00, 0.00),
             1.80,
             Vec(1.00, 1.00, 1.00))
};

static Mtl blinn2SG {
        Mtl (Vec(0.80, 0.80, 0.80),
             Vec(2.00, 2.00, 2.00),
             Vec(1.00, 1.00, 1.00),
             1.00,
             Vec(0.50, 0.50, 0.50)
        )
};

static Mtl initialShadingGroup {
        Mtl (Vec(0.50, 0.50, 0.50),
             Vec(0.00, 0.00, 0.00),
             Vec(1.00, 1.00, 1.00),
             1.00,
             Vec::Zero
        )
};
static Mtl lambert2SG{
        Mtl (Vec(1.00, 0.00, 0.00),
             Vec(0.00, 0.00, 0.00),
             Vec(1.00, 1.00, 1.00),
             1.00,
             Vec::Zero
        )
};
static Mtl lambert3SG{
       Mtl (Vec(0.00, 0.01, 1.00),
            Vec(0.00, 0.00, 0.00),
            Vec(1.00, 1.00, 1.00),
            1.00,
            Vec::Zero
        )
};

Mtl mtl(const char material[20]) {
    if (strcmp(material, "blinn1SG") == 0)
        return blinn1SG;
    else if (strcmp(material, "blinn2SG") == 0)
        return blinn2SG;
    else if (strcmp(material, "initialShadingGroup") == 0) {
        return initialShadingGroup;
    }
    else if (strcmp(material, "lambert2SG") == 0)
        return lambert2SG;
    else if (strcmp(material, "lambert3SG") == 0)
        return lambert3SG;
    else
        return Mtl::Zero;
}

inline bool InFace(const Vec &a, const Vec &b, const Vec &c, const Vec &p) {
    return  ((b - a).cross(c - a).normalization() - (b - a).cross(p - a).normalization()).norm() < 1e-6;
}

struct Face {
	Vec a, b, c, d, n;
    Vec cc;
    Mtl material;
    double maxColor;

	explicit Face(const Vec &a = Vec::Zero, const Vec &b = Vec::Zero,
                  const Vec &c = Vec::Zero, const Vec &d = Vec::Zero,
                  const Vec &n = Vec::Zero, const Mtl &m = Mtl::Zero) :
    a(a), b(b), c(c), d(d), n(n), material(m) {
        maxColor = material.Kd.max();
        cc = material.Kd * (1.0 / maxColor);
    }

    double intersect(const Ray &ray) const {
        double t = (a.dot(n) - ray.origin.dot(n)) / (ray.direction.dot(n));
        Vec point = ray.origin + ray.direction * t;
        const double eps = 1e-4;

        if (t < eps)
            return 1e20;
        else
            return (InFace(a, b, c, point) && InFace(b, c, a, point) && InFace(c, a, b, point))
                   || (InFace(a, c, d, point) && InFace(c, d, a, point) && InFace(d, a, c, point)) ? t : 1e20;
    }

};

static Face faces[11];

inline Face *intersect(const Ray &ray, double &t) {
    t = 1e20;
    Face *ret = nullptr;

    for (Face *f = faces; f != faces + sizeof(faces) / sizeof(Face); ++f) {
        double d = f->intersect(ray);
        if (d < t) {
            t = d;
            ret = f;
        }
    }
    return ret;
}

static Vec radiance(const Ray &ray, int depth, unsigned short *xi) {
    double t;
    Face *obj = intersect(ray, t);

    if(!obj)
        return Vec::Zero;
    else {
        int newDepth = depth + 1;
        bool isMaxDepth = newDepth > 100;
        bool isUseRR = newDepth > 5;
        bool isRR = isUseRR && erand48(xi) < obj->maxColor;


        if (isMaxDepth || (isUseRR && !isRR)){
            return obj->material.Ka;
        }
        else {
            Vec f = isRR ? obj->cc : obj->material.Kd;
            Vec x = ray.origin + ray.direction * t;
            Vec n = obj->n.normalization();
            Vec nl = n.dot(ray.direction) < 0 ? n : n * -1;

            double r1 = 2 * M_PI * erand48(xi);
            double r2 = erand48(xi);
            double r2s = sqrt(r2);

            Vec w = nl;
            Vec wo = fabs(w.x) > .1 ? Vec(0, 1) : Vec(1);
            Vec u = (wo.cross(w)).normalization();
            Vec v = w.cross(u);

            Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 -r2)).normalization();
            return obj->material.Ka + f * radiance(Ray(x, d), newDepth, xi);
        }
    }
}
inline double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

void ReadObj() {
    Vec vertex[100];
    Vec vertex_normal[100];
    Vec a, b, c, d;
    int vertex_index = 0 , faces_index = 0, normal_index = 0, i = 0;
    char line[100];
    char sep[] = " /\n";
    char *p, material[20];

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
            p = strtok(nullptr, sep);
            strcpy(material, p);
        }
        if (line[0] == 'f') {
            strtok(line, sep);
            a.x = atof(strtok(nullptr, sep)) - 1;
            a.y = atof(strtok(nullptr, sep)) - 1;
            a.z = atof(strtok(nullptr, sep)) - 1;
            b.x = atof(strtok(nullptr, sep)) - 1;
            b.y = atof(strtok(nullptr, sep)) - 1;
            b.z = atof(strtok(nullptr, sep)) - 1;
            c.x = atof(strtok(nullptr, sep)) - 1;
            c.y = atof(strtok(nullptr, sep)) - 1;
            c.z = atof(strtok(nullptr, sep)) - 1;
            p = strtok(nullptr, sep);
            if (p == nullptr){
                faces[faces_index] = Face(vertex[(int)a.x], vertex[(int)b.x],
                                          vertex[(int)c.x], vertex[(int)b.x],
                                          vertex_normal[(int)a.z] +
                                          vertex_normal[(int)b.z] +
                                          vertex_normal[(int)c.z], mtl(material));
                continue;
            }
            d.x = atof(p) - 1;
            d.y = atof(strtok(nullptr, sep)) - 1;
            d.z = atof(strtok(nullptr, sep)) - 1;
            faces[faces_index] = Face(vertex[(int)a.x], vertex[(int)b.x],
                                      vertex[(int)c.x], vertex[(int)d.x],
                                      vertex_normal[(int)a.z] +
                                      vertex_normal[(int)b.z] +
                                      vertex_normal[(int)c.z] +
                                      vertex_normal[(int)d.z], mtl(material));
            a.show();
            vertex[(int)a.x].show();
            faces[faces_index].a.show();
            cout << endl;
            faces_index++;
        }
    }
    fclose(read_obj);
}

int main(int argc, char *argv[]) {

    clock_t start = clock();
    const int width = 256;
    const int height = 256;
    const int samples = argc == 2 ? atoi(argv[1]) / 4 : 1;
    const Ray camera(Vec(0, 6, 30), Vec(0, -0.4, -10).normalization());
    const Vec cx(width * .5 / height);
    const Vec cy = (cx.cross(camera.direction)).normalization() * .5;
    auto *color = new Vec[width * height];

    ReadObj();

    for (int y = 0; y < height; y++) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples*4, 100.*y/(height - 1));
        unsigned short xi[3] = {0, 0, y*y*y};
        for (int x = 0; x < width; x++) {
            const int i = (height - y - 1) *  width + x;

            for (int sy = 0; sy < 2; sy++){
                Vec r = Vec::Zero;
                for (int sx = 0; sx < 2; sx++) {
                    for (int s = 0; s < samples; s++) {
                        double r1 = 2 * erand48(xi), dx = r1 < 2 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(xi), dy = r2 < 2 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / height - .5) + camera.direction;

                        r = r + radiance(Ray(camera.origin + d * 15, d.normalization()), 0, xi) * (1.0 / samples);
                    }
                    color[i] = color[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) *.25;
                }
            }
        }
    }
    cout << endl << (float)(clock() - start) / CLOCKS_PER_SEC << endl;

    FILE *draw_pixels = fopen("scene01.ppm", "w");
    fprintf(draw_pixels, "P3\n%d %d \n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(draw_pixels, "%d %d %d\n", toInt(color[i].x), toInt(color[i].y), toInt(color[i].z));
    fclose(draw_pixels);
}
