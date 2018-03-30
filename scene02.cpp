#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>
#include <omp.h>

using namespace std;

// Vec结构 表示空间内的向量，点的坐标以及颜色值。
struct Vec {
    double x, y, z;

    explicit Vec(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

    Vec operator+(const Vec &vec) const {
        return Vec(x + vec.x, y + vec.y, z + vec.z);
    }

    Vec operator-(const Vec &vec) const {
        return Vec(x - vec.x, y - vec.y, z - vec.z);
    }

    Vec operator*(const double b) const {
        return Vec(x * b, y * b, z * b);
    }

    Vec operator/(const Vec &vec) const {
        return Vec(x / vec.x, y / vec.y, z / vec.z);
    }

    friend Vec operator*(const Vec &a, const Vec &b) {
        return Vec(a.x * b.x, a.y * b.y, a.z * b.z);
    }

    double norm() const {// 取模
        return sqrt(x * x + y * y + z * z);
    }

    Vec normalization() const {// 向量单位化
        return *this * (1.0 / sqrt(x * x + y * y + z * z));
    }

    double dot(const Vec &vec) const {// 点乘
        return x * vec.x + y * vec.y + z * vec.z;
    }

    Vec cross(const Vec &vec) const {// 叉乘
        return Vec(y * vec.z - z * vec.y, z * vec.x - x * vec.z, x * vec.y - y * vec.x);
    }

    bool operator==(const Vec &vec) const {

        return (x == vec.x)&&(y == vec.y)&&(z == vec.z);
    }

    double max() const {
        return x > y && y > z ? x : y > z ? y : z;
    }

    void show() {
        cout << x << " " << y << " " << z << endl;
    }
};

// Ray结构 表示射线的原点和方向,方向为单位向量。
struct Ray {
    Vec origin, direction;

    Ray(const Vec &origin, const Vec &direction) : origin(origin), direction(direction.normalization()) {}
};

// Mtl结构 定义材质属性
struct Mtl {
    Vec Kd, Ka, Ks;
    double ns;

    Mtl (const Vec &Kd = Vec(), const Vec &Ka = Vec(), const Vec &Ks = Vec(), double ns = 0) :
            Kd(Kd), Ka(Ka), Ks(Ks), ns(ns) {}
};

//                  漫反射                   环境光照               镜面光           光滑度
static Mtl initialShadingGroup{// 背景墙面
        Mtl (Vec(0.50, 0.50, 0.50), Vec(0.00, 0.00, 0.00), Vec(0.00, 0.00, 0.00), 1)
};
static Mtl lambert2SG{// 光球
        Mtl (Vec(0.00, 0.00, 0.00), Vec(10.0, 10.0, 10.0), Vec(0.00, 0.00, 0.00), 1)
};

// 高光板由上到下
static Mtl mia_material_x_passes1SG{
        Mtl (Vec(0.00, 0.00, 0.00), Vec(0.00, 0.00, 0.00), Vec(1.00, 1.00, 1.00), 4000)
};
static Mtl mia_material_x_passes2SG{
        Mtl (Vec(0.00, 0.00, 0.00), Vec(0.00, 0.00, 0.00), Vec(1.00, 1.00, 1.00), 1500)
};
static Mtl mia_material_x_passes3SG{
        Mtl (Vec(0.00, 0.00, 0.00), Vec(0.00, 0.00, 0.00), Vec(1.00, 1.00, 1.00), 500)
};
static Mtl mia_material_x_passes4SG{
        Mtl (Vec(0.00, 0.00, 0.00), Vec(0.00, 0.00, 0.00), Vec(1.00, 1.00, 1.00), 100)
};

Mtl mtl(const char material[20]) {
    if (strcmp(material, "initialShadingGroup") == 0)
        return initialShadingGroup;
    else if (strcmp(material, "mia_material_x_passes1SG") == 0)
        return mia_material_x_passes1SG;
    else if (strcmp(material, "mia_material_x_passes2SG") == 0)
        return mia_material_x_passes2SG;
    else if (strcmp(material, "mia_material_x_passes3SG") == 0)
        return mia_material_x_passes3SG;
    else if (strcmp(material, "mia_material_x_passes4SG") == 0)
        return mia_material_x_passes4SG;
}

// 判断点p是否在三角形abc内
inline bool InFace(const Vec &a, const Vec &b, const Vec &c, const Vec &p) {
    return  ((b - a).cross(c - a).normalization() - (b - a).cross(p - a).normalization()).norm() < 1e-6;
}

// Face结构 abcd四个点以及材料属性， 法向量n由点abc计算得出。
struct Face {
	Vec a, b, c, d, n;
    Mtl material;

	explicit Face(const Vec &a = Vec(), const Vec &b = Vec(),
                  const Vec &c = Vec(), const Vec &d = Vec(),
                  const Mtl &m = Mtl()) : a(a), b(b), c(c), d(d), material(m) {
        n = (b - a).cross(c - a).normalization();
    }

    // 计算射线与面是否相交，若相交得出距离t。方程：(O + tD - a)n = 0
    double intersect(const Ray &ray) const {
        double t = (a.dot(n) - ray.origin.dot(n)) / (ray.direction.dot(n));
        Vec point = ray.origin + ray.direction * t;
        const double eps = 1e-4;

        if (t < eps)
            return 1e10;
        else
            return (InFace(a, b, c, point) && InFace(b, c, a, point) && InFace(c, a, b, point))
                   || (InFace(a, c, d, point) && InFace(c, d, a, point) && InFace(d, a, c, point)) ? t : 1e10;
    }

};

// Sphere结构 表示球的半径，球心坐标，材料属性
struct Sphere {
    double radius;
    Vec position;
    Mtl material;

    Sphere(double r, const Vec &p, const Mtl &m) : radius(r), position(p), material(m) {}

    // 计算射线与球面是否相交，若相交得出距离t。方程：(O + tD - p)^2 = r^2
    double intersect(const Ray &ray) const {
        Vec o_p = position - ray.origin;
        double b = o_p.dot(ray.direction);
        double delta = b * b - o_p.dot(o_p) + radius * radius;
        double delta_sqrt = sqrt(delta);
        const double eps = 1e-4;

        if(delta < 0)
            return 1e10;
        else{// 返回较小的正解
            return b - delta_sqrt > eps ? b - delta_sqrt : b + delta_sqrt > eps ? b + delta_sqrt : 1e10;
        }
    }
};

static Sphere spheres[] = {
        // 光源 由左到右
        Sphere(0.113, Vec(-3.36951, 6.43046, 3.03575), lambert2SG),
        Sphere(0.331, Vec(-0.01048, 6.43046, 3.03575), lambert2SG),
        Sphere(0.600, Vec( 3.63440, 6.43046, 3.03575), lambert2SG),
        Sphere(0.890, Vec( 7.10389, 6.43046, 3.03575), lambert2SG),
};

static Face faces[26];

// 计算光线与哪个面相交。
inline bool intersect(const Ray &ray, double &t, Mtl &material, Vec &vec, bool &isFace) {
    double inf = t = 1e10;

    for (Face *f = faces; f != faces + sizeof(faces) / sizeof(Face); ++f) {
        double d = f->intersect(ray);
        if (d < t) {
            t = d;
            material = f->material;
            vec = f->n;
            isFace = true;
        }
    }

    for (Sphere *s = spheres; s != spheres + sizeof(spheres) / sizeof(Sphere); ++s) {
        double d = s->intersect(ray);
        if(d < t) {
            t = d;
            material = s->material;
            vec = s->position;
            isFace = false;
        }
    }
    return t < inf;
}

static Vec radiance(const Ray &ray, int depth, unsigned short *xi) {
    double t;           // 距离
    Mtl material;       // 材料属性
    Vec vec;            // 用来计算法向量
    bool isFace;        // 判断是否是平面或球面

    if(!intersect(ray, t, material, vec, isFace))       // 若没有相交返回零。
        return Vec();
    else {
        int newDepth = depth + 1;
        bool isMaxDepth = newDepth > 100;       // 最大光线深度为100。
        bool isUseRR = newDepth > 5;            // 采用俄罗斯轮盘赌停止迭代
        bool isRR = isUseRR && 2 * erand48(xi) < (material.Kd + material.Ks).max();     // 若为真，光线存活。


        if (isMaxDepth || (isUseRR && !isRR)){
            return material.Ka;
        }
        else {
            Vec x = ray.origin + ray.direction * t;
            Vec n = isFace ? vec : (x - vec).normalization();
            Vec nl = n.dot(ray.direction) < 0 ? n : n * -1;

            double r1 = 2 * M_PI * erand48(xi);
            double r2 = erand48(xi);
            double r2s = sqrt(r2);

            // 计算标准正交坐标系
            Vec w = nl;
            Vec wo = fabs(w.x) > .1 ? Vec(0, 1) : Vec(1);
            Vec u = (wo.cross(w)).normalization();
            Vec v = w.cross(u);

            // 在半球面内产生随机均匀分布的射线
            Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normalization();
            Vec r = (u * cos(r1 + M_PI) * r2s + v * sin(r1 + M_PI) * r2s + w * sqrt(1 - r2)).normalization();

            Vec h = (d - ray.direction).normalization();
            double nh = nl.dot(h) > 0 ? nl.dot(h) : 0;

            return material.Ka + (material.Kd + material.Ks * pow(nh, material.ns)) * radiance(Ray(x, d), newDepth, xi);
        }

    }
}

inline double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

// gamma矫正
inline int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}

// 读取obj文件
void ReadObj() {
    Vec vertex[40];
    int vertex_index = 0 , faces_index = 0;
    char line[100];
    char sep[] = " \n";
    char *p, material[20];

    FILE *read_obj = fopen("scene02.obj", "r");
    if(!read_obj)
        cout << "FILE NOT OPEN!" << endl;


    while(fgets(line, 50, read_obj)) {
        if (line[0] == 'v') {       // 坐标点
            if (line[1] == ' '){
                strtok(line, sep);
                vertex[vertex_index].x = atof(strtok(nullptr, sep));
                vertex[vertex_index].y = atof(strtok(nullptr, sep));
                vertex[vertex_index].z = atof(strtok(nullptr, sep));
                vertex_index++;
            }
        }
        if (line[0] == 'u') {       // 材料
            strtok(line, sep);
            p = strtok(nullptr, sep);
            strcpy(material, p);
        }
        if (line[0] == 'f') {       // 面属性
            strtok(line, sep);
            faces[faces_index] = Face(vertex[(int)atof(strtok(nullptr, sep)) - 1],
                                      vertex[(int)atof(strtok(nullptr, sep)) - 1],
                                      vertex[(int)atof(strtok(nullptr, sep)) - 1],
                                      vertex[(int)atof(strtok(nullptr, sep)) - 1],
                                      mtl(material));
            faces_index++;
        }
    }
    fclose(read_obj);
}

int main(int argc, char *argv[]) {

    double start = omp_get_wtime();
    const int width = 800;
    const int height = 600;
    const int samples = argc == 2 ? atoi(argv[1]) : 100;
    const Ray camera(Vec(0, 20, 40), Vec(0, -0.5, -1).normalization());     // 相机位置，视角
    const Vec cx(width * .5 / height);
    const Vec cy = (cx.cross(camera.direction)).normalization() * .5;
    auto *color = new Vec[width * height];

    ReadObj();

#pragma omp parallel for schedule(dynamic, 1)       // openMP

    for (unsigned short y = 0; y < height; y++) {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samples, 100.*y/(height - 1));
        unsigned short xi[3] = {0, 0, y*y*y};
        for (int x = 0; x < width; x++) {
            const int i = (height - y - 1) *  width + x;
            Vec r = Vec();

            for (int s = 0; s < samples; s++) {
                double r1 = 2 * erand48(xi), dx = r1 < 2 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);     // x方向偏移
                double r2 = 2 * erand48(xi), dy = r2 < 2 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);     // y方向偏移
                Vec d = cx * (((.5 + dx) / 2 + x) / width - .5) +
                        cy * (((.5 + dy) / 2 + y) / height - .5) + camera.direction;

                r = r + radiance(Ray(camera.origin + d * 15, d.normalization()), 0, xi) * (1.0 / samples);
            }

            color[i] = color[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));
        }
    }
    cout << endl << (omp_get_wtime() - start) / 60 << " min" << endl;

    FILE *draw_pixels = fopen("ppm/scene02.ppm", "w");
    fprintf(draw_pixels, "P3\n%d %d \n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(draw_pixels, "%d %d %d\n", toInt(color[i].x), toInt(color[i].y), toInt(color[i].z));
    fclose(draw_pixels);
}
