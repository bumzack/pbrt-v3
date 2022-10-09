
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// main/pbrt.cpp*
#include "pbrt.h"
#include "api.h"
#include "parser.h"
#include "parallel.h"
#include <glog/logging.h>
#include "transform.h"
#include "shapes/curve.h"
#include "shapes/loopsubdiv.h"
#include "paramset.h"
#include "geometry.h"

using namespace pbrt;
using namespace std;

#include <iostream>
#include <iomanip>
#include <memory>

#include <iostream>
#include <iomanip>
#include <memory>

void test_curve() {
    std::cout << endl << "================================================================================" << endl;
    cout << "    Curve Shape " << endl;

//    Curve(const Transform *ObjectToWorld, const Transform *WorldToObject,
//    bool reverseOrientation, const std::shared_ptr<CurveCommon> &common,
//    Float uMin, Float uMax)

    const Float obj2world[4][4] = {{1.0, 2.0, 3.0, 4.0},
                                   {1.0, 2.0, 3.0, 4.0},
                                   {1.0, 2.0, 3.0, 4.0},
                                   {1.0, 2.0, 3.0, 4.0}
    };

    const Float w2o[4][4] = {{1.0, 2.0, 3.0, 4.0},
                             {1.0, 2.0, 3.0, 4.0},
                             {1.0, 2.0, 3.0, 4.0},
                             {1.0, 2.0, 3.0, 4.0}
    };

    const pbrt::Transform objectToWorld = pbrt::Transform(obj2world);
    const pbrt::Transform worldToObject = pbrt::Transform(w2o);
    Float umin = 0.0;
    Float umax = 2.0;
    bool reverseOrientation = false;


    CurveType type = CurveType::Ribbon;
    Float w0 = 0.5;
    Float w1 = 1.5;
    const Normal3f norm = Normal3f(Vector3f(1.0, 2.0, 3.0));
    const Point3f cp[4] = {Point3f(0.0f, 1.0f, 1.0f),
                           Point3f(1.0f, 1.0f, 2.0f),
                           Point3f(2.0f, 2.0f, 2.0f),
                           Point3f(3.0f, 3.0f, 3.0f)
    };


    const CurveCommon cc = CurveCommon{cp, w0, w1, type, &norm};
    std::shared_ptr<CurveCommon> common = std::make_shared<CurveCommon>(cc);
//    const Transform *ObjectToWorld, const Transform *WorldToObject,
//    bool reverseOrientation, const std::shared_ptr<CurveCommon> &common,
//    Float uMin, Float uMax


    pbrt::Curve c = Curve(&objectToWorld, &worldToObject, reverseOrientation, common, umin, umax);


//    std::array<Point3f, 4> cp;
//    Bounds3f b = BoundCubicBezier(pstd::MakeConstSpan(cp), 0.f, 1.f);
//    b = Expand(b, 1e-3 * Length(b.Diagonal()));
//    for (Float u = 0; u <= 1.f; u += 1.f / 1024.f) {
//        Point3f p = EvaluateCubicBezier(pstd::MakeConstSpan(cp), u);
//        bool inside = Inside(p, b);
//        if (inside) {
////            cout << "inside   " << p << " @ u = " << u << " not in " << b << endl;
//        } else {
////            cout << "outside   " << p << " @ u = " << u << " not in " << b << endl;
//        }
//    }

//
//    RNG rng;
//    for (int i = 0; i < 1000; ++i) {
//        for (int j = 0; j < 4; ++j)
//            for (int c = 0; c < 3; ++c)
//                cp[j][c] = -5.f + 10.f * rng.Uniform<Float>();
//
//        Bounds3f b = BoundCubicBezier(pstd::MakeConstSpan(cp), 0.f, 1.f);
//        b = Expand(b, 1e-3 * Length(b.Diagonal()));
//        for (Float u = 0; u <= 1.f; u += 1.f / 1024.f) {
//            Point3f p = EvaluateCubicBezier(pstd::MakeConstSpan(cp), u);
//            cout << (Inside(p, b)) << p << " @ u = " << u << " not in " << b << endl;
//        }
//    }
//
    cout << endl << "================================================================================" << endl;

}

void test_clamp() {
    cout << endl << "================================================================================" << endl;
    cout << "clamp     test" << endl << endl;

    Float low = 0.0;
    Float high = 2.0;

    Float values[10] = {-0.000001, 0.0, 0.25, 0.75, 0.5, 1.0, 1.5, 1.75, 2.0, 3.0};

    for (int i = 0; i < 10; i++) {
        Float res = Clamp(values[i], low, high);
        cout << std::fixed << std::setw(11) << std::setprecision(6) << "low " << low << "   high:  " << high
             << "   value " << values[i] << "  clamped to  ==>      " << res << endl;
    }
    cout << endl << "================================================================================" << endl;

}

void test_lerp_f64() {
    cout << endl << "================================================================================" << endl;
    cout << "lerp f64    test" << endl << endl;

    Float min = 0.0;
    Float max = 2.0;

    Float arr[8] = {0.0, 0.25, 0.75, 0.5, 1.0, 1.5, 1.75, 2.0};

    for (int i = 0; i < 8; i++) {
        Float res = Lerp(arr[i], min, max);
        cout << std::fixed << std::setw(11) << std::setprecision(6) << "min " << min << "   max:  " << max
             << " value    " << arr[i] << "  res     " << res << endl;
    }
    cout << endl << "================================================================================" << endl;
}


void test_lerp_vec() {
    cout << endl << "================================================================================" << endl;
    cout << "lerp point3d    test" << endl << endl;

    Point3f a = Point3f(1.0f, 2.0f, 3.0f);
    Point3f b = Point3f(14.0f, -15.0f, 16.0f);

    Float arr[8] = {0.0, 0.25, 0.75, 0.5, 1.0, 1.5, 1.75, 2.0};

    for (int i = 0; i < 8; i++) {
        Point3f res = Lerp(arr[i], a, b);
        cout << std::fixed << std::setw(11) << std::setprecision(6) << "a " << a << "   b:  " << b
             << " value    " << arr[i] << "  res     " << res << endl;
    }
    cout << endl << "================================================================================" << endl;

}


void test_blossom_bezier() {
    cout << endl << "================================================================================" << endl;
    cout << "bloosom bezier test" << endl << endl;

    const Point3f cp[4] = {Point3f(0.0f, 1.0f, 1.0f),
                           Point3f(1.0f, 1.0f, 2.0f),
                           Point3f(2.0f, 2.0f, 2.0f),
                           Point3f(3.0f, 3.0f, 3.0f)
    };

    Float u0 = 0.0;
    Float u1 = 0.25;
    Float u2 = 1.0;

    Point3f f1 = BlossomBezier(cp, u0, u0, u0);
    Point3f f2 = BlossomBezier(cp, u0, u0, u1);
    Point3f f3 = BlossomBezier(cp, u1, u1, u2);
    Point3f f4 = BlossomBezier(cp, u2, u2, u2);

    int i = 1;
    for (Point3f p: cp) {
        cout << "point  " << i << "   x " << p.x << "  y   " << p.y << "   p.z  " << p.z << endl;
        i++;
    }
    cout << "u0  " << u0 << endl;
    cout << "u1  " << u1 << endl;
    cout << "u2  " << u2 << endl;

    cout << "u0 " << u0 << "   u0  " << u0 << "  u0   " << u0 << "   f1       " << f1 << endl;
    cout << "u0 " << u0 << "   u0  " << u0 << "  u1   " << u1 << "   f2       " << f2 << endl;
    cout << "u1 " << u1 << "   u1  " << u1 << "  u2   " << u2 << "   f3       " << f3 << endl;
    cout << "u2 " << u2 << "   u2  " << u2 << "  u2   " << u2 << "   f4       " << f4 << endl;
    cout << endl << "================================================================================" << endl;

}

//
//template<class T>
//std::ostream &
//operator<<(std::ostream &s, const Vector3<T> &v) {
//    s << std::fixed << std::setw(11) << std::setprecision(4) <<v. x << ", " << v.y << " " << v.z;
//    return s;
//}


void test_coordinate_system() {
    cout << endl << "================================================================================" << endl;
    cout << "coordinate_system   " << endl << endl;

    Vector3f v1[4] = {Vector3f(1.0, 0.0, 0.0),
                      Vector3f(1.0, 1.0, 0.0),
                      Vector3f(1.0, 0.0, 1.0),
                      Vector3f(1.0, 1.0, 1.0)};

    Vector3f v2[4] = {Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0)};

    Vector3f v3[4] = {Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0),
                      Vector3f(0.0, 0.0, 0.0)};

    for (int i = 0; i < 4; i++) {
        CoordinateSystem(v1[i], &v2[i], &v3[i]);

        cout << "v1 " << v1[i].x << "     " << v1[i].y << "    " << v1[i].z << endl;
        cout << "v2 " << v2[i].x << "     " << v2[i].y << "    " << v2[i].z << endl;
        cout << "v3 " << v3[i].x << "     " << v3[i].y << "    " << v3[i].z << endl;

        cout << "v1 " << v1[i] << "    v2   " << v2[i] << "  v3    " << v3[i] << endl;

        cout << endl;
    }


    cout << endl << "================================================================================" << endl;

}


void test_look_at() {
    cout << endl << "================================================================================" << endl;
    cout << "   LookAt   " << endl << endl;

    Point3f pos = Point3f(1.0, 2.0, 3.0);
    Point3f look = Point3f(5.0, 6.0, 7.0);
    Vector3f up[3] = {Vector3f(1.0, 0.0, 0.0),
                      Vector3f(0.0, 1.0, 0.0),
                      Vector3f(0.0, 0.0, 1.0)
    };

    cout << "pos:   " << pos.x << "   " << pos.y << "    " << pos.z << endl;
    cout << "look:   " << look.x << "   " << look.y << "    " << look.z << endl;

    for (int i = 0; i < 3; i++) {
        Transform t = LookAt(pos, look, up[i]);

        Matrix4x4 m = t.GetMatrix();
        Matrix4x4 m_inv = t.GetInverseMatrix();

        cout << "up:   " << up[i].x << "   " << up[i].y << "    " << up[i].z << endl;

        cout << "m     :   " << m << endl;
        cout << "m_inv :   " << m_inv << endl;

        cout << endl;

    }


    cout << endl << "================================================================================" << endl;

}


void test_union() {
    cout << endl << "================================================================================" << endl;
    cout << "   Union   " << endl << endl;

    Point3f p0 = Point3f(0.0, 0.0, 0.0);
    Point3f p1 = Point3f(2.0, 3.0, 4.0);

    Point3f p2 = Point3f(0.5, -0.5, -1.5);
    Point3f p3 = Point3f(1.5, 13.0, 3.0);

    Bounds3f bb1 = Bounds3f(p0, p1);
    Bounds3f bb2 = Bounds3f(p2, p3);

    Bounds3f u = Union(bb1, bb2);


    cout << "bb1.min:   " << bb1.pMin.x << "   " << bb1.pMin.y << "    " << bb1.pMin.z << endl;
    cout << "bb1.max:   " << bb1.pMax.x << "   " << bb1.pMax.y << "    " << bb1.pMax.z << endl;

    cout << "bb2.min:   " << bb2.pMin.x << "   " << bb2.pMin.y << "    " << bb2.pMin.z << endl;
    cout << "bb2.max:   " << bb2.pMax.x << "   " << bb2.pMax.y << "    " << bb2.pMax.z << endl;

    cout << "union.min:   " << u.pMin.x << "   " << u.pMin.y << "    " << u.pMin.z << endl;
    cout << "union.max:   " << u.pMax.x << "   " << u.pMax.y << "    " << u.pMax.z << endl;


    cout << endl << "================================================================================" << endl;

}


void test_expand() {
    cout << endl << "================================================================================" << endl;
    cout << "   Expand   " << endl << endl;

    Point3f p0 = Point3f(0.0, 0.0, 0.0);
    Point3f p1 = Point3f(2.0, 3.0, 4.0);

    Bounds3f bb1 = Bounds3f(p0, p1);

    Bounds3f e = Expand(bb1, 0.75);

    cout << "bb1.min:   " << bb1.pMin.x << "   " << bb1.pMin.y << "    " << bb1.pMin.z << endl;
    cout << "bb1.max:   " << bb1.pMax.x << "   " << bb1.pMax.y << "    " << bb1.pMax.z << endl;

    cout << "expand.min:   " << e.pMin.x << "   " << e.pMin.y << "    " << e.pMin.z << endl;
    cout << "expand.max:   " << e.pMax.x << "   " << e.pMax.y << "    " << e.pMax.z << endl;

    cout << endl << "================================================================================" << endl;
}


//void test_ray_bounds() {
//    cout << endl << "================================================================================" << endl;
//    cout << "   rayBounds   " << endl << endl;
//
//    Ray r = Ray(Point3f(1.0, 2.0, 3.0), Vector3f(-2.0, -2.0, -1.5));
//    Float tMax = 1.0;
//    Bounds3f rayBounds(Point3f(0, 0, 0), Point3f(0, 0, r.d.Length() * tMax));
//
//
//    cout << "ray.p:   " << r.o.x << "   " << r.o.y << "    " << r.o.z << endl;
//    cout << "ray.d:   " << r.d.x << "   " << r.d.y << "    " << r.d.z << "    length(r.d)  " << r.d.Length() << endl;
//
//
//    cout << "rayBounds.min:   " << rayBounds.pMin.x << "   " << rayBounds.pMin.y << "    " << rayBounds.pMin.z << endl;
//    cout << "rayBounds.max:   " << rayBounds.pMax.x << "   " << rayBounds.pMax.y << "    " << rayBounds.pMax.z << endl;
//
//    cout << endl << "================================================================================" << endl;
//}


void test_overlaps() {
    cout << endl << "================================================================================" << endl;
    cout << "   Overlaps   " << endl << endl;

    Point3f p0 = Point3f(0.0, 0.0, 0.0);
    Point3f p1 = Point3f(2.0, 3.0, 4.0);

    Point3f p2 = Point3f(0.5, 0.5, 0.5);
    Point3f p3 = Point3f(1.5, 13.0, 3.0);

    Point3f p4 = Point3f(10.5, 10.5, 10.5);
    Point3f p5 = Point3f(1.5, 13.0, 13.0);

    Bounds3f bb1 = Bounds3f(p0, p1);
    Bounds3f bb2 = Bounds3f(p2, p3);
    Bounds3f bb3 = Bounds3f(p4, p5);

    bool overlaps = Overlaps(bb1, bb2);
    bool overlaps2 = Overlaps(bb1, bb3);

    cout << "bb1.min:   " << bb1.pMin.x << "   " << bb1.pMin.y << "    " << bb1.pMin.z << endl;
    cout << "bb1.max:   " << bb1.pMax.x << "   " << bb1.pMax.y << "    " << bb1.pMax.z << endl << endl;

    cout << "bb2.min:   " << bb2.pMin.x << "   " << bb2.pMin.y << "    " << bb2.pMin.z << endl;
    cout << "bb2.max:   " << bb2.pMax.x << "   " << bb2.pMax.y << "    " << bb2.pMax.z << endl << endl;

    cout << "bb3.min:   " << bb3.pMin.x << "   " << bb3.pMin.y << "    " << bb3.pMin.z << endl;
    cout << "bb3.max:   " << bb3.pMax.x << "   " << bb3.pMax.y << "    " << bb3.pMax.z << endl << endl;

    cout << "overlaps:   " << std::boolalpha << overlaps << endl;
    cout << "overlaps2:   " << std::boolalpha << overlaps2 << endl;

    cout << endl << "================================================================================" << endl;
}


void test_intersect_curve(const std::string curve_type) {
    cout << endl << "================================================================================" << endl;
    cout << "intersect curve" << endl << endl;

    const Float obj2world[4][4] = {{1.0, 2.0, 3.0, 1.0},
                                   {2.0, 7.0, 5.0, 4.0},
                                   {3.0, 6.0, 2.0, 1.0},
                                   {4.0, 5.0, 1.0, 1.0}
    };

    const Float w2o[4][4] = {{1.0, 3.0, 9.0, 4.0},
                             {1.0, 4.0, 8.0, 3.0},
                             {1.0, 1.0, 7.0, 4.0},
                             {1.0, 5.0, 6.0, 5.0}
    };

    const pbrt::Transform objectToWorld = pbrt::Transform(obj2world);
    const pbrt::Transform worldToObject = pbrt::Transform(w2o);
    Float umin = 0.0;
    Float umax = 2.0;
    bool reverseOrientation = false;

    CurveType type = CurveType::Ribbon;
    Float w0 = 0.5;
    Float w1 = 1.5;
    const Normal3f norm = Normal3f(Vector3f(1.0, 2.0, 3.0));
    const Point3f cp[4] = {Point3f(0.0f, 1.0f, 1.0f),
                           Point3f(1.0f, 1.0f, 2.0f),
                           Point3f(2.0f, 2.0f, 2.0f),
                           Point3f(3.0f, 3.0f, 3.0f)
    };

    const CurveCommon cc = CurveCommon{cp, w0, w1, type, &norm};
    std::shared_ptr<CurveCommon> common = std::make_shared<CurveCommon>(cc);

    pbrt::Curve c = Curve(&objectToWorld, &worldToObject, reverseOrientation, common, umin, umax);

    ParamSet paramSet = ParamSet();

    auto p1 = Point3f(1.0, 2.0, 3.0);
    std::unique_ptr<Point3f[]> pts(new Point3f[1]);
    //  std::unique_ptr<Point3f[]> testData = std::make_unique<Point3f[]>(1);
    pts[0] = p1;
    paramSet.AddPoint3f("P", std::move(pts), 1);

    std::unique_ptr<std::string[]> basis(new std::string[1]);
    basis[0] = curve_type;
    paramSet.AddString("basis", std::move(basis), 1);

    std::unique_ptr<std::string[]> t(new std::string[1]);
    t[0] = "flat";
    paramSet.AddString("type", std::move(t), 1);

    std::unique_ptr<Normal3f[]> n(new Normal3f[1]);
    n[0] = Normal3f(1.0, 2.0, 3.0);
    // paramSet.AddNormal3f("N", std::move(n), 1);

    std::vector<std::shared_ptr<Shape>> shape = CreateCurveShape(&objectToWorld,
                                                                 &worldToObject,
                                                                 reverseOrientation,
                                                                 paramSet);

    int cnt = 12;
    Ray rays[] = {Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, 0.02, 0.0)),
                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, -0.02, 0.0)),
                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, -0.02, 0.02)),
                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, -0.02, -0.02)),

                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),
                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),
                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),
                  Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, 0.0, 0.0)),

                  Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1.0, 0.0, 0.0)),
                  Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1.0, 0.0, 0.0)),
                  Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1., 0.0, 0.0)),
                  Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1.0, 0.0, 0.0)),
    };

    Float tHit;
    SurfaceInteraction si = SurfaceInteraction();
    bool testAlphaTexture = false;

    for (int i = 0; i < cnt; i++) {
        const Float obj2world[4][4] = {{1.0, 2.0, 3.0, 1.0},
                                       {2.0, 7.0, 5.0, 4.0},
                                       {3.0, 6.0, 2.0, 1.0},
                                       {4.0, 5.0, 1.0, 1.0}
        };

        const Float w2o[4][4] = {{1.0, 3.0, 9.0, 4.0},
                                 {1.0, 4.0, 8.0, 3.0},
                                 {1.0, 1.0, 7.0, 4.0},
                                 {1.0, 5.0, 6.0, 5.0}
        };

        const pbrt::Transform objectToWorld = pbrt::Transform(obj2world);
        const pbrt::Transform worldToObject = pbrt::Transform(w2o);


        std::vector<std::shared_ptr<Shape>> shape = CreateCurveShape(&objectToWorld,
                                                                     &worldToObject,
                                                                     reverseOrientation,
                                                                     paramSet);


        c.Intersect(rays[i], &tHit, &si, testAlphaTexture);

        cout << "intersect_curve " << curve_type << "  ray  point " << rays[i].o << "   direction " << rays[i].d
             << "   i = " << i << endl;

        cout << "thit   " << tHit << endl;
        cout << "si uv  " << si.uv << endl;
        cout << "si dndu " << si.dndu << endl;
        cout << "si dndv " << si.dndv << endl;
        cout << "si dpdu " << si.dpdu << endl;
        cout << "si dpdv " << si.dpdv << endl;
        cout << "si dpdx " << si.dpdx << endl;
        cout << "si dpdy " << si.dpdy << endl;
        cout << "si dvdx " << si.dvdx << endl;
        cout << "si dvdy " << si.dvdy << endl;
        cout << "si faceIndex " << si.faceIndex << endl;
        cout << "si.n " << si.n << endl;
        cout << "si.p " << si.p << endl;
        cout << endl;
    }

    cout << endl << "================================================================================" << endl;
}

void test_loopsubdiv() {
    cout << endl << "================================================================================" << endl;
    cout << "loop sub div" << endl << endl;


    const Float obj2world[4][4] = {{1.0, 2.0, 3.0, 1.0},
                                   {2.0, 7.0, 5.0, 4.0},
                                   {3.0, 6.0, 2.0, 1.0},
                                   {4.0, 5.0, 1.0, 1.0}
    };

    const Float w2o[4][4] = {{1.0, 3.0, 9.0, 4.0},
                             {1.0, 4.0, 8.0, 3.0},
                             {1.0, 1.0, 7.0, 4.0},
                             {1.0, 5.0, 6.0, 5.0}
    };

    const pbrt::Transform objectToWorld = pbrt::Transform(obj2world);
    const pbrt::Transform worldToObject = pbrt::Transform(w2o);

    bool reverseOrientation = false;

    ParamSet paramSet = ParamSet();

    auto p1 = Point3f(1.0, 2.0, 3.0);
    auto p2 = Point3f(2.0, 2.0, 4.0);
    auto p3 = Point3f(3.0, 2.0, 1.0);
    std::unique_ptr<Point3f[]> vertices(new Point3f[3]);
    vertices[0] = p1;
    vertices[1] = p2;
    vertices[2] = p3;
    paramSet.AddPoint3f("P", std::move(vertices), 3);

    std::unique_ptr<int[]> indices(new int[1]);
    indices[0] = 3;
    paramSet.AddInt("indices", std::move(indices), 1);

    std::unique_ptr<int[]> level(new int[1]);
    level[0] = 1;
    paramSet.AddInt("levels", std::move(level), 1);


    std::vector<std::shared_ptr<Shape>> shapes =
            CreateLoopSubdiv(&objectToWorld, &worldToObject, reverseOrientation,
                             paramSet);


    for (auto s: shapes) {

        int cnt = 12;
        Ray rays[] = {Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, 0.02, 0.0)),
                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, -0.02, 0.0)),
                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, -0.02, 0.02)),
                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.0, -0.02, -0.02)),

                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),
                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),
                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),
                      Ray(Point3f(-10.0, 0.0, 0.0), Vector3f(1.02, -0.02, 0.02)),

                      Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1.0, 0.0, 0.0)),
                      Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1.0, 0.0, 0.0)),
                      Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1., 0.0, 0.0)),
                      Ray(Point3f(-10.0, 1.0, 0.0), Vector3f(1.0, 0.0, 0.0)),
        };

        Float tHit;
        SurfaceInteraction si = SurfaceInteraction();
        bool testAlphaTexture = false;

        for (int i = 0; i < cnt; i++) {
            s->Intersect(rays[i], &tHit, &si, testAlphaTexture);

            cout << "loop sub div   ray  point " << rays[i].o << "   direction " << rays[i].d << "   i = " << i << endl;

            cout << "thit   " << tHit << endl;
            cout << "si uv  " << si.uv << endl;
            cout << "si dndu " << si.dndu << endl;
            cout << "si dndv " << si.dndv << endl;
            cout << "si dpdu " << si.dpdu << endl;
            cout << "si dpdv " << si.dpdv << endl;
            cout << "si dpdx " << si.dpdx << endl;
            cout << "si dpdy " << si.dpdy << endl;
            cout << "si dvdx " << si.dvdx << endl;
            cout << "si dvdy " << si.dvdy << endl;
            cout << "si faceIndex " << si.faceIndex << endl;
            cout << "si.n " << si.n << endl;
            cout << "si.p " << si.p << endl;
            cout << endl;
        }
    }

    cout << endl << "================================================================================" << endl;
}


void test_log2() {
    cout << endl << "================================================================================" << endl;
    cout << "loop sub div" << endl << endl;

    srand(time(nullptr));

    auto Log2 = [](float v) -> int {
        if (v < 1) return 0;
        uint32_t bits = FloatToBits(v);
        // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
        // (With an additional add so get round-to-nearest rather than
        // round down.)
        cout << "v = " << v << "  (bits >> 23)   " << (bits >> 23) << "    (1 << 22)   " << (1 << 22) << endl;
        cout << "v = " << v << "  (bits & (1 << 22) ? 1 : 0)   " << (bits & (1 << 22) ? 1 : 0) << endl;
        cout << "v = " << v << "  bits & (1 << 22)   " << (bits & (1 << 22)) << endl;

        if ((bits & (1 << 22))) {
            cout << " (bits & (1 << 22)) = " << (bits & (1 << 22)) << "   == true" << endl;
        } else {
            cout << " (bits & (1 << 22))  =  " << (bits & (1 << 22)) << " == false" << endl;

        }
        return (bits >> 23) - 127 + (bits & (1 << 22) ? 1 : 0);
    };

    // Compute log base 4 by dividing log2 in half.
//
//    for (int i = 0; i < 10; i++) {
//        Float x = rand() % 100 + 1;
//
//        int l2 = Log2(x);
//        cout << std::fixed << std::setw(11) << std::setprecision(6) << "x  " << log2(x) << "   log2(x):  " << l2
//             << endl;
//    }

Float   x = 7.0;
    int l2 = Log2(x);
    cout << std::fixed << std::setw(11) << std::setprecision(6) << " x " << x << "   log2(x):  " << l2 << endl;


    Float   x1 = 46.0;
    int l2a = Log2(x1);
    cout << std::fixed << std::setw(11) << std::setprecision(6) << " x " << x1 << "   log2(x):  " << l2a << endl;

    cout << endl << "================================================================================" << endl;
}


// main program
int main(int argc, char *argv[]) {
    test_lerp_f64();
    test_lerp_vec();
    test_blossom_bezier();
    test_coordinate_system();
    test_look_at();

    test_union();
    test_expand();
    // test_ray_bounds();
    test_overlaps();
    test_clamp();
    test_log2();
//
//    test_intersect_curve("bezier");
//    test_intersect_curve("bspline");
//    test_loopsubdiv();

    return 0;
}
