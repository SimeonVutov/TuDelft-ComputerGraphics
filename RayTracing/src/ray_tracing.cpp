#include "ray_tracing.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
#include <framework/opengl_includes.h>

DISABLE_WARNINGS_PUSH()
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
#include <glm/gtx/string_cast.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>

// This method is based on Barycentric Coordinates
// We express p as the other vertices
// p=v0+β*e0+γ*e1
// equivalently to v = β*e0+γ*e1 where v = p - v0
bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    glm::vec3 e0 = v1 - v0;
    glm::vec3 e1 = v2 - v0;
    glm::vec3 v = p - v0;
    
    // We need to find beta and gamma
    // Since we cannot devide by vectors we will project everything using dotProd
    // So we just "multiply" by vector e0: v*e0 = β*(e0*e0) + γ*(e1*e0) 
    // We do the same with e1
    // Then we get two equation:
    // d20 =  β*d00 + γ*d01
    // d21 =  β*d01 + γ*d11

    float d00 = glm::dot(e0, e0);
    float d01 = glm::dot(e0, e1);
    float d11 = glm::dot(e1, e1);
    float d20 = glm::dot(v, e0);
    float d21 = glm::dot(v, e1);
    
    // We solve the system using determinants
    // We have the matrix d00 d01
    //                    d01 d11 -> det = d00d11 - d01^2
    //                  
    float denom = d00 * d11 - d01 * d01;
    
    // if this determinant is than the triangle is degenerate -> all points lie on a line
    if(denom <= 0.0f) {
        return false;
    }

    //Now we use the Cramer rule to solve the system:
    float beta = (d11 * d20 - d01 * d21) / denom;
    float gamma = (d00 * d21 - d01 * d20) / denom;

    return (beta >= 0.0f) && (gamma >= 0.0f) && (beta + gamma <= 1.0f);
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float denom = glm::dot(ray.direction, plane.normal);
    
    if(denom == 0) {
        return false;
    }

    float t = (plane.D - glm::dot(ray.origin, plane.normal)) / denom;

    if (t > 0.0f && t < ray.t) {
        ray.t = t;
        return true;
    }

    return false;
}

/// Input: the three vertices of the triangle
/// Output: the plane (represented via its normal vector and distance from the origin) defined by the input vertices
Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 A = v1 - v0;
    glm::vec3 B = v1 - v2;

    plane.normal = glm::normalize(glm::cross(A, B));
    plane.D = glm::dot(v0, plane.normal);

    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if ray intersects triangle defined by input vertices, then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray)
{
    Plane plane = trianglePlane(v0, v1, v2);
    Ray temp = ray;

    if(!intersectRayWithPlane(plane, temp)) return false;
    
    glm::vec3 hitPoint = temp.origin + temp.t * temp.direction;

    if(pointInTriangle(v0, v1, v2, plane.normal, hitPoint) && temp.t < ray.t) {
        ray.t = temp.t;
        return true;
    }

    return false;
}

/*
 * INTERSECT RAY WITH GHOSTS
 * Given the vertices of all ghosts and the hulls of all ghosts,
 * check if the ray intersects with any of their triangles
 */
bool intersectRayWithGhosts(const std::span<Vertices> vertices, const std::span<Hull> hulls, Ray &ray) {
    bool hit = false;
    
    for(size_t i = 0; i < hulls.size(); i++) {
        Hull &currentHull = hulls[i];
        Vertices currentGhostVertices = vertices[i];

        for(const auto &segment : currentHull) {
            for(const auto &face : segment) {
                glm::vec3 v0 = currentGhostVertices[face.x];
                glm::vec3 v1 = currentGhostVertices[face.y];
                glm::vec3 v2 = currentGhostVertices[face.z];

                if(intersectRayWithTriangle(v0, v1, v2, ray)) {
                    hit = true;
                }
            }
        }
    }

    return hit;
}

/*
 *
 * TEST SCENE
 * If you want to test your triangle intersections without using the ghost hull,
 * feel free to implement a test scene here and use the UI option test scene.
 * This part serves the purpose of simplifying your debugging and is not required.
 */

TestScene generateTestScene () {
    // implement this properly if you need it
    return {
        .vertices = {{1.0, 1.0, 1.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}},
        .triangles = {{.vertex_indices= {0, 1, 2}, .color={0.5, 0.5, 0.5}}}
    };
}

void drawTestScene (const TestScene &testScene) {
    // implement this properly if you need it
    glBegin(GL_TRIANGLES);
    glColor3fv(glm::value_ptr(testScene.triangles[0].color));
    glVertex3fv(glm::value_ptr(testScene.vertices[testScene.triangles[0].vertex_indices[0]]));
    glVertex3fv(glm::value_ptr(testScene.vertices[testScene.triangles[0].vertex_indices[1]]));
    glVertex3fv(glm::value_ptr(testScene.vertices[testScene.triangles[0].vertex_indices[2]]));
    glEnd();
}

bool intersectTestScene (const TestScene &testScene, Ray &ray) {
    // implement this properly if you need it
    return intersectRayWithTriangle(testScene.vertices[testScene.triangles[0].vertex_indices[0]],
                                    testScene.vertices[testScene.triangles[0].vertex_indices[1]],
                                    testScene.vertices[testScene.triangles[0].vertex_indices[2]],
                                    ray);
}
