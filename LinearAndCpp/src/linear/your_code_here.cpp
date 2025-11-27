// Should trigger CI only in this project
// TEST
#include "linear.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/matrix_transform.hpp>
#include <glm/vec3.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/constants.hpp>
DISABLE_WARNINGS_POP()
#include <iostream>
#include <algorithm>

std::string to_string(const Matrix3& matrix)
{
    return "matrix with columns:\n Col1 " + to_string(matrix.col1) + "\n Col2 " + to_string(matrix.col2) + "\n Col3 " + to_string(matrix.col3);
}

// ==================================
// ========    Exercise 0    ========
// ==================================
// Multiplication of a vector with a scalar
glm::vec3 mul(const glm::vec3& lhs, float rhs)
{
    glm::vec3 result {lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
    return result;
}

// Dot product of two vectors
float dot3(const glm::vec3& lhs, const glm::vec3& rhs)
{
    float result;
    result = lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
    return result;
}

// Cross product of two vectors
glm::vec3 cross3(const glm::vec3& lhs, const glm::vec3& rhs)
{
    glm::vec3 result { lhs.y*rhs.z - rhs.y*lhs.z, lhs.z*rhs.x - rhs.z*lhs.x, lhs.x*rhs.y - rhs.x*lhs.y };
    return result;
}

// Length of a vector
float length(const glm::vec3& lhs)
{
    float result;
    result = glm::sqrt(lhs.x * lhs.x + lhs.y * lhs.y + lhs.z * lhs.z);
    return result;
}

// For the following exercises, a matrix is defined as:
// | m00  m01  m02 |
// | m10  m11  m12 |
// | m20  m21  m22 |
//
// The columns are stored as Vector3's, where:
// col1 = (m00, m10, m20)
// col2 = (m01, m11, m21)
// col3 = (m02, m12, m22)
//
// For a Matrix3 matrix object, these can be accessed via matrix.col1, matrix.col2, matrix.3.

// Matrix multiplication with a scalar
Matrix3 mul(const Matrix3& lhs, float rhs)
{
    Matrix3 result {};
    result.col1 = mul(lhs.col1, rhs);
    result.col2 = mul(lhs.col2, rhs);
    result.col3 = mul(lhs.col3, rhs);
    return result;
}

// Taking the transpose of a matrix means changing it's columns into rows and vice versa.
// Given a matrix of the form:
// | m00  m01  m02 |
// | m10  m11  m12 |
// | m20  m21  m22 |
// its transpose is given as:
// | m00  m10  m20 |
// | m01  m11  m21 |
// | m02  m12  m22 |
Matrix3 transpose(const Matrix3& m)
{
    Matrix3 result {};
    
    result.col1.x = m.col1.x;
    result.col1.y = m.col2.x;
    result.col1.z = m.col3.x;
    
    result.col2.x = m.col1.y;
    result.col2.y = m.col2.y;
    result.col2.z = m.col3.y;

    result.col3.x = m.col1.z;
    result.col3.y = m.col2.z;
    result.col3.z = m.col3.z;
    return result;
}

// The determinant is needed to compute the inverse of a matrix.
// If you need a Linear Algebra refresher, please check out:
// https://www.tudelft.nl/en/eemcs/study/online-education/math-explained/linear-algebra/#c144161
float determinant(const Matrix3& m)
{
    float result;
    result = m.col1.x * m.col2.y * m.col3.z +
             m.col3.x * m.col2.z * m.col1.y +
             m.col1.z * m.col2.x * m.col3.y -

             m.col1.z * m.col2.y * m.col3.x -
             m.col1.y * m.col2.x * m.col3.z -
             m.col1.x * m.col2.z * m.col3.y;
    return result;
}

// Custom function for helping inverse
float det2(float a, float b, float c, float d) {
    return a * d - b * c;
}

// Computing the inverse of the given matrix. If you implemented it correctly then matrix M multiplied
// by its inverse should give the identity matrix). More information on how to compute the inverse of a
// 3x3 matrix can be found here:
// https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
// For this method, you can assume that the given matrix is invertable, i.e., that the determinant is
// not zero. However, you are responsible for checking that all matrices passed to your function in your 
// own code are invertible.
Matrix3 inverse(const Matrix3& matrix)
{
    Matrix3 result {};
    // Hint:
    // You can follow, e.g., the method described here:
    // https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html

    // Step 0: It is probably handy to define a method that computes the determinant of a 2x2 matrix.
    // Step 1: Compute the Matrix of Minors.
    Matrix3 minors;
    minors.col1.x = det2(matrix.col2.y, matrix.col3.y, matrix.col2.z, matrix.col3.z);
    minors.col2.x = det2(matrix.col1.y, matrix.col3.y, matrix.col1.z, matrix.col3.z);
    minors.col3.x = det2(matrix.col1.y, matrix.col2.y, matrix.col1.z, matrix.col2.z);
    
    minors.col1.y = det2(matrix.col2.x, matrix.col3.x, matrix.col2.z, matrix.col3.z);
    minors.col2.y = det2(matrix.col1.x, matrix.col3.x, matrix.col1.z, matrix.col3.z);
    minors.col3.y = det2(matrix.col1.x, matrix.col2.x, matrix.col1.z, matrix.col2.z);

    minors.col1.z = det2(matrix.col2.x, matrix.col3.x, matrix.col2.y, matrix.col3.y);
    minors.col2.z = det2(matrix.col1.x, matrix.col3.x, matrix.col1.y, matrix.col3.y);
    minors.col3.z = det2(matrix.col1.x, matrix.col2.x, matrix.col1.y, matrix.col2.y);
    // Step 2: Compute the Matrix of Cofactors
    minors.col2.x *= -1;
    minors.col1.y *= -1;
    minors.col3.y *= -1;
    minors.col2.z *= -1;
    // Step 3: Adjugate the Matrix (Note that you have the transpose available from above).
    minors = transpose(minors); 
    // Step 4: Multiply by 1/determinant (Note that you have the multiplication available from above).
    float det = determinant(matrix);
    
    if(det == 0.0f) {
        return Matrix3{};
    }

    result = mul(minors, 1.0f / det);
    return result;
}


typedef struct VertexData {
    size_t index;
    float angle;
} VertexData;

bool compareByAngle(const VertexData &a, const VertexData &b) {
    return a.angle < b.angle;
}

// ==================================
// ========    Exercise 1    ========
// ==================================
// Given a set of vertices of a regular or irregular n-Gon (n>2) on a plane in 3D, we want to order them clock- or counterclockwise
// You may assume that the n-Gon is convex.
std::vector<int> orderOfnGonVertices(const std::vector<glm::vec3> nGon)
{
    // You can assume that all vertices of the n-Gon are residing in some plane in 3D, that there are at least 3 vertices in the n-Gon,
    // and that the n-Gon itself is convex.
    // Your solution here, result should have the indices of the vertices from the n-Gon in either counter- or clockwise order.
    
    std::vector<int> result;
    size_t count = nGon.size();

    if(count < 3) return result;

    glm::vec3 center{0.0f};
    for(size_t i = 0; i < count; i++) {
        center += nGon.at(i);
    }
    center /= static_cast<float>(count);

    glm::vec3 vecA = nGon[1] - nGon[0];
    glm::vec3 vecB = nGon[2] - nGon[0];

    glm::vec3 normal = glm::normalize(cross3(vecA, vecB));
    glm::vec3 refDir = glm::normalize(nGon[0] - center);

    std::vector<VertexData> vertexAngles;
    vertexAngles.reserve(count);

    for(size_t i = 0; i < count; i++) {
        glm::vec3 currDir = glm::normalize(nGon[i] - center);
        float dot = dot3(refDir, currDir);

        if(dot > 1.0f) dot = 1.0f;
        else if(dot < -1.0f) dot = -1.0f;

        float angle = std::acos(dot);

        glm::vec3 cross = cross3(refDir, currDir);
        float sign = dot3(cross, normal);

        if(sign < 0.0f) {
            angle = glm::two_pi<float>() - angle;
        }
        
        vertexAngles.push_back({i, angle});
    }

    std::sort(vertexAngles.begin(), vertexAngles.end(), compareByAngle);

    for(const VertexData &v : vertexAngles) {
        result.push_back(static_cast<int>(v.index));
    }
    
    return result;
}


// ==================================
// ========    Exercise 2    ========
// ==================================
// To visualize a given plane, we want to find four points on this plane that span a rectangle.
std::array<glm::vec3, 4> rectangleOnPlane(const Plane& plane)
{
    std::array<glm::vec3, 4> result;

    // Hint:
    // You are given the normal vector and a base point of the plane via the Plane struct.
    // (1) Find a vector that is perpendicular to the normal vector. You can assume that the normal vector is not the zero-vector. 
    //     You can use the dot product you implemented above to check that the vectors are indeed perpendicular.
    // (2) Use the cross-product to compute a second vector that is perpendicular to both the vector found in (1) and the normal 
    //     of the given plane. Again, you can check yourself using the dot product implemented above.
    // (3) As these two vectors you found span the plane, use them to offset the base point in four different directions.
    //     To be able to draw a nice rectangle from them, be sure to order them either clockwise or counter-clockwise.

    // Regarding the visuals: It's really about passing four points that lie on that plane. The rectangles can be tiny, intersect, 
    // be rotated, ... whatever, we don't care. The one thing you should check is that the four points actually lie on the plane 
    // and have the shape of (roughly) a square.
    
    glm::vec3 helper;

    if(std::abs(plane.n.x) < std::abs(plane.n.y)) {
        helper = glm::vec3{1.0f, 0.0f, 0.0f};
    }
    else {
        helper = glm::vec3{0.0f, 1.0f, 0.0f};
    }

    glm::vec3 oneDir = cross3(helper, plane.n);
    oneDir = glm::normalize(oneDir);

    glm::vec3 secDir = cross3(plane.n, oneDir);
    secDir = glm::normalize(secDir);
    
    result[0] = plane.p - oneDir - secDir;
    result[1] = plane.p + oneDir - secDir;
    result[2] = plane.p + oneDir + secDir;
    result[3] = plane.p - oneDir + secDir;
    
    return result;
}

// A solid can be given by the planes that are surrounding it. We want to find the vertices of a solid given by its planes.
std::vector<glm::vec3> verticesFromPlanes(std::span<const Plane> planes)
{
    std::vector<glm::vec3> result;

    // Hint:
    // You are given a set of planes.
    // A triple of these planes can intersect in a point.
    // There are two things to watch out for:
    // - a triple of planes might not have a common intersection.
    // - a triple of planes might have a common intersection point that is not within the solid.

    // Therefore:
    // - Iterate over all triples (p1,p2,p3) of given planes.
    // - Set up a linear system that represents the intersection of the three planes, you can ignore the case where the three planes intersect in a line.
    // - Test whether the system has a solution (you can, e.g., use the determinant for matrices).
    // - If it has a solution, compute it (e.g., via using the inverse matrix).
    // - Assuming that the normal of the plane points to the inside of the solid, test for all other planes whether the computed point of intersection is indeed within the solid.

    size_t n_planes = planes.size();

    for(size_t i = 0; i < n_planes; i++) {
        for(size_t j = i + 1; j < n_planes; j++) {
            for(size_t k = j + 1; k < n_planes; k++) {
                const Plane &p1 = planes[i];
                const Plane &p2 = planes[j];
                const Plane &p3 = planes[k];

                // we need to construct a linear system and solve it to find the intersection
                // A*x = b
                // we start from n*(x - p) = 0 => n*x = n*p
                // Let n*p = b => A * x = b
                Matrix3 A;
                A.col1 = glm::vec3(p1.n.x, p2.n.x, p3.n.x);
                A.col2 = glm::vec3(p1.n.y, p2.n.y, p3.n.y);
                A.col3 = glm::vec3(p1.n.z, p2.n.z, p3.n.z);

                float det = determinant(A);
                if(std::abs(det) < 0.0001f) {
                    continue;
                }
                
                glm::vec3 b(dot3(p1.n, p1.p), dot3(p2.n, p2.p), dot3(p3.n, p3.p));

                Matrix3 inverseA = inverse(A);

                glm::vec3 part1 = mul(inverseA.col1, b.x);
                glm::vec3 part2 = mul(inverseA.col2, b.y);
                glm::vec3 part3 = mul(inverseA.col3, b.z);
                glm::vec3 resultPoint = part1 + part2 + part3;

                bool isInside = true;
                for(const auto &plane : planes) {
                    // n * (x - p) = 0
                    float dist = dot3(plane.n, resultPoint - plane.p);

                    if(std::abs(dist) < 0.0001f) {
                        isInside = true;
                        break;
                    }
                }
                
                if(isInside) {
                    bool duplicate = false;

                    for (const auto &existing : result) {
                        if(length(existing - resultPoint) < 0.0001f) {
                            duplicate = true;
                            break;
                        }
                    }

                    if(!duplicate) {
                        result.push_back(resultPoint);
                    }
                }
            }
        }
    }

    return result;
}

// ==================================
// ========    Exercise 3    ========
// ==================================
// Given a (possibly irregular) polygon by its vertices (at least 3, ordered either clock- or counterclockwise) that all lie on some 2D plane, 
// compute its area by splitting it into triangles and summing up their areas.
// You may assume that the n-Gon is convex.

float areaOfTriangle(const std::array<glm::vec3, 3> triangle)
{
    float result;

    glm::vec3 A = triangle.at(0);
    glm::vec3 B = triangle.at(1);
    glm::vec3 C = triangle.at(2);

    glm::vec3 edge1 = A - B;
    glm::vec3 edge2 = C - B;

    glm::vec3 crossResult = cross3(edge1, edge2);

    return 0.5f * length(crossResult);
}

std::vector<std::array<glm::vec3, 3>> splitNGonIntoTriangles(const std::vector<glm::vec3> nGon)
{
    std::vector<std::array<glm::vec3, 3>> result;
    size_t count = nGon.size();

    if (nGon.size() < 3) {
        return result;
    }

    glm::vec3 center(0.0f);

    for(const auto &v : nGon) {
        center += v;
    }
    center /= static_cast<float>(count);

    for(size_t i = 0; i < count; i++) {
        std::array<glm::vec3, 3> triangle;

        size_t nextIndex = (i + 1) % count;
        triangle[0] = center;
        triangle[1] = nGon[i];
        triangle[2] = nGon[nextIndex];

        result.push_back(triangle);
    }

    return result;
}

float areaOfIrregularNGon(const std::vector<glm::vec3> nGon)
{
    // Split the nGon into triangles
    std::vector<std::array<glm::vec3, 3>> triangles = splitNGonIntoTriangles(nGon);
    // Sum up their areas
    float area = 0.;
    for (const auto& triangle : triangles) {
        area += areaOfTriangle(triangle);
    }
    return area;
}
