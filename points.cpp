/**
 * EXERCISE
Please write a C++ program to solve the following task:
Given a text file containing the locations of a set of points in 2D space, fit a circle through the points using the Modified Least Squares method as described in the attached paper.
Each line of the input file contains the x and y coordinate of one point separated by a tab character.
 */

// 1.0\t1.0
// 1.3\t2.1

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <filesystem>

/**
 * Structures represents a point in 2D space
 */
class Point
{
public:
    double x;
    double y;
};

/**
 * Read the 2D points from input stream and put them in the corresponding structure
 */
class Reader
{
public:
    static std::vector<Point> readPoints(std::istream& input)
    {
        std::vector<Point> points;

        Point point;
        while( input >> point.x >> point.y)
        {
            points.push_back(point);
        }

        return points;
    }
};

/**
 * This class provides an implementation for Modified Least Square algortihm.
 * Please refer to A Few Methods for Fitting Circles to Data by Dale Umbach, Kerry N. Jones (DOI:10.1109/TIM.2003.820472) for more information.
 */
class MLSSolver
{
    double m_a;
    double m_b;
    double m_r;

public:
    MLSSolver(const std::vector<Point>& points)
    {
        if (points.size() < 3)
        {
            throw std::invalid_argument("Number of points (n) should be n >= 3. Found only " + std::to_string(points.size()) + " points.");
        }

        double n = static_cast<double>(points.size());

        double sum_x = 0.0;
        double sum_x_sq = 0.0;
        double sum_x_cube = 0.0;
        double sum_y = 0.0;
        double sum_y_sq = 0.0;
        double sum_y_cube = 0.0;
        double sum_xy = 0.0;
        double sum_xsq_y = 0.0;
        double sum_x_ysq = 0.0;

        for (const auto &[x, y] : points)
        {
            double x_sq = x * x;
            double y_sq = y * y;

            sum_x += x;
            sum_x_sq += x_sq;
            sum_x_cube += x_sq * x;
            sum_y += y;
            sum_y_sq += y_sq;
            sum_y_cube += y_sq * y;
            sum_xy += x * y;
            sum_xsq_y += x_sq * y;
            sum_x_ysq += x * y_sq;
        }

        double A = n * sum_x_sq - sum_x * sum_x;
        double B = n * sum_xy - sum_x * sum_y;
        double C = n * sum_y_sq - sum_y * sum_y;
        double D = 0.5 * (n * sum_x_ysq - sum_x * sum_y_sq + n * sum_x_cube - sum_x * sum_x_sq);
        double E = 0.5 * (n * sum_xsq_y - sum_y * sum_x_sq + n * sum_y_cube - sum_y * sum_y_sq);

        double denom = ( A * C - B * B );
        m_a = ( D * C - B * E ) / denom;
        m_b = ( A * E - B * D ) / denom;
        
        m_r = 0.0;
        for (const auto &[x, y] : points)
        {
            m_r += std::sqrt( (x - m_a) * (x - m_a) + (y - m_b) * (y - m_b) );
        }
        m_r /= n;
    }

    double getA() const
    {
        return m_a;
    }

    double getB() const
    {
        return m_b;
    }

    double getR() const
    {
        return m_r;
    }
};

/// Test using simple circle (A = B = 0, R = 1)
void runTest1()
{
    std::vector<Point> simpleCircle = {
        {1.0, 0.0},
        {0.0, 1.0},
        {-1.0, 0.0},
        {0.0, -1.0}
    };

    MLSSolver solver(simpleCircle);
    assert(std::abs(solver.getA() - 0.0) < 1e-7);
    assert(std::abs(solver.getB() - 0.0) < 1e-7);
    assert(std::abs(solver.getR() - 1.0) < 1e-7);
}

/// Test using shifted circle (A = B = 1, R = 1)
void runTest2()
{
    std::vector<Point> simpleCircle = {
        {2.0, 1.0},
        {1.0, 2.0},
        {0.0, 1.0},
        {1.0, 0.0}
    };

    MLSSolver solver(simpleCircle);
    assert(std::abs(solver.getA() - 1.0) < 1e-7);
    assert(std::abs(solver.getB() - 1.0) < 1e-7);
    assert(std::abs(solver.getR() - 1.0) < 1e-7);
}

/// Test using shifted circle with bigger radius (A = B = 1, R = 2)
void runTest3()
{
    std::vector<Point> simpleCircle = {
        {3.0, 1.0},
        {1.0, 3.0},
        {-1.0, 1.0},
        {1.0, -1.0}
    };

    MLSSolver solver(simpleCircle);
    assert(std::abs(solver.getA() - 1.0) < 1e-7);
    assert(std::abs(solver.getB() - 1.0) < 1e-7);
    assert(std::abs(solver.getR() - 2.0) < 1e-7);
}

/// Test using simple circle (A = B = 0, R = 1) and different points
void runTest4()
{
    double sin_15_degrees = 0.2588190451;
    double cos_15_degrees = 0.96592582628;

    std::vector<Point> simpleCircle = {
        {-sin_15_degrees, cos_15_degrees},
        {0.0, 1.0},
        {sin_15_degrees, cos_15_degrees}
    };

    MLSSolver solver(simpleCircle);
    assert(std::abs(solver.getA() - 0.0) < 1e-7);
    assert(std::abs(solver.getB() - 0.0) < 1e-7);
    assert(std::abs(solver.getR() - 1.0) < 1e-7);
}

/// Test using points between 2 circles sharing same center and expect radius in between these radiuses
void runTest5()
{
    double r1 = 1.0; // inner circle a = 0, b = 0, r = 1
    double r2 = 1.05; // outer circle a = 0, b = 0, r = 1.05

    double sin_45_deg = 0.70710678118;
    double cos_45_deg = 0.70710678118;

    std::vector<Point> simpleCircle = {
        {r1 * 1.0, r1 * 0.0},
        {r2 * sin_45_deg, r2 * cos_45_deg},
        {r1 * 0.0, r1 * 1.0},
        {r2 * -sin_45_deg, r2 * cos_45_deg},
        {r1 * -1.0, r1 * 0.0},
        {r2 * -sin_45_deg, r2 * -cos_45_deg},
        {r1 * 0.0, r1 * -1.0},
        {r2 * sin_45_deg, r2 * -cos_45_deg}
    };

    MLSSolver solver(simpleCircle);
    assert(solver.getR() > r1);
    assert(solver.getR() < r2);
}

int main(int argc, char** argv)
{
#ifdef RUN_TESTS
    runTest1();
    runTest2();
    runTest3();
    runTest4();
    runTest5();
#endif

    if (argc == 1)
    {
        throw std::invalid_argument("File with input data is not provided");
    }

    const char* fname = argv[1];
    if (!std::filesystem::exists(fname))
    {
        throw std::invalid_argument("File " + std::string(fname) + " does not exist");
    }
    
    std::ifstream input(fname);
    std::vector<Point> points = Reader::readPoints(input);
    MLSSolver solver(points);

    std::cout << "Fitted parameters:" << std::endl;
    std::cout << "\ta = " << solver.getA() << std::endl;
    std::cout << "\tb = " << solver.getB() << std::endl;
    std::cout << "\tr = " << solver.getR() << std::endl;
}