#ifndef HELPERSV3_0_0_H
#define HELPERSV3_0_0_H

/***************************************************************************/
/* This header file contains the V3.0.0 API helper functions and classes,  */
/* which are used to speed up calculations. Include this header into your  */
/* contact model project. These functions and methods will only change     */
/* when version changes.                                                   */
/* All functions and declarations are contained purely within this header. */
/***************************************************************************/

#include <istream>
#include <fstream>
#include <math.h>

namespace NApiHelpersV3_0_0
{

    /***************************************************************************/
    /* Forward declarations of classes in this file                            */
    /***************************************************************************/
    class CSimple3DVector;
    class CSimple3x3Matrix;

    /***************************************************************************/
    /* Constant definitions of use in calculations                             */
    /***************************************************************************/

    const double PI = 3.14159265358979323846;
    const double DOUBLE_UNDERFLOW = 1e-25;
    const double REALLY_REALLY_ZERO = 1e-40;

    /***************************************************************************/
    /* Helper functions for double comparison                                  */
    /***************************************************************************/

    /**
     * Checks if a double is zero.  As doubles produce small variations
     * due to rounding this method returns true if
     * -DOUBLE_UNDERFLOW < val < DOUBLE_UNDERFLOW
     * @param val The value to check
     * @return bool True if the value is zero (within the tolerance)
     */
    template <class T>
    bool isZero(const T val)
    {
        if(val < DOUBLE_UNDERFLOW && val > (-DOUBLE_UNDERFLOW))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * Checks if a double is zero.  As doubles produce small variations
     * due to rounding this method returns true if
     * -tolerance < val < tolerance
     * The tolerance defaults to REALLY_REALLY_ZERO which is substantially more
     * limiting than the value used by isZero()
     * @param val The value to check
     * @param tolerance The tolerance to check to
     * @return bool True if the value is zero (within the tolerance)
     */
    template <class T>
    inline bool isReallyReallyZero(const T val, const T tolerance = REALLY_REALLY_ZERO)
    {
        if(val <= tolerance && val >= (-tolerance))
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * Compares 2 doubles.  As doubles produce small variations due to
     * rounding this method allows the difference to be up to
     * REALLY_REALLY_ZERO (or another user supplied value)
     * @param valA The first value to check
     * @param valB The second value to check
     * @param tolerance The tolerance to check to
     * @return bool True if the values are the same (within the tolerance)
     */
    template <class T>
    inline bool areEqual(const T valA, const T valB, const T tolerance = REALLY_REALLY_ZERO)
    {
        if(!isReallyReallyZero(valB, tolerance))
        {
            // valB is not zero
            return ( fabs(valA/valB - 1.0) < 1e-12 );
        }
        else if(!isReallyReallyZero(valA, tolerance))
        {
            // valA is not zero
            return ( fabs(valB/valA - 1.0) < 1e-12 );
        }
        else
        {
            // Both values are effectivly zero so they are the same
            return true;
        }
    }

    /****************************************************************************/
    /* Class definitions for representation of matrices and vectors in 3D space */
    /****************************************************************************/

    /**
     * The CSimple3DVector class represents a vector in 3D space.
     */
    class CSimple3DVector
    {
    public:
        CSimple3DVector() : m_x(0.0), m_y(0.0), m_z(0.0) { }

        CSimple3DVector(const CSimple3DVector& vect);

        CSimple3DVector(double dx, double dy, double dz) : m_x(dx), m_y(dy), m_z(dz) { }

        virtual ~CSimple3DVector() = default;

        inline CSimple3DVector  operator+(const CSimple3DVector& vector) const;
        inline CSimple3DVector& operator+=(const CSimple3DVector& vector);

        inline CSimple3DVector  operator-(const CSimple3DVector& vector) const;
        inline CSimple3DVector  operator-() const;
        inline CSimple3DVector& operator-=(const CSimple3DVector& vector);

        inline CSimple3DVector  operator*(const double multValue) const;
        inline CSimple3DVector  operator*(const CSimple3x3Matrix& mat) const;
        inline CSimple3DVector& operator*=(const double multValue);
        inline CSimple3DVector& operator*=(const CSimple3x3Matrix& mat);

        inline CSimple3DVector  operator/(const double divValue) const;
        inline CSimple3DVector& operator/=(const double divValue);

        inline bool operator==(const CSimple3DVector& vector) const;
        inline bool operator!=(const CSimple3DVector& vector) const;

        inline CSimple3DVector normalise();
        inline CSimple3DVector cross(const CSimple3DVector& vector) const;
        inline double dot(const CSimple3DVector& vector) const;
        inline double length() const { return sqrt(m_x * m_x + m_y * m_y + m_z * m_z); }
        inline double lengthSquared() const { return m_x * m_x + m_y * m_y + m_z * m_z; }

        inline double getX() const { return m_x; }
        inline double getY() const { return m_y; }
        inline double getZ() const { return m_z; }
        inline void setX(double dx) { m_x = dx; }
        inline void setY(double dy) { m_y = dy; }
        inline void setZ(double dz) { m_z = dz; }

        friend std::istream& operator>>(std::istream& stream, CSimple3DVector& point);

    private:
        double m_x, m_y, m_z;
    };

    /**
     * The CSimple3x3Matrix class represents a 3x3 matrix.
     * Various operations to access and manipulate the matrix are provided
     */
    class CSimple3x3Matrix
    {
    public:
        CSimple3x3Matrix();

        CSimple3x3Matrix(const CSimple3x3Matrix& mat);

        CSimple3x3Matrix(double dxx, double dxy, double dxz,
                         double dyx, double dyy, double dyz,
                         double dzx, double dzy, double dzz);

        CSimple3x3Matrix(const CSimple3DVector& vector, double rotangle);

        virtual ~CSimple3x3Matrix() {}

        inline CSimple3x3Matrix  operator+(const CSimple3x3Matrix& mat) const;
        inline CSimple3x3Matrix& operator+=(const CSimple3x3Matrix& mat);

        inline CSimple3x3Matrix  operator-(const CSimple3x3Matrix& mat) const;
        inline CSimple3x3Matrix& operator-=(const CSimple3x3Matrix& mat);

        inline CSimple3x3Matrix  operator-() const;

        inline CSimple3x3Matrix  operator*(const double multValue) const;
        inline CSimple3x3Matrix  operator*(const CSimple3x3Matrix& mat) const;
        inline CSimple3x3Matrix& operator*=(const double multValue);
        inline CSimple3x3Matrix& operator*=(const CSimple3x3Matrix& mat);

        inline CSimple3x3Matrix  operator/(const double divValue) const;
        inline CSimple3x3Matrix& operator/=(const double divValue);

        inline double            det() const;
        inline CSimple3x3Matrix  inv() const;
        inline CSimple3x3Matrix  transpose() const;

        inline double getXX() const { return m_xx; }
        inline double getXY() const { return m_xy; }
        inline double getXZ() const { return m_xz; }
        inline double getYX() const { return m_yx; }
        inline double getYY() const { return m_yy; }
        inline double getYZ() const { return m_yz; }
        inline double getZX() const { return m_zx; }
        inline double getZY() const { return m_zy; }
        inline double getZZ() const { return m_zz; }
        inline void setXX(double dxx) { m_xx = dxx;}
        inline void setXY(double dxy) { m_xy = dxy;}
        inline void setXZ(double dxz) { m_xz = dxz;}
        inline void setYX(double dyx) { m_yx = dyx;}
        inline void setYY(double dyy) { m_yy = dyy;}
        inline void setYZ(double dyz) { m_yz = dyz;}
        inline void setZX(double dzx) { m_zx = dzx;}
        inline void setZY(double dzy) { m_zy = dzy;}
        inline void setZZ(double dzz) { m_zz = dzz;}

    private:
        double m_xx, m_xy, m_xz, m_yx, m_yy, m_yz, m_zx, m_zy, m_zz;
    };

    /***************************************************************************/
    /* Inline Implementations of all CSimple3DVector methods                   */
    /***************************************************************************/

    inline CSimple3DVector::CSimple3DVector(const CSimple3DVector& vect)
    {
        m_x = vect.m_x;
        m_y = vect.m_y;
        m_z = vect.m_z;
    }

    inline CSimple3DVector CSimple3DVector::normalise()
    {
        double lenSquared = lengthSquared();

        if (!isZero(lenSquared)) //only normalise if length isn't zero.
        {
            double dinv = 1.0 / sqrt(lenSquared);
            m_x *= dinv;
            m_y *= dinv;
            m_z *= dinv;
        }

        return *this;
    }

    inline CSimple3DVector CSimple3DVector::operator+(const CSimple3DVector& vector) const
    {
        CSimple3DVector retVal;

        retVal.setX(vector.m_x + m_x);
        retVal.setY(vector.m_y + m_y);
        retVal.setZ(vector.m_z + m_z);

        return retVal;
    }

    inline CSimple3DVector& CSimple3DVector::operator+=(const CSimple3DVector& vector)
    {
        m_x += vector.m_x;
        m_y += vector.m_y;
        m_z += vector.m_z;

        return *this;
    }

    inline CSimple3DVector CSimple3DVector::operator-(const CSimple3DVector& vector) const
    {
        CSimple3DVector retVal;

        retVal.setX(m_x - vector.m_x);
        retVal.setY(m_y - vector.m_y);
        retVal.setZ(m_z - vector.m_z);

        return retVal;
    }

    inline CSimple3DVector& CSimple3DVector::operator-=(const CSimple3DVector& vector)
    {
        m_x -= vector.m_x;
        m_y -= vector.m_y;
        m_z -= vector.m_z;

        return *this;
    }

    inline CSimple3DVector CSimple3DVector::operator-() const
    {
        return CSimple3DVector(-m_x, -m_y, -m_z);
    }

    inline CSimple3DVector CSimple3DVector::operator*(const double multValue) const
    {
        CSimple3DVector retVal;

        retVal.setX(m_x * multValue);
        retVal.setY(m_y* multValue);
        retVal.setZ(m_z * multValue);

        return retVal;
    }

    inline CSimple3DVector CSimple3DVector::operator*(const CSimple3x3Matrix& mat) const
    {
        CSimple3DVector retVect;

        retVect.setX(m_x * mat.getXX() + m_y * mat.getXY() + m_z * mat.getXZ());
        retVect.setY(m_x * mat.getYX() + m_y * mat.getYY() + m_z * mat.getYZ());
        retVect.setZ(m_x * mat.getZX() + m_y * mat.getZY() + m_z * mat.getZZ());

        return retVect;
    }

    inline CSimple3DVector& CSimple3DVector::operator*=(const double multValue)
    {
        m_x *= multValue;
        m_y *= multValue;
        m_z *= multValue;

        return *this;
    }

    inline CSimple3DVector& CSimple3DVector::operator*=(const CSimple3x3Matrix& mat)
    {
        CSimple3DVector tmp(*this); 
        m_x = tmp.m_x * mat.getXX() + tmp.m_y * mat.getXY() + tmp.m_z * mat.getXZ();
        m_y = tmp.m_x * mat.getYX() + tmp.m_y * mat.getYY() + tmp.m_z * mat.getYZ();
        m_z = tmp.m_x * mat.getZX() + tmp.m_y * mat.getZY() + tmp.m_z * mat.getZZ();

        return *this;
    }

    inline CSimple3DVector CSimple3DVector::operator/(const double divValue) const
    {
        CSimple3DVector retVal;

        if(!isZero(divValue))
        {
            retVal.setX(m_x / divValue);
            retVal.setY(m_y / divValue);
            retVal.setZ(m_z / divValue);
        }
        else
        {
            retVal.setX( 9E30);
            retVal.setY( 9E30);
            retVal.setZ( 9E30);
        }

        return retVal;
    }

    inline CSimple3DVector& CSimple3DVector::operator/=(const double divValue)
    {
        if(!isZero(divValue))
        {
            m_x /= divValue;
            m_y /= divValue;
            m_z /= divValue;
        }
        else
        {
            m_x = 9E30;
            m_y = 9E30;
            m_z = 9E30;
        }

        return *this;
    }

    inline bool CSimple3DVector::operator==(const CSimple3DVector& vector) const
    {
        return (areEqual(m_x, vector.m_x)) &&
               (areEqual(m_y, vector.m_y)) &&
               (areEqual(m_z, vector.m_z));
    }

    inline bool CSimple3DVector::operator!=(const CSimple3DVector& vector) const
    {
        return !(*this == vector);
    }

    inline CSimple3DVector CSimple3DVector::cross(const CSimple3DVector& vector) const
    {
        CSimple3DVector retVect;

        retVect.setX(m_y * vector.m_z - m_z * vector.m_y);
        retVect.setY(m_z * vector.m_x - m_x * vector.m_z);
        retVect.setZ(m_x * vector.m_y - m_y * vector.m_x);

        return retVect;
    }

    inline double CSimple3DVector::dot(const CSimple3DVector& vector) const
    {
        return m_x * vector.m_x + m_y * vector.m_y + m_z * vector.m_z;
    }

    inline std::istream& operator>>(std::istream& stream, CSimple3DVector& vector)
    {
        stream >> vector.m_x >> vector.m_y >> vector.m_z;

        return stream;
    }

    /***************************************************************************/
    /* Inline Implementations of all CSimple3x3Matrix methods                  */
    /***************************************************************************/

    inline CSimple3x3Matrix::CSimple3x3Matrix() :
    m_xx(1.0),
    m_xy(0.0),
    m_xz(0.0),
    m_yx(0.0),
    m_yy(1.0),
    m_yz(0.0),
    m_zx(0.0),
    m_zy(0.0),
    m_zz(1.0)
    {
        ;
    }


    inline CSimple3x3Matrix::CSimple3x3Matrix(double dxx, double dxy, double dxz,
                                              double dyx, double dyy, double dyz,
                                              double dzx, double dzy, double dzz)  :
    m_xx(dxx),
    m_xy(dxy),
    m_xz(dxz),
    m_yx(dyx),
    m_yy(dyy),
    m_yz(dyz),
    m_zx(dzx),
    m_zy(dzy),
    m_zz(dzz)
    {
    }

    inline CSimple3x3Matrix::CSimple3x3Matrix(const CSimple3x3Matrix& mat)
    {
        m_xx = mat.m_xx;
        m_xy = mat.m_xy;
        m_xz = mat.m_xz;
        m_yx = mat.m_yx;
        m_yy = mat.m_yy;
        m_yz = mat.m_yz;
        m_zx = mat.m_zx;
        m_zy = mat.m_zy;
        m_zz = mat.m_zz;
    }

    inline CSimple3x3Matrix::CSimple3x3Matrix(const CSimple3DVector& vector, double rotangle)
    {
        double c = cos(rotangle);
        double s = sin(rotangle);

        m_xx = 1 + ( 1 - c ) * ( vector.getX() * vector.getX() - 1 );
        m_xy = -vector.getZ() * s + ( 1 - c ) * vector.getX() * vector.getY();
        m_xz = vector.getY() * s + ( 1 - c ) * vector.getX() * vector.getZ();
        m_yx = vector.getZ() * s + ( 1 - c ) * vector.getX() * vector.getY();
        m_yy = 1 + ( 1 - c ) * ( vector.getY() * vector.getY() - 1 );
        m_yz = -vector.getX() * s + ( 1 - c ) * vector.getY() * vector.getZ();
        m_zx = -vector.getY() * s + ( 1 - c ) * vector.getX() * vector.getZ();
        m_zy = vector.getX() * s + ( 1 - c ) * vector.getY() * vector.getZ();
        m_zz = 1 + ( 1 - c ) * ( vector.getZ() * vector.getZ() - 1 );

    }

    inline CSimple3x3Matrix CSimple3x3Matrix::operator+(const CSimple3x3Matrix& mat) const
    {
        CSimple3x3Matrix retVal;

        retVal.setXX(mat.m_xx + m_xx);
        retVal.setXY(mat.m_xy + m_xy);
        retVal.setXZ(mat.m_xz + m_xz);
        retVal.setYX(mat.m_yx + m_yx);
        retVal.setYY(mat.m_yy + m_yy);
        retVal.setYZ(mat.m_yz + m_yz);
        retVal.setZX(mat.m_zx + m_zx);
        retVal.setZY(mat.m_zy + m_zy);
        retVal.setZZ(mat.m_zz + m_zz);

        return retVal;
    }

    inline CSimple3x3Matrix& CSimple3x3Matrix::operator+=(const CSimple3x3Matrix& mat)
    {
        m_xx += mat.m_xx;
        m_xy += mat.m_xy;
        m_xz += mat.m_xz;
        m_yx += mat.m_yx;
        m_yy += mat.m_yy;
        m_yz += mat.m_yz;
        m_zx += mat.m_zx;
        m_zy += mat.m_zy;
        m_zz += mat.m_zz;

        return *this;
    }

    inline CSimple3x3Matrix CSimple3x3Matrix::operator-(const CSimple3x3Matrix& mat) const
    {
        CSimple3x3Matrix retVal;

        retVal.setXX(m_xx - mat.m_xx);
        retVal.setXY(m_xy - mat.m_xy);
        retVal.setXZ(m_xz - mat.m_xz);
        retVal.setYX(m_yx - mat.m_yx);
        retVal.setYY(m_yy - mat.m_yy);
        retVal.setYZ(m_yz - mat.m_yz);
        retVal.setZX(m_zx - mat.m_zx);
        retVal.setZY(m_zy - mat.m_zy);
        retVal.setZZ(m_zz - mat.m_zz);

        return retVal;
    }

    inline CSimple3x3Matrix& CSimple3x3Matrix::operator-=(const CSimple3x3Matrix& mat)
    {
        m_xx -= mat.m_xx;
        m_xy -= mat.m_xy;
        m_xz -= mat.m_xz;
        m_yx -= mat.m_yx;
        m_yy -= mat.m_yy;
        m_yz -= mat.m_yz;
        m_zx -= mat.m_zx;
        m_zy -= mat.m_zy;
        m_zz -= mat.m_zz;

        return *this;
    }

    inline CSimple3x3Matrix CSimple3x3Matrix::operator*(const double multValue) const
    {
        CSimple3x3Matrix retVal;

        retVal.setXX(m_xx * multValue);
        retVal.setXY(m_xy * multValue);
        retVal.setXZ(m_xz * multValue);
        retVal.setYX(m_yx * multValue);
        retVal.setYY(m_yy * multValue);
        retVal.setYZ(m_yz * multValue);
        retVal.setZX(m_zx * multValue);
        retVal.setZY(m_zy * multValue);
        retVal.setZZ(m_zz * multValue);

        return retVal;
    }

    inline CSimple3x3Matrix CSimple3x3Matrix::operator*(const CSimple3x3Matrix& mat) const
    {
        CSimple3x3Matrix retVal;

        retVal.setXX(m_xx * mat.m_xx + m_xy * mat.m_yx + m_xz * mat.m_zx);
        retVal.setXY(m_xx * mat.m_xy + m_xy * mat.m_yy + m_xz * mat.m_zy);
        retVal.setXZ(m_xx * mat.m_xz + m_xy * mat.m_yz + m_xz * mat.m_zz);
        retVal.setYX(m_yx * mat.m_xx + m_yy * mat.m_yx + m_yz * mat.m_zx);
        retVal.setYY(m_yx * mat.m_xy + m_yy * mat.m_yy + m_yz * mat.m_zy);
        retVal.setYZ(m_yx * mat.m_xz + m_yy * mat.m_yz + m_yz * mat.m_zz);
        retVal.setZX(m_zx * mat.m_xx + m_zy * mat.m_yx + m_zz * mat.m_zx);
        retVal.setZY(m_zx * mat.m_xy + m_zy * mat.m_yy + m_zz * mat.m_zy);
        retVal.setZZ(m_zx * mat.m_xz + m_zy * mat.m_yz + m_zz * mat.m_zz);

        return retVal;
    }

    inline CSimple3x3Matrix& CSimple3x3Matrix::operator*=(const double multValue)
    {
        m_xx *= multValue;
        m_xy *= multValue;
        m_xz *= multValue;
        m_yx *= multValue;
        m_yy *= multValue;
        m_yz *= multValue;
        m_zx *= multValue;
        m_zy *= multValue;
        m_zz *= multValue;

        return *this;
    }

    inline CSimple3x3Matrix& CSimple3x3Matrix::operator*=(const CSimple3x3Matrix& mat)
    {
        CSimple3x3Matrix tmp(*this);
        m_xx = tmp.m_xx * mat.m_xx + tmp.m_xy * mat.m_yx + tmp.m_xz * mat.m_zx;
        m_xy = tmp.m_xx * mat.m_xy + tmp.m_xy * mat.m_yy + tmp.m_xz * mat.m_zy;
        m_xz = tmp.m_xx * mat.m_xz + tmp.m_xy * mat.m_yz + tmp.m_xz * mat.m_zz;
        m_yx = tmp.m_yx * mat.m_xx + tmp.m_yy * mat.m_yx + tmp.m_yz * mat.m_zx;
        m_yy = tmp.m_yx * mat.m_xy + tmp.m_yy * mat.m_yy + tmp.m_yz * mat.m_zy;
        m_yz = tmp.m_yx * mat.m_xz + tmp.m_yy * mat.m_yz + tmp.m_yz * mat.m_zz;
        m_zx = tmp.m_zx * mat.m_xx + tmp.m_zy * mat.m_yx + tmp.m_zz * mat.m_zx;
        m_zy = tmp.m_zx * mat.m_xy + tmp.m_zy * mat.m_yy + tmp.m_zz * mat.m_zy;
        m_zz = tmp.m_zx * mat.m_xz + tmp.m_zy * mat.m_yz + tmp.m_zz * mat.m_zz;

        return *this;
    }

    inline CSimple3x3Matrix CSimple3x3Matrix::operator/(const double divValue) const
    {
        CSimple3x3Matrix retVal;

        if(!isZero(divValue))
        {
            retVal.setXX(m_xx / divValue);
            retVal.setXY(m_xy / divValue);
            retVal.setXZ(m_xz / divValue);
            retVal.setYX(m_yx / divValue);
            retVal.setYY(m_yy / divValue);
            retVal.setYZ(m_yz / divValue);
            retVal.setZX(m_zx / divValue);
            retVal.setZY(m_zy / divValue);
            retVal.setZZ(m_zz / divValue);
        }
        else
        {
            retVal.setXX(9E30);
            retVal.setXY(9E30);
            retVal.setXZ(9E30);
            retVal.setYX(9E30);
            retVal.setYY(9E30);
            retVal.setYZ(9E30);
            retVal.setZX(9E30);
            retVal.setZY(9E30);
            retVal.setZZ(9E30);
        }

        return retVal;
    }

    inline CSimple3x3Matrix& CSimple3x3Matrix::operator/=(const double divValue)
    {
        if(!isZero(divValue))
        {
            m_xx /= divValue;
            m_xy /= divValue;
            m_xz /= divValue;
            m_yx /= divValue;
            m_yy /= divValue;
            m_yz /= divValue;
            m_zx /= divValue;
            m_zy /= divValue;
            m_zz /= divValue;

        }
        else
        {
            m_xx = 9E30;
            m_xy = 9E30;
            m_xz = 9E30;
            m_yx = 9E30;
            m_yy = 9E30;
            m_yz = 9E30;
            m_zx = 9E30;
            m_zy = 9E30;
            m_zz = 9E30;
        }

        return *this;
    }

    inline double CSimple3x3Matrix::det() const
    {
        // Our matrix is organised as:
        //
        // XX XY XZ
        // YX YY YZ
        // ZX ZY ZZ
        //
        // The formula for The determinant of a 3x3 matrix of this form is:
        //
        // XX*YY*ZZ + XY*YZ*ZX + XZ*YX*ZY - XX*YZ*ZY - XY*YX*ZZ - XZ*YY*ZX
        //
        // from:
        //
        // http://en.wikipedia.org/wiki/Determinant#3-by-3_matrices

        return m_xx * m_yy * m_zz +
               m_xy * m_yz * m_zx +
               m_xz * m_yx * m_zy -
               m_xx * m_yz * m_zy -
               m_xy * m_yx * m_zz -
               m_xz * m_yy * m_zx;
    }

    inline CSimple3x3Matrix CSimple3x3Matrix::inv() const
    {
        CSimple3x3Matrix retVal;

        // Our matrix is organised as:
        //
        // XX XY XZ
        // YX YY YZ
        // ZX ZY ZZ
        //
        // The formula for the inverse of a matrix using Cramer's rule is:
        //
        //     Inverse(A) = 1/Determinant(A) * Adjugate(A)
        //
        // Where
        //
        //     Adjugate(A) is The Transpose of The Cofactor matrix of A
        //
        // from:
        //
        // http://en.wikipedia.org/wiki/Cramer's_rule#Finding_inverse_matrix

        // First calculate the determinant
        double ndet = det();

        // If we have a determinant of 0 we can't get the inverse
        // Set all values of the return matrix to be invalid values
        if (isZero(ndet))
        {
            retVal.setXX(9E30);
            retVal.setXY(9E30);
            retVal.setXZ(9E30);
            retVal.setYX(9E30);
            retVal.setYY(9E30);
            retVal.setYZ(9E30);
            retVal.setZX(9E30);
            retVal.setZY(9E30);
            retVal.setZZ(9E30);
            return retVal;
        }

        // Calculate the inverse
        retVal.setXX((m_yy*m_zz - m_yz * m_zy) / ndet);
        retVal.setXY((m_xz*m_zy - m_xy * m_zz) / ndet);
        retVal.setXZ((m_xy*m_yz - m_xz * m_yy) / ndet);
        retVal.setYX((m_yz*m_zx - m_yx * m_zz) / ndet);
        retVal.setYY((m_xx*m_zz - m_xz * m_zx) / ndet);
        retVal.setYZ((m_xz*m_yx - m_xx * m_yz) / ndet);
        retVal.setZX((m_yx*m_zy - m_yy * m_zx) / ndet);
        retVal.setZY((m_xy*m_zx - m_xx * m_zy) / ndet);
        retVal.setZZ((m_xx*m_yy - m_xy * m_yx) / ndet);

        return retVal;
    }

    inline CSimple3x3Matrix CSimple3x3Matrix::transpose() const
    {
        CSimple3x3Matrix retVal;

        // Our matrix is organised as:
        // XX XY XZ
        // YX YY YZ
        // ZX ZY ZZ

        //It needs to be returned as:
        // XX YX ZX
        // XY YY ZY
        // XZ YZ ZZ

        retVal.setXX(m_xx);
        retVal.setXY(m_yx);
        retVal.setXZ(m_zx);
        retVal.setYX(m_xy);
        retVal.setYY(m_yy);
        retVal.setYZ(m_zy);
        retVal.setZX(m_xz);
        retVal.setZY(m_yz);
        retVal.setZZ(m_zz);

        return retVal;
    }

    /** When retrieving triangle nodes, vertex ids will be specified per triangle */
    struct SGeomTriangleNode
    {
        SGeomTriangleNode(unsigned int a, unsigned int b, unsigned int c) :vertId0(a), vertId1(b), vertId2(c) {}
        unsigned int vertId0;
        unsigned int vertId1;
        unsigned int vertId2;
    };

}

#endif // HELPERSV3_0_0_H
