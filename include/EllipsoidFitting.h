#ifndef _LEAST_SQUARES_H_
#define _LEAST_SQUARES_H_

#include "Matrix.h"

class EllipsoidFitting{
    private:
        Matrix<float> P;
        Matrix<float> Theta;
        Matrix<float> r;
        void input(float x, float y, float z);

    public:
        //better to provide Theta
        EllipsoidFitting():P(6,6,1.0f),Theta(6,1,1.0f),r(1,6){};
        EllipsoidFitting(float* theta):P(6,6,1.0f),Theta(6,1,(float*)theta),r(1,6){};
        EllipsoidFitting(float* theta, float* p):P(6,6,(float*)p),Theta(6,1,(float*)theta),r(1,6){};
        void update(float x, float y, float z);
        void getAffineTrans(float& x_c, float& y_c, float& z_c, float& x_r, float& y_r,float& z_r);
};

#endif