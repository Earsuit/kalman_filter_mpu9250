#include <Arduino.h>
#include <Wire.h>
#include <math.h>
#include "EllipsoidFitting.h"

void EllipsoidFitting::input(float x, float y, float z){
    r[0][0] = y*y;
    r[0][1] = z*z;
    r[0][2] = x;
    r[0][3] = y;
    r[0][4] = z;
    r[0][5] = 1;
}

void EllipsoidFitting::update(float x, float y, float z){
    input(x,y,z);

    Matrix<float> L = (P*Matrix<float>::transpose(r))/(1.0f+r*P*Matrix<float>::transpose(r));

    Theta = Theta + L*(-x*x - r*Theta);

    P = P - L*r*P;
}

void EllipsoidFitting::getAffineTrans(float& x_c, float& y_c, float& z_c, float& x_r, float& y_r,float& z_r){
    x_c = Theta[2][0]/-2;
    y_c = Theta[3][0]/(-2*Theta[0][0]);
    z_c = Theta[4][0]/(-2*Theta[1][0]);
    x_r = sqrt(x_c*x_c+Theta[0][0]*y_c*y_c+Theta[1][0]*z_c*z_c-Theta[5][0]);
    y_r = sqrt(x_r*x_r/Theta[0][0]);
    z_r = sqrt(x_r*x_r/Theta[1][0]);
}   