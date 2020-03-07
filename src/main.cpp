#include <Arduino.h>
#include <math.h>
#include <Wire.h>
#include "Matrix.h"
#include "EllipsoidFitting.h"

#define WRITE 0x00
#define READ 0x01

//MPU9250
#define MPU9250_AD 0x68
#define FIFO_EN_AD 0x23
#define PWR_MGMT_1_AD 0x6B
#define ACCEL_XOUT_H_AD 0x3B
#define GYRO_XOUT_H_AD 0x43
#define EXT_SENS_DATA_00_AD 0x49
#define ACCEL_CONFIG_1_AD 0x1C
#define ACCEL_CONFIG_2_AD 0x1D
#define GYRO_CONFIG_AD 0x1B
#define CONFIG_AD 0x1A
#define I2C_MST_CTRL_AD 0x24
#define I2C_SLV0_ADDR_AD 0x25
#define I2C_SLV0_REG_AD 0x26
#define I2C_SLV0_CTRL_AD 0x27
#define INT_BYPASS_CONFIG_AD 0x37
#define USER_CTRL_AD 0x6A
#define ACCEL_SENS 8192.0f
#define GYRO_SENS 32.8f

//Magnetometer
#define MAG_AD 0xC
#define WIA_AD 0x00
#define INFO 0x01
#define STATUS_1_AD 0x02
#define HXL_AD 0x03    //X-axis measurement data lower 8bit
#define HXH_AD 0x04    //X-axis measurement data higher 8bit
#define HYL_AD 0x05    //Y-axis measurement data lower 8bit
#define HYH_AD 0x06    //Y-axis measurement data higher 8bit
#define HZL_AD 0x07    //Z-axis measurement data lower 8bit
#define HZH_AD 0x08    //Z-axis measurement data higher 8bit
#define STATUS_2_AD 0x09
#define CNTL1_AD 0x0A   //control 1
#define CNTL2_AD 0x0B   //control 2
#define ASTC_AD 0x0C    //Self-Test Control
#define TS1_AD 0x0D    //test 1
#define TS2_AD 0x0E   //test 2
#define I2CDIS_AD 0x0F    //I2C disable
#define ASAX_AD 0x10    //Magnetic sensor X-axis sensitivity adjustment value
#define ASAY_AD 0x11    //Magnetic sensor Y-axis sensitivity adjustment value
#define ASAZ_AD 0x12    //Magnetic sensor Z-axis sensitivity adjustment value
#define MAGNE_SENS 6.67f
#define SCALE 0.1499f  // 4912/32760 uT/tick

#define Rad2Deg 57.2958f
#define SELF_AD 0x18
#define I2C_FREQ 400

#define gyroOffsetX -0.26552559f
#define gyroOffsetY -0.253429884f
#define gyroOffsetZ -1.423399411f

#define accel_xc 0.020324488769860f
#define accel_yc 0.012416583502088f
#define accel_zc -0.036144032958744f
#define accel_xr 0.998546693908833f
#define accel_yr 1.005679687468392f
#define accel_zr 1.018845762774923f

// #define DEBUG

volatile static float accelX,accelY,accelZ,gyroX,gyroY,gyroZ,asax,asay,asaz;
volatile static int16_t magneX, magneY,magneZ;
static float x_c,y_c,z_c,x_r,y_r,z_r;

TwoWire i2c(1,I2C_FAST_MODE);

//Kalman Filter
uint32_t prev = 0;
float dt = 0.01;  

Matrix<float> accel(3,1);
Matrix<float> gyro(3,1);
Matrix<float> magnet(3,1);

Matrix<float> omega(4,4);
Matrix<float> A(4,4);
Matrix<float> P_pri(4,4);
Matrix<float> P_post(4,4);
Matrix<float> I_4(4,4,'I');
Matrix<float> q_pri(4,1);
Matrix<float> q_post(4,1);
Matrix<float> q_cor(4,1);
Matrix<float> H(3,4);
Matrix<float> K(4,3);
Matrix<float> h(3,1);

float arrayQ[4][4] = {{0.001,0,0,0},{0,0.001,0,0},{0,0,0.001,0},{0,0,0,0.0001}};
Matrix<float> Q(4,4,(float*)arrayQ);

float arrayR[3][3] = {{1,0,0},{0,1,0},{0,0,0.01}};
Matrix<float> R(3,3,(float*)arrayR);

//initial values for recursive Least-Squares
float magnetrometerTheta[6][1] = {{0.8940},{1.0781},{41.4295},{-77.9046},{-8.6748},{185.6022}};
float magnetrometerP[6][6] = {{0.000000001061542f,0.000000000544315f,-0.000000004173639f,-0.000000090257378f,-0.000000004275027f,0.000000719572396f},
                                {0.000000000544315f,0.000000002204715f,-0.000000003104380f,-0.000000050155517f,-0.000000017238803f,-0.000000293813843f},
                                {-0.000000004173639f,-0.000000003104380f,0.000000555330195f,0.000000435892319f,0.000000061887221f,0.000008429485196f},
                                {-0.000000090257378f,-0.000000050155517f,0.000000435892319f,0.000008178476387f,0.000000336868927f,-0.000076492015500f},
                                {-0.000000004275027f,-0.000000017238803f,0.000000061887221f,0.000000336868927f,0.000001054258614f,0.000002739961276f},
                                {0.000000719572396f,-0.000000293813843f,0.000008429485196f,-0.000076492015500f,0.000002739961276f,0.002015553457120f}};
//ellipsoid fitting 
EllipsoidFitting magnetrometer((float*)magnetrometerTheta,(float*)magnetrometerP);

//z axis pointing upwards at the beginning


void MPU9250Setup(){
    i2c.beginTransmission(MPU9250_AD);
    i2c.write(PWR_MGMT_1_AD);
    i2c.write(0x01); //set the clock reference to X axis gyroscope to get a better accuracy
    i2c.endTransmission();

    i2c.beginTransmission(MPU9250_AD);
    i2c.write(ACCEL_CONFIG_1_AD);
    i2c.write(0x08); //set the accel scale to 4g
    i2c.endTransmission();

    i2c.beginTransmission(MPU9250_AD);
    i2c.write(ACCEL_CONFIG_2_AD);
    i2c.write(0x03);     //turn on the internal low-pass filter for accel with 44.8Hz bandwidth
    i2c.endTransmission();

    i2c.beginTransmission(MPU9250_AD);
    i2c.write(GYRO_CONFIG_AD);
    i2c.write(0x10); //set the gyro scale to 1000 degree/s and FCHOICE_B on
    i2c.endTransmission();

    // turn on the internal low-pass filter for gyro with 41Hz bandwidth
    i2c.beginTransmission(MPU9250_AD);
    i2c.write(CONFIG_AD);
    i2c.write(0x03);
    i2c.endTransmission();

    /*
        disable the I2C Master I/F module; pins ES_DA and ES_SCL are isolated
        from pins SDA/SDI and SCL/ SCLK.
    */
    i2c.beginTransmission(MPU9250_AD);
    i2c.write(USER_CTRL_AD);
    i2c.write(0x00);
    i2c.endTransmission();

    /*
        When asserted, the i2c_master interface pins(ES_CL and ES_DA) will go
        into bypass mode when the i2c master interface is disabled.
    */
    i2c.beginTransmission(MPU9250_AD);
    i2c.write(INT_BYPASS_CONFIG_AD);
    i2c.write(0x02);
    i2c.endTransmission();

    // setup the Magnetometer to fuse ROM access mode to get the Sensitivity
    // Adjustment values and 16-bit output
    i2c.beginTransmission(MAG_AD);
    i2c.write(CNTL1_AD);
    i2c.write(0x1F);
    i2c.endTransmission();

    //wait for the mode changes
    delay(100);

    //read the Sensitivit Adjustment values
    i2c.beginTransmission(MAG_AD);
    i2c.write(ASAX_AD);
    i2c.endTransmission(false);
    i2c.requestFrom(MAG_AD,3);
    asax = (i2c.read()-128)*0.5/128+1;
    asay = (i2c.read()-128)*0.5/128+1;
    asaz = (i2c.read()-128)*0.5/128+1;

    //reset the Magnetometer to power down mode
    i2c.beginTransmission(MAG_AD);
    i2c.write(CNTL1_AD);
    i2c.write(0x00);
    i2c.endTransmission();

    //wait for the mode changes
    delay(100);

    //set the Magnetometer to continuous mode 2(100Hz) and 16-bit output
    i2c.beginTransmission(MAG_AD);
    i2c.write(CNTL1_AD);
    i2c.write(0x16);
    i2c.endTransmission();

    //wait for the mode changes
    delay(100);
}

void readAccel(){
    //read the accelerate
    i2c.beginTransmission(MPU9250_AD);
    i2c.write(ACCEL_XOUT_H_AD);
    i2c.endTransmission();  
    i2c.requestFrom(MPU9250_AD,6);
    accelX = int16_t((i2c.read()<<8) | i2c.read());
    accelY = int16_t((i2c.read()<<8) | i2c.read());
    accelZ = int16_t((i2c.read()<<8) | i2c.read());
}

void readGyro(){
    //read the gyro
    i2c.beginTransmission(MPU9250_AD);
    i2c.write(GYRO_XOUT_H_AD);
    i2c.endTransmission();   
    i2c.requestFrom(MPU9250_AD,6);
    gyroX = int16_t(i2c.read()<<8) | i2c.read();
    gyroY = int16_t(i2c.read()<<8) | i2c.read();
    gyroZ = int16_t(i2c.read()<<8) | i2c.read();
}

void readMagnetometer(){
    i2c.beginTransmission(MAG_AD);
    i2c.write(STATUS_1_AD);
    i2c.endTransmission();   
    i2c.requestFrom(MAG_AD,1);
    if((i2c.read() & 0x01) == 0x01){
        i2c.beginTransmission(MAG_AD);
        i2c.write(HXL_AD);
        i2c.endTransmission();   
        i2c.requestFrom(MAG_AD,7);

        unsigned char buffer[7] = {0};
        for(int i=0;i<7;i++){
            buffer[i] = i2c.read();
        }

        if(!(buffer[6] & 0x8)){
          magneX = int16_t(buffer[0]) | int16_t(buffer[1]<<8);
          magneY = int16_t(buffer[2]) | int16_t(buffer[3]<<8);
          magneZ = int16_t(buffer[4]) | int16_t(buffer[5]<<8);
        }
    }
}

inline void composeOmega(){
	omega[0][1] = -gyro[0][0];
	omega[0][2] = -gyro[1][0];
	omega[0][3] = -gyro[2][0];

	omega[1][0] = gyro[0][0];
	omega[1][2] = gyro[2][0];
	omega[1][3] = -gyro[1][0];
	
	omega[2][0] = gyro[1][0];
	omega[2][1] = -gyro[2][0];
	omega[2][3] = gyro[0][0];

	omega[3][0] = gyro[2][0];
	omega[3][1] = gyro[1][0];
	omega[3][2] = -gyro[0][0];
}

inline void composeJacobian_1(Matrix<float>& jacb, Matrix<float>& q){
	jacb[0][0] = -2*q[2][0];
	jacb[0][1] = 2*q[3][0];
	jacb[0][2] = -2*q[0][0];
	jacb[0][3] = 2*q[1][0];

	jacb[1][0] = 2*q[1][0];
	jacb[1][1] = 2*q[0][0];
	jacb[1][2] = 2*q[3][0];
	jacb[1][3] = 2*q[2][0];

	jacb[2][0] = 2*q[0][0];
	jacb[2][1] = -2*q[1][0];
	jacb[2][2] = -2*q[2][0];
	jacb[2][3] = 2*q[3][0];
}

inline void composeJacobian_2(Matrix<float>& jacb, Matrix<float>& q){
    jacb[0][0] = 2*q[3][0];
    jacb[0][1] = 2*q[2][0];
    jacb[0][2] = 2*q[1][0];
    jacb[0][3] = 2*q[0][0];

    jacb[1][0] = 2*q[0][0];
    jacb[1][1] = -2*q[1][0];
    jacb[1][2] = -2*q[2][0];
    jacb[1][3] = -2*q[3][0];

    jacb[2][0] = -2*q[1][0];
    jacb[2][1] = -2*q[0][0];
    jacb[2][2] = 2*q[3][0];
    jacb[2][3] = 2*q[2][0];
}

//The estimated gravity
inline void composeGravity(Matrix<float>& q){
    h[0][0] = 2*q[1][0]*q[3][0] - 2*q[0][0]*q[2][0];
    h[1][0] = 2*q[0][0]*q[1][0] + 2*q[2][0]*q[3][0];
    h[2][0] = q[0][0]*q[0][0] - q[1][0]*q[1][0] - q[2][0]*q[2][0]+q[3][0]*q[3][0];
}

//The estimated magnetic field
inline void composeMagnet(Matrix<float>& q){    
    h[0][0] = 2*q[1][0]*q[2][0] + 2*q[0][0]*q[3][0];
    h[1][0] = q[0][0]*q[0][0] - q[1][0]*q[1][0] - q[2][0]*q[2][0] - q[3][0]*q[3][0];
    h[2][0] = 2*q[2][0]*q[3][0] - 2*q[0][0]*q[1][0];
}

void setup() {
    i2c.begin();
    MPU9250Setup();
    delay(1000);

    //The initial q should be (1,0,0,0)
    //means no rotation
    q_post[0][0] = 1;
}

void loop(){
    uint32_t now = micros();
    if(now - prev >= 10000){    //it takes 5ms to complete on a 72MHz processor
        prev = micros();

        readAccel();
        readGyro();
        readMagnetometer();

        //calibrate acceleration
        accel[0][0] = (accelX/ACCEL_SENS - accel_xc)/accel_xr;
        accel[1][0] = (accelY/ACCEL_SENS - accel_yc)/accel_yr;
        accel[2][0] = (accelZ/ACCEL_SENS - accel_zc)/accel_zr;

        //calibrate gyro
        gyro[0][0] = gyroX/GYRO_SENS - gyroOffsetX;
        gyro[1][0] = gyroY/GYRO_SENS - gyroOffsetY;
        gyro[2][0] = gyroZ/GYRO_SENS - gyroOffsetZ;

        //convert unit
        gyro = gyro/Rad2Deg;

        //Recursive Least-Squares for online identification
        magnetrometer.update(magneX*asax*SCALE,magneY*asay*SCALE,magneZ*asaz*SCALE);
        magnetrometer.getAffineTrans(x_c,y_c,z_c,x_r,y_r,z_r);

        magnet[0][0] = (magneY*asay*SCALE-y_c)/y_r;
        magnet[1][0] = (magneX*asax*SCALE-x_c)/x_r;
        magnet[2][0] = -(magneZ*asaz*SCALE-z_c)/z_r;

        
		composeOmega();
		A = I_4 + 0.5f*omega*dt;
		q_pri = A*q_post;
		P_pri = A*P_post*Matrix<float>::transpose(A) + Q;

        //stage 1
		composeJacobian_1(H, q_pri);
        K = P_pri*Matrix<float>::transpose(H)*(H*P_pri*Matrix<float>::transpose(H) + R);
        composeGravity(q_pri);
        q_cor = K*(accel - h);
        q_cor[3][0] = 0;
        q_post = q_cor + q_pri;
        P_post = (I_4 - K*H)*P_pri;

        //stage 2
        composeJacobian_2(H,q_pri);
        K = P_pri*Matrix<float>::transpose(H)*(H*P_pri*Matrix<float>::transpose(H)+R);
        composeMagnet(q_pri);
        q_cor = K*(magnet - h);
        q_cor[1][0] = 0;
        q_cor[2][0] = 0;
        q_post = q_cor + q_post;
        P_post = (I_4 - K*H)*P_post;

#ifdef DEBUG  //use Serial0, through st-link
        Matrix<float> X(3,1);
        X[0][0] = atan2f(2*(q_post[2][0]*q_post[3][0]+q_post[0][0]*q_post[1][0]), q_post[0][0]*q_post[0][0]-q_post[1][0]*q_post[1][0]-q_post[2][0]*q_post[2][0]+q_post[3][0]*q_post[3][0]);
        X[1][0] = asinf(2*(q_post[0][0]*q_post[2][0]-q_post[1][0]*q_post[3][0])/(q_post[0][0]*q_post[0][0]+q_post[1][0]*q_post[1][0]+q_post[2][0]*q_post[2][0]+q_post[3][0]*q_post[3][0]));
        X[2][0] = atan2f(2*(q_post[1][0]*q_post[2][0]+q_post[0][0]*q_post[3][0]), q_post[0][0]*q_post[0][0]+q_post[1][0]*q_post[1][0]-q_post[2][0]*q_post[2][0]-q_post[3][0]*q_post[3][0]);

        X = X*Rad2Deg;

        Serial.print(X[0][0]);
        Serial.print(" ");
        Serial.print(X[1][0]);
        Serial.print(" ");
        Serial.println(X[2][0]);
#else  
        //for unity3D
        Serial.write(85);
        Serial.write( (byte*)&q_post[0][0], 4);
        Serial.write( (byte*)&q_post[1][0], 4);
        Serial.write( (byte*)&q_post[2][0], 4);
        Serial.write( (byte*)&q_post[3][0], 4);
        Serial.write(86);        
#endif
    }
}

