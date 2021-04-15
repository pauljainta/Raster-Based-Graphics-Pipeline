#include<stdio.h>
#include<iostream>
#include<stack>
using namespace std;

struct point
{
    double x,y,z;
};

struct vector
{
    double x,y,z;
};


struct matrix
{
    double arr[4][4];
};

stack<matrix> matrixStack;


struct vector crossMultiply(struct vector v1,struct vector v2)
{
    struct vector crossProductVector;

    crossProductVector.x = v1.y * v2.z - v1.z * v2.y;
    crossProductVector.y = v1.z * v2.x - v1.x * v2.z;
    crossProductVector.z = v1.x * v2.y - v1.y * v2.x;

    return crossProductVector; 


}


double dotMultiply(struct vector v1, struct vector v2)
{

    double product = 0;
    product = product + (v1.x*v1.x)+(v1.y*v2.y)+(v1.z*v2.z);
    return product;
}


struct vector rotateVector()
{


}


matrix generateIdentityMatrix()
{
    matrix I;
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            if(i==j)
            {

                I.arr[i][j]=1;
            }
            else
            {
                I.arr[i][j]=0;
            }
        }

    }
    return I;


}

matrix generateTranslationMatrix(double tx,double ty,double tz)
{
    matrix T=generateIdentityMatrix();
    T.arr[0][3]=tx;
    T.arr[1][3]=ty;
    T.arr[2][3]=tz;

    return T;



}
matrix generateScalingMatrix(double sx,double sy,double sz)
{
    matrix S=generateIdentityMatrix();
    S.arr[0][0]=sx;
    S.arr[1][1]=sy;
    S.arr[2][2]=sz;

    return S;



}

matrix generateRotationMatrix(double rotation_angle,double ax,double ay,double az)
{
    matrix S=generateIdentityMatrix();
    // S.arr[0][0]=sx;
    // S.arr[1][1]=sy;
    // S.arr[2][2]=sz;

    return S;



}




int main()
{
    matrix I=generateIdentityMatrix();
    matrixStack.push(I);

    //testing translation matrix

    matrix T=generateTranslationMatrix(1,2.3,4);

     for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
           printf("%lf ",T.arr[i][j]);
        }
        printf("\n");

    }

    printf("\n\n");

      //testing scaling matrix

    matrix S=generateScalingMatrix(1,2.3,4);

     for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
           printf("%lf ",S.arr[i][j]);
        }
        printf("\n");

    }






    return 0;
}

