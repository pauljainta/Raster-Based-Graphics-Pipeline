#include<stdio.h>
#include<iostream>
#include<stack>
#include<math.h>
using namespace std;

struct point
{
    double x,y,z;
};

struct vector
{
    double x,y,z;
    vector(){

    }
    vector(double xx,double yy,double zz)
    {
        x=xx;
        y=yy;
        z=zz;
    }
};


struct matrix
{
    double arr[4][4];
};

stack<matrix> matrixStack;


struct vector normalizeVector(struct vector v)
{
    vector returnVec;
    double vectorValue=sqrt(pow(v.x,2)+pow(v.y,2)+pow(v.z,2));
    returnVec.x=v.x/vectorValue;
    returnVec.y=v.y/vectorValue;
    returnVec.z=v.z/vectorValue;
    return returnVec;




}



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


struct vector rotateVector(struct vector  axis,struct vector v,double rotation_angle)
{
    vector returnVec;
    double dotProduct=dotMultiply(axis,v);
    vector crossProduct=crossMultiply(v,axis);
    returnVec.x=axis.x*cos(rotation_angle)+crossProduct.x*sin(rotation_angle)+(v.x)*dotProduct*(1-cos(rotation_angle));
    returnVec.y=axis.y*cos(rotation_angle)+crossProduct.y*sin(rotation_angle)+(v.y)*dotProduct*(1-cos(rotation_angle));
    returnVec.z=axis.z*cos(rotation_angle)+crossProduct.z*sin(rotation_angle)+(v.z)*dotProduct*(1-cos(rotation_angle));
    return returnVec; 

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
    vector a(ax,ay,az);
    vector normlized_a=normalizeVector(a);
    matrix R=generateIdentityMatrix();
    vector i_vector(1,0,0),j_vector(0,1,0),k_vector(0,0,1);
    vector c1=rotateVector(i_vector,normlized_a,rotation_angle);
    vector c2=rotateVector(j_vector,normlized_a,rotation_angle);
    vector c3=rotateVector(k_vector,normlized_a,rotation_angle);
    R.arr[0][0]=c1.x;
    R.arr[1][0]=c1.y;
    R.arr[2][0]=c1.z;

    R.arr[0][1]=c2.x;
    R.arr[1][1]=c2.y;
    R.arr[2][1]=c2.z;

    R.arr[0][2]=c3.x;
    R.arr[1][2]=c3.y;
    R.arr[2][2]=c3.z;

    return R;



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

    //testing Rotation matirx
     matrix R=generateRotationMatrix(90,1,2.3,4);

     for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
           printf("%lf ",R.arr[i][j]);
        }
        printf("\n");

    }






    return 0;
}

