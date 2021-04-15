#include<stdio.h>
#include<iostream>
#include<stack>
#include<math.h>
#include<fstream>
#include<string>
#include<cstdlib>
using namespace std;

#define pi (2*acos(0.0))

ifstream in("scene.txt");
ofstream stage1Output("stage1.txt");


void skipNFileLines(int N)
{

    string s;
    s.reserve(30);

   //skip N lines
    for(int i = 0; i < N; ++i)
        getline(in, s);
   

}


struct point
{
    double x,y,z;
    point(){

    }
    point(double xx,double yy,double zz)
    {
        x=xx;
        y=yy;
        z=zz;
    }


};

struct point_matrix
{
    double arr[4][1];
     point_matrix(){

    }
    point_matrix(struct point p)
    {
      arr[0][0]=p.x;
      arr[1][0]=p.y;
      arr[2][0]=p.z;
      arr[3][0]=1.0;
    }


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

void printPoint(struct point p)
{
    cout<<"POINT"<<endl;
    cout<<p.x<<" ";
    cout<<p.y<<" ";
    cout<<p.z<<" ";

    cout<<endl;

}

void printMatirx(struct matrix m)
{
    cout<<"Matrix"<<endl;

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            cout<<m.arr[i][j]<<" ";
        }
        cout<<endl;
    }

}

void printPointMatirx(struct point_matrix m)
{
    cout<<"POINTMatrix"<<endl;
    

    for(int i=0;i<4;i++)
    {
       // //for(int j=0;j<4;j++)
      //  {
            cout<<m.arr[i][0]<<" ";
       // }
       // cout<<endl;
    }
    cout<<endl;


}

void printVector( struct vector v)
{
    cout<<"Vector"<<endl;
    cout<<v.x<<" ";
    cout<<v.y<<" ";
    cout<<v.z<<" ";

    cout<<endl;

}

struct point pointFromPointMatrix(struct point_matrix p)
{
    struct point returnPoint;
    returnPoint.x=p.arr[0][0];
    returnPoint.y=p.arr[1][0];
    returnPoint.z=p.arr[2][0];
    return returnPoint;

}



stack<matrix> matrixStack;
stack<int> stackSizeStoreStack;


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
    product = product + (v1.x*v2.x)+(v1.y*v2.y)+(v1.z*v2.z);
    return product;
}


struct vector rotateVector(struct vector  axis,struct vector v,double rotation_angle)
{
    vector returnVec;
    double dotProduct=dotMultiply(axis,v);
   // cout<<"dot "<<dotProduct<<endl;
    vector crossProduct=crossMultiply(v,axis);
 //   cout<<"cross P"<<endl;
   // printVector(crossProduct);

  //  cout<<"sxis vector"<<endl;
  //  printVector(axis);

   // cout<<"a vector normalized"<<endl;
  //  printVector(v);



   // cout<<"angle "<<cos(rotation_angle)<<endl;

    double s_theta=sin(rotation_angle);
    double c_theta=sqrt(1-pow(s_theta,2));
    double one_minus_c_theta=1-c_theta;

    returnVec.x=axis.x*c_theta+crossProduct.x*s_theta+(dotProduct)*v.x*one_minus_c_theta;
    returnVec.y=axis.y*c_theta+crossProduct.y*s_theta+(dotProduct)*v.y*one_minus_c_theta;
    returnVec.z=axis.z*c_theta+crossProduct.z*s_theta+(dotProduct)*v.z*one_minus_c_theta;
   // cout<<"returnVec Inside RotateVec"<<endl;
  //  printVector(returnVec);
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
    rotation_angle = rotation_angle * (pi/180);
    vector a(ax,ay,az);
    vector normlized_a=normalizeVector(a);
    matrix R=generateIdentityMatrix();
    vector i_vector(1,0,0),j_vector(0,1,0),k_vector(0,0,1);
    vector c1=rotateVector(i_vector,normlized_a,rotation_angle);
   // cout<<"C1"<<endl;
    //printVector(c1);
    vector c2=rotateVector(j_vector,normlized_a,rotation_angle);
   // cout<<"C2"<<endl;
    //printVector(c2);
    vector c3=rotateVector(k_vector,normlized_a,rotation_angle);
  //  cout<<"C3"<<endl;
    //printVector(c3);
    R.arr[0][0]=c1.x;
    R.arr[0][1]=c2.x;
    R.arr[0][2]=c3.x;

    R.arr[1][0]=c1.y;
    R.arr[1][1]=c2.y;
    R.arr[1][2]=c3.y;

    R.arr[2][0]=c1.z;
    R.arr[2][1]=c2.z;
    R.arr[2][2]=c3.z;

   // cout<<"inside genR"<<endl;

  //  printMatirx(R);

    return R;



}

struct point_matrix pointTranformationMatrixMultiply(struct matrix m1,struct point_matrix m2)
{
   // printMatirx(m1);
    struct point_matrix returnMatrix;
     for(int i=0; i<4; ++i)
     {
         for(int j=0; j<1; ++j)
         {
            returnMatrix.arr[i][j] = 0;
            for(int k=0; k<4; ++k) 
            {
               returnMatrix.arr[i][j]+=m1.arr[i][k]*m2.arr[k][j];
            }
         }
     }

     return returnMatrix;
}


struct matrix matrixMultiply(struct matrix m1,struct  matrix m2)
{
    struct matrix returnMatrix;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            returnMatrix.arr[i][j] = 0;
            for (int k = 0; k < 4; k++)
                returnMatrix.arr[i][j] += m1.arr[i][k] * m2.arr[k][j];
        }
    }
    return returnMatrix;

}


struct point transformPoint(struct point p,struct matrix m)
{
    struct point_matrix p_matrix(p);
   // printPointMatirx(p_matrix);
    struct point_matrix transformedPointMatrix=pointTranformationMatrixMultiply(m,p_matrix);
    //printPointMatirx(transformedPointMatrix);
    struct point transformedPoint=pointFromPointMatrix(transformedPointMatrix);

    //printPoint(transformedPoint);

    return transformedPoint;

}






int main()
{
   
    matrix I=generateIdentityMatrix();
    matrixStack.push(I);
    string command;
    skipNFileLines(4);





    while (true)
   {
        getline(in,command);
        if(command=="triangle")
        {
          //  stage1Output<<"bal1"<<endl;
            string fileLine;

            struct point trianglePoints[3],transformedTrianglePoints[3];
            for(int i=0;i<3;i++)
            {

                getline(in,fileLine);
                size_t pos = 0;
                pos = fileLine.find(" ");
                double x = stod(fileLine.substr(0,pos));
                fileLine.erase(0, pos + 1); 

                pos = fileLine.find(" ");
                double y = stod(fileLine.substr(0,pos));
                fileLine.erase(0, pos + 1); 

                pos = fileLine.find(" ");
                double z = stod(fileLine.substr(0,pos));

                trianglePoints[i].x=x;
                trianglePoints[i].y=y;
                trianglePoints[i].z=z;
            


            }



            for(int i=0;i<3;i++)
            {
               /// struct matrix m=matrixStack.top();
               // printMatirx(m);

                
          
                transformedTrianglePoints[i]=transformPoint(trianglePoints[i],matrixStack.top());
                stage1Output<<transformedTrianglePoints[i].x<<" ";
                stage1Output<<transformedTrianglePoints[i].y<<" ";
                stage1Output<<transformedTrianglePoints[i].z<<" ";
                stage1Output<<endl;
              
                 

            }
            stage1Output<<endl;


        }

        else if(command=="translate")
        {
           // stage1Output<<"bal2"<<endl;
            string fileLine;

            struct point translatePoints[3],transformedTranslatePoints[3];
            
     

            getline(in,fileLine);
            size_t pos = 0;
            pos = fileLine.find(" ");
            double x = stod(fileLine.substr(0,pos));
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double y = stod(fileLine.substr(0,pos));
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double z = stod(fileLine.substr(0,pos));

            struct matrix T=generateTranslationMatrix(x,y,z);
            matrixStack.push(matrixMultiply(matrixStack.top(),T));


        }

        else if(command=="scale")
        {
           // stage1Output<<"bal3"<<endl;
            string fileLine;

            struct point scalePoints[3],transformedScalePoints[3];
    

            getline(in,fileLine);
            size_t pos = 0;
            pos = fileLine.find(" ");
            double x = stod(fileLine.substr(0,pos));
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double y = stod(fileLine.substr(0,pos));
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double z = stod(fileLine.substr(0,pos));
            struct matrix S=generateScalingMatrix(x,y,z);
            matrixStack.push(matrixMultiply(matrixStack.top(),S));


        }


        else if(command=="rotate")
        {
           // stage1Output<<"bal4"<<endl;
            string fileLine;

        
            getline(in,fileLine);
            size_t pos = 0;
            pos = fileLine.find(" ");
            double angle = stod(fileLine.substr(0,pos));
           // cout<<angle<<endl;
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double x = stod(fileLine.substr(0,pos));
           // cout<<x<<endl;
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double y = stod(fileLine.substr(0,pos));
           // cout<<y<<endl;
            fileLine.erase(0, pos + 1); 

            pos = fileLine.find(" ");
            double z = stod(fileLine.substr(0,pos));
          //  cout<<z<<endl;

           // stage1Output<<"bal5"<<endl;

            struct matrix R=generateRotationMatrix(angle,x,y,z);
           // cout<<"Rotation matrix";
           // printMatirx(R);
            matrixStack.push(matrixMultiply(matrixStack.top(),R));



        }

        else if(command=="push")
        {
          //  stage1Output<<"bal push"<<endl;
            stackSizeStoreStack.push(matrixStack.size());   
        }

        else if(command=="pop")
        {
           // stage1Output<<"bal pop"<<endl;
            int popped_length=stackSizeStoreStack.top();
            stackSizeStoreStack.pop();
            while(matrixStack.size()!=popped_length)
            {
                matrixStack.pop();
            }

        }

        else if(command=="end")
        {
            break;
        }



   }


    return 0;
}

