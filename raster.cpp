#include<iostream>
//#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <bits/stdc++.h>
#include<string.h>
#include<list>

#define pi (2*acos(0.0))

using namespace std;

class point
{
public:

    double x;
    double y;
    double z;
    double w;
};

struct point_vector
{
    double x,y,z;
};

class triangle
{
public:

    point first_point;
    point second_point;
    point third_point;
};

class matrix
{
public:

    double square[4][4];
};

class simply_push
{
public:

    int stack_size_when_I_was_done;
};

list<triangle> triangles;
list<triangle> triangles_viewed;
list<triangle> triangles_projected;
stack<matrix> stack_for_matrices;
stack<simply_push> stack_for_push;

struct point_vector i_unit;
struct point_vector j_unit;
struct point_vector k_unit;

double first_13_value[13];

void initialize_axes()
{
    i_unit.x = 1;
    i_unit.y = 0;
    i_unit.z = 0;

    j_unit.x = 0;
    j_unit.y = 1;
    j_unit.z = 0;

    k_unit.x = 0;
    k_unit.y = 0;
    k_unit.z = 1;
}

double dotProduct(double vect_A[], double vect_B[])
{
    double product = 0;
    for (int i = 0; i < 3; i++)
        product = product + vect_A[i] * vect_B[i];
    return product;
}

point_vector crossProduct(double vect_A[] , double vect_B[])
{
    point_vector cross_P;

    cross_P.x = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
    cross_P.y = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
    cross_P.z = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];

    return cross_P;
}

point_vector rodrigues_rotation(point_vector vect , point_vector axes , double rotation_angle)
{
    //cout<<rotation_angle<<endl;

    point_vector after_rotation;

    double sin_value = sin(rotation_angle);
    //double cos_value = cos(rotation_angle);
    double cos_value = sqrt(1 - (sin_value*sin_value));
    double one_minus_cos = 1 - cos_value;

    //cout<<pi/2<<endl;
    //cout<<cos_value;

    double vect_v[] = {vect.x , vect.y , vect.z};
    double vect_k[] = {axes.x , axes.y , axes.z};

    point_vector k_cross_v = crossProduct(vect_k , vect_v);

    double k_dot_v = dotProduct(vect_k , vect_v);

    after_rotation.x = vect.x*cos_value + k_cross_v.x*sin_value + axes.x*k_dot_v*one_minus_cos;
    after_rotation.y = vect.y*cos_value + k_cross_v.y*sin_value + axes.y*k_dot_v*one_minus_cos;
    after_rotation.z = vect.z*cos_value + k_cross_v.z*sin_value + axes.z*k_dot_v*one_minus_cos;

    return after_rotation;

}

matrix matrix_matrix_multiplication(matrix first , matrix second)
{
    matrix after_multiplication;
    int i,j,k;

    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            after_multiplication.square[i][j] = 0;
        }
    }

    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            for(k=0 ; k<4 ; k++)
            {
                after_multiplication.square[i][j] += first.square[i][k] * second.square[k][j];
            }
        }
    }

    return after_multiplication;
}

point matrix_point_multiplication(matrix mat , point pot)
{
    point after_transformation;
    int i,j,k;

    after_transformation.x = 0;
    after_transformation.y = 0;
    after_transformation.z = 0;
    after_transformation.w = 0;

    double point_container[4] = {pot.x , pot.y , pot.z , pot.w};
    double point_after_transformation[4] = {0 , 0 , 0 , 0};

    for(i=0 ; i<4 ; i++)
    {
        for(k=0 ; k<4 ; k++)
        {
            point_after_transformation[i] += mat.square[i][k] * point_container[k];
        }
    }

    after_transformation.x = point_after_transformation[0];
    after_transformation.y = point_after_transformation[1];
    after_transformation.z = point_after_transformation[2];
    after_transformation.w = point_after_transformation[3];

    return after_transformation;
}

void modeling_transformation()
{
    int i,j;
    double identity_array[4][4] = {{1,0,0,0} , {0,1,0,0} , {0,0,1,0} , {0,0,0,1}};
    matrix identity_matrix;
    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            identity_matrix.square[i][j] = identity_array[i][j];
        }
    }
    stack_for_matrices.push(identity_matrix);

    freopen("scene.txt" , "r" , stdin);
    double skip;
    string command;
    for(i=0 ; i<13 ; i++)
    {
       cin>>skip;
       //cout<<skip<<endl;
       first_13_value[i] = skip;
    }

    while(true)
    {
        cin>>command;

        if(command.compare("triangle") == 0){
            //cout<<"bingo"<<endl;
            triangle a_triangle;
            for(i=0 ; i<3 ; i++)
            {
                point a_point;
                cin>>a_point.x;
                cin>>a_point.y;
                cin>>a_point.z;
                a_point.w = 1;

                matrix current_top = stack_for_matrices.top();
                point after_transformation;
                after_transformation = matrix_point_multiplication(current_top , a_point);

                if(i == 0)
                {
                    a_triangle.first_point.x = after_transformation.x;
                    a_triangle.first_point.y = after_transformation.y;
                    a_triangle.first_point.z = after_transformation.z;
                    a_triangle.first_point.w = after_transformation.w;
                }
                else if(i == 1)
                {
                    a_triangle.second_point.x = after_transformation.x;
                    a_triangle.second_point.y = after_transformation.y;
                    a_triangle.second_point.z = after_transformation.z;
                    a_triangle.second_point.w = after_transformation.w;
                }
                else
                {
                    a_triangle.third_point.x = after_transformation.x;
                    a_triangle.third_point.y = after_transformation.y;
                    a_triangle.third_point.z = after_transformation.z;
                    a_triangle.third_point.w = after_transformation.w;
                }
            }
            triangles.push_back(a_triangle);
        }
        else if(command.compare("translate") == 0){
            double tx,ty,tz;
            cin>>tx;
            cin>>ty;
            cin>>tz;

            matrix translation_matrix;
            double translation_array[4][4] = {{1,0,0,tx} , {0,1,0,ty} , {0,0,1,tz} , {0,0,0,1}};
            for(i=0 ; i<4 ; i++)
            {
                for(j=0 ; j<4 ; j++)
                {
                    translation_matrix.square[i][j] = translation_array[i][j];
                }
            }

            matrix current_top = stack_for_matrices.top();
            matrix after_translation;
            //after_translation = matrix_matrix_multiplication(translation_matrix , current_top);
            after_translation = matrix_matrix_multiplication(current_top , translation_matrix);

            stack_for_matrices.push(after_translation);
        }
        else if(command.compare("scale") == 0){
            double sx,sy,sz;
            cin>>sx;
            cin>>sy;
            cin>>sz;

            matrix scaler_matrix;
            double scaler_array[4][4] = {{sx,0,0,0} , {0,sy,0,0} , {0,0,sz,0} , {0,0,0,1}};
            for(i=0 ; i<4 ; i++)
            {
                for(j=0 ; j<4 ; j++)
                {
                    scaler_matrix.square[i][j] = scaler_array[i][j];
                }
            }

            matrix current_top = stack_for_matrices.top();
            matrix after_scaling;
            //after_scaling = matrix_matrix_multiplication(scaler_matrix , current_top);
            after_scaling = matrix_matrix_multiplication(current_top , scaler_matrix);

            stack_for_matrices.push(after_scaling);
        }
        else if(command.compare("rotate") == 0){
            double angle,ax,ay,az;
            cin>>angle;
            cin>>ax;
            cin>>ay;
            cin>>az;

            //cout<<ax<<ay<<az<<angle;

            double rotation_angle = angle * (pi/180);
            //double rotation_angle = angle;
            double axes_value = sqrt(ax*ax + ay*ay + az*az);
            point_vector normalized_axes;
            normalized_axes.x = ax / axes_value;
            normalized_axes.y = ay / axes_value;
            normalized_axes.z = az / axes_value;

            //cout<<normalized_axes.x<<normalized_axes.y<<normalized_axes.z;
            //cout<<rotation_angle;

            point_vector c1;
            point_vector c2;
            point_vector c3;
            c1 = rodrigues_rotation(i_unit , normalized_axes , rotation_angle);
            c2 = rodrigues_rotation(j_unit , normalized_axes , rotation_angle);
            c3 = rodrigues_rotation(k_unit , normalized_axes , rotation_angle);

            //cout<<c3.x<<c3.y<<c3.z;

            matrix rotation_matrix;
            double rotation_array[4][4] = {{c1.x,c2.x,c3.x,0} , {c1.y,c2.y,c3.y,0} , {c1.z,c2.z,c3.z,0} , {0,0,0,1}};
            for(i=0 ; i<4 ; i++)
            {
                for(j=0 ; j<4 ; j++)
                {
                    rotation_matrix.square[i][j] = rotation_array[i][j];
                }
            }

            matrix current_top = stack_for_matrices.top();
            matrix after_rotation;
            //after_rotation = matrix_matrix_multiplication(rotation_matrix , current_top);
            after_rotation = matrix_matrix_multiplication(current_top , rotation_matrix);

            stack_for_matrices.push(after_rotation);
        }
        else if(command.compare("push") == 0){
            int stack_size = stack_for_matrices.size();
            simply_push just_a_push;
            just_a_push.stack_size_when_I_was_done = stack_size;
            stack_for_push.push(just_a_push);
        }
        else if(command.compare("pop") == 0){
            if(stack_for_push.size() != 0)
            {
                simply_push last_push = stack_for_push.top();
                int items_to_be_popped = stack_for_matrices.size() - last_push.stack_size_when_I_was_done;
                for(i=0 ; i<items_to_be_popped ; i++)
                {
                    stack_for_matrices.pop();
                }
                stack_for_push.pop();
            }
            else
            {
                cout<<"Why POP ? No PUSH was done..."<<endl;
            }
        }
        else if(command.compare("end") == 0){
            cout<<"done"<<endl;
            break;
        }
    }

    fclose(stdin);
}

void view_transformation()
{
    point_vector eye;
    point_vector look;
    point_vector up;

    //freopen("scene.txt" , "r" , stdin);
    eye.x = first_13_value[0];
    eye.y = first_13_value[1];
    eye.z = first_13_value[2];
    look.x = first_13_value[3];
    look.y = first_13_value[4];
    look.z = first_13_value[5];
    up.x = first_13_value[6];
    up.y = first_13_value[7];
    up.z = first_13_value[8];

    point_vector l;
    point_vector u;
    point_vector r;

    //cout<<look.x<<look.y<<look.z<<"??"<<endl;
    //cout<<eye.x<<eye.y<<eye.z<<"??"<<endl;

    l.x = look.x - eye.x;
    l.y = look.y - eye.y;
    l.z = look.z - eye.z;

    //cout<<l.x<<l.y<<l.z<<"??"<<endl;

    double l_value = sqrt(l.x*l.x + l.y*l.y + l.z*l.z);
    l.x = l.x / l_value;
    l.y = l.y / l_value;
    l.z = l.z / l_value;

    //cout<<l.x<<l.y<<l.z<<"??"<<endl;

    double l_vect[] = {l.x , l.y , l.z};
    double up_vect[] = {up.x , up.y , up.z};
    r = crossProduct(l_vect , up_vect);

    double r_value = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);
    r.x = r.x / r_value;
    r.y = r.y / r_value;
    r.z = r.z / r_value;

    double r_vect[] = {r.x , r.y , r.z};
    u = crossProduct(r_vect , l_vect);

    matrix T;
    matrix R;
    double translation_array[4][4] = {{1,0,0,-eye.x} , {0,1,0,-eye.y} , {0,0,1,-eye.z} , {0,0,0,1}};
    double rotation_array[4][4] = {{r.x,r.y,r.z,0} , {u.x,u.y,u.z,0} , {-l.x,-l.y,-l.z,0} , {0,0,0,1}};
    int i,j;
    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            T.square[i][j] = translation_array[i][j];
            R.square[i][j] = rotation_array[i][j];
        }
    }
    matrix V;
    V = matrix_matrix_multiplication(R , T);

    list <triangle> :: iterator it;
    for(it = triangles.begin() ; it != triangles.end() ; ++it)
    {
        triangle now = *it;
        triangle triangle_viewed;
        triangle_viewed.first_point = matrix_point_multiplication(V , now.first_point);
        triangle_viewed.second_point = matrix_point_multiplication(V , now.second_point);
        triangle_viewed.third_point = matrix_point_multiplication(V , now.third_point);
        triangles_viewed.push_back(triangle_viewed);
    }
}

void projection_transformation()
{
    double fovY = first_13_value[9];
    double aspectRatio = first_13_value[10];
    double near = first_13_value[11];
    double far = first_13_value[12];

    double fovX = fovY * aspectRatio;
    double t = near * tan((fovY/2) * (pi/180));
    double r = near * tan((fovX/2) * (pi/180));

    double not_interested1 = -(far+near)/(far-near);
    double not_interested2 = -(2*far*near)/(far-near);
    double projection_array[4][4] = {{near/r,0,0,0} , {0,near/t,0,0} , {0,0,not_interested1,not_interested2} , {0,0,-1,0}};
    matrix P;
    int i,j;
    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            P.square[i][j] = projection_array[i][j];
        }
    }

    list <triangle> :: iterator it;
    for(it = triangles_viewed.begin() ; it != triangles_viewed.end() ; ++it)
    {
        triangle now = *it;
        triangle triangle_projected;
        triangle_projected.first_point = matrix_point_multiplication(P , now.first_point);
        triangle_projected.second_point = matrix_point_multiplication(P , now.second_point);
        triangle_projected.third_point = matrix_point_multiplication(P , now.third_point);

        triangle_projected.first_point.x = triangle_projected.first_point.x / triangle_projected.first_point.w;
        triangle_projected.first_point.y = triangle_projected.first_point.y / triangle_projected.first_point.w;
        triangle_projected.first_point.z = triangle_projected.first_point.z / triangle_projected.first_point.w;
        triangle_projected.first_point.w = triangle_projected.first_point.w / triangle_projected.first_point.w;

        triangle_projected.second_point.x = triangle_projected.second_point.x / triangle_projected.second_point.w;
        triangle_projected.second_point.y = triangle_projected.second_point.y / triangle_projected.second_point.w;
        triangle_projected.second_point.z = triangle_projected.second_point.z / triangle_projected.second_point.w;
        triangle_projected.second_point.w = triangle_projected.second_point.w / triangle_projected.second_point.w;

        triangle_projected.third_point.x = triangle_projected.third_point.x / triangle_projected.third_point.w;
        triangle_projected.third_point.y = triangle_projected.third_point.y / triangle_projected.third_point.w;
        triangle_projected.third_point.z = triangle_projected.third_point.z / triangle_projected.third_point.w;
        triangle_projected.third_point.w = triangle_projected.third_point.w / triangle_projected.third_point.w;

        triangles_projected.push_back(triangle_projected);
    }
}

void print_triangles()
{
    list <triangle> :: iterator it;
    for(it = triangles.begin() ; it != triangles.end() ; ++it)
    {
        triangle now = *it;
        cout<<now.first_point.x<<" "<<now.first_point.y<<" "<<now.first_point.z<<endl;
        cout<<now.second_point.x<<" "<<now.second_point.y<<" "<<now.second_point.z<<endl;
        cout<<now.third_point.x<<" "<<now.third_point.y<<" "<<now.third_point.z<<endl;

        cout<<endl;
    }
}

void print_triangles2()
{
    list <triangle> :: iterator it;
    for(it = triangles_viewed.begin() ; it != triangles_viewed.end() ; ++it)
    {
        triangle now = *it;
        cout<<now.first_point.x<<" "<<now.first_point.y<<" "<<now.first_point.z<<endl;
        cout<<now.second_point.x<<" "<<now.second_point.y<<" "<<now.second_point.z<<endl;
        cout<<now.third_point.x<<" "<<now.third_point.y<<" "<<now.third_point.z<<endl;

        cout<<endl;
    }
}

void print_triangles3()
{
    list <triangle> :: iterator it;
    for(it = triangles_projected.begin() ; it != triangles_projected.end() ; ++it)
    {
        triangle now = *it;
        cout<<now.first_point.x<<" "<<now.first_point.y<<" "<<now.first_point.z<<endl;
        cout<<now.second_point.x<<" "<<now.second_point.y<<" "<<now.second_point.z<<endl;
        cout<<now.third_point.x<<" "<<now.third_point.y<<" "<<now.third_point.z<<endl;

        cout<<endl;
    }
}

int main()
{
    /*freopen("scene.txt" , "r" , stdin);
    double a;
    string s;
    int i;
    for(i=0 ; i<13 ; i++)
    {
       cin>>a;
    }
    cin>>s;
    if(s.compare("triangle") == 0){
        cout<<"bingo"<<endl;
    }
    int arr[4] = {1,5,7,8};
    cout<<arr[3]<<endl;*/

    /*int i , j;
    matrix a;
    matrix b;
    double arra[4][4] = {{5,7,9,10} , {2,3,3,8} , {8,10,2,3} , {3,3,4,8}};
    double arrb[4][4] = {{3,10,12,18} , {12,1,4,9} , {9,10,12,2} , {3,12,4,10}};
    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            a.square[i][j] = arra[i][j];
            b.square[i][j] = arrb[i][j];
        }
    }*/
    /*matrix c;
    c = matrix_matrix_multiplication(a , b);
    for(i=0 ; i<4 ; i++)
    {
        for(j=0 ; j<4 ; j++)
        {
            cout<<c.square[i][j]<<" ";
        }
        cout<<endl;
    }*/
    /*point p;
    p.x = 4;
    p.y = 2;
    p.z = 3;
    p.w = 1;

    point q;
    q = matrix_point_multiplication(a , p);
    cout<<q.x<<" "<<q.y<<" "<<q.z<<" "<<q.w;*/

    initialize_axes();
    modeling_transformation();
    view_transformation();
    projection_transformation();
    print_triangles3();
    //cout<<stack_for_matrices.top().square[3][2];
}
