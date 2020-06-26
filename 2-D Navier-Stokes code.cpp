//2-D Navier-Stokes(Cavity Flow)
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

//Input set
const int nx=41;   //mesh point in x coordinate
const int ny=41;   //mesh point in y coordinate
const int nt=700;  //number of timesteps
const int nit=50;  //poisson equation iteration

//Each point per distance
const double dx=2.0/(nx-1);
const double dy=2.0/(ny-1);

//Navier-Stokes parameter
const double rho=1;
const double nu=0.1;
const double dt=0.001;

//Initialization Function
void in(double *u, double *un, double *v, double *vn, double *p, double *pn,double *b, int nx, int ny){
    int i, j;
    for(i = 0; i < nx ; i++){
        for(j = 0; j < ny ; j++){
            u[j*nx+i] = 0.0;
            un[j*nx+i] = 0.0;
            v[j*nx+i] = 0.0;
            vn[j*nx+i] = 0.0;
            p[j*nx+i] = 0.0;
            pn[j*nx+i] = 0.0;
            b[j*nx+i] = 0.0;
        }
    }
    for(i = 0 ; i < nx ;i++){
        u[(nx-1)*nx + i] = 1.0;
        un[(nx-1)*nx + i] = 1.0;
    }

cout<<"Initialization"<<endl;
}

//Build_up_b Function
void build_up_b(double *b, double *u,double *v){
   int i, j;
    for(i = 1; i < nx - 1 ; i++){
        for(j = 1; j < ny - 1 ; j++){
            b[j*nx+i] = (rho * ( 1.0/dt * ((u[j*nx + i+1] - u[j*nx + i-1]) / (2 * dx) + (v[(j+1)*nx+i] - v[(j-1)*nx+i]) / (2 * dy)) - ((u[j*nx+i+1] - u[j*nx+i-1]) / (2*dx)) * ((u[j*nx+i+1] - u[j*nx+i-1]) / (2*dx)) - 2 * ((u[(j+1)*nx+i] - u[(j-1)*nx+i]) / (2*dy) * (v[j*nx+i+1] - v[j*nx + i-1]) / (2*dx)) - ((v[(j+1)*nx+i] - v[(j-1)*nx+i]) / (2*dy)) * ((v[(j+1)*nx+i] - v[(j-1)*nx+i]) / (2*dy)) ));
        }
    }
}

//Pressure_poisson Function
void pressure_poisson(double *p,double *pn, double *b){
   int i, j;
    for(int k = 0; k < nit ; k++){
        for(i = 0; i < nx ; i++){
            for(j = 0; j < ny ; j++){
                pn[j*nx+i] = p[j*nx+i];
            }
        }
        for(i = 1; i < nx - 1 ; i++){
            for(j = 1; j < ny - 1; j++){
                p[j*nx+i] = (((pn[j*nx+i+1] + pn[j*nx+i-1]) * dy * dy + (pn[(j+1)*nx+i] + pn[(j-1)*nx+i]) * dx * dx)/ (2 * (dx * dx + dy * dy)) - dx * dx * dy * dy * b[j*nx+i] * rho / (2 * (dx *dx + dy * dy)));
            }
        }
        for(int m = 1 ; m < ny - 1 ; m++){
            p[m*nx+nx-1] = p[m*nx+nx-2]; 
            p[m*nx+0]    = p[m*nx+1];    
        }
        for(int n = 0 ; n < nx ; n++){
            p[0*nx+n] = p[1*nx+n]; 
            p[(nx-1)*nx+n] = 0.0;  
        }
    }
}

//Cavity_flow Function
void cavity_flow(double *u, double *un, double *v, double *vn, double *p, double *pn, double *b){
 int i, j;
    for(i = 0; i < nx ; i++){
        for(j = 0; j < ny ; j++){
            un[j*nx+i] = u[j*nx+i];
            vn[j*nx+i] = v[j*nx+i];
        }
    }
    
    build_up_b(b, u, v);
    pressure_poisson(p, pn, b);
    
    for(i = 1; i < nx - 1 ; i++){
        for(j = 1; j < ny - 1 ; j++){
            u[j*nx+i] = (un[j*nx+i] - un[j*nx+i] * dt / dx * (un[j*nx+i] - un[j*nx+i-1]) - vn[j*nx+i] * dt / dy * (un[j*nx+i] - un[(j-1)*nx+i]) - dt / (2 * rho * dx) * (p[j*nx+i+1] - p[j*nx+i-1]) + nu * (dt / (dx * dx) * (un[j*nx+i+1] - 2*un[j*nx+i] + un[j*nx+i-1]) + dt / (dy * dy) * (un[(j+1)*nx+i] - 2*un[j*nx+i] + un[(j-1)*nx+i])));
            v[j*nx+i] = ( vn[j*nx+i] - un[j*nx+i] * dt / dx * (vn[j*nx+i] - vn[j*nx+i-1]) - vn[j*nx+i] * dt / dy * (vn[j*nx+i] - vn[(j-1)*nx+i]) - dt / (2 * rho * dy) * (p[(j+1)*nx+i] - p[(j-1)*nx+i]) + nu * (dt / (dx * dx) * (vn[j*nx+i+1] - 2*vn[j*nx+i] + vn[j*nx+i-1]) + dt / (dy * dy) * (vn[(j+1)*nx+i] - 2*vn[j*nx+i] + vn[(j-1)*nx+i])));
        }
    }
    //Boundary
    for(i = 0 ; i < nx ; i++){
        u[0*nx+i] = 0.0;
        u[(nx-1)*nx+i] = 1.0;
        v[0*nx+i] = 0.0;
        v[(nx-1)*nx+i] = 0.0;
    }
    //Boundary
    for(j = 1; j < ny - 1 ; j++){
        u[j*nx+0] = 0.0;
        u[j*nx+nx-1] = 0.0;
        v[j*nx+0] = 0.0;
        v[j*nx+nx-1] = 0.0;
    }
}

//Main Function
int main(void){
  double *u,*un,*v,*vn,*p,*pn,*b;
  int size= nx*ny*sizeof(double);

  u=(double*)malloc(size);
  un=(double*)malloc(size);
  v=(double*)malloc(size);
  vn=(double*)malloc(size);
  p=(double*)malloc(size);
  pn=(double*)malloc(size);
  b=(double*)malloc(size);

  in(u,un,v,vn,p,pn,b,nx,ny);
  for(int i=1;i<nt+1;i++){
    cavity_flow(u,un,v,vn,p,pn,b);
  }

//Write in CSV
  ofstream outFile;
  outFile.open("data.csv",ios::out);
  for(int i=0;i<nx*ny;i++){
    outFile<<*(u+i)<<','<<*(v+i)<<','<<*(p+i)<<endl;
    }
  outFile.close();
  cout<<"Finish"<<endl;

  free(u);
  free(un);
  free(v);
  free(vn);
  free(p);
  free(pn);
  free(b);
return 0;
}

// #Using Python to generate figure.
//import numpy
//import pandas
//from matplotlib import pyplot, cm
//from mpl_toolkits.mplot3d import Axes3D
//%matplotlib inline
//nx = 41
//ny = 41
//nt = 700
//nit = 50
//c = 1
//dx = 2 / (nx - 1)
//dy = 2 / (ny - 1)
//x = numpy.linspace(0, 2, nx)
//y = numpy.linspace(0, 2, ny)
//X, Y = numpy.meshgrid(x, y)
//rho = 1
//nu = .1
//dt = .001
//u = numpy.zeros((ny, nx))
//v = numpy.zeros((ny, nx))
//p = numpy.zeros((ny, nx)) 
//b = numpy.zeros((ny, nx))

//#Load csv data
//data=pandas.read_csv("data.csv")
//#u=data.copy
//uu=data.u.values
//vv=data.v.values
//pp=data.p.values

//#Tranfrom data to u,v,p
//for i in range(1,ny):
//    for j in range(1,nx):
//        u[i,j]=uu[ny*(i-1)+j]
//        v[i,j]=vv[ny*(i-1)+j]
//        p[i,j]=pp[ny*(i-1)+j]

//fig = pyplot.figure(figsize=(11,7), dpi=100)
//# plotting the pressure field as a contour
//pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)  
//pyplot.colorbar()
//# plotting the pressure field outlines
//pyplot.contour(X, Y, p, cmap=cm.viridis)  
//# plotting velocity field
//pyplot.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
//#label
//pyplot.xlabel('X')
//pyplot.ylabel('Y');