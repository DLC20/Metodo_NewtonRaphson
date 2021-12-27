#include <math.h>
#include<iostream>
#define nqq 2
#define nqqq 3
using namespace std;

double jmf1, jmf2, jmf3;
  
  int inversa23(double mat[nqqq][nqqq],double FF,double FFF, double FFFF){
  
  int m=nqqq;
    int i, j, z, h, d;
  long double  a[m][m*2];
  

  long double matriz[m][m*2], dif, mult, f;
  for (i=0; i<m; i++) {
    for (j=0; j<m*2; j++) {
      matriz[i][j]=0;
    }
  }
  j=m;
  for (i=0; i<m; i++) {
    matriz[i][j]=1;
    j++;
  }

  for (i=0; i<m; i++) {
    for (j=0; j<m; j++) {
       matriz[i][j]=mat[i][j];   
    }
  }
  
  double cro=((matriz[0][0]*matriz[1][1]*matriz[2][2])+(matriz[0][1]*matriz[1][2]*matriz[2][0])+(matriz[0][2]*matriz[1][0]*matriz[2][1]));
  double cru=((matriz[2][2]*matriz[1][0]*matriz[0][1])+(matriz[2][1]*matriz[1][2]*matriz[0][0])+(matriz[2][0]*matriz[1][1]*matriz[0][2]));

  f=cro-cru;

  if(f>0 || f<0){

  for (z=0; z<m; z++) {
    dif=matriz[z][z];
      
    //Cambio de renglon
    if(dif==0){
      for(d=0;d<m*2;d++){
        a[0][d]=matriz[1][d];
        matriz[1][d]=matriz[0][d];
        matriz[0][d]=a[0][d];  
      }
        
    }
    
    dif=matriz[z][z];
  
    for (h=0; h<m*2; h++) {
      matriz[z][h]=matriz[z][h]/dif;
    }
    for (i=0; i<m; i++) {
      if (i!=z){
        mult= -matriz[i][z];
        for (j=0; j<m*2; j++) {
          matriz[i][j]=matriz[z][j]*mult + matriz[i][j];
        }
      }
    }
  }

  }else{
    cout<<"\n\tNo se puede calcular la inversa de la matriz Jacobiana ya que el determinante es 0\n";
    return 1;
  }
  
  jmf1=matriz[0][3]*FF+matriz[0][4]*FFF+matriz[0][5]*FFFF;
  jmf2=matriz[1][3]*FF+matriz[1][4]*FFF+matriz[1][5]*FFFF;
  jmf3=matriz[2][3]*FF+matriz[2][4]*FFF+matriz[2][5]*FFFF; 
  return 0;
  }
  
int inversa2(double mat[nqq][nqq],double FF,double FFF){
  
  int m=nqq;
    int i, j, z, h;
  double  a[m][m*2];
  double matriz[m][m*2], dif, mult, f;
  for (i=0; i<m; i++) {
    for (j=0; j<m*2; j++) {
      matriz[i][j]=0;
    }
  }
  j=m;
  for (i=0; i<m; i++) {
    matriz[i][j]=1;
    j++;
  }

  for (i=0; i<m; i++) {
    for (j=0; j<m; j++) {
       matriz[i][j]=mat[i][j];   
    }
  }
  
  double cro=(matriz[0][0]*matriz[1][1]);
  double cru=(matriz[1][0]*matriz[0][1]);

  f=cro-cru;

  if(f>0 || f<0){

  for (z=0; z<m; z++) {
    dif=matriz[z][z];
    
    //Cambio de renglon    
    if(dif==0){
      for(int d=0;d<m*2;d++){
        a[0][d]=matriz[1][d];
        matriz[1][d]=matriz[0][d];
        matriz[0][d]=a[0][d];  
      }
        
    }
    
    dif=matriz[z][z];  
    
    for (h=0; h<m*2; h++) {
      matriz[z][h]=matriz[z][h]/dif;
    }
    for (i=0; i<m; i++) {
      if (i!=z){
        mult= -matriz[i][z];
        for (j=0; j<m*2; j++) {
          matriz[i][j]=matriz[z][j]*mult + matriz[i][j];
        }
      }
    }
  }

  }else{
    cout<<"\n\tNo se puede calcular la inversa de la matriz Jacobiana, el determinante es 0\n";
    return 1;
  }
  
  jmf1=matriz[0][2]*FF+matriz[0][3]*FFF;
  jmf2=matriz[1][2]*FF+matriz[1][3]*FFF;
  return 0;
  }

void newton(int op){
  
  int it;
  long double tol, Er, Ex, Ey, nx, ny, f1x, f1y, f2x, f2y, F1, F2, x, y, z, Ez, nz, f1z, f2z, f3x, f3y, f3z, F3; 
    
  if(op>2){
    cout<<"\n\tEScribe los Puntos Iniciales (x, y, z)\n";
    cout<<"x: ";
    cin>>x;
    cout<<"y: ";
    cin>>y;  
    cout<<"z: ";
    cin>>z;
  }else{
    cout<<"\n\n\tEscribe los Puntos Iniciales (x, y)\n";
    cout<<"\tx: ";
    cin>>x;
    cout<<"\ty: ";
    cin>>y;  
  }
  
  cout<<"\n\tDar el numero Maximo de Iteraciones:";
  cin>>it;
  cout<<"\n\tDar el valor de Tolerancia:";
  cin>>tol;
  cout<<"\n";
  
  for(int freno=0;freno<it;freno++){
    cout<<"\tIteracion no."<<freno+1<<"___________________________________";
    
    switch(op){
      case 1:
        f1x=2*x+y;
        f1y=x;
        f2x=3*(y*y);
        f2y=6*x*y+1;
        F1=((x*x)+(x*y))-10;
        F2=(y+(3*x*(y*y)))-50;
        break;
      case 2:
        f1x=2*x;
        f1y=2*y;
        f2x=-exp(x);
        f2y=-2;
        F1=(pow(x, 2))+(pow(y, 2))-9;
        F2=-exp(x)-(2*y)-3;
        break;
      case 3:
        f1x=(4*x)-4;
        f1y=2*y;
        f1z=(6*z)+6;
        f2x=2*x;
        f2y=(2*y)-2;
        f2z=4*z;
        f3x=(6*x)-12;
        f3y=2*y;
        f3z=(-6*z);
        F1=((2*(pow(x,2)))-(4*x)+(pow(y,2))+(3*(pow(z,2)))+(6*z)+2);
        F2=((pow(x,2))+(pow(y,2))-(2*y)+(2*(pow(z,2)))-5);
        F3=((3*(pow(x,2)))-(12*x)+(pow(y,2))-(3*(pow(z,2)))+8);
        break;
      case 4:
        f1x=(2*x)-4;
        f1y=2*y;
        f1z=0;
        f2x=(2*x)-1;
        f2y=-12;
        f2z=0;
        f3x=(6*x)-12;
        f3y=2*y;
        f3z=(-6*z);
        F1=(pow(x,2)-(4*x)+(pow(y,2)));
        F2=((pow(x,2))-x-(12*y)+1);
        F3=((3*(pow(x,2)))-(12*x)+(pow(y,2))-(3*(pow(z,2)))+8);
        break;      
    }  
    
    if(op>2){
      double Jac[nqqq][nqqq];
      Jac[0][0]=f1x;
      Jac[0][1]=f1y;
      Jac[0][2]=f1z;
      Jac[1][0]=f2x;
      Jac[1][1]=f2y;
      Jac[1][2]=f2z;
      Jac[2][0]=f3x;
      Jac[2][1]=f3y;
      Jac[2][2]=f3z;  
      
      if(inversa23(Jac,F1,F2,F3)==1){
      freno=it;
      }else{
        
      nx=x-jmf1;
      ny=y-jmf2;
      nz=z-jmf3;
      
      cout<<"\n\tPunto x= "<<x;
      cout<<"\n\tPunto y= "<<y;
      cout<<"\n\tPunto z= "<<z;
      
      cout<<"\n\t\tMatriz Jacobiana:\n";
      for (int i=0; i<nqqq; i++)
      {
        for (int j=0; j<nqqq; j++)
        {
          cout<<"\t|"<<Jac[i][j]<<"|";
          cout<<"   ";
        }
      cout<<"\n";
      }
      cout<<"\n\t\tf(x)";
      cout<<"\n\tf1(x,y,z)= "<<F1;
      cout<<"\n\tf2(x,y,z)= "<<F2; 
      cout<<"\n\tf3(x,y,z)= "<<F3;
      
      cout<<"\n\t\tInversa de la matriz Jacobiana:\t";
      cout<<"\n\tf1= "<<jmf1;
      cout<<"\n\tf2 = "<<jmf2;
      cout<<"\n\tf3 ="<<jmf3;  
      
      long double abs, absy, noy, nox, absz, noz;
      
      if(it>1){
            
      if((nx-x)<0){
        abs=(nx-x)*-1;
      }else{
        abs=nx-x;
      }
      if(nx<0){
        long double nox=nx*-1;
        Ex=abs/nox;
      }else{
        nox=nx;
        Ex=abs/nox;
      }

      if((ny-y)<0){
                
      long double absy=(ny-y)*-1;
          
      }else{
        absy=ny-y;
      }
      if(ny<0){
        long double noy=ny*-1;
        Ey=absy/noy;
      }else{
        noy=ny;
        Ey=absy/noy;
      }
            
      if((nz-z)<0){
    
      long double absz=(nz-z)*-1;
        
      }else{
        absz=nz-z;
      }
      if(nz<0){
        long double noz=nz*-1;
        Ez=absz/noz;
      }else{
        noz=nz;
        Ez=absz/noz;
      }
  
      if(Ex>Ey){
        Er=Ex;
      }else{
        Er=Ey;
      }if(Er>Ez){
        Er=Er;
      }else{
        Er=Ez;
      }

cout<<"\n\tEr= "<<Er;
  
      if(Er<tol){
        cout<<"\n\tConverge en el punto ("<<x<<", "<<y<<", "<<z<<") ";
        cout<<"de la iteracion: "<<freno+1<<"\n\n\n";
        freno=it;      
      }else{
        cout<<"\n\tSe deben hacer mas iteraciones, el metodo aun no converge\n\n\n";
      }  
      
      }
      
      x=nx;
      y=ny;
      z=nz;  
      }
      
    }else{
    double Jac[nqq][nqq];
    Jac[0][0]=f1x;
    Jac[0][1]=f1y;
    Jac[1][0]=f2x;
    Jac[1][1]=f2y;
        
    if(inversa2(Jac,F1,F2)==1){
      freno=it;
    }else{
      

    nx=x-jmf1;
    ny=y-jmf2;
    
    cout<<"\n\tPunto x= "<<x;
    cout<<"\n\tPunto y= "<<y;
      
    cout<<"\n\t\t[Matriz Jacobiana]\n";
    for (int i=0; i<nqq; i++)
    {
      for (int j=0; j<nqq; j++)
      {
        cout<<"\t|"<<Jac[i][j]<<"|";
        cout<<"   ";
      }
    cout<<"\n";
    }
    cout<<"\n\t\t[f(x)]";
    cout<<"\n\tf1(x,y)= "<<F1;
    cout<<"\n\tf2(x,y)= "<<F2; 

        cout<<"\n\n\t\t[Inversa de la Matriz Jacobiana]";
    cout<<"\n\tf1= "<<jmf1;
    cout<<"\n\tf2= "<<jmf2;
          
    long double abs, absy, noy, nox;
    
    if(it>1){
    
    
    if((nx-x)<0){
      abs=(nx-x)*-1;
    }else{
      abs=nx-x;
    }
    if(nx<0){
      long double nox=nx*-1;
      Ex=abs/nox;
    }else{
      nox=nx;
      Ex=abs/nox;
    }

    if((ny-y)<0){
      
    long double absy=(ny-y)*-1;
        
    }else{
      absy=ny-y;
    }
    if(ny<0){
      long double noy=ny*-1;
      Ey=absy/noy;
    }else{
      noy=ny;
      Ey=absy/noy;
    }
      
    if(Ex>Ey){
      Er=Ex;
    }else{
      Er=Ey;
    }
    
    cout<<"\n\t[Error= "<<Er<<"]";
    
    if(Er<tol){
        
      cout<<"\n\tEl Sistema converge en el punto ("<<x<<", "<<y<<") ";
      cout<<"de la iteracion: "<<freno+1<<"\n\n\n";
      freno=it;
    }else{
      cout<<"\n\tSe deben hacer mas iteraciones, el metodo aun no converge\n\n\n";
    }  
    
    }
    
    x=nx;
    y=ny;
      
    }  
    }
  
  }
}

void SE(int op){
  int resp;
  
  do{
    switch(op){
    case 1:
      cout<<"\n\tElegiste el Sistema (1) \n\tf1(x, y)= x^2 + xy - 10 = 0\n\tf2(x, y)= y + 3xy^2 - 50 = 0";
      newton(op);
      break;
    case 2:
      cout<<"\n\tElegiste el Sistema (2)\n\tf1(x, y)= x^2 + y^2 - 9 = 0\n\tf2(x, y)= -e^x - 2y - 3 = 0";
      newton(op);
      break;
    case 3:
      cout<<"\n\tElegiste el Sistema (3)\n\tf1(x, y, z)= 2x^2 - 4x + y^2 + 3z^2 + 6z + 2 = 0\n\tf2(x, y, z)= x^2 + y^2 -2y + 2z^2 -5 = 0\n\tf3(x, y, z)= 3x^2 - 12x + y^2 -3z^2 + 8= 0";
      newton(op);
      break;
    case 4:
      cout<<"\n\tElegiste el Sistema (4)\n\tf1(x, y, z)= x^2 + 4x + y^2 = 0\n\tf2(x, y, z)= x^2 -x -12y + 1 = 0\n\tf3(x, y, z)= 3x^2 -12x + y^2 -3z^2 + 8 = 0";
      newton(op);
      break;  
    
    }  
    cout<<"\n\tDesea probar con otros puntos? (1=Si, 2=No)";
    cin>>resp;  
  }while(resp==1);
}
/*void caratula(){
  cout<<"\tMETODO DE NEWTON-RAPHSON\n\n";
  cout<<"\n\t\tIntegrantes: \n ";
  cout<<"\n\tHernández Cabrera Carlos Joaquín";
  cout<<"\n\tRomo Rangel Edgar Otniel";
  cout<<"\n\tOrtega Lobato Dulce Ivonne";
  cout<<"\n\tPresione Enter, para entrar al Menu...";
  getchar();
}*/

int main(){
  
  
  /*caratula();*/
  
  
  int opcion, resp;  
  do{
    system("cls");
    cout<<"\tMETODO DE NEWTON-RAPHSON\n\n";
    cout<<"\t\tMenu\n";
    cout<<"\n\t1.Sistema de Ecuaciones no Lineales\n\tf1(x, y)= x^2 + xy - 10 = 0\n\tf2(x, y)= y + 3xy^2 - 50 = 0";
    cout<<"\n\n\t2.Sistema de Ecuaciones no Lineales\n\tf1(x, y)= x^2 + y^2 - 9 = 0\n\tf2(x, y)= -e^x - 2y - 3 = 0";
    cout<<"\n\n\t3.Sistema de Ecuaciones no Lineales\n\tf1(x, y, z)= 2x^2 - 4x + y^2 + 3z^2 + 6z + 2 = 0\n\tf2(x, y, z)= x^2 + y^2 -2y + 2z^2 -5 = 0\n\tf3(x, y, z)= 3x^2 - 12x + y^2 -3z^2 + 8= 0";
    cout<<"\n\n\t4.Sistema de Ecuaciones no Lineales\n\tf1(x, y, z)= x^2 + 4x + y^2 = 0\n\tf2(x, y, z)= x^2 -x -12y + 1 = 0\n\tf3(x, y, z)= 3x^2 -12x + y^2 -3z^2 + 8 = 0";
    cout<<"\n\n\t5. Salir";
    cout<<"\n\tPresione el numero del Sistema de Ecuaciones que desea resolver [1-5]: ";
    cin>>opcion;
    
    
    if(opcion<5){
      SE(opcion);  
      cout<<"\tResolver otro sistema (1=Si, 2=No)";
      cin>>resp;  
    }
    if(resp==2||opcion==5){
      cout<<"\n\tGracias por utilizar nuestro programa";
      cout<<"\n\t\tIntegrantes: \n ";
      cout<<"\n\tHernández Cabrera Carlos Joaquín";
      cout<<"\n\tRomo Rangel Edgar Otniel";
      cout<<"\n\tOrtega Lobato Dulce Ivonne";
      resp=2;
    }
  }while(resp==1);  
}
