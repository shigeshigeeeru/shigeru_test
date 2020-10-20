//reference /Desktop/RS232/PIDemulator 
#include<stdio.h>
#include<cmath>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

#define XSIZE 0
#define YSIZE 0
#define ZSIZE 28
#define STEP 30000
//#define TIME 0.1
#define INTERVAL 0.01

#define INITIAL_TEMP 25.0
#define TARGET_TEMP -10.0

#define TRUE 1
#define FALSE 0

#define N_REF 5
#define PID_STEP 100
double K_P;
double K_I;
double K_D;
//#define K_P 0.01
//#define K_I 0
//#define K_D 0.1


//std::string material = "alminum";
int runnumber=0;
double x_length = 0.05;//m
double y_length = 0.05;//m
double z_length = 0.04;//m
double A = x_length * y_length; 
double row = A/z_length;//mm
//double h = z_length/ZSIZE

double density_Al =2688;// kg/m^3
double hinetsu_Al =0.905;// kJ/(kg K)
double lambda = 237.0;//(W/(m K))
double Cv = 1000.0 * density_Al *hinetsu_Al;

double V=0.0;
double T[ZSIZE+1][STEP];
//double x[XSIZE+1][STEP];
//double y[YSIZE+1][STEP];
//double z[ZSIZE+1][STEP];
//double position[SIZE][SIZE][SIZE];
//define 
/*
class Aluminum{
private: 
  const double m_density =2688;// kg/m^3
  const double m_specific_heat =0.905;// kJ/(kg K)
  const double m_conductivity = 237.0;//(W/(m K))
  double m_Cv = 1000.0 * denisity *Specific_heat;
public:
  Aluminum();
  virtual ~Aluminum();

  const double& Den(){return m_density;}
  const double& Spe(){return m_specific_heat;}
  const double& Con(){return m_conductivity;}
   
}

Aluminum::~Aluminum(){
  std::cout << "Aluminum was deleted" << std::endl;
}
*/
/*
function(double Tx_pre ,double Ty_pre ,double Tz_pre , double Tx_post ,double Ty_post ,double Tz_post ,double T){
    double Q;

    //developing...
	
}
*/


void inverse(int t){
  int n=ZSIZE;
  double h=z_length/double(n);
  double d=INTERVAL;
  double alpha=lambda/density_Al/hinetsu_Al/1000;
  double g = 2.0 *h*h/d/alpha + 2.0;
  //std::cout << "h=" << h << std::endl;
  //cin >> n;
  double a[n-1][n];
  double b[n-1];
  for(int i=0;i<n-1;i++){
    for(int j=0;j<=n-1;j++){
      a[i][j]=0;
      if(i==j){
	 a[i][j]=g;
      }
      if(i==j+1 or i==j-1){
	 a[i][j] = -1;
      }
      //std::cout << i << j <<" "<< a[i][j] << std::endl;
    }
  }

  for(int i=0;i<n-1;i++){
    b[i]=T[i+2][t-1] + (g-4)*T[i+1][t-1] + T[i][t-1]; 
    //std::cout << i << " "<< b[i] << std::endl;
  
  }

  for(int i=1;i<n-2;i++){ 
    a[i][n-1]=b[i];
    //std::cout << i << " "<< b[i] << std::endl;  
  }
  a[0][n-1]= b[0] + T[0][t];
  a[n-2][n-1] = b[n-2] + T[n][t];

  

  //( 00  01 = 02);
  //( 10  11 = 12);
 
  for(int i=0;i<n-1;i++){
    int maxLine =i;
    for(int j=i;j<n-1;j++) if(std::abs(a[maxLine][i])<std::abs(a[j][i])) maxLine=j;
    for(int j=i;j<=n-1;j++) std::swap(a[i][j],a[maxLine][j]);

    double beginNum = a[i][i];
    if(!beginNum){
      std::cout << "can not calc" << std::endl;
      //return 0;
    }

    for(int j=i;j<=n-1;j++) a[i][j] /=beginNum;

    for(int j=i+1;j<n-1;j++){
      double beginDelete = a[j][i];
      for(int k=i;k<=n-1;k++) a[j][k] -= beginDelete*a[i][k];
    }
  }
   
  for(int i=n-2;i>=0;i--){
    for(int j=i+1;j<n-1;j++){
      a[i][n-1] -= a[j][n-1]*a[i][j];
      //for(int k=0;k<n-1;k++) a[i][k] -= a[j][k]*a[i][j];
    }
  }
  for(int i=1;i<n;i++){
    T[i][t]=a[i-1][n-1];
    //std::cout << a[i-1][n-2] <<std::endl;
  }
}

double PeltierCalculation(double T_pre){
    double T_h=INITIAL_TEMP;//heatsink
    double T_next;//cool
    double zebeck=0.8; 
    //double zebeck=0.000211;
    //double R=1.08497e-05;//resistance
    double R=27.5;//resistance
    double conduct=1.659;//Kp
    //double conduct=160.9;//Kp
    double mc=0.05472;//not correct ,need to research 
    //T_c = T_h - (zebeck*T_h + V/2.0) * V/R/conduct;
    //T_next = T_pre - ( zebeck*V*V*(zebeck*T_h + V/2.0)/R/R/conduct + V*V/R)*INTERVAL/mc;  
    //T_next = T_pre +( (-1)*zebeck*(T_h + T_pre)*V/R + 2.0*conduct * (T_h - T_pre) )* INTERVAL/mc;  
    //T_next = T_pre + (zebeck*(T_h - T_pre)*V/R + V*V/R )* INTERVAL/mc;  
    T_next = T_pre + ( (-1)*zebeck*(T_pre+273)*V/R + V*V/R/2.0 + conduct * (T_h - T_pre) )* INTERVAL/mc;  
    //T_next = T_pre + (zebeck*(T_pre+273)*V/R + V*V/R/2.0 - conduct * (T_h - T_pre) )* INTERVAL/mc;  
     
    //std::cout << T_next << std::endl;
    return T_next; 

}


void PID(int z){
    //double e = T-TARGET_TEMP;
    //double proportional; 
    double P,I,D;
    double sum_x,sum_y,sum_xy,sum_x2;    
    double proportional;
    double integral;
    double derivative;
    double deltaTime; 
    double LimitVoltage=10.0;

    if(runnumber > N_REF * PID_STEP){
      for(int i=0;i<N_REF;i++){
	integral += T[z][runnumber - i * PID_STEP];
	
    	sum_xy += INTERVAL * i * PID_STEP * (T[z][runnumber - (N_REF - i - 1) * PID_STEP] - TARGET_TEMP);
	sum_x += INTERVAL * i * PID_STEP;
	sum_y += T[z][runnumber - (N_REF - i - 1) * PID_STEP] - TARGET_TEMP;
	sum_x2 += pow(INTERVAL * i * PID_STEP ,2);
      }
      derivative=(N_REF *sum_xy - sum_x * sum_y ) / (N_REF * sum_x2 -pow(sum_x ,2)); 
      
      proportional = T[z][runnumber] - TARGET_TEMP;
      P = K_P * proportional;
      I = K_I * integral;
      D = K_D * derivative;
 
      V += P + I + D;

      if(V > LimitVoltage){
	V = LimitVoltage;
	//std::cout<<"== over LimitVoltage ==" << std::endl;
      }

      if(V < LimitVoltage*(-1)){
	V = LimitVoltage*(-1);
	//std::cout<<"== under LimitVoltage ==" << std::endl;
      }

      //std::cout << "Voltage=" << T << std::endl;
      if(runnumber>STEP) std::cout << " " << V << std::endl;
  }  
}

int main(){
//std::ofstream ofs("data/ThermalSimulation.txt");
std::ofstream ofs;
for(int p=0;p<5;p++){
for(int d=0;d<5;d++){
  std::string str_p = std::to_string(p);
  std::string str_d = std::to_string(d);
  std::string file = "data/ThermalSimulation_P" + str_p + "_D" + str_d + ".txt";
  std::ofstream ofs(file.c_str());

  K_P = 0.01*p;
  K_I = 0.0;
  K_D = 0.01*d; 
  int i,j,k;
  //intialize//
	for(k=0;k<ZSIZE+1;k++){
	  /*
	  if(i==0 or i==XSIZE){}
	  if(j==0 or j==YSIZE){}
	  if(k==0){}
	  if(k==ZSIZE){}
	  */
	  T[k][0] = INITIAL_TEMP;
	}

	//T[0][0] = TARGET_TEMP;
      	for(runnumber=1;runnumber<STEP;runnumber++){	
	  T[0][runnumber] = PeltierCalculation(T[0][runnumber-1]);
 	  //T[0][runnumber] = T[0][runnumber-1]; 
	  T[ZSIZE][runnumber]=T[ZSIZE-1][runnumber-1];
	  inverse(runnumber);// PROBLEM HERE!!!
	  //std::cout << runnumber << " "<< T[ZSIZE][runnumber] << std::endl;  


	  //ThermalCalculation(runnumber);    
	  if(runnumber%PID_STEP == 0) PID(ZSIZE);
	ofs << runnumber << " " <<  T[ZSIZE][runnumber] << " " << V << std::endl; 
	}

  std::cout << "===Simulation Completed===" <<std::endl;
  std::cout << "PID final result" << std::endl;
  std::cout << "TARGET_TEMP = " << TARGET_TEMP << std::endl;
  std::cout << "CURRENT_TEMP = " << T[ZSIZE][STEP-1] << " #circC" << std::endl;
}
}
  return 0;
}

