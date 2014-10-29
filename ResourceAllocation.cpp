#include <iostream>
#include <fstream>
#include <cmath>
#include <queue>
#include <vector>
#include <string>
using namespace std;

const int len = 14542;
const float target = 0.7;
const double mmin = 1e-20;

int maptocell[len];

void get_init(const double (&P)[len], const double (&Tchs)[len], double (&x)[len]){//produce the optimum point in Real Number Space
    double c1=0.0,c2=0.0;
    for (int i = 0; i < len; i++){
        c1+=P[i]/target-Tchs[i];
        c2+=P[i]*P[i];
    }
    c1*=16.0;
    double lamda=c1/c2;
    for (int i = 0; i < len; i++){
        x[i]=-1*P[i]*P[i]*lamda/128.0+(P[i]/target-Tchs[i])/8.0;
        //fprintf(fp,"%lf %lf %lf\n", P[i], Tchs[i], x[i]);
    }
}

template <typename T>
T sum(T a[len]){
        T s=0;
        for (int i = 0; i < len; i++)
            s+=a[i];
        return  s;
}

template<typename T>
double f(const T x[len], const double P[len], const double Tchs[len]){
    double sum = 0.0,t = 0.0;
    for (int i = 0; i < len; i++){
        //t = x[i]-(P[i]/target-Tchs[i])/8.;
        //sum+=64.*t*t/(P[i]*P[i]);
        t = P[i]/(Tchs[i]+8*x[i])-target;
        sum+=t*t; 
    }
    return sum;
}

template<typename T>
double f_x(const T x, const double P, const double Tchs){
       double t = P/(Tchs+8*x)-target;
       return t*t;
}

double partial_f(int x, double P, double Tchs){
      return 128.0*x/(P*P)-16/(P*target)+16*Tchs/(P*P);
}

struct node {
       double dfx,dx;
       int idx;
       node(int a, double b, double c):idx(a),dfx(b),dx(c){}
};

bool operator < (const node &a, const node &b){        //big root heap
     return a.dfx < b.dfx;
}

void greedy2(const double (&P)[len], const double (&Tchs)[len], double (&c)[len], int (&x)[len], int day){
        get_init(P, Tchs, c);
        cout << sum(c) << endl;
        cout << f(c,P,Tchs) << endl;
        bool h[len];
        memset(h,0,sizeof(h));
        priority_queue<node> incq,decq;
        for (int i = 0; i < len; i++) {
            if(c[i]>=0){
                    if(abs(Tchs[i]+8*(int(c[i]-mmin)+1))>0.0001)
                        incq.push(node(i,f_x(c[i],P[i],Tchs[i])-f_x(int(c[i]-mmin)+1,P[i],Tchs[i]),int(c[i]-mmin)+1-c[i]));
                    if(abs(Tchs[i]+8*(int(c[i])))>0.0001)
                       decq.push(node(i,f_x(c[i],P[i],Tchs[i])-f_x(int(c[i]),P[i],Tchs[i]),c[i]-int(c[i])));
            }
            else{
                 if(abs(Tchs[i]+8*(int(c[i])))>0.0001)
                      incq.push(node(i,f_x(c[i],P[i],Tchs[i])-f_x(int(c[i]),P[i],Tchs[i]),int(c[i])-c[i]));
                 if(abs(Tchs[i]+8*(int(c[i]+mmin)-1))>0.0001)
                      decq.push(node(i,f_x(c[i],P[i],Tchs[i])-f_x(int(c[i]+mmin)-1,P[i],Tchs[i]),c[i]-(int(c[i]+mmin)-1)));
            }
        }
        
        
        double resi=0.;
        if(incq.top().dfx>decq.top().dfx){
                          resi = incq.top().dx;
                          h[incq.top().idx] = 1;
                          x[incq.top().idx] = c[incq.top().idx]>=0?int(c[incq.top().idx]-mmin)+1:int(c[incq.top().idx]);
                          incq.pop();
                          for (int i = 0; i < len-1; i++){
                              node t(0,0,0);
                              if(resi>0.&&!decq.empty()){
                                          while(!decq.empty()){
                                                 t = decq.top();
                                                 if(h[t.idx]==0) break;
                                                 decq.pop();
                                          }
                                          
                                          if(decq.empty()){
                                                           i--;
                                                           continue;
                                          }
                                          
                                          resi -= t.dx;
                                          x[t.idx] = c[t.idx]>=0?int(c[t.idx]):int(c[t.idx]+mmin)-1;
                                          h[t.idx] = 1;
                                          decq.pop();
                              }
                              else{
                                   while(!incq.empty()){
                                            t = incq.top();
                                            if(h[t.idx]==0) break;
                                            incq.pop();
                                   }
                                   
                                   if(incq.empty()){
                                                    i--;
                                                    continue;
                                   }
                                   cout << resi << endl;
                                   resi += t.dx;
                                   cout << resi << endl;
                                   x[t.idx] = c[t.idx]>=0?int(c[t.idx]-mmin)+1:int(c[t.idx]);
                                   h[t.idx] = 1;
                                   incq.pop();
                              }
                          }
        }
        else{
             resi -= decq.top().dx;
             h[decq.top().idx] = 1;
             x[decq.top().idx] = c[decq.top().idx]>=0?int(c[decq.top().idx]):int(c[decq.top().idx]+mmin)-1;
             decq.pop();
             for (int i = 0; i < len-1; i++){
                 node t(0,0,0);
                 if(resi>0.&&!decq.empty()){
                                          while(!decq.empty()){
                                                 t = decq.top();
                                                 if(h[t.idx]==0) break;
                                                 decq.pop();
                                          }
                                          
                                          if(decq.empty()){
                                                           i--;
                                                           continue;
                                          }
                                          cout << resi << endl;
                                          resi -= t.dx;
                                          cout << resi << endl;
                                          x[t.idx] = c[t.idx]>=0?int(c[t.idx]):int(c[t.idx]+mmin)-1;
                                          h[t.idx] = 1;
                                          decq.pop();
                 }
                 else{
                                   while(!incq.empty()){
                                            t = incq.top();
                                            if(h[t.idx]==0) break;
                                            incq.pop();
                                   }
                                   
                                   if(incq.empty()){
                                                    i--;
                                                    continue;
                                   }
                                   cout << resi << endl;
                                   resi += t.dx;
                                   cout << resi << endl;
                                   x[t.idx] = c[t.idx]>=0?int(c[t.idx]-mmin)+1:int(c[t.idx]);
                                   h[t.idx] = 1;
                                   incq.pop();
                 }
             }
        }
        
        cout << "x: " <<sum(x) << endl;
        cout << resi << endl;
        cout << f(x,P,Tchs) << endl;
        
        FILE *fp;
        fp = fopen("./greedy2_data/result.txt","a");
        fprintf(fp, "%d %lf\n", day, f(x,P,Tchs));
        fclose(fp);
        
        
        char filename[50],*path="./greedy2_data/Tch",*filetype = ".txt";
        
        int plen=strlen(path);
        for (int i = 0; i < plen; i++)    filename[i] = path[i];
        if(day/10) {
                   filename[plen] = char(day/10+'0');
                   filename[plen+1] = char(day%10+'0');
                   int i = 2;
                   for(; i <6; i++)
                         filename[plen+i] = filetype[i-2];
                   filename[plen+i] = '\0';
        }
        else{
             filename[plen] = char(day%10+'0');
             int i = 1;
             for (; i < 5; i++)
                 filename[plen+i] = filetype[i-1];
             filename[plen+i] = '\0';
        }
        
        
        fp = fopen(filename,"w");
        //fprintf(fp, "%lf\n", f(x,P,Tchs));
        for (int i = 0; i < len-1; i++)
            fprintf(fp, "%d %lf\n", maptocell[i], x[i]*8+Tchs[i]);
        fprintf(fp, "%d %lf", maptocell[len-1], x[len-1]*8+Tchs[len-1]);
        fclose(fp);
}

int main(){
    double P[len],Tchs[len],y,s=0.,c[len];
    int temp[len],x[len];
    
    FILE *fp;
    
    char filename[50],*filetype = ".txt", *path="./greedy2_data/P", filename_TCH[50], *TCH_path="./greedy2_data/Tch";
    int plen = strlen(path),tlen = strlen(TCH_path);
    for (int i = 0; i < plen; i++)
        filename[i] = path[i];
    for (int i = 0; i < tlen; i++)
        filename_TCH[i] = TCH_path[i];
    for (int i = 1; i <= 31; i++){
        if(i/10){
                 filename[plen] = filename_TCH[tlen] = char(i/10+'0');
                 filename[plen+1] = filename_TCH[tlen+1] = char(i%10+'0');
                 for(int i = 2; i <6; i++)
                         filename[plen+i] = filename_TCH[tlen+i] = filetype[i-2];
        }
        else{
             filename[plen] = filename_TCH[tlen] = char(i+'0');
             for (int i = 1; i < 5; i++)
                 filename[plen+i] = filename_TCH[tlen+i] = filetype[i-1];
        } 
        cout << filename << " " << filename_TCH << endl;
        
        fp = fopen(filename,"r");
        if(fp == NULL) exit(1);
        int readline = 0;
        while(1){
                 fscanf(fp,"%d %lf", &maptocell[readline], &P[readline]);
                 readline++;
                 if(feof(fp))
                             break;
        }
        //cout << readline << endl;
        fclose(fp);
        
        fp = fopen(filename_TCH, "r");
        if(fp==NULL) exit(1);
        readline = 0;
        while(1){
                 fscanf(fp, "%d %lf", &maptocell[readline], &Tchs[readline]);
                 readline++;
                 if(feof(fp))
                             break;
        }
        if(i==1)
        greedy2(P,Tchs,c,x,i+1);
    }
    
    system("pause");
}
