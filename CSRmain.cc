#include <iostream>
using namespace std;

#include "mathTools.h"
#include "CSRmatrix.h"
#include "gnuplot.h"



int main(){

    srand(time(NULL));
    
    CSRmatrix M5;
    M5.generate(5,5,40);
    //~ M5.infoCSR(GET_VARIABLE_NAME(M5));   
    //~ M5.spy(GET_VARIABLE_NAME(M5));
    M5.save("M5.dat");
    
    CSRmatrix M6;
    //~ M6 = M5.transpose();
    //~ M6.generate(7,4,40);
    //~ M6.infoCSR(GET_VARIABLE_NAME(M6));   
    M6.read("M5.dat");
    
    if ( M5 == M6 ) {
        cout << "M5 = M6" << endl;
    } else {
        cout << "M5 != M6" << endl;
    }
    
    if ( M5 != M6 ) {
        cout << "M5 != M6" << endl;
    } else {
        cout << "M5 = M6" << endl;
    }
    
    //~ M5 *= M6;
    //~ M5.infoCSR(GET_VARIABLE_NAME(M5));   
    CSRmatrix M7;
    M7 = M5*M6;
    //~ M7.infoCSR(GET_VARIABLE_NAME(M7));   
    //~ M5.infoCSR(GET_VARIABLE_NAME(M5));   
    //~ M6.infoCSR(GET_VARIABLE_NAME(M6));   
    if ( M5 == M7 ) {
        cout << "M5 = M7" << endl;
    } else {
        cout << "M5 != M7" << endl;
    }
      
    if ( M5 != M7 ) {
        cout << "M5 != M7" << endl;
    } else {
        cout << "M5 = M7" << endl;
    }
    
    M7.spy("M7");
    
    
    //~ vector<double> v1({1.0, 3.0, 2.0, 0.0, 1.0});
    //~ vector<double> v2;
     
     //~ v2 = M5 * v1;
     //~ cout << "v1 =" << endl << "[ ";
     //~ for (int i(0); i<5; ++i) {
         //~ cout << "  " << v1[i] << endl;
     //~ }
     //~ cout << "]" << endl;
     //~ M5.infoCSR(GET_VARIABLE_NAME(M5));
     //~ cout << M5;
     
     
     //~ cout << "v2 =" << endl << "[";
     //~ for (int i(0); i<5; ++i) {
         //~ cout << " " << v2[i] << endl << " ";
     //~ }
     //~ cout << "]" << endl;
    
    //~ CSRmatrix M6;
    //~ M6.generate(20000,20000,1.5,false);
    //~ M6.infoCSR(GET_VARIABLE_NAME(M6));    
    //~ M6.save("M27.dat");
    //~ start = chrono::high_resolution_clock::now(); 
    //~ M6.read("M6.dat");
    //~ finish = chrono::high_resolution_clock::now();    
    //~ elapsed = finish - start;
    //~ cout << "Elapsed time: " << elapsed.count() << " s\n";
        
    //~ M5.spy(GET_VARIABLE_NAME(M5));
    //~ M6.spy(GET_VARIABLE_NAME(M6));
    
    //~ CSRmatrix M3;
    //~ M3.generate(20000,20000,1.2,false);   
    //~ M3.save("M37.dat");
    
    //~ M5 += M6;
    //~ M5.infoCSR(GET_VARIABLE_NAME(M5));        
    //~ M5.spy(GET_VARIABLE_NAME(M5));

    //~ start = chrono::high_resolution_clock::now(); 
    //~ CSRmatrix M7;
    //~ M7 = M5+M6;    
    //~ finish = chrono::high_resolution_clock::now();    
    //~ elapsed = finish - start;
    //~ cout << "Elapsed time: " << elapsed.count() << " s\n";

    
    
    //~ auto finish = chrono::high_resolution_clock::now();
    
    //~ chrono::duration<double> elapsed = finish - start;
    //~ cout << "Elapsed time: " << elapsed.count() << " s\n";
    
    //~ M7.infoCSR(GET_VARIABLE_NAME(M7));        
    //~ M7.spy(GET_VARIABLE_NAME(M7));
    

    
    return 0;
}
