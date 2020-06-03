#ifndef GNUPLOT_H
#define GNUPLOT_H
#include <iostream>
#include <cstring>
using namespace std;

#define GET_VARIABLE_NAME(Variable) (void(Variable), #Variable)


class gnuplot {
    
    public:
    
        // constructeur
        gnuplot();
        
        // destructeur
        ~gnuplot();
    
    
        
        void operator() (const string& command);
    
    protected:
    
        FILE* gnuplotpipe;
};


#endif


void plotSpyInFile(int, int, int, string);
void plotSpyGUI(int, int, int, string);
