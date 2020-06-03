#include <iostream>    // std::cout, std::fixed
#include <iomanip>     // std::setprecision
#include <limits>      // std::numeric_limits
#include <cstdlib>     // srand, rand 
#include <ctime>       // time 
#include <vector>
#include <cstring>
#include <fstream>
#include <cmath>
#include "mathTools.h"

using namespace std;



class CSRmatrix{

    public :

        // Constructeur
        CSRmatrix();
        CSRmatrix(int nbRows, int nbCols, int nnz);
        CSRmatrix(CSRmatrix const&);
        
        
        // Destructeur
        ~CSRmatrix();
        
        // Getter
        int getRows() const;
        int getCols() const;
        int getNNZ() const;
        //~ bool getSymmetric() const;
        
        // setter
        
        
        // methodes
        void infoCSR(string const& name = "") const;
        ostream& affiche(ostream&) const;
        void spy(string const& name = "") const;
        CSRmatrix& generate(int const, int const, double const);
        CSRmatrix transpose() const;
        void save(string const&) const;
        CSRmatrix& read(string const&);
        
        CSRmatrix& operator+=(CSRmatrix const&);
        CSRmatrix& operator-=(CSRmatrix const&);
        CSRmatrix& operator*=(CSRmatrix const&);
        bool operator==(CSRmatrix const&) const;
        vector<double> operator*(vector<double> const&) const;
        
        
    protected :                
        
        int nbRows;        // number of rows,
        int nbCols;        // number of columns,
        int nnz;           // number of none-zeros elements,
        //~ bool symmetric;    // If the matrix is symmetric, only the low triangular part is stocked
        vector <int> Ai;   // vector of integer contening the row coordinate of the first element of each row
        vector <int> Aj;   // vector of integer contening the column coordinate of all element  
        vector <double> A; // vector of double contening the value of all element
        
        

    private :



};


// Prototype externe
const CSRmatrix operator+(CSRmatrix , CSRmatrix const&);    
const CSRmatrix operator-(CSRmatrix , CSRmatrix const&);    
const CSRmatrix operator*(CSRmatrix , CSRmatrix const&);    
bool operator!=(CSRmatrix const&, CSRmatrix const&);
ostream& operator<<(ostream&, CSRmatrix const&);
