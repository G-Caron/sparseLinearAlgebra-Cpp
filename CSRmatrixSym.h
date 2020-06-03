#include "CSRmatrix.h"


class CSRmatrixSym : public CSRmatrix {

    public :
    
        // Constructeur
        using CSRmatrix :: CSRmatrix;
        /* Grace à la norma C++11, il est possible d'utiliser les constructeurs de la classe parente,
         * mais la sous-classe NE DOIT PAS avoir d'arguments supplémentaires car il ne seront pas
         * construits. A utiliser avec précautions. On a donc construit implicitement les constructeurs :
         *   CSRmatrixSym();
         *   CSRmatrixSym(int nbRows, int nbCols, int nnz);
         *   CSRmatrixSym(CSRmatrix const&);
         * A verifier: si le constructeur de copie est bon, car la norme c++ veut qu'il soit appelé
         * explicitement. */ 
        
        // methodes specialisées
        //~ void infoCSR(string const& name = "") const;
        //~ ostream& affiche(ostream&) const;
        //~ void spy(string const& name = "") const;
        //~ CSRmatrix& generate(int const, int const, double const);
        //~ CSRmatrix transpose() const;
        //~ void save(string const&) const;
        //~ CSRmatrix& read(string const&);
        
        //~ CSRmatrix& operator+=(CSRmatrix const&);
        //~ CSRmatrix& operator-=(CSRmatrix const&);
        //~ CSRmatrix& operator*=(CSRmatrix const&);
        //~ bool operator==(CSRmatrix const&) const;
        //~ vector<double> operator*(vector<double> const&) const;
    
    protected :
    
    private :

}; 
