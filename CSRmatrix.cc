//~ #include <chrono>
#include "CSRmatrix.h"
#include "gnuplot.h"
/* TODO :
 *   - test la taille produit matrice-vecteur
 *   - Optimisation de la taille des vecteurs temporaire
 *   - compression
 *   - test unitaire
 *   - vallgring callgrind
 *   - comment assigner les vecteurs soit meme, il faut des setter!!!
 *   - produit vecteur*matrice
 *   - isSymemetric, isEmpty, diag, 
 *   - Reflexion : faut il creer une classe CSRmatrixSymmetrique heritant de CSRmatrix 
 *   - creation des tableaux (utilisation des template <T> intger, double precision, complex)
 *   - modification de la dimension de la matrice (dans ce cas faut-il réinitialiser?) Prevoir les differentes possibilites
 *   - ones, zeros, det, diag, transpose, triu, tril (sont ils des constructeurs? sortes de constructeurs par copies?)
 *   - compression.
 *   - isSymetric ? 
 *   - if symmetric alors suppression de la partie triangulaire sup
 *   - matrice inverse (lapack, mumps)
 *   - insert 
 *   - extraire
 *   - Surcharge d'operateur (interne ou externe) : ^ / dans la classe ?
 *   -                                            : >=, >, <, <=, 
 *   -                                            : -a, ++, --, /=
 *   -                                            : = (par defaut, copie de surface mais il peut etre interessant de la supprimer pour les tres grosses classes)
 * 
 */ 




// ======================================================================
//                           CONSTRUCTEUR 
// ======================================================================


// constructeur par defaut par defaut
CSRmatrix::CSRmatrix() 
    : nbRows(0), nbCols(0), nnz(0)
    {}

 
// constructeur par defaut 
CSRmatrix::CSRmatrix(int nbRows, int nbCols, int nnz)
    : nbRows(nbRows), nbCols(nbCols), nnz(nnz), Ai(nbRows+1,0), Aj(nnz,0), A(nnz,0)
    {}

// constructeur par copie 
CSRmatrix::CSRmatrix(CSRmatrix const& M)
    : nbRows(M.nbRows), nbCols(M.nbCols), nnz(M.nnz) , Ai(M.Ai), Aj(M.Aj), A(M.A)
    {}

// ======================================================================
//                           DESTRUCTEUR 
// ======================================================================

CSRmatrix::~CSRmatrix() {
        Ai.clear();
        Aj.clear();
        A.clear();
}


// ======================================================================
//                              GETTER 
// ======================================================================

int CSRmatrix::getRows() const { return nbRows;}
int CSRmatrix::getCols() const { return nbCols;}
int CSRmatrix::getNNZ() const { return nnz;}



// ======================================================================
//                              METHODES 
// ======================================================================

void CSRmatrix::infoCSR(string const& name) const {
    

    if ( nbRows != 0 and nbCols != 0 and not Ai.empty() ) {      

        // Matrice
        cout << "Matrice "<< name << " ";
        //~ if (symmetric) {
            //~ cout << "symmetrique ";
        //~ } else {
            //~ cout << "non symmetrique ";
        //~ }
         cout << "de " << nbRows << " lignes, "
             << nbCols << " colonnes et de " 
             << nnz << " elements non nuls." << endl;
         
        // Densite        
        cout << "Elle a une densité de ";
        //~ if (symmetric) {
            //~ int diagElmt(0);
                
            //~ for (int i(0); i<nbRows; ++i){
                //~ if ( Aj[Ai[i+1]-1] == i ) {
                    //~ ++diagElmt;
                //~ }                
            //~ }            
            //~ cout << (100.0*(2*nnz-diagElmt))/(nbRows*nbCols);
        //~ } else {
            cout << (100.0*nnz)/(nbRows*nbCols);
        //~ }
        cout << "%." << endl;
        
        // Affichage des tableaux
        if (nnz < 300 ) {
            
            // Ai
            cout << "Ai = [";
            for (int k(0); k<=nbRows; ++k){
                cout <<  Ai[k] << " ";
            }
            cout << "]" << endl;
            
            // Aj
            cout << "Aj = [";
            for (int k(0); k<nnz; ++k){
                cout <<  Aj[k] << " ";
            }
            cout << "]" << endl;
            
            // A
            cout << "A = [";
            for (int k(0); k<nnz; ++k){
                cout <<  A[k] << " ";
            }
            cout << "]" << endl;
                       
        }       
        
         
        // affiche le nombre de chiffre, identique à std::setprecision(10)
        cout.precision(10);        
        
        // Affichage des tableaux
        if (nbRows <= 20 and nbCols <= 10) {
            
            vector<double> v(nbCols,0);
            
            cout << name << " = " << endl;
            
            for (int i(0); i<nbRows; ++i){
                
                v.assign(nbCols, 0);

                cout << "[ " ;
                for ( int j(Ai[i]); j<Ai[i+1]; ++j) {
                    v[Aj[j]] = A[j];
                }
                
                for (int j(0); j<nbCols; ++j) {
                    //~ cout <<   std::setprecision (std::numeric_limits<double>::digits10-6 ) << A[k] << " ";
                    //~ cout <<  fixed << A[k] << " "; // nombre de chiffre après la décimale fixee
                    //~ cout <<  scientific << A[k] << " "; // ecriture scientique
                    //~ cout <<  defaultfloat << A[k] << " "; // nombre max de chiffre fixee (avant et apres) 
                    cout <<  std::scientific << v[j] << "   ";
                }
                cout << "]" << endl;
                
            }
                       
        }        
        
        
        cout << endl;
        
    } else {
        cout << "Matrice " << name << " vide." << endl << endl;
    }
    

    
}


//***********************************************************************
ostream& CSRmatrix::affiche(ostream& out) const{
    

    if ( nbRows != 0 and nbCols != 0 and not Ai.empty() ) { 
        // affiche le nombre de chiffre, identique à std::setprecision(10)
        out.precision(10);        
        
        // Affichage la matrice plein si elle n'est pas trop grande, sinon en csr 
        if (nbRows <= 20 and nbCols <= 10) {
            
            vector<double> v(nbCols,0);            
            
            for (int i(0); i<nbRows; ++i){
                
                v.assign(nbCols, 0);
                for ( int j(Ai[i]); j<Ai[i+1]; ++j) {
                    v[Aj[j]] = A[j];
                }

                out << "[ " ;                
                for (int j(0); j<nbCols; ++j) {
                    out <<  std::scientific << v[j] << "   ";
                }
                out << "]" << endl;
                
            }
                       
        }  else {
            
            // Ai
            out << "Ai = [";
            for (int k(0); k<=nbRows; ++k){
                out <<  Ai[k] << " ";
            }
            out << "]" << endl;
            
            // Aj
            out << "Aj = [";
            for (int k(0); k<nnz; ++k){
                out <<  Aj[k] << " ";
            }
            out << "]" << endl;
            
            // A
            out << "A = [";
            for (int k(0); k<nnz; ++k){
                out <<  A[k] << " ";
            }
            out << "]" << endl;
                       
        }       
        
         

        
    } else {
        out << "Matrice vide." << endl << endl;
    }   
    
    return out;       
}


//***********************************************************************
void CSRmatrix::spy(string const& name) const {
    
    vector <int> tabi(nnz,0);
    
    // 1- Creation de Ai
    for (int i(0); i<nbRows; ++i) {
        for (int j(Ai[i]); j<Ai[i+1]; ++j){
            tabi[j] = i;
        }
    }
    
    // 2- Creation du fichier de data
    FILE * file;
    file = fopen ("spy.dat","w");
    
    if (file != NULL ) {
        for (int k(0); k<nnz; ++k){
            fprintf(file, "%d \t%d\n",Aj[k],tabi[k]);           
        }
        fclose(file);
    }   

    // 3- Suppression du tableau Ai
    tabi.clear();


    // 4- Affichage du spy    
    //~ plotSpyInFile( nnz, nbRows, nbCols, name); 
    plotSpyGUI( nnz, nbRows, nbCols, name); 
    

    // 5- Suppression du fichier de data
    if( remove( "spy.dat" ) != 0 ){
        perror( "Error deleting file" );    
    }
    
}


//***********************************************************************
CSRmatrix& CSRmatrix::generate(int const n, int const m, double const densite) {
    
    
    //~ if (sym) {
        //~ nnz = int( (densite * (n*m)) / 200); 
    //~ } else {
        nnz = int( (densite * (n*m)) / 100); 
    //~ }
    nbRows = n;
    nbCols = m;
    
    Ai.resize(nbRows+1);
    Aj.resize(nnz);
    A.resize(nnz);

                
    int l(0);

    // remplisage de Ai
    for (int i(0); i<nnz; ++i){
        
        // genere un nombre aleatoire entre 0 et nbRows tant que la ligne n'est pas remplie
        do {
            l=( rand() % ( (nbRows-1) - 0 + 1)) + 0;
        } while ( (Ai[l+1]-Ai[l])>=nbRows );
        
        // incremente le tableau Ai de l'indice l+1 à la fin 
        for (int j(l+1); j <= nbRows; ++j) {
            ++Ai[j];
        } 
                    
    }        
    
    // remplisage de Aj et A
    bool test(false);
    for (int i(0); i<nbRows; ++i){            
        for ( int j(Ai[i]); j<Ai[i+1]; ++j){
            
            // genere un nombre aleatoire entre 0 et nbColss tant que le nombre n'est pas deja apparue
            do {
                l=( rand() % ( (nbCols-1) - 0 + 1)) + 0;
                test = false;
                for (int k(j-1); k>=Ai[i]; --k){
                    if ( l == Aj[k]) {
                        test = true;
                    }
                }
            } while ( test);                
            Aj[j] = l;
            
            // genere un nombre entier aleatoire entre -10000 et 10000
            //~ A[j] = ( rand() % ( 10000 + 10000 + 1)) - 10000;
            // genere un nombre entier aleatoire entre -10 et 10
            A[j] = ( rand() % ( 10 + 10 + 1)) - 10;            
            // genere un nombre reel aleatoire entre -10000/1.333 et 10000/1.333
            //~ A[j] = (( rand() % ( 10000 + 10000 + 1)) - 10000)/1.333;
        }
        
        // tri la ligne des Aj dans l'ordre croissant
        tri_rapide(Aj, Ai[i], Ai[i+1]-1);
                
    }
        
    //~ if (sym) {        
    //~ int r1(0), r2(0);
    //~ bool test(false);
    //~ int cpt(0);
    //~ // remplisage de Ai, Aj et A
    //~ for (int i(0); i<nnz; ++i){            
        //~ do {   
            //~ // choisit aleatoirement une poisition (r1,r2)             
            //~ r1=( rand() % ( (nbRows-1) - 0 + 1)) + 0;
            //~ r2=( rand() % ( (nbRows-1) - 0 + 1)) + 0;
            //~ if ( r1<r2 ) {
                //~ int temp(r1);
                //~ r1 = r2;
                //~ r2 = temp;                
            //~ }
            
            //~ // verifie que la position n'existe pas
            //~ test = false;
            //~ for (int k(Ai[r1]); k<Ai[r1+1]; ++k){
                //~ if ( r2 == Aj[k]) {
                    //~ test = true;
                //~ }
            //~ }
        //~ } while ( test );  
             

        //~ // insere l'element dans le tableau Aj à la bonne place
        //~ if ( Ai[r1+1] < cpt ) {
            //~ for (int k(cpt); k>Ai[r1+1]; --k) {
                //~ Aj[k]= Aj[k-1];
            //~ }
            //~ Aj[Ai[r1+1]] = r2;
        
        //~ } else {
            //~ Aj[cpt] = r2;
        //~ }
        //~ ++cpt;
        


        //~ // incremente le tableau Ai de l'indice r1+1 à la fin 
        //~ for (int j(r1+1); j <= nbRows; ++j) {
            //~ ++Ai[j];
        //~ } 
             
                    
        //~ // genere un nombre aleatoire entre -10000 et 10000
        //~ A[i] = ( rand() % ( 10000 + 10000 + 1)) - 10000;
    //~ }        
    
    //~ // tri la ligne des Aj dans l'ordre croissant
    //~ for (int i(0); i<nbRows; ++i){                        
        //~ tri_rapide(Aj, Ai[i], Ai[i+1]-1);                    
    //~ }
        

    
    return *this;
}


//***********************************************************************
CSRmatrix CSRmatrix::transpose() const {
    
    CSRmatrix M(nbCols,nbRows,nnz);
    vector<int> AiTemp(nbCols+1,0);
    //~ vector<int> AiTemp2(Ai);
    //~ vector<int> AjTemp(nnz,0);
    //~ vector<double> ATemp(nnz,0);
    
    // Calcul le nb de nnz par ligne de la transposée,
    for ( int i(0); i < nnz; ++i) {
        ++AiTemp[ Aj[i]+1 ];
    }

    // Construction de Ai
    for ( int i(0); i < nbCols; ++i) {
        AiTemp[i+1] += AiTemp[i];
    }
    //~ Ai = AiTemp;
    M.Ai = AiTemp;


    // Construction de Aj, A
    int m(0), n(0);
    for ( int i(0); i < nbRows; ++i) {
        for ( int j(Ai[i]); j < Ai[i+1]; ++j) {
            
            m = Aj[j];
            n = AiTemp[m];
            
            //~ AjTemp[ n ] = i;
            //~ ATemp[ n ]  = A[j];
            M.Aj[ n ] = i;
            M.A[ n ]  = A[j];
            
            ++AiTemp[m];
        }
    }
    
    //~ int temp(nbRows);
    //~ nbRows = nbCols;
    //~ nbCols = temp;
    //~ A = ATemp;
    //~ Aj = AjTemp;
    
    
    return M;
}


//***********************************************************************
void CSRmatrix::save(string const& name) const{
    
    // ecriture via <stdio>, commande C, a priori plus rapide sur une majorite de compilo
    //~ FILE* file;
    //~ file = fopen(name.c_str(), "w");
    
    //~ fprintf(file, "%d\n", nbRows);
    //~ fprintf(file, "%d\n", nbCols);
    //~ fprintf(file, "%d\n", nnz);
    
    //~ for (int i(0); i<=nbRows;++i) {
        //~ fprintf(file, "%d\t", Ai[i]);
    //~ }
    //~ fprintf(file, "\n");
    
    //~ for (int i(0); i<nnz;++i) {
        //~ fprintf(file, "%d\t", Aj[i]);
    //~ }    
    //~ fprintf(file, "\n");
    
    //~ for (int i(0); i<nnz;++i) {
        //~ fprintf(file, "%lf\t", A[i]);
    //~ }
    //~ fclose(file);
    
    // ecriture via <fstream>, commande c++, plus safe que fopen
    ofstream file(name);
    file << nbRows << endl;
    file << nbCols << endl;
    file << nnz << endl;
    
   for (int i(0); i<=nbRows;++i) {
        file << Ai[i] << " ";
    }
    file <<endl;
    
    for (int i(0); i<nnz;++i) {
        file <<  Aj[i] << " ";
    }    
    file <<endl;
    
    file.precision(17);
    for (int i(0); i<nnz;++i) {
        file <<  A[i] << " ";
    }
    
    file.close();
    
}


//***********************************************************************
CSRmatrix& CSRmatrix::read(string const& name) {
    
    // lecture via <stdio>, commande C, a priori plus rapide sur une majorite de compilo
    //~ FILE* file;
    //~ file = fopen(name.c_str(), "r");
    
    //~ fscanf(file,"%d", &nbRows);
    //~ fscanf(file,"%d", &nbCols);
    //~ fscanf(file,"%d", &nnz);
    
    //~ Ai.resize(nbRows+1);
    //~ Aj.resize(nnz);
    //~ A.resize(nnz);
    
    //~ for (int i(0); i<=nbRows;++i) {
        //~ fscanf(file, "%d", &Ai[i]);
    //~ }
    
    //~ for (int i(0); i<nnz;++i) {
        //~ fscanf(file, "%d", &Aj[i]);
    //~ }    
    
    //~ for (int i(0); i<nnz;++i) {
        //~ fscanf(file, "%lf", &(A[i]));
    //~ }
    //~ fclose(file);
    
    
    // lecture via <fstream>, commande c++, plus safe que fopen
    ifstream file(name);
    
    file >> nbRows;
    file >> nbCols;
    file >> nnz;

    Ai.resize(nbRows+1);
    Aj.resize(nnz);
    A.resize(nnz);
    
   for (int i(0); i<=nbRows;++i) {
        file >> Ai[i];
    }
    
    for (int i(0); i<nnz;++i) {
        file >> Aj[i];
    }    
    
    for (int i(0); i<nnz;++i) {
        file >> A[i];
    }
    file.close();
    
    return *this;
}


//***********************************************************************
CSRmatrix& CSRmatrix::operator +=(CSRmatrix const& M) {
    
    vector<int> AiTemp(Ai);    
    vector<int> AjTemp(nnz+M.nnz);
    vector<double> ATemp(nnz+M.nnz);
    
    if ( nbRows != M.nbRows or nbCols != M.nbCols ) {
        cout << "Erreur dans l'addition, les dimensions ne correspondent pas." << endl;
        exit (EXIT_FAILURE);
    }
                  
    int cpt1(0), cpt2(0), cpt3(0);
    int nnz1(0), nnz2(0);
    int mj1(0), mj2(0);

    nnz = 0;
    Ai.assign(nbRows+1, 0);

    // Parcours les tailleSysteme des 2 matrices ligne après ligne
    for ( int i(0); i<nbRows; ++i) {

        nnz1 +=  AiTemp[i+1] - AiTemp[i];
        nnz2 +=  M.Ai[i+1] - M.Ai[i];
        Ai[i+1] = Ai[i];

        // tant que le nb de nnz de chaque ligne n'a pas été atteint faire
        while ( cpt3 < (nnz1 + nnz2) ) {

            if (cpt2 > nnz2) {
                mj2 = 0;                // n'importequelle valeur suffit
                mj1 =  Aj[cpt1];
            } else if (cpt1 > nnz1) {
                mj1 = 0;                // n'importequelle valeur suffit
                mj2 =  M.Aj[cpt2];
            } else { 
                mj1 =  Aj[cpt1];
                mj2 =  M.Aj[cpt2];
            }
            
            // CAS 1, nnz dans M1 et 0 dans M2
            if ( ( mj1 < mj2 and  cpt1 < nnz1 ) or cpt2 >= nnz2 ) {

                AjTemp[nnz] = Aj[cpt1];
                ATemp[nnz]  = A[cpt1];
                ++nnz;
                ++cpt1;
                ++cpt3;

            // CAS 2, 0 dans M1 et nnz dans M2
            } else if  ( ( mj1 >  mj2 and  cpt2 < nnz2  ) or cpt1 >= nnz1 ) {

                AjTemp[nnz] = M.Aj[cpt2];
                ATemp[nnz]  = M.A[cpt2];
                ++nnz;
                ++cpt2;
                ++cpt3;

            // CAS 3, nnz dans M1 ET nnz dans M2
            } else {

                AjTemp[nnz] = Aj[cpt1];
                ATemp[nnz]  = A[cpt1] + M.A[cpt2];
                ++nnz;
                ++cpt1;
                ++cpt2;
                cpt3 += 2;
            }
            
            // update Ai table
            ++Ai[i+1];
            
        }
    }
    
    Aj = { AjTemp.begin(), AjTemp.begin()+nnz};
    A  = { ATemp.begin(), ATemp.begin()+nnz};
        
    return *this;
}


//***********************************************************************
CSRmatrix& CSRmatrix::operator -=(CSRmatrix const& M) {
    
    vector<int> AiTemp(Ai);    
    vector<int> AjTemp(nnz+M.nnz);
    vector<double> ATemp(nnz+M.nnz);
    
    if ( nbRows != M.nbRows or nbCols != M.nbCols ) {
        cout << "Erreur dans la soustraction, les dimensions ne correspondent pas." << endl;
        exit (EXIT_FAILURE);
    }
                  
    int cpt1(0), cpt2(0), cpt3(0);
    int nnz1(0), nnz2(0);
    int mj1(0), mj2(0);

    nnz = 0;
    Ai.assign(nbRows+1, 0);

    // Parcours les tailleSysteme des 2 matrices ligne après ligne
    for ( int i(0); i<nbRows; ++i) {

        nnz1 +=  AiTemp[i+1] - AiTemp[i];
        nnz2 +=  M.Ai[i+1] - M.Ai[i];
        Ai[i+1] = Ai[i];

        // tant que le nb de nnz de chaque ligne n'a pas été atteint faire
        while ( cpt3 < (nnz1 + nnz2) ) {

            if (cpt2 > nnz2) {
                mj2 = 0;                // n'importequelle valeur suffit
                mj1 =  Aj[cpt1];
            } else if (cpt1 > nnz1) {
                mj1 = 0;                // n'importequelle valeur suffit
                mj2 =  M.Aj[cpt2];
            } else { 
                mj1 =  Aj[cpt1];
                mj2 =  M.Aj[cpt2];
            }
            
            // CAS 1, nnz dans M1 et 0 dans M2
            if ( ( mj1 < mj2 and  cpt1 < nnz1 ) or cpt2 >= nnz2 ) {

                AjTemp[nnz] = Aj[cpt1];
                ATemp[nnz]  = A[cpt1];
                ++nnz;
                ++cpt1;
                ++cpt3;

            // CAS 2, 0 dans M1 et nnz dans M2
            } else if  ( ( mj1 >  mj2 and  cpt2 < nnz2  ) or cpt1 >= nnz1 ) {

                AjTemp[nnz] = M.Aj[cpt2];
                ATemp[nnz]  = -M.A[cpt2];
                ++nnz;
                ++cpt2;
                ++cpt3;

            // CAS 3, nnz dans M1 ET nnz dans M2
            } else {

                AjTemp[nnz] = Aj[cpt1];
                ATemp[nnz]  = A[cpt1] - M.A[cpt2];
                ++nnz;
                ++cpt1;
                ++cpt2;
                cpt3 += 2;
            }
            
            // update Ai table
            ++Ai[i+1];
            
        }
    }
    
    Aj = { AjTemp.begin(), AjTemp.begin()+nnz};
    A  = { ATemp.begin(), ATemp.begin()+nnz};
        
    return *this;
}


//***********************************************************************
CSRmatrix& CSRmatrix::operator *=(CSRmatrix const& M) {

    vector<int> AiTemp(Ai);    
    int siz(nnz*M.nnz);
    if (siz>50000000) {
        siz = 50000000;
    }
    vector<int> AjTemp(siz);
    vector<double> ATemp(siz);
    vector<double> C(M.nbCols);

    if ( nbCols != M.nbRows ) {
        cout << "Erreur dans la multiplication, les dimensions ne correspondent pas." << endl;
        exit (EXIT_FAILURE);
    }
                  
    nnz = 0;
    Ai.assign(nbRows+1, 0);
    
    int row(0);
    
    for (int i(0); i<nbRows; ++i) {
        
        
        Ai[i+1] = Ai[i];
        C.assign(M.nbCols,0.0);
        
        // calcul de la ligne i de la matrice resultante
        for (int j(AiTemp[i]); j<AiTemp[i+1]; j++) {
        
            row = Aj[j];
            for (int k(M.Ai[row]); k<M.Ai[row+1]; ++k) {
                C[M.Aj[k]] += ( A[j] * M.A[k] );
            } 
        
        }
        
        
        // parcours de la ligne pour stocker les elements non nuls
        for (int j(0); j<M.nbCols; ++j) {
            if ( C[j] != 0.0 ) {                
                AjTemp[nnz] = j;
                ATemp[nnz] = C[j];
                ++nnz;
                ++Ai[i+1];
            }
        } 
        
    }
 
    nbCols = M.nbCols;
    Aj = { AjTemp.begin(), AjTemp.begin()+nnz};
    A  = { ATemp.begin(), ATemp.begin()+nnz};

    return *this;
}


//***********************************************************************
bool CSRmatrix::operator==(CSRmatrix const& M) const {
    
    // test rapide
    if ( nbRows != M.nbRows or nbCols != M.nbCols or nnz != M.nnz ) {
        return false;
    }
    
    if ( Ai != M.Ai  or Aj != M.Aj  ) { //or A != M.A
        return false;
    }
    
    double eps(1e-16);
    for (int i(0); i<nnz; ++i) {
        if ( fabs(A[i] - M.A[i]) > eps ) {
            return false;
        }
    }
        
    return true;
}


//***********************************************************************
vector<double> CSRmatrix::operator*(vector<double> const& v1) const {
    vector<double> v2(v1.size(),0.0);
    
    for (int i(0); i<nbRows; ++i) {
        
        for (int j(Ai[i]); j<Ai[i+1]; ++j) {
            
            v2[i] += ( A[j] * v1[Aj[j]] );
            
        }
    
    }
    
    return v2;
}


// ======================================================================
//                          METHODES EXTERNES 
// ======================================================================

const CSRmatrix operator+(CSRmatrix M1, CSRmatrix const& M2) {    
    return M1+=M2;     
}


//***********************************************************************
const CSRmatrix operator-(CSRmatrix M1, CSRmatrix const& M2) {    
    return M1-=M2;     
}


//***********************************************************************
const CSRmatrix operator*(CSRmatrix M1, CSRmatrix const& M2) {    
    return M1*=M2;     
}


//***********************************************************************
bool operator!=(CSRmatrix const& M1, CSRmatrix const& M2) {
    return !(M1 == M2);
}


//***********************************************************************
ostream& operator<<(ostream& out, CSRmatrix const& M) {
    M.affiche(out);
    return out;
}
