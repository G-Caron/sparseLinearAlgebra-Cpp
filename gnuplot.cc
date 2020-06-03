#include "gnuplot.h"


gnuplot::gnuplot() {
    gnuplotpipe = popen("gnuplot -persist","w");
    if (!gnuplotpipe) {
        cerr << "gnuplot not find!" << endl;
    }    
}

gnuplot::~gnuplot() {
    fprintf(gnuplotpipe, "exit\n");
    pclose(gnuplotpipe);
    //~ delete gnuplotpipe;
}

void gnuplot::operator() (const string& command) {
    
    fprintf(gnuplotpipe,"%s\n",command.c_str());
    fflush(gnuplotpipe);    
}


void plotSpyInFile(int nnz, int nbRows, int nbCols, string name) {
    
    string message = "";
    gnuplot g;

    message = "set xlabel \"nnz = " + to_string(nnz) + "\" ";
    g(message);

    message = "set xrange [" + to_string(-1) + ":" + to_string(nbCols) + "]";
    g(message);

    message = "set yrange [" + to_string(nbRows) + ":" + to_string(-1) + "] reverse";
    g(message);

    g("set size ratio -1");
    
    g("set colorsequence classic");

    g("set term postscript enhanced eps color font \"Palatino,15\" ");

    message = "set output \"spy_" + name + ".ps" + "\" ";
    g(message);

    //~ if (symmetric) {
        //~ message = "set title \"Sparsity pattern of the symmetric matrix " + name + "\" font \"Palatino-Bold,17\"";
    //~ } else {
        message = "set title \"Sparsity pattern of the matrix " + name + "\" font \"Palatino-Bold,17\"";        
    //~ }
    g(message);

    //~ if (symmetric) {
        //~ message = "set label \"Only the lower triangular half\" at " + to_string(40*nbCols/100) + "," + to_string(10*nbCols/100);
        //~ g(message);
        //~ message = "set label \"of the matrix is stored\" at " + to_string(50*nbCols/100) + "," + to_string(20*nbCols/100);
        //~ g(message);
    //~ }

    g("plot \'spy.dat\' u 1:2 w p lc 1 notitle");
    
    //~ g.~gnuplot();
}


void plotSpyGUI(int nnz, int nbRows, int nbCols, string name) {
    
    string message = "";
    gnuplot g;

    message = "set xlabel \"nnz = " + to_string(nnz) + "\" ";
    g(message);

    message = "set xrange [" + to_string(-1) + ":" + to_string(nbCols) + "]";
    g(message);

    message = "set yrange [" + to_string(nbRows) + ":" + to_string(-1) + "] reverse";
    g(message);

    g("set size ratio -1");
    
    g("set colorsequence classic");


    //~ if (symmetric) {
        //~ message = "set title \"Sparsity pattern of the symmetric matrix " + name + "\" font \"Palatino-Bold,17\"";
    //~ } else {
        message = "set title \"Sparsity pattern of the matrix " + name + "\" font \"Palatino-Bold,17\"";        
    //~ }
    g(message);

    //~ if (symmetric) {
        //~ message = "set label \"Only the lower triangular half\" at " + to_string(40*nbCols/100) + "," + to_string(10*nbCols/100);
        //~ g(message);
        //~ message = "set label \"of the matrix is stored\" at " + to_string(50*nbCols/100) + "," + to_string(20*nbCols/100);
        //~ g(message);
    //~ }
    
    g("plot \'spy.dat\' u 1:2 w p lc 1 notitle");
    
    //~ g.~gnuplot();
}
