#include "mathTools.h"


void tri_rapide(vector<int>& tableau, int premier, int dernier) {
    int pivot;
    if (premier<dernier) {
        pivot = (rand() % (dernier - premier + 1)) + premier;
        pivot = partitionner(tableau, premier, dernier, pivot);
        tri_rapide(tableau, premier, pivot-1);
        tri_rapide(tableau, pivot+1, dernier);
    }
}

int partitionner(vector<int>& tableau, int premier, int dernier, int pivot){
    
   // échange le pivot avec le dernier du tableau , le pivot devient le dernier du tableau
    int temp = tableau[dernier];
    tableau[dernier] = tableau[pivot];
    tableau[pivot] = temp;
    
    
    int j = premier;
    for (int i=premier; i<dernier; i++){
        if ( tableau[i] <= tableau[dernier] ) {
            temp = tableau[i];
            tableau[i] = tableau[j];
            tableau[j] = temp;
            j++;
        }
    }
    temp = tableau[dernier];
    tableau[dernier] = tableau[j];
    tableau[j] = temp;
    return j;
}
