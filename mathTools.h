#include <iostream>
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <vector>
using namespace std;



void tri_rapide(vector<int>& tableau, int premier, int dernier) ;
int partitionner(vector<int>& tableau, int premier, int dernier, int pivot);
