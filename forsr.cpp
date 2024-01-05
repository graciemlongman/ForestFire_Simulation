#include <omp.h>
#include <time.h> 
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <cmath>
#include <random>

////////////////////////////////////////////////////////////////
//           initialise and fill random grid with trees       //
//                      trees in a circle?                    //
////////////////////////////////////////////////////////////////

std::vector < std::vector < int > > init_forest(int Ni, int Nj, double p)
{
    //********
    // Paramaters:
    //
    // Inputs:
    //      int Ni, size of vector in x direction
    //      int Nj, size of vector in y direction
    //      double p, probability p of tree
    //
    // Outputs:
    //      Vector of Vectors, with top row of trees set to burning
    //********

    //init grid called forest
    std::vector<std::vector<int> > forest(Ni, std::vector<int>(Nj,0));

    //randomly fill grid
    #pragma omp parallel
    {
        //fixed seed for reproducible behvaiour
        int seed = 42 + omp_get_thread_num();

        #pragma omp for
        for (int i=0; i < Ni; ++i){
            for (int j=0; j < Nj; ++j){
                
                //generate rand floating number
                std::mt19937 gen(seed);
                std::uniform_real_distribution<double> dis(0.0, 1.0);
                double rand = dis(gen);

                //if the random number is less than our probability p, fill site with a live tree
                if (rand / RAND_MAX < p){
                    forest[i][j] = 1;
                }
            }
        }
    }

    //set top row to burning tree i.e when at positions of j=0 cell state = 2
    for (int i=0; i<Ni; ++i){
        if (forest[i][0] = 1){
            forest[i][0] =2;
        }
    }

    return forest;
}


//////////////////////////////////////////////////////////////////
//                        forest fire                           //
//////////////////////////////////////////////////////////////////

bool burn_forest(int Ni, int Nj, std::vector<std::vector<int>>& old_forest, std::vector<std::vector<int>>& new_forest)
{
    //********
    // Paramaters:
    //scheduling?
    // Inputs:
    //      int Ni, size of vector in x direction
    //      int Nj, size of vector in y direction
    //      Vector of vectors old_forest, state of forest at timestep t-1
    //      Vector of vectors new_forest, state of forest at current timestep
    //
    // Outputs:
    //      Void output

    //assume forest is no longer burning unless told otherwise
    bool burning = false;

    #pragma omp parallel for

    //loop through the forest
    for (int i=0; i<Ni; ++i){
        for (int j=0; j<Nj; ++j){

            //living trees next to the burning tree catches fire
            if (old_forest[i][j] == 2){

                //North
                if (i!= 0){
                    if (old_forest[i-1][j] == 1){
                            new_forest[i-1][j] == 2;
                    }
                }

                //South
                if (i!=Ni-1){
                    if (old_forest[i+1][j] == 1){
                        new_forest[i+1][j] == 2;
                    }
                }

                //West
                if (j!=0){
                    if(old_forest[i][j-1] == 1){
                        new_forest[i][j-1] == 2;
                    }
                }

                //East
                if (j!=Nj-1){
                    if (old_forest[i][j+1] == 1){
                        new_forest[i][j+1] == 2;
                    }
                }

                //The burning tree dies
                new_forest[i][j] == 3;

                //we had at least one burning tree
                burning = true;

            }

        }
    }
    
    //all other trees remain in the same state so do nothing else 

    return burning;
}

int main(int argc, char **argv){

    //read in command line variables and initialise forest
    int Ni = atoi(argv[1]);  // the size of the grid in the x-dimension
    int Nj = atoi(argv[2]);  // the size of the grid in the y-dimension
    double p = atof(argv[3]); //probabilty p of fire spreading

    std::vector < std::vector < int > > old_forest = init_forest(Ni, Nj, p);
    std::vector < std::vector < int > > new_forest(&Ni, &Nj);

    //start time
    double start = omp_get_wtime();

    //run the programme - while the forest is burning - then end
    bool burning = true;
    int nsteps = 0;

    while (burning){

        burning = burn_forest(Ni, Nj, old_forest, new_forest);

        if (burning){
            nsteps++;
        }

    }

    //end timer 
    double end = omp_get_wtime();

    std::cout << "Forest burnt out in " << nsteps << " iterations" << std::endl;
    std::cout << "Total time taken = "  << end - start << " s" << std::endl;

    return EXIT_SUCCESS;

}