#include <omp.h>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <random>


////////////////////////////////////////////////////////////////
//                      write to grid_file                    //
////////////////////////////////////////////////////////////////

void write_grid( std::ofstream& grid_file, int Ni, int Nj, std::vector < std:: vector < int >> grid){
    
    //write size of image
    grid_file << Ni << " " << Nj << std::endl;

    //write out the grid
    for (int i=0; i<Ni; ++i){
        for (int j=0; j<Nj; ++j){
            grid_file << i << " " << j << " " << grid[i][j] << std::endl;
        }
    }
}



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
        //fixed seed for reproducible behaviour
        unsigned int omp_seed = 42 + omp_get_thread_num();
        std::mt19937 gen(omp_seed);

        #pragma omp for schedule(dynamic)
        for (int i=0; i < Ni; ++i){
            for (int j=0; j < Nj; ++j){
                
                //generate rand floating number
                std::uniform_real_distribution<double> dis(0.0, 1.0);
                double rand = dis(gen);

                //if the random number is less than our probability p, fill site with a live tree
                if (rand < p){
                    forest[i][j] = 1;
                }
            }
        }
    }

    //set top row to burning tree i.e when at positions of j=0 cell state = 2
    for (int i=0; i<Ni; ++i){
        if (forest[0][i] == 1){
            forest[0][i] = 2;
        }
    }

    return forest;
}



//////////////////////////////////////////////////////////////////
//                        forest fire                           //
//////////////////////////////////////////////////////////////////

bool burn_forest(int Ni, int Nj, char wind, std::vector<std::vector<int>>& old_forest, std::vector<std::vector<int>>& new_forest)
{
    //********
    // Paramaters:
    // Inputs:
    //      int Ni, size of vector in x direction
    //      int Nj, size of vector in y direction
    //      bool wind, whether wind is part of the sim or not
    //      Vector of vectors old_forest, state of forest at timestep t-1
    //      Vector of vectors new_forest, state of forest at current timestep
    //
    // Outputs:
    //      Void output

    //assume forest is no longer burning unless told otherwise
    bool burning = false;


    #pragma omp parallel for schedule(dynamic) reduction(||:burning)
    for (int i=0; i<Ni; ++i){
        for (int j=0; j<Nj; ++j){

            //if tree is burning, otherwise ignore
            if (old_forest[i][j] == 2){

                //North
                if (i>=2 && old_forest[i-2][j] == 1 && wind == 'N'){
                    new_forest[i-2][j] = 2;
                }
                            
                if (i!=0 && old_forest[i-1][j] == 1){       
                    new_forest[i-1][j] = 2;
                }

                //South
                if (i<=Ni-3 && old_forest[i+2][j] == 1 && wind == 'S'){
                    new_forest[i+2][j] = 2;
                }

                if (i<=Ni-2 && old_forest[i+1][j] == 1) {
                    new_forest[i+1][j] = 2;
                }

                //West
                if (j>=2 && old_forest[i][j-2] == 1 && wind == 'W'){
                    new_forest[i][j-2] = 2;
                }

                if (j>=1 && old_forest[i][j-1] == 1) {
                    new_forest[i][j-1] = 2;
                }

                //East
                if (j<=Nj-3 && old_forest[i][j+2] == 1 && wind == 'E'){
                    new_forest[i][j+2] = 2;
                }

                if (j <= Nj-2 && old_forest[i][j+1] == 1){
                    new_forest[i][j+1] = 2;
                }

                //The burning tree dies
                new_forest[i][j] = 3;

                //we had at least one burning tree
                #pragma omp critical  //avoid race conditions
                burning = true;
                

            }

        }
    }
    
    //all other trees remain in the same state so do nothing else 

    return burning;
}



/////////////////////////////////////////////////////////////////////
//                        get num of trees                         //
/////////////////////////////////////////////////////////////////////

int get_trees(int Ni, int Nj, std::vector < std::vector < int >> grid){
    //counter
    int num_trees = 0;
    
    #pragma omp parallel for reduction (+:num_trees)
     
     //loop through grid
     for (int i=0; i<Ni; ++i){
        for (int j=0; j<Nj; ++j){
            if (grid[i][j]>0){
                ++num_trees;
            }
        }
     }

    return num_trees;
}



/////////////////////////////////////////////////////////////////////
//                             main                                //
/////////////////////////////////////////////////////////////////////

int main(int argc, char **argv){

    //read in command line variables and initialise forest
    int Ni = atoi(argv[1]);  // the size of the grid in the y-dimension
    int Nj = atoi(argv[2]);  // the size of the grid in the x-dimension
    double p = atof(argv[3]); // probabilty p of fire spreading
    char wind = *(argv[4]); // direction of wind, any character but NSEW if no wind wanted

    std::string write_str = argv[5]; //whether to write out the grid
    bool write = (write_str == "true");


    std::vector < std::vector < int > > old_forest = init_forest(Ni, Nj, p);
    std::vector < std::vector < int > > new_forest = old_forest;

    //start time
    double start = omp_get_wtime();

    //open the grid file
    std::ofstream grid_file("forest_wind.dat");

    //write the initial grid
    if (write){
        write_grid(grid_file, Ni, Nj, old_forest);
    }

    //run the programme - while the forest is burning - then end
    bool burning = true;
    int nsteps = 0;

    while (burning){

        burning = burn_forest(Ni, Nj, wind, old_forest, new_forest);
        # pragma omp barrier

        // to continue loop, old forest is reassigned to new forest
        old_forest = new_forest;

        if (burning){
            nsteps++;
            if (write){
                write_grid(grid_file, Ni, Nj, new_forest);
            }
        }

    }

    bool bottom = false;
    // see if the fire reached the bottom of the grid after sim
    for (int j = 0; j < Nj; ++j){
        if (new_forest[Ni-1][j] == 3){
            bottom = true;
        }
    }

    //end timer 
    double end = omp_get_wtime();

    std::cout << omp_get_max_threads() << "," << nsteps << "," << bottom << "," << end - start << std::endl;

    return EXIT_SUCCESS;

}
