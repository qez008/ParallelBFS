/******************************************************************************
 *
 *    This program contains necessary support for running sequential and 
 *    parallel graph algorithms. 
 *
 *    The outline of the program is as follows:
 *    Information about the graphs to use and the setup for the experiments 
 *    is first read from the file given as argument (See "data_file").
 *
 *    Next each graph is read in from an external file and processed. For each 
 *    graph the program sets up the necessary data structures and then calls 
 *    the different graph algorithms using this graph. Sequential algorithm
 *    are executed first and then parallel ones. For each algorithm the program 
 *    performs timings and verifies the result. These results are first printed 
 *    to the screen and also stored to a file.
 *
 *******************************************************************************/

#include "driver.h"	    // Include system files and define variables
#include "cFiles.h"     // Include c-files used for different algorithms

int main(int argc, char *argv[]) {

    /* The following code is used for binding the threads to different cores */
    if (numa_available() < 0) {
        printf("No NUMA support available on this system.\n");
        exit(1);
    }
    numa_set_interleave_mask(numa_all_nodes_ptr);

    // Primary graph data structure, using compressed neighbor lists.
    // ver stores pointers into the edges list. Note that numbering of both vertices and edges starts from 1

    int *ver;          // Pointers to edge lists
    int *edges;        // Edge lists
    double *weight;    // Corresponding edge weights

    // Pointers to 5 lists, each of length n=|V|. These are used in the graph algorithms
    int *p,*p1,*p2,*p3,*p4;

    // List of Edges in the order they were read from file
    // The edge data structure is defined in driver.h

    edge *e;           // List of all edges, no weight, one struct for each edge

    // OpenMP locks, not needed in INF236
    omp_lock_t *nlocks;
    nlocks = (omp_lock_t *) malloc((max_n)*sizeof(omp_lock_t));

    // Routines for performing IO and for testing solutions. The code is 
    // given at the end of the file
    //
    int read_graph(),get_input(),allocate_memory(),read_general_graph();
    int read_bin_graph(),nrComponents();
    void store_value(),store_p_value();

    int i,j;          // Loop indices
    double mt1,mt2;   // Used for timing of each algorithm

    // Boolean variables for specifying which algorithms to run 
    // Currently there are three available algorithms.
    // Only the algorithms where the boolean value is set to true are run
    //
    int BFS               = true;    /* Sequential BFS */
    int PBFS              = false;   /* Parallel BFS */
    int ABFS              = true;    /* Alternative Parallel BFS */

    // Read input and check that it is in order
    if (!get_input(argc,argv,&n_graphs,&n_runs,&n_conf,&conf,&name))	
        return(false);

    // Prepare output file
    if (!prepare_output(&wf,n_conf,conf)) {				
        return(false);
    }

    // Allocate memory
    if (!allocate_memory(n_graphs,&n,&m,&timer,&cost,&p_timer,&nrC,n_conf)) {		
        return(false);
    }
    int *c;
    c  = (int *) malloc(sizeof(int)*100);

    // Allocating memory for sequential graph algorithms
    if (!allocate_graph_memory(&weight,&e,&ver,&edges,max_n,&p,&p1,&p2,&p3,&p4)) {			
        return(false);
    }

    // Run the main loop, processing one graph at a time, first sequentially and then in parallel
    //
    for(i=0;i<n_graphs;i++) {

        printf("\n");
        printf("************************************************\n");
        printf("*  Graph %2d: %20s             *\n",i,name[i]);


        // Reading the i'th graph

        //
        if (!read_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e,&weight,p,p1,p2)) {
            printf("Problem reading graph %s \n",name[i]);
            return(false);
        }


        /*

           if (!read_general_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e,&weight,&we,&swe,p,p1)) {
           printf("Problem reading graph %s \n",name[i]);
           return(false);
           }


           if (!read_bin_graph(&(n[i]),&(m[i]),name[i],&ver,&edges,&e,&weight,p1,p2)) {
           printf("Problem reading graph %s \n",name[i]);
           return(false);
           }
           */


        printf("*  |V| = %9d, |E| = %10d           *\n",n[i],m[i]);
        printf("************************************************\n");

        //    printf("Done reading graph \n");


        // Run sequential code 

        printf("\n************************************************\n");
        printf("***           Sequential algorithms          ***\n");
        printf("************************************************\n");

        int l,nrCS;


        // Compute BFS starting from vertex 1
        // Upon return p will contain parent pointers and p1 shortest distance from vertex 1

        if (BFS) {
            for(l=0;l<n_runs;l++) {
                mt1 = omp_get_wtime();      // Get current time
                sbfs(n[i],ver,edges,p,p1,p2,p3);

                mt2 = omp_get_wtime();      // Get current time
                if ((timer[0][i] < 0) || (mt2-mt1 < timer[0][i])) // Capture best run
                    timer[0][i] = mt2-mt1;
            }


            printf("************************************************\n");
            printf("*  Sequential BFS                              *\n");
            printf("*  Time = %8.6lf seconds                     *\n",timer[0][i]);
            printf("************************************************\n");
        }




        // Allocate space for locks

        //    omp_lock_t *nlocks;
        //    nlocks = (omp_lock_t *) malloc((n[i]+1)*sizeof(omp_lock_t));

        printf("************************************************\n");
        printf("***            Parallel algorithms           ***\n");
        printf("************************************************\n");

        // Run parallel algorithm
        // One iteration for each thread configuration
        for(j=0;j<n_conf;j++) {
            printf("************************************************\n");
            printf("*       Configuration %2d, using %2d threads     *\n",j,conf[j]);
            printf("************************************************\n");

            omp_set_num_threads(conf[j]);  // Set number of threads in this configuration

#pragma omp parallel 
            {
                int k;
                int threads = omp_get_num_threads();
                int my_id = omp_get_thread_num();


                if (threads != conf[j]) {
                    printf("***** Did not get the expected numer of threads! Wanted %d, got %d *****\n",conf[i],threads);
                }

                // Now we are ready to run the parallel algorithms
                // Run as many times as required, only store timings for best run

#pragma omp barrier
                for(k=0;k<n_runs;k++) {   // Run all experiments n_runs times, select best run time as final

                    if (PBFS) {   // Parallel BFS, distance values will be returned in p2
#pragma omp barrier
#pragma omp master
                        { mt1 = omp_get_wtime(); }
                        pbfs(n[i],ver,edges,p,p2,p3,p4);

#pragma omp master
                        { mt2 = omp_get_wtime();
                            if ((p_timer[0][i][j] < 0) || (mt2-mt1 < p_timer[0][i][j]))
                                p_timer[0][i][j] = mt2-mt1;

                        } // master
                    }
                } // End of k iterations over the same configuration

                // Verifying solution, check that distance values are the same as the sequential bfs
#pragma omp master
                if (PBFS)
                    verify_bfs(n[i],p1,p2);

#pragma omp barrier

                // Print results to screen

#pragma omp master
                {
                    if (PBFS) {
                        printf("************************************************\n");
                        printf("*  Parallel BFS                                *\n");
                        printf("*  Time = %8.6lf seconds                     *\n",p_timer[0][i][j]);
                        printf("*  Speedup = %4.2lf, Relative speedup = %4.2lf     *\n",timer[0][i]/p_timer[0][i][j],p_timer[0][i][0]/p_timer[0][i][j]);
                        printf("************************************************\n");
                    }
                } // master


#pragma omp barrier
                for(k=0;k<n_runs;k++) {   // Run all experiments n_runs times, select best run time as final

                    if (ABFS) {   // Alternative parallel BFS
#pragma omp barrier
#pragma omp master
                        { mt1 = omp_get_wtime(); }
                        abfs(n[i],ver,edges,p,p2,p3,p4);

#pragma omp master
                        { mt2 = omp_get_wtime();
                            if ((p_timer[1][i][j] < 0) || (mt2-mt1 < p_timer[1][i][j]))
                                p_timer[1][i][j] = mt2-mt1;

                        } // master
                    }
                } // End of k iterations over the same configuration

                // Verifying solution
#pragma omp master
                if (ABFS)   
                    verify_bfs(n[i],p1,p2);

#pragma omp barrier

                // Print results to screen

#pragma omp master
                {
                    if (ABFS) {
                        printf("************************************************\n");
                        printf("*  Alternative parallel BFS                    *\n");
                        printf("*  Time = %8.6lf seconds                     *\n",p_timer[1][i][j]);
                        printf("*  Speedup = %4.2lf, Relative speedup = %4.2lf     *\n",timer[0][i]/p_timer[1][i][j],p_timer[1][i][0]/p_timer[1][i][j]);
                        printf("************************************************\n");
                    }
                } // master

#pragma omp barrier

            }  // End of parallel section

        }  // End of different thread configurations

    } // Loop over graphs


    // Print results to file


    char *tc;         // Name of current graph without the .mtx extension

    fprintf(wf,"name = {");       // Print the names of the graphs
    for(i=0;i<n_graphs;i++)  {
        tc = strtok(name[i],"."); 		// Remove .mtx extension from name
        fprintf(wf,"'%s' ",tc);
        if (i != n_graphs-1)
            fprintf(wf,",");
    }
    fprintf(wf,"};\n");

    fprintf(wf,"n = [");       // Print the number of vertices of the graphs
    for(i=0;i<n_graphs;i++)  {
        fprintf(wf,"%d ",n[i]);
    }
    fprintf(wf,"];\n");

    fprintf(wf,"nz = [");       // Print the number of non-zeros of the graphs
    for(i=0;i<n_graphs;i++)  {
        fprintf(wf,"%d ",m[i]);
    }
    fprintf(wf,"];\n");

    // Sequential timings
    // ******************
    // timer[0][]	Sequential BFS
    //
    // Parallel timings, all timings gives the best run
    // ****************
    // p_timer[0][i][conf]	Parallel BFS
    //

    if (BFS) 
        store_value(wf,"TimeSequentialBFS",timer[0],n_graphs);    // Print timing results for sequential BFS

    // Print out parallel values

    if (PBFS) 
        store_p_value(wf,"TimeParallelBFS",p_timer[0],n_graphs,n_conf); // Print timing results for parallel BFS

    if (ABFS) 
        store_p_value(wf,"TimeAParallelBFS",p_timer[1],n_graphs,n_conf); // Print timing results for alt parallel BFS

    fclose(wf);

    printf("************************************************\n");
    printf("***               Normal exit                ***\n");
    printf("************************************************\n");
}

int read_bin_graph(int *n,int *m,char *f_name,int **ver,int **edges,edge **e,double **weight,int *count,int *start) {

    int i;
    int num_edges;
    FILE *fp;
    int x,y;

    fp = fopen(f_name,"r");
    if (fp == NULL) {
        printf("Could not open file %s \n",f_name);
        return(false);
    }

    fread(n,sizeof(int),1,fp);
    fread(m,sizeof(int),1,fp);
    printf("n=%d m=%d \n",*n,*m);


    // Start counting of degrees by setting counters to zero
    for(i=1;i<=*n;i++) {
        count[i] = 0;
    }

    fread(*e,(*m)*(2*sizeof(int)),1,fp);
    num_edges = 0;

    for(i=0;i<*m;i++) {

        x = (*e)[num_edges].x;
        y = (*e)[num_edges].y;

        count[x]++; // Get vertex degrees
        count[y]++;

        num_edges++;
        if (num_edges > *m) {
            printf("Have set num_edges to %d while i=%d \n",num_edges,i);
            return;
        }
    }

    *m = num_edges;  // Make sure m is the correct number of edges

    // Find starting positions in edge list for each vertex
    start[1] = 0;
    (*ver)[1] = 0;
    for(i=2;i<=*n+1;i++) {
        start[i] = start[i-1]+count[i-1];
        (*ver)[i] = start[i];
    }

    // Place edges in edge lists, once for each endpoint

    for(i=0;i<*m;i++) {
        x = (*e)[i].x;
        y = (*e)[i].y;

        if ((x == 0) || (y == 0)) {
            printf("edge %d: %d and %d are neighbors, numbering starts from 1! \n",i,x,y);
            return(false);
        }

        (*edges)[start[x]] = y;
        //(*weight)[start[x]] = v;
        start[x]++;

        (*edges)[start[y]] = x;
        // (*weight)[start[y]] = v;
        start[y]++;
    }

    fclose(fp);
    //  free(start);
    //  free(count);

    return(true);
}

int read_graph(int *n,int *m,char *f_name,int **ver,int **edges,edge **e,
        double **weight,int *count, int *start,int *where) {

    // Note reuse of p as count and p1 as start

    FILE *fp;
    long int i;
    int x,y,z;
    double v,w;
    MM_typecode matcode;
    int cols;

    int mxdeg = 0;
    //int *count;
    //int *start;

    // printf("Opening file %s \n",f_name);

    // First make sure graph is of the right type

    fp = fopen(f_name,"r");
    if (fp == NULL) {
        printf("Could not open file %s \n",f_name);
        return(false);
    }

    if (mm_read_banner(fp, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return(false);
    }

    if ((mm_read_mtx_crd_size(fp, &cols, n, m)) !=0) {
        printf("Could not read size of graph.\n");
        return(false);
    }

    if (*m > max_m) {
        printf("Graph is too large. Asking for %d edges, max is %d \n",*m,max_m);
        return(false);
    }
    if (*n > max_n) {
        printf("Graph is too large. Asking for %d vertices, max is %d \n",*n,max_n);
        return(false);
    }

    if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) || !mm_is_sparse(matcode) || !mm_is_symmetric(matcode)) {
        printf("The program can only read files that are sparse symmmetric matrices in coordinate format! \n");
        return(false);
    }


    // Start counting of degrees by setting counters to zero
    for(i=1;i<=*n;i++) {
        count[i] = 0;
    }

    int num_edges = 0;

    // printf("Number of possible edges is %d \n",*m);

    // Read inn the edges
    if (mm_is_real(matcode)) {
        printf("Real matrix \n");
        //printf("Real matrix, Starting to read %d edges \n",*m);
        srand48(time(NULL));
        for(i=0;i<*m;i++) {
            // if ((i % 1000000) == 0)
            // printf("Have read %d edges \n",i);
            fscanf(fp,"%d %d %lf",&x,&y,&v); // Use this line if there is exactly one double weight 
            //  fscanf(fp,"%d %d %lf %lf",&x,&y,&v,&w); // Use this line for complex weights
            if (x != y) { // Avoid self-edges

                (*e)[num_edges].x = x; // Store edges
                (*e)[num_edges].y = y;

                //        int intv = (int) fabs(v);
                //        (*we)[num_edges].w = (double) intv; 

                count[x]++; // Get vertex degrees
                count[y]++;

                num_edges++;
                if (num_edges > *m) {
                    printf("Have set num_edges to %d while i=%d \n",num_edges,i);
                    return;
                }
            }
        }
    }
    else if (mm_is_integer(matcode)) {
        printf("Integer matrix, Starting to read %d edges \n",*m);
        srand48(time(NULL));
        for(i=0;i<*m;i++) {
            fscanf(fp,"%d %d %lf",&x,&y,&v); // Use this line if there is exactly one double weight 

            if (x != y) { // Avoid self-edges

                (*e)[num_edges].x = x; // Store edges
                (*e)[num_edges].y = y;

                count[x]++; // Get vertex degrees
                count[y]++;

                num_edges++;
                if (num_edges > *m) {
                    printf("Have set num_edges to %d while i=%d \n",num_edges,i);
                    return;
                }
            }
        }
    }
    else
    {          // Symbolic matrix
        srand48(time(NULL));
        //printf("Symbolic matrix \n");
        for(i=0;i<*m;i++) {
            fscanf(fp,"%d %d",&x,&y);
            if (x != y) { // Avoid self-edges

                (*e)[num_edges].x = x; // Store edges
                (*e)[num_edges].y = y;

                count[x]++; // Get vertex degrees
                count[y]++;

                num_edges++;
            }
        }
    }
    //printf("original edges %d, used edges %d \n",*m,num_edges);
    printf("*  |V| = %d, |E| = %d \n",*n,num_edges);
    printf("************************************************\n");

    *m = num_edges;  // Make sure m is the correct number of edges

    // Find starting positions in edge list for each vertex
    start[1] = 0;
    (*ver)[1] = 0;
    for(i=2;i<=((*n)+1);i++) {
        start[i] = start[i-1]+count[i-1];
        (*ver)[i] = start[i];
    }

    // Place edges in edge lists, once for each endpoint

    for(i=0;i<*m;i++) {
        x = (*e)[i].x;
        y = (*e)[i].y;
        //v = (*e)[i].w;

        if ((x == 0) || (y == 0)) {
            printf("edge %d: %d and %d are neighbors, numbering starts from 1! \n",i,x,y);
            return(false);
        }

        // where[start[x]] = start[y];
        // where[start[y]] = start[x];

        (*edges)[start[x]] = y;
        //(*weight)[start[x]] = v;
        start[x]++;

        (*edges)[start[y]] = x;
        //(*weight)[start[y]] = v;
        start[y]++;

    }

    int maxnode;

    for(i=1;i<=*n;i++) {
        if (((*ver)[i+1] - (*ver)[i]) > mxdeg) {
            mxdeg = (*ver)[i+1] - (*ver)[i];
            maxnode = i;
        }
    }
    /*
       printf("Vertex %d has maximum degree %d \n",maxnode,mxdeg);
       printf("Average degree is %d \n",((*m) * 2)/(*n));
       */


    //  printf("Done read_graph, freeing up memory \n");
    fclose(fp);
    // free(count);
    // free(start);
    return(true);
}

int read_general_graph(int *n,int *m,char *f_name,int **ver,int **edges,edge **e,
        double **weight,int *count, int *start) {

    // Note reuse of p as count and p1 as start

    FILE *fp;
    long int i;
    int x,y,z;
    double v,w;
    MM_typecode matcode;
    int cols;

    //int *count;
    //int *start;

    // printf("Opening file %s \n",f_name);

    // First make sure graph is of the right type

    fp = fopen(f_name,"r");
    if (fp == NULL) {
        printf("Could not open file %s \n",f_name);
        return(false);
    }

    if (mm_read_banner(fp, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        return(false);
    }

    if ((mm_read_mtx_crd_size(fp, &cols, n, m)) !=0) {
        printf("Could not read size of graph.\n");
        return(false);
    }

    // If this is an unsymmetric matrix then we only keep the square part.
    if (cols < *n)
        *n = cols;

    /*
       if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) || !mm_is_sparse(matcode) || !mm_is_symmetric(matcode)) {
       printf("The program can only read files that are sparse symmmetric matrices in coordinate format! \n");
       return(false);
       }
       */


    // Start counting of degrees by setting counters to zero
    for(i=1;i<=*n;i++) {
        count[i] = 0;
    }

    int num_edges = 0;

    //  printf("Number of possible edges is %d \n",*m);

    // Read inn the edges
    if (mm_is_real(matcode)) {
        //    printf("Starting to read %d edges \n",*m);
        for(i=0;i<*m;i++) {
            fscanf(fp,"%d %d %lf",&x,&y,&v); // Use this line if there is exactly one double weight 
            // fscanf(fp,"%d %d %lf %lf",&x,&y,&v,&w); // Use this line for complex weights
            if ((x < y) && (x <= *n) && (y <= *n)) { // Only use edges in square part
                (*e)[num_edges].x = x; // Store edges
                (*e)[num_edges].y = y;

                count[x]++; // Get vertex degrees
                count[y]++;

                num_edges++;
                if (num_edges > *m) {
                    printf("Have set num_edges to %d while i=%d \n",num_edges,i);
                    return;
                }
            }
        }
    }
    else {          // Symbolic matrix
        printf("Trouble ahead, the code now assumes weighted graphs \n");
        for(i=0;i<*m;i++) {
            fscanf(fp,"%d %d",&x,&y);
            if (x != y) { // Avoid self-edges
                (*e)[num_edges].x = x; // Store edges
                (*e)[num_edges].y = y;

                count[x]++; // Get vertex degrees
                count[y]++;

                num_edges++;
            }
        }
    }
    //   printf("original edges %d, used edges %d \n",*m,num_edges);
    printf("*  |V| = %d, |E| = %d \n",*n,num_edges);
    printf("************************************************\n");

    *m = num_edges;  // Make sure m is the correct number of edges

    // Find starting positions in edge list for each vertex
    start[1] = 0;
    (*ver)[1] = 0;
    for(i=2;i<=*n+1;i++) {
        start[i] = start[i-1]+count[i-1];
        (*ver)[i] = start[i];
    }

    // Place edges in edge lists, once for each endpoint

    for(i=0;i<*m;i++) {
        x = (*e)[i].x;
        y = (*e)[i].y;
        if ((x == 0) || (y == 0)) {
            printf("edge %d: %d and %d are neighbors, numbering starts from 1! \n",i,x,y);
            return(false);
        }

        (*edges)[start[x]] = y;
        //(*weight)[start[x]] = v;
        start[x]++;

        (*edges)[start[y]] = x;
        //(*weight)[start[y]] = v;
        start[y]++;

    }
    //  printf("Done read_graph, freeing up memory \n");
    fclose(fp);
    // free(count);
    // free(start);
    return(true);
}

// Get input. This is given as a file name to the executable. 
// This file should contain the following information (per line):
// Number of input graphs
// Number of times each configuration is run, the program reports the best running time out of these
// Number of configurations, followed by the number of threads in each configuration
// One line giving the complete name of each graph file that is to be run

int get_input(int argc, char *argv[],int *n_graphs,int *n_runs,int *n_conf,int **conf,char ***name) {

    FILE *rf;	// File pointer
    int i;

    if (argc != 2) {
        printf("Give data file name as first input parameter!\n");
        return(false);
    }

    printf("************************************************\n");
    printf("*       Reading setup from %10s          *\n",argv[1]);
    printf("************************************************\n");

    // Opening file containing file names for graphs

    rf = fopen(argv[1],"r");
    if (rf == NULL) {
        printf("Cannot open data file: %s \n",argv[1]);
        return(false);
    }


    // Get number of graphs to read
    fscanf(rf,"%d",n_graphs);

    *name =  malloc(sizeof(char *)*(*n_graphs));  // Allocate one pointer for each graph 

    if (*name == NULL) {
        printf("Unable to allocate space for names of %d graphs\n",max_graphs);
        return(0);
    }

    // Get number of runs per configuration
    fscanf(rf,"%d",n_runs);

    // Get number of thread configurations
    fscanf(rf,"%d",n_conf);
    *conf = (int *) malloc(sizeof(int)*(*n_conf));  // Allocate space for configurations


    if (*conf == NULL) {
        printf("Unable to allocate memory for %d different thread configurations \n",n_conf);
        return(0);
    }

    // Get the different configurations
    for(i=0;i< *n_conf;i++) {
        fscanf(rf,"%d",&(*conf)[i]);
    }

    // Get the different file names
    for(i=0;i< *n_graphs;i++) {
        (*name)[i] = (char *) malloc(sizeof(char)*100);
        if ((*name)[i] == NULL) {
            printf("Unable to allocate memory for graph name %d \n",i);
            return(0);
        }
        // Read name of graph
        fscanf(rf,"%s",(*name)[i]);
    }

    fclose(rf);
    return(true);
}


int allocate_memory(int size,int **n,int **m,double ***timer,double ***cost,double ****p_timer,int ****nrC,int n_conf) {


    *n = (int *) malloc(sizeof(int)*size);  // List holding the number of vertices in each graph

    if (*n == NULL) {
        printf("Unable to allocate memory for n[] in allocate_memory() \n");
        return(false);
    }

    *m = (int *) malloc(sizeof(int)*size);  // List holding the number of edges in each graph
    if (*m == NULL) {
        printf("Unable to allocate memory for m[] in allocate_memory() \n");
        return(false);
    }

    *timer = malloc(sizeof(double *)*(max_experiment));  // Allocate one pointer for each sequential experiment
    if (*timer == NULL) {
        printf("Unable to allocate memory for timer in allocate_memory() \n");
        return(false);
    }

    *cost = malloc(sizeof(double *)*(max_experiment));  // Allocate one pointer for each sequential experiment
    if (*timer == NULL) {
        printf("Unable to allocate memory for cost in allocate_memory() \n");
        return(false);
    }

    // For each experiment, allocate one number for each graph
    int i;
    for(i=0;i<max_experiment;i++) {
        (*cost)[i] = (double *) malloc(sizeof(double)*size);
        if ((*cost)[i] == NULL) {
            printf("Unable to allocate memory for cost %d in allocate_memory() \n",i);
            return(false);
        }
        (*timer)[i] = (double *) malloc(sizeof(double)*size);
        if ((*timer)[i] == NULL) {
            printf("Unable to allocate memory for timer %d in allocate_memory() \n",i);
            return(false);
        }
        int j;
        for(j=0;j<size;j++)          // Set each timer to -1
            (*timer)[i][j] = -1.0;
    }

    *p_timer = malloc(sizeof(double **)*(max_experiment));  // Allocate one pointer for each parallel experiment
    *nrC     = malloc(sizeof(int **)*(max_experiment));  // Allocate one pointer for each parallel experiment

    if (*p_timer == NULL) {
        printf("Unable to allocate memory for p_timer in allocate_memory() \n");
        return(false);
    }
    if (*nrC== NULL) {
        printf("Unable to allocate memory for nrC in allocate_memory() \n");
        return(false);
    }
    // For each experiment, allocate "size" pointers for each graph
    for(i=0;i<max_experiment;i++) {
        (*p_timer)[i] = (double **) malloc(sizeof(double *)*size);
        (*nrC)[i]  = (int **) malloc(sizeof(int *)*size);
        if ((*p_timer)[i] == NULL) {
            printf("Unable to allocate memory for p_timer %d in allocate_memory() \n",i);
            return(false);
        }
        if ((*nrC)[i] == NULL) {
            printf("Unable to allocate memory for nrC %d in allocate_memory() \n",i);
            return(false);
        }
        int j;
        for(j=0;j<size;j++) {
            (*p_timer)[i][j] = (double *)  malloc(sizeof(double)*n_conf);
            (*nrC)[i][j] = (int *)  malloc(sizeof(int)*n_conf);
            if ((*p_timer)[i][j] == NULL) {
                printf("Unable to allocate memory for p_timer %d in allocate_memory(),%d \n",i,j);
                return(false);
            }
            if ((*nrC)[i][j] == NULL) {
                printf("Unable to allocate memory for nrC %d,%d in allocate memory \n",i,j);
                return(false);
            }
            int k;
            for(k=0;k<n_conf;k++) {
                (*p_timer)[i][j][k] = -1.0;
                (*nrC)[i][j][k] = -1;
            }
        }
    }
    return(true);
}


int prepare_output(FILE **wf,int n_conf,int conf[]) {

    // Open file for writing of results
    *wf = fopen("results.m","w");
    if (*wf == NULL) {
        printf("Unable to open results.m for writing \n");
        return(false);
    }

    // print the thread configurations to file
    fprintf(*wf,"x = [");
    int i;
    for(i=0;i<n_conf;i++) 
        fprintf(*wf,"%d ",conf[i]);
    fprintf(*wf,"];\n");

    return(true);
}


// Compute a random ordering of the vertices

int random_order(int n,int *order) {

    int l;

    // Compute a random ordering of the vertices
    for(l=1;l<=n;l++) {
        order[l] = l;
    }
    for(l=1;l<n;l++) {
        long int x = random() % (long int) (n-l+1);
        x++;
        int tmp = order[n-l+1];
        order[n-l+1] = order[x];
        order[x] = tmp;
    }
    return(true);
}

// Allocating memory specifically for this graph

int allocate_graph_memory(double **weight,edge **e,int **ver,int **edges,int n,int **p,int **p1,int **p2,int **p3, int **p4) {		


    // 'weight' contains the weight of each edge, stored in the same way as the edge-lists
    *weight= (double *) malloc(2*sizeof(double)*(max_m));
    if (*weight== NULL) {
        printf("Unable to allocate space for weight-array in allocate_graph_memory() \n");
        return(false);
    }
    // The raw edge list is stored in e, i.e. in the same order as it was read in

    *e = (edge *) malloc(sizeof(edge)*(max_m));
    if (*e == NULL) {
        printf("Unable to allocate space for e-array in allocate_graph_memory() \n");
        return(false);
    }

    *ver = (int *) malloc(sizeof(int)*(2+n));
    if (*ver == NULL) {
        printf("Unable to allocate space for vertex-array in allocate_graph_memory() \n");
        return(false);
    }

    *edges = (int *) malloc(2*sizeof(int)*(max_m));
    if (*edges == NULL) {
        printf("Unable to allocate space for edges-array in allocate_graph_memory \n");
        return(false);
    }

    *p = (int *) malloc(sizeof(int)*(n+2));	// Storage for parent pointers 
    if (*p == NULL) {
        printf("Unable to allocate space for p[] in allocate_graph_memory() \n");
        return(false);
    }

    *p1 = (int *) malloc(sizeof(int)*(n+2));	// Extra storage of length n+2
    if (*p1 == NULL) {
        printf("Unable to allocate space for p1[] in allocate_graph_memory() \n");
        return(false);
    }

    *p2 = (int *) malloc(sizeof(int)*(n+2));	// Extra storage of length n+2
    if (*p2 == NULL) {
        printf("Unable to allocate space for p2[] in allocate_graph_memory() \n");
        return(false);
    }

    *p3 = (int *) malloc(sizeof(int)*(n+2));	// Extra storage of length n+2
    if (*p3 == NULL) {
        printf("Unable to allocate space for p3[] in allocate_graph_memory() \n");
        return(false);
    }

    *p4 = (int *) malloc(sizeof(int)*(n+2));	// Extra storage of length n+2
    if (*p4 == NULL) {
        printf("Unable to allocate space for p4[] in allocate_graph_memory() \n");
        return(false);
    }

    return(true);
}


void store_p_value(FILE *wf,char *s,double **values,int n,int n_conf) {
    int i,j;
    fprintf(wf,"%s = [",s); 
    for(i=0;i<n;i++) {
        for(j=0;j<n_conf;j++)  {
            fprintf(wf,"%lf ",values[i][j]);
        }
    }
    fprintf(wf,"];\n");
}

void store_value(FILE *wf,char *s,double *values,int n) {
    int i;
    fprintf(wf,"%s = [",s); 
    for(i=0;i<n;i++)  {
        fprintf(wf,"%lf ",values[i]);
    }
    fprintf(wf,"];\n");
}

int numberComponents(int n,int *p) {

    uint nrC = 0;
    uint i;
    for(i=1;i<=n;i++)
        if (p[i] == i)
            nrC++;
    return nrC;
}


int verify_bfs(int n,int *p1,int *p2) {

    int i;

    for(i=1;i<=n;i++) {
        if (p1[i] != p2[i]) {
            printf("Sequential: p[%d] = %d, Parallel: p[%d] = %d \n",i,p1[i],i,p2[i]);
            return(false);
        }
    }
    return(true);
}
