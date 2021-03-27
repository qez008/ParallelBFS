// Parallel Breadth First Search
// -----------------------------
// Berforms a BFS starting from vertex 1
// The parent of each vertex in the BFS tree along with its distance from the starting
// vertex is computed.
//
// The algorithm should first perform some rounds of sequential BFS before starting a parallel
// execution. In the parallel part each thread should be allocated a part of the vertices from the
// last round of the sequential algorithm. Any discovered vertices in the parallel part should 
// remain with the thread that discovered them. This continues until the entire graph has been
// explored.
//
//
// Parameters:
// n     : number of vertices
// ver   : ver[i] points to the start of the neighbor list of vertex i in edges
// edges : lists of neighbors of each vertex, each edge is listed in both direction
// p     : array of length n used for parent pointers
// dist  : array of length n used for distance from starting vertex
// S     : array of length n used for maintaining queue of vertices to be processed, only used in the 
//         sequential part. 
// T     : array of length n where n >> number of threads. 
//
// Note that the vertices are numbered from 1 to n (inclusive). Thus there is
// no vertex 0.
void ibfs(int n, int *ver, int *edges, int *p, int *dist, int *S, const int offset);
void abfs(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T) 
{
    // write code here

    int i, j;         // Loop indices
    int v, w;         // Pointers to vertices
    int *temp;        // Temporary pointer

    const int offset = 2;

    // reqular bfs. code from sbfs exectued by a single thread.
    // when desired depth is reach split the search between threads.

    int depth = 0;
    int break_depth = 10;

    int *single_T;
    single_T = (int*) malloc(n * sizeof(int));

#pragma omp single

    {
        printf("single bfs performed by thread %d\n", omp_get_thread_num());


        for(i = 1; i <= n; i++) {   // Set that every node is unvisited
            p[i] = -1;              // Using -1 to mark that a vertex is unvisited
            dist[i] = -1;
        }

        p[1] = 1;        // Set the parent of starting vertex to itself
        dist[1] = 0;     // Set the distance from the starting vertex to itself
        S[offset] = 1;   // Add the starting vertex to S

        S[0] = 1;        // Number of vertices in S
        single_T[0] = 0;        // Number of vertices in T

        while (S[0] != 0) {                            // Loop until all vertices have been discovered
            depth++;
            if (depth >= break_depth) {
                //printf("reached break depth of %d\n", break_depth);
                break;
            }

            for (i = 0; i < S[0]; i++) {               // Loop over vertices in S
                v = S[offset + i];                      // Grab next vertex v in S
                for (j = ver[v]; j < ver[v + 1]; j++) {   // Go through the neighbors of v
                    w = edges[j];               // Get next neighbor w of v
                    if (p[w] == -1) {           // Check if w is undiscovered
                        p[w] = v;               // Set v as the parent of w
                        dist[w] = dist[v] + 1;  // Set distance of w 
                        // Add w to T and increase number of vertices discovered 
                        single_T[offset + single_T[0]++] = w; 
                    }
                }  // End loop over neighbors of v
            }  // End loop of vertices in S
            temp = S;  // Swap S and T
            S = single_T;
            single_T = temp;
            single_T[0] = 0;
        }
        
        printf("single bfs done\n");

        // move elements to shared T because pointer to S has switched
        for (i = 0; i < offset + S[0]; i++) {
            T[i] = S[i];
        }
    }

#pragma omp barrier

    //free(single_T);

    ibfs(n, ver, edges, p, dist, T, offset);
}


// synced individual bfs
void ibfs(int n, int *ver, int *edges, int *p, int *dist, int *S, const int offset)
{
    int i, j, k;
    int v, w;
    int *temp;

#pragma omp barrier

    int *local_S = (int*) malloc(n * sizeof(int));
    int *local_T = (int*) malloc(n * sizeof(int));

    const int tid = omp_get_thread_num();
    const int num_t = omp_get_num_threads();

    int num_r, num_w;

    // split vertices in S evenly between threads

    for (i = tid; i < S[0]; i += num_t) {
        local_S[num_r++] = S[offset + i];
    }

#pragma omp critical
    {
        printf("thread %d: ", tid);
        for (i = 0; i < num_r; i++)
            printf("%d ", local_S[i]);
        printf("\n");
    }

#pragma omp barrier

    k = 0;
    while (S[0] != 0) {
        //printf("thread %d depth: %d\n", tid, k++);

        for (i = 0; i < num_r; i++) {
            k++;
            v = local_S[i];
            for (j = ver[v]; j < ver[v+1]; j++) {
                w = edges[j];
                if (p[w] == -1) {
                    p[w] = v;
                    dist[w] = dist[v] + 1;
                    local_T[num_w++] = w;
                } 
            }
        }

        temp = local_S;
        local_S = local_T;
        local_T = temp;
        num_r = num_w;
        num_w = 0;

#pragma omp barrier
#pragma omp single
        { S[0] = 0; }
#pragma omp barrier
#pragma omp critical
        { S[0] += num_r; }
#pragma omp barrier

    }

    free(local_S);
    free(local_T);
    printf("thread %d done, explored %d vertices\n", tid, k);

#pragma omp barrier
}

