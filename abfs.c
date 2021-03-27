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

#pragma omp single

    {
        for(i = 1; i <= n; i++) {   // Set that every node is unvisited
            p[i] = -1;              // Using -1 to mark that a vertex is unvisited
            dist[i] = -1;
        }

        p[1] = 1;        // Set the parent of starting vertex to itself
        dist[1] = 0;     // Set the distance from the starting vertex to itself
        S[offset] = 1;   // Add the starting vertex to S

        S[0] = 1;        // Number of vertices in S
        T[0] = 0;        // Number of vertices in T

        while (S[0] != 0) {                            // Loop until all vertices have been discovered
            depth++;
            if (depth > 10) break;

            for (i = 0; i < S[0]; i++) {               // Loop over vertices in S
                v = S[offset + i];                      // Grab next vertex v in S
                for (j = ver[v]; j < ver[v + 1]; j++) {   // Go through the neighbors of v
                    w = edges[j];               // Get next neighbor w of v
                    if (p[w] == -1) {           // Check if w is undiscovered
                        p[w] = v;               // Set v as the parent of w
                        dist[w] = dist[v] + 1;  // Set distance of w 
                        T[offset + T[0]++] = w; // Add w to T and increase number of vertices discovered 
                    }
                }  // End loop over neighbors of v
            }  // End loop of vertices in S
            temp = S;  // Swap S and T
            S = T;
            T = temp;
            T[0] = 0;
        }
    }



#pragma omp single nowait
    { printf("splitting search, S size: %d\n", S[0]); }

    int *local_S = (int*) malloc(n * sizeof(int));
    int *local_T = (int*) malloc(n * sizeof(int));

    const int tid = omp_get_thread_num();
    const int num_t = omp_get_num_threads();

    int num_r, num_w;

    // split vertices in S evenly between threads

    for (i = tid; i < S[0]; i += num_t) {
        local_S[num_r++] = S[offset + i];
    }
#pragma omp barrier

    // explore local vertices

    // first two indices in S are used for a custom barrier

#pragma omp single
    {
        T[0] = num_t;   // first lock
        T[1] = 0;       // first counter
        T[2] = true;

        S[0] = num_t;   // second lock 
        S[1] = 0;       // second counter
        S[2] = true;
    }

    while (num_r != 0) {

        // entry point
#pragma omp critical
        {
            S[1]++;
            printf("tid %d hit entry. (%d/%d)\n", tid, S[1], S[0]);
            if (S[1] == S[0]) {     // when the last thread arrives
                T[1] = 0;           // reset the other counter 
                T[2] = true;        // raise other barrier
                S[2] = false;       // release this barrier
            }
        }
        while(S[2]); // wait

        printf("tid %d passed entry\n", tid);

        for (i = 0; i < num_r; i++) {
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



        // exit point
#pragma omp critical
        {
            printf("tid %d hit exit. (%d/%d)\n", tid, T[1], T[0]);
            T[1]++;
            if (T[1] == T[0]) {
                S[1] = 0;       // reset first counter
                S[2] = true;    // raise first barrier
                T[2] = false;   // relase this barrier
            }
        }
        while(T[2]); // wait


    }
    // last barrier
#pragma omp critical
    {
        printf("tid %d hit final barrier\n", tid);
        S[0]--;
        T[0]--;
        if (S[1] == S[0]) {
            T[1] = 0;           // last thread to arrive resets the ohter counter
            T[2] = true;
            S[2] = false;
            printf("tid %d done\n", tid);
        }
    }

    free(local_S);
    free(local_T);

}

