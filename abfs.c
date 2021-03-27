
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
void individual_bfs(int n, int *ver, int *edges, int *p, int *dist, int *local_S, int size);
void abfs(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T) 
{
    int i, j, k;
    int v, w;

    int local_counter;
    int *local_T;

    const int thread = omp_get_thread_num();
    const int num_threads = omp_get_num_threads();

    const int offset = 2;   // points to the starting index in S and T 
                            // (shared values are stored before this)

    int depth = 1;

    local_T = (int*) malloc(n * sizeof(int));

#pragma omp single
    {
        for (i = 1; i <= n; i++) {
            p[i] = -1;
            dist[i] = -1;
        }

        p[1] = 1;
        dist[1] = 0;
        S[offset] = 1;

        S[0] = 1;
        S[1] = 0;
    }


    while (S[0] != 0) {

        depth++;
        local_counter = 0;

        if (depth == 10) {  // once search reaches desired depth, break out of the
            break;          // while loop and contiune the search individually
        }

        // move all neighbors of vertices in S to local Ts in parallel

#pragma omp for
        for (i = 0; i < S[0]; i++) {
            v = S[offset + i];
            for (j = ver[v]; j < ver[v+1]; j++) {
                w = edges[j];
                if (p[w] == -1) {
                    p[w] = v;
                    dist[w] = dist[v] + 1;
                    local_T[local_counter++] = w;
                }
            }
        }

        T[offset + thread] = local_counter;

#pragma omp barrier

        // find the next size of S

#pragma omp single
        { 
            S[0] = 0; 
            for (i = 0; i < num_threads; i++) {
                S[0] += T[offset + i];
            }

            // prefix sum last index for each thread

            for (i = 1; i < num_threads; i++) {
                T[offset + i] += T[offset + i - 1];
            }

        }

        // each thread moves elements from local Ts to S


        k = T[offset + thread] - local_counter;     // find start index for this thread
        for (i = 0; i < local_counter; i++) {
            S[offset + k++] = local_T[i];
        }

        // reset local size in T

        T[offset + thread] = 0;
#pragma omp barrier

    } // end while


    if (S[0] != 0) {    // S[0] is only zero if search finished before it reached split depth

#pragma omp single
        { printf("broke loop at depth: %d, S size: %d\n", depth, S[0]); }

        // share S evenly between threads

        local_counter = 0;
        for (i = thread; i < S[0]; i += num_threads) {
            local_T[local_counter++] = S[i];
        }

        // each thread exhausts its local T

        individual_bfs(n, ver, edges, p, dist, local_T, local_counter);
    }
}


void individual_bfs(int n, int *ver, int *edges, int *p, int *dist, int *local_S, int size)
{
    int i, j, k;
    int v, w;
    int *local_T;
    int *temp;
    int writes;
    
    local_T = (int*) malloc(n * sizeof(int));
    
    while (size != 0) {
        for (i = 0; i < size; i++) {
            v = local_S[i];
            k = dist[v] + 1;
            for (j = ver[v]; j < ver[v+1]; j++) {
                w = edges[j];
                if (p[w] == -1) {
                    p[w] = v;
                    dist[w] = k;
                    local_T[writes++] = w;
                }
            }
        }
        temp = local_S;
        local_S = local_T;
        local_T = temp;

        size = writes;
        writes = 0;
    }

    free(local_S);
    free(local_T);
}
