// Parallel Breadth First Search
// -----------------------------
// Berforms a BFS starting from vertex 1
// The parent of each vertex in the BFS tree along with its distance from the starting
// vertex is computed.
//
// The algorithm should gather all discovered vertices from round i, so that they can be 
// distributed among the threads before the search in round i+1.
//
// Parameters:
// n     : number of vertices
// ver   : ver[i] points to the start of the neighbor list of vertex i in edges
// edges : lists of neighbors of each vertex, each edge is listed in both direction
// p     : array of length n used for parent pointers
// dist  : array of length n used for distance from starting vertex
// S     : array of length n used for maintaining queue of vertices to be processed 
// T     : array of length n where n >> number of threads. 
//
// Note that the vertices are numbered from 1 to n (inclusive). Thus there is
// no vertex 0.
void pbfs(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T) 
{
    int i, j, k;
    int v, w;
    int *temp;

    int local_counter;
    int *local_T;

    local_T = (int*) malloc(n * sizeof(int));

    const int thread = omp_get_thread_num();
    const int num_threads = omp_get_num_threads();

    // points to the starting index in S and T (shared values are stored before this)
    const int offset = 2;

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

        local_counter = 0;

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


        k = T[offset + thread] - local_counter;
        for (i = 0; i < local_counter; i++) {
            S[offset + k++] = local_T[i];
        }

        // reset local counter in T

        T[offset + thread] = 0;

#pragma omp barrier

    } // end while


    free(local_T);

#pragma omp barrier

}

