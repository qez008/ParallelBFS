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


void alt1(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T);
void alt2(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T);

void show(int* arr, int size);


void pbfs(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T) 
{
    alt2(n, ver, edges, p, dist, S, T);
    return;

    // Write code here
    int debug = 0;

    int iteration;

    int i, j, k;
    int v, w;
    int *temp;
    int t, ts;

    t = omp_get_thread_num();
    ts = omp_get_num_threads();

    int local_counter;
    int *local_T;

    int offset;

    offset = 2; // number of indices in T used for shared variables

    local_T = (int*) malloc(n * sizeof(int));

#pragma omp single

    {
        for (i = 1; i <= n; i++) {
            p[i] = -1;
            dist[i] = -1;
        }

        p[1] = 1;
        dist[1] = 0;


        T[0] = 1;
        T[1] = 0;

        S[0] = 1;
    }


#pragma omp barrier

    while (T[0] != 0) {

        iteration++;
        local_counter = 0;

        // divide vertices between threads

#pragma omp for
        for (i = 0; i < T[0]; i++ ) {
            v = S[i];

            // explore v's edges

            if (debug) printf("thread %d exploring %d \n", t, v);

            for (j = ver[v]; j < ver[v+1]; j++) {
                w = edges[j];

                // found vertex w

                if (p[w] == -1) {
                    p[w] = v;
                    dist[w] = dist[v] +1;

                    // add w to local T

                    local_T[local_counter++] = w;
                }
            }

            // store local counter in T + offset

            T[i + offset] = local_counter;

#pragma omp critical
            {
                if (debug) {
                    printf("thread %d found: ", t);
                    show(local_T, local_counter);
                }
            }

        } // end parallel for

        // wait for all threads to finish

#pragma omp barrier

        // all vertieces in S have been moved to local Ts

        // do the prefix sums to find out where each thread may place elements
        // and find the total size of the next S

#pragma omp single
        {
            // reset size of S and number of elements found
            T[0] = 0;
            T[1] = 0;

            // sum all local counters. this will be the next size of S
            for (i = 0; i < ts; i++) {
                T[0] += T[i + offset];
            }

            // find last write index for each thread
            for (i = 1; i < ts; i++) {
                T[i + offset] += T[i - 1 + offset];    
            }
        }

#pragma omp barrier


        // each thread moves its elements to S from its local T

        // printf("thread %d counter: %d\n", t, local_counter);

#pragma omp for
        for (i = 0; i < ts; i++) {

            k = T[i + offset] - local_counter;

            for (j = 0; j < local_counter; j++) {
                S[k++] = local_T[j];
            }

            // reset counter in T

            T[i + offset] = 0;

        } // end parallel for

#pragma omp barrier
#pragma omp single
        {
            if (debug) {
                printf("queue: ");
                show(S, T[0]);
                printf("---\n");
            }

        }
#pragma omp barrier
        if (debug && iteration == 10) {
#pragma omp single
            { printf("returned after %d iterations\n", iteration); }
            return;
        }
    }

    free(local_T);

#pragma omp barrier

}



void alt2(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T) 
{
    int i, j, k, t;
    int v, w;
    int *temp;

    int threads;
    threads = omp_get_num_threads();

    int offset;

    offset = 2;

    int local_counter;
    int *local_T;

    local_T = (int*) malloc(n * sizeof(int));

#pragma omp single
    {
        printf("--alt2--\n");

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
            T[offset + omp_get_thread_num()] = local_counter;
        }

#pragma omp barrier

#pragma omp single
        {
            // find the next size of S

            T[1] = 0;
            for (t = 0; t < threads; t++) {
                T[1] += T[offset + t];
            }

            S[0] = T[1];

            // prefix sum last index for each thread

            for (t = 1; t < threads; t++) {
                T[offset + t] += T[offset + t - 1];
            }

        }

        // finally each thread moves elements from local T to S

#pragma omp for
        for (t = 0; t < threads; t++) {

            // find start index for thread t

            k = T[offset + t] - local_counter;
            for (i = 0; i < local_counter; i++) {
                S[offset + k++] = local_T[i];
            }

            // reset local counter in T

            T[offset + t] = 0;
        }
#pragma omp barrier

    } // end while


    free(local_T);
}



// regular bfs exectued by first arriving thread 
void alt1(int n, int *ver, int *edges, int *p, int *dist, int *S, int *T) 
{
    int i, j, k;
    int v, w;
    int *temp;

    int offset;

    offset = 2;

    int local_counter;
    int *local_T;

    local_T = (int*) malloc(n * sizeof(int));

#pragma omp single
    {
        printf("--alt1--\n");
        printf("thread %d resetting\n", omp_get_thread_num());

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

#pragma omp single
    {
        printf("thread %d running bfs\n", omp_get_thread_num());

        while (S[0] != 0) {

            local_counter = 0;

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

            T[1] = local_counter;
            for (i = 0; i < local_counter; i++) {
                T[offset + i] = local_T[i];
            }

            T[0] = T[1];
            T[1] = 0;
            temp = S;
            S = T;
            T = temp;
        }
    }
}


void show(int *arr, int size)
{
    int i;
    for (i = 0; i < size; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}
