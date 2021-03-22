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

void pbfs(int n,int *ver,int *edges,int *p,int *dist,int *S,int *T) 
{
    // Write code here
    int i, j;
    int v, w;
    int num_r, num_w;
    int *temp;

    // set all parents and distances to -1
    for (i = 1; i <= n; i++) {
        p[i] = -1;
        dist[i] = -1;
    }

    p[1] = 1;
    dist[1] = 0;
    S[0] = 1;

    num_r = 1;
    num_w = 0;

    while (num_r != 0) {
        for (i = 0; i < num_r; i++) {
            v = S[i];
            for (j = ver[v]; j < ver[v+1]; j++) {
                w = edges[j];
                if (p[w] == -1) {
                    p[w] = v;
                    dist[w] = dist[v] +1;
                    T[num_w++] = w;
                }
            }
        }
        temp = S;
        S = T;
        T = temp;
        num_r = num_w;
        num_w = 0;
    }
}

