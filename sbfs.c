
// Sequential Breadth First Search
// -------------------------------
// Berforms a BFS starting from vertex 1
// The parent of each vertex along with its distance from the starting
// vertex is computed.
//
// Parameters:
// n     : number of vertices
// ver   : ver[i] points to the start of the neighbor list of vertex i in edges
// edges : lists of neighbors of each vertex, each edge is listed in both direction
// p     : array of length n used for parent pointers
// dist  : array of length n used for distance from starting vertex
// S     : array of length n used for maintaining queue of vertices to be processed 
// T     : array of length n used for maintaining queue of newly discovered vertices  
//
// Note that the vertices are numbered from 1 to n (inclusive). Thus there is
// no vertex 0.

void sbfs(int n,int *ver,int *edges,int *p,int *dist,int *S,int *T) 
{

    int i,j;          // Loop indices
    int v,w;          // Pointers to vertices
    int num_r,num_w;  // Number of vertices in S and T, respectively
    int *temp;        // Temporary pointer

    for(i=1;i<=n;i++) {   // Set that every node is unvisited
        p[i] = -1;        // Using -1 to mark that a vertex is unvisited
        dist[i] = -1;
    }

    p[1] = 1;        // Set the parent of starting vertex to itself
    dist[1] = 0;     // Set the distance from the starting vertex to itself
    S[0] = 1;        // Add the starting vertex to S

    num_r = 1;       // Number of vertices in S
    num_w = 0;       // Number of vertices in T

    while (num_r != 0) {                    // Loop until all vertices have been discovered
        for(i=0;i<num_r;i++) {              // Loop over vertices in S
            v = S[i];                       // Grab next vertex v in S
            for(j=ver[v];j<ver[v+1];j++) {  // Go through the neighbors of v
                w = edges[j];               // Get next neighbor w of v
                if (p[w] == -1) {           // Check if w is undiscovered
                    p[w] = v;               // Set v as the parent of w
                    dist[w] = dist[v]+1;    // Set distance of w 
                    T[num_w++] = w;         // Add w to T and increase number of vertices discovered 
                }
            }  // End loop over neighbors of v
        }  // End loop of vertices in S
        temp = S;  // Swap S and T
        S = T;
        T = temp;
        num_r = num_w; // Set number of elements in S
        num_w = 0;     // Set T as empty
    } //  End loop over entire graph
}
