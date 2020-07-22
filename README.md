# STATEMENT
The author (Hai-Jun Zhou) of this code is against any form of terrorism. This
code should only be used for pure academic purposes and for practical purposes
that improve the welfare of humanity. The author strongly condemn any intention
of adapting this code for destructive purposes.

# DESCRIPTION:
Belief-propagation guided decimation as a solver for the maximum feedback
vertex set problem on a directed simple graph. An arc [i,j] is a directed
edge from vertex i to vertex j. If two vertices i and j are connected by
bi-directional edges, namely there are two arcs [i, j] and [j, i] between
these two vertices, the interaction between vertex i and vertex j should be
treated specially:

> If h_i >0 then it is required that h_j=0; 
> If h_j>0 then it is required that h_i=0. 
> If there is a self-arc [i,i] from a vertex i to itself, this vertex
> should be in the empty state (h_i=0).

The present algorithm applies both to digraphs that are free of bi-directional
edges and to digraphs with bi-directional edges. Multiple arcs from a vertex i
to a vertex j are considered only once.

At each round of the decimation process, some of the remaining active vertices
are fixed to be un-occupied (i.e., moved to the DFVS set) sequentially. After
each vertex is fixed to be un-occupied, the remaining active subgraph is
simplified.

The messages received by a central vertex j are divided into two groups, the
messages q_{i->j}^{h_i} from in-comming arcs [i,j] and the messages
q_{k->j}^{h_k} from out-going arcs [j,k].

The version v04 is identical to v03, the only change is in the output format.

The version v03 differs from v02 in that if the number of active vertices is
less than MaxHeight, the actual maxHeight is set ot be the active vertices
number.

The version v02 differs from v01 in (1) the treatment of vertex structures
(two types of vertices are considered, one without bi-directional edges, one
with bi-directional edges), and (2) a final refinement process is added.

# REFERENCES
The program is applicable on a single graph instance. The theoretical details
of this BPD algorithm can be found at: 

> Hai-Jun Zhou, "A spin glass approach to the directed feedback vertex set problem",
> Journal of Statistical Mechanics: Theory and Experiment, 2016, 073303,
> doi:10.1088/1742-5468/2016/07/073303
> (see also: https://arxiv.org/abs/1604.00873 )

# COMPILE
To generate the executive file, you can simply compile as

`c++ -O3 -o dfvsbpd.exe DFVSbpdV04.cpp`

!!! please make sure some of the key parameters, such as input graph name,
number of edges, number of vertices, output file names, are appropriately
specified in the DFVSbpdV04.cpp file.

# USE
After you successfully compiled the code, you can then run the algorithm as

`dfvsbpd.exe`

# HISTORY:
- 08.05.2016: copied from DFVSbpdv03.cpp into DFVSbpdV04.cpp.
- 27.12.2015: copied from DFVSbpdv03MPI.cpp to DFVSbpdv03.cpp (basically not
change to the algorithm, only the random number generator is different).
The DFVScheck procedure is slightly optimized.
- 19.12.2013: copied from DFVSbpdv02MPI.cpp into DFVSbpdv03MPI.cpp 
- 17.12.2013: copied from DFVSbpdv01MPI.cpp into DFVSbpdv02MPI.cpp
- 10.12.2013: copied from DFVSbpdv00MPI.cpp into DFVSbpdv01MPI.cpp

# PROGRAMMER
- Hai-Jun Zhou
- Institute of Theoretical Physics, Chinese Academy of Sciences
- Zhong-Guan-Cun East Road 55, Beijing 100190
- email: zhouhj@itp.ac.cn
- website: http://home.itp.ac.cn/~zhouhj

