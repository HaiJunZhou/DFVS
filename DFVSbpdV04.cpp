/*
*******************************************************************************
DFVSbpdV04.cpp

DESCRIPTION:
Belief-propagation guided decimation as a solver for the maximum feedback
vertex set problem on a directed simple graph. An arc [i,j] is a directed
edge from vertex i to vertex j. If two vertices i and j are connected by
bi-directional edges, namely there are two arcs [i, j] and [j, i] between
these two vertices, the interaction between vertex i and vertex j should be
treated specially:
if h_i >0 then it is required that h_j=0; if h_j>0 then it is required that
h_i=0. If there is a self-arc [i,i] from a vertex i to itself, this vertex
should be in the empty state (h_i=0).

The program is applicable on a single graph instance. The theoretical details
of this BPD algorithm can be found at:
Hai-Jun Zhou, 
"A spin glass approach to the directed feedback vertex set problem",
Journal of Statistical Mechanics: Theory and Experiment, 2016, 073303,
doi:10.1088/1742-5468/2016/07/073303
(see also: https://arxiv.org/abs/1604.00873 )

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

LOG:
08.05.2016: copied from DFVSbpdv03.cpp into DFVSbpdV04.cpp.
27.12.2015: copied from DFVSbpdv03MPI.cpp to DFVSbpdv03.cpp (basically not
change to the algorithm, only the random number generator is different).
The DFVScheck procedure is slightly optimized.
copied from DFVSbpdv02MPI.cpp into DFVSbpdv03MPI.cpp on 19.12.2013
copied from DFVSbpdv01MPI.cpp into DFVSbpdv02MPI.cpp on 17.12.2013
copied from DFVSbpdv00MPI.cpp into DFVSbpdv01MPI.cpp on 10.12.2013.

PROGRAMMER:
Hai-Jun Zhou
Institute of Theoretical Physics, Chinese Academy of Sciences
Zhong-Guan-Cun East Road 55, Beijing 100190
zhouhj@itp.ac.cn
http://power.itp.ac.cn/~zhouhj/
*******************************************************************************
*/

#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <valarray>

//#include "/Users/zhouhj/Programs/zhjrandom.h"
#include "zhjrandom.h"

using namespace std;

//                                   random real uniformly distributed in [0,1)
double u01prn(ZHJRANDOMv3 *rd)
{
  return rd->rdflt();
}


//    IntPair describes an arc [i,j] from vertex i (first) to vertex j (second)
struct IntPair
{
  int first;           
  int second;
  IntPair(void);
  IntPair(int, int);
};

IntPair::IntPair(void)
{
  first = 0;
  second = 0;
  return ;
}

IntPair::IntPair(int a, int b)
{
  first = a;
  second = b;
  return ;
}

bool operator<(IntPair a, IntPair b)
{
  if(a.first < b.first) return true;
  else if(a.first > b.first) return false;
  else {                                                    //a.first = b.first
    if(a.second<b.second) return true;
    else return false; }
}

bool operator==(IntPair a, IntPair b)
{
  return (a.first == b.first) && (a.second == b.second);
}

//vertex structure
struct vertexstruct
{  
  int index;                                //index of vertex, positive integer
  int in_degree;                                     //number of inputing edges
  int active_in_degree;                      //number of active parent vertices
  int out_degree;                                      //number of output edges
  int active_out_degree;                      //number of active child vertices
  int inout_degree;                            //number of bi-directional edges
  int active_inout_degree;              //number of active bi-directional edges
  double q_0;                                               //empty probability
  bool active; /* true (not being removed from graph); false (has been deleted
		  from graph). an inactive vertex might be occupied. */
  bool occupied;                  //vertex being occupied, not belonging to FVS
  struct messagestruct *pim_ptr; //start position of input message from parents
  struct messagestruct *cim_ptr;//start position of input message from children
  struct messagestruct *pcim_ptr; /* start position of input message from
				     bi-directional edges */
  struct simplifiedmessagestruct *opm_ptr; /* start position of output message
					      to parents */
  struct simplifiedmessagestruct *ocm_ptr; /* start position of output message
					      to children */
  struct simplifiedmessagestruct *opcm_ptr; /* start position of output message
					       to bi-directional edges */
  vertexstruct(void);
};

vertexstruct::vertexstruct(void)
{
  index             = 0;
  in_degree         = 0;
  active_in_degree  = 0;
  out_degree        = 0;
  active_out_degree = 0;
  inout_degree      = 0;
  active_inout_degree=0;
  q_0               = 0;
  active            = true;                 //initially all vertices are active
  occupied          = true;               //initially all vertices are occupied
  pim_ptr=0;
  cim_ptr=0;
  pcim_ptr=0;
  opm_ptr=0;
  ocm_ptr=0;
  opcm_ptr=0;
  return;
}

/* structure of input message to a vertex. this meassage includes a pointer
   to the neighboring (parent or child) vertex, a pointer to the cavity
   probability profile of the neighboring vertex. */
struct messagestruct
{
  struct vertexstruct *v_ptr; /* pointer to the other vertex which sends the
				 message */
  double *q_cavity; /* pointer to q_{j->k}^{0}, the 1st element of
		       q_{j->k}^{h_j}. The cavity probability distribution
		       q_{j->k} is represented by a sequence of double-
		       precision numbers. The length of this sequence is
		       MaxHeight+1. */
  messagestruct(void);
};

messagestruct::messagestruct(void)
{
  v_ptr    =0;
  q_cavity =0;
  return ;
}

struct simplifiedmessagestruct
{
  double *q_cavity; /* pointer to q_{j->k}^{0}, the 1st element of
		       q_{j->k}^{h_j}. The cavity probability distribution
		       q_{j->k} is represented by a sequence of double-
		       precision numbers. The length of this sequence is
		       MaxHeight+1. */
  simplifiedmessagestruct(void);
};

simplifiedmessagestruct::simplifiedmessagestruct(void) {
  q_cavity =0;
  return ;
}

//the data structure for the directed feedback vertex set problem
class DFVS
{
public:
  //  DFVS(trng::lcg64_shift & );
  DFVS(ZHJRANDOMv3* );
  void SetX(double);                             //set re-weighting parameter X
  void SetH(int);                                 //set value of maximal height
  void SetDamping(double);                                 //set damping factor
  bool Digraph(string&, string&, int, bool);           //read graph
  void Initialization(void);                              //initialize messages
  bool BeliefPropagation(double, int);                    //population updating
  void Thermodynamics(string&);                //thermodynamic quantities
  void VertexRanking(void);                     //select vertices to be removed
  int  Fix0(void );         //fix some variables to be un-occupied and simplify
  void Simplify(int);                     //simplify the graph by removing leaf
  bool CheckDFVS(string&); /* check whether the occupation pattern
				    corresponds to a FVS */
  bool RefineDFVS(string&); /* check whether the occupation pattern
				     corresponds to a FVS, and if yes, then
				     refine the DFVS pattern. */
  
private:
  int VertexNumber;                                  //total number of vertices
  int ActiveVertexNumber;                     //total number of active vertices
  int ActiveVertexNumberOld;  /* old value of total number of active vertices,
				 needed in updating the set Permutation */
  int EdgeNumber;                                       //total number of edges
  int OneArrowEdgeNumber; /* total number of one-directional edges between
			     vertices */
  int TwoArrowEdgeNumber;     /* total number of bi-directional edges
				 between vertices */
  int SelfArcNumber; //number of self arcs [i,i]
  int MultArcNumber; //number of multiple arcs [i, j] between the same two vtx
  
  int Max_in_degree;                    //maximal vertex in-degree in the graph
  int Max_out_degree;                  //maximal vertex out-degree in the graph
  int Max_inout_degree;              //maximal vertex inout-degree in the graph

  int MaxHeight;                                       //maximal allowed height
  int HeightOld;                                     //old maximal Height value
  int HeightNew;                                     //new maximal height value
  
  double X;                                            //re-weighting parameter
  double Weight0;                                                    //=exp(-X)
  
  double WeightOld;       //in message updating, weight of old value of message
  double WeightNew;       //in message updating, weight of new value of message
  
  float FixFraction;
  /* fraction of active vertices to be fixed to un-occupied in each round of
     the belief-propagation decimation process */
  int MaxFixNumber;             //maximal number of fixed vertices in one round
  /*
    trng::lcg64_shift PRNG;                           //random number generator
    trng::uniform01_dist<> u01prn;                             //[0,1) uniform
  */
  ZHJRANDOMv3* PRNG;
  
  set<int> SelfArcVertices; //vertices with self-arc attached
  
  valarray<vertexstruct> Vertex;                          //the set of vertices
  valarray<int> Permutation; /* permutation of vertex order. This array is used
				in random sequential updating. */
  
  valarray<double> Message;
  /* messages q_{i->j} and q_{j->i} for each arc [i,j]. the messages are listed
     in the following order:
     q_{i1->j1}, q_{i1<-j1} (for arc [i1,j1]),
     q_{i2->j2}, q_{i2<-j2} (for arc [i2,j2]),
     ....
     q_{iM->jM}, q_{iM<-jM} (for arc [iM,jM]). */
  
  valarray<messagestruct> MessageFromParent;
  /* Stores local graph information. parental connection to vertex. message
     from a parent vertex i of a central vertex j, namely the cavity
     probability profile p_{i->j} on an arc [i,j]. the probability profile is
     composed of (MaxHeight+1) elements, namely p_{i->j}^{h_i}, with h=0, 1,
     ..., MaxHeight. */
  valarray<simplifiedmessagestruct> MessageToParent;
  //stores local graph information. parental connection to vertex

  valarray<messagestruct> MessageFromChild;
  /* Stores local graph information. child connection to vertex. message from
     a child vertex k of a central vertex j, namely the cavity probability
     profile p_{k->j} on an arc [j,k]. The probability profile is composed of
     (MaxHeight+1) elements, namely p_{k->j}^{h_k}, with h=0, 1, ...,
     MaxHeight.*/
  valarray<simplifiedmessagestruct> MessageToChild;
  //stores local graph inforamtion. child connection to vertex

  valarray<messagestruct> MessageFromBDE; 
  /* Message from a neighboring vertex k of a central vertex j. There is a
     bidirectional edge (BDE) between vertex k and vertex j. */
  valarray<simplifiedmessagestruct> MessageToBDE;
  /* message to a neighboring vertex k of a central vertex j. There is a
     bidirectional edge (BDE) between vertex k and vertex j. */
  
  valarray<double> NewMessageToParent;     //newly generated messages to parent
  valarray<double> NewMessageToChild;       //newly generated messages to child
  valarray<double> NewMessageToBDE;
  /* newly generated messages to a neighboring vertex k of a central vertex j.
     There is a bi-directional edge (BDE) between a vertex k and a vertex j. */
  valarray<double> NormalizationOPM;
  valarray<double> NormalizationOCM;
  valarray<double> NormalizationBDE;
  
  valarray<double> VectorIn;
  valarray<double> VectorOut;
  valarray<double> VectorBDE;
  valarray<bool>   IsActiveParent;
  valarray<bool>   IsActiveChild;
  valarray<bool>   IsActiveBDE;
  //auxiliary arrays needed for message updating
  
  valarray<int> CandidateVertices;     //list of candidate vertices to be fixed
  valarray<int> CandidateSize;          /* number of candidate vertices in each
					   empty_prob range */
  
  int MinRange; /* minimal range of empty_prob of fixed vertices. MinRange is
		   the maximal value of r0 such that
		   sum_{r>=r0} CandidateSize[r] >= MaxFixNumber.  */
  
  void UpdateMessage(struct vertexstruct *, double&);
  void UpdateMessage(struct vertexstruct *);
};


int main(int argc, char ** argv)
{
  int rdseed=63456792;
  ZHJRANDOMv3 PRNG(rdseed);
  for(int i=0; i<10000119; ++i)
    PRNG.rdflt();
  DFVS system( &PRNG );
  
  //maximal allowed height
  int MaxHeight=200;

  //reweighting factor (inverse temperature)
  double X=50;

  //damping factor (weight of new message)
  double WeightNew=0.99e0;
  
  system.SetH(MaxHeight);
  system.SetX( X );
  system.SetDamping(WeightNew);
  
  //          specify the total number of vertices and the total number of arcs
  int VertexNumber = 10000;
  int ArcNumber    = 100000;
  //                                                   specify the digraph name
  string graphfile="RRRa10n10k.g00";
  string thermfile="RRRa10n10k.g00bp";
  string DFVSfile="RRRa10n10kg00.DFVS00";

  bool succeed=system.Digraph(graphfile, thermfile, ArcNumber, true);
  if(succeed==false)
    return -1;
  
  system.Initialization();
  
  float epsilon = 1.0e-15;
  int eiterations = 100;
  system.BeliefPropagation(epsilon, eiterations);
  system.Thermodynamics(thermfile);
  
  epsilon = 1.0e-7 ;
  int iterations=10;  
  int ActiveVertexNumber=0;
  do
    {
      system.BeliefPropagation(epsilon, iterations);
      system.VertexRanking();
      ActiveVertexNumber=system.Fix0();
      cout<<ActiveVertexNumber<<endl;
    }
  while(ActiveVertexNumber > 0);
  cout<<"DFVS refine.\n";
  if(system.RefineDFVS(DFVSfile)==true)
    return 1;
  else
    {
      cerr<<"Not a feedback vertex set.\n";
      return -1;
    }
}

DFVS::DFVS(ZHJRANDOMv3 *rd)
{
  PRNG = rd;                                          //random number generator
  FixFraction = 0.001;  //fraction of active vertices to be fixed in each round
  MaxHeight = 0;               //MaxHeight needs to be set through SetH()
  HeightOld=MaxHeight;
  HeightNew=MaxHeight;
  return ;
}

//set the value of maximal height
void DFVS::SetH(int dval)
{
  MaxHeight = dval;
  HeightOld = MaxHeight;
  HeightNew = MaxHeight;
  return;
}


//set the value of X
void DFVS::SetX(double xval)
{
  X = xval;                                            //re-weighting parameter
  Weight0=exp(-X);
  return ;
}

//set the value of DampingFactor
void DFVS::SetDamping(double nw)
{
  WeightNew=nw;                            //weight of newly calculated message
  if(WeightNew >=1.0e0) WeightNew=1.0e0;
  else if(WeightNew <=0.0001e0) WeightNew=0.0001e0;
  WeightOld=1.0e0-WeightNew;                            //weight of old message
  return;
}

//Read digraph from file
bool DFVS::Digraph
(
  string& gname,                                                 //digraph
  string& ofname,                  //output some statistics of the digraph
 int enumber,                                  //number of arcs to be read into
 bool signal                       //whether output some statistics  
 ) {
  if(MaxHeight<=0) {   //MaxHeight not yet been set
    cerr<<"MaxHeight not yet set.\n";
    return false; }

  ifstream graphf( gname.c_str() );  if( !graphf.good() ) {
    cerr<<"Graph probably non-existant.\n";
    return false;  }

  graphf >>VertexNumber
	 >>EdgeNumber;                //total number of arcs in the input graph
  if(EdgeNumber<enumber) {
    cerr<<"No so many edges in the graph.\n";
    graphf.close();
    return false;  }
  EdgeNumber=enumber;                   //only the first enumber edges are read

  try { Vertex.resize(VertexNumber+1); } //vertex index starts from 1
  catch(bad_alloc) {
    cerr<<"Vertex construction failed.\n";
    graphf.close();
    return false;  }

  //graph reading
  OneArrowEdgeNumber=0;
  TwoArrowEdgeNumber=0;
  SelfArcNumber=0;
  MultArcNumber=0;

  set<IntPair> OneArrowEdgeSet;
  set<IntPair> TwoArrowEdgeSet;
  SelfArcVertices.clear();

  bool succeed=true;
  for(int eindex=0; eindex<EdgeNumber && succeed; ++eindex)  {
    int v1, v2;
    graphf >>v1>>v2;
    if(v1 <=0 || v1 >VertexNumber || v2<=0 || v2 >VertexNumber) {
      cerr<<"Graph incorrect at line "<<eindex+1<<endl;
      succeed=false; }
    else {
      if( OneArrowEdgeSet.find(IntPair(v1,v2) ) != OneArrowEdgeSet.end() )
	++MultArcNumber;
      else { //not a multiple arc
	if(v1 == v2) { //self arc
	  ++SelfArcNumber;
	  Vertex[v1].active = false;   //inactive due to self-arc
	  Vertex[v1].occupied = false; //belongs to DFVS
	  SelfArcVertices.insert(v1); }
	else if( OneArrowEdgeSet.find(IntPair(v2,v1) ) 
		 != OneArrowEdgeSet.end() ) {  //bi-directional edge
	  ++TwoArrowEdgeNumber;
	  --OneArrowEdgeNumber;
	  TwoArrowEdgeSet.insert( IntPair(v1, v2) );
	  --(Vertex[v2].out_degree);
	  --(Vertex[v1].in_degree);
	  ++(Vertex[v1].inout_degree);
	  ++(Vertex[v2].inout_degree); }
	else { //treated as one-directional edge for the moment
	  ++OneArrowEdgeNumber;
	  ++(Vertex[v1].out_degree);
	  ++(Vertex[v2].in_degree); }
	OneArrowEdgeSet.insert( IntPair(v1, v2) ); } } }
  graphf.close();
  if(succeed==false) return false;

  /* OneArrowEdgeSet contains at the moment all the (non-redundant) arcs
     of the digraph. It therefore needs to be polished */
  for(set<IntPair>::const_iterator sci=TwoArrowEdgeSet.begin();
      sci != TwoArrowEdgeSet.end(); ++sci) {
    int v1= sci->first;
    int v2= sci->second;
    OneArrowEdgeSet.erase( IntPair(v1, v2) );
    OneArrowEdgeSet.erase( IntPair(v2, v1) ); }
  for(set<int>::const_iterator sci=SelfArcVertices.begin();
      sci != SelfArcVertices.end(); ++sci) {
    int v1 = *sci;
    OneArrowEdgeSet.erase( IntPair(v1, v1) ); }

  if(signal) { //output basic statistics of digraph
    ofstream output(ofname.c_str() );
    output<<gname<<" has "<<VertexNumber<<" vertices, and the first "<<
      EdgeNumber<<" arcs are read.\n"
	  <<"# multi-arcs = "<<MultArcNumber<<endl
	  <<"# self-arcs = "<<SelfArcNumber<<endl
	  <<"# oppositely directed arc pairs = "<<TwoArrowEdgeNumber<<endl
	  <<"# normal (one-directional) arcs = "<<OneArrowEdgeNumber<<endl;
    output.close(); }

  try { 
    Message.resize( 2*(TwoArrowEdgeNumber+OneArrowEdgeNumber)*(MaxHeight+1) );
  } catch(bad_alloc) {
    cerr<<"Message construction failed.\n";
    return false; }

  try { MessageFromParent.resize( OneArrowEdgeNumber ); } catch(bad_alloc) {
    cerr<<"MessageFromParent construction failed.\n";
    return false; }
  try { MessageToParent.resize( OneArrowEdgeNumber ); } catch(bad_alloc) {
    cerr<<"MessageToParent construction failed.\n";
    return false;  }

  try { MessageFromChild.resize( OneArrowEdgeNumber ); } catch(bad_alloc) {
    cerr<<"MessageFromChild construction failed.\n";
    return false;  }
  try { MessageToChild.resize( OneArrowEdgeNumber ); } catch(bad_alloc) {
    cerr<<"MessageToChild construction failed.\n";
    return false; }

  try { MessageFromBDE.resize( 2*TwoArrowEdgeNumber ); } catch(bad_alloc) {
    cerr<<"MessageFromBDE construction failed.\n";
    return false; }
  try { MessageToBDE.resize( 2*TwoArrowEdgeNumber ); } catch(bad_alloc) {
    cerr<<"MessageToBDE construction failed.\n";
    return false; }

  int index_in_edge=0;                                      //index of in-edges
  int index_out_edge=0;                                    //index of out-edges
  int index_inout_edge=0;                                //index of inout-edges

  Max_in_degree=0;                                          //maximal in-degree
  Max_out_degree=0;                                        //maximal out-degree
  Max_inout_degree=0;                                    //maximal inout-degree

  struct vertexstruct *v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v) {
    v_ptr->index=v;
    v_ptr->active_in_degree = v_ptr->in_degree;
    v_ptr->active_out_degree = v_ptr->out_degree;
    v_ptr->active_inout_degree = v_ptr->inout_degree;
    if(v_ptr->in_degree>Max_in_degree) Max_in_degree=v_ptr->in_degree;
    if(v_ptr->out_degree>Max_out_degree) Max_out_degree=v_ptr->out_degree;
    if(v_ptr->inout_degree>Max_inout_degree)
      Max_inout_degree=v_ptr->inout_degree;
    v_ptr->pim_ptr = &MessageFromParent[index_in_edge];
    v_ptr->opm_ptr = &MessageToParent[index_in_edge];
    index_in_edge += v_ptr->in_degree;
    v_ptr->in_degree=0;
    v_ptr->cim_ptr = &MessageFromChild[index_out_edge];
    v_ptr->ocm_ptr = &MessageToChild[index_out_edge];
    index_out_edge += v_ptr->out_degree;
    v_ptr->out_degree=0;
    v_ptr->pcim_ptr = &MessageFromBDE[index_inout_edge];
    v_ptr->opcm_ptr = &MessageToBDE[index_inout_edge];
    index_inout_edge += v_ptr->inout_degree;
    v_ptr->inout_degree=0;
    ++v_ptr; }

  int index_message=0;    //notice that q_{i->j} and q_{j->i} are two messages
  for(set<IntPair>::const_iterator sci = OneArrowEdgeSet.begin();
      sci != OneArrowEdgeSet.end(); ++sci) {
    struct vertexstruct *v1_ptr=&Vertex[ sci->first ];
    struct vertexstruct *v2_ptr=&Vertex[ sci->second ];
    //construction for vertex v1 (=sci->first)
    struct messagestruct *im_ptr=v1_ptr->cim_ptr;
    im_ptr[v1_ptr->out_degree].v_ptr=v2_ptr;
    im_ptr[v1_ptr->out_degree].q_cavity
      =&Message[(index_message+1)*(MaxHeight+1)];
    struct simplifiedmessagestruct *om_ptr=v1_ptr->ocm_ptr;
    om_ptr[v1_ptr->out_degree].q_cavity=&Message[index_message*(MaxHeight+1)];
    v1_ptr->out_degree += 1;
    //construction for vertex v2 (=sci->second)
    im_ptr=v2_ptr->pim_ptr;
    im_ptr[v2_ptr->in_degree].v_ptr=v1_ptr;
    im_ptr[v2_ptr->in_degree].q_cavity=&Message[index_message*(MaxHeight+1)];
    om_ptr=v2_ptr->opm_ptr;
    om_ptr[v2_ptr->in_degree].q_cavity
      =&Message[(index_message+1)*(MaxHeight+1)];
    v2_ptr->in_degree += 1;
    index_message += 2; }
  OneArrowEdgeSet.clear();

  for(set<IntPair>::const_iterator sci = TwoArrowEdgeSet.begin();
      sci != TwoArrowEdgeSet.end(); ++sci) {
    struct vertexstruct *v1_ptr=&Vertex[ sci->first ];
    struct vertexstruct *v2_ptr=&Vertex[ sci->second ];
    //construction for vertex v1 (=sci->first)
    struct messagestruct *im_ptr=v1_ptr->pcim_ptr;
    im_ptr[v1_ptr->inout_degree].v_ptr=v2_ptr;
    im_ptr[v1_ptr->inout_degree].q_cavity
      =&Message[(index_message+1)*(MaxHeight+1)];
    struct simplifiedmessagestruct *om_ptr=v1_ptr->opcm_ptr;
    om_ptr[v1_ptr->inout_degree].q_cavity=
      &Message[index_message*(MaxHeight+1)];
    v1_ptr->inout_degree += 1;
    //construction for vertex v2 (=sci->second)
    im_ptr=v2_ptr->pcim_ptr;
    im_ptr[v2_ptr->inout_degree].v_ptr=v1_ptr;
    im_ptr[v2_ptr->inout_degree].q_cavity=
      &Message[index_message*(MaxHeight+1)];
    om_ptr=v2_ptr->opcm_ptr;
    om_ptr[v2_ptr->inout_degree].q_cavity
      =&Message[(index_message+1)*(MaxHeight+1)];
    v2_ptr->inout_degree += 1;
    index_message += 2; }
  TwoArrowEdgeSet.clear();

  try { Permutation.resize(VertexNumber); } catch(bad_alloc) {
    cerr<<"Permutation construction failed.\n";
    return false; }   //vertex order permutation needed in random sequential BP

  //update the permutation array and the active degrees of each vertex
  ActiveVertexNumber=0;
  for(int v=1; v<=VertexNumber; ++v) {
    if(Vertex[v].active) Permutation[ActiveVertexNumber++ ] = v;
    else {                  //vertex inactive (and unoccupied) due to self-arc
      struct messagestruct *im_ptr=Vertex[v].pim_ptr;
      for(int i=0; i<Vertex[v].in_degree; ++i) {
	--(im_ptr->v_ptr->active_out_degree);
	++im_ptr; }
      im_ptr=Vertex[v].cim_ptr;
      for(int k=0; k<Vertex[v].out_degree; ++k) {
	--(im_ptr->v_ptr->active_in_degree);
	++im_ptr; }
      im_ptr=Vertex[v].pcim_ptr;
      for(int l=0; l<Vertex[v].inout_degree; ++l) {
	--(im_ptr->v_ptr->active_inout_degree);
	++im_ptr; } } }
  ActiveVertexNumberOld=ActiveVertexNumber;

  //recursively removing vertices whose in-degree or out-degree is zero
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v) {
    if(v_ptr->active) {
      if( v_ptr->active_inout_degree ==0 && 
	  (v_ptr->active_in_degree ==0 || v_ptr->active_out_degree ==0) ) {
	v_ptr->active = false;
	v_ptr->occupied = true;                                //being occupied
	--ActiveVertexNumber;
	Simplify(v_ptr->index); } }
    ++v_ptr; }

  if(ActiveVertexNumber==0) {
    cerr<<"The graph has no directed cycles. Done! \n";
    return false; }

  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0) MaxFixNumber = 1;

  //the empty probability q_0 is divided into bins of width 0.01
  try { CandidateVertices.resize( MaxFixNumber * 101 ); } catch( bad_alloc ) {
    cerr<<"CandidateVertices construction failed.\n";
    return false; }
  CandidateSize.resize(101); 

  MinRange=0;

  //accumulative values from parent vertices to a central vertex
  try { VectorIn.resize(Max_in_degree); } catch(bad_alloc) {
    cerr<<"VectorIn construction failed.\n";
    return false; }
  //accumulative values from child vertices to a central vertex
  try { VectorOut.resize(Max_out_degree); } catch(bad_alloc) {
    cerr<<"VectorOut construction failed.\n";
    return false; }
  //accumulative values from BDE connected vertices to a central vertex
  try { VectorBDE.resize( Max_inout_degree); } catch(bad_alloc) {
    cerr<<"VectorBDE construction failed.\n";
    return false; }

  try { NewMessageToParent.resize((MaxHeight + 1)*Max_in_degree); }
  catch(bad_alloc) {
    cerr<<"NewMessageToParent construction failed.\n";
    return false; }
  try { NormalizationOPM.resize(Max_in_degree); } catch(bad_alloc) {
    cerr<<"NormalizationOPM construction failed.\n";
    return false; }

  try { NewMessageToChild.resize((MaxHeight + 1)*Max_out_degree); }
  catch(bad_alloc) {
    cerr<<"NewMessageToChild construction failed.\n";
    return false; }
  try { NormalizationOCM.resize(Max_out_degree); } catch(bad_alloc) {
    cerr<<"NormalizationOCM construction failed.\n";
    return false; }

  try { NewMessageToBDE.resize((MaxHeight + 1)*Max_inout_degree); }
  catch(bad_alloc) {
    cerr<<"NewMessageToBDE construction failed.\n";
    return false; }
  try { NormalizationBDE.resize(Max_inout_degree); } catch(bad_alloc) {
    cerr<<"NormalizationBDE construction failed.\n";
    return false; }

  try { IsActiveParent.resize(Max_in_degree); } catch(bad_alloc) {
    cerr<<"IsActiveParent construction failed.\n";
    return false; }
  try { IsActiveChild.resize(Max_out_degree); } catch(bad_alloc) {
    cerr<<"IsActiveChild construction failed.\n";
    return false; }
  try { IsActiveBDE.resize(Max_inout_degree); } catch(bad_alloc) {
    cerr<<"IsActiveBDE construction failed.\n";
    return false; }

  return true;
}

//simplify after fixing vertex
void DFVS::Simplify
(
 int vtx          //vertex that is newly fixed
 )
{
  queue<int> VertexInQue;                              //vertices to be treated
  VertexInQue.push(vtx);
  while( !VertexInQue.empty() )
    {
      struct vertexstruct *v_ptr=&Vertex[ VertexInQue.front() ];
      VertexInQue.pop();
      //first consider all child vertices of the central vertex j
      struct messagestruct *im_ptr=v_ptr->cim_ptr;
      for(int k=0; k< v_ptr->out_degree; ++k)
	{ 
	  if( (--(im_ptr->v_ptr->active_in_degree)) == 0 )
	    {
	      if( im_ptr->v_ptr->active && 
		  im_ptr->v_ptr->active_inout_degree == 0 )
		{
		  im_ptr->v_ptr->active   = false;
		  im_ptr->v_ptr->occupied = true;
		  --ActiveVertexNumber;
		  VertexInQue.push(im_ptr->v_ptr->index);
		}
	    }
	  ++im_ptr;
	}
      //then consider all parent vertices of the central vertex j
      im_ptr=v_ptr->pim_ptr;
      for(int i=0; i<v_ptr->in_degree; ++i)
	{
	  if( (--(im_ptr->v_ptr->active_out_degree)) ==0 )
	    {
	      if(im_ptr->v_ptr->active &&
		 im_ptr->v_ptr->active_inout_degree==0 )
		{
		  im_ptr->v_ptr->active   = false;
		  im_ptr->v_ptr->occupied = true;
		  --ActiveVertexNumber;
		  VertexInQue.push(im_ptr->v_ptr->index);
		}
	    }
	  ++im_ptr;
	}

      // then consider all the bi-directionally connected neighboring vertices
      if(v_ptr->inout_degree>0)
	{
	  im_ptr=v_ptr->pcim_ptr;
	  for(int l=0; l<v_ptr->inout_degree; ++l)
	    {
	      if( (--(im_ptr->v_ptr->active_inout_degree)) ==0 )
		{
		  if(im_ptr->v_ptr->active &&             
		     ( im_ptr->v_ptr->active_in_degree==0 ||
		       im_ptr->v_ptr->active_out_degree ==0 ) )
		    {
		      im_ptr->v_ptr->active   = false;
		      im_ptr->v_ptr->occupied = true;
		      --ActiveVertexNumber;
		      VertexInQue.push(im_ptr->v_ptr->index);
		    }
		}
	      ++im_ptr;
	    }
	}
    }
  return;
}

//initialize messages
void DFVS::Initialization(void) {
  double uniformprob=1.0e0/(1.0e0+MaxHeight);
  Message=uniformprob;
  return ;
}

/* Belief propagation */
bool DFVS::BeliefPropagation
(
 double error,                                          //convergence criterion
 int count                                       //maximal number of iterations
 )
{
  if(HeightNew>ActiveVertexNumber) HeightNew=ActiveVertexNumber;
  /* first update the list of active vertices */
  for(int index=ActiveVertexNumberOld-1; index >=0; --index)
    {
      if(Vertex[ Permutation[index] ].active == false)
	{
	  struct vertexstruct *v_ptr = &Vertex[ Permutation[index ] ];
	  --ActiveVertexNumberOld;
	  Permutation[index]=Permutation[ActiveVertexNumberOld];
	  /* updating the output messages of the nearest neighbors of the newly
	     inactivated vertex. */
	  struct messagestruct *im_ptr = v_ptr->pim_ptr; //parent input message
	  for(int i=0; i< v_ptr->in_degree; ++i)
	    {                                   //parent vertex is denoted as i
	      if(im_ptr->v_ptr->active) UpdateMessage( im_ptr->v_ptr );
	      ++im_ptr; }
	  im_ptr=v_ptr->cim_ptr;
	  for(int k=0; k < v_ptr->out_degree; ++k)
	    {                                    //child vertex is denoted as k
	      if(im_ptr->v_ptr->active) UpdateMessage( im_ptr->v_ptr );
	      ++im_ptr;
	    }
	  if(v_ptr->active_inout_degree >0)
	    {
	      im_ptr=v_ptr->pcim_ptr;
	      for(int l=0; l< v_ptr->inout_degree; ++l)
		{                                   //BDE neighbor denoted as l
		  if(im_ptr->v_ptr->active)
		    UpdateMessage( im_ptr->v_ptr );
		  ++im_ptr;
		}
	    }
	}
    }
  
  double max_error;           //maximal difference of messages in one iteration
  int iter=1;
  do {
    max_error=0;
    for(int quant=ActiveVertexNumber; quant>=1; --quant)
      {
	int iii = static_cast<int>(quant * u01prn(PRNG) );
	int v = Permutation[iii];
	Permutation[iii] = Permutation[quant-1];
	Permutation[quant-1] = v;
	UpdateMessage( &Vertex[v], max_error );
      }
    if(HeightOld>HeightNew) HeightOld=HeightNew;
    if( iter % 10 == 0)
      cerr<<' '<<max_error;
    ++iter;
  } while ( max_error>error && iter<=count );
  if(max_error<=error)
    {
      cerr<<' '<<max_error<<" :-)\n";
      return true;
    }
  else
    {
      cerr<<' '<<max_error<<" :-(\n";
      return false;
    }
}

/* update the output messages from a vertex. this central vertex j is assumed
   to be active and has at least one active parent vertex and one active child
   vertex, namely v_ptr->active == true &&  (v_ptr->active_in_degree +
   v_ptr->active_inout_degree) >= 1 && (v_ptr->active_out_degree +
   v_ptr->active_inout_degree) >= 1 (if the central vertex does not satisfy
   these properties, it should have already been claimed as v_ptr->active = 
   false in the graph simplification process). */
void DFVS::UpdateMessage
(
 struct vertexstruct *v_ptr,     //vertex whose outgoing messages to be updated
 double& max_prob_diff   //max difference between new and old outgoing messages
 )
{
  int in_degree  =v_ptr->in_degree;
  int out_degree =v_ptr->out_degree;
  valarray<double> in_value(in_degree+1);                    //auxilliary array
  valarray<double> out_value(out_degree+1);                  //auxilliary array
  struct messagestruct *im_ptr = v_ptr->pim_ptr;         //parent input message
  for(int i=0; i<in_degree; ++i) {       //parent vertex is denoted as i,i2,...
    if(im_ptr->v_ptr->active) {
      IsActiveParent[i]=true;
      NewMessageToParent[i*(MaxHeight+1)]=Weight0;           //q_{j->i}^{h_j=0}
      NormalizationOPM[i]=Weight0;
      if(HeightOld==HeightNew) VectorIn[i]=im_ptr->q_cavity[0];
      else {
        double zeroval=im_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
          zeroval += im_ptr->q_cavity[he];
        VectorIn[i]=zeroval; } }
    else IsActiveParent[i]=false;
    ++im_ptr; }
  im_ptr = v_ptr->cim_ptr;                                //child input message
  for(int k=0; k<out_degree; ++k) {       //child vertex is denoted as k,k2,...
    if(im_ptr->v_ptr->active) {
      IsActiveChild[k]=true;
      NewMessageToChild[k*(MaxHeight+1)]=Weight0;            //q_{j->k}^{h_j=0}
      NormalizationOCM[k]=Weight0;
      if(HeightOld==HeightNew) VectorOut[k]=im_ptr->q_cavity[0];
      else {
        double zeroval=im_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
          zeroval += im_ptr->q_cavity[he];
        VectorOut[k]= zeroval; } }
    else IsActiveChild[k]=false;
    ++im_ptr; }
  if(v_ptr->active_inout_degree == 0) {        //no active bi-directional edges
    for(int h=1; h<=HeightNew; ++h) {
      in_value=1.0e0;
      im_ptr=v_ptr->pim_ptr;                             //parent input message
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  in_value[in_degree] *= VectorIn[i];  //all the arcs [i,j] to vertex j
	  for(int i2=0; i2<in_degree; ++i2) {
	    if(IsActiveParent[i2] && i2 != i) in_value[i2] *= VectorIn[i]; }
	  VectorIn[i] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i])
	  NewMessageToParent[i*(MaxHeight+1)+h] =in_value[i]; }
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k])
	  NewMessageToChild[k*(MaxHeight+1)+h]=in_value[in_degree]; } }
    for(int h=HeightNew; h>=1; --h) {
      out_value=1.0e0;
      im_ptr=v_ptr->cim_ptr;
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  out_value[out_degree] *= VectorOut[k];
	  for(int k2=0; k2<out_degree; ++k2) {
	    if(IsActiveChild[k2] && k2 != k) out_value[k2] *= VectorOut[k]; }
	  VectorOut[k] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  int position=k*(MaxHeight+1)+h;
	  NewMessageToChild[position] *= out_value[k];
	  NormalizationOCM[k] += NewMessageToChild[position]; } }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  int position=i*(MaxHeight+1)+h;
	  NewMessageToParent[position] *= out_value[out_degree];
	  NormalizationOPM[i] += NewMessageToParent[position]; } } } }
  else { //v_ptr->active_inout_degree > 0,  having active bi-directional edges
    int inout_degree=v_ptr->inout_degree;
    valarray<double> inout_value(1.0e0, inout_degree+1);     //auxilliary array
    im_ptr=v_ptr->pcim_ptr;                                 //BDE input message
    for(int l=0; l<inout_degree; ++l) { //bi-dir neighbor denoted as l, l2, ...
      if(im_ptr->v_ptr->active) {
	IsActiveBDE[l]=true;
	NewMessageToBDE[l*(MaxHeight+1)]=Weight0;
	NormalizationBDE[l]=Weight0;
        if(HeightOld==HeightNew) VectorBDE[l]=im_ptr->q_cavity[0];
        else {
	  double zeroval=im_ptr->q_cavity[0];
	  for(int he=HeightOld; he>HeightNew; --he)
	    zeroval += im_ptr->q_cavity[he];
	  VectorBDE[l]=zeroval; } }
      else IsActiveBDE[l]=false;
      ++im_ptr; }
    for(int l=0; l<inout_degree; ++l) {
      if(IsActiveBDE[l]) {
	inout_value[inout_degree] *= VectorBDE[l];
	for(int l2=0; l2<inout_degree; ++l2) {
	  if(IsActiveBDE[l2] && l2 != l) inout_value[l2] *= VectorBDE[l]; } } }
    for(int h=1; h<=HeightNew; ++h) {
      in_value=1.0e0;
      im_ptr=v_ptr->pim_ptr;                             //parent input message
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  in_value[in_degree] *= VectorIn[i];  //all the arcs [i,j] to vertex j
	  for(int i2=0; i2<in_degree; ++i2) {
	    if(IsActiveParent[i2] && i2 != i) in_value[i2] *= VectorIn[i]; }
	  VectorIn[i] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i])
	  NewMessageToParent[i*(MaxHeight+1)+h]
	    =in_value[i]*inout_value[inout_degree]; }
      double wvalue=in_value[in_degree]*inout_value[inout_degree];
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) NewMessageToChild[k*(MaxHeight+1)+h]=wvalue; }
      for(int l=0; l<inout_degree; ++l) {
	if(IsActiveBDE[l])
	  NewMessageToBDE[l*(MaxHeight+1)+h] 
	    = in_value[in_degree]*inout_value[l]; } }
    for(int h=HeightNew; h>=1; --h) {
      out_value=1.0e0;
      im_ptr=v_ptr->cim_ptr;
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  out_value[out_degree] *= VectorOut[k];
	  for(int k2=0; k2<out_degree; ++k2) {
	    if(IsActiveChild[k2] && k2 != k) out_value[k2] *= VectorOut[k]; }
	  VectorOut[k] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  int position=k*(MaxHeight+1)+h;
	  NewMessageToChild[position] *= out_value[k];
	  NormalizationOCM[k] += NewMessageToChild[position]; } }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  int position=i*(MaxHeight+1)+h;
	  NewMessageToParent[position] *= out_value[out_degree];
	  NormalizationOPM[i] += NewMessageToParent[position]; } }
      for(int l=0; l<inout_degree; ++l) {
	if(IsActiveBDE[l]) {
	  int position = l * (MaxHeight+1)+h;
	  NewMessageToBDE[position ] *= out_value[out_degree];
	  NormalizationBDE[l] += NewMessageToBDE[position]; } } }
    struct simplifiedmessagestruct *opcm_ptr = v_ptr->opcm_ptr;
    for(int l=0; l<inout_degree; ++l) {
      if(IsActiveBDE[l]) {
	int position=l*(MaxHeight+1);
	double diff=0;
        if(HeightOld>HeightNew) {
	  double zeroval=opcm_ptr->q_cavity[0];
	  for(int he=HeightOld; he>HeightNew; --he)
	    zeroval += opcm_ptr->q_cavity[he];
	  opcm_ptr->q_cavity[0]=zeroval; }
	for(int h=0; h<=HeightNew; ++h) {
	  NewMessageToBDE[position] /= NormalizationBDE[l];
	  diff += abs(NewMessageToBDE[position] - opcm_ptr->q_cavity[h]);
	  opcm_ptr->q_cavity[h] = WeightNew*NewMessageToBDE[position]
	    +WeightOld * (opcm_ptr->q_cavity[h]);
	  ++position; }
	if(diff>max_prob_diff) max_prob_diff=diff; }
      ++opcm_ptr; } }
  struct simplifiedmessagestruct *om_ptr = v_ptr->opm_ptr;
  for(int i=0; i<in_degree; ++i) {
    if(IsActiveParent[i]) {
      int position=i*(MaxHeight+1);
      double diff=0;
      if(HeightOld>HeightNew) {
        double zeroval=om_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
	  zeroval += om_ptr->q_cavity[he];
        om_ptr->q_cavity[0]=zeroval; }
      for(int h=0; h<=HeightNew; ++h) {
	NewMessageToParent[position] /= NormalizationOPM[i];
	diff += abs(NewMessageToParent[position] - om_ptr->q_cavity[h]);
	om_ptr->q_cavity[h] = WeightNew*NewMessageToParent[position]
	  +WeightOld * (om_ptr->q_cavity[h]);
	++position; }
      if(diff>max_prob_diff) max_prob_diff=diff; }
    ++om_ptr; }
  om_ptr = v_ptr->ocm_ptr;
  for(int k=0; k<out_degree; ++k) {
    if(IsActiveChild[k]) {
      int position=k*(MaxHeight+1);
      double diff=0;
      if(HeightOld>HeightNew) {
        double zeroval=om_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
	  zeroval += om_ptr->q_cavity[he];
        om_ptr->q_cavity[0]=zeroval; }
      for(int h=0; h<=HeightNew; ++h) {
	NewMessageToChild[position] /= NormalizationOCM[k];
	diff += abs(NewMessageToChild[position] - om_ptr->q_cavity[h]);
	om_ptr->q_cavity[h] = WeightNew*NewMessageToChild[position]
	  +WeightOld * (om_ptr->q_cavity[h]);
	++position; }
      if(diff>max_prob_diff) max_prob_diff=diff; }
    ++om_ptr; }
  return;
}

/* update the output messages from a vertex. this central vertex j is assumed
   to be active and has at least one active parent vertex and one active child
   vertex, namely v_ptr->active == true &&  (v_ptr->active_in_degree +
   v_ptr->active_inout_degree) >= 1 && (v_ptr->active_out_degree +
   v_ptr->active_inout_degree) >= 1 (if the central vertex does not satisfy
   these properties, it should have already been claimed as v_ptr->active = 
   false in the graph simplification process). */
void DFVS::UpdateMessage
(
 struct vertexstruct *v_ptr     //vertex whose outgoing messages to be updated
 )
{
  int in_degree  =v_ptr->in_degree;
  int out_degree =v_ptr->out_degree;
  valarray<double> in_value(in_degree+1);                    //auxilliary array
  valarray<double> out_value(out_degree+1);                  //auxilliary array

  struct messagestruct *im_ptr = v_ptr->pim_ptr;         //parent input message
  for(int i=0; i<in_degree; ++i) {       //parent vertex is denoted as i,i2,...
    if(im_ptr->v_ptr->active) {
      IsActiveParent[i]=true;
      NewMessageToParent[i*(MaxHeight+1)]=Weight0;           //q_{j->i}^{h_j=0}
      NormalizationOPM[i]=Weight0;
      if(HeightOld==HeightNew) VectorIn[i]=im_ptr->q_cavity[0];
      else {
        double zeroval=im_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
	  zeroval += im_ptr->q_cavity[he];
        VectorIn[i]=zeroval; } }
    else IsActiveParent[i]=false;
    ++im_ptr; }
  im_ptr = v_ptr->cim_ptr;                                //child input message
  for(int k=0; k<out_degree; ++k) {       //child vertex is denoted as k,k2,...
    if(im_ptr->v_ptr->active) {
      IsActiveChild[k]=true;
      NewMessageToChild[k*(MaxHeight+1)]=Weight0;            //q_{j->k}^{h_j=0}
      NormalizationOCM[k]=Weight0;
      if(HeightOld==HeightNew) VectorOut[k]=im_ptr->q_cavity[0];
      else {
        double zeroval=im_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
	  zeroval += im_ptr->q_cavity[he];
        VectorOut[k]= zeroval; } }
    else IsActiveChild[k]=false;
    ++im_ptr; }
  if(v_ptr->active_inout_degree == 0) {        //no active bi-directional edges
    for(int h=1; h<=HeightNew; ++h) {
      in_value=1.0e0;
      im_ptr=v_ptr->pim_ptr;                             //parent input message
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  in_value[in_degree] *= VectorIn[i];  //all the arcs [i,j] to vertex j
	  for(int i2=0; i2<in_degree; ++i2) {
	    if(IsActiveParent[i2] && i2 != i) in_value[i2] *= VectorIn[i]; }
	  VectorIn[i] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i])
	  NewMessageToParent[i*(MaxHeight+1)+h] =in_value[i]; }
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k])
	  NewMessageToChild[k*(MaxHeight+1)+h]=in_value[in_degree]; } }
    for(int h=HeightNew; h>=1; --h) {
      out_value=1.0e0;
      im_ptr=v_ptr->cim_ptr;
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  out_value[out_degree] *= VectorOut[k];
	  for(int k2=0; k2<out_degree; ++k2) {
	    if(IsActiveChild[k2] && k2 != k) out_value[k2] *= VectorOut[k]; }
	  VectorOut[k] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  int position=k*(MaxHeight+1)+h;
	  NewMessageToChild[position] *= out_value[k];
	  NormalizationOCM[k] += NewMessageToChild[position]; } }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  int position=i*(MaxHeight+1)+h;
	  NewMessageToParent[position] *= out_value[out_degree];
	  NormalizationOPM[i] += NewMessageToParent[position]; } } } }
  else { //v_ptr->active_inout_degree > 0,  having active bi-directional edges
    int inout_degree=v_ptr->inout_degree;
    valarray<double> inout_value(1.0e0, inout_degree+1);     //auxilliary array
    im_ptr=v_ptr->pcim_ptr;                                 //BDE input message
    for(int l=0; l<inout_degree; ++l) { //bi-dir neighbor denoted as l, l2, ...
      if(im_ptr->v_ptr->active) {
	IsActiveBDE[l]=true;
	NewMessageToBDE[l*(MaxHeight+1)]=Weight0;
	NormalizationBDE[l]=Weight0;
	if(HeightOld==HeightNew) VectorBDE[l]=im_ptr->q_cavity[0];
        else { 
	  double zeroval=im_ptr->q_cavity[0];
	  for(int he=HeightOld; he>HeightNew; --he)
	    zeroval += im_ptr->q_cavity[he];
	  VectorBDE[l]=zeroval; } }
      else IsActiveBDE[l]=false;
      ++im_ptr; }
    for(int l=0; l<inout_degree; ++l) {
      if(IsActiveBDE[l]) {
	inout_value[inout_degree] *= VectorBDE[l];
	for(int l2=0; l2<inout_degree; ++l2) {
	  if(IsActiveBDE[l2] && l2 != l) inout_value[l2] *= VectorBDE[l]; } } }
    for(int h=1; h<=HeightNew; ++h) {
      in_value=1.0e0;
      im_ptr=v_ptr->pim_ptr;                             //parent input message
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  in_value[in_degree] *= VectorIn[i];  //all the arcs [i,j] to vertex j
	  for(int i2=0; i2<in_degree; ++i2) {
	    if(IsActiveParent[i2] && i2 != i) in_value[i2] *= VectorIn[i]; }
	  VectorIn[i] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i])
	  NewMessageToParent[i*(MaxHeight+1)+h]
	    =in_value[i]*inout_value[inout_degree]; }
      double wvalue=in_value[in_degree]*inout_value[inout_degree];
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) NewMessageToChild[k*(MaxHeight+1)+h]=wvalue; }
      for(int l=0; l<inout_degree; ++l) {
	if(IsActiveBDE[l])
	  NewMessageToBDE[l*(MaxHeight+1)+h] 
	    = in_value[in_degree]*inout_value[l]; } }
    for(int h=HeightNew; h>=1; --h) {
      out_value=1.0e0;
      im_ptr=v_ptr->cim_ptr;
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  out_value[out_degree] *= VectorOut[k];
	  for(int k2=0; k2<out_degree; ++k2) {
	    if(IsActiveChild[k2] && k2 != k) out_value[k2] *= VectorOut[k]; }
	  VectorOut[k] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  int position=k*(MaxHeight+1)+h;
	  NewMessageToChild[position] *= out_value[k];
	  NormalizationOCM[k] += NewMessageToChild[position]; } }
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  int position=i*(MaxHeight+1)+h;
	  NewMessageToParent[position] *= out_value[out_degree];
	  NormalizationOPM[i] += NewMessageToParent[position]; } }
      for(int l=0; l<inout_degree; ++l) {
	if(IsActiveBDE[l]) {
	  int position = l * (MaxHeight+1)+h;
	  NewMessageToBDE[position ] *= out_value[out_degree];
	  NormalizationBDE[l] += NewMessageToBDE[position]; } } }
    struct simplifiedmessagestruct *opcm_ptr = v_ptr->opcm_ptr;
    for(int l=0; l<inout_degree; ++l) {
      if(IsActiveBDE[l]) {
	int position=l*(MaxHeight+1);
	if(HeightOld>HeightNew) {
	  double zeroval=opcm_ptr->q_cavity[0];
	  for(int he=HeightOld; he>HeightNew; --he)
	    zeroval += opcm_ptr->q_cavity[he];
	  opcm_ptr->q_cavity[0]=zeroval; }
	for(int h=0; h<=HeightNew; ++h) {
	  NewMessageToBDE[position] /= NormalizationBDE[l];
	  opcm_ptr->q_cavity[h] = WeightNew*NewMessageToBDE[position]
	    +WeightOld * (opcm_ptr->q_cavity[h]);
	  ++position; } }
      ++opcm_ptr; } }
  struct simplifiedmessagestruct *om_ptr = v_ptr->opm_ptr;
  for(int i=0; i<in_degree; ++i) {
    if(IsActiveParent[i]) {
      int position=i*(MaxHeight+1);
      if(HeightOld>HeightNew) {
        double zeroval=om_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
	  zeroval += om_ptr->q_cavity[he];
        om_ptr->q_cavity[0]=zeroval; }
      for(int h=0; h<=HeightNew; ++h) {
	NewMessageToParent[position] /= NormalizationOPM[i];
	om_ptr->q_cavity[h] = WeightNew*NewMessageToParent[position]
	  +WeightOld * (om_ptr->q_cavity[h]);
	++position; } }
    ++om_ptr; }
  om_ptr = v_ptr->ocm_ptr;
  for(int k=0; k<out_degree; ++k) {
    if(IsActiveChild[k]) {
      int position=k*(MaxHeight+1);
      if(HeightOld>HeightNew) {
        double zeroval=om_ptr->q_cavity[0];
        for(int he=HeightOld; he>HeightNew; --he)
	  zeroval += om_ptr->q_cavity[he];
        om_ptr->q_cavity[0]=zeroval; }
      for(int h=0; h<=HeightNew; ++h) {
	NewMessageToChild[position] /= NormalizationOCM[k];
	om_ptr->q_cavity[h] = WeightNew*NewMessageToChild[position]
	  +WeightOld * (om_ptr->q_cavity[h]);
	++position; } }
    ++om_ptr; }
  return;
}

/* thermodynamic quantities */
void DFVS::Thermodynamics
(
  string& bpfname                //output file of thermodynamic quantities
 ) {
  if(ActiveVertexNumber == 0) return ;
  MaxFixNumber=static_cast<int>( ActiveVertexNumber*FixFraction );
  if(MaxFixNumber==0) MaxFixNumber = 1;
  CandidateSize=0; /* number of candidate vertices in each q_0 bin of width
		      0.01. this value is initialized to be zero. */
  MinRange=0; /* the maximal value of r0 such that sum_{r>=r0} CandidateSize[r]
		 >= MaxFixNumber is satisfied */
  double phi_vtx=0;                       //vertex contribution to free entropy
  double phi_edge=0;                        //edge contribution to free entropy
  double phi_BDE =0;                          //BDE contribution to free energy
  double rho_vtx=0;                                 //vertex occupation density
  valarray<double> in_value(MaxHeight+1);                    //auxilliary array
  valarray<double> normEdge(Max_in_degree);      //used in calculating phi_edge
  for(int v_index=0; v_index<ActiveVertexNumber; ++v_index) {
    struct vertexstruct *v_ptr = &Vertex[ Permutation[v_index] ];
    /*  the central vertex is assumed to have at least one active parent and
	one active child. if it has zero or only one active parent/child, it
	it should have been removed from the active subgraph */
    int in_degree  =v_ptr->in_degree;
    int out_degree =v_ptr->out_degree;
    int inout_degree = v_ptr->inout_degree;
    double norm=Weight0;
    struct messagestruct *im_ptr = v_ptr->pim_ptr;       //parent input message
    struct simplifiedmessagestruct *om_ptr
      = v_ptr->opm_ptr;                              //output message to parent
    for(int i=0; i<in_degree; ++i) {     //parent vertex is denoted as i,i2,...
      if(im_ptr->v_ptr->active) {
	IsActiveParent[i]=true;
	normEdge[i]=om_ptr->q_cavity[0];
	VectorIn[i]=im_ptr->q_cavity[0]; }
      else IsActiveParent[i]=false;
      ++im_ptr;
      ++om_ptr; }
    im_ptr = v_ptr->cim_ptr;                              //child input message
    for(int k=0; k<out_degree; ++k) {     //child vertex is denoted as k,k2,...
      if(im_ptr->v_ptr->active) {
	IsActiveChild[k]=true;
	VectorOut[k]= im_ptr->q_cavity[0]; }
      else IsActiveChild[k]=false;
      ++im_ptr; }
    double inout_value=1.0e0;
    if(v_ptr->active_inout_degree>0) { 
      //having active bi-directional connected neighbors
      im_ptr=v_ptr->pcim_ptr;                          //bi-dir edge in-message
      om_ptr=v_ptr->opcm_ptr;                         //bi-dir edge out-message
      for(int l=0; l<inout_degree; ++l) {        //bi-dir neighbor denoted as l
	if(im_ptr->v_ptr->active) {
	  inout_value *= im_ptr->q_cavity[0];
	  if(im_ptr->v_ptr->index > v_ptr->index) 
	    phi_BDE += log(im_ptr->q_cavity[0] + (1.0e0 - im_ptr->q_cavity[0]) 
			   * om_ptr->q_cavity[0] ); }
	++im_ptr;
	++om_ptr;  }  }
    for(int h=1; h<=HeightNew; ++h) {
      in_value[h]=1.0e0;
      im_ptr=v_ptr->pim_ptr;                             //parent input message
      om_ptr=v_ptr->opm_ptr;                         //output message to parent
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  in_value[h] *= VectorIn[i];
	  normEdge[i] += VectorIn[i]* om_ptr->q_cavity[h];
	  VectorIn[i] += im_ptr->q_cavity[h]; }
	++im_ptr;
	++om_ptr;  } }
    for(int h=HeightNew; h>=1; --h) {
      double out_value=1.0e0;
      im_ptr=v_ptr->cim_ptr;
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  out_value *= VectorOut[k];
	  VectorOut[k] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      norm += inout_value*in_value[h]*out_value; }
    for(int i=0; i<in_degree; ++i) {
      if(IsActiveParent[i]) phi_edge += log(normEdge[i]);  }
    v_ptr->q_0 = Weight0/norm;
    phi_vtx += X+log(norm);
    rho_vtx += 1.0e0-v_ptr->q_0;
    int rrr = static_cast<int>(v_ptr->q_0 * 100);
    if(rrr>=MinRange) {
      if(CandidateSize[rrr]<MaxFixNumber) {
	CandidateVertices[rrr*MaxFixNumber + CandidateSize[rrr] ] 
	  = v_ptr->index;
	++CandidateSize[rrr]; }
      else MinRange=rrr; } }
  phi_vtx  /= ActiveVertexNumber;
  phi_edge /= ActiveVertexNumber;
  phi_BDE  /= ActiveVertexNumber;
  rho_vtx  /= ActiveVertexNumber;
  double phi = (phi_vtx - phi_edge-phi_BDE)/X;
  ofstream output(bpfname.c_str(), ios_base::app );
  output<< X << '\t'
	<< ActiveVertexNumber << '\t'
	<< rho_vtx <<'\t'
	<< phi <<'\t'
	<< X * (phi-rho_vtx) << endl;
  output.close();
  return ;
}



/* calculate the unoccupation probability of each vertex and rank the vertices
   accordingly */
void DFVS::VertexRanking(void) {
  if(ActiveVertexNumber == 0) return;
  MaxFixNumber=static_cast<int>(ActiveVertexNumber*FixFraction);
  if(MaxFixNumber==0) MaxFixNumber = 1;
  CandidateSize=0; /* number of candidate vertices in each q_0 bin of width
		      0.01. this value is initialized to be zero. */
  MinRange=0; /* the maximal value of r0 such that sum_{r>=r0} CandidateSize[r]
		 >= MaxFixNumber is satisfied */
  valarray<double> in_value(MaxHeight+1);                     //auxiliary array
  for(int v_index=0; v_index < ActiveVertexNumber; ++v_index) {
    struct vertexstruct *v_ptr = &Vertex[ Permutation[v_index] ];
    /*  the central vertex is assumed to have at least one active parent and
	one active child. if it has zero or only one active parent/child, it
	it should have been removed from the active subgraph */
    int in_degree  =v_ptr->in_degree;
    int out_degree =v_ptr->out_degree;
    double norm=Weight0;
    struct messagestruct *im_ptr = v_ptr->pim_ptr;       //parent input message
    for(int i=0; i<in_degree; ++i) {     //parent vertex is denoted as i,i2,...
      if(im_ptr->v_ptr->active) {
	IsActiveParent[i]=true;
	VectorIn[i]=im_ptr->q_cavity[0];  }
      else IsActiveParent[i]=false;
      ++im_ptr; }
    im_ptr = v_ptr->cim_ptr;                              //child input message
    for(int k=0; k<out_degree; ++k) {     //child vertex is denoted as k,k2,...
      if(im_ptr->v_ptr->active) {
	IsActiveChild[k]=true;
	VectorOut[k]= im_ptr->q_cavity[0];  }
      else IsActiveChild[k]=false;
      ++im_ptr; }
    double inout_value=1.0e0;
    if(v_ptr->active_inout_degree >0) { //with active BDE neighbors
      im_ptr=v_ptr->pcim_ptr;                          //bi-dir edge in-message
      for(int l=0; l<v_ptr->inout_degree; ++l) { //bi-dir neighbor denoted as l
	if(im_ptr->v_ptr->active) inout_value *= im_ptr->q_cavity[0];
	++im_ptr; } }
    for(int h=1; h<=HeightNew; ++h) {
      in_value[h]=1.0e0;
      im_ptr=v_ptr->pim_ptr;                             //parent input message
      for(int i=0; i<in_degree; ++i) {
	if(IsActiveParent[i]) {
	  in_value[h] *= VectorIn[i];
	  VectorIn[i] += im_ptr->q_cavity[h]; }
	++im_ptr; } }
    for(int h=HeightNew; h>=1; --h) {
      double out_value=1.0e0;
      im_ptr=v_ptr->cim_ptr;
      for(int k=0; k<out_degree; ++k) {
	if(IsActiveChild[k]) {
	  out_value *= VectorOut[k];
	  VectorOut[k] += im_ptr->q_cavity[h]; }
	++im_ptr; }
      norm += inout_value* in_value[h] * out_value; }
    v_ptr->q_0 = Weight0/norm;
    int rrr = static_cast<int>(v_ptr->q_0 * 100);
    if(rrr>=MinRange) {
      if(CandidateSize[rrr]<MaxFixNumber) {
	CandidateVertices[rrr*MaxFixNumber + CandidateSize[rrr] ] 
	  = v_ptr->index;
	++CandidateSize[rrr]; }
      else MinRange=rrr; } }
  return ;
}



//externally fixing some variables to be empty and simplify the system
int DFVS::Fix0(void) {
  int rank = 100; /* q_0 is distributed to 101 bins [0,0.01), [0.01,0.02), ...,
		     [0.99,1), [1,1] */
  int num_fixed_empty = 0;  //number of vertices externally fixed to unoccupied
  double mean_emptyprob = 0;        // ... and mean q_0 value of these vertices
  int num_examined_vertices = 0;        //number of examined candidate vertices
  while(num_examined_vertices < MaxFixNumber && ActiveVertexNumber > 0) {
    int csize=CandidateSize[rank];
    if(csize>0) {
      int *i_ptr = &CandidateVertices[rank * MaxFixNumber];
      if((num_examined_vertices + csize) <= MaxFixNumber) {
	num_examined_vertices += csize;
	for(int s=0; s<csize; ++s) {
	  int vtx = *i_ptr;
	  if(Vertex[vtx].active) {
	    ++num_fixed_empty;
	    Vertex[vtx].active = false;
	    Vertex[vtx].occupied = false;
	    --ActiveVertexNumber;
	    mean_emptyprob += Vertex[vtx].q_0;
	    Simplify(vtx); }
	  ++i_ptr; } }
      else {
	while(num_examined_vertices < MaxFixNumber) {
	  ++num_examined_vertices;
	  int vtx = *i_ptr;
	  if(Vertex[vtx].active) {
	    Vertex[vtx].active = false;
	    Vertex[vtx].occupied = false;
	    ++num_fixed_empty;
	    --ActiveVertexNumber;
	    mean_emptyprob += Vertex[vtx].q_0;
	    Simplify(vtx) ; }
	  ++i_ptr; } } }
    --rank; }
  return ActiveVertexNumber;
}

/* check whether the final occupation pattern corresponds to a DFVS */
bool DFVS::CheckDFVS
(
  string& dfvsfname    //the occupation pattern 0: belongs to DFVS; 1: not
 )
{
  int numOccupy =0;                            //number of vertices not in DFVS
  int numEmpty =0;                                //number of vertices in DFVS
  int numActiveEdge =0;           //number of arcs within the active subgraph
  queue<int> ParentFreeVertices;                     //vertices without parents
  bool noconflict=true;
  struct vertexstruct *v_ptr = &Vertex[1];
  for(int v=1; v<=VertexNumber && noconflict; ++v)
    {
      if(v_ptr->occupied  == false) ++numEmpty;          //vertex belong to DFVS
      else
	{                                                 //vertex not deleted
	  ++numOccupy;
	  int active_in_degree=0;
	  struct messagestruct *pim_ptr = v_ptr->pim_ptr;//pure parent vertices
	  for(int i=0; i<v_ptr->in_degree; ++i)
	    {
	      if(pim_ptr->v_ptr->occupied) {
		++active_in_degree;
		++numActiveEdge;
	      }
	      ++pim_ptr;
	    }
	  v_ptr->active_in_degree=active_in_degree; //number of active parents
	  if(active_in_degree == 0) ParentFreeVertices.push(v_ptr->index);
	  if(v_ptr->inout_degree>0)
	    {
	      struct messagestruct *pcim_ptr=v_ptr->pcim_ptr;
	      //bi-directional edges
	      for(int i=0; i<v_ptr->inout_degree && noconflict; ++i)
		{
		  if(pcim_ptr->v_ptr->occupied) noconflict=false;
		  ++pcim_ptr;
		}
	    }
	}
      ++v_ptr;
    }
  if(noconflict==false)
    return false; //directed cycle due to bi-directional edges
  
  /* if there is an arc connecting a vertex to itself, the involved vertex
     should have been set to be unoccupied when the digraph was read into.
     we therefore did not check for this special case. */
  
  int ResidueSize=numOccupy;
  while( !ParentFreeVertices.empty() )
    {
      v_ptr = &Vertex[ ParentFreeVertices.front() ];
      ParentFreeVertices.pop();
      --ResidueSize;
      struct messagestruct *cim_ptr = v_ptr->cim_ptr;
      for(int k=0; k<v_ptr->out_degree; ++k)
	{
	  if(cim_ptr->v_ptr->occupied)
	    {
	      if( (--(cim_ptr->v_ptr->active_in_degree)) == 0)
		ParentFreeVertices.push(cim_ptr->v_ptr->index);
	    }
	  ++cim_ptr;
	}
    }
  
  if(ResidueSize >0)
    {
      cerr<<"Not a proper DFVS. The size of the residue cyclic subgraph is "
	  <<ResidueSize<<endl;
      cerr.flush();
      return false;
    }
  
  ofstream optfile(dfvsfname.c_str() );
  optfile << numOccupy<<'\t'<<numEmpty<<'\t'<<numActiveEdge<<endl<<endl;
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v)
    {
      optfile << v_ptr->index <<'\t' << (v_ptr->occupied? 1 : 0)<<endl;
      ++v_ptr;
    }
  optfile.close();
  return true;
}

/* check whether the final occupation pattern corresponds to a DFVS and, if
   yes then refine the DFVS pattern. An attempt is made to remove each of the
   vertices of the constructed DFVS set. If this removal does not cause
   conflict then this removal is accepted and the DFVS set size decreases by
   one. */
bool DFVS::RefineDFVS
(
  string& dfvsfname  //the occupation pattern.  0: belongs to DFVS; 1: not
 )
{
  bool noconflict=true;
  int numOccupy =0;                            //number of vertices not in DFVS
  int numEmpty =0;                                 //number of vertices in DFVS
  int numActiveEdge =0;             //number of arcs within the active subgraph
  queue<int> ParentFreeVertices;                     //vertices without parents
  struct vertexstruct *v_ptr = &Vertex[1];
  for(int v=1; v<=VertexNumber && noconflict; ++v)
    {
      if(v_ptr->occupied  == false) ++numEmpty;        //vertex belong to DFVS
      else
	{                                                //vertex not deleted
	  ++numOccupy;
	  v_ptr->active_in_degree = 0;
	  struct messagestruct *pim_ptr = v_ptr->pim_ptr;//pure parent vertices
	  for(int i=0; i<v_ptr->in_degree; ++i)
	    {
	      if(pim_ptr->v_ptr->occupied)
		{
		  ++(v_ptr->active_in_degree);
		  ++numActiveEdge; }
	      ++pim_ptr;
	    }
	  if(v_ptr->active_in_degree == 0)
	    ParentFreeVertices.push(v_ptr->index);
	  if(v_ptr->inout_degree >0)
	    {
	      struct messagestruct *pcim_ptr=v_ptr->pcim_ptr; 
	      //bi-directional edges
	      for(int l=0; l<v_ptr->inout_degree && noconflict; ++l)
		{
		  if(pcim_ptr->v_ptr->occupied) noconflict=false;
		  ++pcim_ptr;
		}
	    }
	}
      ++v_ptr;
    }
  if(noconflict==false)
    {
      cerr<<"Not a proper DFVS. Directed cycle due to bi-directional edges.\n";
      return false;
    }                //directed cycle due to bi-directional edges
  
  /* if there is an arc connecting a vertex to itself, the involved vertex
     should have been set to be unoccupied when the digraph was read into.
     we therefore did not check for this special case. */
  
  int numOccupyNew=numOccupy;
  int numEmptyNew=numEmpty;
  int numActiveEdgeNew=numActiveEdge;
  valarray<int> DFVS(numEmpty);
  v_ptr = &Vertex[1];
  numEmpty=0;
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->occupied == false)
	DFVS[numEmpty++ ] = v_ptr->index;
      ++v_ptr;
    }
  int ResidueSize=numOccupy;
  while( !ParentFreeVertices.empty() )
    {
      v_ptr = &Vertex[ ParentFreeVertices.front() ];
      ParentFreeVertices.pop();
      --ResidueSize;
      struct messagestruct *cim_ptr = v_ptr->cim_ptr;
      for(int k=0; k<v_ptr->out_degree; ++k)
	{
	  if(cim_ptr->v_ptr->occupied) {
	    if( (--(cim_ptr->v_ptr->active_in_degree)) == 0)
	      ParentFreeVertices.push(cim_ptr->v_ptr->index);
	  }
	  ++cim_ptr;
	}
    }
  if(ResidueSize >0)
    {
      cerr<<"Not a proper DFVS. The size of the residue cyclic subgraph is "
	  <<ResidueSize<<endl;
      return false;
    }
  cout << numOccupy  <<'\t'
       << numEmpty  <<'\t'
       << numActiveEdge
       << endl;
  //                        check if each vertex of the DFVS set can be deleted
  int reduced=1;                            //where the size of DFVS is reduced
  while(reduced>0)  //repeat refinement if the DFVS was reduced in ealier round
    {
      reduced=0;
      for(int quant =numEmptyNew; quant >=1; --quant)
	{
	  int iii = static_cast<int>( quant * u01prn(PRNG) );
	  int v = DFVS[iii];
	  DFVS[iii] = DFVS[quant-1];
	  DFVS[quant-1]=v;
	  if(SelfArcVertices.find(v) != SelfArcVertices.end() )
	    continue;
	  noconflict=true;
	  v_ptr=&Vertex[v];
	  if(v_ptr->inout_degree >0)
	    {
	      struct messagestruct *pcim_ptr = v_ptr->pcim_ptr;
	      for(int l=0; l<v_ptr->inout_degree && noconflict; ++l)
		{
		  if(pcim_ptr->v_ptr->occupied) noconflict=false;
		  ++pcim_ptr;
		}
	    }
	  if(noconflict==false)
	    continue;
	  set<int> oparents;         //the set of occupied parent vertices of v
	  struct messagestruct *im_ptr = v_ptr->pim_ptr;
	  for(int i=0; i<v_ptr->in_degree; ++i)
	    {
	      if(im_ptr->v_ptr->occupied)
		oparents.insert(im_ptr->v_ptr->index);
	      ++im_ptr;
	    }
	  stack<int> candidates;
	  set<int> odescendents;                         //occupied descendents
	  candidates.push(v);
	  while( !candidates.empty() && noconflict)
	    {
	      int v2=candidates.top();
	      candidates.pop();
	      im_ptr=Vertex[ v2 ].cim_ptr;
	      for(int k=0; k<Vertex[v2].out_degree && noconflict; ++k)
		{
		  if(im_ptr->v_ptr->occupied)
		    {
		      int v3 = im_ptr->v_ptr->index;
		      if(odescendents.find(v3) == odescendents.end() )
			{
			  odescendents.insert(v3);
			  candidates.push(v3);
			  if(oparents.find(v3) != oparents.end() )
			    noconflict = false;
			}
		    }
		  ++im_ptr;
		}
	    }
	  if(noconflict==false)
	    continue;
	  Vertex[v].occupied=true;
	  ++numOccupyNew;
	  --numEmptyNew;
	  ++reduced;
	  DFVS[quant-1]=DFVS[numEmptyNew];
	  DFVS[numEmptyNew]=v;
	  im_ptr=v_ptr->pim_ptr;
	  for(int i=0; i<v_ptr->in_degree; ++i)
	    {
	      if(im_ptr->v_ptr->occupied)
		++numActiveEdgeNew;
	      ++im_ptr;
	    }
	  im_ptr=v_ptr->cim_ptr;
	  for(int k=0; k<v_ptr->out_degree; ++k)
	    {
	      if(im_ptr->v_ptr->occupied)
		++numActiveEdgeNew;
	      ++im_ptr;
	    }
	}
      cout<<"reduced: "<<reduced<<endl;
    }
  //check if the DFVS set after refinement is still a true DFVS set
  noconflict=true;
  ResidueSize=0;
  while( !ParentFreeVertices.empty() ) ParentFreeVertices.pop();
  v_ptr = &Vertex[1];
  for(int v=1; v<=VertexNumber && noconflict; ++v)
    {
      if(v_ptr->occupied)
	{
	  ++ResidueSize;
	  v_ptr->active_in_degree = 0;
	  struct messagestruct *pim_ptr= v_ptr->pim_ptr; //pure parent vertices
	  for(int i=0; i<v_ptr->in_degree; ++i)
	    {
	      if(pim_ptr->v_ptr->occupied) ++(v_ptr->active_in_degree);
	      ++pim_ptr;
	    }
	  if(v_ptr->active_in_degree == 0)
	    ParentFreeVertices.push(v_ptr->index);
	  if(v_ptr->inout_degree >0)
	    {
	      struct messagestruct *pcim_ptr=v_ptr->pcim_ptr;
	      //bi-directional edges
	      for(int l=0; l<v_ptr->inout_degree && noconflict; ++l)
		{
		  if(pcim_ptr->v_ptr->occupied)
		    noconflict=false;
		  ++pcim_ptr;
		}
	    }
	}
      ++v_ptr;
    }
  if(noconflict==false)
    {
      cerr<<
	"Problem in refinement: Directed cycle due to bi-directional edges.\n";
      return false;
    }
  while( !ParentFreeVertices.empty() )
    {
      v_ptr = &Vertex[ ParentFreeVertices.front() ];
      ParentFreeVertices.pop();
      --ResidueSize;
      struct messagestruct *cim_ptr = v_ptr->cim_ptr;
      for(int k=0; k<v_ptr->out_degree; ++k)
	{
	  if(cim_ptr->v_ptr->occupied)
	    {
	      if( (--(cim_ptr->v_ptr->active_in_degree)) == 0)
		ParentFreeVertices.push(cim_ptr->v_ptr->index);
	    }
	  ++cim_ptr;
	}
    }
  if(ResidueSize >0)
    {
      cerr<<"Problem in refinment. The size of the residue cyclic subgraph is "
	  <<ResidueSize<<endl;
      return false;
    }
  ofstream optfile(dfvsfname.c_str() );
  optfile<<numEmptyNew
	 <<endl<<endl;
  v_ptr=&Vertex[1];
  for(int v=1; v<=VertexNumber; ++v)
    {
      if(v_ptr->occupied == false)
	optfile << v_ptr->index<<endl;
      ++v_ptr;
    }
  optfile.close();
  return true;
}
