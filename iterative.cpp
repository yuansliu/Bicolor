#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <stack>
#include <getopt.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include <gatb/gatb_core.hpp>

using namespace std;

#define MAX_READ_LEN 500000
#define MAX_PATH_LEN 1000
#define DEF_MAX_BRANCH 200
#define DEF_ERROR_RATE 0.40
#define DEF_TRIALS 5
#define DEF_THREADS 0
#define MIN(a,b) ((a) < (b) ? (a) : (b)) // min macro
#define MAX(a,b) ((a) > (b) ? (a) : (b)) // max macro

// alignment scores
#define ALIGN_MATCH 1
#define ALIGN_MISMATCH -3
#define ALIGN_INDEL -2

float max_error_rate  = DEF_ERROR_RATE;
int max_branch = DEF_MAX_BRANCH;
int max_trials  = DEF_TRIALS;
int strict_mode = 0;

void reverse(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    switch(source[len-i-1]) {
    case 'a':
      target[i] = 't';
      break;
    case 'c':
      target[i] = 'g';
      break;
    case 'g':
      target[i] = 'c';
      break;
    case 't':
      target[i] = 'a';
      break;
    case 'n':
      target[i] = 'n';
      break;
    case 'A':
      target[i] = 'T';
      break;
    case 'C':
      target[i] = 'G';
      break;
    case 'G':
      target[i] = 'C';
      break;
    case 'T':
      target[i] = 'A';
      break;
    case 'N':
      target[i] = 'N';
      break;
    }
  }
  target[len] = '\0';
}

class stack_element {
public:
	Node  node;
	int pos;
	char  nt;

	stack_element(Node _node, int _pos, char _nt):node(_node), pos(_pos), nt(_nt) {}
};


int extend( Graph graph, int kmer_len, Node begin, char *extended_path, int *max_b, char *read_part, int part_len, int *dp, int *best_ed, int *best_part_len) {
	std::stack <stack_element *>  nodes;
	char path[MAX_PATH_LEN + 1];
	int best_len = -1;
	*best_ed  = max_error_rate * (part_len + kmer_len) + 1;
	*best_part_len  = 0;

	/* store the neighbors of begin node in the stack */
	// Graph::Vector<Edge> neighbors = graph.successors<Edge>( begin );
	GraphVector<Edge> neighbors = graph.successorsEdge (begin);

	for ( int i = neighbors.size() - 1; i >= 0; i-- ) {
		nodes.push( new stack_element( neighbors[i].to, 1, ascii( neighbors[i].nt ) ) );
	}

	while ( nodes.size() > 0 && *max_b > 0 ) {
		stack_element *current = nodes.top();
		nodes.pop();
		Node  cnode = current->node;
		int pos = current->pos;
		path[pos - 1] = current->nt;
		delete current;

		/* compute the current row of the DP matrix and save its minimum */
		dp[pos * (MAX_PATH_LEN + 1) + 0] = pos;
		/* dp[pos][0] = pos; */
		int min = MAX_PATH_LEN;
		int mink  = -1;
		for ( int k = 1; k <= part_len; k++ ) {
			if ( path[pos - 1] == read_part[k - 1] )
			{
				dp[pos * (MAX_PATH_LEN + 1) + k] = MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k - 1], MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k] + 1, dp[pos * (MAX_PATH_LEN + 1) + k - 1] + 1 ) );
			} else {
				dp[pos * (MAX_PATH_LEN + 1) + k] = MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k - 1] + 1, MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k] + 1, dp[pos * (MAX_PATH_LEN + 1) + k - 1] + 1 ) );
			}
			if ( dp[pos * (MAX_PATH_LEN + 1) + k] < min ) {
				min = dp[pos * (MAX_PATH_LEN + 1) + k];
				mink  = k;
			}
		}

		/* Check if this path can eventually be under the error threshold */
		if ( min <= max_error_rate * (part_len + kmer_len) + 1 ) {
			/* Is it currently under the error threshold? */
			if ( (mink + kmer_len) * max_error_rate >= min ) {
				/* Is it better than what we have found before? */
				if ( mink > *best_part_len || (mink == *best_part_len && min < *best_ed) ) {
					best_len  = pos;
					*best_part_len  = mink;
					*best_ed  = min;
					memcpy( extended_path, path, pos );
				}
			}
			// Graph::Vector<Edge> neighbors = graph.successors<Edge>( cnode );
			GraphVector<Edge> neighbors = graph.successorsEdge (cnode);
			/*      if (neighbors.size() == 1) { */
			for ( int i = neighbors.size() - 1; i >= 0; i-- ) {
				Edge e = neighbors[i];
				nodes.push( new stack_element( e.to, pos + 1, ascii( e.nt ) ) );
			}
			if ( neighbors.size() == 0 ) {
				(*max_b)--;
			}
		} else {
		  (*max_b)--;
		}
	}

	/* Empty the stack and free memory */
	while ( nodes.size() > 0 ) {
		stack_element *current = nodes.top();
		nodes.pop();

		delete current;
	}

	/* Check for strict mode */
	if ( strict_mode && *max_b <= 0 ) {
		return(-1);
	}

	/* Compute the best correcting extension */
	int align_dp1[MAX_PATH_LEN];
	int align_dp2[MAX_PATH_LEN];

	int max = 0, maxi = 0, maxj = 0;

	/* Initialize the first row */
	for ( int i = 0; i <= best_len; i++ ) {
		align_dp1[i] = i * ALIGN_INDEL;
	}

	for ( int j = 1; j <= *best_part_len; j++ ) {
		align_dp2[0] = j * ALIGN_INDEL;
		for ( int i = 1; i <= best_len; i++ ) {
			if ( read_part[j - 1] == extended_path[i - 1] ) {/* match */
				align_dp2[i] = MAX( align_dp1[i - 1] + ALIGN_MATCH, MAX( align_dp1[i] + ALIGN_INDEL, align_dp2[i - 1] + ALIGN_INDEL ) );
			} else {                                        /* mismatch */
				align_dp2[i] = MAX( align_dp1[i - 1] + ALIGN_MISMATCH, MAX( align_dp1[i] + ALIGN_INDEL, align_dp2[i - 1] + ALIGN_INDEL ) );
			}
			if ( align_dp2[i] >= max ) {
				max = align_dp2[i];
				maxi  = i;
				maxj  = j;
			}
		}

		if ( j + 1 <= *best_part_len ) {
			j++;
			align_dp1[0] = j * ALIGN_INDEL;
			for ( int i = 1; i <= best_len; i++ ) {
				if ( read_part[j - 1] == extended_path[i - 1] ) {/* match */
					align_dp1[i] = MAX( align_dp2[i - 1] + ALIGN_MATCH, MAX( align_dp2[i] + ALIGN_INDEL, align_dp1[i - 1] + ALIGN_INDEL ) );
				} else {/* mismatch */
					align_dp1[i] = MAX( align_dp2[i - 1] + ALIGN_MISMATCH, MAX( align_dp2[i] + ALIGN_INDEL, align_dp1[i - 1] + ALIGN_INDEL ) );
				}
				if ( align_dp1[i] >= max ) {
					max = align_dp1[i];
					maxi = i;
					maxj = j;
				}
			}
		}
	}

	best_len  = maxi;
	*best_part_len  = maxj;

	extended_path[best_len] = '\0';
	return best_len;
}

int best_path( Graph graph, Node begin, Node end, char *best_path, int *max_b, char *read_part, int part_len, int *dp, int *best_ed ) {
  std::stack <stack_element *>  nodes;
  char        path[MAX_PATH_LEN + 1];
  *best_ed = max_error_rate * part_len + 1;
  int best_len = -1;

  /* store the neighbors of begin node in the stack */
  GraphVector<Edge> neighbors = graph.successorsEdge(begin);

  for ( int i = neighbors.size() - 1; i >= 0; i-- )
  {
    nodes.push( new stack_element( neighbors[i].to, 1, ascii( neighbors[i].nt ) ) );
  }

  while ( nodes.size() > 0 && *max_b > 0 )
  {
    stack_element *current = nodes.top();
    nodes.pop();
    Node  cnode = current->node;
    int pos = current->pos;
    path[pos - 1] = current->nt;

    delete current;

    if ( pos >= MAX_PATH_LEN )
    {
      std::cout << "Too long path searched for" << std::endl;
      exit( EXIT_FAILURE );
    }

    dp[pos * (MAX_PATH_LEN + 1) + 0] = pos;
    int min = MAX_PATH_LEN;
    for ( int k = 1; k <= part_len; k++ )
    {
      if ( path[pos - 1] == read_part[k - 1] )
      {
        dp[pos * (MAX_PATH_LEN + 1) + k] = MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k - 1], MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k] + 1, dp[pos * (MAX_PATH_LEN + 1) + k - 1] + 1 ) );
      } else {
        dp[pos * (MAX_PATH_LEN + 1) + k] = MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k - 1] + 1, MIN( dp[(pos - 1) * (MAX_PATH_LEN + 1) + k] + 1, dp[pos * (MAX_PATH_LEN + 1) + k - 1] + 1 ) );
      }
      if ( dp[pos * (MAX_PATH_LEN + 1) + k] < min )
        min = dp[pos * (MAX_PATH_LEN + 1) + k];
    }

    /*
     * Can the current path eventually be extended to be better that
     * the best we have found so far?
     */
    if ( min < *best_ed )
    {
      /* Have we reached the end node? */
      if ( cnode == end )
      {
        (*max_b)--;
        /* int ed = dp[pos][part_len]; */
        int ed = dp[pos * (MAX_PATH_LEN + 1) + part_len];
        /* Is this path better than the previous one we have found? */
        if ( ed < *best_ed )
        {
          *best_ed = ed;
          memcpy( best_path, path, pos );
          best_len = pos;
        }
      } else {
        // Graph::Vector<Edge> neighbors = graph.successors<Edge>( cnode );
        GraphVector<Edge> neighbors = graph.successorsEdge (cnode);
        for ( int i = neighbors.size() - 1; i >= 0; i-- )
        {
          Edge e = neighbors[i];
          nodes.push( new stack_element( e.to, pos + 1, ascii( e.nt ) ) );
        }
        if ( neighbors.size() == 0 )
        {
          (*max_b)--;
        }
      }
    } else {
      (*max_b)--;
    }
  }

  /* Empty the stack and free memory */
  while ( nodes.size() > 0 )
  {
    stack_element *current = nodes.top();
    nodes.pop();

    delete current;
  }

  if ( strict_mode && *max_b <= 0 )
  {
    return(-1);
  } else {
    return(best_len);
  }
}

void copy_lower_case(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    target[i] = tolower(source[i]);
  }
  // target[len] = '\0';
}

void copy_upper_case(char *target, char *source, int len) {
  for(int i = 0; i < len; i++) {
    target[i] = toupper(source[i]);
  }
  // target[len] = '\0';
}

int correct_one_read(Sequence seq, char *corrected, Graph graph, int kmer_len) {
	int corrected_len = 0;
    char *read = seq.getDataBuffer();
    int read_len = seq.getDataSize();

    read[read_len] = '\0';
    
    // Kmer<>::ModelDirect model(kmer_len);
    // Kmer<>::ModelDirect::Iterator itKmer(model);
    // itKmer.setData (seq.getData());

    // We iterate the kmers.
    // First enumerate all positions with solid k-mers
    int pos = 0;
    int num_solid = 0;
    int *solid = new int[MAX_READ_LEN];

	int idx = 0;

	char *kmers0 = new char [kmer_len + 1], *kmers1 = new char [kmer_len + 1];
	kmers0[kmer_len] = '\0'; kmers1[kmer_len] = '\0';

	for (idx = 0; idx < read_len - kmer_len + 1; idx++) {
		for (int i = 0, j = idx; i < kmer_len; i++, j++) {
			kmers0[i] = read[j];
		}
		for (int i = 0, j = idx+1; i < kmer_len; i++, j++) {
			kmers1[i] = read[j];
		}

		Node node0 = graph.buildNode(kmers0);
		Node node1 = graph.buildNode(kmers1);
		if (graph.contains(node0) && graph.contains(node1) && graph.indegree(node0) > 0 && graph.outdegree(node0) > 0 && graph.indegree(node1) > 0 && graph.outdegree(node1) > 0) {
			solid[num_solid++] = idx;
		}
	}

 //    for (itKmer.first(); !itKmer.isDone(); itKmer.next())   {  
	//   Node node = graph.buildNode(Data((char *)(model.toString(itKmer->value()).c_str())));
 //      // The graph has false positives -> filtering out nodes without incoming or outgoing edges helps
 //      if (graph.contains(node) && graph.indegree(node) >= 1 && graph.outdegree(node)>=1) {
	// solid[num_solid++] = pos;
 //      }
 //      pos++;
 //    }

    // Do we need this? (May be something strange with seq.getDataSize()?)
    read_len=strlen(read);

    if (num_solid == 0) {
      // No solid k-mers -> just copy the read as it is
      copy_lower_case(&corrected[corrected_len], &read[0], read_len);
      corrected_len += read_len;
      corrected[corrected_len] = '\0';
      delete [] solid;
      return 0;
    }

    int *dp;
    char *path = new char[MAX_PATH_LEN+1];
    char *path2 = new char[MAX_PATH_LEN+1];
    char *buffer = new char[MAX_READ_LEN+1];
    
    // int dp[MAX_PATH_LEN+1][MAX_PATH_LEN+1];

    // Allocate memory for dp table
    dp = new int[(MAX_PATH_LEN+1)*(MAX_PATH_LEN+1)];
    // Initialize DP matrix
    for(int k=0; k <= (1+max_error_rate)*read_len+1; k++) {
      // dp[0][k] = k;
      dp[0*(MAX_PATH_LEN+1)+k] = k;
    }

    if (solid[0] > 0) {
      // There is a head to correct
      if ((1.0+max_error_rate)*solid[0] < MAX_PATH_LEN) {

	// Attempt to correct the head of the read
	Node n = graph.buildNode(Data(&read[solid[0]]));
	n = graph.reverse(n);
	reverse(buffer, &read[0], solid[0]);

	int ed;
	int best_part_len;

	int max_b = max_branch;
	int len = extend(graph, kmer_len, n, path, &max_b, buffer, solid[0], dp, &ed, &best_part_len);

	if (len >= 0) {
	  path[len] = '\0';
	  copy_lower_case(&corrected[corrected_len], &read[0], solid[0]-best_part_len);
	  corrected_len += solid[0]-best_part_len;

	  reverse(buffer, path, len);
	  copy_upper_case(&corrected[corrected_len], buffer, len);
	  copy_upper_case(&corrected[corrected_len+len], &read[solid[0]], kmer_len);
	  corrected_len += len+kmer_len;

	} else {
	  // No correction of the head was found
	  copy_lower_case(&corrected[corrected_len], &read[0], solid[0]);
	  copy_upper_case(&corrected[corrected_len+solid[0]], &read[solid[0]], kmer_len);
	  corrected_len += (solid[0]+kmer_len);

	}
      } else {
	// The head is too long: no correction was attempted (memory constraints)
	copy_lower_case(&corrected[corrected_len], &read[0], solid[0]);
	copy_upper_case(&corrected[corrected_len+solid[0]], &read[solid[0]], kmer_len);
	corrected_len += (solid[0]+kmer_len);
      }
    } else {
      // The first k-mer is solid
      copy_upper_case(&corrected[corrected_len], &read[solid[0]], kmer_len);
      corrected_len += kmer_len;
    }

    ////////////////////////////////////////////////////////
    // INNER REGION CORRECTION

    if (num_solid >= 2) {
      // There is an intermediate part to correct

      // Data structure for a correction for the region between two solid k-mers
      struct Path {
	int ed;      // Edit distance between the region and the corrected region
	int len;     // Length of the corrected region
	char *str;   // Corrected region
      };
      
      // We use the graph from boost library
      // Path graph: describes the corrections that we have found between the solid k-mers
      //   - Nodes: solid k-mers: identified by their position in the solid-array
      //   - Edges: Found paths between the solid k-mers in dbg

      // Edge in the path graph
      typedef std::pair<int, int> Edge;
      // Path graph: no properties attached to the nodes, Path struct (above) gives the properties attached to the edges
      typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS, boost::no_property, Path > path_graph_t;
      typedef boost::graph_traits < path_graph_t >::vertex_descriptor vertex_descriptor;
      typedef boost::graph_traits < path_graph_t >::edge_descriptor edge_descriptor;
      
      // Reserve memory for edges and edge properties
      Edge *paths;
      struct Path *path_prop;
      paths = (Edge *)MemoryAllocatorStdlib::singleton().malloc(MAX_READ_LEN*max_trials*sizeof(Edge));
      path_prop = (struct Path *)MemoryAllocatorStdlib::singleton().malloc(MAX_READ_LEN*max_trials*sizeof(struct Path));

      int num_paths=0;

      for (int i = 0; i < num_solid-1;i++) {
	int found = 0;
	
	// Number of correction paths we have attempted to find for this k-mer
	int trials = 0;
	for(int j = i+1; j < num_solid && trials < max_trials; j++) {
	  if (solid[j] == solid[i] + (j-i) ) { // run of solid k-mers (all adjacent)
	    trials++;
	    if (num_paths >= MAX_READ_LEN*max_trials) {
	      std::cout << "Not enough memory allocated" << std::endl;
	      exit(EXIT_FAILURE);
	    }
	    paths[num_paths] = std::make_pair(i,j);
	    path_prop[num_paths].ed = 0;
	    path_prop[num_paths].len = solid[j]-solid[i];
	    path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((solid[j]-solid[i]+1)*sizeof(char));
	    copy_upper_case(path_prop[num_paths].str, &read[solid[i]+kmer_len], solid[j]-solid[i]);
	    num_paths++;
	    found++;
	  } else if (solid[j] - solid[i] < (1.0-max_error_rate)*kmer_len) { // solid k-mers are within the same k-mer: too near from eachother
	    continue;
	  } else if ((1.0+max_error_rate)*(solid[j] - solid[i]) >= MAX_PATH_LEN) { // solid k-mers are too far from eachother: not enough memory for the DP matrix
	    break;
	  } else { // solid k-mers are neither too near, nor to far: try find a path
	    trials++;

	    Node ni = graph.buildNode(Data(&read[solid[i]]));
	    Node nj = graph.buildNode(Data(&read[solid[j]]));

	    strncpy(buffer, &read[solid[i]+kmer_len], solid[j]-solid[i]);
	    buffer[solid[j]-solid[i]] = '\0';

	    int max_b = max_branch;
	    int ed = MAX_PATH_LEN;                       // computed edit distance value


	    int len = best_path(graph, ni, nj, path, &max_b, buffer, solid[j]-solid[i], dp, &ed);
	    if (len >= 0) {
	      path[len] = '\0';
	      if (num_paths >= MAX_READ_LEN*max_trials) {
		std::cout << "Not enough memory allocated" << std::endl;
		exit(EXIT_FAILURE);
	      }
	      paths[num_paths] = std::make_pair(i,j); // insert a new edge in the path graph
	      // record the prop of the found path
	      path_prop[num_paths].ed= ed;
	      path_prop[num_paths].len = len;
	      path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((len+1)*sizeof(char));
	      copy_upper_case(path_prop[num_paths].str, path, len); // copy the path to be saved with the path
	      num_paths++;
	      found++;

	    } else { // no found path
	    }
	  }
	} // End of trials for k-mer at position i
	if (found == 0) {
	  if (solid[i+1]-solid[i] > kmer_len) {
	    // Attempt to extend from the left side up to half the gap
	    Node n = graph.buildNode(Data(&read[solid[i]]));
	    int l = (solid[i+1]-solid[i]-kmer_len)/2;

	    if ((1.0+max_error_rate)*(l+kmer_len) >= MAX_PATH_LEN-1) {
	      l = (int)(MAX_PATH_LEN/(1.0+max_error_rate))-2-kmer_len;
	    }

	    strncpy(buffer, &read[solid[i]+kmer_len], l);
	    buffer[l] = '\0';

	    int ed = 0;
	    int best_part_len = 0;
	    int max_b = max_branch;
	    int len = extend(graph, kmer_len, n, path, &max_b, buffer, l, dp, &ed, &best_part_len);

	    // Attempt to extend the gap from the right up to half the gap
	    Node n2 = graph.reverse(graph.buildNode(Data(&read[solid[i+1]])));

	    // If the gap is of odd length add the extra base to this extension
	    if ((solid[i+1]-solid[i]-kmer_len) % 2 == 1)
	      l++;

	    reverse(buffer, &read[solid[i+1]-l], l);
	    buffer[l] = '\0';

	    int ed2 = 0;
	    int best_part_len2 = 0;
	    int max_b2 = max_branch;
	    int len2 = extend(graph, kmer_len, n2, path2, &max_b2, buffer, l, dp, &ed2, &best_part_len2);

	    if (num_paths >= MAX_READ_LEN*max_trials) {
	      std::cout << "Not enough memory allocated" << std::endl;
	      exit(EXIT_FAILURE);
	    }

	    if (len < 0) {
	      ed = 0;
	      len = 0;
	      best_part_len = 0;
	    }
	    if (len2 < 0) {
	      ed2 = 0;
	      len2 = 0;
	      best_part_len2 = 0;
	    }

	    int uncorrected_len = solid[i+1]-solid[i]-kmer_len - best_part_len - best_part_len2;
	    if (uncorrected_len < 0) {
	      std::cout << "Overlapping gap extensions!" << std::endl;
	      exit(EXIT_FAILURE);
	    }

	    paths[num_paths] = std::make_pair(i,i+1);
	    path_prop[num_paths].ed = ed+ed2+uncorrected_len;
	    path_prop[num_paths].len = len+uncorrected_len+len2+kmer_len;
	    path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((path_prop[num_paths].len+1)*sizeof(char));

	    if (len > 0) {
	      copy_upper_case(path_prop[num_paths].str, path, len); // Copy the left extension
	    }	      
	    if (uncorrected_len > 0) {
	      copy_lower_case(&path_prop[num_paths].str[len], &read[solid[i] + kmer_len + best_part_len], uncorrected_len); // Copy the uncorrected part
	    }
	    if (len2 > 0) {
	      reverse(buffer, path2, len2);
	      copy_upper_case(&path_prop[num_paths].str[len+uncorrected_len], buffer, len2); // Copy the right extension
	    }
	    copy_upper_case(&path_prop[num_paths].str[len+uncorrected_len+len2], &read[solid[i+1]], kmer_len); // Copy the right k-mer

	    num_paths++;

	  } else {
	    // Add a dummy edge if the kmers overlap so that at least
	    // one path from 1st to last solid k-mer is found
	    if (num_paths >= MAX_READ_LEN*max_trials) {
	      std::cout << "Not enough memory allocated" << std::endl;
	      exit(EXIT_FAILURE);
	    }
	    paths[num_paths] = std::make_pair(i,i+1);
	    path_prop[num_paths].ed = solid[i+1]-solid[i];
	    path_prop[num_paths].len = solid[i+1]-solid[i];
	    path_prop[num_paths].str = (char *)MemoryAllocatorStdlib::singleton().malloc((solid[i+1]-solid[i]+1)*sizeof(char));
	    if (solid[i+1]-solid[i]-kmer_len > 0) {
	      copy_lower_case(path_prop[num_paths].str, &read[solid[i]+kmer_len], solid[i+1]-solid[i]-kmer_len);
	      copy_upper_case(&path_prop[num_paths].str[solid[i+1]-solid[i]-kmer_len], &read[solid[i+1]], kmer_len);
	    } else {
	      copy_upper_case(path_prop[num_paths].str, &read[solid[i]+kmer_len], solid[i+1]-solid[i]);
	    }
	    num_paths++;
	  }
	}
      }

      path_graph_t path_graph(num_solid);

      for (int i = 0; i < num_paths; i++) {
	boost::add_edge(paths[i].first, paths[i].second, path_prop[i], path_graph);
      }

      std::vector<vertex_descriptor> p(num_vertices(path_graph));
      std::vector<int> d(num_vertices(path_graph));
      vertex_descriptor s = vertex(0, path_graph);

      // Find the shortest paths from the first solid k-mer
      boost::dijkstra_shortest_paths(path_graph, s,
				     boost::weight_map(get(&Path::ed, path_graph)).
				     predecessor_map(boost::make_iterator_property_map(p.begin(), get(boost::vertex_index, path_graph))).
				     distance_map(boost::make_iterator_property_map(d.begin(), get(boost::vertex_index, path_graph))));

      // Backtrack the shortest path from the last solid k-mer to the first
      int *successor = new int[MAX_READ_LEN];
      int curr = num_solid-1;

      while (curr != 0) {
	if (curr < 0 || curr >= num_solid) {
	  std::cout << "Curr: " << curr << std::endl;
	  exit(EXIT_FAILURE);
	}
	successor[p[curr]] = curr;
	curr = p[curr];
      }

      // Traverse the shortest path and correct the read
      while(curr != num_solid-1) {
	boost::graph_traits < path_graph_t>::edge_descriptor  e = boost::edge(curr, successor[curr], path_graph).first;
	if (path_graph[e].len < 0 || corrected_len + path_graph[e].len >= MAX_READ_LEN) {
	  std::cout << "Invalid path length in the path graph" << std::endl;
	  exit(EXIT_FAILURE);
	}
	strncpy(&corrected[corrected_len], path_graph[e].str, path_graph[e].len);
	corrected_len += path_graph[e].len;
	curr = successor[curr];
      }

      // Free the memory of the path graph
      for(int i = 0; i < num_paths; i++) {
	MemoryAllocatorStdlib::singleton().free(path_prop[i].str);
      }
      
      MemoryAllocatorStdlib::singleton().free(paths);
      MemoryAllocatorStdlib::singleton().free(path_prop);
      delete [] successor;
    }

    // Attempt to correct the tail of the read
    if (solid[num_solid-1] < read_len-kmer_len) {
      if ((1.0+max_error_rate)*(read_len-solid[num_solid-1]-kmer_len) < MAX_PATH_LEN) {
	Node n = graph.buildNode(Data(&read[solid[num_solid-1]]));

	strncpy(buffer, &read[solid[num_solid-1]+kmer_len], read_len-solid[num_solid-1]-kmer_len);
	buffer[read_len-solid[num_solid-1]-kmer_len] = '\0';

	int ed;
	int best_part_len;

	int max_b = max_branch;
	int len = extend(graph, kmer_len, n, path, &max_b, buffer, read_len-solid[num_solid-1]-kmer_len, dp, &ed, &best_part_len);
	if (len >= 0) {
	  path[len] = '\0';
	  copy_upper_case(&corrected[corrected_len], path, len);
	  corrected_len += len;
	  if (read_len - solid[num_solid-1] - kmer_len > best_part_len) {
	    copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len+best_part_len],
			    read_len - solid[num_solid-1]-kmer_len-best_part_len);
	    corrected_len += read_len - solid[num_solid-1]-kmer_len-best_part_len;
	  }

	} else {
	  copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len], read_len - solid[num_solid-1]-kmer_len);
	  corrected_len += (read_len - solid[num_solid-1]-kmer_len);

	}
      } else {
	copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len], read_len - solid[num_solid-1]-kmer_len);
	corrected_len += (read_len - solid[num_solid-1]-kmer_len);
      }
    // } else {
    //   copy_lower_case(&corrected[corrected_len], &read[solid[num_solid-1]+kmer_len], read_len - solid[num_solid-1]-kmer_len);
    //   corrected_len += (read_len - solid[num_solid-1]-kmer_len);
    }

    delete [] solid;
    delete [] dp;
    delete [] path;
    delete [] path2;
    delete [] buffer;

// cerr << "corrected_len: " << corrected_len << endl;

    corrected[corrected_len] = '\0';

// cerr << "-------\n";
// cerr << corrected << endl;
// cerr << "-------\n";

    return 1;
}

struct NODE {
	int cnt;
	char ch;
	NODE():cnt(0),ch('-') {}
};

void split_kmer(char* str, vector<int>& kmer_vec) {
	const char *sep = ",";
	char *p = strtok(str, sep);
	while(p) {
		kmer_vec.push_back(atoi(p));
		p = strtok(NULL, sep);
	}
}

bool cmpCnt(const NODE &a, const NODE &b) {
	return a.cnt > b.cnt;
}

int main(int argc, char const *argv[]) {
	// long.fasta short.fasta kmer_size foler
	// string long_reads_file = "DATA/_pacbio-test.fa", short_reads_files = "DATA/ill-test-5K-1.fa,DATA/ill-test-5K-2.fa", output_file = "test.fasta";
	// string long_reads_file = "long.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	// string long_reads_file = "long_only_one.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	// string long_reads_file = "long_three.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	string long_reads_file = argv[1], short_reads_files = argv[2];
	// string long_reads_file = "long_one_no_solid.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	int nCores = DEF_THREADS;

	// vector<int> kmer_vec;
	int kmer_vec_len = 0;
	// char kmerstr[] = "13,15,17,19,21,23";
	// char *kmerstr = new char[100];
	// sprintf(kmerstr, "%s", argv[3]);
	// split_kmer(kmerstr, kmer_vec);
	// kmer_vec_len = kmer_vec.size();
	int kmer_size;
	kmer_size = atoi(argv[3]);
	kmer_vec_len = atoi(argv[4]);

	char *folder = new char[100];
	sprintf(folder, argv[5]);
	// cerr << folder << endl;

	cerr << "noisy long reads: " << argv[1] << "; " << "short reads: " << argv[2]  << "; " << "initial length of kmer: " << argv[3] << "; " << "iterative rounds: " << argv[4] << "; " << "temporary file: " << argv[5] << endl;

	Graph *graph = new Graph[kmer_vec_len];

	IBank *ibk = Bank::open(short_reads_files);
	LOCAL(ibk);
	
	// cerr << "kmer_size: " << kmer_size << " ...\n";

	try {
		for (int i = 0; i < kmer_vec_len; i ++) {
			// kmer_size = kmer_vec[i];
			graph[i] = Graph::create(ibk, (const char *) "-kmer-size %d -abundance-min %d -bloom cache -debloom original -debloom-impl basic -nb-cores %d -verbose 0", kmer_size+2*i, 3, nCores);
			
			// cerr << "completed graph[" << kmer_size+2*i << "]\n";
			// graph[1] = Graph::create(ibk, (const char *) "-kmer-size %d -abundance-min %d -bloom cache -debloom original -debloom-impl basic -nb-cores %d -verbose 0", kmer_size+2, 3, nCores);
		}
	} catch (Exception& e){
		cerr << "EXCEPTION: " << e.getMessage() << "\n";
		return 0;
	}

	IBank *ptrIbk = NULL;

	try {
		ptrIbk = Bank::open(long_reads_file);
	} catch (gatb::core::system::Exception& bankPBExc) {
		cout << "Error message PacBio bank:\n" << bankPBExc.getMessage() << "\n";
		return 0;
    }
    Iterator<Sequence> *longSeqIt = ptrIbk->iterator();

	cerr << "begin correcting long reads...\n";

    //versioin 4
    IDispatcher::Status status = Dispatcher(nCores).iterate(longSeqIt, [&](const Sequence &seq) {
    	char *read = new char[MAX_READ_LEN];
    	char *rev_read = new char[MAX_READ_LEN];
    	char *result_str = new char[MAX_READ_LEN];
    	char *result_str1 = new char[MAX_READ_LEN];
    	int result = 0, temp_str_len;
    	// cerr<<seq.toString()<<"\n";
    	// toUpperCaseStr(result_str, seq.getDataBuffer());
    	copy_upper_case(read, seq.getDataBuffer(), seq.getDataSize());
    	read[seq.getDataSize()] = '\0';
    	// reverse(rev_read, read, strlen(read));

    	char *temp_file = new char[100];
    	sprintf(temp_file, "%s/%d.txt", folder, seq.getIndex());
    	// cerr << seq.getIndex() <<"\n";
    	
    	FILE *out = fopen(temp_file, "a");//append
    	//forward

    	strcpy(result_str1, read);

    	for (int i = 0; i < kmer_vec_len; i ++) {
    		// reverse(result_str, result_str1, strlen(result_str1));
	    	Sequence seq0(result_str1);
	    	seq0._comment = seq.getComment();
	    	result = correct_one_read(seq0, result_str, graph[i], kmer_size+2*i);

			temp_str_len = strlen(result_str);
	    	copy_upper_case(result_str1, result_str, temp_str_len);
	    	result_str1[temp_str_len] = '\0';    		

    	}
    	fprintf(out, ">%s-%d\n%s\n", argv[3], result, result_str1);
    	
    	fclose(out);
    	delete[] result_str;
    	delete[] result_str1;
    	delete[] rev_read;

    	delete[] read;
    });
    
    delete[] folder;

	return 0;
}