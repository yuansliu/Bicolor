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

void split_kmer(char* str, vector<int>& kmer_vec) {
	const char *sep = ",";
	char *p = strtok(str, sep);
	while(p) {
		kmer_vec.push_back(atoi(p));
		p = strtok(NULL, sep);
	}
}

#define MAX_LONG_READ_LEN 500000
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

struct NODE {
	int cnt;
	char ch;
	NODE():cnt(0),ch('-') {}
};

bool cmpCnt(const NODE &a, const NODE &b) {
	return a.cnt > b.cnt;
}

int main(int argc, char const *argv[]) {
	// string long_reads_file = "DATA/_pacbio-test.fa", short_reads_files = "DATA/ill-test-5K-1.fa,DATA/ill-test-5K-2.fa", output_file = "test.fasta";
	// string long_reads_file = "long.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	// string long_reads_file = "long_only_one.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	// string long_reads_file = "long_three.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	string long_reads_file = argv[1], output_file = argv[3];
	// string long_reads_file = "long_one_no_solid.fasta", short_reads_files = "ERR022075_smp.fasta", output_file = "test.fasta";
	int nCores = DEF_THREADS;
	
    char *folder = new char[100];
	
	sprintf(folder, argv[2]);

	cerr << folder << endl;

	IBank *ptrIbk = NULL;

	try {
		ptrIbk = Bank::open(long_reads_file);
	} catch (gatb::core::system::Exception& bankPBExc) {
		cout << "Error message PacBio bank:\n" << bankPBExc.getMessage() << "\n";
		return 0;
    }
    Iterator<Sequence> *longSeqIt = ptrIbk->iterator();
    BankFasta output(output_file);
    // cerr<<"xxx\n";

	cerr << "begin to combine long reads...\n";

    ISynchronizer *sync = System::thread().newSynchronizer();

    IDispatcher::Status status = Dispatcher(nCores).iterate(longSeqIt, [&](const Sequence &seq) {
    	char *read = new char[MAX_LONG_READ_LEN];

    	char *align_file = new char[100];
    	sprintf(align_file, "%s/_%d.txt", folder, seq.getIndex());

    	BankFasta bfile(align_file);
    	BankFasta::Iterator bfile_it (bfile);

    	bool flag = false;
    	vector<string> aligned;
    	int aa, bb;

    	for (bfile_it.first(); !bfile_it.isDone(); bfile_it.next()) {
    		// cerr << bfile_it->getComment() << endl;
    		sscanf(bfile_it->getComment().c_str(),  "%d-%d", &aa, &bb);
    		if (bb == 0) {
    			if (!flag) {
    				flag = true;
    				aligned.push_back(bfile_it->toString());
    			}
    		} else {
	    		aligned.push_back(bfile_it->toString());
	    		// cerr << bfile_it->toString() << endl;
    		}
    	}
    	
    	int aligned_size = aligned.size();

    	if (aligned_size == 1) {
    		strcpy(read, aligned[0].c_str());
    	} else {
    		int final_len = aligned[0].length();

    		int read_len = 0;
    		//version 1
			/*for (int i = 0; i < final_len; i++) {
				int ca = 0, ct = 0, cg = 0, cc = 0, cn = 0;
				for (int j = 0; j < aligned_size; j++) {
					if (aligned[j][i]=='A') ca++; else
					if (aligned[j][i]=='T') ct++; else
					if (aligned[j][i]=='G') cg++; else
					if (aligned[j][i]=='C') cc++; else
					if (aligned[j][i]=='-') cn++;
				}
				if (ca == ct && ct == cg && cg == cc ) {
					continue;
				}
				int max_atgc = MAX(MAX(MAX(ca, ct), cg), cc);
				int max_count = MAX(max_atgc, cn);
				if( max_count == ca) read[read_len++] = 'A'; else
				if( max_count == ct) read[read_len++] = 'T'; else
				if( max_count == cg) read[read_len++] = 'G'; else
				if( max_count == cc) read[read_len++] = 'C'; 
				// else {
				// 	if( max_atgc == ca) read[read_len++] = 'A'; else
	   //  			if( max_atgc == ct) read[read_len++] = 'T'; else
	   //  			if( max_atgc == cg) read[read_len++] = 'G'; else
	   //  			if( max_atgc == cc) read[read_len++] = 'C';
				// }
			}*/

			//version 2
    		for (int i = 0; i < final_len; i++) {
		    	NODE *node_fre = new NODE[5];
		    	node_fre[0].ch = 'A'; node_fre[1].ch = 'T'; node_fre[2].ch = 'G'; node_fre[3].ch = 'C';

    			for (int j = 0; j < aligned_size; j++) {
    				if (aligned[j][i]=='A') node_fre[0].cnt++; else
    				if (aligned[j][i]=='T') node_fre[1].cnt++; else
    				if (aligned[j][i]=='G') node_fre[2].cnt++; else
    				if (aligned[j][i]=='C') node_fre[3].cnt++; else
    				if (aligned[j][i]=='-') node_fre[4].cnt++;
    			}
    			sort(node_fre, node_fre+5, cmpCnt);

				if (node_fre[0].ch == '-') {
					if (node_fre[1].cnt > 0) {
	    				read[read_len++] = node_fre[1].ch;
					}
    			} else {
    				read[read_len++] = node_fre[0].ch;
    			}	    	
    			delete[] node_fre;
    		}

			//version 3
    		/*for (int i = 0; i < final_len; i++) {
		    	NODE *node_fre = new NODE[5];
		    	NODE *node_sol = new NODE[4];
		    	node_fre[0].ch = 'A'; node_fre[1].ch = 'T'; node_fre[2].ch = 'G'; node_fre[3].ch = 'C';
		    	node_sol[0].ch = 'A'; node_sol[1].ch = 'T'; node_sol[2].ch = 'G'; node_sol[3].ch = 'C';

    			for (int j = 0; j < aligned_size; j++) {
    				if (aligned[j][i]=='A') node_fre[0].cnt++; else
    				if (aligned[j][i]=='T') node_fre[1].cnt++; else
    				if (aligned[j][i]=='G') node_fre[2].cnt++; else
    				if (aligned[j][i]=='C') node_fre[3].cnt++; else
    				if (aligned[j][i]=='-') node_fre[4].cnt++;
    			}
    			sort(node_fre, node_fre+5, cmpCnt);

    			if ( read_len < kmer_vec[kmer_vec_len-1]) {
	    			if (node_fre[0].ch == '-') {
	    				read[read_len++] = node_fre[1].ch;
	    			} else {
	    				read[read_len++] = node_fre[0].ch;
	    			}
	    		} else {
	    			//count the number of solid
	    			char *temp_str = new char[kmer_vec[kmer_vec_len]+10];
	    			for (int k = 0; k < 4; k++) {
	    				read[read_len] = node_sol[k].ch;
	    				for (int idx = 0; idx < kmer_vec_len; idx++) {
    						copyStr(temp_str, &read[read_len - kmer_vec[idx] + 1], kmer_vec[idx]);
    						Node node = graph[idx].buildNode(temp_str);
							if (graph[idx].contains(node) && graph[idx].indegree(node) > 0 && graph[idx].outdegree(node) > 0) {
								node_sol[k].cnt ++;
							}
    					}
	    			}
	    			sort(node_sol, node_sol+4, cmpCnt);
	    			// if (node_sol[0].cnt >= kmer_vec_len/2) {
	    			if (node_sol[0].cnt > 0) {
	    				read[read_len++] = node_sol[0].ch;
	    			} else {
	    				if (node_fre[0].ch == '-') {
		    				read[read_len++] = node_fre[1].ch;
		    			} else 
		    			if (node_fre[0].cnt >= kmer_vec_len/2) {
		    				read[read_len++] = node_fre[0].ch;
		    			}	
	    			}
	    			delete[] temp_str;
	    		}

    			delete[] node_fre;
    			delete[] node_sol;
    		}*/

    		read[read_len] = '\0';
    	}
    	//find the final result form _%s(temp_file); save to read
    	// sprintf(cmd, "rm -f %s\n", align_file);
    	// system(cmd);
    	// delete[] temp_file;
    	delete[] align_file;
    	
    	Sequence final_seq(read);
    	final_seq._comment = seq.getComment();
    	{
    		LocalSynchronizer local(sync);
    		output.insert(final_seq);
    	}
    	delete[] read;
    });

    output.flush();
    delete sync;

	return 0;
}