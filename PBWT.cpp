#if defined MEMORY_MAP
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
using namespace boost::iostreams;
#endif
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <limits>
#include <map>
#include <cassert>
#include <cstring>
#include <algorithm>

using namespace std;

int M, N, L, W, G;  // # of sequences, # of sites, length, width, gap size

struct SparseTable {
	int N;
	vector<vector<int>> table;

	SparseTable(int _N) {
		this->N = _N;
		table = vector<vector<int>>(N, vector<int>((int)log2(N) + 1));
	}
	
	void init() {
		for (int j = 1; (1 << j) <= N; ++j) {
			for (int i = 0; (i + (1 << j) - 1) < N; ++i) {
				table[i][j] = max(table[i][j - 1], table[i + (1 << (j - 1))][j - 1]);
			}
		}
	}

	int query(int l, int r) { //inclusive range
		int j = (int)log2(r - l + 1);
		return max(table[l][j], table[r - (1 << j) + 1][j]);
	}
};

struct State {
	int f_mini, f_maxi, r_mini, r_maxi;
	vector<int> IDs, zero, one;

	State() {
		f_mini = r_mini = numeric_limits<int>::max();
		f_maxi = r_maxi = -1;
		zero = vector<int>(G);
		one = vector<int>(G);
	}
};

int p1 = 0;
vector<vector<int>> gap; // stores haplotype data in the gap

int getGap(int g, int idx) {
	int jump = G - g;
	return gap[(p1 + G - jump) % G][idx];
}

int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0); cin.tie(0);

	string writeTo = string(argv[2]);
	#if defined MEMORY_MAP
	mapped_file_source file(writeTo + ".rpbwt");
	stream<mapped_file_source> backward(file, ios::binary);
	#else
	ifstream backward(writeTo + ".rpbwt");
	#endif
	ifstream in(argv[1]), sites(writeTo + ".sites"), meta(writeTo + ".meta");
	ofstream results(writeTo + ".results"), resultIDs(writeTo + ".IDs"), resultMI(writeTo + ".MI");

	int checkpoint = atoi(argv[3]);
	L = atoi(argv[4]), G = atoi(argv[5]);

	// retrieve M and N from meta file
	meta >> M >> N;

	// input chromosome site positions
	vector<int> positions(N);
	for (int i = 0; i < N; ++i) sites >> positions[i];

	// skip meta-info lines and get header line
	string header;
	while (getline(in, header)) {
		if ((int)header.size() < 2 || header[0] != '#' || header[1] != '#') break;
	}

	// input sample IDs
	vector<string> ID(M);
	stringstream ss(header);
	for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
	for (int i = 0; i < M; ++i) {
		getline(ss, ID[i], '\t');
	}

	vector<int> pre(M), div(M); // prefix and divergence array
	iota(pre.begin(), pre.end(), 0);
	vector<int> a(M), b(M), d(M), e(M);
	char s[2 * M + 5000]; // assumes fixed fields take up less than 5000 characters
	SparseTable forwardSparse(M), backwardSparse(M);
	vector<int> idx(M); // idx[i] = index of sample i in the reverse positional prefix array
	vector<int> backwardPre(M), block(M), blockSize(M + 1); // block[i] = block ID of sample i in the reverse PBWT
	gap = vector<vector<int>>(G, vector<int>(M));
	int fp = -1, bp = 0; // pointers for gemomic distance

	// initialize the gap
	for (int site = 0; site < G; ++site) {
		in.getline(s, 2 * M + 5000);
		// skip fixed fields
		int offset = 0, tabs = 0; // offset = position in "s" of first sequence - points to the first character after 9 tabs
		while (tabs < 9) {
			if (s[offset] == '\t') ++tabs;
			++offset;
		}
		for (int i = 0; i < M; ++i) {
			assert(s[offset + (i / 2) * 4 + 1] == '|'); // sanity check
			gap[p1][i] = (s[offset + 2 * i] == '0' ? 0 : 1);
		}
		p1 = (p1 + 1) % G;
	}

	for (int site = 0; site + G < N; ++site) {
		if (site != 0) {
			int rsite = (N - 1) - site - G; // index of the corresponding reverse site
			backward.seekg((long long)rsite * M * 8);

			// update genomic pointers
			while (positions[site - 1] - (fp != -1 ? positions[fp] : 0) >= L) ++fp; // points to first position < L
			while ((bp < N ? positions[bp] : numeric_limits<int>::max()) - positions[site + G] < L) ++bp; // points to the first position >= L
			int rsite_bp = (N - 1) - bp;

			// initialize backward sparse table, idx, backwardPre, and block
			int start = -1, id = 0;
			for (int i = 0; i < M; ++i) {
				backward.read((char*)&backwardPre[i], sizeof backwardPre[i]);
				idx[backwardPre[i]] = i;
				int rDiv; backward.read((char*)&rDiv, sizeof rDiv);
				backwardSparse.table[i][0] = rDiv;

				if (rDiv > rsite_bp) {
					for (int j = start; j < i && j != -1; ++j) block[backwardPre[j]] = id;
					blockSize[id] = i - start;
					++id;
					start = i;
				}
			}
			// special case where a matching block extends up to the final haplotype
			for (int j = start; j < M; ++j)	block[backwardPre[j]] = id;
			blockSize[id] = M - start;
			
			// initialize forward sparse table
			for (int i = 0; i < M; ++i) forwardSparse.table[i][0] = div[i];

			// build both sparse tables
			forwardSparse.init();
			backwardSparse.init();
			
			// block finding algorithm
			double MI = 0; // mutual information
			start = 0;
			for (int i = 1; i < M; ++i) {
				if (div[i] >= fp) {
					// process reverse blocks in the forward block
					map<int, State> mp;
					for (int j = start; j < i; ++j) {
						id = block[pre[j]];
						if (id == -1) continue;
						mp[id].f_mini = min(mp[id].f_mini, j);
						mp[id].f_maxi = max(mp[id].f_maxi, j);
						mp[id].r_mini = min(mp[id].r_mini, idx[pre[j]]);
						mp[id].r_maxi = max(mp[id].r_maxi, idx[pre[j]]);
						mp[id].IDs.push_back(pre[j]);
						for (int k = 0; k < G; ++k) {
							if (getGap(k, pre[j]) == 0) ++mp[id].zero[k];
							else ++mp[id].one[k];
						}
					}

					// go through all candidates
					for (const auto& p: mp) {
						State state = p.second;
						// compute MI
						double pxy = (double)((int)state.IDs.size()) / M;
						double px = (double)(i - start) / M;
						double py = (double)blockSize[p.first] / M;
						MI += pxy * log2(pxy / px / py);

						if ((int)state.IDs.size() < W) continue; // width too small
						bool flag = false;
						for (int j = 0; j < G; ++j) {
							if (min(state.zero[j], state.one[j]) == 0) flag = true; // no mismatch at this site in the gap
						}
						if (flag) continue;

						int fL = (site - 1) - forwardSparse.query(state.f_mini + 1, state.f_maxi) + 1; // length of forward block
						int rL = rsite - backwardSparse.query(state.r_mini + 1, state.r_maxi) + 1; // length of reverse block

						if (fL >= L && rL >= L) {
							results << site << ' ' << positions[site] << ' ' << fL << ' ' << rL << ' ' << positions[site - fL] << ' ' << positions[site + G + rL - 1] << ' ' << (int)state.IDs.size() << '\n';
							for (int j = 0; j < (int)state.IDs.size(); ++j) resultIDs << ID[state.IDs[j]] << ' ';
							resultIDs << '\n';
						}
					}
					start = i;
				}
			}
			// special case where a matching block extends up to the final haplotype
			// process reverse blocks in the forward block
			map<int, State> mp;
			for (int j = start; j < M; ++j) {
				id = block[pre[j]];
				if (id == -1) continue;
				mp[id].f_mini = min(mp[id].f_mini, j);
				mp[id].f_maxi = max(mp[id].f_maxi, j);
				mp[id].r_mini = min(mp[id].r_mini, idx[pre[j]]);
				mp[id].r_maxi = max(mp[id].r_maxi, idx[pre[j]]);
				mp[id].IDs.push_back(pre[j]);
				for (int k = 0; k < G; ++k) {
					if (getGap(k, pre[j]) == 0) ++mp[id].zero[k];
					else ++mp[id].one[k];
				}
			}
			// go through all candidates
			for (const auto& p: mp) {
				State state = p.second;
				// compute MI
				double pxy = (double)((int)state.IDs.size()) / M;
				double px = (double)(M - start) / M;
				double py = (double)blockSize[p.first] / M;
				MI += pxy * log2(pxy / px / py);

				if ((int)state.IDs.size() < W) continue; // width too small
				bool flag = false;
				for (int j = 0; j < G; ++j) {
					if (min(state.zero[j], state.one[j]) == 0) flag = true; // no mismatch at this site in the gap
				}
				if (flag) continue;

				int fL = (site - 1) - forwardSparse.query(state.f_mini + 1, state.f_maxi) + 1; // length of forward block
				int rL = rsite - backwardSparse.query(state.r_mini + 1, state.r_maxi) + 1; // length of reverse block

				if (fL >= L && rL >= L) {
					results << site << ' ' << positions[site] << ' ' << fL << ' ' << rL << ' ' << positions[site - fL] << ' ' << positions[site + G + rL - 1] << ' ' << (int)state.IDs.size() << '\n';
					for (int j = 0; j < (int)state.IDs.size(); ++j) resultIDs << ID[state.IDs[j]] << ' ';
					resultIDs << '\n';
				}
			}

			resultMI << positions[site] << ' ' << MI << '\n'; 
		}


		// pbwt algorithm
		int u = 0, v = 0, p = site + 1, q = site + 1;
		for (int i = 0; i < M; ++i) {
			int id = pre[i];
			if (div[i] > p) p = div[i];
			if (div[i] > q) q = div[i];
			if (getGap(0, id) == 0) {
				a[u] = id;
				d[u] = p;
				++u;
				p = 0;
			}
			else {
				b[v] = id;
				e[v] = q;
				++v;
				q = 0;
			}
		}
		for (int i = 0; i < u; ++i) {
			pre[i] = a[i];
			div[i] = d[i];
		}
		for (int i = 0; i < v; ++i) {
			pre[u + i] = b[i];
			div[u + i] = e[i];
		}

		in.getline(s, 2 * M + 5000);
		// skip fixed fields
		int offset = 0, tabs = 0; // offset = position in "s" of first sequence - points to the first character after 9 tabs
		while (tabs < 9) {
			if (s[offset] == '\t') ++tabs;
			++offset;
		}

		for (int i = 0; i < M; ++i) {
			assert(s[offset + (i / 2) * 4 + 1] == '|'); // sanity check
			gap[p1][i] = (s[offset + 2 * i] == '0' ? 0 : 1);
		}
		p1 = (p1 + 1) % G;

		if (site % checkpoint == 0) cout << "Checkpoint " << site << endl;
	}

	in.close();
	sites.close();
	meta.close();
	results.close();
	resultIDs.close();
	resultMI.close();

	return 0;
}
