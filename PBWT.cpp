#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <numeric>
#include <map>
#include <cassert>
#include <cstring>
#include <algorithm>

using namespace std;

int M, N, L, W, G;  // # of sequences, # of sites, length, width, gap size

struct SparseTable {
	static const int B = 30;
	int N, blocks;
	vector<int> v, mask;
	vector<vector<int>> table;

	SparseTable(vector<int> _v) {
		this->v = _v;
		this->N = (int)v.size();
		blocks = N / B;
		mask = vector<int>(N);	
		table = vector<vector<int>>(blocks, vector<int>(msb(blocks) + 1));

		int cur = 0; // sliding mask
		for (int i = 0; i < N; ++i) {
			cur = (cur << 1) & ((1 << B) - 1);
			while (cur > 0 && max(v[i], v[i - msb(lsb(cur))]) == v[i]) cur ^= lsb(cur); 
			cur |= 1;
			mask[i] = cur;
		}

		for (int i = 0; i < blocks; ++i) table[i][0] = mini_query(B * i + B - 1);
		for (int j = 1; (1 << j) <= blocks; ++j) {
			for (int i = 0; i + (1 << j) - 1 < blocks; ++i) {
				table[i][j] = max(table[i][j - 1], table[i + (1 << (j - 1))][j - 1]);
			}
		}
	}

	// least significant set bit
	int lsb(int num) {return num & -num;}

	// index of most significant set bit
	int msb(int num) {return __builtin_clz(1) - __builtin_clz(num);}

	int mini_query(int r, int len = B) {
		return v[r - msb(mask[r] & ((1 << len) - 1))];
	}

	int query(int l, int r) {
		if (r - l + 1 <= B) return mini_query(r, r - l + 1);
		int ret = max(mini_query(l + B - 1), mini_query(r));
		int blockL = l / B + 1, blockR = r / B - 1;
		if (blockL <= blockR) {
			int j = msb(blockR - blockL + 1);
			ret = max({ret, table[blockL][j], table[blockR - (1 << j) + 1][j]});
		}
		return ret;
	}
};

struct VCFReader {
	ifstream vcf;
	int G, M, p1 = 0;
	vector<vector<int>> gap; // stores haplotype data in the gap
	vector<string> ID;

	VCFReader(string file, int _G, int _M) {
		vcf = ifstream(file);
		G = max(_G, 1), M = _M;
		gap = vector<vector<int>>(G, vector<int>(M));
		ID = vector<string>(M);
		preprocess();
		initGap();
	}
	
	// gets haplotype IDs and moves input stream pointer to start of raw data
	void preprocess() { 
		// skip meta-info lines and get header line
		string header;
		while (getline(vcf, header)) {
			if ((int)header.size() < 2 || header[0] != '#' || header[1] != '#') break;
		}

		// input sample IDs
		stringstream ss(header);
		for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
		for (int i = 0; i < M / 2; ++i) {
			getline(ss, ID[2 * i], '\t');
			ID[2 * i + 1] = ID[2 * i] + "-1";
			ID[2 * i] += "-0";
		}
	}

	// initializes sliding window for the gap
	void initGap() {
		for (int i = 0; i < G; ++i) nextSite();
	}

	// reads the next site in the VCF file
	void nextSite() {
		char s[2 * M + 5000]; // assumes fixed fields take up less than 5000 characters
		vcf.getline(s, 2 * M + 5000);
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

	int getGap(int g, int idx) {
		int jump = G - g;
		return gap[(p1 + G - jump) % G][idx];
	}

	void close() {vcf.close();}
};


void countingSort(vector<vector<int>>& v, int idx) {
	vector<vector<vector<int>>> table(M + 1);
	for (int i = 0; i < M; ++i) {
		table[v[i][idx]].push_back(v[i]);
	}
	int p = 0;
	for (int i = 0; i <= M; ++i) {
		for (int j = 0; j < (int)table[i].size(); ++j) {
			v[p++] = table[i][j];
		}
	}
}

void processBlock(vector<vector<int>>& link, int start, int end, vector<int>& idx, vector<int>& rIdx, SparseTable& forwardSparse, SparseTable& backwardSparse, int site, int rsite, vector<int>& positions, vector<string>& ID, ofstream& blocks, ofstream& blockIDs, double& MI, vector<int>& blockSize, vector<int>& rBlockSize, VCFReader& vcf) { // [start, end)
	// compute MI
	double pxy = (double)(end - start) / M;
	double px = (double)blockSize[link[start][1]] / M;
	double py = (double)rBlockSize[link[start][2]] / M;
	MI += pxy * log2(pxy / px / py);
	
	if (end - start < W) return; // width too small
	
	int f_mini = M - 1, f_maxi = 0, r_mini = M - 1, r_maxi = 0;
	vector<int> zero(G), one(G);
	for (int j = start; j < end; ++j) {
		int id = link[j][0];
		f_mini = min(f_mini, idx[id]);
		f_maxi = max(f_maxi, idx[id]);
		r_mini = min(r_mini, rIdx[id]);
		r_maxi = max(r_maxi, rIdx[id]);
		for (int k = 0; k < G; ++k) {
			if (vcf.getGap(k, id) == 0) ++zero[k];
			else ++one[k];
		}
	}

	// check if every site in the gap has a mismatch
	bool flag = false;
	for (int k = 0; k < G; ++k) {
		if (min(zero[k], one[k]) == 0) flag = true; // no mismatch at this site in the gap
	}
	if (flag) return;

	int fL = (site - 1) - forwardSparse.query(f_mini + 1, f_maxi) + 1; // length of forward block
	int rL = rsite - backwardSparse.query(r_mini + 1, r_maxi) + 1; // length of reverse block

	blocks << site << ' ' << positions[site] << ' ' << fL << ' ' << rL << ' ' << positions[site - fL] << ' ' << positions[site + G + rL - 1] << ' ' << (end - start) << '\n';
	for (int j = start; j < end; ++j) blockIDs << ID[link[j][0]] << ' ';
	blockIDs << '\n';
}

int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0); cin.tie(0);

	string writeTo = string(argv[2]);
	ifstream backward(writeTo + ".rpbwt"), sites(writeTo + ".sites"), meta(writeTo + ".meta");
	ofstream blocks(writeTo + ".blocks"), blockIDs(writeTo + ".IDs"), resultMI(writeTo + ".MI");

	int checkpoint = atoi(argv[3]);
	L = atoi(argv[4]), W = atoi(argv[5]), G = atoi(argv[6]);

	// retrieve M and N from meta file
	meta >> M >> N;

	VCFReader vcf(string(argv[1]), G, M);

	// input chromosome site positions
	vector<int> positions(N);
	for (int i = 0; i < N; ++i) sites >> positions[i];

	vector<int> pre(M), div(M), backwardPre(M), backwardDiv(M); // prefix and divergence arrays for forward and backward PBWT
	iota(pre.begin(), pre.end(), 0);
	vector<int> a(M), b(M), d(M), e(M);
	vector<int> idx(M), rIdx(M); // idx[i] = index of sample i in the positional prefix array; r = reverse
	vector<int> block(M), blockSize(M + 1), rBlock(M), rBlockSize(M + 1); // block[i] = block ID of sample i in the reverse PBWT; block IDs go from [1, M]

	for (int site = 0; site + G < N; ++site) {
		if (site != 0) {
			int rsite = (N - 1) - site - G; // index of the corresponding reverse site
			backward.seekg((long long)rsite * M * 8);

			// initialize rIdx, backwardPre, rBlock, and rBlockSize
			int start = -1, id = 0;
			for (int i = 0; i < M; ++i) {
				backward.read((char*)&backwardPre[i], sizeof backwardPre[i]);
				rIdx[backwardPre[i]] = i;
				int rDiv; backward.read((char*)&rDiv, sizeof rDiv);
				backwardDiv[i] = rDiv;
				rDiv = (N - 1) - rDiv; // get forward index for position comparision

				if ((string(argv[7]) == "0" && positions[rDiv] < positions[site + (G - 1)] + L) || (string(argv[7]) == "1" && rDiv < site + (G - 1) + L)) {
					for (int j = start; j < i && j != -1; ++j) rBlock[backwardPre[j]] = id;
					rBlockSize[id] = i - start;
					++id;
					start = i;
				}
			}
			// special case where a matching block extends up to the final haplotype
			for (int j = start; j < M; ++j)	rBlock[backwardPre[j]] = id;
			rBlockSize[id] = M - start;

			// initialize idx, block, and blockSize
			start = -1, id = 0;
			for (int i = 0; i < M; ++i) {
				idx[pre[i]] = i;
				if ((string(argv[7]) == "0" && positions[div[i]] > positions[site] - L) || (string(argv[7]) == "1" && div[i] > site - L)) {
					for (int j = start; j < i && j != -1; ++j) block[pre[j]] = id;
					blockSize[id] = i - start;
					++id;
					start = i;
				}
			}
			// special case where a matching block extends up to the final haplotype
			for (int j = start; j < M; ++j) block[pre[j]] = id;
			blockSize[id] = M - start;

			SparseTable forwardSparse(div), backwardSparse(backwardDiv); // build sparse tables

			// Algorithm 2 - block matching
			vector<vector<int>> link(M, vector<int>(3)); // [sample ID, forward block ID, reverse block ID]
			for (int i = 0; i < M; ++i) {
				link[i][0] = i, link[i][1] = block[i], link[i][2] = rBlock[i];
			}
			// radix sort
			countingSort(link, 2);
			countingSort(link, 1);

			double MI = 0; // mutual information
			start = 0;
			for (int i = 1; i < M; ++i) {
				if (link[i][1] != link[i - 1][1] || link[i][2] != link[i - 1][2]) {
					processBlock(link, start, i, idx, rIdx, forwardSparse, backwardSparse, site, rsite, positions, vcf.ID, blocks, blockIDs, MI, blockSize, rBlockSize, vcf);
					start = i;	
				}
			}
			processBlock(link, start, M, idx, rIdx, forwardSparse, backwardSparse, site, rsite, positions, vcf.ID, blocks, blockIDs, MI, blockSize, rBlockSize, vcf);

			resultMI << positions[site] << ' ' << MI << '\n'; 
		}

		// pbwt algorithm
		int u = 0, v = 0, p = site + 1, q = site + 1;
		for (int i = 0; i < M; ++i) {
			int id = pre[i];
			if (div[i] > p) p = div[i];
			if (div[i] > q) q = div[i];
			if (vcf.getGap(0, id) == 0) {
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

		vcf.nextSite(); // move input stream pointer

		if (site % checkpoint == 0) cout << "Checkpoint " << site << endl;
	}

	vcf.close();
	sites.close();
	meta.close();
	blocks.close();
	blockIDs.close();
	resultMI.close();

	return 0;
}
