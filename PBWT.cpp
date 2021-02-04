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

int M, N, L; // # of sequences, # of sites, L

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
	ofstream results(writeTo + ".results");

	int checkpoint = atoi(argv[3]);
	L = atoi(argv[4]);

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

	vector<int> pre(M), div(M); // prefix and divergence array
	iota(pre.begin(), pre.end(), 0);
	vector<int> a(M), b(M), d(M), e(M);
	char s[2 * M + 5000]; // assumes fixed fields take up less than 5000 characters
	vector<int> idx(M); // idx[i] = index of sample i in the reverse positional prefix array
	vector<int> backwardPre(M), block(M), blockSize(M + 1); // block[i] = block ID of sample i in the reverse PBWT
	int fp = -1, bp = 0; // pointers for gemomic distance

	for (int site = 0; site + 1 < N; ++site) {
		in.getline(s, 2 * M + 5000);
		// skip fixed fields
		int offset = 0, tabs = 0; // offset = position in "s" of first sequence - points to the first character after 9 tabs
		while (tabs < 9) {
			if (s[offset] == '\t') ++tabs;
			++offset;
		}

		if (site != 0) {
			int rsite = (N - 1) - site - 1; // index of the corresponding reverse site
			backward.seekg((long long)rsite * M * 8);

			// update genomic pointers
			while (positions[site - 1] - (fp != -1 ? positions[fp] : 0) >= L) ++fp; // points to first position < L
			while ((bp < N ? positions[bp] : numeric_limits<int>::max()) - positions[site + 1] < L) ++bp; // points to the first position >= L
			int rsite_bp = (N - 1) - bp;

			// initialize backward sparse table, idx, backwardPre, and block
			int start = -1, id = 0;
			for (int i = 0; i < M; ++i) {
				backward.read((char*)&backwardPre[i], sizeof backwardPre[i]);
				idx[backwardPre[i]] = i;
				int rDiv; backward.read((char*)&rDiv, sizeof rDiv);

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
			
			
			// block finding algorithm
			double MI = 0; // mutual information
			start = 0;
			for (int i = 1; i < M; ++i) {
				if (div[i] >= fp) {
					// process reverse blocks in the forward block
					map<int, vector<int>> mp;
					for (int j = start; j < i; ++j) {
						id = block[pre[j]];
						mp[id].push_back(pre[j]);
					}

					// go through all blocks and compute MI
					for (const auto& p: mp) {
						double pxy = (double)(p.second.size()) / M;
						double px = (double)(i - start) / M;
						double py = (double)blockSize[p.first] / M;
						MI += pxy * log2(pxy / px / py);

					}
					start = i;
				}
			}
			// special case where a matching block extends up to the final haplotype
			// process reverse blocks in the forward block
			map<int, vector<int>> mp;
			for (int j = start; j < M; ++j) {
				id = block[pre[j]];
				mp[id].push_back(pre[j]);
			}

			// go through all candidates and compute MI
			for (const auto& p: mp) {
				double pxy = (double)(p.second.size()) / M;
				double px = (double)(M - start) / M;
				double py = (double)blockSize[p.first] / M;
				MI += pxy * log2(pxy / px / py);
			}

			results << positions[site] << ' ' << MI << '\n'; 
		}


		// pbwt algorithm
		int u = 0, v = 0, p = site + 1, q = site + 1;
		for (int i = 0; i < M; ++i) {
			int id = pre[i];
			if (div[i] > p) p = div[i];
			if (div[i] > q) q = div[i];
			assert(s[offset + (id / 2) * 4 + 1] == '|'); // sanity check
			if (s[offset + (id / 2) * 4 + (id % 2) * 2] == '0') {
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

		if (site % checkpoint == 0) cout << "Checkpoint " << site << endl;
	}

	in.close();
	sites.close();
	meta.close();
	results.close();

	return 0;
}
