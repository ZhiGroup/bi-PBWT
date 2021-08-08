#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <istream>
#include <cassert>
#include <fstream>

using namespace std;

int M, W;
long double L;

// skip meta-info lines and header line
string skipMeta(ifstream& in) {
	string s;
	while (getline(in, s)) {
		if ((int)s.size() < 2 || s[0] != '#' || s[1] != '#') break;
	}

	return s;	
}

// find M using vcf's info field - resets input stream pointer when finished
void getM(ifstream& in) {
	int startPos = in.tellg(); 
	string s; getline(in, s);
	int tabs = 0;
	for (int i = 0; i < (int)s.size(); ++i) {
		if (s[i] == '\t') ++tabs;
		if (tabs >= 9 && s[i] == '|') ++M;
	}
	in.seekg(startPos);
	M *= 2;
}

int main(int argc, char* argv[]) {
	ios_base::sync_with_stdio(0); cin.tie(0);

	ifstream in(argv[1]), rMap(argv[5]);
	ofstream out(string(argv[2]) + ".blocks");
	L = atof(argv[3]), W = atoi(argv[4]);
	
	string header = skipMeta(in);
	getM(in);

	// input sample IDs
	vector<string> ID(M);
	stringstream ss(header);
	for (int i = 0; i < 9; ++i) getline(ss, ID[0], '\t'); // skip fixed columns, assumes 9 columns (FORMAT column) 
	for (int i = 0; i < M / 2; ++i) {
		getline(ss, ID[2 * i], '\t');
		ID[2 * i + 1] = ID[2 * i] + "-1";
		ID[2 * i] += "-0";
	}

	int site = 0;
	vector<int> pre(M), div(M);
	for (int i = 0; i < M; ++i) pre[i] = i;
	vector<int> a(M), b(M), d(M), e(M);
	vector<int> positions;
	vector<long double> map;
	char s[2 * M + 5000];

	while (in.getline(s, 2 * M + 5000)) {
		// skip fixed fields
		string position;
		int offset = 0, tabs = 0;
		while (tabs < 9) {
			if (tabs == 1 && s[offset] != '\t') position += s[offset];
			if (s[offset] == '\t') ++tabs;
			++offset;
		}
		positions.push_back(stoi(position));

		int trash;
		long double genetic_pos;
		rMap >> trash >> genetic_pos;
		map.push_back(genetic_pos);

		if (site != 0) {
			int start = 0, maxi = 0; 
			for (int i = 1; i < M; ++i) {
				if (map[div[i]] > map[site] - L) {
					if (i - start >= W) {
						out << positions[maxi] << ' ' << positions[site - 1] << ' ' << map[maxi] << ' ' << map[site - 1] << ' ' << (i - start);
						for (int j = start; j < i; ++j) {
							out << ' ' << ID[pre[j]];
						}
						out << '\n';
					}
					start = i;
					maxi = 0;
				}
				else maxi = max(maxi, div[i]);
			}
			if (M - start >= W) {
				out << positions[maxi] << ' ' << positions[site - 1] << ' ' << map[maxi] << ' ' << map[site - 1] << ' ' << (M - start);
				for (int j = start; j < M; ++j) {
					out << ' ' << ID[pre[j]];
				}
				out << '\n';
			}
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

		++site;
	}

	in.close();
	out.close();
}
