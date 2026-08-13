module;
#include <chrono> // fix for GCC 15.2.1
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
export module htslib_wrapper;

import <filesystem>;
import <iostream>;
import <string>;
import <cassert>;

using std::chrono::high_resolution_clock, std::chrono::duration_cast, std::chrono::milliseconds;
using std::filesystem::path, std::filesystem::exists, std::filesystem::last_write_time, std::filesystem::remove;
using std::cout, std::cerr, std::endl, std::flush;
using std::string, std::to_string;

namespace htslib_wrapper {
	export {
		using ::hts_pos_t; // sequence index
	}

	/*
	 * read-only fasta file, possibly bgzipped
	 */
	export class fasta_file {
		private:
			faidx_t *idx;
			hts_pos_t n;

		public:
			fasta_file() = delete;

			/*
			 * index and open a given (gzipped) FASTA file
			 * NB: we follow the convention of checking and generating ".fai"/".gzi" indices
			 */
			fasta_file(const string &fa_path)
			{
				path fastap(fa_path);
				path fastaindex(fa_path + ".fai");
				path gzip_index(fa_path + ".gzi");
				bool compressed = false;
				{
					BGZF *bgzf = bgzf_open(fastap.c_str(), "r");
					if (bgzf->is_compressed) {
						compressed = true;
					}
					bgzf_close(bgzf);
				}

				if (!exists(fastaindex) or (compressed and !exists(gzip_index))) {
					auto start = high_resolution_clock::now();
					cerr << "[htslib_wrapper] index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " or \"" + gzip_index.string() + "\"" : "") << " not found, generating the index" << ((compressed) ? "es" : "") << "..." << flush;
					if (fai_build3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str()) == -1) {
						cerr << "[htslib_wrapper] error: failed to create index" << endl;
						exit(1);
					}
					auto stop = high_resolution_clock::now();
					auto duration = duration_cast<milliseconds>(stop - start);
					cerr << " done" << " (" + to_string(duration.count()) + "ms)" << endl;
				} else if (last_write_time(fastaindex) < last_write_time(fastap) or
						(compressed and last_write_time(gzip_index) < last_write_time(fastap))) {
					auto start = high_resolution_clock::now();
					cerr << "[htslib_wrapper] index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " or \"" + gzip_index.string() + "\" are" : " is") << " older than MSA, regenerating the index" << ((compressed) ? "es" : "") << "..." << flush;
					remove(fastaindex);
					if (compressed and exists(gzip_index))
						remove(gzip_index);
					if (fai_build3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str()) == -1) {
						cerr << "[htslib_wrapper] error: failed to create index" << endl;
						exit(1);
					}
					auto stop = high_resolution_clock::now();
					auto duration = duration_cast<milliseconds>(stop - start);
					cerr << " done" << " (" + to_string(duration.count()) + "ms)" << endl;
				} else {
					cerr << "[htslib_wrapper] index" << ((compressed) ? "es " : " ") << fastaindex << ((compressed) ? " and \"" + gzip_index.string() + "\"" : "") << " found" << endl;
				}

				assert(exists(fastaindex) and (!compressed or exists(gzip_index)));
				if (!(idx = fai_load3(fastap.c_str(), fastaindex.c_str(), gzip_index.c_str(), FAI_NONE))) {
					cerr << "[htslib_wrapper] error: failed to create index" << endl;
					exit(1);
				}

				n = faidx_nseq(idx);
			}

			hts_pos_t nseq() const
			{
				return n;
			}

			string id(hts_pos_t i) const
			{
				assert(i >= 0 and i < n);
				return string(faidx_iseq(idx, i));
			}

			hts_pos_t len(hts_pos_t i) const
			{
				assert(i >= 0 and i < n);
				return faidx_seq_len64(idx, faidx_iseq(idx, i));
			}

			string fetch(hts_pos_t i) const
			{
				assert(i >= 0 and i < n);
				hts_pos_t out_len;
				char *buffer = faidx_fetch_seq64(idx, faidx_iseq(idx, i), 0, faidx_seq_len64(idx, faidx_iseq(idx, i)) - 1, &out_len);
				assert(out_len == faidx_seq_len64(idx, faidx_iseq(idx, i)));
				const string res(buffer);
				assert((hts_pos_t)res.length() == out_len);
				free(buffer);
				return res;
			}

			~fasta_file() {
				fai_destroy(idx);
			}
	};
}
