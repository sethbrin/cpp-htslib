#include "easehts/genome_loc.h"
#include "easehts/sam_sequence_dictionary.h"

#include <string>
#include <vector>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

namespace ncic {
namespace easehts {


const GenomeLoc GenomeLoc::kUnmapped = GenomeLoc();
const GenomeLoc GenomeLoc::kWholeGenome = GenomeLoc("all");

bool GenomeLocParser::ValidateGenomeLoc(const std::string& contig, int start,
                                        int stop, bool must_be_on_reference) const {
  if (validation_level_ == ValidationLevel::NONE) {
    return true;
  } else {
    if (stop < start) {
      WARN(utils::StringFormatCStr("The stop position %d is less than start %d in contig %s",
                               stop, start, contig.c_str()));
      return false;
    }
    if (must_be_on_reference &&
        !ref_header_.HasSequence(contig)) {
      WARN(utils::StringFormatCStr("The reference has no contig %s\n", contig.c_str()));
      return false;
    }

    if (must_be_on_reference) {
      if (start < 1) {
        WARN(utils::StringFormatCStr("The start position %d is less than 1\n", start));
        return false;
      }
      if (stop < 1) {
        WARN(utils::StringFormatCStr("The stop position %d is less than 1\n", stop));
        return false;
      }
      int contig_size = ref_header_.GetSequenceLen(contig);
      if (start > contig_size || stop > contig_size) {
        WARN(utils::StringFormatCStr("The genome loc coordinates %d-%d exceed the contig size (%d)\n", start, stop, contig_size));
        return false;
      }

    }
  }
  return true;
}
int GenomeLocParser::GetContigIndex(const std::string& contig) const {
  return ref_header_.GetSequenceId(contig);
}

std::vector<GenomeLoc> IntervalUtils::IntervalFileToList(
    const GenomeLocParser& gl_parser,
    const std::string& filename) {
#define _read_token(_p) (_p); for (; *(_p) && *(_p) != '\t'; ++(_p)) {};
#define _read_token_aux(_p) (_p); for (; *(_p) && *(_p) != '\t'; ++(_p)); *(_p)++ = 0
#define _parse_err(cond, msg) do { if ((cond) && hts_verbose >= 1) { fprintf(stderr, "[E::%s] " msg "\n", __func__); goto err_ret; } } while (0)
#define _parse_warn(cond, msg) if ((cond) && hts_verbose >= 2) fprintf(stderr, "[W::%s] " msg "\n", __func__)

  std::vector<GenomeLoc> genome_locs;
  samFile* fp = sam_open(filename.c_str(), "r");
  SAMSequenceDictionary header = ReadHeader(fp);

  std::string name;
  std::string contig;
  int start;
  int stop;
  bool negative;
  // the fileds:
  // sequence\tstart\tend\t(+/-)\tname
  do {
    char *p = fp->line.s, *q;
    // sequence
    q = _read_token(p);
    _parse_warn(p - q < 1, "empty sequence");
    contig = std::string(q, p-q);

    // start
    start = ::strtol(p, &p, 10);
    if (*p++ != '\t') goto err_ret; // malformated start

    // end
    stop = ::strtol(p, &p, 10);
    if (*p++ != '\t') goto err_ret; // malformated end

    // +/-
    negative = true;
    if (*p == '+') {
      negative = false;
    } else if (*p == '-') {
      negative = true;
    } else {
      goto err_ret;
    }
    if (*(++p) != '\t') goto err_ret; // malformated end

    // name
    q = _read_token_aux(p);
    _parse_warn(p - q < 1, "empty sequence");
    name = std::string(q, p-q);

    _parse_warn(!header.HasSequence(contig), "Ignore interval for unknown reference");

    if (gl_parser.ValidateGenomeLoc(contig, start, stop, true)) {
      genome_locs.emplace_back(contig, gl_parser.GetContigIndex(contig), start, stop);
    } else {
      _parse_warn(true, "Ignore invalid genome loc");
    }
    continue;

err_ret:
    fprintf(stderr, "[W::%s] malformated line \n", __func__);
  } while (::hts_getline(fp, KS_SEP_LINE, &fp->line) >= 0);

#undef _read_token
#undef _read_token_aux
#undef _parse_err
#undef _parse_warn

  ::hts_close(fp);
  return genome_locs;
}


SAMSequenceDictionary IntervalUtils::ReadHeader(htsFile* fp) {
  return SAMSequenceDictionary::ReadHeader(fp);
}

}// easehts
} // ncic
