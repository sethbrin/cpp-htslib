#include "reference_sequence.h"

#include <set>
#include <string>

namespace ncic {
namespace easehts {

AbstractFastaSequenceFile::AbstractFastaSequenceFile(const std::string& file_path) :
  file_path_(file_path) {
  std::string dictionary_name = AbstractFastaSequenceFile::FindSequenceDictionaryFileName(file_path_);
  if (!dictionary_name.empty()) {
    samFile* fp = sam_open(dictionary_name.c_str(), "r");
    ref_header_ = SAMSequenceDictionary::ReadHeader(fp);
  }
}

std::string AbstractFastaSequenceFile::FindSequenceDictionaryFileName(const std::string& file_path) {
  bool filetype_supported = false;
  std::string dictionary_name;
  // the name with ext, such as file_path: a.fasta
  // then dictionary_name: a.dict, dictionary_name_ext: a.fasta.dict
  // first check dictionary_name exists
  std::string dictionary_name_ext;
  for (auto& ext : AbstractFastaSequenceFile::kFastaExtensions) {
    if (utils::EndsWith(file_path, ext)) {
      dictionary_name_ext = file_path + ".dict";
      dictionary_name = file_path.substr(0, file_path.size() - ext.size()) + ".dict";
      filetype_supported = true;
      break;
    }
  }

  if (!filetype_supported) {
    ERROR(utils::StringFormatCStr("File is not a supported references file type: %s", file_path.c_str()));
  } else {
    if (utils::FileExists(dictionary_name)) {
      return dictionary_name;
    } else {
      return utils::FileExists(dictionary_name_ext) ? dictionary_name_ext : std::string();
    }
  }
}

const std::set<std::string> AbstractFastaSequenceFile::kFastaExtensions = {
  ".fasta",
  ".fasta.gz",
  ".fa",
  ".fa.gz",
  ".fna",
  ".fna.gz",
  ".txt",
  ".txt.gz"
};

} // easehts
} //ncic
