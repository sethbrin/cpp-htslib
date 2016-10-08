//
// Created by zp on 8/24/16.
//

#ifndef EASEHTSLIB_SAM_TEXT_HEADER_CODEC_H
#define EASEHTSLIB_SAM_TEXT_HEADER_CODEC_H

#include "sam_file_header.h"

#include <string>

#include <htslib/sam.h>
#include <map>

namespace ncic {

namespace easehts {

// For the detail meaning, see https://samtools.github.io/hts-specs/SAMv1.pdf
typedef enum HeaderRecordType {
  HD, SQ, RG, PG, CO, NUL // invalid header type
} HeaderRecordType;

class ParsedHeaderLine;
// sam file header codec
class SAMTextHeaderCodec {
 public:
  SAMFileHeaderPtr Parse(htsFile *fp);

  /**
   * Takes a header line as a String and converts it into a HeaderRecordType, and a map of key:value strings.
   * If the line does not contain a recognized HeaderRecordType, then the line is considered invalid, and will
   * not have any key:value pairs.
   */
  class ParsedHeaderLine {
   public:
    explicit ParsedHeaderLine(const std::string& line);

    const HeaderRecordType& GetHeaderRecordType() const{
      return header_type_;
    }

    bool HasKey(const std::string& key) const {
      if (key_value_pairs_.find(key) != key_value_pairs_.end()) {
        return true;
      } else {
        return false;
      }
    }

    const std::string& GetValue(const std::string& key) const {
      if (!HasKey(key)) {
        fprintf(stderr, "You should first check if has key\n");
      }

      return key_value_pairs_.at(key);
    }

    void RemoveKey(const std::string& key) {
      key_value_pairs_.erase(key);
    }

    const std::map<std::string, std::string>& GetKeyValuePairs() {
      return key_value_pairs_;
    };
   private:
    static std::map<std::string, HeaderRecordType> InitHeaderTypeMap() {
      std::map<std::string, HeaderRecordType> header_type_map;

      header_type_map["@HD"] = HD;
      header_type_map["@SQ"] = SQ;
      header_type_map["@RG"] = RG;
      header_type_map["@PG"] = PG;
      header_type_map["@CO"] = CO;

      return header_type_map;
    };

    /**
     * parse the header type of the line
     */
    void ParseHeaderType(const std::string& type);

   private:
    const static std::map<std::string, HeaderRecordType> HEADER_TYPE_MAP;

    // Store the line attributes
    std::map<std::string, std::string> key_value_pairs_;
    HeaderRecordType header_type_;
  };

 private:
  /** read next header line, if has next, return true and set to line,
   * otherwise return false
   *
   * @param
   *    fp: the sam file
   *    line: the header line
   *
   *  @return
   *     true if has next header line
   */
  bool AdvanceLine(htsFile *fp, std::string* line);

  /**
  * parse SQ line, it will remove the parsed_header_line attribute if needed
  */
  void ParseSQLine(const SAMFileHeaderPtr& sam_header,
                   ParsedHeaderLine *parsed_header_line);

 private:
  static const char HEADER_LINE_START = '@';
  static const char FIELD_SEPARATOR_CHAR = '\t';
  static const char TAG_KEY_VALUE_SEPARATOR_CHAR = ':';

  static const std::string FIELD_SEPARATOR;
  static const std::string TAG_KEY_VALUE_SEPARATOR;
};

} // easehts
} // ncic
#endif //EASEHTSLIB_SAM_TEXT_HEADER_CODEC_H
