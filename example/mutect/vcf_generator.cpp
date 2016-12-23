#include "vcf_generator.h"

#include <easehts/genotype.h>

namespace ncic {
namespace mutect {

const std::string VCFGenerator::kTumorName = "MG225_tumor";
const std::string VCFGenerator::kNormalName = "MG225_normal";

void VCFGenerator::WriteHeader(
    const easehts::SAMSequenceDictionary& dict,
    const std::string assemby) {
  // use default version string set in htslib
  easehts::VCFHeader* header = writer_.GetHeader();
  header->SetVersion(header->GetVersion());

  // Add Filter header
  header->AddHeaderLine(
      "##FILTER=<ID=PASS,Description=\"Accept as a confident"
      " somatic mutation\">");
  header->AddHeaderLine(
      "##FILTER=<ID=REJECT,Description=\"Rejected as a confident"
      " somatic mutation\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths"
      " for the ref and alt alleles in the order listed\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=BQ,Number=A,Type=Float,Description=\"Average base quality"
      " for reads supporting alleles\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read "
      "depth (reads with MQ=255 or with bad mates are filtered)\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=FA,Number=A,Type=Float,Description=\"Allele fraction "
      "of the alternate allele with regard to reference\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,"
      "Description=\"Genotype Quality\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, "
      "Phred-scaled likelihoods for genotypes as defined in the "
      "VCF specification\">");
  header->AddHeaderLine(
      "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status "
      "relative to non-adjacent Normal,0=wildtype,1=germline,"
      "2=somatic,3=LOH,4=post-transcriptional modification,5=unknown\">");
  header->AddHeaderLine(
      "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP Membership\">");
  header->AddHeaderLine(
      "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total "
      "Mapping Quality Zero Reads\">");
  header->AddHeaderLine(
      "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">");
  header->AddHeaderLine(
      "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant "
      "type, can be SNP, INS or DEL\">");
  header->AddContigs(dict.GetRawHeader(), assemby);
  header->AddSample(kTumorName);
  header->AddSample(kNormalName);
  header->FinishSample(); // update inner

  writer_.WriteHeader();
}

//void VCFGenerator::Add(const CandidateMutation& candidate) {
//  easehts::VCFHeader* header = writer_.GetHeader();
//  easehts::VariantContext vc;
//  vc.SetLocation(candidate.location);
//
//  easehts::GenotypeBuilder tumor_genetype(
//      candidate.tumor_sample_name,
//      {easehts::Allele::Create(candidate.ref_allele, true),
//      easehts::Allele::Create(candidate.alt_allele)});
//
//  tumor_genetype
//    .AD({candidate.initial_tumor_ref_counts,
//        candidate.initial_tumor_alt_counts})
//    .Attribute("FA", std::to_string(candidate.tumor_F))
//    .DP(candidate.initial_tumor_read_depth);
//
//  easehts::GenotypeBuilder normal_genetype(
//      candidate.normal_sample_name,
//      {easehts::Allele::Create(candidate.ref_allele, true)});
//
//  normal_genetype
//    .AD({candidate.initial_normal_ref_counts,
//        candidate.initial_normal_alt_counts})
//    .Attribute("FA", std::to_string(candidate.normal_F))
//    .DP(candidate.initial_normal_read_depth)
//    .Attribute(easehts::VCFConstants::kRmsBaseQualityKey, ".");
//
//
//  if (candidate.dbsnp_VC != nullptr) {
//    vc.SetId(candidate.dbsnp_VC->GetId(), header->GetRawHeader());
//    vc.UpdateFlagInfo(easehts::VCFConstants::kDbsnpKey,
//                      header->GetRawHeader());
//  }
//
//  if (candidate.initial_tumor_alt_counts > 0) {
//    vc.UpdateInfo(easehts::VCFConstants::kRmsBaseQualityKey,
//                     &candidate.initial_tumor_alt_counts, 1,
//                     header->GetRawHeader());
//  }
//
//  if (!candidate.rejected) {
//    vc.Filter("PASS", header->GetRawHeader());
//    vc.UpdateFlagInfo("SOMATIC", header->GetRawHeader());
//    vc.UpdateStringInfo("VT", "SNP", header->GetRawHeader());
//
//    tumor_genetype.Attribute("SS", "2");
//    normal_genetype.Attribute("SS", "0");
//  } else {
//    vc.Filter("REJECT", header->GetRawHeader());
//  }
//
//  writer_.Add(vc);
//}

void VCFGenerator::Add(const CandidateMutation& candidate) {
  easehts::VCFHeader* header = writer_.GetHeader();
  std::unique_ptr<easehts::VariantContext> vc(new easehts::VariantContext());
  vc->SetLocation(candidate.location);
  std::string ref_alt;
  ref_alt.push_back(candidate.ref_allele);
  ref_alt.push_back(',');
  ref_alt.push_back(candidate.alt_allele);
  vc->UpdateAlleleStr(ref_alt, header->GetRawHeader());

  int sample_cnt = 2;
  // two samples
  std::unique_ptr<int[]> tmpi(
      new int[2 * sample_cnt]);

  // Add Format
  // Add genetype
  tmpi[0] = easehts::VariantContext::GenotypePhaseInt(0);
  tmpi[1] = INT32_VECTOR_MISSING;
  tmpi[2] = easehts::VariantContext::GenotypePhaseInt(0);
  tmpi[3] = easehts::VariantContext::GenotypePhaseInt(1);
  vc->UpdateGenetype(tmpi.get(), 2*sample_cnt, header->GetRawHeader());

  // Add AD
  tmpi[0] = candidate.initial_normal_ref_counts;
  tmpi[1] = candidate.initial_normal_alt_counts;
  tmpi[2] = candidate.initial_tumor_ref_counts;
  tmpi[3] = candidate.initial_tumor_alt_counts;
  vc->UpdateFormat("AD", tmpi.get(), 2 * sample_cnt, header->GetRawHeader());

  // Add BQ
  if (candidate.initial_tumor_alt_counts > 0) {
    tmpi[0] = INT32_MISSING;
    tmpi[1] = candidate.initial_tumor_alt_quality_sum /
      candidate.initial_tumor_alt_counts;
    vc->UpdateFormat("BQ", tmpi.get(), sample_cnt, header->GetRawHeader());
  }


  // Add DP
  tmpi[0] = candidate.initial_normal_read_depth;
  tmpi[1] = candidate.initial_tumor_read_depth;
  vc->UpdateFormat("DP", tmpi.get(), sample_cnt, header->GetRawHeader());

  std::unique_ptr<float[]> tmpf(
      new float[2 * sample_cnt]);
  // Add FA
  tmpf[0] = candidate.normal_F;
  tmpf[1] = candidate.tumor_F;
  vc->UpdateFormat("FA", tmpf.get(), sample_cnt, header->GetRawHeader());

  if (candidate.dbsnp_VC != nullptr) {
    vc->SetId(candidate.dbsnp_VC->GetId(), header->GetRawHeader());
    vc->UpdateFlagInfo(easehts::VCFConstants::kDbsnpKey,
                      header->GetRawHeader());
  }


  if (!candidate.rejected) {
    vc->Filter("PASS", header->GetRawHeader());
    vc->UpdateFlagInfo("SOMATIC", header->GetRawHeader());
    vc->UpdateStringInfo("VT", "SNP", header->GetRawHeader());

    // Add SS
    tmpi[0] = 0;
    tmpi[1] = 2;
    vc->UpdateFormat("SS", tmpi.get(), sample_cnt, header->GetRawHeader());
  } else {
    vc->Filter("REJECT", header->GetRawHeader());
  }

  writer_.Add(vc);
}


} // mutect
} // ncic
