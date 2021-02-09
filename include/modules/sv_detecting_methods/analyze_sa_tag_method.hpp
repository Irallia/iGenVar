#pragma once

#include "structures/aligned_segment.hpp"   // for struct AlignedSegment

/*! \brief Splits a string by a given delimiter and stores substrings in a given container.
 *
 * \param str   string to split
 * \param cont  container for the splitted substrings
 * \param delim delimiter
 */
template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ');

/*! \brief Parse the SA tag from the SAM/BAM alignment of a chimeric/split-aligned read. Build
 *         [aligned_segments](\ref AlignedSegment), one for each alignment segment of the read.
 *
 * \param sa_string         "SA" tag string
 * \param aligned_segments  vector of [aligned_segments](\ref AlignedSegment).
 *
 * \details The SA tag describes the alignments of a chimeric read and is like a small SAM within a SAM file:
 *          "SA:Z:(rname,pos,strand,CIGAR,mapQ,NM;)+"
 *          Each element (in parentheses) represents one alignment segment of the chimeric alignment formatted as
 *          a colon-delimited list.
 *          We add all segments to our candidate list `aligned_segments` and examine them in the following function
 *          `analyze_aligned_segments()`.
 *          For more information about this tag, see the
 *          ([Map Optional Fields Specification](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)).
 */
void retrieve_aligned_segments(std::string sa_string, std::vector<AlignedSegment> & aligned_segments);

/*! \brief Build junctions out of aligned_segments.
 *
 * \param aligned_segments  vector of [aligned_segments](\ref AlignedSegment).
 * \param junctions         vector for storing junctions
 * \param read_name         QNAME field of the SAM/BAM file
 */
void analyze_aligned_segments(const std::vector<AlignedSegment> & aligned_segments,
                              std::vector<Junction> & junctions,
                              const std::string & read_name);