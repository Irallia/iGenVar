#include <gtest/gtest.h>

#include <array>
#include <vector>
#include <deque>
#include <forward_list>
#include <list>
#include <stack>
#include <queue>
#include <set>
#include <map>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"    // for the split read method
#include "modules/sv_detection_methods/analyze_sa_tag_method.hpp"   // for the cigar string method

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

/* -------- detection methods tests -------- */

// TODO (irallia): implement test cases

TEST(junction_detection, cigar_string_method_simple)
{
    std::string const read_name = "read021";
    std::string const chromosome = "chr1";
    int32_t const query_start_pos = 1;
    std::vector<seqan3::cigar> cigar_string = {{5, 'S'_cigar_operation}, {15, 'M'_cigar_operation},
                                               {6, 'D'_cigar_operation}, {30, 'M'_cigar_operation}}; //5S15M6D24M
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};
    std::vector<Junction> junctions_res{};
    uint64_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 16, strand::forward};
    Breakend new_breakend_2 {chromosome, 22, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name}};

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
    }
}

TEST(junction_detection, split_read_method_simple)
{
    std::string const read_name = "read021";
    seqan3::sam_flag const flag{0u};
    std::string const chromosome = "chr1";
    int32_t const pos = 80;
    uint8_t const mapq = 60;
    std::vector<seqan3::cigar> cigar_string = {{44, 'M'_cigar_operation}, {6, 'S'_cigar_operation}}; //44M6S
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};
    std::string const sa_tag = "chr1,108,+,44S6M,60,0;";
    std::vector<Junction> junctions_res{};
    uint64_t const min_var_length = 4;
    analyze_sa_tag(read_name, flag, chromosome, pos, mapq, cigar_string, seq, sa_tag, junctions_res);

    Breakend new_breakend_1 {chromosome, 108, strand::reverse};
    Breakend new_breakend_2 {chromosome, 124, strand::reverse};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name}};


    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
        // For debugging #include <seqan3/core/debug_stream.hpp> and use:
        // seqan3::debug_stream << "-----------------------------------------------------------------------------------\n"
        //                      << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1()) << ": \n"
        //                      << junctions_expected_res[i].get_mate1() << " == " << junctions_res[i].get_mate1() << "\n"
        //                      << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2()) << ": \n"
        //                      << junctions_expected_res[i].get_mate2() << " == " << junctions_res[i].get_mate2() << "\n";
    }
}

// TODO (irallia): these tests should be implemented when the associated detection methods are implemented

// TEST(junction_detection, read_pairs_method_simple)
// {
//     testing::internal::CaptureStdout();

//     std::vector<Junction> resulting_junctions{};
//     analyze_read_pairs();

//     std::vector<Junction> expected_junctions{};
//     EXPECT_EQ(expected_junctions, resulting_junctions);
// }

// TEST(junction_detection, read_depth_method_simple)
// {
//     testing::internal::CaptureStdout();

//     std::vector<Junction> resulting_junctions{};
//     analyze_read_depth();

//     std::vector<Junction> expected_junctions{};
//     EXPECT_EQ(expected_junctions, resulting_junctions);
// }

/* ---- Tests of single functions: ---- */

/* analyze_sa_tag_method.cpp/hpp */

TEST(junction_detection, split_string)
{
    std::vector<std::string> test_strings{"a,a,a,aa", "b;  b;b;  bbb", "c c cccc"};

    // Container std::vector
    {
        {
            std::vector<std::string> resulting_strings{};
            split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
            std::vector<std::string> expected{"a", "a", "a", "aa"};
            EXPECT_EQ(expected, resulting_strings);
        }
        {
            std::vector<std::string> resulting_strings{};
            split_string(test_strings[1], resulting_strings, ';');      // split "b;  b;b;  bbb" with ';'
            std::vector<std::string> expected{ "b", "  b", "b", "  bbb" };
            EXPECT_EQ(expected, resulting_strings);
        }
        {
            std::vector<std::string> resulting_strings1{};
            split_string(test_strings[2], resulting_strings1, ' ');     // split "c c cccc"     with ' '
            std::vector<std::string> expected{"c", "c", "cccc"};
            EXPECT_EQ(expected, resulting_strings1);

            std::vector<std::string> resulting_strings2{};
            split_string(test_strings[2], resulting_strings2);          // split "c c cccc"     with default ' '
            EXPECT_EQ(resulting_strings1, resulting_strings2);
        }
    }
    // Other containers ..., std::list, std::stack, std::set, ...
    {
        // {
        //     std::array<std::string, 4> resulting_strings{};
        //     split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
        //     std::array<std::string, 4> expected{"a", "a", "a", "aa"};
        //     EXPECT_EQ(expected, resulting_strings);
        // }
        // {
        //     std::deque<std::string> resulting_strings{};
        //     split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
        //     std::deque<std::string> expected{"a", "a", "a", "aa"};
        //     EXPECT_EQ(expected, resulting_strings);
        // }
        // {
        //     std::list<std::string> resulting_strings{};
        //     split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
        //     std::list<std::string> expected{"a", "a", "a", "aa"};
        //     EXPECT_EQ(expected, resulting_strings);
        // }
        // {
        //     std::stack<std::string> resulting_strings{};
        //     split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
        //     std::stack<std::string> expected;
        //     expected.push("a");
        //     expected.push("a");
        //     expected.push("a");
        //     expected.push("aa");
        //     // std::stack<std::string> expected{"a", "a", "a", "aa"};
        //     EXPECT_EQ(expected, resulting_strings);
        // }
        // {
        //     std::queue<std::string> resulting_strings{};
        //     split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
        //     std::queue<std::string> expected{};
        //     expected.push("a");
        //     expected.push("a");
        //     expected.push("a");
        //     expected.push("aa");
        //     EXPECT_EQ(expected, resulting_strings);
        // }
        // {
        //     std::set<std::string> resulting_strings{};
        //     split_string(test_strings[0], resulting_strings, ',');      // split "a,a,a,aa"     with ','
        //     std::set<std::string> expected{"a", "a", "a", "aa"};
        //     EXPECT_EQ(expected, resulting_strings);
        // }
    }
}

// TEST(junction_detection, retrieve_aligned_segments)
// {

// }

// TEST(junction_detection, analyze_aligned_segments)
// {

// }

// TEST(junction_detection, analyze_cigar)
// {

// }