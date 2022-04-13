#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_parse_phylip.h"
#include "helper_methods.h"
#include "gubbins.h"

START_TEST (check_gubbins_no_recombinations)
{
	remove("../tests/data/no_recombinations.tre");
	cp("../tests/data/no_recombinations.tre", "../tests/data/no_recombinations.original.tre");
	run_gubbins("../tests/data/no_recombinations.aln.vcf", "../tests/data/no_recombinations.tre","../tests/data/no_recombinations.aln.snp_sites.aln",3,"../tests/data/no_recombinations.aln.snp_sites.aln",100,10000,0.05,1,0);
	ck_assert(file_exists("../tests/data/no_recombinations.tre.tab") == 1);
	ck_assert(file_exists("../tests/data/no_recombinations.tre.vcf") == 1);
	ck_assert(file_exists("../tests/data/no_recombinations.tre.phylip") == 1);
	ck_assert(file_exists("../tests/data/no_recombinations.tre.snp_sites.aln") == 1);
	ck_assert(file_exists("../tests/data/no_recombinations.tre.stats") == 1);
	ck_assert(file_exists("../tests/data/no_recombinations.tre.gff") == 1);
	
  ck_assert(number_of_recombinations_in_file("../tests/data/no_recombinations.tre.tab") == 0);
  ck_assert(compare_files("../tests/data/no_recombinations.tre","../tests/data/no_recombinations.expected.tre") == 1);
  ck_assert(compare_files("../tests/data/no_recombinations.tre.branch_snps.tab","../tests/data/no_recombinations.tre.branch_snps.expected.tab") == 1);

  remove("../tests/data/no_recombinations.tre");
	remove("../tests/data/no_recombinations.tre.tab");
	remove("../tests/data/no_recombinations.tre.vcf");
	remove("../tests/data/no_recombinations.tre.phylip");
	remove("../tests/data/no_recombinations.tre.stats");
	remove("../tests/data/no_recombinations.tre.gff");
	remove("../tests/data/no_recombinations.tre.snp_sites.aln");
	remove("../tests/data/no_recombinations.tre.branch_snps.tab");
}
END_TEST

START_TEST (check_gubbins_one_recombination)
{
	remove("../tests/data/one_recombination.tre");
	cp("../tests/data/one_recombination.tre", "../tests/data/one_recombination.original.tre");
	run_gubbins("../tests/data/one_recombination.aln.vcf", "../tests/data/one_recombination.tre","../tests/data/one_recombination.aln.snp_sites.aln",3,"../tests/data/one_recombination.aln.snp_sites.aln",100,10000,0.05,1,0);
	ck_assert(file_exists("../tests/data/one_recombination.tre.tab") == 1);
	ck_assert(file_exists("../tests/data/one_recombination.tre.vcf") == 1);
	ck_assert(file_exists("../tests/data/one_recombination.tre.phylip") == 1);
	ck_assert(file_exists("../tests/data/one_recombination.tre.stats") == 1);
	ck_assert(file_exists("../tests/data/one_recombination.tre.gff") == 1);
	ck_assert(file_exists("../tests/data/one_recombination.tre.snp_sites.aln") == 1);
	
	ck_assert(number_of_recombinations_in_file("../tests/data/one_recombination.tre.tab") == 1);
	ck_assert(compare_files("../tests/data/one_recombination.tre.vcf","../tests/data/one_recombination.expected.vcf") == 1);
	ck_assert(compare_files("../tests/data/one_recombination.tre.stats","../tests/data/one_recombination.expected.stats") == 1);
	ck_assert(compare_files("../tests/data/one_recombination.tre.branch_snps.tab","../tests/data/one_recombination.tre.branch_snps.expected.tab") == 1);
	ck_assert(compare_files("../tests/data/one_recombination.tre.tab","../tests/data/one_recombination.tre.expected.tab") == 1);
	ck_assert(compare_files("../tests/data/one_recombination.tre.gff","../tests/data/one_recombination.tre.expected.gff") == 1);
  
  remove("../tests/data/one_recombination.tre");
	remove("../tests/data/one_recombination.tre.tab");
	remove("../tests/data/one_recombination.tre.vcf");
	remove("../tests/data/one_recombination.tre.phylip");
	remove("../tests/data/one_recombination.tre.stats");
	remove("../tests/data/one_recombination.tre.gff");
	remove("../tests/data/one_recombination.tre.snp_sites.aln");
	remove("../tests/data/one_recombination.tre.branch_snps.tab");
}
END_TEST

START_TEST (check_gubbins_multiple_recombinations)
{
	remove("../tests/data/multiple_recombinations.tre");
//	cp("../tests/data/multiple_recombinations.tre", "../tests/data/multiple_recombinations.original.tre");
    cp("../tests/data/multiple_recombinations.original.tre", "../tests/data/multiple_recombinations.tre");

	run_gubbins("../tests/data/multiple_recombinations.aln.vcf",
                "../tests/data/multiple_recombinations.tre",
                "../tests/data/multiple_recombinations.jar.aln",
                3,
                "../tests/data/multiple_recombinations.aln",
                30,100,0.05,1,0);

	ck_assert(file_exists("../tests/data/multiple_recombinations.tre.tab") == 1);
	ck_assert(file_exists("../tests/data/multiple_recombinations.tre.vcf") == 1);
	ck_assert(file_exists("../tests/data/multiple_recombinations.tre.phylip") == 1);
	ck_assert(file_exists("../tests/data/multiple_recombinations.tre.stats") == 1);
	ck_assert(file_exists("../tests/data/multiple_recombinations.tre.gff") == 1);
	ck_assert(file_exists("../tests/data/multiple_recombinations.tre.snp_sites.aln") == 1);

	ck_assert(number_of_recombinations_in_file("../tests/data/multiple_recombinations.tre.tab") == 5);
    ck_assert(compare_files("../tests/data/multiple_recombinations.tre","../tests/data/multiple_recombinations.expected.tre") == 1);
    ck_assert(compare_files("../tests/data/multiple_recombinations.tre.branch_snps.tab","../tests/data/multiple_recombinations.tre.branch_snps.expected.tab") == 1);

	remove("../tests/data/multiple_recombinations.tre");
	remove("../tests/data/multiple_recombinations.tre.tab");
	remove("../tests/data/multiple_recombinations.tre.vcf");
	remove("../tests/data/multiple_recombinations.tre.phylip");
	remove("../tests/data/multiple_recombinations.tre.stats");
	remove("../tests/data/multiple_recombinations.tre.gff");
	remove("../tests/data/multiple_recombinations.tre.snp_sites.aln");
	remove("../tests/data/multiple_recombinations.tre.branch_snps.tab");
}
END_TEST

START_TEST (check_recombination_at_root)
{
	remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1");
	cp("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1", "../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.original.tre");

    
	run_gubbins("../tests/data/recombination_at_root/recombination_at_root.aln.gaps.vcf",
                "../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1",
                "../tests/data/recombination_at_root/recombination_at_root.aln.gaps.snp_sites.aln",
                3,
                "../tests/data/recombination_at_root/recombination_at_root.aln",
                100,10000,0.05,1,0);

    ck_assert(compare_files("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.tab","../tests/data/recombination_at_root/expected_RAxML_result.recombination_at_root.iteration_1.tab") == 1);

    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.vcf")             == 1);
    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.tab")             == 1);
    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.stats")           == 1);
    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.snp_sites.aln")   == 1);
    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.phylip")          == 1);
    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.gff")             == 1);
    ck_assert(file_exists("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.branch_snps.tab") == 1);
    
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.vcf");
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.tab");
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.stats");
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.snp_sites.aln");
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.phylip");
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.gff");
    remove("../tests/data/recombination_at_root/RAxML_result.recombination_at_root.iteration_1.branch_snps.tab");

}
END_TEST



Suite * run_gubbins_suite(void)
{
  Suite *s = suite_create ("Checking the gubbins functionality");
  TCase *tc_gubbins = tcase_create ("check_gubbins_produces_files");
  tcase_add_test (tc_gubbins, check_gubbins_no_recombinations);
  //tcase_add_test (tc_gubbins, check_gubbins_one_recombination);
  tcase_add_test (tc_gubbins, check_gubbins_multiple_recombinations);
  tcase_add_test (tc_gubbins, check_recombination_at_root);
  suite_add_tcase (s, tc_gubbins);
  return s;
}



