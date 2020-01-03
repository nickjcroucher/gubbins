# Change Log

## [Unreleased](https://github.com/sanger-pathogens/gubbins/tree/HEAD)

[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.3.4...HEAD)

**Implemented enhancements:**

- Retaining recombinant free non-polymorphic sites [\#208](https://github.com/sanger-pathogens/gubbins/issues/208)

**Fixed bugs:**

- Segmentation fault with some input files [\#170](https://github.com/sanger-pathogens/gubbins/issues/170)
- ltmain.sh should be removed? [\#157](https://github.com/sanger-pathogens/gubbins/issues/157)

**Closed issues:**

- Gubbins not generating final\_tree.tre [\#253](https://github.com/sanger-pathogens/gubbins/issues/253)
- Can gubbins take the fasta file produced by snp-sites as inputs? [\#250](https://github.com/sanger-pathogens/gubbins/issues/250)
- rho/theta--\> total recombination events? [\#247](https://github.com/sanger-pathogens/gubbins/issues/247)
- Gubbins output: internal\_5 / internal\_6 node? [\#246](https://github.com/sanger-pathogens/gubbins/issues/246)
- recombination\_predictions.gff regions are overlapping and duplicated ? [\#240](https://github.com/sanger-pathogens/gubbins/issues/240)
- Gubbins crashed, please ensure you have enough free memory [\#239](https://github.com/sanger-pathogens/gubbins/issues/239)
- GFF file output columns [\#234](https://github.com/sanger-pathogens/gubbins/issues/234)
- what does "Bases in Recombinations" in the stats file refer to? [\#232](https://github.com/sanger-pathogens/gubbins/issues/232)
- Make the FTP urls actual hyperlinks [\#229](https://github.com/sanger-pathogens/gubbins/issues/229)
- Failing test on 2.1.0 [\#177](https://github.com/sanger-pathogens/gubbins/issues/177)

**Merged pull requests:**

- Rename exec in docker [\#260](https://github.com/sanger-pathogens/gubbins/pull/260) ([puethe](https://github.com/puethe))
- Mark as unmaintained [\#259](https://github.com/sanger-pathogens/gubbins/pull/259) ([puethe](https://github.com/puethe))
- 627258 codecov [\#244](https://github.com/sanger-pathogens/gubbins/pull/244) ([ssjunnebo](https://github.com/ssjunnebo))
- Add changelog [\#242](https://github.com/sanger-pathogens/gubbins/pull/242) ([ssjunnebo](https://github.com/ssjunnebo))
- Rt574328 readme [\#241](https://github.com/sanger-pathogens/gubbins/pull/241) ([ssjunnebo](https://github.com/ssjunnebo))
- Remove INSTALL file [\#238](https://github.com/sanger-pathogens/gubbins/pull/238) ([ssjunnebo](https://github.com/ssjunnebo))
- 574328 edit readme [\#237](https://github.com/sanger-pathogens/gubbins/pull/237) ([ssjunnebo](https://github.com/ssjunnebo))
- update changelog [\#228](https://github.com/sanger-pathogens/gubbins/pull/228) ([puethe](https://github.com/puethe))

## [v2.3.4](https://github.com/sanger-pathogens/gubbins/tree/v2.3.4) (2018-07-30)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.3.3...v2.3.4)

**Merged pull requests:**

- dynamically allocate memory to avoid stack overflow [\#227](https://github.com/sanger-pathogens/gubbins/pull/227) ([puethe](https://github.com/puethe))

## [v2.3.3](https://github.com/sanger-pathogens/gubbins/tree/v2.3.3) (2018-07-26)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.3.2...v2.3.3)

**Closed issues:**

- ' INSTALL'  doesn't link to INSTALL.md [\#222](https://github.com/sanger-pathogens/gubbins/issues/222)

**Merged pull requests:**

- Dockerfile independent from debian package [\#226](https://github.com/sanger-pathogens/gubbins/pull/226) ([puethe](https://github.com/puethe))
- universal link [\#224](https://github.com/sanger-pathogens/gubbins/pull/224) ([puethe](https://github.com/puethe))
- link to install.md [\#223](https://github.com/sanger-pathogens/gubbins/pull/223) ([puethe](https://github.com/puethe))

## [v2.3.2](https://github.com/sanger-pathogens/gubbins/tree/v2.3.2) (2018-06-15)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.3.1...v2.3.2)

**Implemented enhancements:**

- Direct output to a directory [\#209](https://github.com/sanger-pathogens/gubbins/issues/209)

**Fixed bugs:**

- Warning while parsing tree: non-numeric label N38 for internal node  [\#154](https://github.com/sanger-pathogens/gubbins/issues/154)

**Closed issues:**

- Is there minimum sample size to use gubbins? [\#213](https://github.com/sanger-pathogens/gubbins/issues/213)
- gubbins\_drawer graph [\#212](https://github.com/sanger-pathogens/gubbins/issues/212)
- Running example files [\#211](https://github.com/sanger-pathogens/gubbins/issues/211)
- TypeError: reroot\_at\_midpoint\(\) got an unexpected keyword argument 'update\_splits' [\#205](https://github.com/sanger-pathogens/gubbins/issues/205)
- x [\#204](https://github.com/sanger-pathogens/gubbins/issues/204)
- Input format is not clear [\#185](https://github.com/sanger-pathogens/gubbins/issues/185)
- Inconsistencies in branch stats file [\#181](https://github.com/sanger-pathogens/gubbins/issues/181)
- error while loading shared libraries: libgubbins.so.0 [\#156](https://github.com/sanger-pathogens/gubbins/issues/156)

**Merged pull requests:**

- Replaced the script gubbins\_drawer.py with an empty script displaying… [\#221](https://github.com/sanger-pathogens/gubbins/pull/221) ([puethe](https://github.com/puethe))
- Link check to lsubunit [\#220](https://github.com/sanger-pathogens/gubbins/pull/220) ([puethe](https://github.com/puethe))
- corrected typo: 'fastree' -\> 'fasttree' [\#219](https://github.com/sanger-pathogens/gubbins/pull/219) ([puethe](https://github.com/puethe))
- update instructions for Docker [\#217](https://github.com/sanger-pathogens/gubbins/pull/217) ([ssjunnebo](https://github.com/ssjunnebo))
- added choices to `--tree\_builder` and other args [\#207](https://github.com/sanger-pathogens/gubbins/pull/207) ([schultzm](https://github.com/schultzm))
- ignore ltmain.sh for language statistics [\#206](https://github.com/sanger-pathogens/gubbins/pull/206) ([ssjunnebo](https://github.com/ssjunnebo))
- update manual [\#203](https://github.com/sanger-pathogens/gubbins/pull/203) ([andrewjpage](https://github.com/andrewjpage))

## [v2.3.1](https://github.com/sanger-pathogens/gubbins/tree/v2.3.1) (2017-09-22)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.3.0...v2.3.1)

**Fixed bugs:**

- gubbins\_drawer.py help improvement [\#183](https://github.com/sanger-pathogens/gubbins/issues/183)

**Merged pull requests:**

- Internal nodes [\#202](https://github.com/sanger-pathogens/gubbins/pull/202) ([andrewjpage](https://github.com/andrewjpage))

## [v2.3.0](https://github.com/sanger-pathogens/gubbins/tree/v2.3.0) (2017-09-21)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.2.3...v2.3.0)

**Fixed bugs:**

- Some options have default listed twice in --help [\#199](https://github.com/sanger-pathogens/gubbins/issues/199)

**Merged pull requests:**

- Gubbins drawer help text [\#201](https://github.com/sanger-pathogens/gubbins/pull/201) ([andrewjpage](https://github.com/andrewjpage))

## [v2.2.3](https://github.com/sanger-pathogens/gubbins/tree/v2.2.3) (2017-09-15)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.2.2...v2.2.3)

**Fixed bugs:**

- Error running gubbins\_drawer.py with v. 2.2.0 [\#184](https://github.com/sanger-pathogens/gubbins/issues/184)

**Closed issues:**

- Support raxml-AVX2 if available [\#196](https://github.com/sanger-pathogens/gubbins/issues/196)
- brew install gubbins errors - no formula for pillow [\#194](https://github.com/sanger-pathogens/gubbins/issues/194)
- gubbins would not run [\#187](https://github.com/sanger-pathogens/gubbins/issues/187)

**Merged pull requests:**

- Fix default extended help text. [\#200](https://github.com/sanger-pathogens/gubbins/pull/200) ([andrewjpage](https://github.com/andrewjpage))
- allow for AVX2 [\#197](https://github.com/sanger-pathogens/gubbins/pull/197) ([andrewjpage](https://github.com/andrewjpage))

## [v2.2.2](https://github.com/sanger-pathogens/gubbins/tree/v2.2.2) (2017-09-04)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.2.1...v2.2.2)

**Closed issues:**

- Forgot to make a 2.2.1 release? [\#195](https://github.com/sanger-pathogens/gubbins/issues/195)

## [v2.2.1](https://github.com/sanger-pathogens/gubbins/tree/v2.2.1) (2017-08-02)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.10...v2.2.1)

**Closed issues:**

- Gubbins - Ns in output alignment [\#192](https://github.com/sanger-pathogens/gubbins/issues/192)
- from gubbins import common                                     ImportError: cannot import name common [\#191](https://github.com/sanger-pathogens/gubbins/issues/191)
- Losing taxa? [\#189](https://github.com/sanger-pathogens/gubbins/issues/189)

**Merged pull requests:**

- updated manual [\#193](https://github.com/sanger-pathogens/gubbins/pull/193) ([ssjunnebo](https://github.com/ssjunnebo))
- Manual in docx format [\#188](https://github.com/sanger-pathogens/gubbins/pull/188) ([ssjunnebo](https://github.com/ssjunnebo))
- Fixed biopython sub\_features bug [\#186](https://github.com/sanger-pathogens/gubbins/pull/186) ([ssjunnebo](https://github.com/ssjunnebo))

## [v1.4.10](https://github.com/sanger-pathogens/gubbins/tree/v1.4.10) (2016-10-31)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.2.0...v1.4.10)

**Closed issues:**

- removing identical isolates [\#179](https://github.com/sanger-pathogens/gubbins/issues/179)

## [v2.2.0](https://github.com/sanger-pathogens/gubbins/tree/v2.2.0) (2016-10-31)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.1.0...v2.2.0)

**Closed issues:**

- Linuxbrew installs 1.4.7 still [\#176](https://github.com/sanger-pathogens/gubbins/issues/176)
- Can't find RAxML in path [\#174](https://github.com/sanger-pathogens/gubbins/issues/174)
- Source install doesn't install python wrappers [\#173](https://github.com/sanger-pathogens/gubbins/issues/173)

**Merged pull requests:**

- dont filter out identical sequences by default [\#180](https://github.com/sanger-pathogens/gubbins/pull/180) ([andrewjpage](https://github.com/andrewjpage))

## [v2.1.0](https://github.com/sanger-pathogens/gubbins/tree/v2.1.0) (2016-07-22)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v2.0.0...v2.1.0)

**Fixed bugs:**

- slow fastml [\#167](https://github.com/sanger-pathogens/gubbins/issues/167)

**Merged pull requests:**

- Use GTRCAT by default with RAxML [\#175](https://github.com/sanger-pathogens/gubbins/pull/175) ([andrewjpage](https://github.com/andrewjpage))

## [v2.0.0](https://github.com/sanger-pathogens/gubbins/tree/v2.0.0) (2016-07-15)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.9...v2.0.0)

**Closed issues:**

- .embl output file not opening in Artemis [\#165](https://github.com/sanger-pathogens/gubbins/issues/165)

**Merged pull requests:**

- Raxml reconstruction [\#172](https://github.com/sanger-pathogens/gubbins/pull/172) ([andrewjpage](https://github.com/andrewjpage))
- Raxml reconstruction [\#171](https://github.com/sanger-pathogens/gubbins/pull/171) ([andrewjpage](https://github.com/andrewjpage))

## [v1.4.9](https://github.com/sanger-pathogens/gubbins/tree/v1.4.9) (2016-04-15)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.8...v1.4.9)

**Merged pull requests:**

- Duplicate sequences [\#166](https://github.com/sanger-pathogens/gubbins/pull/166) ([andrewjpage](https://github.com/andrewjpage))

## [v1.4.8](https://github.com/sanger-pathogens/gubbins/tree/v1.4.8) (2016-04-13)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.7...v1.4.8)

**Closed issues:**

- numpy issues [\#159](https://github.com/sanger-pathogens/gubbins/issues/159)

**Merged pull requests:**

- Convert the proc output from bytes to strings [\#164](https://github.com/sanger-pathogens/gubbins/pull/164) ([aaronk](https://github.com/aaronk))

## [v1.4.7](https://github.com/sanger-pathogens/gubbins/tree/v1.4.7) (2016-02-29)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.6...v1.4.7)

**Merged pull requests:**

- midpoint rerooting by default [\#163](https://github.com/sanger-pathogens/gubbins/pull/163) ([andrewjpage](https://github.com/andrewjpage))
- Fix allocation from stack bug. [\#161](https://github.com/sanger-pathogens/gubbins/pull/161) ([jeromekelleher](https://github.com/jeromekelleher))

## [v1.4.6](https://github.com/sanger-pathogens/gubbins/tree/v1.4.6) (2016-02-29)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.5...v1.4.6)

## [v1.4.5](https://github.com/sanger-pathogens/gubbins/tree/v1.4.5) (2016-01-13)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v.1.4.5...v1.4.5)

## [v.1.4.5](https://github.com/sanger-pathogens/gubbins/tree/v.1.4.5) (2016-01-13)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.4...v.1.4.5)

**Merged pull requests:**

- remove install-sh [\#160](https://github.com/sanger-pathogens/gubbins/pull/160) ([andrewjpage](https://github.com/andrewjpage))
- Debian bug 807150 [\#158](https://github.com/sanger-pathogens/gubbins/pull/158) ([andrewjpage](https://github.com/andrewjpage))

## [v1.4.4](https://github.com/sanger-pathogens/gubbins/tree/v1.4.4) (2015-12-09)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.3...v1.4.4)

## [v1.4.3](https://github.com/sanger-pathogens/gubbins/tree/v1.4.3) (2015-11-19)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.2...v1.4.3)

**Closed issues:**

- multithreading [\#151](https://github.com/sanger-pathogens/gubbins/issues/151)

## [v1.4.2](https://github.com/sanger-pathogens/gubbins/tree/v1.4.2) (2015-09-02)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.1...v1.4.2)

**Merged pull requests:**

- check for input files which are too long [\#153](https://github.com/sanger-pathogens/gubbins/pull/153) ([andrewjpage](https://github.com/andrewjpage))
- update install docs [\#150](https://github.com/sanger-pathogens/gubbins/pull/150) ([andrewjpage](https://github.com/andrewjpage))

## [v1.4.1](https://github.com/sanger-pathogens/gubbins/tree/v1.4.1) (2015-07-07)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.4.0...v1.4.1)

## [v1.4.0](https://github.com/sanger-pathogens/gubbins/tree/v1.4.0) (2015-07-01)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.3.4...v1.4.0)

**Closed issues:**

- "Binary package - For 64 bit Linux" links to nowhere [\#142](https://github.com/sanger-pathogens/gubbins/issues/142)
- INSTALL file not obvious [\#129](https://github.com/sanger-pathogens/gubbins/issues/129)
- segmentation fault [\#49](https://github.com/sanger-pathogens/gubbins/issues/49)

**Merged pull requests:**

- Convert python to 3 [\#149](https://github.com/sanger-pathogens/gubbins/pull/149) ([andrewjpage](https://github.com/andrewjpage))
- update raxml version for travis [\#148](https://github.com/sanger-pathogens/gubbins/pull/148) ([andrewjpage](https://github.com/andrewjpage))
- Add build status to README [\#147](https://github.com/sanger-pathogens/gubbins/pull/147) ([bewt85](https://github.com/bewt85))
- Add TravisCI support [\#146](https://github.com/sanger-pathogens/gubbins/pull/146) ([bewt85](https://github.com/bewt85))
- Reference URL updated to final published version not early access [\#145](https://github.com/sanger-pathogens/gubbins/pull/145) ([aslett1](https://github.com/aslett1))
- Manually merge Aidan Delaneys pull request [\#144](https://github.com/sanger-pathogens/gubbins/pull/144) ([andrewjpage](https://github.com/andrewjpage))
- Detect blocks at end of sequence [\#143](https://github.com/sanger-pathogens/gubbins/pull/143) ([andrewjpage](https://github.com/andrewjpage))

## [v1.3.4](https://github.com/sanger-pathogens/gubbins/tree/v1.3.4) (2015-05-18)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.3.3...v1.3.4)

**Closed issues:**

- No ./configure in github - need to run autoconf? [\#130](https://github.com/sanger-pathogens/gubbins/issues/130)

## [v1.3.3](https://github.com/sanger-pathogens/gubbins/tree/v1.3.3) (2015-04-16)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.3.1...v1.3.3)

**Merged pull requests:**

- Window size self error [\#140](https://github.com/sanger-pathogens/gubbins/pull/140) ([andrewjpage](https://github.com/andrewjpage))
- Window size self error [\#139](https://github.com/sanger-pathogens/gubbins/pull/139) ([andrewjpage](https://github.com/andrewjpage))

## [v1.3.1](https://github.com/sanger-pathogens/gubbins/tree/v1.3.1) (2015-04-16)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.3.0...v1.3.1)

**Merged pull requests:**

- version 1.3.1 [\#138](https://github.com/sanger-pathogens/gubbins/pull/138) ([andrewjpage](https://github.com/andrewjpage))

## [v1.3.0](https://github.com/sanger-pathogens/gubbins/tree/v1.3.0) (2015-04-16)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.2.4...v1.3.0)

**Merged pull requests:**

- Recalculate genome length [\#137](https://github.com/sanger-pathogens/gubbins/pull/137) ([andrewjpage](https://github.com/andrewjpage))
- Downstream recombinations [\#136](https://github.com/sanger-pathogens/gubbins/pull/136) ([andrewjpage](https://github.com/andrewjpage))

## [v1.2.4](https://github.com/sanger-pathogens/gubbins/tree/v1.2.4) (2015-03-31)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.2.3...v1.2.4)

**Merged pull requests:**

- output version in python script [\#135](https://github.com/sanger-pathogens/gubbins/pull/135) ([andrewjpage](https://github.com/andrewjpage))

## [v1.2.3](https://github.com/sanger-pathogens/gubbins/tree/v1.2.3) (2015-03-30)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.2.2...v1.2.3)

**Merged pull requests:**

- include last element of array [\#134](https://github.com/sanger-pathogens/gubbins/pull/134) ([andrewjpage](https://github.com/andrewjpage))

## [v1.2.2](https://github.com/sanger-pathogens/gubbins/tree/v1.2.2) (2015-03-30)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.2.1...v1.2.2)

## [v1.2.1](https://github.com/sanger-pathogens/gubbins/tree/v1.2.1) (2015-03-30)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.2.0...v1.2.1)

**Merged pull requests:**

- Correct vcf format and EMBL format [\#133](https://github.com/sanger-pathogens/gubbins/pull/133) ([andrewjpage](https://github.com/andrewjpage))
- debian changelog [\#132](https://github.com/sanger-pathogens/gubbins/pull/132) ([andrewjpage](https://github.com/andrewjpage))
- Dont fill in gaps in parent with bases from child [\#131](https://github.com/sanger-pathogens/gubbins/pull/131) ([andrewjpage](https://github.com/andrewjpage))

## [v1.2.0](https://github.com/sanger-pathogens/gubbins/tree/v1.2.0) (2015-03-06)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.1.2...v1.2.0)

## [v1.1.2](https://github.com/sanger-pathogens/gubbins/tree/v1.1.2) (2015-03-06)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.1.1...v1.1.2)

**Merged pull requests:**

- Update robustness [\#128](https://github.com/sanger-pathogens/gubbins/pull/128) ([andrewjpage](https://github.com/andrewjpage))
- message about invalid outgroups [\#127](https://github.com/sanger-pathogens/gubbins/pull/127) ([andrewjpage](https://github.com/andrewjpage))

## [v1.1.1](https://github.com/sanger-pathogens/gubbins/tree/v1.1.1) (2015-01-23)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.1.0...v1.1.1)

**Merged pull requests:**

- Clade as outgroup [\#126](https://github.com/sanger-pathogens/gubbins/pull/126) ([andrewjpage](https://github.com/andrewjpage))

## [v1.1.0](https://github.com/sanger-pathogens/gubbins/tree/v1.1.0) (2015-01-23)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.12...v1.1.0)

## [v1.0.12](https://github.com/sanger-pathogens/gubbins/tree/v1.0.12) (2015-01-20)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.11...v1.0.12)

**Merged pull requests:**

- test data [\#125](https://github.com/sanger-pathogens/gubbins/pull/125) ([andrewjpage](https://github.com/andrewjpage))
- Remove internal nodes from final tree [\#124](https://github.com/sanger-pathogens/gubbins/pull/124) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.11](https://github.com/sanger-pathogens/gubbins/tree/v1.0.11) (2015-01-19)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.10...v1.0.11)

**Merged pull requests:**

- check fasta files valid and check for intermediate raxml files [\#123](https://github.com/sanger-pathogens/gubbins/pull/123) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.10](https://github.com/sanger-pathogens/gubbins/tree/v1.0.10) (2015-01-19)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.9...v1.0.10)

## [v1.0.9](https://github.com/sanger-pathogens/gubbins/tree/v1.0.9) (2015-01-06)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/1.0.9...v1.0.9)

## [1.0.9](https://github.com/sanger-pathogens/gubbins/tree/1.0.9) (2015-01-06)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.8...1.0.9)

**Merged pull requests:**

- Dont use cpu info to find the flags [\#122](https://github.com/sanger-pathogens/gubbins/pull/122) ([andrewjpage](https://github.com/andrewjpage))
- class namespace [\#121](https://github.com/sanger-pathogens/gubbins/pull/121) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.8](https://github.com/sanger-pathogens/gubbins/tree/v1.0.8) (2014-12-08)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.7...v1.0.8)

**Merged pull requests:**

- reduced functionality but more stable cpu info lookup [\#120](https://github.com/sanger-pathogens/gubbins/pull/120) ([andrewjpage](https://github.com/andrewjpage))
- Rho theta [\#119](https://github.com/sanger-pathogens/gubbins/pull/119) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.7](https://github.com/sanger-pathogens/gubbins/tree/v1.0.7) (2014-12-08)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.6...v1.0.7)

**Merged pull requests:**

- allow for variation and exclude zero track length [\#118](https://github.com/sanger-pathogens/gubbins/pull/118) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.6](https://github.com/sanger-pathogens/gubbins/tree/v1.0.6) (2014-12-03)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.5...v1.0.6)

## [v1.0.5](https://github.com/sanger-pathogens/gubbins/tree/v1.0.5) (2014-11-19)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.2...v1.0.5)

**Merged pull requests:**

- fasttree executables with titlecase and lowercase [\#117](https://github.com/sanger-pathogens/gubbins/pull/117) ([andrewjpage](https://github.com/andrewjpage))
- warning message when changing number of threads [\#116](https://github.com/sanger-pathogens/gubbins/pull/116) ([andrewjpage](https://github.com/andrewjpage))
- Change raxml exec search order [\#115](https://github.com/sanger-pathogens/gubbins/pull/115) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.2](https://github.com/sanger-pathogens/gubbins/tree/v1.0.2) (2014-10-16)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0.1...v1.0.2)

**Merged pull requests:**

- rename convergence methods and allow for more [\#114](https://github.com/sanger-pathogens/gubbins/pull/114) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0.1](https://github.com/sanger-pathogens/gubbins/tree/v1.0.1) (2014-10-14)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v1.0...v1.0.1)

**Merged pull requests:**

- unstable test fixed [\#113](https://github.com/sanger-pathogens/gubbins/pull/113) ([andrewjpage](https://github.com/andrewjpage))

## [v1.0](https://github.com/sanger-pathogens/gubbins/tree/v1.0) (2014-10-14)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.6...v1.0)

**Merged pull requests:**

- Converging recombinations [\#112](https://github.com/sanger-pathogens/gubbins/pull/112) ([andrewjpage](https://github.com/andrewjpage))

## [0.6](https://github.com/sanger-pathogens/gubbins/tree/0.6) (2014-09-24)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.7...0.6)

## [0.7](https://github.com/sanger-pathogens/gubbins/tree/0.7) (2014-09-24)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v0.6...0.7)

## [v0.6](https://github.com/sanger-pathogens/gubbins/tree/v0.6) (2014-09-09)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v0.5...v0.6)

**Merged pull requests:**

- Gubbins drawer can use recombination embl file [\#111](https://github.com/sanger-pathogens/gubbins/pull/111) ([andrewjpage](https://github.com/andrewjpage))

## [v0.5](https://github.com/sanger-pathogens/gubbins/tree/v0.5) (2014-09-09)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v0.4.1...v0.5)

**Merged pull requests:**

- thread support for raxml [\#110](https://github.com/sanger-pathogens/gubbins/pull/110) ([andrewjpage](https://github.com/andrewjpage))

## [v0.4.1](https://github.com/sanger-pathogens/gubbins/tree/v0.4.1) (2014-09-08)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v0.4...v0.4.1)

**Merged pull requests:**

- Check for lowercase fasttree exec [\#109](https://github.com/sanger-pathogens/gubbins/pull/109) ([andrewjpage](https://github.com/andrewjpage))

## [v0.4](https://github.com/sanger-pathogens/gubbins/tree/v0.4) (2014-09-05)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/v0.3...v0.4)

**Merged pull requests:**

- pairwise rename output files [\#108](https://github.com/sanger-pathogens/gubbins/pull/108) ([andrewjpage](https://github.com/andrewjpage))
- Update stats headers [\#107](https://github.com/sanger-pathogens/gubbins/pull/107) ([andrewjpage](https://github.com/andrewjpage))

## [v0.3](https://github.com/sanger-pathogens/gubbins/tree/v0.3) (2014-09-04)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.2...v0.3)

## [0.2](https://github.com/sanger-pathogens/gubbins/tree/0.2) (2014-09-02)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.1.11...0.2)

**Merged pull requests:**

- Catch too few sequences [\#106](https://github.com/sanger-pathogens/gubbins/pull/106) ([andrewjpage](https://github.com/andrewjpage))

## [0.1.11](https://github.com/sanger-pathogens/gubbins/tree/0.1.11) (2014-08-22)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.1.10...0.1.11)

**Merged pull requests:**

- version in configure [\#105](https://github.com/sanger-pathogens/gubbins/pull/105) ([andrewjpage](https://github.com/andrewjpage))
- Support for bipartition trees as input [\#104](https://github.com/sanger-pathogens/gubbins/pull/104) ([andrewjpage](https://github.com/andrewjpage))

## [0.1.10](https://github.com/sanger-pathogens/gubbins/tree/0.1.10) (2014-08-22)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.1.9...0.1.10)

**Merged pull requests:**

- print list of excluded sequences [\#103](https://github.com/sanger-pathogens/gubbins/pull/103) ([andrewjpage](https://github.com/andrewjpage))

## [0.1.9](https://github.com/sanger-pathogens/gubbins/tree/0.1.9) (2014-08-22)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.1.8...0.1.9)

**Merged pull requests:**

- print errors for invalid input files [\#102](https://github.com/sanger-pathogens/gubbins/pull/102) ([andrewjpage](https://github.com/andrewjpage))

## [0.1.8](https://github.com/sanger-pathogens/gubbins/tree/0.1.8) (2014-08-15)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.1.7...0.1.8)

**Merged pull requests:**

- Fasttree exec [\#101](https://github.com/sanger-pathogens/gubbins/pull/101) ([andrewjpage](https://github.com/andrewjpage))

## [0.1.7](https://github.com/sanger-pathogens/gubbins/tree/0.1.7) (2014-08-14)
[Full Changelog](https://github.com/sanger-pathogens/gubbins/compare/0.1.6...0.1.7)

**Merged pull requests:**

- Build for trusty [\#100](https://github.com/sanger-pathogens/gubbins/pull/100) ([andrewjpage](https://github.com/andrewjpage))

## [0.1.6](https://github.com/sanger-pathogens/gubbins/tree/0.1.6) (2014-08-12)
**Merged pull requests:**

- Changes from Aidan Delaney [\#99](https://github.com/sanger-pathogens/gubbins/pull/99) ([andrewjpage](https://github.com/andrewjpage))
- Stats missing snps [\#98](https://github.com/sanger-pathogens/gubbins/pull/98) ([andrewjpage](https://github.com/andrewjpage))
- update no support option in fastml [\#97](https://github.com/sanger-pathogens/gubbins/pull/97) ([andrewjpage](https://github.com/andrewjpage))
- move blocks inwards [\#96](https://github.com/sanger-pathogens/gubbins/pull/96) ([andrewjpage](https://github.com/andrewjpage))
- define is starting tree method as static [\#95](https://github.com/sanger-pathogens/gubbins/pull/95) ([andrewjpage](https://github.com/andrewjpage))
- starting tree needed to come from args [\#94](https://github.com/sanger-pathogens/gubbins/pull/94) ([andrewjpage](https://github.com/andrewjpage))
- Python tests [\#93](https://github.com/sanger-pathogens/gubbins/pull/93) ([andrewjpage](https://github.com/andrewjpage))
- remove duplicate move [\#92](https://github.com/sanger-pathogens/gubbins/pull/92) ([andrewjpage](https://github.com/andrewjpage))
- dont move in twice [\#91](https://github.com/sanger-pathogens/gubbins/pull/91) ([andrewjpage](https://github.com/andrewjpage))
- typo in method, too many args [\#90](https://github.com/sanger-pathogens/gubbins/pull/90) ([andrewjpage](https://github.com/andrewjpage))
- move end block in correct number of snps [\#89](https://github.com/sanger-pathogens/gubbins/pull/89) ([andrewjpage](https://github.com/andrewjpage))
- Span gaps [\#88](https://github.com/sanger-pathogens/gubbins/pull/88) ([andrewjpage](https://github.com/andrewjpage))
- translate snps and gaps to original genome space coords [\#87](https://github.com/sanger-pathogens/gubbins/pull/87) ([andrewjpage](https://github.com/andrewjpage))
- extend blocks over gaps [\#86](https://github.com/sanger-pathogens/gubbins/pull/86) ([andrewjpage](https://github.com/andrewjpage))
- dont call methods and variables the same thing [\#85](https://github.com/sanger-pathogens/gubbins/pull/85) ([andrewjpage](https://github.com/andrewjpage))
- remove latest tree symlink if it exists, even if its broken [\#84](https://github.com/sanger-pathogens/gubbins/pull/84) ([andrewjpage](https://github.com/andrewjpage))
- delete intermediate raxml files [\#83](https://github.com/sanger-pathogens/gubbins/pull/83) ([andrewjpage](https://github.com/andrewjpage))
- Eliminate memory leaks and uninitalised reads/writes [\#82](https://github.com/sanger-pathogens/gubbins/pull/82) ([andrewjpage](https://github.com/andrewjpage))
- GPL [\#81](https://github.com/sanger-pathogens/gubbins/pull/81) ([andrewjpage](https://github.com/andrewjpage))
- reduce memory usage and fix seg fault [\#80](https://github.com/sanger-pathogens/gubbins/pull/80) ([andrewjpage](https://github.com/andrewjpage))
- move math lib linking to end [\#79](https://github.com/sanger-pathogens/gubbins/pull/79) ([andrewjpage](https://github.com/andrewjpage))
- move zlib linking to end [\#78](https://github.com/sanger-pathogens/gubbins/pull/78) ([andrewjpage](https://github.com/andrewjpage))
- swap output suppression [\#77](https://github.com/sanger-pathogens/gubbins/pull/77) ([andrewjpage](https://github.com/andrewjpage))
- memory leaks [\#76](https://github.com/sanger-pathogens/gubbins/pull/76) ([andrewjpage](https://github.com/andrewjpage))
- update window creation [\#75](https://github.com/sanger-pathogens/gubbins/pull/75) ([andrewjpage](https://github.com/andrewjpage))
- part4 log10 [\#74](https://github.com/sanger-pathogens/gubbins/pull/74) ([andrewjpage](https://github.com/andrewjpage))
- Improving likelihood [\#73](https://github.com/sanger-pathogens/gubbins/pull/73) ([andrewjpage](https://github.com/andrewjpage))
- Improving likelihood [\#72](https://github.com/sanger-pathogens/gubbins/pull/72) ([andrewjpage](https://github.com/andrewjpage))
- reroot on outgroup should keep internal nodes [\#71](https://github.com/sanger-pathogens/gubbins/pull/71) ([andrewjpage](https://github.com/andrewjpage))
- Use bases for gaps for internal nodes after fastml [\#70](https://github.com/sanger-pathogens/gubbins/pull/70) ([andrewjpage](https://github.com/andrewjpage))
- propagate unambiguous bases up the tree [\#69](https://github.com/sanger-pathogens/gubbins/pull/69) ([andrewjpage](https://github.com/andrewjpage))
- 	 move inwards one at a time, remove gap extending, Bonferroni correction, sliding window 1 [\#68](https://github.com/sanger-pathogens/gubbins/pull/68) ([andrewjpage](https://github.com/andrewjpage))
- new window shifting method and bug with detecting snps [\#67](https://github.com/sanger-pathogens/gubbins/pull/67) ([andrewjpage](https://github.com/andrewjpage))
- dont move window past size of genome [\#66](https://github.com/sanger-pathogens/gubbins/pull/66) ([andrewjpage](https://github.com/andrewjpage))
- move window in to first snp [\#65](https://github.com/sanger-pathogens/gubbins/pull/65) ([andrewjpage](https://github.com/andrewjpage))
- Updates to coords other minor bugs [\#64](https://github.com/sanger-pathogens/gubbins/pull/64) ([andrewjpage](https://github.com/andrewjpage))
- sliding window at half minimum window size [\#63](https://github.com/sanger-pathogens/gubbins/pull/63) ([andrewjpage](https://github.com/andrewjpage))
- swap bases in snp tab file and cleanup filenames [\#62](https://github.com/sanger-pathogens/gubbins/pull/62) ([andrewjpage](https://github.com/andrewjpage))
- Gaps can be in ref as well as alt [\#61](https://github.com/sanger-pathogens/gubbins/pull/61) ([andrewjpage](https://github.com/andrewjpage))
- pairwise [\#60](https://github.com/sanger-pathogens/gubbins/pull/60) ([andrewjpage](https://github.com/andrewjpage))
- remove filtered sequences from starting tree [\#59](https://github.com/sanger-pathogens/gubbins/pull/59) ([andrewjpage](https://github.com/andrewjpage))
- missing comma [\#58](https://github.com/sanger-pathogens/gubbins/pull/58) ([andrewjpage](https://github.com/andrewjpage))
- unquote tree nodes for fastml and check if too much is filtered out [\#57](https://github.com/sanger-pathogens/gubbins/pull/57) ([andrewjpage](https://github.com/andrewjpage))
- update test files [\#56](https://github.com/sanger-pathogens/gubbins/pull/56) ([andrewjpage](https://github.com/andrewjpage))
- Fastml no gaps [\#55](https://github.com/sanger-pathogens/gubbins/pull/55) ([andrewjpage](https://github.com/andrewjpage))
- readme [\#54](https://github.com/sanger-pathogens/gubbins/pull/54) ([andrewjpage](https://github.com/andrewjpage))
- snp branch output and prefiltering aln [\#53](https://github.com/sanger-pathogens/gubbins/pull/53) ([andrewjpage](https://github.com/andrewjpage))
- make drawer work with empty tracks [\#52](https://github.com/sanger-pathogens/gubbins/pull/52) ([andrewjpage](https://github.com/andrewjpage))
- strip off input params when testing if exec is in path [\#51](https://github.com/sanger-pathogens/gubbins/pull/51) ([andrewjpage](https://github.com/andrewjpage))
- Fastml integration [\#50](https://github.com/sanger-pathogens/gubbins/pull/50) ([andrewjpage](https://github.com/andrewjpage))
- versioning [\#48](https://github.com/sanger-pathogens/gubbins/pull/48) ([andrewjpage](https://github.com/andrewjpage))
- Cleanup intermediate files [\#47](https://github.com/sanger-pathogens/gubbins/pull/47) ([andrewjpage](https://github.com/andrewjpage))
- treat N in the same way as gaps when snp searching [\#46](https://github.com/sanger-pathogens/gubbins/pull/46) ([andrewjpage](https://github.com/andrewjpage))
- Treat N as gap when extending over gaps at block boundries [\#45](https://github.com/sanger-pathogens/gubbins/pull/45) ([andrewjpage](https://github.com/andrewjpage))
- use same fasttree params in hybrid [\#44](https://github.com/sanger-pathogens/gubbins/pull/44) ([andrewjpage](https://github.com/andrewjpage))
- vcf coords start from 1 [\#43](https://github.com/sanger-pathogens/gubbins/pull/43) ([andrewjpage](https://github.com/andrewjpage))
- unquote strings [\#42](https://github.com/sanger-pathogens/gubbins/pull/42) ([andrewjpage](https://github.com/andrewjpage))
- default preserve underscores to true [\#41](https://github.com/sanger-pathogens/gubbins/pull/41) ([andrewjpage](https://github.com/andrewjpage))
- preserve underscores [\#40](https://github.com/sanger-pathogens/gubbins/pull/40) ([andrewjpage](https://github.com/andrewjpage))
- Reroot around edge [\#39](https://github.com/sanger-pathogens/gubbins/pull/39) ([andrewjpage](https://github.com/andrewjpage))
- remove calls to reading beginnign of line [\#38](https://github.com/sanger-pathogens/gubbins/pull/38) ([andrewjpage](https://github.com/andrewjpage))
- Read line memory leak [\#37](https://github.com/sanger-pathogens/gubbins/pull/37) ([andrewjpage](https://github.com/andrewjpage))
- increase read line memory [\#36](https://github.com/sanger-pathogens/gubbins/pull/36) ([andrewjpage](https://github.com/andrewjpage))
- resize read line buffer [\#35](https://github.com/sanger-pathogens/gubbins/pull/35) ([andrewjpage](https://github.com/andrewjpage))
- Copyfile [\#34](https://github.com/sanger-pathogens/gubbins/pull/34) ([andrewjpage](https://github.com/andrewjpage))
- increase read buffer size [\#33](https://github.com/sanger-pathogens/gubbins/pull/33) ([andrewjpage](https://github.com/andrewjpage))
- semi colon at end of tree [\#32](https://github.com/sanger-pathogens/gubbins/pull/32) ([andrewjpage](https://github.com/andrewjpage))
- speedup gubbins [\#31](https://github.com/sanger-pathogens/gubbins/pull/31) ([andrewjpage](https://github.com/andrewjpage))
- scale distances by num snps [\#30](https://github.com/sanger-pathogens/gubbins/pull/30) ([andrewjpage](https://github.com/andrewjpage))
- providing a starting tree skips first tree building step [\#29](https://github.com/sanger-pathogens/gubbins/pull/29) ([andrewjpage](https://github.com/andrewjpage))
- calc bases in recombs [\#28](https://github.com/sanger-pathogens/gubbins/pull/28) ([andrewjpage](https://github.com/andrewjpage))
- improved stats [\#27](https://github.com/sanger-pathogens/gubbins/pull/27) ([andrewjpage](https://github.com/andrewjpage))
- calc snp stats and parameter to set min snps [\#26](https://github.com/sanger-pathogens/gubbins/pull/26) ([andrewjpage](https://github.com/andrewjpage))
- Get all gaps [\#25](https://github.com/sanger-pathogens/gubbins/pull/25) ([andrewjpage](https://github.com/andrewjpage))
- increase window size to 20k [\#24](https://github.com/sanger-pathogens/gubbins/pull/24) ([andrewjpage](https://github.com/andrewjpage))
- Span gaps [\#23](https://github.com/sanger-pathogens/gubbins/pull/23) ([andrewjpage](https://github.com/andrewjpage))
- fasttree needs to use output of gubbins [\#22](https://github.com/sanger-pathogens/gubbins/pull/22) ([andrewjpage](https://github.com/andrewjpage))
- extend blocks over gaps [\#21](https://github.com/sanger-pathogens/gubbins/pull/21) ([andrewjpage](https://github.com/andrewjpage))
- raxml exec check [\#20](https://github.com/sanger-pathogens/gubbins/pull/20) ([andrewjpage](https://github.com/andrewjpage))
- Merge blocks straddling gaps [\#19](https://github.com/sanger-pathogens/gubbins/pull/19) ([andrewjpage](https://github.com/andrewjpage))
- fix remove [\#18](https://github.com/sanger-pathogens/gubbins/pull/18) ([andrewjpage](https://github.com/andrewjpage))
- speedup phylip file writing second try [\#17](https://github.com/sanger-pathogens/gubbins/pull/17) ([andrewjpage](https://github.com/andrewjpage))
- speedup phylip file writing [\#16](https://github.com/sanger-pathogens/gubbins/pull/16) ([andrewjpage](https://github.com/andrewjpage))
- provide alignment filename not full path for fastree  [\#15](https://github.com/sanger-pathogens/gubbins/pull/15) ([andrewjpage](https://github.com/andrewjpage))
- remove hardcoded path to exec [\#14](https://github.com/sanger-pathogens/gubbins/pull/14) ([andrewjpage](https://github.com/andrewjpage))
- rerun autoreconf [\#13](https://github.com/sanger-pathogens/gubbins/pull/13) ([andrewjpage](https://github.com/andrewjpage))
- midpoint rerooting and pairwise comparison [\#12](https://github.com/sanger-pathogens/gubbins/pull/12) ([andrewjpage](https://github.com/andrewjpage))
- Use get opt long [\#11](https://github.com/sanger-pathogens/gubbins/pull/11) ([andrewjpage](https://github.com/andrewjpage))
- dont point to gubbins in same dir [\#10](https://github.com/sanger-pathogens/gubbins/pull/10) ([andrewjpage](https://github.com/andrewjpage))
- dont print out rerooted tree to stdout [\#9](https://github.com/sanger-pathogens/gubbins/pull/9) ([andrewjpage](https://github.com/andrewjpage))
- Hybrid tree building [\#8](https://github.com/sanger-pathogens/gubbins/pull/8) ([andrewjpage](https://github.com/andrewjpage))
- Rewrite perl as python [\#7](https://github.com/sanger-pathogens/gubbins/pull/7) ([andrewjpage](https://github.com/andrewjpage))
- send in most recent phylip file to raxml [\#6](https://github.com/sanger-pathogens/gubbins/pull/6) ([andrewjpage](https://github.com/andrewjpage))
- run gubbins original snp sites [\#5](https://github.com/sanger-pathogens/gubbins/pull/5) ([andrewjpage](https://github.com/andrewjpage))
- run gubbins with original snp sites [\#4](https://github.com/sanger-pathogens/gubbins/pull/4) ([andrewjpage](https://github.com/andrewjpage))
- Gubbins overall tests [\#3](https://github.com/sanger-pathogens/gubbins/pull/3) ([andrewjpage](https://github.com/andrewjpage))
- Fix memory leak and properly remove recombinations [\#2](https://github.com/sanger-pathogens/gubbins/pull/2) ([andrewjpage](https://github.com/andrewjpage))
- dont ignore internal nodes [\#1](https://github.com/sanger-pathogens/gubbins/pull/1) ([andrewjpage](https://github.com/andrewjpage))



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*