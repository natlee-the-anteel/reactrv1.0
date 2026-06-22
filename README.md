[![DOI](https://zenodo.org/badge/1129249052.svg)](https://doi.org/10.5281/zenodo.18306541)
DOWNLOAD AND INSTALL INSTRUCTIONS:
---------------------------------------------
please contact natlee@nuevaschool.org for a copy of current citations, if help is needed, or want to reach out.
## 1. System Preparation
Before installing the software, you must prepare the MacOS environment to handle developer tools and older bioinformatics binaries. you must also download java. 
        
    xcode-select --install
    softwareupdate --install-rosetta --agree-to-license

## 2. Installation of Package Manager (Miniforge)
We use Miniforge to manage most of the software dependencies automatically.

    curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh -o miniforge.sh
    bash miniforge.sh

## 3. Installation of reactr
run the following
  
    git clone https://github.com/natlee-the-anteel/reactrv1.0.git
    cd reactrv1.0
    export PATH=/path/to/reactrv1.0/bin/:$PATH

## 3. Startup the environment and activate it
run the following

    conda env create --file environment.yaml
    conda activate reactr

## 5. Download the non-conda dependencies
1. for Flashfry
    
        mkdir -p preset
        cd preset
        wget https://github.com/mckennalab/FlashFry/releases/download/1.15/FlashFry-assembly-1.15.jar
        mv FlashFry-assembly-1.15.jar FlashFry.jar

2. for MCScanX

        wget https://github.com/wyp1125/MCScanX/archive/refs/heads/master.zip -O MCScanX.zip
        unzip MCScanX.zip
        cd MCScanX-master
        make
        mv MCScanX ../
        mv duplicate_gene_classifier ../
        cd ..
        rm -rf MCScanX-master MCScanX.zip

3. for pfam database

        mkdir -p pfam
        cd pfam
        wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        gunzip Pfam-A.hmm.gz
        cd ../.. 
        hmmpress -f preset/pfam/Pfam-A.hmm
4. install deeploc2 from DTU (requires academic liscence, use the readme, can be long to install, so if needed, remove rule from main command)
-----------------------------------------------------
Usage Instructions

steps for users
if you want to load new genomes: (1) edit the taxonomy IDs in the config.yaml (2) delete the folder "reactr/data" if applicable, (3) run 
    
    snakemake -s LoadDatasets.smk --cores 8 --rerun-incomplete --forceall -p

if you are content with the current genomes or don't have one yet: (1) edit in desired protein fasta(s) from the base genome, in the top of the config.yaml, (2) check if it matches the base genome and ideally, the same sequenced version, (3) delete the folder "reactr/output" if applicable, (4) run

    snakemake -s MainPipeline.smk --cores all --rerun-incomplete --forceall -p
    
wait. this should take a few minutes max, though it scales with the number of proteins you query. next, once it says it's complete, then run 

    snakemake -s MainPipeline.smk --cores all --rerun-incomplete --forceall -p
------------------------------------------------------
Format Instructions
The only things that you really should need to edit (unless you're directly manipulating to code), is just the config.yaml. Specifically, just the taxonids (they're ncbi ids, they autodownload all the necesary stuff if you simply run the loaddatasets.smk rule (see above), and the query). Below is an example format: 

        query_contents: |
          >sp|Q9LQT8.1|GAI_ARATH RecName: Full=DELLA protein GAI; AltName: Full=GRAS family protein 3; Short=AtGRAS-3; AltName: Full=Gibberellic acid-insensitive mutant protein; AltName: Full=Restoration of growth on ammonia protein 2
          MKRDHHHHHHQDKKTMMMNEEDDGNGMDELLAVLGYKVRSSEMADVAQKLEQLEVMMSNVQEDDLSQLAT
          ETVHYNPAELYTWLDSMLTDLNPPSSNAEYDLKAIPGDAILNQFAIDSASSSNQGGGGDTYTTNKRLKCS
          NGVVETTTATAESTRHVVLVDSQENGVRLVHALLACAEAVQKENLTVAEALVKQIGFLAVSQIGAMRKVA
          TYFAEALARRIYRLSPSQSPIDHSLSDTLQMHFYETCPYLKFAHFTANQAILEAFQGKKRVHVIDFSMSQ
          GLQWPALMQALALRPGGPPVFRLTGIGPPAPDNFDYLHEVGCKLAHLAEAIHVEFEYRGFVANTLADLDA
          SMLELRPSEIESVAVNSVFELHKLLGRPGAIDKVLGVVNQIKPEIFTVVEQESNHNSPIFLDRFTESLHY
          YSTLFDSLEGVPSGQDKVMSEVYLGKQICNVVACDGPDRVERHETLSQWRNRFGSAGFAAAHIGSNAFKQ
          ASMLLALFNGGEGYRVEESDGCLMLGWHTRPLIATSAWKLSTN
          >sp|Q9SLH3.1|RGA_ARATH RecName: Full=DELLA protein RGA; AltName: Full=GAI-related sequence; AltName: Full=GRAS family protein 10; Short=AtGRAS-10; AltName: Full=Repressor on the ga1-3 mutant; AltName: Full=Restoration of growth on ammonia protein 1
          MKRDHHQFQGRLSNHGTSSSSSSISKDKMMMVKKEEDGGGNMDDELLAVLGYKVRSSEMAEVALKLEQLE
          TMMSNVQEDGLSHLATDTVHYNPSELYSWLDNMLSELNPPPLPASSNGLDPVLPSPEICGFPASDYDLKV
          IPGNAIYQFPAIDSSSSSNNQNKRLKSCSSPDSMVTSTSTGTQIGGVIGTTVTTTTTTTTAAGESTRSVI
          LVDSQENGVRLVHALMACAEAIQQNNLTLAEALVKQIGCLAVSQAGAMRKVATYFAEALARRIYRLSPPQ
          NQIDHCLSDTLQMHFYETCPYLKFAHFTANQAILEAFEGKKRVHVIDFSMNQGLQWPALMQALALREGGP
          PTFRLTGIGPPAPDNSDHLHEVGCKLAQLAEAIHVEFEYRGFVANSLADLDASMLELRPSDTEAVAVNSV
          FELHKLLGRPGGIEKVLGVVKQIKPVIFTVVEQESNHNGPVFLDRFTESLHYYSTLFDSLEGVPNSQDKV
          MSEVYLGKQICNLVACEGPDRVERHETLSQWGNRFGSSGLAPAHLGSNAFKQASMLLSVFNSGQGYRVEE
          SNGCLMLGWHTRPLITTSAWKLSTAAY 

As you can see you only need to manipulate the line below the bar (|), and make sure that it startings with a "<", then continues with the name (the name format really does not matter, we clean your query's header to turn into the ncbi header in the output files). After that, you have the query in amino acid fasta format, all in caps with all the letters. 

------------------------------------------------------ 
Notes

Note 0.

        Q: why do i have to run it twice? A: we have wildcards based on domains detected, the first one is to identify them,
        and the second is to do all the domain_sorted rules (i.e. meme, iqtree)

Note 1.

    its highly recommended, for accuracy, that you put a gene FAMILY's fastas into the the true_query file
    rather than just a singular gene. however, it's not the end of the world if you cant; we automatically
    do a blastp on your arabidopsis query against the arabidopsis genome, just to find similar arabdiposis genes
    when building trees/meme/msa --to note, we get rid of all duplicates
    it's also highly recommeneded, for clarity, that your gene headers for true_query.fasta are the common gene ID
    i.e. AT2G45160 rather than something like NP_182041.1
Note 2.

    LoadDatasets.smk takes some time (at least many minutes), depending on the annotation level
    and size of the genomes. Example, avocado target + arabidopsis base = ~20 minutes. 
    (mostly became of gmap databse loading). you should expect MainPipeline.smk to be a lot faster (a few minutes max), but
    this also depends on the number of query sequences you upload
    everything else should generally run decently fast.
Note 3.

    sure. all you need to do is comment out what you dont want in
    the top of mainpipeline.smk rule all: input: etc (begins right around line 30). 
    do it line by line, but be aware that if that thing is demanded in line after it, it will still run
Note 4.

    we only select one sequenced version for simplicity, but be aware that it might not be the correct strain or correct
    sequenced dataset that you may want (i.e. you might get West Indian avocado instead of Hass avocado when you enter in 3435)
    you can manually upload your data, but make sure the paths are updated
    and the folder format is wildcard_constraints. Also, the purpose of the project is mainly to automate plant comparitivive
    genomics and gene family characterization studies (it can be run on ANY genome, it's just more common and reasonable to be
        running this on plants, esp. non-model, little annotated, but commerically valuable crops.
------------------------------------------------------------------------------------------

Acknowledgements
The author conducted calls with Dr. Zhiyong Wang, Dr. Jeffrey Groh, Shane Brubaker, and Dr. Lindy Jense for consultation on infrastructure changes, though this release does not necessarily reflect their views, nor does their inclusion constitute an endorsement of this publication's content. The author would also like to thank Dr. Morton Nielson at the Technical University of Denmark for granting access to Deeploc (subcellular localization prediction). 

------------------------------------------------------------------------------------------

Citations

Sequence Alignment & Search (BLAST, DIAMOND, FASTA)
Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2
Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18(4), 366–368. https://doi.org/10.1038/s41592-021-01101-x
Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, Article 421. https://doi.org/10.1186/1471-2105-10-421
Pearson, W. R., & Lipman, D. J. (1988). Improved tools for biological sequence comparison. Proceedings of the National Academy of Sciences, 85(8), 2444–2448. https://doi.org/10.1073/pnas.85.8.2444
Multiple Sequence Alignment & Phylogenetics (MUSCLE, IQ-TREE)
Edgar, R. C. (2004). MUSCLE: Multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research, 32(5), 1792–1797. https://doi.org/10.1093/nar/gkh340
Edgar, R. C. (2021). MUSCLE v5 enables improved estimates of phylogenetic tree confidence by ensemble bootstrapping. bioRxiv. https://doi.org/10.1101/2021.06.20.449169
Minh, B. Q., et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution, 37(5), 1530–1534. https://doi.org/10.1093/molbev/msaa015
Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268–274. https://doi.org/10.1093/molbev/msu300
Python Libraries & Core Infrastructure (Biopython, Pandas, Matplotlib)
Chapman, B., & Chang, J. (2000). Biopython: Python tools for computational biology. ACM SIGBIO Newsletter, 20(2), 15–19.
Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & de Hoon, M. J. L. (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423. https://doi.org/10.1093/bioinformatics/btp163
Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9(3), 90–95. https://doi.org/10.1109/MCSE.2007.55
McKinney, W. (2010). Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference, 56–61. https://doi.org/10.25080/Majora-92bf1922-00a
Genomic File Manipulation (gffutils, pybedtools, SAMtools)
Dale, R. K. (2024). gffutils: GFF and GTF file manipulation and interconversion (Version 0.13) [Computer software]. https://github.com/daler/gffutils
Dale, R. K., Pedersen, B. S., & Quinlan, A. R. (2011). Pybedtools: A flexible Python library for manipulating genomic datasets and annotations. Bioinformatics, 27(24), 3423–3424. https://doi.org/10.1093/bioinformatics/btr539
Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), Article giab008. https://doi.org/10.1093/gigascience/giab008
Genomics Utilities & Short Read Mapping (GMAP, GSNAP)
Wu, T. D., & Nacu, S. (2010). Fast and SNP-tolerant detection of complex variants and splicing in short reads. Bioinformatics, 26(7), 873–881. https://doi.org/10.1093/bioinformatics/btq057
Wu, T. D., & Watanabe, C. K. (2005). GMAP: A genomic mapping and alignment program for mRNA and EST sequences. Bioinformatics, 21(9), 1859–1875. https://doi.org/10.1093/bioinformatics/bti310
Databases & Ontologies (Pfam, RefSeq, Sequence Ontology)
Bateman, A., et al. (2004). The Pfam protein families database. Nucleic Acids Research, 32(Database issue), D138–D141. https://doi.org/10.1093/nar/gkh121
Eilbeck, K., et al. (2005). The Sequence Ontology: A tool for the unification of genome annotations. Genome Biology, 6, R44. https://doi.org/10.1186/gb-2005-6-5-r44
Mistry, J., et al. (2021). Pfam: The protein families database in 2021. Nucleic Acids Research, 49(D1), D412–D419. https://doi.org/10.1093/nar/gkaa913
O'Leary, N. A., Wright, M. W., Brister, J. R., Ciufo, S., Haddad, D., McVeigh, R., ... & Pruitt, K. D. (2016). Reference sequence (RefSeq) database at NCBI: Current status, taxonomic expansion, and functional annotation. Nucleic Acids Research, 44(D1), D733–D745. https://doi.org/10.1093/nar/gkv1189
Evolutionary Analysis & Codon Modeling (PAL2NAL, KaKs_Calculator, MCScanX)
Nei, M., & Gojobori, T. (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Molecular Biology and Evolution, 3(5), 418–426. https://doi.org/10.1093/oxfordjournals.molbev.a040410
Suyama, M., Torrents, D., & Bork, P. (2006). PAL2NAL: Robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research, 34(Web Server issue), W609–W612. https://doi.org/10.1093/nar/gkl315
Wang, D., Zhang, Y., Zhang, Z., Zhu, J., & Yu, J. (2010). KaKs_Calculator 2.0: A toolkit incorporating gamma-series methods and sliding window strategies. Genomics Proteomics Bioinformatics, 8(1), 77–80. https://doi.org/10.1016/S1672-0229(10)60008-3
Wang, Y., Tang, H., DeBarry, J. D., Tan, X., Li, J., Wang, X., ... & Paterson, A H. (2012). MCScanX: A toolkit for detection and evolutionary analysis of gene synteny and collinearity. Nucleic Acids Research, 40(7), e49. https://doi.org/10.1093/nar/gkr1293
Zhang, Z., Li, J., Zhao, X. Q., Wang, J., Wong, G. K., & Yu, J. (2006). KaKs_Calculator: Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics, 4(4), 259–263. https://doi.org/10.1016/S1672-0229(07)60007-2
Workflow Management (Snakemake)
Köster, J., & Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520–2522. https://doi.org/10.1093/bioinformatics/bts480
Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., ... & Köster, J. (2021). Sustainable data analysis with Snakemake. F1000Research, 10, Article 33. https://doi.org/10.12688/f1000research.29032.2
Visualization Tools (pyGenomeViz, DNA Features Viewer)
Yamada, K. D. (2022). pyGenomeViz: A genome visualization python package for comparative genomics. Bioinformatics, 38(19), 4623–4625. https://doi.org/10.1093/bioinformatics/btac488
Zulkower, V., & Rosser, S. (2020). DNA Features Viewer: A sequence annotation formatting and plotting library for Python. Bioinformatics, 36(15), 4350–4352. https://doi.org/10.1093/bioinformatics/btaa404
CRISPR Target Design (FlashFry, Hsu et al., Doench et al.)
Doench, J. G., et al. (2014). Rational design of highly active sgRNAs for CRISPR-Cas9–mediated gene inactivation. Nature Biotechnology, 32, 1262–1267. https://doi.org/10.1038/nbt.3026
Hsu, P. D., et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. Nature Biotechnology, 31, 827–832. https://doi.org/10.1038/nbt.2647
McKenna, A., & Shendure, J. (2018). FlashFry: A fast and flexible tool for large-scale CRISPR target design. BMC Biology, 16, Article 74. https://doi.org/10.1186/s12915-018-0545-0
Motif Discovery & Hidden Markov Models (MEME Suite, HMMER)
Bailey, T. L., & Elkan, C. (1994). Fitting a mixture model by expectation maximization to discover motifs in biopolymers. Proceedings of the International Conference on Intelligent Systems for Molecular Biology, 2, 28–36. PMID: 7584402
Bailey, T. L., Johnson, J., Grant, C. E., & Noble, W. S. (2015). The MEME Suite. Nucleic Acids Research, 43(W1), W39–W49. https://doi.org/10.1093/nar/gkv416
Eddy, S. R. (2011). Accelerated Profile HMM Searches. PLOS Computational Biology, 7(10), e1002195. https://doi.org/10.1371/journal.pcbi.1002195
Primer Design & General Toolkits (Primer3, EMBOSS)
Koressaar, T., & Remm, M. (2007). Enhancements and modifications of primer design program Primer3. Bioinformatics, 23(10), 1289–1291. https://doi.org/10.1093/bioinformatics/btm091
Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics, 16(6), 276–277. https://doi.org/10.1016/S0168-9525(00)02024-2
Untergasser, A., et al. (2012). Primer3—new capabilities and interfaces. Nucleic Acids Research, 40(15), e115. https://doi.org/10.1093/nar/gks596

