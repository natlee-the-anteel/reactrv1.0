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

4. install deeploc2 from DTU (requires academic license, use the readme, can be long to install, so if needed, remove rule from main command). More installation details can be found here: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/
-----------------------------------------------------
Usage Instructions

steps for users
if you want to load new genomes: (1) edit the taxonomy IDs in the config.yaml and the asssembly accession number if needed. (2) delete the folder "reactr/data" if applicable, (3) run 
    
    snakemake -s LoadDatasets.smk --cores 8 --rerun-incomplete --forceall -p

if you are content with the current genomes or don't have one yet: (1) edit in desired protein fasta(s) from the base genome, in the top of the config.yaml, (2) check if it matches the base genome and ideally, the same sequenced version, (3) delete the folder "reactr/output" if applicable, (4) run

    snakemake -s MainPipeline.smk --cores all --rerun-incomplete --forceall -p
    
wait. this should take a few minutes max, though it scales with the number of proteins you query. next, once it says it's complete, then run 

    snakemake -s MainPipeline.smk --cores all --rerun-incomplete --forceall -p

Format Instructions
-----------------------------------------------------
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

Q&A
------------------------------------------
1. We have wildcards based on domains detected, the first one is to identify them,and the second is to do all the domain_sorted rules (i.e. meme, iqtree) its highly recommended, for accuracy, that you put a gene FAMILY's fastas into the the true_query file rather than just a singular gene.
2. Make sure that you know what assembly accession you're using (go to NCBI for more details, example for arabidopsis: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=3702. Enter the Gen Bank ID please (RefSeq ID's preferable whenever (even though they can often be misannotated with tons of unresolved isoforms)). Sidenote: we do have rules built in place to remove isoforms, but it's not 100% perfect yet. Bad annotations WILL cause problems with the pipeline as we're merely analyzing that data.
3. We automatically do a blastp on your arabidopsis query against the arabidopsis genome, just to find similar arabdiposis genes when building trees/meme/msa --to note, we get rid of all duplicates
4. LoadDatasets.smk takes some time (at least many minutes), depending on the annotation level and size of the genomes. This may last up to a few hours for certain large genomes (due to tools that build larger databases). Example: avocado target + arabidopsis base = ~20 minutes. See snakemake validation logs (/reactr_validation_examples/{example_case/}
5. You should expect MainPipeline.smk to be a lot faster (a few minutes max), but this also depends on the number of query sequences you upload (scales really high). For small gene families (<10), expect a few minutes. However, for larger ones, it make more time and may crash your computer depending on your infrastructure capabilites.
6. If you don't want to run all the analysis tools, all you need to do is comment out what you dont want in the top of MainPipeline.smk rule all: input. Ddo it line by line, but be aware that if that thing is demanded in line after it, it will still run we allow only select one sequenced version for simplicity.
7. You can manually upload your data, but make sure the paths are updated and the folder format is wildcard_constraints.
8. The purpose of the project is mainly to automate plant comparitivive genomics and gene family characterization studies (it can be run on ANY genome, it's just more common and reasonable to be running this on plants, esp. non-model, little annotated, but commerically valuable crops.

Further directions
------------------------------------------------------------------------------------------
Please reach out if you or your lab has access to HPCs. We're looking to use the HPCs to run this pipeline en masse and create a larger database that has the characterization profiles of as many species as possible, therby decreasing the need for more independent papers that analyze singular genes with the tools we already integrate. Also haven't conducted as many tests with larger gene families (i.e. 100+ genes per family) yet, due to limiations in personal computing power (the author has a macbook pro). Next, we're looking to increase compatibility with other OS systems (specifically Linux then Windows), which will likely be done through some Docker testing. Finally, a fully drafted manuscript that comprehensively describes this program is in the works

Validation/Case studies Walkthough
--------------------------------------------------------------------------------------------
In each subfolder of /reactr_validation_examples, you can find the output table (similar structure to what you'll get if you run your own sequences and species in it). You'll see stuff like trees and protein properties, motifs and chromosomal localization, etc. These things are all sorted by domain (Pfam database). Addtionally, you'll see three snakemake logs, which detail the 3 commands and their outputs in your terminal (i.e. the step progress )when you run the program. 

I did validation on the DELLA blueberry (Vaccinium darrowii), HOOKLESS tomato (Solanum lycopersicum), and PHY rice gene families . These were generally smaller (due to computational limits stated earlier) gene families. The number of genes, their identity, chromosomal localization, and phylogenetic trees strongly agree with previous literature. The papers can be found here.

        DELLA (Zhou et al., 2024) https://pmc.ncbi.nlm.nih.gov/articles/PMC11360860/
        HOOKLESS (Chaabouni et al., 2016) https://www.sciencedirect.com/science/article/abs/pii/S0176161716300931
        PHY (Takano et al., 2009) https://www.pnas.org/doi/10.1073/pnas.0907378106

Sidenote: there were some discrepancies (relating to the isoforms), in which additional splice variants/isoforms are included during BLAST searches. However, further examinations with the phylogenetic trees makes it obvious to researchers that those are isoforms (we do have some counter measures by removing identical sequences, but small errors i.e. small variations in the C-terminus, do pass through). The speed of all of them was remarkably quick compared to doing it by hand and also was highly accurate. 

Acknowledgements
--------------------------------------------------------------
The author conducted consultation calls with Dr. Zhiyong Wang, Dr. Jeffrey Groh, Shane Brubaker, and Dr. Lindy Jensen for advice on infrastructure changes, though this release does not necessarily reflect their views, nor does their inclusion constitute an endorsement of this code's content. The author would also like to thank Dr. Morton Nielson at the Technical University of Denmark for granting access to Deeploc (subcellular localization prediction). 

Citations:
https://drive.google.com/file/d/1yXfE7mirTI3_arohewrS54e4YhNeUgOE/view?usp=sharing
