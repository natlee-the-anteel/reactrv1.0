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
MISC Notes

We have wildcards based on domains detected, the first one is to identify them,and the second is to do all the domain_sorted rules (i.e. meme, iqtree) its highly recommended, for accuracy, that you put a gene FAMILY's fastas into the the true_query file rather than just a singular gene. however, it's not the end of the world if you cant; we automatically do a blastp on your arabidopsis query against the arabidopsis genome, just to find similar arabdiposis genes when building trees/meme/msa --to note, we get rid of all duplicates it's also highly recommeneded, for clarity, that your gene headers for true_query.fasta are the common gene ID i.e. AT2G45160 rather than something like NP_182041.1 LoadDatasets.smk takes some time (at least many minutes), depending on the annotation level and size of the genomes. Example, avocado target + arabidopsis base = ~20 minutes. (mostly became of gmap databse loading). you should expect MainPipeline.smk to be a lot faster (a few minutes max), but this also depends on the number of query sequences you upload everything else should generally run decently fast. If you don't want to run all the analysis tools, all you need to do is comment out what you dont want in the top of mainpipeline.smk rule all: input: etc (begins right around line 30). do it line by line, but be aware that if that thing is demanded in line after it, it will still run we allow only select one sequenced version for simplicity, make sure that you know what assembly accession you're using (go to NCBI for more details, example for arabidopsis: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=3702. Enter the Gen Bank ID please). You can manually upload your data, but make sure the paths are updated and the folder format is wildcard_constraints. Also, the purpose of the project is mainly to automate plant comparitivive genomics and gene family characterization studies (it can be run on ANY genome, it's just more common and reasonable to be running this on plants, esp. non-model, little annotated, but commerically valuable crops.
------------------------------------------------------------------------------------------

Acknowledgements:

The author conducted consultation calls with Dr. Zhiyong Wang, Dr. Jeffrey Groh, Shane Brubaker, and Dr. Lindy Jensen for advice on infrastructure changes, though this release does not necessarily reflect their views, nor does their inclusion constitute an endorsement of this code's content. The author would also like to thank Dr. Morton Nielson at the Technical University of Denmark for granting access to Deeploc (subcellular localization prediction). 

Citations:
https://drive.google.com/file/d/1yXfE7mirTI3_arohewrS54e4YhNeUgOE/view?usp=sharing
