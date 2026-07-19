[![DOI](https://zenodo.org/badge/1129249052.svg)](https://doi.org/10.5281/zenodo.18306541)
DOWNLOAD AND INSTALL INSTRUCTIONS:
---------------------------------------------

## Quick install (macOS and Linux)

    git clone https://github.com/natlee-the-anteel/reactrv1.0.git
    cd reactrv1.0
    bash setup.sh

This handles Miniforge installation, environment creation, PATH setup, and 
downloading/building the non-conda dependencies (FlashFry, MCScanX, Pfam-A). 
It's safe to re-run if a step fails partway — it skips anything already 
completed.

Supported platforms: macOS (Apple Silicon) and Linux (x86_64). 
Windows users: run via WSL2 (Windows Subsystem for Linux), which behaves 
as Linux for this pipeline.

## Manual installation (if you'd prefer to run each step yourself, or 
## setup.sh doesn't work for your system)

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

OPTIONAL for subcelllular localization prediction:
Install deeploc2 from DTU (requires academic license). More installation details can be found here: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/
If you chose to proceed with this, you can uncomment out lines 12 and 34 on the MainPipeline.smk

Usage Instructions
-----------------------------------------------------

Steps for users
If you want to load new genomes: (1) edit the taxonomy IDs in the config.yaml and the asssembly accession number if needed. (2) delete the folder "reactr/data" if applicable, (3) run 
    
    snakemake -s LoadDatasets.smk --cores 8 

If you are content with the current genomes or don't have one yet: (1) edit in desired protein fasta(s) from the base genome, in the top of the config.yaml, (2) check if it matches the base genome and ideally, the same sequenced version, (3) delete the folder "reactr/output" if applicable, (4) run

    snakemake -s MainPipeline.smk --cores all 
    
After the snakemake job finishes (shows x out of x jobs complete in terminal, or the .smakemake/logs), you can run it once again.The reason why you have to rerun it is because this is meant to be a way for you to get an estimate of the jobs and time that it will take for the following command to run. Usually, the first MainPipeline run will be a bit more static in terms of amount of jobs it runs. However, the second iteration which has the brunt of the analysis (domain sorted analysis), will scale with the number of domains and proteins identified. It's suggested that you take a look at output/hmmer/domains to do a quick check to make sure that everything's reasonable and ready to proceed. 

    snakemake -s MainPipeline.smk --cores all 

Detailed Capabilites
-----------------------------------------------------
**List of current sequential modules:** MCSCANX (synteny, collinear blocks); pair duplication type identification; Ka/Ks calculations when applicable; double blast calculation to retrieve orthologs of base species in target species and rescue additional homologs of base species; DIAMOND of target and base proteomes; conserved domain identification and sorting; protein motif detection among protein families; multiple sequence alignment and phylogenetic tree generation of gene families of target and base proteomes; gene annotation extraction; mRNA extraction; gene structure visualization, such as exons and introns; chromosomal localization and visualization; protein physiochemical property calculation, such as Pi/Mw; protein promoter detection and motif analysis of promoter; DeepLoc subcellular localization; PCR primer generation for qPCR, validation, and cloning; Primer Search off-target PCR scanning; FlashFry CRISPR gRNA generation and off-target scanning.

**Currently developing:** Gene expression fetching per target ortholog via NCBI GEO (Gene Expression Omnibus); Cis-regulatory element identification with promoter motifs; 3D protein visualization; protein-protein interaction prediction


Format Instructions
-----------------------------------------------------
The only things that you really should need to edit (unless you're directly manipulating to code), is just the config.yaml. Specifically, just the taxonids and assembly IDs (they're ncbi ids, they autodownload all the necesary stuff if you simply run the loaddatasets.smk rule (see above), and the query). Below is an example format: 

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

Notes
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
See Detailed Capabilites for in-development modules. 

One of the big aims of this project is to run this en mass. Due to personal computational limits, we're looking to collaborate with any labs that may have any access to HPCs to run this pipeline en masse and create a larger database that has standard characterization profiles of as many species as possible, therby decreasing the need for more independent papers that analyze singular genes with the tools we already integrate. Please reach out to natlee@nuevaschool.org if you or your lab is interested.

Validation/Case studies Walkthough
--------------------------------------------------------------------------------------------
See the REACTR Validation Summary Table (a pdf in the main REACTR directory) for a summary with % homologs identified, Robinson Foulds distance matrices, and false isoform rates. Raw and example outputs can be found here: https://dataverse.harvard.edu/dataverse/reactr

Acknowledgements
--------------------------------------------------------------
The author conducted consultation calls with Dr. Zhiyong Wang, Dr. Jeffrey Groh, Shane Brubaker, and Dr. Lindy Jensen for advice on infrastructure changes, though this release does not necessarily reflect their views, nor does their inclusion constitute an endorsement of this code's content. The author would also like to thank Dr. Morton Nielson at the Technical University of Denmark for granting access to Deeploc (subcellular localization prediction). 
