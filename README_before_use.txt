=================================================================
   IMPORTANT – READ ME **BEFORE** RUNNING ANY PART OF THIS TOOL
=================================================================

CLUES2Companion can be obtained in two ways:

  1. Clone the repository with:
       git clone https://github.com/<your_repo>/Clues2Companion.git
     or download the precompiled release (.zip/.tar.gz) from GitHub.

  2. Download the release on github:
       download the precompiled release (.zip/.tar.gz) from GitHub.

Dependencies:

All the dependencies requires are installable via Conda:
       conda env create -f Clues2Companion.yml
       conda activate Clues2Companion

This will provide all required Python libraries (cyvcf2, numpy, pandas, 
matplotlib, adjustText, etc.) and make CLUES2Companion ready to use.

-----------------------------------------------------------------
FOLDER STRUCTURE
-----------------------------------------------------------------
The script relies on *relative paths*. Please keep the following folders 
exactly where they are, next to this script:

  • Relate/ or Relate-Linux/    (binary + helper scripts)
  • Singer/ or Singer-Linux/    (binary executables)
  • CLUES2/                     (inference.py, RelateToCLUES.py, etc.)
  • required_files/             (ancestral FASTA, recomb maps, masks, etc.)
  • phase1/, phase2/            (produced automatically by previous steps)

Moving or renaming **any** of these directories – or this script itself –
will break Phase 1, Phase 2 and Phase 3.

If you prefer to use your own versions of Relate or Singer, make sure
the folders are renamed exactly as above (Relate/Relate-Linux, 
Singer/Singer-Linux) so that the pipeline can locate them.

-----------------------------------------------------------------
   Happy dating!
-----------------------------------------------------------------
File:  README_before_use.txt
Place: same directory as CLUES2Companion_*.sh
=================================================================
