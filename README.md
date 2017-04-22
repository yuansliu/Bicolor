# Bicolor

Bi-level error correction for PacBio long reads

## Compile
	make all
	
## Run
	1. sh runCorr.sh
	2. sh runCombine.sh

## Tips
1. GATB should be installed. Please see http://gatb-core.gforge.inria.fr. After installing, the variable *GATB* in the file *Makefile* is the path to search GATB.

2. If errors occur when compiling GATB, the problem may be the version of G++.

3. MUSCLE download from http://drive5.com/muscle. Even we provide an executable file of MUSCLE, it may not support your machine.

4. When runs *sh runCombine.sh*, please make sure that *sh runCorr.sh* must finished. As *muscle* runs in multiple cores with command "&".

5. Before running the codes, please check long read file, short read file and output file in *runCorr.sh* or *runCombine.sh*.

---
**Declaration**<br />
Some parts of the codes are changed from LoRDEC and GATB library.

## Status
Submitted to *InCob 2017*.

### Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>

