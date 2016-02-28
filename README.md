# Algorithms In Bioinformatics
To run the algorithms execute the file "algorithmsInBioinformatics.py" in the folder source/bin.

## Parameters

### Help
  -h, --help            

  Show this help message and exit

### Algorithms

  -a {nw,gotoh,nw3,fengDoolittle,sumOfPairs,upgma,wpgma,nussinov}, --algorithm {nw,gotoh,nw3,fengDoolittle,sumOfPairs,upgma,wpgma,nussinov}

  Define which algorithm should be executed. Options
  are: 'nw' for the algorithm of Needleman and Wunsch,
  'gotoh' for the algorithm of Osamu Gotoh, 'nw3' for
  the Needleman-Wunsch algorithm with three sequences,
  'fengDoolittle' for the heuristic multiple sequence
  alignment algorithm by Da-Fei Feng and Russell F.
  Doolittle,'sumOfPairs' for the scoring of a multiple
  sequence alignment by Humberto Carrillo and David
  Lipman.'upgma' or 'wpgma' is a clustering method to
  generate phylogenetic trees, 'nussinov' for the RNA
  secondary structure prediction algorithm by Ruth
  Nussinov.

### Input file

  -f INPUTFILE, --inputFile INPUTFILE
                        Define the file in which the input sequences are
                        defined. It have to be in fasta-format.

### Output file

  -o OUTPUTFILE, --outputFile OUTPUTFILE
  
  Define in which file the output should be written. If
  not defined, it is written to "outputFile.fas" in the
  local directory.

### Weight function

  -w WEIGHTFUNCTION, --weightFunction WEIGHTFUNCTION
  
  Name of a weight function definde in class
  PairwiseAligmentHelper.

### Gap costs    

  -gc GAPCOST, --gapCost GAPCOST
  
  Name of a gap function definde in class PairwiseAligmentHelper.

### Number of solutions     

  --numberOfSolutions NUMBEROFSOLUTIONS

  Define the number of optimal solutions the Needleman-Wunsch algorithm should compute.

### Output format    

  --outputFormat {graphML,newickTree}

  Define the output format of the output file. This function is only parsed if you choose 'upgma' or 'wpgma' as an algorithm. Default is Newick tree.

### similarity score   

  --similarityScore SIMILARITYSCORE

  Name of a similarity score defined in class PairwiseAligmentHelper.

## Support

If you are having issues, please let me know. Mail adress: wolffj[at]informatik[dot]uni-freiburg[dot]de