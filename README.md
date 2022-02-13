## Requirements

install the required modules for run the code

Run `pip install -r requirements.txt`

## Generate Homologous Gene Sequences

Run `python getHomologousGeneSequences.py` 
This command will generate each Homologous Gene Sequence and write them in the homologousGeneSequences folder.
There will be seperate fasta file for each protein in the protein_set.

## Generate Phylogenetic Trees
To generate the Phylogenetic Trees using the UPGMA clustering method, DNA sequence alignments are required.

1.To generate DNA sequence alignments run following command in the main directory. 

`python clustalo.py --email <email_address> --stype dna --sequence <gene_sequence_file_name>`
A new file with .clustal_num extension will be created in the main directory with sequence alignments.

For more details refer https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Clustal+Omega+Help+and+Documentation

2.To generate the Phylogenetic Trees run following comant in the main directory

`python simple_phylogeny.py --email <email_address> --clustering UPGMA --sequence <path_to_the_sequnce_alignment_file_generated_in_above_step>`
A new file with .ph will be created in the main directory with Phylogenetic Trees.

For more details refer https://www.ebi.ac.uk/seqdb/confluence/display/JDSAT/Simple+Phylogeny+Help+and+Documentation#SimplePhylogenyHelpandDocumentation-WebServices

To visualize the Phylogenetic Tree upload the generated .ph file to this website. https://icytree.org/

## Get Robinson Distances

Run `python Robinson_Distance.py` 
This command will generate robinson_distances.txt file in the project folder which include robinson distances. This command also generate phylogenetic trees for each protein as text files.

