import pandas as pd
from Bio import SeqIO

# extract common_bacteria_set


def getCommonSet(proteins, preteinTableFilePath):
    allSheets = pd.read_excel(preteinTableFilePath,
                              sheet_name=None, usecols="C,D,K")

    commonBacteriaSet = []

    for name, sheet in allSheets.items():
        count = 0
        for protein in proteins:
            for i, row in sheet.iterrows():
                if(protein in row['Protein name']):
                    count += 1
                    break

        if(count == 4):
            commonBacteriaSet.append(name)

    return commonBacteriaSet


# extract postions for each protein in bacterias
def getPositions(preteinTableFilePath, commonBacteriaSet):

    positions = {
        "amino acid permease": {},
        "LysR family transcriptional regulator": {},
        "helix-turn-helix domain-containing protein": {},
        "efflux transporter outer membrane subunit": {}
    }

    for protein in positions:
        for bacteria in commonBacteriaSet:
            sheet = pd.read_excel(preteinTableFilePath,
                                  sheet_name=bacteria, usecols="C,D,K")
            for i, row in sheet.iterrows():
                if(protein in row['Protein name']):
                    start = int(row['Start'])
                    stop = int(row['Stop'])
                    positions[protein][bacteria] = {
                        'Start': start, 'Stop': stop}
                    break

    return positions

# extract gene sequences


def extractGeneSequences(positions):
    for protein in positions:
        for bacteria in positions[protein]:
            filePath = "Data/Sequences/"+bacteria.strip()+".fasta"
            start = positions[protein][bacteria]["Start"]
            stop = positions[protein][bacteria]["Stop"]
            with open(filePath) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequence = record.seq
                    positions[protein][bacteria] = sequence[start-1:stop]
    geneSequences = positions
    return geneSequences

# write gene sequences to fasta files


def writeGeneSequences(geneSequences):
    for protein in geneSequences:
        new_fasta = []
        for bacteria in geneSequences[protein]:
            new_fasta.append('>%s\n%s' %
                             (bacteria, geneSequences[protein][bacteria]))
        outputFile = "homologousGeneSequences/" + \
            '_'.join(protein.split(' '))+".fasta"
        with open(outputFile, 'w') as f:
            f.write('\n'.join(new_fasta))

    return


proteins = ["amino acid permease", "LysR family transcriptional regulator",
            "helix-turn-helix domain-containing protein", "efflux transporter outer membrane subunit"]
preteinTableFilePath = 'Data/protein_tables.xlsx'

print('Getting the common bacteria set...')
commonBacteriaSet = getCommonSet(proteins, preteinTableFilePath)
print('Found common bacteria set\n')
print(commonBacteriaSet)
print('\nGetting positions of gene sequenses...')
positions = getPositions(preteinTableFilePath, commonBacteriaSet)
geneSequences = extractGeneSequences(positions)
print('Writing Gene sequences to fasta files...')
writeGeneSequences(geneSequences)
print('Finished Writing Gene sequences!')
