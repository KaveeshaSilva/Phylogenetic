import os
from tabulate import tabulate
from ete3 import PhyloTree
import sys

proteinList = ["amino acid permease", "LysR family transcriptional regulator",
               "helix-turn-helix domain-containing protein", "efflux transporter outer membrane subunit"]
robinsonDistanceFile = "robinson_distances.txt"


def calculateRobinsonDistance(trees):
    robinsonDistances = []
    for i in range(len(trees)-1):
        tempList = []
        tempList.append(proteinList[i])
        for k in range(0, i+1):
            tempList.append("-")
        for j in range(i+1, len(trees)):
            tempList.append(trees[i].robinson_foulds(trees[j])[0])
        robinsonDistances.append(tempList)
    return robinsonDistances


def treeGeneration(file):
    file = open(file, "r")
    treeData = file.read()
    tree = PhyloTree(treeData)
    file.close()
    return tree


def printTree(tree, protein, treeFile):
    sys.stdout = open(treeFile, "w")
    print("====== " + protein + " ======")
    print(tree)


def printRobinsonDistance(rf_distancees, robinsonDistanceFile):

    sys.stdout = open(robinsonDistanceFile, "a")
    print("=========="+" Robinson distance "+"==========")
    print(tabulate(rf_distancees, headers=[
          " ", proteinList[0], proteinList[1], proteinList[2], proteinList[3]]))


def main():
    con = sys.stdout
    trees = []
    for protein in proteinList:
        folder = "results\\" + protein
        for f in os.listdir(folder):
            if protein in f:
                tree = treeGeneration(folder + "\\" + f)
                trees.append(tree)
                treeFile = "phylogeneticTrees\\" + protein + "\\" + \
                    protein + " phylogenetic tree" + ".txt"
                printTree(tree, protein, treeFile)

    printRobinsonDistance(calculateRobinsonDistance(
        trees), robinsonDistanceFile)
    sys.stdout.close()


if __name__ == "__main__":
    main()
