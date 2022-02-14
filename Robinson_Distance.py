import os
from tabulate import tabulate
from ete3 import PhyloTree
import sys

proteinList = ["amino acid permease", "LysR family transcriptional regulator",
               "helix-turn-helix domain-containing protein", "efflux transporter outer membrane subunit"]
robinsonDistanceFile = "robinson_distances.txt"


def calculateRobinsonDistance(trees):
    robinsonDistances = []
    maximumRobinsonDistances = []
    for i in range(len(trees)-1):
        tempList1 = []
        tempList1.append(proteinList[i])
        tempList2 = []
        tempList2.append(proteinList[i])
        for k in range(0, i+1):
            tempList1.append("-")
            tempList2.append("-")
        for j in range(i+1, len(trees)):
            tempList1.append(trees[i].robinson_foulds(trees[j])[0])
            tempList2.append(trees[i].robinson_foulds(trees[j])[1])
        robinsonDistances.append(tempList1)
        maximumRobinsonDistances.append(tempList2)
    return robinsonDistances, maximumRobinsonDistances


def treeGeneration(file):
    file = open(file, "r")
    treeData = file.read()
    tree = PhyloTree(treeData)
    file.close()
    return tree


def printTree(tree, protein, treeFile):
    sys.stdout = open(treeFile, "a")
    print("====== " + protein + " ======")
    print(tree)


def printRobinsonDistance(robinsonDistances, robinsonDistanceFile):
    robinsonDistances_, maximumRobinsonDistance = robinsonDistances
    sys.stdout = open(robinsonDistanceFile, "a")
    print("=========="+" Robinson distances "+"==========")
    print("\n"*2)
    print(tabulate(robinsonDistances_, headers=[
          " ", proteinList[0], proteinList[1], proteinList[2], proteinList[3]]))
    print("\n"*6)
    print("=========="+" Maximum Robinson distances "+"==========")
    print("\n"*2)
    print(tabulate(maximumRobinsonDistance, headers=[
          " ", proteinList[0], proteinList[1], proteinList[2], proteinList[3]]))


def main():
    con = sys.stdout
    trees = []

    folder = "phylogenetic_trees\\"
    files = os.listdir(os.getcwd()+'\\phylogenetic_trees')

    for file in files:
        tree = treeGeneration(folder + file.split('\\')[-1])
        trees.append(tree)
        protein = (file.split('\\')[-1]).split('.')[0]
        treeFile = "phylogenetic_tree_structures\\" + \
            protein + "_phylogenetic_tree" + ".txt"
        printTree(tree, protein, treeFile)

    printRobinsonDistance(calculateRobinsonDistance(
        trees), robinsonDistanceFile)
    sys.stdout.close()


if __name__ == "__main__":
    folder = 'phylogenetic_tree_structures'
    path = os.path.join(os.getcwd(), folder)
    if(not os.path.isdir(path)):
        os.mkdir(path)
    main()
