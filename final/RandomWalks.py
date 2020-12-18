from random import choice, shuffle, random
from collections import deque
import time

#x = read length
#y = length of repeat
#z = number of repeats

#Genome Generation
def genome_gen(genomeLength, bases, stringRepeats, z):
    genome = ""
    count = 0

    zRepeats =  int(genomeLength/z - z)

    while(genomeLength > 0):
        if(count == zRepeats):
            genome += stringRepeats
            genomeLength -= len(stringRepeats)
            count = 0
            continue
        else:
            inp = choice(bases)
            genome += inp
            genomeLength -= len(inp)
            count = count + 1

    return genome

#Repeat Generator
def repeat_gen(length, y, z):
    repeat = genome_gen(y, ['A','C','T','G'], "", -1)
    print("Here is the repeat: " +  repeat + "\n\n")
    return genome_gen(length, ['A', 'C', 'T', 'G'], repeat, z)

def kmerList(gene, x):
    k_mers = []
    for m in range(len(gene)):
        for n in range(m + 1, len(gene) + 1):
            if(len(gene[m:n]) == x):
                k_mers.append(gene[m:n])

    return k_mers

#Create De Bruijn
def createDeBruijn(input):
    edges = input
    graph = {}

    for edge in edges:
        head = edge[:len(edge)-1]
        tail = edge[1:]
        if head in graph:
            graph[head].append(tail)
        else:
            graph[head]=[tail]

    for value in graph.values():
        value.sort()
    return (graph)

def genomeToString(eulCycle):
    genomeString = eulCycle[0]
    for m in range(1, len(eulCycle)):
        genomeString += eulCycle[m][-1]
    return genomeString

#Create Genome from De Bruijn
def createGenomeDeBruijn(graph, firstN, lastN):
    head = lastN[1:]
    tail = firstN[:len(firstN) - 1]
    if head in graph:
        graph[head].append(tail)
    else:
        graph[head]=[tail]

    current = tail
    eulCycle = deque()
    numberOfEdges = sum(map(len, graph.values()))
    while(numberOfEdges > 0):
        possibleChoices = graph[current]

        while possibleChoices:
            eulCycle.append(current)
            numberOfEdges -= 1
            current = possibleChoices.pop()
            possibleChoices = graph.get(current, None)

        if numberOfEdges == 0:
            break

        rotate = 0
        for current in eulCycle:
            if graph[current]:
                break
            rotate += 1

        eulCycle.rotate(-rotate)
    eulCycle.rotate(-eulCycle.index(tail))

    return genomeToString(list(eulCycle))

#Assemble Genome
def genomeAssembly(x, y, z, repeat):
    if(not repeat):
        genome = genome_gen(1000, ['A','C','T','G'], "", -1) #human genome with 1,000 nucleotides
    else:
        genome = repeat_gen(1000, y, z)
    print("Genome: " + genome + "\n\n")
    kmer_List = kmerList(genome, x)
    count = 0
    timeout = time.time() + 60
    originalGenome = ""
    while(genome != originalGenome):
        if(time.time() > timeout):
            print("One moment...")
            break
        graph = createDeBruijn(kmer_List)
        for m in graph.values():
            shuffle(m)
        originalGenome = (createGenomeDeBruijn(graph, kmer_List[0], kmer_List[-1]))
        count += 1
        print("Genome Reconstructed: " + originalGenome + "\n\n")

    print("Random Walks Count: " + str(count))
    if(genome == originalGenome):
        print("Congrats!")
    else:
        print("Error!")

#X needs tail be bigger than 25
def main():
    genomeAssembly(25, 20, 10, True)


if __name__ == "__main__":
    main()
