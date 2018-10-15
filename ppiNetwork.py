#Written by Aidan Fike
#October 14, 2018

import collections
from numpy.random import randint
import networkx as nx
import matplotlib.pyplot as plt

#Class to represent the protein protein interaction network. Able to read in
#files, test cluster coefficients, degree distributions, shortest paths and
#majority vote algorithms
class PPI():
    #Initialize by reading in confidences and true labels
    def __init__(self):
        self.ppi = nx.Graph()
        self.confidences = {}
        self.trueLabels = {}
        
        self.readInNetwork("confidence.ppi.txt")
        self.readInTrueLabels('labels.txt')

    #Read in the edges and confidences of the graph and use them to construct
    #the entire network. Save the confidences to the dictionary
    #"self.confidences" using the pair of nodes in a tuple as the keys
    def readInNetwork(self, filename):
        with open(filename, 'r') as network:
            for edge in network:
                data = edge.split()
                self.ppi.add_edge(data[0], data[1])
                self.confidences[(data[0], data[1])] = float(data[2])
                self.confidences[(data[1], data[0])] = float(data[2])

    #Read in the true labels of all the proteins in the network and save them
    #to the dictionary trueLabels
    def readInTrueLabels(self, filename):
        with open(filename, 'r') as file:
            for line in file:
                data = line.split() 
                currLabels = []
                protein = data[0]
                for word in data:
                    if word != protein:
                        currLabels.append(int(word))

                self.trueLabels[protein] = currLabels
     

    #Visualize the graph to determine it was correctly constructed
    def printGraph(self):
        nx.draw(self.ppi, with_labels=True, font_weight='bold')
        plt.savefig("ppi.png")

    #Calculate the degree distribution of the ppi network and print the
    #histogram to a png
    def degreeDist(self):
        degree_sequence = sorted([d for n, d in self.ppi.degree()], reverse=True)  
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = zip(*degreeCount.items())

        plt.bar(deg, cnt, width=0.80, color='b')

        plt.title("Degree Histogram")
        plt.ylabel("Count")
        plt.xlabel("Degree")

        plt.savefig("degreeDist.png") 

    #Write the clustering coefficients for each of the proteins in the graph 
    #and write the results to the file clusteringCoeffs.txt
    def clusteringCoeff(self):
        coeffs = nx.clustering(self.ppi)
        with open('clusteringCoeffs.txt', 'w') as file:
            for x in coeffs:
                file.write(x + " " + str(coeffs[x]) + "\n")

    #Find all of the 3-cliques in the graph and write them into the file 
    #triangles.dat.
    def cliqueFinding(self):
        with open('triangles.dat', 'w') as file:
            numTri = 0
            for clique in nx.enumerate_all_cliques(self.ppi):
                if (len(clique) == 3):
                    numTri += 1
                    file.write(str(clique) + "\n")
                    
                if (len(clique) == 4):
                    break

            print "numTri = " + str(numTri)
    
    #Find the lengths of the shortest paths between every node of a randomized 
    #subgraph made up of 1000 nodes from the original ppi network (created in 
    #self.chooseRandom). Write the distribution of these shortest paths in a 
    #histogram, shortestPaths.png
    #
    # Params: sampleSize - The number of nodes that the subgraph is made up of
    # Return: void
    def shortestPaths(self, sampleSize):

        #Create the subgraph based with the nodes based on the randomized 
        #sampling of .chooseRandom()
        subGraph = nx.Graph()
        nodes = self.chooseRandomProteins(sampleSize)
        subGraph.add_nodes_from(nodes)

        #Go through each of the nodes of the subgraph and add an edge between
        #them if a corresponding edge appears in the full ppi network.
        for node in nodes:
            for neighbor in self.ppi.neighbors(node):
                for otherNode in nodes:
                    if neighbor == otherNode:
                        if self.ppi.has_edge(node, neighbor) or \
                                self.ppi.has_edge(neighbor, node):
                            subGraph.add_edge(node, neighbor)

        #Find all of the shortest path lengths of the graph
        shortestPaths = dict(nx.all_pairs_shortest_path_length(subGraph))

        #Add all of the path lengths into a histogram
        pathLengths = []
        for x,v in shortestPaths.items():
            for y in v:
                if x != y:
                    pathLengths.append(shortestPaths[x][y])
        pathLengths = sorted(pathLengths)
        degreeCount = collections.Counter(pathLengths)
        deg, cnt = zip(*degreeCount.items())

        #Print data about the shortest paths discovered
        print "Number of paths", len(pathLengths)
        print "Largest path length: ", pathLengths[len(pathLengths) - 1]

        #Create a graph to display the histogram in shortestPaths.png
        fig, ax = plt.subplots()
        plt.bar(deg, cnt, width=0.60, color='b')
        
        plt.title("Shortest Paths Histogram")
        plt.ylabel("Count")
        plt.xlabel("Shortest Path Length")

        plt.savefig("shortestPaths.png") 
        
    #Choose 1000 random nodes from  the ppi network and add them to a list
    #
    #Params: sampleSize - The number of nodes that the subgraph is made up of
    #Return: A list of 1000 random proteins from the ppi graph
    def chooseRandomProteins(self, sampleSize):
        randNodes = []
        zeros = []
        nonZeros = []

        #Create a list of self.ppi.number_of_nodes() random numbers from 0 to 4
        rand = randint(0, int(self.ppi.number_of_nodes() / sampleSize), \
                                        self.ppi.number_of_nodes())

        #Find all of the zeros in the random list above and add their indicies
        #a list of all zeros. All other values (1-4) have their indicies added
        #to the list of nonZeros
        for i, x in enumerate(rand):
            if x==0:
                zeros.append(i)
            else:
                nonZeros.append(i)

        #If, randomly, more than a sampleSize number of zeros are found, delete
        #some indicies until there are exactly a sampleSize number of indicies
        #Similarly, if randomly less than a sampleSize number of zeros are
        #found, add indicies until exactly a sampleSize number of indicies 
        #is created
        if (len(zeros) > sampleSize):
            for i in range(len(zeros) - sampleSize):
                randomIndex = randint(0, len(zeros))
                del zeros[randomIndex]
        elif (len(zeros) < sampleSize):
            for i in range(sampleSize - len(zeros)):
                randomIndex = randint(0, len(nonZeros))
                zeros.append(nonZeros[randomIndex])
                 
        #Collect protein names corresponding to the random array of indicies 
        #and add them to a list of nodes meant to represent a randomized
        #subgraph
        nodes = list(self.ppi.nodes())
        for x in zeros: 
            randNodes.append(nodes[x])

        return randNodes
 
    #Guess the label of a network protein using a majority vote algorithm. This
    #Algorithm looks at the labels of the neighbors around the given protein
    #and find the most important label based on the number of times a label
    #appears in neighbors
    #
    # Params: protein - the protein who the user desires the guessed label of
    # Return: The integer label of the desired protein
    def voteSingleProtein(self, protein):
        votes = {}
        for neighbor in self.ppi.neighbors(protein):
            for label in self.trueLabels[neighbor]:
                if label in votes.keys():
                    votes[label] += 1
                else: 
                    votes[label] = 1 

        highestLabel = ['abc', 0]
        for guess, weight in votes.items():
            if (weight > highestLabel[1]):
                highestLabel = [guess, weight]
            elif (weight == highestLabel[1]):
                if (guess < highestLabel[0]): 
                    highestLabel = [guess, weight]



    #Guess the label of a network protein using a majority vote algorithm. This
    #Algorithm looks at the labels of the neighbors around the given protein
    #and find the most important label based on the confidences that the given
    #protein has with its neighbors
    #
    # Params: protein - the protein who the user desires the guessed label of
    # Return: The integer label of the desired protein
    def voteSingleProteinImproved(self, protein):
        votes = {}
        for neighbor in self.ppi.neighbors(protein):
            for label in self.trueLabels[neighbor]:
                if label in votes.keys():
                    votes[label] += self.confidences[(neighbor, protein)]
                else: 
                    votes[label] = self.confidences[(neighbor, protein)]

        highestLabel = ['abc', 0]
        for guess, weight in votes.items():
            if (weight > highestLabel[1]):
                highestLabel = [guess, weight]
            elif (weight == highestLabel[1]):
                if (guess < highestLabel[0]): 
                    highestLabel = [guess, weight]

        return highestLabel[0]
       
    #Test the accuracy of the majority vote algorithms using leave one out
    #cross validation. 
    #
    #Params: improved - boolean representing which majority vote algorithm to
    #                   test. True indicates the improved majority vote
    #                   algorithm, false indicates the non-confidence based
    #                   algorithm
    #Return: prints the found accuracy
    def testMajorityVote(self, improved):
        correctGuesses = 0
        for node in self.ppi.nodes:
            if (improved):
                guessedLabel = self.voteSingleProteinImproved(node) 
            else: 
                guessedLabel = self.voteSingleProtein(node) 
            for truelabel in self.trueLabels[node]:
                if(truelabel == guessedLabel):
                    correctGuesses += 1

        if (improved):
            print "Improved majority vote Accuracy: " + \
                    str(float(correctGuesses) / float(len(self.ppi.nodes)))
        else: 
            print "Non-confidence based majority vote Accuracy: " + \
                    str(float(correctGuesses) / float(len(self.ppi.nodes)))
    

#Quick main to test the methods created above
network = PPI()
#network.printGraph()
#network.degreeDist()
#network.clusteringCoeff()
#network.cliqueFinding()
network.shortestPaths(1000)
#network.voteSingleProtein('YDR143C')
#network.voteSingleProteinImproved('YDR143C')
#network.testMajorityVote(True)
#network.testMajorityVote(False)
