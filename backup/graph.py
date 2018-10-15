import collections
import networkx as nx
import matplotlib.pyplot as plt

def getDegree(node):
    return self.ppi.degree(node)

class PPI():
    def __init__(self):
        self.ppi = nx.Graph()
        self.confidences = {}
        
        self.readIn("confidence.ppi.txt")

    def readIn(self, filename):
        with open(filename, 'r') as network:
            for edge in network:
                data = edge.split()
                self.ppi.add_edge(data[0], data[1])
                self.confidences[(data[0], data[1])] = data[2]

    def printGraph(self):
        nx.draw(self.ppi, with_labels=True, font_weight='bold')
        plt.savefig("ppi.png")

    def degreeDist(self):
        degree_sequence = sorted([d for n, d in self.ppi.degree()], reverse=True)  
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = zip(*degreeCount.items())

        plt.bar(deg, cnt, width=0.80, color='b')

        plt.title("Degree Histogram")
        plt.ylabel("Count")
        plt.xlabel("Degree")

        plt.savefig("degreeDist.png") 

    def clusteringCoeff(self):
        coeffs = nx.clustering(self.ppi)
        with open('clusteringCoeffs.txt', 'w') as file:
            for x in coeffs:
                file.write(x + " " + str(coeffs[x]) + "\n")

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

    def shortestPaths():


                 

network = PPI()
#network.printGraph()
#network.degreeDist()
#network.clusteringCoeff()
#network.cliqueFinding()
