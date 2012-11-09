"""
    Graph Generator Class

    Developed by:
        Summa Lab, Dept. of Computer Science
        University of New Orleans
        2000 Lakeshore Dr.
        New Orleans, LA 70148

    Authors:
        Jonathan Redmann
        Christopher Summa, PhD.

    Developed from:
        June to November 2012

    Version: Alpha

Code herein copyright (c) 2012, the University of New Orleans, written by David L. Mobley, Shuai Liu, Jonathan Redmann, and Christopher Summa.
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
   * Neither the name of the University of New Orleans nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

try:

    import sys
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    from operator import itemgetter
    import copy 

except ImportError as err:

    sys.stderr.write(str(err))
    raise err


class GraphGenerator4:
    """Generates a networkx graph object that is optimized for calculation planning"""


    def __init__(self, titles, scores, secondaryScores, netCharges, similarityScoresLimit, maxPathLength, knownCompoundsByNameList, debug=True):
        """Constructor contains initialization as well a primary processing calls to generate an optimized graph"""

        # Check to see if the data in the input file is valid
        if not len(titles) == scores.shape[0] == scores.shape[1]:
            raise InputMismatchError(len(titles), scores.shape[0], scores.shape[1])

        # Titles of compounds
        self.titles = titles

        # Net Charges of each compound
        self.netCharges = netCharges

        # Set the number of nodes to the length of the list of compound titles
        self.totalNumberOfCompounds = len(titles)

        # Minimum similarity score limit
        self.similarityScoresLimit = similarityScoresLimit

        # Limit on the path length from any unknown compound to some known compound
        self.maxPathLength = maxPathLength

        # List of nodes representing the known compounds from knownCompoundsByNameList
        self.knownCompoundsByNodeList = map(self.titles.index, knownCompoundsByNameList)

        # List of nodes representing the compounds with unknown free binding energies
        self.unknownCompoundsByNodeList = [compound for compound in range(self.totalNumberOfCompounds) if compound not in self.knownCompoundsByNodeList]

        # Convert scores ndarray from strings to floats
        self.scoresAsFloatsArray = scores.astype(np.float)

        # Trim allPairsWeights to upper triangle
        self.scoresAsFloatsArray = np.triu(self.scoresAsFloatsArray, 1)

        # Convert secondaryScores ndarray from strings to floats
        self.secondaryScoresAsFloatsArray = secondaryScores.astype(np.float)

        # Generate a list of subgraphs from the similarity scores in allPairsWeights
        self.initialSubgraphList = self.generateInitialSubgraphList()

        # A list of lists of the edge weights for each subgraph
        self.subgraphScoresLists = self.generateSubgraphScoresLists(self.initialSubgraphList)

        # Elimintates from each subgraph those edges whose weights are less than the hard limit
        self.removeEdgesBelowHardLimit()

        # Make a new master list of subgraphs now that there may be more disconnected components
        self.workingSubgraphsList = self.generateWorkingSubgraphsList()

        # Make a new master sorted list of edge weights for each subgraph now that there may be new subgraphs
        self.workingSubgraphScoresLists = self.generateSubgraphScoresLists(self.workingSubgraphsList)

        # List of nodes that are not in the cycle cover for a given subgraph
        self.nonCycleNodesSet = set()

        # Remove edges, whose removal does not violate constraints, from the subgraphs, 
        # starting with lowest similarity score first
        self.minimizeEdges()

        # Make a copy of the resultingSubgraphsList for use in connectGraphComponents_brute_force2()
        self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)

        # Merge all of the subgraphs into a single networkx graph object
        self.resultGraph = self.mergeAllSubgraphs()

        # Make a copy of the resultGraph to be used in connectGraphComponents_brute_force2()
        self.copyResultGraph = self.resultGraph.copy()

        # Contains a list of edges added to resultGraph in connectGraphComponents_brute_force()
        # This will be used in connectGraphComponents_brute_force2()
        self.edgesAddedInFirstTreePass = []

        # Connect individual components of the resultGraph
	self.connectSubgraphs();

    # End __init__ def

    def connectSubgraphs(self):
        """Adds edges to the resultGraph to connect as many components of
           the final graph as possible"""

        connectSuccess = self.connectGraphComponents_brute_force()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force()

        connectSuccess = self.connectGraphComponents_brute_force_2()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force_2()

    # End connectSubgraphs def



    def generateInitialSubgraphList(self):
        """Generates a list of subgraphs where all nodes of the same net charge are in one subgraph"""

        compoundsGraph = nx.Graph()

        for i in range(self.totalNumberOfCompounds):

            isKnown = False

            if i in self.knownCompoundsByNodeList: isKnown = True

            compoundsGraph.add_node(i, title=self.titles[i], netCharge=self.netCharges[i], known=isKnown)

            for j in range(self.totalNumberOfCompounds):

                if self.scoresAsFloatsArray[i][j] > 0.0:

                    compoundsGraph.add_edge(i, j, weight=self.scoresAsFloatsArray[i][j] )

        initialSubgraphList = nx.connected_component_subgraphs(compoundsGraph)

        return initialSubgraphList

    # End generateChargeSubgraphList def



    def generateSubgraphScoresLists(self, subgraphList):
        """Generate a list of lists where each inner list is the weights of each edge in
           a given subgraph in the subgraphList, sorted from lowest to highest"""

        subgraphScoresLists = []

        for subgraph in subgraphList:

            weightsDictionary = nx.get_edge_attributes(subgraph, 'weight')
            subgraphWeightsList = [(edge[0], edge[1], weightsDictionary[edge]) for edge in weightsDictionary.iterkeys()]
            subgraphWeightsList.sort(key = lambda entry: entry[2])
            subgraphScoresLists.append(subgraphWeightsList)

        return subgraphScoresLists

    # End generateSubgraphScoresLists def



    def removeEdgesBelowHardLimit(self):
        """Remove edges below hard limit from each subGraph and from each weightsList"""

        for subgraph in self.initialSubgraphList:

            weightsList = self.subgraphScoresLists[self.initialSubgraphList.index(subgraph)]
            index = 0

            for edge in weightsList:

                if edge[2] < self.similarityScoresLimit:

                    subgraph.remove_edge(edge[0],edge[1])
                    index = weightsList.index(edge)

            del weightsList[:index + 1]

    # End removeEdgesBelowHardLimit def



    def generateWorkingSubgraphsList(self):
        """Make a new master list of subgraphs now that there may be more disconnected components"""

        workingSubgraphsList = []

        for subgraph in self.initialSubgraphList:

            newSubgraphList = nx.connected_component_subgraphs(subgraph)

            for newSubgraph in newSubgraphList:

                workingSubgraphsList.append(newSubgraph)

        return workingSubgraphsList

    # End generateWorkingSubgraphsList def



    def minimizeEdges(self):
        """Minimize edges in each subgraph while ensuring constraints are met"""

        #edgesRemovedCount = 0

        for subgraph in self.workingSubgraphsList:

            localKnownCompoundsList = self.getLocalKnownCompounds(subgraph)
            localUnknownCompoundsList = self.getLocalUnknownCompounds(subgraph)
            weightsList = self.workingSubgraphScoresLists[self.workingSubgraphsList.index(subgraph)]
            self.nonCycleNodesSet = self.findNonCyclicNodes(subgraph, localKnownCompoundsList)
            numberOfComponents = nx.number_connected_components(subgraph)

            if len(subgraph.edges()) > 2:   # Graphs must have at least 3 edges to be minimzed

                for edge in weightsList:

                    subgraph.remove_edge(edge[0], edge[1])

                    if self.checkConstraints(subgraph, localKnownCompoundsList, localUnknownCompoundsList, numberOfComponents) == False:

                        subgraph.add_edge(edge[0], edge[1], weight = edge[2])

    # End minimizeEdges def



    def getLocalKnownCompounds(self, subgraph):
        """Create a list of compounds that are known for a given subgraph from the master list of known compounds, list may be empty"""

        localKnownCompounds = [compound for compound in subgraph.nodes() if compound in self.knownCompoundsByNodeList]

        return localKnownCompounds

    # End getLocalKnownCompounds def



    def getLocalUnknownCompounds(self, subgraph):
        """Create a list of compounds that are unknown for a given subgraph from the master list of unknown compounds, list may be empty"""

        localUnknownCompounds = [compound for compound in subgraph.nodes() if compound in self.unknownCompoundsByNodeList]

        return localUnknownCompounds

    # End getLocalUnknownCompounds def



    def checkConstraints(self, subgraph, localKnownCompoundsList, localUnknownCompoundsList, numComp ):
        """Determine if the given subgraph still meets the constraints"""

        constraintsMet = True

        if not self.remainsConnected(subgraph, numComp): constraintsMet = False

	if constraintsMet :

        	if not self.checkCycleCovering(subgraph, localKnownCompoundsList): constraintsMet = False

	if constraintsMet :

        	if localKnownCompoundsList:

            		if not self.checkMaxPathLength(subgraph, localKnownCompoundsList, localUnknownCompoundsList): constraintsMet = False

        	else:

            		if not self.checkMaxDistance(subgraph): constraintsMet = False

        return constraintsMet

    # End checkConstraints def



    def remainsConnected(self, subgraph, numComponents):
        """Determines if the subgraph is still connected after removing an edge"""

        isConnected = False

        if numComponents == nx.number_connected_components(subgraph): isConnected = True

        return isConnected

    # End remainsConnected def



    def findNonCyclicNodes(self, subgraph, localKnownCompoundsList):
        """Finds all nodes in the subgraph that are not in cycles"""

        missingNodesSet = set()
        cycleNodes = []
        cycleList = nx.cycle_basis(subgraph)
        cycleNodes = [node for cycle in cycleList for node in cycle]
        missingNodesSet = set([node for node in subgraph.nodes() if node not in cycleNodes])

        return missingNodesSet

    # End findNonCyclicNodes def



    def checkCycleCovering(self, subgraph, localKnownCompoundsList):
        """Checks if the subgraph has a cycle covering returns boolean"""

        hasCovering = False

        if(not self.findNonCyclicNodes(subgraph, localKnownCompoundsList).difference(self.nonCycleNodesSet)): hasCovering = True

        return hasCovering

    # End checkCycleCovering def



    def checkMaxPathLength(self, subgraph, localKnownCompoundsList, localUnknownCompoundsList):
        """Check to see if the compound graph has paths from all unknowns to some known within the maximum path length"""

        withinMaxPathLength = True
        pathLengthDictionary = nx.all_pairs_shortest_path_length(subgraph)

        if len(pathLengthDictionary.keys()) > 0:

            for compound in localUnknownCompoundsList:

                neighborsDictionary = pathLengthDictionary[compound]
                testLength = lambda x: x <= self.maxPathLength
                nodesToTest = [neighborsDictionary[x] for x in neighborsDictionary if x in localKnownCompoundsList]

                if not filter(testLength, nodesToTest): withinMaxPathLength = False

        return withinMaxPathLength

    # End checkMaxPathLength def



    def checkMaxDistance(self, subgraph):
        """Check to see if the graph has paths from all compounds to all other compounds within a specified limit"""

        withinMaxDistance = True

        if not nx.is_connected(subgraph): print nx.info(subgraph)

        for node in subgraph:

            eccentricity = nx.eccentricity(subgraph, node)

            if eccentricity > self.maxPathLength: withinMaxDistance = False

        return withinMaxDistance

    # End checkMaxDistance def



    def mergeAllSubgraphs(self):
        """Combine networkx graph objects representing each of the subgraphs into
           a single networkx graph object with the former subgraphs being components""" 

        finalGraph = nx.Graph()

        for subgraph in self.workingSubgraphsList:

            finalGraph = nx.union(finalGraph, subgraph)

        return finalGraph

    # End mergeAllSubgraphs def

            

    def connectGraphComponents_brute_force(self):
        """Adds edges to the resultGraph to connect all components that can be connected,
           only one edge is added per component, to form a tree like structure between
           the different components of the resultGraph"""

        self.workingSubgraphsList = nx.connected_component_subgraphs(self.resultGraph)

        if len(self.workingSubgraphsList) == 1:

            return False

        edgesToCheck = []
        edgesToCheckAdditionalInfo = []
        numzeros = 0

        for i in xrange(0,len(self.workingSubgraphsList)-1):

            nodesOfI = self.workingSubgraphsList[i].nodes()

            for j in xrange(i+1,len(self.workingSubgraphsList)):

                nodesOfJ = self.workingSubgraphsList[j].nodes()

                for k in xrange(0,len(nodesOfI)):

                    for l in xrange(0,len(nodesOfJ)):

                        """produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights push this edge into possibleEdgeList """

                        weight = max(self.secondaryScoresAsFloatsArray[nodesOfI[k]][nodesOfJ[l]], self.secondaryScoresAsFloatsArray[nodesOfJ[l]][nodesOfI[k]])

                        if weight > 0.0 :

                          edgesToCheck.append((nodesOfI[k], nodesOfJ[l], weight))
                          edgesToCheckAdditionalInfo.append((nodesOfI[k], nodesOfJ[l], weight, i, j))

                        else :

                          numzeros = numzeros + 1

        if len(edgesToCheck) > 0:

            sortedList = sorted(edgesToCheck, key = itemgetter(2), reverse=True)
            sortedListAdditionalInfo = sorted(edgesToCheckAdditionalInfo, key = itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            edgeToAddAdditionalInfo = sortedListAdditionalInfo[0]
            self.edgesAddedInFirstTreePass.append(edgeToAdd)
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], weight=edgeToAdd[2])
            self.workingSubgraphsList = nx.connected_component_subgraphs(self.resultGraph)

            return True

        else:

            return False

    # End connectGraphComponents_brute_force def



    def connectGraphComponents_brute_force_2(self):
        """Adds a second edge between each of the (former) components of the
           resultGraph to try to provide cycles between (former) components"""

        if len(self.resultingSubgraphsList) == 1:

            return False

        edgesToCheck = []

        for i in xrange(0,len(self.resultingSubgraphsList)-1):

            nodesOfI = self.resultingSubgraphsList[i].nodes()

            for j in xrange(i+1,len(self.resultingSubgraphsList)):

                nodesOfJ = self.resultingSubgraphsList[j].nodes()

                for k in xrange(0,len(nodesOfI)):

                    for l in xrange(0,len(nodesOfJ)):

                        """produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights push this edge into possibleEdgeList """

                        weight = max(self.secondaryScoresAsFloatsArray[nodesOfI[k]][nodesOfJ[l]], self.secondaryScoresAsFloatsArray[nodesOfJ[l]][nodesOfI[k]])

                        if (weight > 0.0):

                            edgesToCheck.append((nodesOfI[k], nodesOfJ[l], weight))

        finalEdgesToCheck = [edge for edge in edgesToCheck if edge not in self.edgesAddedInFirstTreePass]

        if len(finalEdgesToCheck) > 0:

            sortedList = sorted(finalEdgesToCheck, key = itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], weight=edgeToAdd[2])
            self.copyResultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], weight=edgeToAdd[2])
            self.resultingSubgraphsList = nx.connected_component_subgraphs(self.copyResultGraph)

            return True

        else:

            return False

    # End connectGraphComponents_brute_force_2 def



###############################################################

# METHODS FOR OUTPUT AND TO GET DATA ABOUT THE GRAPH

###############################################################



    def generateDotFileForGraph(self, graph, counter):
        """Produces the dot file for a given subgraph"""

        try:

            m = nx.drawing.write_dot.__module__

        except Execption as err:

            sys.stderr.write("Pygraphviz or Pydot were not found.")
            sys.stderr.write("Cannot generate dot files.")
            raise err

        nx.write_dot(graph, "Graph" + str(counter) + ".dot")

    # End generateDotFileForGraph def



    def generateDotFiles(self, offset):
        """Generate Dot Files for all subgraphs in the graph"""

        counter = offset

        for subgraph in self.resultingSubgraphsList:

            self.generateDotFileForGraph(subgraph, counter)

            counter += 1

    # End generateDotFiles def



    def getMasterEdgeList(self):
        """Returns a list of all edges that need to be calculated"""

        masterEdgeList = []

        for subgraph in self.resultingSubgraphsList:

            subgraphEdgeList = []
            titlesDictionary = nx.get_node_attributes(subgraph, 'title')

            for edge in nx.edges(subgraph):

                subgraphEdgeList.append((titlesDictionary[edge[0]],titlesDictionary[edge[1]]))

            masterEdgeList.append(subgraphEdgeList)

        return masterEdgeList

    # End getEdgeList def



    def generateMasterEdgeListFile(self):
        """Creates a file containing the edge list for the graph"""

        try:

            edgeListFile = open('Master Edge List.txt','w')

        except IOError as err:

            sys.stderr.write("Cannot Create Master Edge List File")
            raise err

        masterEdgeList = self.getMasterEdgeList()

        for list in masterEdgeList:

            print >> edgeListFile,"\n"

            for edge in list:

                print >> edgeListFile, edge

    # End generateEdgeListFile def

                    

    def getGraphObject(self):
        """Resturns the networkx object that is the resultGraph"""

        return self.resultGraph

    # End getGraphObject() def

# End Class GraphGenerator Definition



class InputMismatchError(Exception):
    """Custom Exception class for class GraphGenerator when input data does not match."""

    def __init__(self, titlesLength, scoresLengthX, scoresLengthY):

        self.__errorMessage = "Number of titles, scores x dimension, y dimension do not match; values are: "
        self.__value =  ErrorMessage + str(titlesLength) + str(scoresLengthX) + str(scoresLengthY)

    # End __init__ def



    def __str__(self):

        return repr(self.value)

    # End __str__ def

# End Class InputMismatchError definition



