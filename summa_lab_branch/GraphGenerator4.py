"""
Graph Genererator Class:

Developed by:
 Summa Lab, Dept. of Computer Science, University of New Orleans
 2000 Lakeshore Dr.
 New Orleans, LA 70148

Developed from:
 June to November 2012

Authors:
 Jonathan Redmann,
 Christopher Summa, PhD.

"""

try:

    import sys
    import numpy as np
    import networkx as nx
    import matplotlib.pyplot as plt
    import itertools as itr
    from operator import itemgetter
    import copy

except ImportError as err:

    sys.stderr.write(str(err))
    raise err

###############################################################
# CLASS GraphGenerator DEFINITION
###############################################################

class GraphGenerator4:

    def __init__(self, scores, secondaryScores, similarityScoresLimit, maxPathLength, titles, ids, knownCompoundsByNameList, debug=True):

        # Check to see if the data in the input file is valid
        if not len(titles) == scores.shape[0] == scores.shape[1]:
            raise InputMismatchError(len(titles), scores.shape[0], scores.shape[1])

        ##############################
        # Initialize Local Variables
        ##############################

        # Titles of compounds
        self.titles = titles

        # Net Charges of each compound
        # self.netCharges = netCharges

        # Set the number of nodes to the length of the list of compound titles
        self.totalNumberOfCompounds = len(titles)

        # Minimum similarity score limit
        self.similarityScoresLimit = similarityScoresLimit

        # Limit on the path length from any unknown compound to some known compound
        self.maxPathLength = maxPathLength

        # List of nodes representing the known compounds from knownCompoundsByNameList
        self.knownCompoundsByNodeList = map(self.titles.index, knownCompoundsByNameList)

        # List of nodes representing the compounds with unknown free binding energies
        self.unknownCompoundsByNodeList = [compound for compound in xrange(self.totalNumberOfCompounds) if compound not in self.knownCompoundsByNodeList]


        ###############################
        # Process Primary Scores Array
        ###############################

        # Convert scores ndarray from strings to floats
        self.scoresAsFloatsArray = scores.astype(np.float)

        # Trim allPairsWeights to upper triangle
        self.scoresAsFloatsArray = np.triu(self.scoresAsFloatsArray, 1)

        ################################
        # Process Secondary Scores Array
        ################################

        # Convert scores ndarray from strings to floats
        self.secondaryScoresAsFloatsArray = secondaryScores.astype(np.float)


        ################################
        # Process Graph
        ################################

        # Generate a list of subgraphs from the similarity scores in allPairsWeights
        self.initialSubgraphList = self.generateInitialSubgraphList()

        # A list of lists of the edge weights for each subgraph
        self.subgraphScoresLists = self.generateSubgraphScoresLists(self.initialSubgraphList)

        # Elimintates from each subgraph those edges whose weights are less than the hard limit
        self.removeEdgesBelowHardLimit()

        # Make a new master list of subgraphs now that there may be more disconnected components
        self.workingSubgraphsList = self.generateWorkingSubgraphsList()

        # List of nodes that are not in the cycle cover for a given subgraph
        self.nonCycleNodesSet = set()

        # Make a new sorted list of edge weights for each subgraph now that there may be new subgraphs
        self.workingSubgraphScoresLists = self.generateSubgraphScoresLists(self.workingSubgraphsList)

        # Remove edges, whose removal does not violate constraints, from the subgraphs,
        # starting with lowest similarity score first
        self.minimizeEdges()

        # Collect together disjoint subgraphs of like charge into subgraphs
        self.resultingSubgraphsList = copy.deepcopy(self.workingSubgraphsList)

        # Combine seperate subgraphs into a single resulting graph
        self.resultGraph = self.mergeAllSubgraphs()
        #print "Number of Components before bruteforce: "
        #print str(nx.number_connected_components(self.resultGraph))
        #print "Number of edges before brute force: "
        #print self.resultGraph.number_of_edges()
        # Make a copy of the resulting graph for later processing in connectResultingComponents()
        self.copyResultGraph = self.resultGraph.copy()

        # Holds list of edges that were added in the connect components phase
        self.edgesAddedInFirstTreePass = []

        # Add edges to the resultingGraph to connect its components
        # self.connectResultingComponents()
        #self.edgeFile = open('edgeFile.dat', 'w')
        self.connectSubgraphs()

        self.ids = ids

        self.labels = []
        for index in xrange(self.totalNumberOfCompounds):
            self.labels.append((index, self.ids[index]))
        self.labelsDict = {key:value for (key,value) in self.labels}
        self.finalResultGraph = nx.relabel_nodes(self.resultGraph, self.labelsDict)

    # End __init__ def


    ###############################################################
    # CLASS GraphGenerator METHODS FOR GRAPH PROCESSING
    ###############################################################


    def generateInitialSubgraphList(self):
        """Generates a list of subgraphs where all nodes of the same net charge are in one subgraph"""

        compoundsGraph = nx.Graph()

        for i in xrange(self.totalNumberOfCompounds):

            isKnown = False

            if i in self.knownCompoundsByNodeList: isKnown = True

            compoundsGraph.add_node(i, title=self.titles[i], known=isKnown)

            for j in xrange(self.totalNumberOfCompounds):

                if self.scoresAsFloatsArray[i][j] > 0.0:

                    compoundsGraph.add_edge(i, j, similarity=self.scoresAsFloatsArray[i][j] )
        
        #print "Number of edges in the initial graph: "
        #print compoundsGraph.number_of_edges()
        initialSubgraphList = nx.connected_component_subgraphs(compoundsGraph)


        return initialSubgraphList

    # End generateChargeSubgraphList def



    def generateSubgraphScoresLists(self, subgraphList):
        """Generate a list of lists where each inner list is the weights of each edge in
           a given subgraph in the subgraphList, sorted from lowest to highest"""

        subgraphScoresLists = []

        for subgraph in subgraphList:

            weightsDictionary = nx.get_edge_attributes(subgraph, 'similarity')

            subgraphWeightsList = [(edge[0], edge[1], weightsDictionary[edge]) for edge in weightsDictionary.iterkeys()]

            subgraphWeightsList.sort(key = lambda entry: entry[2])

            subgraphScoresLists.append(subgraphWeightsList)


        return subgraphScoresLists

    # End generateSubgraphScoresLists def



    def removeEdgesBelowHardLimit(self):
        """Remove edges below hard limit from each subGraph and from each weightsList"""
        totalEdges = 0
        for subgraph in self.initialSubgraphList:

            weightsList = self.subgraphScoresLists[self.initialSubgraphList.index(subgraph)]

            index = 0

            for edge in weightsList:

                if edge[2] < self.similarityScoresLimit:

                    subgraph.remove_edge(edge[0],edge[1])

                    index = weightsList.index(edge)

            del weightsList[:index + 1]
            totalEdges = totalEdges + subgraph.number_of_edges()
        #print "Number of Edges After Removing Below Hard Limit: " 
        #print  totalEdges
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

        for subgraph in self.workingSubgraphsList:

            localKnownCompoundsList = self.getLocalKnownCompounds(subgraph)

            localUnknownCompoundsList = self.getLocalUnknownCompounds(subgraph)

            weightsList = self.workingSubgraphScoresLists[self.workingSubgraphsList.index(subgraph)]

            self.nonCycleNodesSet = self.findNonCyclicNodes(subgraph, localKnownCompoundsList)

            numberOfComponents = nx.number_connected_components(subgraph)
            #counter = 0
            if len(subgraph.edges()) > 2:   # Graphs must have at least 3 edges to be minimzed

                for edge in weightsList:

                    subgraph.remove_edge(edge[0], edge[1])

                    if self.checkConstraints(subgraph, localKnownCompoundsList, localUnknownCompoundsList, numberOfComponents) == False:

                        subgraph.add_edge(edge[0], edge[1], similarity = edge[2])
                    #else:
                        #counter = counter + 1
                        #print "Edge " + str(counter) + " removed" 

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
        """Determine if the subgraph remains connected after an edge has been removed"""

        isConnected = False

        if numComponents == nx.number_connected_components(subgraph): isConnected = True

        return isConnected

    # End remainsConnected def



    def findNonCyclicNodes(self, subgraph, localKnownCompoundsList):
        """Generates a list of any nodes of the subgraph that are not in cycles"""

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

        for node in subgraph:

            eccentricity = nx.eccentricity(subgraph, node)

            if eccentricity > self.maxPathLength: withinMaxDistance = False

        return withinMaxDistance

    # End checkMaxDistance def



    def mergeAllSubgraphs(self):
        """Generates a single networkx graph object from the subgraphs that have been processed"""

        finalGraph = nx.Graph()

        for subgraph in self.workingSubgraphsList:

            finalGraph = nx.union(finalGraph, subgraph)

        return finalGraph

    # End mergeAllSubgraohs def


    """
    def connectResultingComponents(self):
        \"""Adds edges to the resultGraph to connect its componenets\"""

        connectSuccess = self.connectGraphComponents_brute_force()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force()

        connectSuccess = self.connectGraphComponents_brute_force_2()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force_2()

    # End connectResultingComponents def
    """


    """
    def generateEdgesToCheck(self, iteration):
        \"""Generate a list of edges to check for connecting components of resultGraph\"""

        subgraphsList = []

        if iteration == 1: subgraphsList = self.workingSubgraphsList
        elif iteration == 2: subgraphsList = self.resultingSubgraphsList
        else: print "Bad input on generateEdgesToCheck"

        edgesToCheck = []

        for i in xrange(0,len(self.workingSubgraphsList)-1):

            nodesOfI = self.workingSubgraphsList[i].nodes()

            for j in xrange(i+1,len(self.workingSubgraphsList)):

                nodesOfJ = self.workingSubgraphsList[j].nodes()

                for k in xrange(0,len(nodesOfI)):

                    for l in xrange(0,len(nodesOfJ)):

                        \"""produce an edge from nodesOfI[k] and nodesofJ[l] if nonzero weights

                           push this edge into possibleEdgeList \"""

                        weight = max(self.secondaryScoresAsFloatsArray[nodesOfI[k]][nodesOfJ[l]], self.secondaryScoresAsFloatsArray[nodesOfJ[l]][nodesOfI[k]])

                        if weight > 0.0 :

                          edgesToCheck.append((nodesOfI[k], nodesOfJ[l], weight))

        return edgesToCheck
    """

    def connectSubgraphs(self):
        """
        Adds edges to the resultGraph to connect as many components of the final graph
        as possible
        """

        connectSuccess = self.connectGraphComponents_brute_force()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force()

        connectSuccess = self.connectGraphComponents_brute_force_2()

        while (connectSuccess) :

            connectSuccess = self.connectGraphComponents_brute_force_2()

    # End connectSubgraphs def

    """
    def connectGraphComponents_brute_force(self):

        self.workingSubgraphsList = nx.connected_component_subgraphs(self.resultGraph)

        if len(self.workingSubgraphsList) == 1:

            return False

        edgesToCheck = self.generateEdgesToCheck(1)

        if len(edgesToCheck) > 0:

            sortedList = sorted(edgesToCheck, key = itemgetter(2), reverse=True)

            edgeToAdd = sortedList[0]

            self.edgesAddedInFirstTreePass.append(edgeToAdd)

            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], weight=edgeToAdd[2])

            self.workingSubgraphsList = nx.connected_component_subgraphs(self.resultGraph)

            return True

        else:

            return False

    # End connectGraphComponents_brute_force def



    def connectGraphComponents_brute_force_2(self):

        if len(self.resultingSubgraphsList) == 1:

            return False

        edgesToCheck = self.generateEdgesToCheck(2)

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

    # End connectGraphComponents_brute_force_2
    """

    def connectGraphComponents_brute_force(self):
        """
        Adds edges to the resultGraph to connect all components that can be connected,
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

                        similarity = max(self.secondaryScoresAsFloatsArray[nodesOfI[k]][nodesOfJ[l]], self.secondaryScoresAsFloatsArray[nodesOfJ[l]][nodesOfI[k]])

                        if similarity > 0.0 :

                          edgesToCheck.append((nodesOfI[k], nodesOfJ[l], similarity))
                          edgesToCheckAdditionalInfo.append((nodesOfI[k], nodesOfJ[l], similarity, i, j))

                        else :

                          numzeros = numzeros + 1

        if len(edgesToCheck) > 0:

            sortedList = sorted(edgesToCheck, key = itemgetter(2), reverse=True)
            sortedListAdditionalInfo = sorted(edgesToCheckAdditionalInfo, key = itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            #self.edgeFile.write("\n" + str(edgeToAdd))
            edgeToAddAdditionalInfo = sortedListAdditionalInfo[0]
            self.edgesAddedInFirstTreePass.append(edgeToAdd)
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2])
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

                        similarity = max(self.secondaryScoresAsFloatsArray[nodesOfI[k]][nodesOfJ[l]], self.secondaryScoresAsFloatsArray[nodesOfJ[l]][nodesOfI[k]])

                        if (similarity > 0.0):

                            edgesToCheck.append((nodesOfI[k], nodesOfJ[l], similarity))

        finalEdgesToCheck = [edge for edge in edgesToCheck if edge not in self.edgesAddedInFirstTreePass]

        if len(finalEdgesToCheck) > 0:

            sortedList = sorted(finalEdgesToCheck, key = itemgetter(2), reverse=True)
            edgeToAdd = sortedList[0]
            #self.edgeFile.write("\n" + str(edgeToAdd))
            self.resultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2])
            self.copyResultGraph.add_edge(edgeToAdd[0], edgeToAdd[1], similarity=edgeToAdd[2])
            self.resultingSubgraphsList = nx.connected_component_subgraphs(self.copyResultGraph)

            return True

        else:

            return False

    # End connectGraphComponents_brute_force_2 def



###############################################################
# CLASS GraphGenerator METHODS FOR OUTPUT
###############################################################


    def getGraphObject(self):

        return self.finalResultGraph

    # End getGraphObject def


# End Class GraphGenerator Definition


###############################################################
# CLASS InputMismatchError USED BY CLASS GraphGenerator
###############################################################


class InputMismatchError(Exception):
    """Custom Exception class for class GraphGenerator when input data does not match."""


    def __init__(self, titlesLength, scoresLengthX, scoresLengthY):

        self.__errorMessage = "Number of titles, scores x dimension, y dimension do not match; values are: "
        self.__value =  ErrorMessage + str(titlesLength) + str(scoresLengthX) + str(scoresLengthY)


    def __str__(self):

        return repr(self.value)


# End Class InputMismatchError definition
