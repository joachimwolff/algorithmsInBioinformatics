#! /usr/bin/python
# Copyright 2015 Joachim Wolff
# Programming Course: Algorithms in Bioinformatics
# Tutors: Robert Kleinkauf, Omer Alkhnbashi
# Winter semester 2014/2015
#
# Chair of Bioinformatics
# Department of Computer Science
# Faculty of Engineering
# Albert-Ludwig-University Freiburg im Breisgau

from helper import MathHelper


class UpgmaWpgma():
    """Upgma/Wpgma is a clustering method to generate phylogenetic trees. """

    def __init__(self, distance_dictionary, node_count, upgma_wpgma=True, sequence_size_mapping={}):
        """To initalize a object of this class, please define the following:
                distance_dictionary:    A dictionary with the distance between two sequences.
                                        Should have the form \"Key0 key1\":distance. The key0 and key1 have to be integers.
                node_count:             The number of sequences.
                upgma_wpgma:            If True, the upgma weighting is used, if False, wpgma.
                sequence_size_mapping:  Only necessary if wpgma is executed. It defines the size of each sequence.
                                        Should have the form: \"Key:len(sequence)\""""
        self.distance_dictionary = distance_dictionary
        self.mapping = {}
        self.node_count = node_count
        self.number_of_nodes = node_count
        self.upgma_wpgma = upgma_wpgma
        self.sequence_size_mapping = sequence_size_mapping
        self.edge_weight = {}

    def compute_clustering(self):
        """This function computes the clustering to get the phylogenetic tree."""
        computation_is_done = False
        j = 0
        while not computation_is_done:
            j += 1
            minimum_cluster = self.compute_minimal_distance()
            nodes = minimum_cluster[0].split(" ")
            if len(nodes) > 1:
                self.mapping[minimum_cluster[0]] = self.node_count
                self.compute_edge_weight(minimum_cluster[1], nodes)

                if minimum_cluster[0] in self.distance_dictionary:
                    del self.distance_dictionary[minimum_cluster[0]]

                for i in range(0, self.node_count + 1):
                    key_value_0 = nodes[0] + " " + str(i)
                    key_value_1 = nodes[1] + " " + str(i)
                    key_value = self.key_in_dictionary(key_value_0, key_value_1)
                    if key_value[0] != "":
                        key_for_new_cluster_distance = str(i) + " " + str(self.node_count)
                        self.distance_dictionary[key_for_new_cluster_distance] = self.compute_new_distance(
                            self.distance_dictionary[key_value[0]], self.distance_dictionary[key_value[1]], nodes[0],
                            nodes[1])
                        # try:
                        # except:
                        #     "something wring"
                            # "something wring"
                        if not self.upgma_wpgma:
                            self.sequence_size_mapping[self.node_count] = self.sequence_size_mapping[int(nodes[0])] + \
                                                                          self.sequence_size_mapping[int(nodes[1])]
                        del self.distance_dictionary[key_value[0]]
                        del self.distance_dictionary[key_value[1]]
                self.node_count += 1
            else:
                computation_is_done = True

    def key_in_dictionary(self, key_value_0, key_value_1):
        """Returns True if the given keys are in the distance dictionary, False otherwise.
            key_value_0: The first key value.
            key_value_1: The second key value."""
        for i in range(0, 4):
            if key_value_0 in self.distance_dictionary and key_value_1 in self.distance_dictionary:
                return [key_value_0, key_value_1]
            elif key_value_0[::-1] in self.distance_dictionary and key_value_1 in self.distance_dictionary:
                return [key_value_0[::-1], key_value_1]
            elif key_value_0 in self.distance_dictionary and key_value_1[::-1] in self.distance_dictionary:
                return [key_value_0, key_value_1[::-1]]
            elif key_value_0[::-1] in self.distance_dictionary and key_value_1[::-1] in self.distance_dictionary:
                return [key_value_0[::-1], key_value_1[::-1]]
            else:
                return ["", ""]


    def compute_minimal_distance(self):
        """Returns the next two clusters for merging."""
        minimum = ["", MathHelper.Inf]
        for i in self.distance_dictionary:
            if minimum[1] > self.distance_dictionary[i]:
                minimum[0] = i
                minimum[1] = self.distance_dictionary[i]
        return minimum

    def compute_new_distance(self, distance_a_x, distance_b_x, index_a, index_b):
        """Returns the new distance between the new merged cluster and an other cluster.
            distance_a_x:   The old distance between cluster a and x.
            distance_b_x:   The old distance between cluster b and x.
            index_a:        The index of a.
            index_b:        The index of b."""
        if self.upgma_wpgma:
            return self.upgma_distance(distance_a_x, distance_b_x)
        else:
            return self.wpgma_distance(distance_a_x, distance_b_x, self.sequence_size_mapping[int(index_a)],
                                       self.sequence_size_mapping[int(index_b)])

    def upgma_distance(self, distance_a_x, distance_b_x):
        """Returns the upgma-distance between the new merged cluster a and an other cluster x.
            distance_a_x:   The old distance between cluster a and x.
            distance_b_x:   The old distance between cluster b and x."""
        return (distance_a_x + distance_b_x) / 2

    def wpgma_distance(self, distance_a_x, distance_b_x, length_of_a, length_of_b):
        """Returns the wpgma-distance between the new merged cluster a and an other cluster x.
            distance_a_x:   The old distance between cluster a and x.
            distance_b_x:   The old distance between cluster b and x.
            length_of_a:        The index of a.
            length_of_b:        The index of b."""
        return (length_of_a * distance_a_x + length_of_b * distance_b_x) / (length_of_a + length_of_b)

    def compute_edge_weight(self, weight, nodes):
        """This method computes the new edge weight for a new cluster.
            weight: The edge weight equal to the distance of the to merged clusters.
            nodes:  A list containing the indices of the two merged clusters."""
        node0= int(nodes[0])
        node1 = int(nodes[1])
        if node0 < self.number_of_nodes and node1 < self.number_of_nodes:
            # self.edge_weight[self.node_count] = 1
            self.edge_weight[self.node_count] = [weight / float(2), weight / float(2)]
        elif node0 < self.number_of_nodes:
            weightToLeafs = self.edge_weight[node1][1]
            self.edge_weight[self.node_count] = [weight / float(2) - weightToLeafs, weight / float(2)]
        elif node1 < self.number_of_nodes:
            weightToLeafs = self.edge_weight[node0][1]
            self.edge_weight[self.node_count] = [weight / float(2), weight / float(2) - weightToLeafs]
        else:
            weightToLeafs = self.edge_weight[node0][1]
            weightToLeafs1 = self.edge_weight[node1][1]
            self.edge_weight[self.node_count] = [weight / float(2) - weightToLeafs, weight / float(2) - weightToLeafs1]


    def get_newick_tree(self, with_edge_weights=False):
        """Returns the computed cluster in the Newick tree format.
            with_edge_weights:  If True, edge weights are part of the output, if False, not."""
        # expectedValue = {"2 3": 5, "0 1": 7, "4 5": 6, "6 7": 8}
        newick_dictionary = dict([[v, k] for k, v in self.mapping.items()])
        if with_edge_weights:
            for i in newick_dictionary:
                if i in self.edge_weight:
                    nodesWithWeights = newick_dictionary[i].split(" ")
                    nodesWithWeights[0] = nodesWithWeights[0].strip(" ")
                    nodesWithWeights[0] += ":" + str(self.edge_weight[i][1])
                    nodesWithWeights[1] = nodesWithWeights[1].strip(" ")
                    nodesWithWeights[1] += ":" + str(self.edge_weight[i][0])
                    newick_dictionary[i] = nodesWithWeights[0] + " " + nodesWithWeights[1]
            self.mapping = dict([[v, k] for k, v in newick_dictionary.items()])
        for i in self.mapping:
            index = -1
            leading_sequence = True
            for j in newick_dictionary:
                string_to_find = " " + str(self.mapping[i]) + ""
                if newick_dictionary[j].find(string_to_find) != -1:
                    index = j
                    leading_sequence = False
                    break
                string_to_find = str(self.mapping[i]) + " "
                if newick_dictionary[j].find(string_to_find) != -1:
                    index = j
                    leading_sequence = True
                    break
                if with_edge_weights:
                    string_to_find = str(self.mapping[i]) + ":"
                else:
                    string_to_find = str(self.mapping[i]) + ","
                if newick_dictionary[j].find(string_to_find) != -1:
                    index = j
                    leading_sequence = True
                    break
                string_to_find = "," + str(self.mapping[i])
                if newick_dictionary[j].find(string_to_find) != -1:
                    index = j
                    leading_sequence = False
                    break

            if index != -1:
                if leading_sequence:
                    stringToReplace = "(" + newick_dictionary[int(string_to_find.strip().strip(",").strip(":"))].replace(" ",
                                                                                                              ",") + "):"
                else:
                    stringToReplace = ",(" + newick_dictionary[int(string_to_find.strip().strip(",").strip(":"))].replace(" ",
                                                                                                               ",") + ")"
                newick_dictionary[index] = newick_dictionary[index].replace(string_to_find, stringToReplace).replace(
                    ",,", ",")
                del newick_dictionary[int(string_to_find.strip().strip(",").strip(":"))]

        for i in newick_dictionary:
            return "(" + newick_dictionary[i] + ")"
