import networkx as nx
import pandas as pd

# Create two sample graphs
# G1 = nx.cycle_graph(5)
# G2 = nx.path_graph(5)


def computeGED(el1, el2, labels1=None, labels2=None, size1=None, size2=None):
    # print("Computing GED")
    el1 = pd.DataFrame(el1)
    el2 = pd.DataFrame(el2)

    G1 = nx.Graph()
    for index, row in el1.iterrows():
        G1.add_edge(row['V1'], row['V2'], weight=row['weight'])
    for index, row in labels1.iterrows():
        G1.nodes[row['id']]['label'] = row['label']
    for index, row in size1.iterrows():
        G1.nodes[row['id']]['size'] = row['size']
    # print(G1)
    # for node in G1.nodes(data=True):
        # print(f"Node {node[0]}: Label = {node[1].get('label', 'None')}, Size = {node[1].get('size', 1)}")

    G2 = nx.Graph()
    for index, row in el2.iterrows():
        G2.add_edge(row['V1'], row['V2'], weight=row['weight'])
    for index, row in labels2.iterrows():
        G2.nodes[row['id']]['label'] = row['label']
    for index, row in size2.iterrows():
        G2.nodes[row['id']]['size'] = row['size']
    # print("\n", G2)
    # for node in G2.nodes(data=True):
        # print(f"Node {node[0]}: Label = {node[1].get('label', 'None')}, Size = {node[1].get('size', 1)}")

    ged = nx.optimize_graph_edit_distance(G1, G2, node_subst_cost=node_subst_cost, node_del_cost=node_del_cost,
                                          node_ins_cost=node_ins_cost,
                                          edge_subst_cost=edge_subst_cost, edge_del_cost=edge_del_cost,
                                          edge_ins_cost=edge_ins_cost)

    return next(ged)


"""Function node_subst_cost overrides node_match if specified. 
If neither node_match nor node_subst_cost are specified then default node substitution cost of 0 is used 
(node attributes are not considered during matching).
If node_del_cost is not specified then default node deletion cost of 1 is used. 
If node_ins_cost is not specified then default node insertion cost of 1 is used."""


def node_subst_cost(u, v):
    s1 = u.get('size', 1)
    s2 = v.get('size', 1)
    # print("node subst cost: ", abs(s1 - s2))
    return abs(s1 - s2)
    # com_u = u.get('community', 'None')
    # com_v = v.get('community', 'None')
    # if com_u == com_v:
        # label_u = u.get('label', 'None')
        # label_v = v.get('label', 'None')
        # return 1 if label_u != label_v else 0
    # else:
        # return 2


def node_del_cost(u):
    # print("node deletion cost: ", u.get('size', 1))
    return u.get('size', 1)


def node_ins_cost(v):
    # print("node insertion cost: ", v.get('size', 1))
    return v.get('size', 1)


"""Function edge_subst_cost overrides edge_match if specified. 
If neither edge_match nor edge_subst_cost are specified then default edge substitution cost of 0 is used 
(edge attributes are not considered during matching).
If edge_del_cost is not specified then default edge deletion cost of 1 is used. 
If edge_ins_cost is not specified then default edge insertion cost of 1 is used."""


def edge_subst_cost(e1, e2):
    weight1 = e1.get('weight', 1)
    weight2 = e2.get('weight', 1)
    # print("edge subst cost: ", abs(weight1 - weight2))
    return abs(weight1 - weight2)


def edge_del_cost(e1):
    # print("edge del cost: ", e1.get('weight', 1))
    return e1.get('weight', 1)


def edge_ins_cost(e2):
    # print("edge ins cost: ", e2.get('weight', 1))
    return e2.get('weight', 1)
    # print(e2)
    # u2 = e2.get('end1', None)
    # v2 = e2.get('end2', None)
    # com_u2 = u2.get('community', None)
    # com_v2 = v2.get('community', None)
    # return 2 if com_u2 != com_v2 else 1


# ged = nx.graph_edit_distance(G1, G2, node_subst_cost=node_subst_cost, edge_subst_cost=edge_subst_cost)
# print(list(ged))

"""if labels1:
        for node, label in labels1.items():
            G1.nodes[node]['label'] = label"""

"""    if community1:
        for node, com in community1.items():
            G1.nodes[node]['community'] = com"""
