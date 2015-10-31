#Developed by Debasmit Das
#Code forked from https://github.com/keiichishima/pcalg.git 
#This pertains to Problem 2 of HW 2


import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import logging
import sys


from itertools import combinations, permutations

#Libraries for G-Square Test
from gsq.ci_tests import ci_test_bin, ci_test_dis
from gsq.gsq_testdata import bin_data, dis_data

_logger = logging.getLogger(__name__)

def _create_complete_graph(node_ids):
	"""Create a complete graph from the list of node ids.
	Args:
		node_ids: a list of node ids
	Returns:
		An undirected graph (as a networkx.Graph)
	"""
	g = nx.Graph()
	g.add_nodes_from(node_ids)
	for (i, j) in combinations(node_ids, 2):
		g.add_edge(i, j)
		pass
	return g


def estimate_skeleton(indep_test_func, data_matrix, alpha, **kwargs):
	#Estimate a skeleton graph from the sample (Step 1 of IC Algorithm)
	node_ids = range(data_matrix.shape[1])
	g = _create_complete_graph(node_ids)

	node_size = data_matrix.shape[1]
	sep_set = [[set() for i in range(node_size)] for j in range(node_size)]
	l = 0
	while True:
		cont = False
		for (i, j) in permutations(node_ids, 2):
			adj_i = list(g.neighbors(i))
			if j not in adj_i:
			    continue
			else:
			    adj_i.remove(j)
			    pass
			if len(adj_i) >= l:
			    _logger.debug('testing %s and %s' % (i,j))
			    _logger.debug('neighbors of %s are %s' % (i, str(adj_i)))
			    if len(adj_i) < l:
			        continue
			    for k in combinations(adj_i, l):
			        _logger.debug('indep prob of %s and %s with subset %s'
			                      % (i, j, str(k)))
			        p_val = ci_test_bin(dm, i, j, set(k))
			                                
			        _logger.debug('p_val is %s' % str(p_val))
			        if p_val > alpha:
			            if g.has_edge(i, j):
			                _logger.debug('p: remove edge (%s, %s)' % (i, j))
			                g.remove_edge(i, j)
			                pass
			            sep_set[i][j] |= set(k)
			            sep_set[j][i] |= set(k)
			            break
			        pass
			    cont = True
			    pass
			pass
		l += 1
		if cont is False:
			break
		if ('max_reach' in kwargs) and (l > kwargs['max_reach']):
			break
		pass
	return (g, sep_set)

def estimate_cpdag(skel_graph, sep_set):
    """Estimate a CPDAG from the skeleton graph and separation sets
    returned by the estimate_skeleton() function.
    Args:
        skel_graph: A skeleton graph (an undirected networkx.Graph).
        sep_set: An 2D-array of separation set.
            The contents look like something like below.
                sep_set[i][j] = set([k, l, m])
    Returns:
        An estimated DAG.
    """
    dag = skel_graph.to_directed()
    node_ids = skel_graph.nodes()
    for (i, j) in combinations(node_ids, 2):
        adj_i = set(dag.successors(i))
        if j in adj_i:
            continue
        adj_j = set(dag.successors(j))
        if i in adj_j:
            continue
        common_k = adj_i & adj_j
        for k in common_k:
            if k not in sep_set[i][j]:
                if dag.has_edge(k, i):
                    _logger.debug('S: remove edge (%s, %s)' % (k, i))
                    dag.remove_edge(k, i)
                    pass
                if dag.has_edge(k, j):
                    _logger.debug('S: remove edge (%s, %s)' % (k, j))
                    dag.remove_edge(k, j)
                    pass
                pass
            pass
        pass

    def _has_both_edges(dag, i, j):
        return dag.has_edge(i, j) and dag.has_edge(j, i)

    def _has_any_edge(dag, i, j):
        return dag.has_edge(i, j) or dag.has_edge(j, i)

    def _has_one_edge(dag, i, j):
        return ((dag.has_edge(i, j) and (not dag.has_edge(j, i))) or
                (not dag.has_edge(i, j)) and dag.has_edge(j, i))

    def _has_no_edge(dag, i, j):
        return (not dag.has_edge(i, j)) and (not dag.has_edge(j, i))

    # For all the combination of nodes i and j, apply the following
    # rules.
    for (i, j) in combinations(node_ids, 2):
        # Rule 1: For each pair of non-adjacent nodes a and b with a common
        #neighbor c, if the link between a and c has an arrowhead into c and if
        #the link between c and b has no arrowhead into c, then add an
        #arrowhead on the link between c and b pointing to b and mark that
        #link to obtain c (star) b
        # such that k and j are nonadjacent.
        #
        # Check if i-j.
        if _has_both_edges(dag, i, j):
            # Look all the predecessors of i.
            for k in dag.predecessors(i):
                # Skip if there is an arrow i->k.
                if dag.has_edge(i, k):
                    continue
                # Skip if k and j are adjacent.
                if _has_any_edge(dag, k, j):
                    continue
                # Make i-j into i->j
                _logger.debug('R1: remove edge (%s, %s)' % (j, i))
                dag.remove_edge(j, i)
                dag.remove_edge(i,j)
                dag.add_edge(i,j,weight=0)
                break
            pass

        # Rule 2:  If a and b are adjacent and there is a directed path (composed
        #strictly of marked links) from a to b, then add an arrowhead pointing
        #toward b on the link between a and b
        # i->k->j.
        #
        # Check if i-j.
        if _has_both_edges(dag, i, j):
            # Find nodes k where k is i->k.
            succs_i = set()
            for k in dag.successors(i):
                if not dag.has_edge(k, i):
                    succs_i.add(k)
                    pass
                pass
            # Find nodes j where j is k->j.
            preds_j = set()
            for k in dag.predecessors(j):
                if not dag.has_edge(j, k):
                    preds_j.add(k)
                    pass
                pass
            # Check if there is any node k where i->k->j.
            if len(succs_i & preds_j) > 0:
                # Make i-j into i->j
                _logger.debug('R2: remove edge (%s, %s)' % (j, i))
                dag.remove_edge(j, i)
                dag.remove_edge(i,j)
                dag.add_edge(i,j,weight=0)
                break
            pass

        # Rule 3: Orient i-j into i->j whenever there are two chains
        # i-k->j and i-l->j such that k and l are nonadjacent.
        #
        # Check if i-j.
        """if _has_both_edges(dag, i, j):
            # Find nodes k where i-k.
            adj_i = set()
            for k in dag.successors(i):
                if dag.has_edge(k, i):
                    adj_i.add(k)
                    pass
                pass
            # For all the pairs of nodes in adj_i,
            for (k, l) in combinations(adj_i, 2):
                # Skip if k and l are adjacent.
                if _has_any_edge(dag, k, l):
                    continue
                # Skip if not k->j.
                if dag.has_edge(j, k) or (not dag.has_edge(k, j)):
                    continue
                # Skip if not l->j.
                if dag.has_edge(j, l) or (not dag.has_edge(l, j)):
                    continue
                # Make i-j into i->j.
                _logger.debug('R3: remove edge (%s, %s)' % (j, i))
                dag.remove_edge(j, i)
                dag.remove_edge(i,j)
                dag.add_edge(i,j,weight=0)
                break
            pass"""

        # Rule 4: Orient i-j into i->j whenever there are two chains
        # i-k->l and k->l->j such that k and j are nonadjacent.
        #
        # However, this rule is not necessary when the PC-algorithm
        # is used to estimate a DAG.
        pass

    return dag



if __name__ == '__main__':
	import matplotlib.pyplot as plt
	import networkx as nx
	import numpy as np
	np.seterr(divide='ignore', invalid='ignore')

	from gsq.ci_tests import ci_test_bin, ci_test_dis
	from gsq.gsq_testdata import bin_data, dis_data

	# ch = logging.StreamHandler()
	# ch.setLevel(logging.DEBUG)
	# _logger.setLevel(logging.DEBUG)
	# _logger.addHandler(ch)

	n=1000 #Change the number of samples
	px=0.25
	p=0.2#Change the parameter for the disturbances

	#Generate the random samples  
	x0=np.random.binomial(1,px,n)
	y0=np.random.binomial(1,px,n)
	ex0=np.random.binomial(1,p,n)
	ey0=np.random.binomial(1,p,n)
	ex1=np.random.binomial(1,p,n)
	ey1=np.random.binomial(1,p,n)
	ex2=np.random.binomial(1,p,n)
	ey2=np.random.binomial(1,p,n)
	ex3=np.random.binomial(1,p,n)
	ey3=np.random.binomial(1,p,n)
	ex4=np.random.binomial(1,p,n)
	ey4=np.random.binomial(1,p,n)
	ex5=np.random.binomial(1,p,n)
	ey5=np.random.binomial(1,p,n)

	
	#Now, evaluating all the other samples
	x1=np.logical_xor(np.logical_or(x0,y0),ex0).astype(int)
	y1=np.logical_xor(np.logical_or(x0,y0),ey0).astype(int)
	x2=np.logical_xor(np.logical_or(x1,y1),ex1).astype(int)
	y2=np.logical_xor(np.logical_or(x1,y1),ey1).astype(int)
	x3=np.logical_xor(np.logical_or(x2,y2),ex2).astype(int)
	y3=np.logical_xor(np.logical_or(x2,y2),ey2).astype(int)
	x4=np.logical_xor(np.logical_or(x3,y3),ex3).astype(int)
	y4=np.logical_xor(np.logical_or(x3,y3),ey3).astype(int)
	x5=np.logical_xor(np.logical_or(x3,y3),ex4).astype(int)
	y5=np.logical_xor(np.logical_or(x3,y3),ey4).astype(int)

	#Now we have to create the Data Matrix to which
	#we use G-square test to find the conditional independencies
	dm=np.column_stack((x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5))

	

	(g, sep_set) = estimate_skeleton(indep_test_func=ci_test_bin,
		                             data_matrix=dm,
		                             alpha=0.01)

	g = estimate_cpdag(skel_graph=g, sep_set=sep_set)	
	labels={}
	labels[0]=r'$x0$'
	labels[1]=r'$y0$'
	labels[2]=r'$x1$'
	labels[3]=r'$y1$'
	labels[4]=r'$x2$'
	labels[5]=r'$y2$'
	labels[6]=r'$x3$'
	labels[7]=r'$y3$'
	labels[8]=r'$x4$'
	labels[9]=r'$y4$'
	labels[10]=r'$x5$'
	labels[11]=r'$y5$'
	pos=nx.shell_layout(g)	
	#nx.draw_networkx_labels(g,pos,labels)
	nx.draw_networkx(g,pos,arrows='true',with_labels='true',labels=labels)
	plt.show()

	
	

	
	
