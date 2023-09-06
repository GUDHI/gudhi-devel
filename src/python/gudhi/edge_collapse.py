from tqdm import tqdm

def _collapse_edge_list(edges, num:int=0, full:bool=False, strong:bool=False, progress:bool=False):
	"""
	Given an edge list defining a 1 critical 2 parameter 1 dimensional simplicial complex, simplificates this filtered simplicial complex, using filtration-domination's edge collapser.
	"""
	from filtration_domination import remove_strongly_filtration_dominated, remove_filtration_dominated
	n = len(edges)
	if full:
		num = 100
	with tqdm(range(num), total=num, desc="Removing edges", disable=not(progress)) as I:
		for i in I:
			if strong:
				edges = remove_strongly_filtration_dominated(edges) # nogil ?
			else:
				edges = remove_filtration_dominated(edges)
			# Prevents doing useless collapses
			if len(edges) >= n:
				if full and strong:
					strong = False
					n = len(edges)
					# n = edges.size() # len(edges)
				else : 
					break
			else:
				n = len(edges)
				# n = edges.size()
	return edges

