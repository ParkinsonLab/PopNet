'''
Created on Jul 30, 2015

@author: javi
'''

inf = float('inf')
import collections as col
import re
import string
Edge = col.namedtuple('Edge', 'start, end, cost')


def parseToGraph(filepath):
    '''
    parse aggregate matrix into edges and nodes
    '''
    
    with open(filepath, 'r') as input:
        data = input.read()
        
    nodes = []
    edges = []
    
    lines = re.split("\n", data)[:-1]
    
    samples = re.split('\t', lines[0])[1:-1]
    
    #make a list of names and edges from weights
    for sample in samples:
        parent_a, parent_b, name = findName(sample)
        nodes.append(name)
        if '-' not in parent_a:
            edges.append((parent_a, name))
            edges.append((parent_b, name))
    
    #make the edge distance matrix
    matrix = {}
    for line, node in zip(lines[1:], nodes):
        matrix[node] = {}
        line_split = [float(x) for x in re.split('\t', line)[1:-1]]
        highest = max(line_split)
                
        for ind, e in enumerate(line_split):
            matrix[node][nodes[ind]] = 1 - (e / highest)
        
    weighted_edges = [(parent, child, matrix[parent][child]) for parent, child in edges]
    
    return weighted_edges, matrix
            

def findName(raw_name):
    
    result = tuple(re.split('_', string.lstrip(raw_name)))
    return result

def getDistance(matrix, path):
    
    if len(path) == 0:
        return inf
    
    previous = path[0]
    total = 0
    while path:
        step = path.pop()
        total += matrix[previous][step]
        previous = step 
    
    return total

class GraphScore():
    
    def __init__(self, edges):
        self.nodes = set(reduce(lambda x, y: x + y, [[edge[0], edge[1]] for edge in edges]))
        self.edges = edges
        
    def dijkstra(self, source, dest):
        source = findName(source)[2]
        dest = findName(dest)[2]
        
        assert source in self.nodes
        dist = {node: inf for node in self.nodes}
        previous = {node: None for node in self.nodes}
        dist[source] = 0
        q = self.nodes.copy()
        neighbours = {node: set() for node in self.nodes}
        for start, end, cost in self.edges:
            neighbours[start].add((end, cost))
            
        while q:
            u = min(q, key = lambda node: dist[node])
            q.remove(u)
            if dist[u] == inf or u == dest:
                break
            for v, cost in neighbours[u]:
                alt = dist[u] + cost
                if alt < dist[v]:
                    dist[v] = alt
                    previous[v] = u
        
        dist = 0
        s, u = col.deque(), dest
        while previous[u]:
            s.appendleft(u)
            u = previous[u]
        s.appendleft(u)
        s.appendleft(source)
        return s


if __name__ == '__main__':
    import GroupComposition as gc
    
    directory = '/data/new/javi/toxo/simulations2/matrix/'
    filename = 'cytoscape/countMatrices/aggregate.txt'
    groupname = 'groups.txt'
#     outputname = ''
    
    filepath = directory + filename
#     outpath = directory + outputname
    grouppath = directory + groupname
    
    groups = gc.loadGroups(grouppath, None)
    del groups['NONE']
    
    edges, matrix = parseToGraph(filepath)
    graph = GraphScore(edges)
    
    print('Starting test')
    
    for group in groups:
        self_strains = groups[group]
        other_strains = reduce(lambda x, y: x + y, [groups[other] for other in groups if other is not group])        
        for strain in self_strains:
            self_dists = []
            other_dists = []
            for self_strain in self_strains:
                self_dists.append(getDistance(matrix, graph.dijkstra(strain, self_strain)))
            for other_strain in other_strains:
                other_dists.append(getDistance(matrix, graph.dijkstra(strain, other_strain)))
        
            if not (max(self_dists) <= min(other_dists)):
                print('Failed: {0}. Max self = {1}, Min Other = {2}'.format(strain, str(max(self_dists)), str(min(other_dists))))
            else:
                print('Passed: {0}. Max self = {1}, Min Other = {2}'.format(strain, str(max(self_dists)), str(min(other_dists))))
    
    print('GraphScore Test Complete.')
    