'''
Created on May 15, 2014

@author: javi
'''
'''helper class for the SimilarityContigs script.
Same functions as an iterator except it has a skip function'''

class ContigIterator:

    
    def next(self):
        self.index += 1
        return self.iterable[self.index]
    
    def skip(self, number):
        self.index += number
    
    def __init__(self, param):
        self.iterable = param
        self.index = -1