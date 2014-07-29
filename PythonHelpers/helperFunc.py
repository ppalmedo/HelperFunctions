'''
helper Func: Created July 2014
    helper Func contains helper functions for Python.  Developed primarily for bioinformatic purposes.

@author: ppalmedo
'''

# function to take two lists (i.e. sample or genes) and return indices for them
def findOLIndices(listA,listB):
    indexList = [(listA.index(g), listB.index(g)) for g in listA if g in listB]
    iA, iB = map(list, zip(*indexList))
    return iA,iB