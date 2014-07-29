'''
Created on Jul 29, 2014

@author: ppalmedo
'''

def findOLIndices(listA,listB):
    indexList = [(listA.index(g), listB.index(g)) for g in listA if g in listB]
    iA, iB = map(list, zip(*indexList))
    return iA,iB