# -*- coding: utf-8 -*-
"""
Created on Thu May 18 09:33:10 2017

@author: Moa
"""

class Node(object):
    """
    A node class for the taxonomic plot.
    """
    def __init__(self,parent,name,level):
        self.parent=parent
        if parent!=None:
            self.type=parent.type
        else:
            self.type=None
        self.name=name
        self.level=level
        self.start=0
        self.stop=0
        self.children=[]
        self.counts=0
        self.indexes=[]
        self.serie=None
        self.max=0
                    
class Node_y(object):
    """
    A (smaller) node class for  the year plot.
    """
    def __init__(self,parent,name,level):
        self.parent=parent
        self.name=name
        self.level=level
        self.children=[]