import os, sys
import numpy as np

class E():
    def __init__(self, i, name, aux_dir):
        self.I = i
        self.name = name
        self.aux_dir = aux_dir

    def _evaluate_dijkstra_endocardial_activation(self):
        print('Taking as input root node locations...')
        with open(self.aux_dir+'root_node_idxs.txt', 'r') as f:
            data = f.readlines()
        root_node_idxs = data.split().astype(int)
