'''
 Copyright: Copyright (c) 2018
 Created: 2019-7-25
 Author: Zhang_Licheng
 Title: tran mesh in data file to graph
 All rights reserved
'''
import re
import os
from sys import argv,exit

def mesh2graph(mesh_name):

    field_num = 2

    meshfile  = open(mesh_name,  mode='r')

    # skip coor
    node_num = int(meshfile.readline().rstrip('\n'))
    for i in range(node_num+2):
        line = next(meshfile)


    # skip constraint
    line = next(meshfile)
    for i in range(field_num):
        constraint_num = int(next(meshfile).split()[0])
        for j in range(constraint_num+1):
            line = next(meshfile)

        constraint_num = int(next(meshfile).split()[0])
        for j in range(constraint_num+1):
            line = next(meshfile)

    line = next(meshfile)
    line = next(meshfile)
    

    # read mesh
    mesh_list = []

    # for field in field_mesh:
        # for mesh in field:

    elem_num, node_add1 = map(int, line.rstrip('\n').split())

    for i in range(elem_num):

        line = next(meshfile)
        mesh_list.append(line.lstrip().rstrip('\n').split()[1:node_add1])

    meshfile.close()


    # write mesh file metis could read
    metis_mesh_file_name = mesh_name.split('.')[0] + '.mesh'
    meshfile = open(metis_mesh_file_name, mode='w')
    meshfile.write(str(elem_num) + '\n')
    for elem in mesh_list:
        meshfile.write(' '.join(elem) + '\n')
    meshfile.close()

    # generate graph file by 'm2gmetis', need to set it to env
    metis_graph_file_name = mesh_name.split('.')[0] + '.graph'

    command = f"m2gmetis -gtype=nodal {metis_mesh_file_name} {metis_graph_file_name} "
    #command = 'm2gmetis -gtype=nodal {} {}'.format(metis_mesh_file_name, metis_graph_file_name)

    os.system(command)

def main(argvs=None):
    if argvs is None:
        argvs = argv

    mesh2graph(argvs[1])

if __name__ == "__main__":
    exit(main())