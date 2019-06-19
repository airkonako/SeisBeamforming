# -*- coding: utf-8 -*-
"""
Partisioning MPI Domain adaptive to fault geometry
for Okubo et al., 2018 real fault geometry
28 December 2017
"""

def adaptiveMPIdomain_2D(plotfig=False):

    import matplotlib.pyplot as plt
    import random
    import numpy as np
    import sys    
    import math
    import time

    #------Model Values-----#
    meshname = './kaikoura_mesh/mesh_kaikoura.inp'
    fname = 'MPIDomains.input'


    domainsize = [-30e3-100, -30e3-100, 60e3+100, 60e3+100] #Domain size [xmin ymin xmax ymax] which includes buffer space
    np = 64 #128 #Number of mpi domain: 2^n is preferable
    num_of_tri_par_cell_alpha = 1.1 #threshold of maximum num of tri in a mpi cell defined by num_of_tri_par_cell_alpha * avg_tri
    maxbuffer = 1.0 #200 #maxbuffer in MPIDomain context for HOSS
    max_itr = 500 #maximum iteration number for partitioning
    neighbor_minimum_distance =  50 #domain coordinates will merge if two points are within this distance
    #-----------------------#

    def Get_angle(x, y):
        """
        return angle of two vectors x and y
        """
        import numpy as np

        dot_xy = np.dot(x, y)
        norm_x = np.linalg.norm(x)
        norm_y = np.linalg.norm(y)
        cos = dot_xy / (norm_x*norm_y)
        theta = np.degrees(np.arccos(cos))
        
        return theta

    def Get_distance(x1, y1, x2, y2):
        """return distance between two points"""
        return math.sqrt( ((x1-x2)**2)+((y1-y2)**2) )


    class c_mpicell:
        def __init__(self):
            self.xmin = ""
            self.ymin = ""
            self.zmin = 0
            self.xmax = ""
            self.ymax = ""
            self.zmax = 0
            self.num_tri_in_cell = ""
            self.centre_tri_coord_x = []
            self.centre_tri_coord_y = []


    #output file
    print('#---------------------------#')
    print('AdaptiveMPIDomainForHOSS_v1.1')
    print('#---------------------------#')

    #read .inp file
    node_x =  []
    node_y = []
    centre_x_all = []
    centre_y_all = []
    tri = []

    f = open(meshname)
    line = f.readline()

    while line:
        if "NODE" in line:
            while not line.startswith("**"):
                #read node
                if not line.startswith("*"):
                    temp = [float(x) for x in line.split(',')]
                    #print(["%12.8f" % i for i in temp])
                    node_x.append(temp[1])
                    node_y.append(temp[2])
                    line = f.readline()
                else:
                    line = f.readline()


        elif "ELEMENT" in line:
            while not line.startswith("**"):
                tri_temp = []
                #read node
                if not line.startswith("*"):
                    temp = [int(x) for x in line.split(',')]

                    tri_temp.append((temp[1]-1))
                    tri_temp.append((temp[2]-1))
                    tri_temp.append((temp[3]-1))
                    tri.append(tri_temp)
                    #print(tri_temp)
                    line = f.readline()
                else:
                    line = f.readline()

        line = f.readline()
    f.close

    num_tri_all = len(tri)
    avg_tri = round(num_tri_all/np) #avarage of number of tri in one mpi cell

    #derive element distribution by center of triangle
    for i in range(num_tri_all):
        i1 = tri[i][0]
        i2 = tri[i][1]
        i3 = tri[i][2]

        centre_x_all.append((node_x[i1] + node_x[i2] + node_x[i3])/3.0)
        centre_y_all.append((node_y[i1] + node_y[i2] + node_y[i3])/3.0)

    #plt.plot(centre_x_all, centre_y_all, 'k.')
    #plt.show()

    #MPI domain partition algorithm
    """
    while num of mpi cell = np
        1. select mpi cell
        2. evaluate if the num of tri in the cell is more than avg_tri
            -yes:
                devide the domain into two which has same num of tri
            -no:
                finish dividing the mpi cell
    """

    np_itr = 0
    C_mpicell = []    
    C_mpicell_temp = c_mpicell()
    
    #initialize mpi cell by whole domain
    C_mpicell_temp.xmin = domainsize[0]
    C_mpicell_temp.ymin = domainsize[1]
    C_mpicell_temp.zmin = 0
    C_mpicell_temp.xmax = domainsize[2]
    C_mpicell_temp.ymax = domainsize[3]
    C_mpicell_temp.zmax = 0
    C_mpicell_temp.num_tri_in_cell = num_tri_all
    C_mpicell_temp.centre_tri_coord_x = centre_x_all
    C_mpicell_temp.centre_tri_coord_y = centre_y_all
    C_mpicell.append(C_mpicell_temp)
    
    num_mpicell = len(C_mpicell)
    
    #step 0:registar all points for MPI Domains
    DomainCoords = []
    DomainCoords.append([C_mpicell_temp.xmin, C_mpicell_temp.ymin])
    DomainCoords.append([C_mpicell_temp.xmin, C_mpicell_temp.ymax])
    DomainCoords.append([C_mpicell_temp.xmax, C_mpicell_temp.ymin])
    DomainCoords.append([C_mpicell_temp.xmax, C_mpicell_temp.ymax])

    while num_mpicell < np:
        print('num of mpi cell:%d'%num_mpicell)

        #step 1:select mpi cell
        for i in range(num_mpicell):
            #print(i)
            #print(C_mpicell[0].xmin)
            xmin = C_mpicell[i].xmin
            ymin = C_mpicell[i].ymin
            xmax = C_mpicell[i].xmax
            ymax = C_mpicell[i].ymax
            num_tri_in_cell = C_mpicell[i].num_tri_in_cell
            centre_tri_coord_x = C_mpicell[i].centre_tri_coord_x
            centre_tri_coord_y = C_mpicell[i].centre_tri_coord_y

            #step 2:evaluate if the num of tri in the cell is more than avg_tri
            if num_tri_in_cell > avg_tri * num_of_tri_par_cell_alpha:
                #devide the domain into two which has same num of tri
                #consider the direction of partitioning

                angle = Get_angle([xmax - xmin, ymax - ymin],[1, 0])
                
                if angle >0 and angle <= 45.0:
                    direction_of_partition = 'vertical'
                elif angle < 90.0:
                    direction_of_partition = 'horizontal'
                else:
                    print('mpi cell partitioning error.')
                    sys.exit(1) 
                #Monte Carlo approach to find partitioning position
                itr = 0
                R = num_tri_all ** 2
                if direction_of_partition == 'vertical':
                    
                    while itr < max_itr:
                        trial_partition_position = random.uniform(xmin, xmax)
                        trial_mpicell1 = [xmin, ymin, trial_partition_position, ymax]
                        trial_mpicell2 = [trial_partition_position, ymin, xmax, ymax]

                        #count number of tri in new cell
                        num_tri_in_newcell_1 = 0
                        num_tri_in_newcell_2 = 0

                        for j in range(num_tri_in_cell):
                            x_temp = centre_tri_coord_x[j]
                            y_temp = centre_tri_coord_y[j]
                            if trial_mpicell1[0] <= x_temp <= trial_mpicell1[2] and trial_mpicell1[1] <= y_temp <= trial_mpicell1[3]:
                                num_tri_in_newcell_1 = num_tri_in_newcell_1 + 1
                            else:
                                num_tri_in_newcell_2 = num_tri_in_newcell_2 + 1

                        #misfit function R to minimize
                        R_trial = (num_tri_in_newcell_1 - num_tri_in_newcell_2) ** 2
                        
                        if R_trial < R:
                            new_partition_position = trial_partition_position
                            R = R_trial

                        itr = itr + 1


                    #Check the neighbors of new points and merge if they are closer than neighbor_minimum_distance

                    for k in range(len(DomainCoords)):

                        x_test = DomainCoords[k][0]
                        y_test = DomainCoords[k][1]

                        if Get_distance(x_test, y_test, new_partition_position,ymin) < neighbor_minimum_distance \
                            or Get_distance(x_test, y_test, new_partition_position,ymax) < neighbor_minimum_distance:
                            #merge to pre-existing points
                            new_partition_position = x_test
                            print('merge points')


                    #renew C_mpicell

                    centre_tri_coord_newcell1_x = []
                    centre_tri_coord_newcell1_y = []
                    centre_tri_coord_newcell2_x = []
                    centre_tri_coord_newcell2_y = []

                    new_mpicell1 = [xmin, ymin, new_partition_position, ymax]
                    new_mpicell2 = [new_partition_position, ymin, xmax, ymax]

                    DomainCoords.append([new_partition_position, ymin])
                    DomainCoords.append([new_partition_position, ymax])


                    num_tri_in_newcell_1 = 0
                    num_tri_in_newcell_2 = 0

                    for j in range(num_tri_in_cell):
                        x_temp = centre_tri_coord_x[j]
                        y_temp = centre_tri_coord_y[j]
                        if new_mpicell1[0] <= x_temp <= new_mpicell1[2] and new_mpicell1[1] <= y_temp <= new_mpicell1[3]:
                            num_tri_in_newcell_1 = num_tri_in_newcell_1 + 1
                            centre_tri_coord_newcell1_x.append(x_temp)
                            centre_tri_coord_newcell1_y.append(y_temp)
                        else:
                            num_tri_in_newcell_2 = num_tri_in_newcell_2 + 1
                            centre_tri_coord_newcell2_x.append(x_temp)
                            centre_tri_coord_newcell2_y.append(y_temp)

                    del C_mpicell[i]

                    C_mpicell_temp_1 = c_mpicell()
    
                    C_mpicell_temp_1.xmin = new_mpicell1[0]
                    C_mpicell_temp_1.ymin = new_mpicell1[1]
                    C_mpicell_temp_1.zmin = 0
                    C_mpicell_temp_1.xmax = new_mpicell1[2]
                    C_mpicell_temp_1.ymax = new_mpicell1[3]
                    C_mpicell_temp_1.zmax = 0
                    C_mpicell_temp_1.num_tri_in_cell =  num_tri_in_newcell_1
                    C_mpicell_temp_1.centre_tri_coord_x = centre_tri_coord_newcell1_x
                    C_mpicell_temp_1.centre_tri_coord_y = centre_tri_coord_newcell1_y
                    C_mpicell.append(C_mpicell_temp_1)

                    C_mpicell_temp_2 = c_mpicell()
    
                    C_mpicell_temp_2.xmin = new_mpicell2[0]
                    C_mpicell_temp_2.ymin = new_mpicell2[1]
                    C_mpicell_temp_2.zmin = 0
                    C_mpicell_temp_2.xmax = new_mpicell2[2]
                    C_mpicell_temp_2.ymax = new_mpicell2[3]
                    C_mpicell_temp_2.zmax = 0
                    C_mpicell_temp_2.num_tri_in_cell =  num_tri_in_newcell_2
                    C_mpicell_temp_2.centre_tri_coord_x = centre_tri_coord_newcell2_x
                    C_mpicell_temp_2.centre_tri_coord_y = centre_tri_coord_newcell2_y
                    C_mpicell.append(C_mpicell_temp_2)

                else:
                    while itr < max_itr:
                        trial_partition_position = random.uniform(ymin, ymax)
                        trial_mpicell1 = [xmin, ymin, xmax, trial_partition_position]
                        trial_mpicell2 = [xmin, trial_partition_position, xmax, ymax]

                        #count number of tri in new cell
                        num_tri_in_newcell_1 = 0
                        num_tri_in_newcell_2 = 0

                        for j in range(num_tri_in_cell):
                            x_temp = centre_tri_coord_x[j]
                            y_temp = centre_tri_coord_y[j]
                            if trial_mpicell1[0] <= x_temp <= trial_mpicell1[2] and trial_mpicell1[1] <= y_temp <= trial_mpicell1[3]:
                                num_tri_in_newcell_1 = num_tri_in_newcell_1 + 1
                            else:
                                num_tri_in_newcell_2 = num_tri_in_newcell_2 + 1

                        #misfit function R to minimize
                        R_trial = (num_tri_in_newcell_1 - num_tri_in_newcell_2) ** 2
                        
                        if R_trial < R:
                            new_partition_position = trial_partition_position
                            R = R_trial
                        itr = itr + 1

                    for k in range(len(DomainCoords)):

                        x_test = DomainCoords[k][0]
                        y_test = DomainCoords[k][1]

                        if Get_distance(x_test, y_test, xmin, new_partition_position) < neighbor_minimum_distance \
                            or Get_distance(x_test, y_test, xmax, new_partition_position) < neighbor_minimum_distance:
                            #merge to pre-existing points
                            new_partition_position = y_test
                            print('merge points')


                    #renew C_mpicell

                    centre_tri_coord_newcell1_x = []
                    centre_tri_coord_newcell1_y = []
                    centre_tri_coord_newcell2_x = []
                    centre_tri_coord_newcell2_y = []

                    new_mpicell1 = [xmin, ymin, xmax, new_partition_position]
                    new_mpicell2 = [xmin, new_partition_position, xmax, ymax]

                    DomainCoords.append([xmin, new_partition_position])
                    DomainCoords.append([xmax, new_partition_position])

                    num_tri_in_newcell_1 = 0
                    num_tri_in_newcell_2 = 0

                    for j in range(num_tri_in_cell):
                        x_temp = centre_tri_coord_x[j]
                        y_temp = centre_tri_coord_y[j]
                        if new_mpicell1[0] <= x_temp <= new_mpicell1[2] and new_mpicell1[1] <= y_temp <= new_mpicell1[3]:
                            num_tri_in_newcell_1 = num_tri_in_newcell_1 + 1
                            centre_tri_coord_newcell1_x.append(x_temp)
                            centre_tri_coord_newcell1_y.append(y_temp)
                        else:
                            num_tri_in_newcell_2 = num_tri_in_newcell_2 + 1
                            centre_tri_coord_newcell2_x.append(x_temp)
                            centre_tri_coord_newcell2_y.append(y_temp)

                    del C_mpicell[i]

                    C_mpicell_temp_1 = c_mpicell()

                    C_mpicell_temp_1.xmin = new_mpicell1[0]
                    C_mpicell_temp_1.ymin = new_mpicell1[1]
                    C_mpicell_temp_1.zmin = 0
                    C_mpicell_temp_1.xmax = new_mpicell1[2]
                    C_mpicell_temp_1.ymax = new_mpicell1[3]
                    C_mpicell_temp_1.zmax = 0
                    C_mpicell_temp_1.num_tri_in_cell =  num_tri_in_newcell_1
                    C_mpicell_temp_1.centre_tri_coord_x = centre_tri_coord_newcell1_x
                    C_mpicell_temp_1.centre_tri_coord_y = centre_tri_coord_newcell1_y
                    C_mpicell.append(C_mpicell_temp_1)

                    C_mpicell_temp_2 = c_mpicell()

                    C_mpicell_temp_2.xmin = new_mpicell2[0]
                    C_mpicell_temp_2.ymin = new_mpicell2[1]
                    C_mpicell_temp_2.zmin = 0
                    C_mpicell_temp_2.xmax = new_mpicell2[2]
                    C_mpicell_temp_2.ymax = new_mpicell2[3]
                    C_mpicell_temp_2.zmax = 0
                    C_mpicell_temp_2.num_tri_in_cell =  num_tri_in_newcell_2
                    C_mpicell_temp_2.centre_tri_coord_x = centre_tri_coord_newcell2_x
                    C_mpicell_temp_2.centre_tri_coord_y = centre_tri_coord_newcell2_y
                    C_mpicell.append(C_mpicell_temp_2)

                if len(C_mpicell) == np:
                    break

            else:
                continue

            #print(DomainCoords)


        num_mpicell = len(C_mpicell)



    #Output MPI Domain and summary of result
    fo_mpi = open(fname,'w');

    fo_mpi.write('<MPI Domains>\n');
    fo_mpi.write('    ifAutoSetBufferZone =              0\n');
    fo_mpi.write('    Number Of Interval Steps  =       30\n');
    fo_mpi.write('    Maximum Buffer Zone Size = %12.8e\n\n' % maxbuffer);
    fo_mpi.write('    Number Of MPI Domains  =       %d\n' % np)

    fo_mpi.write('    Domain Minimum Coordinate X =\n');

    for i in range(np):
        fo_mpi.write('      %12.8e' %  C_mpicell[i].xmin);
        if (i+1) % 6 == 0:
            fo_mpi.write('\n');

    fo_mpi.write('\n');
    fo_mpi.write('    Domain Minimum Coordinate Y =\n');

    for i in range(np):
        fo_mpi.write('      %12.8e' %  C_mpicell[i].ymin);
        if (i+1) % 6 == 0:
            fo_mpi.write('\n');

    fo_mpi.write('\n');
    fo_mpi.write('    Domain Minimum Coordinate Z =\n');

    for i in range(np):
        fo_mpi.write('      %12.8e' %  0);
        if (i+1) % 6 == 0:
            fo_mpi.write('\n');

    fo_mpi.write('\n');
    fo_mpi.write('    Domain Maximum Coordinate X =\n');

    for i in range(np):
        fo_mpi.write('      %12.8e' %  C_mpicell[i].xmax);
        if (i+1) % 6 == 0:
            fo_mpi.write('\n');
    
    fo_mpi.write('\n');
    fo_mpi.write('    Domain Maximum Coordinate Y =\n');

    for i in range(np):
        fo_mpi.write('      %12.8e' %  C_mpicell[i].ymax);
        if (i+1) % 6 == 0:
            fo_mpi.write('\n');

    fo_mpi.write('\n');
    fo_mpi.write('    Domain Maximum Coordinate Z =\n');

    for i in range(np):
        fo_mpi.write('      %12.8e' %  0);
        if (i+1) % 6 == 0:
            fo_mpi.write('\n');

    fo_mpi.write('\n');
    fo_mpi.write('</MPI Domains>');
    fo_mpi.close();



    numoftri_list = [C_mpicell[n].num_tri_in_cell for n in range(len(C_mpicell))]
    maxnumoftri = max(numoftri_list)
    minnumoftri = min(numoftri_list)
    print('Average number of triangle par mpi cell: %d'%avg_tri)
    print('Mamimum number of triangle par mpi cell: %d'%maxnumoftri)
    print('Minimum number of triangle par mpi cell: %d'%minnumoftri)
    print('MPIDomains.input has been output.')

    if plotfig:
        #plot center point of element and mpi domain
        fig = plt.figure(figsize=(16, 12), dpi=80)
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.85)
        ax.set_title('Average number of triangle par mpi cell:%d\nMamimum number of triangle par mpi cell: %d\nMinimum number of triangle par mpi cell: %d'%(avg_tri,maxnumoftri,minnumoftri))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.axis('equal')
        ax.ticklabel_format(style='sci', scilimits=(-2,2))
        ax.plot(centre_x_all, centre_y_all, 'k.')
        for i in range(np):
            xmin = C_mpicell[i].xmin
            ymin = C_mpicell[i].ymin
            xmax = C_mpicell[i].xmax
            ymax = C_mpicell[i].ymax
            ax.plot([xmin, xmax], [ymin, ymin], 'b-')
            ax.plot([xmax, xmax], [ymin, ymax], 'b-')
            ax.plot([xmax, xmin], [ymax, ymax], 'b-')
            ax.plot([xmin, xmin], [ymax, ymin], 'b-')

        plt.savefig('./MPIDomains.png', dpi=80, facecolor='w', edgecolor='w')
        plt.show()

if __name__ == '__main__':

    plotfig = True
    adaptiveMPIdomain_2D(plotfig)
