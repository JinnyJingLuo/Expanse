import os
from shutil import copyfile
import re
from numpy import dot

#unit_v1 = [0.70710678118654746172, -0.70710678118654746172, 0.0000]
#unit_v2 = [0.57735026918962584208, 0.57735026918962584208, 0.57735026918962584208]
#unit_v3 = cross(unit_v1,unit_v2) = [-0.408248290463863, -0.408248290463863, 0.816496580927726]

# This is the path where all the files are stored.
#strr = '10'
length_str = '200' 
stress_str = '100'
path = 'Cases/'+length_str+'_'+stress_str+'_FRsource/'
for strr_num in range(0, 50):
    strr=str(strr_num+1)
    if strr_num < 9:
        file_r = open(path+'/OUTPUT/restart/rs000'+strr+'.data','r')
    else:
        if strr_num < 99:
            file_r = open(path+'OUTPUT/restart/rs00'+strr+'.data','r')
        else:
            file_r = open(path +'OUTPUT/restart/rs0'+strr+'.data','r')
    file_w1 = open(path+length_str+'_'+stress_str+'_'+strr+'.txt','w')
    file_w2 = open(path+length_str+'_'+stress_str+'_links_'+strr+'.txt','w')
    file_w3 = open(path+length_str+'_'+stress_str+'_burgersv_'+strr+'.txt','w')
    
    lines=file_r.readlines()
    if len(lines) > 35:
        node_count = [int(s) for s in lines[16].split() if s.isdigit()]
        #print (node_count[0])
        rs = [[0 for x in range(3)] for y in range(node_count[0]+1)]
        line_num = 35
        for i in range(1,node_count[0]+1):
            l0 = lines[line_num].replace("0,", "")
            numbers_str = l0.split()
            numbers_float = [float(x) for x in numbers_str]
            rs[i] = [numbers_float[1],numbers_float[2],numbers_float[3]]
            
            file_w1.write(str(int(numbers_float[0])))
            file_w1.write(" ")
            file_w1.write(str(rs[i][0]))
            file_w1.write(" ")
            file_w1.write(str(rs[i][1]))
            file_w1.write(" ")
            file_w1.write(str(rs[i][2]))
            file_w1.write('\n')

            for j in range(1, int(numbers_float[4])+1):
                l0 = lines[line_num+2*int(j)].replace("0,", "")
                numbers_str = l0.split()
                numbers_float1 = [float(x) for x in numbers_str]
                # print numbers_float1[0], numbers_float1[1]
                if (numbers_float1[1]<0):
                    file_w2.write(str(int(numbers_float[0])))
                    file_w2.write(" ")
                    file_w2.write(str(int(numbers_float1[0])))
                    file_w2.write('\n')
                    file_w3.write(str((numbers_float1[1])))
                    file_w3.write(" ")
                    file_w3.write(str((numbers_float1[2])))
                    file_w3.write(" ")
                    file_w3.write(str((numbers_float1[3])))
                    file_w3.write('\n')
            line_num = line_num+2+2*int(numbers_float[4])
    file_r.close()
    file_w1.close()
    file_w2.close()
    file_w3.close()
