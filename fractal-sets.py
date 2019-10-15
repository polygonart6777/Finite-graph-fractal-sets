import numpy
from numpy import linalg #Need for Matrix operations
from numpy.linalg import inv #Need to be able to find the inverse of a matrix
from numpy.linalg import matrix_power #Need to be able to find the powers of a
import matplotlib as mpl #Need for plotting
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt
import itertools  #Need for finding all possible combinations of pairs from two lists.


#Putting this here because I will be using each of these for several of the projection functions
A = numpy.array([[1,1,1],[1,0,0],[0,1,0]]);	
#Find the eigenvalues and eigvenvectors
u,v=numpy.linalg.eig(A)
#Combine the real parts of two eigenvectors to define the real matrix similar to A, so that we can work in R^3
real_evector = numpy.transpose(v)[0].real;
im_evector_real_part = numpy.transpose(v)[1].real;
im_evector_im_part = numpy.transpose(v)[1].imag;
real_matrix = [real_evector,im_evector_real_part,im_evector_im_part];

def forward_paths_v1(path_length):
#All forward paths of length, path_length, on the tribonacci graph starting at vertext 1.
	paths = [['edge1'],['edge2'],['edge3']]; 											   	
	while (len(paths[0])<path_length):		   						
#Apply the for loop until paths are the desired length.
		paths1=[];				   						
		paths2=[];
		paths3=[];
		for j in range(len(paths)):
# Run through all of the paths of length j	
# For each path, check the last element and append only the edges that are allowed to follow that last element.																
# There are three possible cases for the last element of each path, corresponding to the the three vertices of the graph.
			if paths[j][-1] in {'edge1','edge2','edge3'}:
				#Need to keep this recorded in each loop in order to make sure all paths are of the same size.
				k = paths[j]; 
				paths1.extend([k+['edge1'],k+['edge4']]);		
			elif paths[j][-1]=='edge4':	
				l = paths[j];				
				paths2.extend([l+['edge2'],l+['edge5']])
			elif paths[j][-1]=='edge5':
				m = paths[j];				
				paths3.append(m+['edge3'])
			
			else:									
				break
		paths=paths1+paths2+paths3;
	return paths


def backward_paths_v1(path_length):
#All forward paths of length, path_length, on the tribonacci graph starting at vertext 1.
	paths = [['edge4'],['edge1']]; 											   	
	while (len(paths[0])<path_length):		   						
#Apply the for loop until paths are the desired length.
		paths1=[];				   						
		paths2=[];
		paths3=[];
		for j in range(len(paths)):
# Run through all of the paths of length j	
# For each path, check the last element and append only the edges that are allowed to follow that last element.																
# There are three possible cases for the last element of each path, corresponding to the the three vertices of the graph.
			if paths[j][-1] in {'edge1','edge4'}:
				k=paths[j];
				paths1.extend([k+['edge1'],k+['edge2'],k+['edge3']]);		
			elif paths[j][-1] in {'edge2', 'edge5'}:
				l= paths[j];					
				paths2.append(l+['edge4'])
			elif paths[j][-1]=='edge3':	
				m= paths[j];	
				paths3.append(m+['edge5'])
			else:									
				break
		paths=paths1+paths2+paths3;
	return paths

def relabel(paths, label_edge1, label_edge2, label_edge3, label_edge4, label_edge5):
#Relabels each of the edges of the paths, according to a chosen label on the edges e1 to e5.
		labellings={'edge1': label_edge1,'edge2': label_edge2,'edge3':label_edge3, 'edge4':label_edge4, 'edge5':label_edge5}																			
		label_list=[[labellings[element] for element in path] for path in paths]	
		return label_list																								



def matrix_sum_forward(labeled_paths):
#Converts each forward path into a vector. 
#For each path, we apply the matrix A to each labels of the path according to the index of the label within the path.  
#Then, we sums up all of these up.    
	A_inverse= numpy.linalg.inv(A);  																																													
	summed=[-sum([matrix_power(A_inverse,i+1).dot(labeled_paths[j][i]) for i in range(len(labeled_paths[j]))]) for j in range(len(labeled_paths))]																																	
	return summed

def matrix_sum_backward(labeled_paths):	
#Converts each backward path into a vector. 
#For each path, we apply the matrix A_inverse to each labels of the path according to the index of the label within the path.  
#Then, we sums up all of these up. 																																																											
	summed=[sum([matrix_power(A,i).dot(labeled_paths[j][i]) for i in range(len(labeled_paths[j]))]) for j in range(len(labeled_paths))]																																	
	return summed												


def projection_expandline(vectors):
#This function projects each of the vector onto the expanding line.		
#By solving a system of equations that give the desired coefficient on the expanding eigenvector.					\		
	proj_list=[numpy.linalg.solve(numpy.transpose(real_matrix),vector)[0]*real_evector for vector in vectors]																	
	return proj_list

	
def projection_contplane(vectors):
#This function projects each of the vector onto the contracting plane. 
#By solving a system of equations that give the desired coefficeint on the expanding eigenvector, and so by using different, 
#the coefficients on the contracting eigenvectors as well. 									
	proj_list=[vector-numpy.linalg.solve(numpy.transpose(real_matrix),vector)[0]*real_evector for vector in vectors]																		
	return proj_list



def projection_xyplane(list):
#Projects vectors from the contracting plane onto the xy-plane.
#rewriting using coordinates [1,0,0] and [0,0,1].  
#I chose to take basis vectors that have the second component as zero, by changing up one of the scalars (\alpha \beta where v= \alpha v1 + \beta v2)						
	proj_list=[];
	for i in range(len(list)):
		solution=numpy.linalg.solve(numpy.transpose(real_matrix),list[i]);
		x=solution[1];
		y=solution[2]*im_evector_im_part[1]/im_evector_real_part[1];
		proj_list.append([x,y])
	return proj_list


#Gives the set of points on the expanding line.
#intervals=projection_expandline(matrix_sum_forward(relabel(forward_paths_v1(7),[0,0,0],[0,0,0],[1,0,0],[0,0,0],[1,0,0])));

#Gives the set of points in the contracting plane.
#rauzy=projection_contplane(matrix_sum_backward(relabel(backward_paths_v1(7),[0,0,0],[0,0,0],[1,0,0],[0,0,0],[1,0,0])));

#Gives the set of points from the contracting plane onto the xy plane.
rauzy_xyplane= projection_xyplane(projection_contplane(matrix_sum_backward(relabel(backward_paths_v1(20),[0,0,1],[0,1,0],[1,0,0],[1,0,0],[0,1,0]))))


#set_combo = list(itertools.product(intervals, rauzy))
#set_sum = [sum(item) for item in set_combo]


def extract_xyz(list):
#This function breaks apart a list of 3d vectors into their x, y and z components and outputs the list of x's,y's,and z's.
	xlist, ylist, zlist = zip(*[list[i] for i in range(len(list))])
	return xlist,ylist,zlist


def extract_xy(list):
#This function breaks apart a list of 2d vectors into their x, y components and outputs the list of x's and y's.
	xlist, ylist= zip(*[list[i] for i in range(len(list))])
	return xlist,ylist


#Break apart x's, y's, and z's for plotting.
#xlist_expand,ylist_expand,zlist_expand =extract_xyz(intervals)
#xlist_cont,ylist_cont,zlist_cont =extract_xyz(rauzy)
#x_list_sum,y_list_sum,z_list_sum = extract_xyz(set_sum)
xlist_plane_cont,ylist_plane_cont = extract_xy(rauzy_xyplane)


#Do the plotting in a single call.
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(xlist_expand,
#           ylist_expand,
#           zlist_expand,marker='.',s=.1)
#plt.show()


#Do the plotting in a single call.
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(xlist_cont,
#           ylist_cont,
#           zlist_cont,marker='.',s=.1)
#plt.show()


#Do the plotting in a single call.
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.scatter(x_list_sum,
#           y_list_sum,
#           z_list_sum,marker='.',s=.1)
#plt.show()


#2d plotting 
x_min=min(min(xlist_plane_cont),min(ylist_plane_cont));
y_min=max(max(xlist_plane_cont),max(ylist_plane_cont));

plt.xlim(x_min, y_min)
plt.ylim(x_min, y_min)
plt.scatter(xlist_plane_cont, ylist_plane_cont, marker = '.',s=.05)
plt.show()





###########Testing



#print(relabel([['e1','e2','e3'],['e1','e4','e5']],[0,0,0],[0,0,0],[0,0,0],[1,0,0],[1,0,0]))	
#print(matrix_sum([[[1,0,0],[0,0,0],[1,0,0]],[[1,0,0],[0,1,0],[1,0,0]]]))
#print(trib_forward_paths_v1(5))
#print(relabel(trib_forward_paths_v1(5),[0,0,0],[0,0,0],[0,0,0],[1,0,0],[1,0,0]))
#print(matrix_sum(relabel(trib_forward_paths_v1(5),[0,0,0],[0,0,0],[0,0,0],[1,0,0],[1,0,0])))