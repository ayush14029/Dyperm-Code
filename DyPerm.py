import networkx as nx
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
import operator
from collections import defaultdict
import time
G=nx.Graph()
global comm_list
comm_list=[]

def permanence(G,u,c, comm_list):						#tested, correct
	e_max=0
	numerator=0
	neighbors_u=G.neighbors(u)
	internal_neighbors=len((set(comm_list[c]) & set(neighbors_u)))
	i_neigh=(set(comm_list[c]) & set(neighbors_u))
	d_u=len(neighbors_u)
	for comm in comm_list[0:c]:
		e_c=len(list(set(neighbors_u) & set(comm)))
		if e_c>e_max:
			e_max=e_c
	for comm in comm_list[c+1:]:
		e_c=len(list(set(neighbors_u) & set(comm)))
		if e_c>e_max:
			e_max=e_c
	if e_max==0 and d_u!=0:
		perm=internal_neighbors/float(d_u)
	elif e_max==0 and d_u==0:
		perm=0
	else:
		for node_1 in comm_list[c]:
			for node_2 in comm_list[c]:
				if node_1!=u and node_2!=u and (node_1 in i_neigh) and (node_2 in i_neigh) and (G.has_edge(node_1,node_2) or G.has_edge(node_2,node_1)):
					numerator+=1
		denominator= internal_neighbors*(internal_neighbors-1)/2.0
		if denominator==0:
			denominator=1
		c_in=numerator/float(denominator)
		perm=(internal_neighbors/float(d_u))*(1/float(e_max))-1+c_in
	return perm

def comm_node(G,u,comm_list):				#tested, correct
	for i in range(0,len(comm_list)):
		if u in comm_list[i]:
				return i
	return -1

def perm_comm(G,c,comm_list):				#tested, correct
	perm=0
	nodes=comm_list[c]
	if len(nodes)>0:
		for node in nodes:
			if G.has_node(node):
				perm+=permanence(G,node,c,comm_list)
		perm/=float(len(nodes))
	return perm

def edge_addition(G,u,v, comm_list):
	if comm_node(G,u,comm_list)==-1:
		G.add_node(u)
		comm_list.append([u])
	if comm_node(G,v,comm_list)==-1:
		G.add_node(v)
		comm_list.append([v])
	if comm_node(G,u,comm_list)!=comm_node(G,v,comm_list):
		#check if 'u' moves
		visited={}
		c1=comm_node(G,u,comm_list)
		c2=comm_node(G,v,comm_list)
		c1_perm_old=perm_comm(G,c1,comm_list)
		c2_perm_old=perm_comm(G,c2,comm_list)
		queue=[]
		temp_comm_list=comm_list[:]
		G.add_edge(u,v)
		queue.append(u)
		visited[queue[0]]=0
		while len(queue)>0:
			c=comm_node(G,queue[0],temp_comm_list)
			evaluated=[]
			for vis in visited.keys():
				if visited[vis]==0 and vis!=v:
					p_1=permanence(G,vis,c,temp_comm_list)
					temp_comm_list_new=temp_comm_list[:]
					temp_comm_list_new[c]=list(set(temp_comm_list_new[c])-set([vis]))
					temp_comm_list_new[comm_node(G,v,temp_comm_list_new)]=temp_comm_list_new[comm_node(G,v,temp_comm_list_new)]+[vis]
					p_2=permanence(G,vis,comm_node(G,v,temp_comm_list_new),temp_comm_list_new)
					visited[vis]=1
					if p_2>p_1:
						temp_comm_list[c]=list(set(temp_comm_list[c])-set([vis]))
						temp_comm_list[comm_node(G,v,temp_comm_list)]=temp_comm_list[comm_node(G,v,temp_comm_list)]+[vis]
					else:
						evaluated.append(vis)
			que_add=[]
			for q in queue:
				if q not in evaluated:
					neigh=G.neighbors(q)
					removal=[]
					for n in neigh:
						if (comm_node(G,n,temp_comm_list_new)!=c) or (n in visited) or (n==v):
							removal.append(n)
					neigh=list(set(neigh)-set(removal))
					visited.update({key: 0 for key in neigh})
					que_add+=neigh
			queue=que_add
		diff_c1=c1_perm_old-perm_comm(G,c1,temp_comm_list)
		visited={}
		c1=comm_node(G,u,comm_list)
		c2=comm_node(G,v,comm_list)
		c1_perm_old=perm_comm(G,c1,comm_list)
		c2_perm_old=perm_comm(G,c2,comm_list)
		queue=[]
		temp_comm_list2=comm_list[:]
		queue.append(v)
		visited[queue[0]]=0
		while len(queue)>0:
			c=comm_node(G,queue[0],temp_comm_list2)
			evaluated=[]
			for vis in visited.keys():
				if visited[vis]==0 and vis!=u:
					p_1=permanence(G,vis,c,temp_comm_list2)
					temp_comm_list_new=temp_comm_list2[:]
					temp_comm_list_new[c]=list(set(temp_comm_list_new[c])-set([vis]))
					temp_comm_list_new[comm_node(G,u,temp_comm_list_new)]=temp_comm_list_new[comm_node(G,u,temp_comm_list_new)]+[vis]
					p_2=permanence(G,vis,comm_node(G,u,temp_comm_list_new),temp_comm_list_new)
					visited[vis]=1
					if p_2>p_1:
						temp_comm_list2[c]=list(set(temp_comm_list2[c])-set([vis]))
						temp_comm_list2[comm_node(G,u,temp_comm_list2)]=temp_comm_list2[comm_node(G,u,temp_comm_list2)]+[vis]
					else:
						evaluated.append(vis)

			que_add=[]
			for q in queue:
				if q not in evaluated:
					neigh=G.neighbors(q)
					removal=[]
					for n in neigh:
						if (comm_node(G,n,temp_comm_list2)!=c) or (n in visited) or (n==u):
							removal.append(n)
					neigh=list(set(neigh)-set(removal))
					visited.update({key: 0 for key in neigh})
					que_add+=neigh
			queue=que_add
		diff_c2=c2_perm_old-perm_comm(G,c2,temp_comm_list2)
		#retain community structure having greater difference
		if diff_c1>diff_c2:
			comm_list=temp_comm_list
		else:
			comm_list=temp_comm_list2
	else:
		G.add_edge(u,v)
	return G, comm_list



def edge_deletion(G,u,v,comm_list):							#tested, correct
	if G.has_edge(u,v) and G.degree(u)==1 and G.degree(v)==1:
		c1=comm_node(G,u,comm_list)
		comm_list=comm_list[0:c1]+comm_list[(c1+1):]
		comm_list.append([u])
		comm_list.append([v])
	elif G.has_edge(u,v) and comm_node(G,u,comm_list)==comm_node(G,v,comm_list) and G.degree(v)==1 and comm_node(G,u,comm_list)!=-1:
		index=comm_node(G,v,comm_list)
		l=comm_list[index]
		l.remove(v)
		if len(l)>0:
			comm_list[index]=l
		comm_list.append([v])
	elif G.has_edge(u,v) and comm_node(G,u,comm_list)==comm_node(G,v,comm_list) and G.degree(u)==1 and comm_node(G,u,comm_list)!=-1:
		index=comm_node(G,u,comm_list)
		l=comm_list[index]
		l.remove(u)
		if len(l)>0:
			comm_list[index]=l
		comm_list.append([u])
	elif G.has_edge(u,v) and comm_node(G, u, comm_list)==comm_node(G,v,comm_list):
		visited={}
		c1=comm_node(G,u,comm_list)
		#c2=comm_node(G,v,comm_list)
		c1_perm_old=perm_comm(G,c1,comm_list)
		queue=[]
		temp_comm_list=comm_list[:]
		G.remove_edge(u,v)
		queue.append(u)
		visited[queue[0]]=0
		while len(queue)>0:
			c=comm_node(G,queue[0],temp_comm_list)
			evaluated=[]
			for vis in visited.keys():
				if visited[vis]==0 and vis!=v:
					temp_comm_list=[]
					temp_comm_list.append(vis)
					visited[vis]=1
					evaluated.append(vis)
			que_add=[]
			for q in queue:
				if q not in evaluated:
					neigh=G.neighbors(q)
					removal=[]
					for n in neigh:
						if (comm_node(G,n,comm_list)!=c) or (n in visited) or (n==v):
							removal.append(n)
					neigh=list(set(neigh)-set(removal))
					visited.update({key: 0 for key in neigh})
					que_add+=neigh
			queue=que_add
		#check if 'v' moves
		visited={}
		queue=[]
		temp_comm_list2=comm_list[:]
		queue.append(v)
		visited[queue[0]]=0
		while len(queue)>0:
			c=comm_node(G,queue[0],temp_comm_list2)
			evaluated=[]
			for vis in visited.keys():
				if visited[vis]==0 and vis!=u:
					temp_comm_list2=[]
					temp_comm_list2.append(vis)
					visited[vis]=1
					evaluated.append(vis)

			que_add=[]
			for q in queue:
				if q not in evaluated:
					neigh=G.neighbors(q)
					removal=[]
					for n in neigh:
						if (comm_node(G,n,comm_list)!=c) or (n in visited) or (n==u):
							removal.append(n)
					neigh=list(set(neigh)-set(removal))
					visited.update({key: 0 for key in neigh})
					que_add+=neigh
			queue=que_add
		temp_comm_list=[temp_comm_list]
		temp_comm_list2=[temp_comm_list2]
		if perm_comm(G,0,temp_comm_list2)+perm_comm(G,0,temp_comm_list)> perm_comm(G,c,comm_list):
			comm_list.pop(c)
			comm_list.append(temp_comm_list[0])
			comm_list.append(temp_comm_list2[0])

	elif G.has_edge(u,v):
		G.remove_edge(u,v)
	return G, comm_list


def node_deletion(G,u,comm_list):					#tested, correct
	neighbors_u=G.neighbors(u)
	for v in neighbors_u:
		G,comm_list=edge_deletion(G,u,v,comm_list)
	c=comm_node(G,u,comm_list)
	comm_list=comm_list[0:c]+comm_list[(c+1):]
	G.remove_node(u)
	return G,comm_list

def node_addition(G,u,comm_list):					#tested, correct
	neighbors_u=G.neighbors(u)
	for v in neighbors_u:
		G,comm_list=edge_addition(G,u,v,comm_list)
	c=comm_node(G,u,comm_list)
	comm_list=comm_list[0:c]+comm_list[(c+1):]
	G.add_node(u)
	return G,comm_list

def nmi(comm_new,comm_list):
	actual=[]
	predicted=[]
	for i in range(0,len(comm_new)):
		for j in range(0,len(comm_new[i])):
			actual.append(i)
	for i in range(0,len(comm_list)):
		for j in range(0,len(comm_list[i])):
			predicted.append(i)
	return normalized_mutual_info_score(actual, predicted)
def ari(comm_new,comm_list):
	actual=[]
	predicted=[]
	for i in range(0,len(comm_new)):
		for j in range(0,len(comm_new[i])):
			actual.append(i)
	for i in range(0,len(comm_list)):
		for j in range(0,len(comm_list[i])):
			predicted.append(i)
	return adjusted_rand_score(actual, predicted)

def str_to_int(x):
	return [[int(v) for v in line.split()] for line in x]
#----------main-----------------
edges_added = []
edges_removed = []
nodes_added = []
nodes_removed = []

edge_file='switch.t01.edges'
with open(edge_file,'r') as f:
	edge_list=f.readlines()
	for edge in edge_list:
		edge=edge.split()
		G.add_node(int(edge[0]))
		G.add_node(int(edge[1]))
		G.add_edge(int(edge[0]),int(edge[1]))
G=G.to_undirected()
comm_file='switch.t01.comm'
with open(comm_file,'r') as f:
	comm_list=f.readlines()
	comm_list=str_to_int(comm_list)

start=time.time()
for i in range(2,21):
	comm_new_file=open('output_new_'+str(i)+'.txt','r')
	if i<10:
		edge_list_old_file=open('switch.t0'+str(i-1)+'.edges','r')
		edge_list_old=edge_list_old_file.readlines()
		edge_list_new_file=open('switch.t0'+str(i)+'.edges','r')
		edge_list_new=edge_list_new_file.readlines()
		comm_new=comm_new_file.readlines()
	elif i==10:
		edge_list_old_file=open('switch.t09.edges','r')
		edge_list_old=edge_list_old_file.readlines()
		edge_list_new_file=open('switch.t10.edges','r')
		edge_list_new=edge_list_new_file.readlines()
		comm_new=comm_new_file.readlines()
	else:
		edge_list_old_file=open('switch.t'+str(i-1)+'.edges','r')
		edge_list_old=edge_list_old_file.readlines()
		edge_list_new_file=open('switch.t'+str(i)+'.edges','r')
		edge_list_new=edge_list_new_file.readlines()
		comm_new=comm_new_file.readlines()
	comm_new=str_to_int(comm_new)
	current_snapshot = defaultdict(set)
	previous_snapshot = defaultdict(set)

	previous_nodes = set()
	current_nodes = set()
	for line in edge_list_old:
 		temp = line.strip().split()
 		temp = map(int,temp)
 		previous_snapshot[temp[0]].add(temp[1])
 		previous_nodes.add(temp[0])
 		previous_nodes.add(temp[1])
 	
 	for line in edge_list_new:
 		temp = line.strip().split()
 		temp = map(int,temp)

 		current_snapshot[temp[0]].add(temp[1])
 		current_nodes.add(temp[0])
 		current_nodes.add(temp[1])
 	

 	total_nodes = previous_nodes.union(current_nodes)

 	for n in nodes_added:
 		gr,com_lis = node_addition(G,n,comm_list)

 	for n in nodes_removed:
 		gr, com_lis = node_deletion(G,n,comm_list)


 	for nodes in total_nodes:
 		edgeinprevious = previous_snapshot[nodes]
 		edgeincurrent = current_snapshot[nodes]

 		edgesadded = edgeincurrent - edgeinprevious
 		edgesremoved = edgeinprevious - edgeincurrent



 		for e in edgesadded:
 			edges_added.append([nodes,e])

 		for e in edgesremoved:
 			edges_removed.append([nodes,e])
	for edge in edges_added:
		if edge[0]!=edge[1]:
			G,comm_list=edge_addition(G,edge[0],edge[1],comm_list)
	for edge in edges_removed:
		if edge[0]!=edge[1]:
			G,comm_list=edge_deletion(G,edge[0],edge[1],comm_list)
	actual=[]
	baseline=[]
	for j in range(len(comm_new)):
		for c in comm_new[j]:
			flag=False
			for k in range(len(comm_list)):
				if c in comm_list[k] and flag==False:
					flag=True 
					actual.append(j)
					baseline.append(k)
					break
	print 'nmi', normalized_mutual_info_score(actual, baseline)
	print 'ari', adjusted_rand_score(actual, baseline)
end=time.time()
print 'time', (end-start),'seconds'
