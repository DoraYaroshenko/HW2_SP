import math
import sys
import numpy as np
import pandas as pd
import mykmeanssp

def euclidian_distance(v1, v2):
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

def checkInput(K,iter,eps,N):   
    if(K <= 1 or K>=N): 
        print("Incorrect number of clusters!")
        sys.exit(1)
    if(iter>=1000 or iter<=1):
            print("Incorrect maximum iteration!")
            sys.exit(1)
    if(eps<0):
        print("Invalid epsilon!")
        sys.exit(1)

def getInputVariables():
    if len(sys.argv)>6 or len(sys.argv)<5:
        print("An Error has Occured")
        sys.exit(1)
    K=-1
    iter=300
    eps=0
    file_name_1 = ""
    file_name_2 = ""
    try:
        K = float(sys.argv[1])
        if not K.is_integer():
            raise ValueError
        K = int(K)
    except:
        print("Incorrect number of clusters!")
        sys.exit(1)
    iter_flag = len(sys.argv) > 5
    i=2
    if iter_flag:
        i+=1
        try:
            iter = float(sys.argv[2])
            if not iter.is_integer():
                raise ValueError
            iter = int(iter)
        except:
            print("Incorrect maximum iteration!")
            sys.exit(1)
    try:
        eps = float(sys.argv[i])
    except:
        print("Invalid epsilon!")
        sys.exit(1)
    file_name_1 = sys.argv[i+1]
    file_name_2 = sys.argv[i+2]
    return K,iter,eps,file_name_1,file_name_2

def getObservations(file_name_1,file_name_2):
    vectors1_df = pd.read_csv(file_name_1,header=None)
    vectors2_df = pd.read_csv(file_name_2,header=None)
    columns_names = [f"col {i}" for i in range(len(vectors1_df.columns))]
    columns_names[0] = "key"
    vectors1_df.columns=columns_names
    vectors2_df.columns=columns_names
    vectors = vectors1_df.merge(vectors2_df,how='left', on='key')
    vectors = vectors.set_index('key')
    vectors.sort_values(by='key', ascending=True, inplace=True)
    columns_names = [f"coordinate {i}" for i in range(len(vectors.columns))]
    vectors.columns=columns_names
    return vectors

# def dist_from_nearest(point, centroids):


def initCentroids(vectors,points_array,K):
    points = np.copy(points_array)
    np.random.seed(1234)
    index_chosen = np.random.choice(len(points_array))
    centroid = points[index_chosen]
    centroids = np.zeros((K,len(points[0])))
    centroids_indexes = []
    centroids[0]=centroid
    centroids_indexes.append(vectors.index[index_chosen])
    points = np.delete(points,index_chosen,axis=0)
    # print(points)
    for i in range(1,K):
        distances = []
        for p in points:
            # print(p)
            # print(centroids)
            point_distances = [euclidian_distance(centroid,p) for centroid in centroids[:i]]
            dist = min(point_distances)
            distances.append(dist)
        probabilities=[distance/sum(distances) for distance in distances]
        index_chosen=np.random.choice(len(points),p=probabilities)
        centroids[i]=points[index_chosen]
        centroids_indexes.append(vectors.index[np.where((points_array==centroids[i]).all(axis=1))[0][0]])
        points = np.delete(points,index_chosen, axis=0)
    return centroids,centroids_indexes

        

if _name_ == "_main_":
    K,iter,eps,file_name_1,file_name_2 = getInputVariables()
    vectors = getObservations(file_name_1,file_name_2)
    # print(vectors)
    N = len(vectors)
    checkInput(K,iter,eps,N)
    points_array=vectors.to_numpy()
    centroids,centroids_indexes = initCentroids(vectors,points_array, K)
    print(",".join(str(int(v)) for v in centroids_indexes))
    final_centroids = mykmeanssp.fit(centroids, points_array, iter, eps)