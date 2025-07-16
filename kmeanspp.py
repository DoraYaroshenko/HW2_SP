import math
import sys
import numpy as np
import pandas as pd
import mykmeanspp

def euclidian_distance(v1, v2):
    sum = 0
    for i in range(len(v1)):
       sub = v1[i] - v2[i]
       sum += sub ** 2
    return math.sqrt(sum)

def checkInput(K,iter,eps):
    if not K.is_integer():
        print("Invalid number of clusters!")
        sys.exit(1)
    K=int(K)   
    if(K <= 1): 
        print("Invalid number of clusters!")
        sys.exit(1)
    if iter!=300 and not iter.is_integer():
        print("Invalid maximum iteration!")
        sys.exit(1)
    iter=int(iter)
    if(iter>=1000 or iter<=1):
        print("Invalid maximum iteration!")
        sys.exit(1)
    if(eps<0):
        print("Invalid epsilon!")
        sys.exit(1)
    return K,iter

def getInputVariables():
    if len(sys.argv)>6 or len(sys.argv)<5:
        print("An Error Has Occurred")
        sys.exit(1)
    K=-1
    iter=300
    eps=0
    file_name_1 = ""
    file_name_2 = ""
    try:
        K = float(sys.argv[1])
    except:
        print("Invalid number of clusters!")
        sys.exit(1)
    iter_flag = len(sys.argv) > 5
    i=2
    if iter_flag:
        i+=1
        try:
            iter = float(sys.argv[2])
        except:
            print("Invalid maximum iteration!")
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
    columns_names1 = [f"col {i}" for i in range(len(vectors1_df.columns))]
    columns_names2 = [f"col {i}" for i in range(len(vectors2_df.columns))]
    columns_names1[0] = "key"
    columns_names2[0] = "key"
    vectors1_df.columns=columns_names1
    vectors2_df.columns=columns_names2
    vectors = vectors1_df.merge(vectors2_df,how='inner', on='key')
    vectors = vectors.set_index('key')
    vectors.sort_values(by='key', ascending=True, inplace=True)
    columns_names = [f"coordinate {i}" for i in range(len(vectors.columns))]
    vectors.columns=columns_names
    return vectors

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
    for i in range(1,K):
        distances = []
        for p in points:
            point_distances = [euclidian_distance(centroid,p) for centroid in centroids[:i]]
            dist = min(point_distances)
            distances.append(dist)
        probabilities=[distance/sum(distances) for distance in distances]
        index_chosen=np.random.choice(len(points),p=probabilities)
        centroids[i]=points[index_chosen]
        centroids_indexes.append(vectors.index[np.where((points_array==centroids[i]).all(axis=1))[0][0]])
        points = np.delete(points,index_chosen, axis=0)
    return centroids,centroids_indexes

        

if __name__ == "__main__":
    K,iter,eps,file_name_1,file_name_2 = getInputVariables()
    K,iter = checkInput(K,iter,eps)
    vectors = getObservations(file_name_1,file_name_2)
    N = len(vectors)
    if(K >= N): 
        print("Invalid number of clusters!")
        sys.exit(1)
    points_array=vectors.to_numpy().tolist()
    centroids,centroids_indexes = initCentroids(vectors,points_array, K)
    print(",".join(str(int(v)) for v in centroids_indexes))
    centroids = centroids.tolist()
    final_centroids = mykmeanspp.fit(centroids, points_array, iter, eps)
    if final_centroids is None:
        print("An Error Has Occurred")
        sys.exit(1)
    for centroid in final_centroids:
        print(",".join('{:.4f}'.format(coordinate) for coordinate in centroid))