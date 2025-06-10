#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define eps 0.001

typedef struct
{
    double *coordinates;
    int dimension;
} vector;

typedef struct
{
    vector *centroid;
    vector *members;
    int num_of_members;
} cluster;

typedef struct
{
    vector *all_vectors;
    int num_vectors;
} all_vecs;

int checkArg(char *str);
int checkConvergence(vector *v1, vector *v2);
void assignVectorToCluster(vector *v, cluster *clus);
void updateCentroid(cluster *clus);
double distance(vector *v1, vector *v2);
vector *sumVectors(vector *vectors, int num_of_vecs);
vector *mulByScalar(vector *v, double scalar);
void emptyCluster(cluster *clus);
cluster *initiateClusters(all_vecs *all_vectors, int num_of_clusters);
cluster *iterateAlgorithm(cluster *cluster_array, all_vecs *all_vectors, int K, int N, int iters);
all_vecs getInput();
void errorHandling();
void printOutput(cluster *clus, int K);
void freeMemory(cluster *clus, all_vecs *all_vectors, int K, int N);
void printVector(vector *vec);

int checkArg(char *str){
    int i;
    int num_of_chars=0;

    while(str[num_of_chars]!='\0'){
        num_of_chars++;
    }

    for(i=0;i<num_of_chars;i++){
        if((((int)str[i])<48 || ((int)str[i])>57)&&(int)str[i]!=46){
            return(0);
        }
    }
    return(1);
}

void printVector(vector *vec)
{
    int i;
    for (i = 0; i < vec->dimension; i++)
    {
        if (i == vec->dimension - 1)
        {
            printf("%.4f", (vec->coordinates)[i]);
        }
        else
            printf("%.4f,", (vec->coordinates)[i]);
    }
    printf("\n");
}

int checkConvergence(vector *v1, vector *v2)
{
    return distance(v1, v2) < eps;
}

void assignVectorToCluster(vector *v, cluster *clus)
{
    clus->members = (vector *)realloc(clus->members, ((clus->num_of_members) + 1) * sizeof(vector));
    clus->num_of_members++;
    if (clus->members == NULL)
    {
        errorHandling();
    }
    clus->members[clus->num_of_members - 1] = *v;
}

void updateCentroid(cluster *clus)
{
    vector *old_centroid = clus->centroid;
    double scalar;
    vector *sum_vector;
    scalar = (double)1 / (clus->num_of_members);
    sum_vector = sumVectors(clus->members, clus->num_of_members);
    clus->centroid = mulByScalar(sum_vector, scalar);
    free(sum_vector->coordinates);
    free(sum_vector);
    free(old_centroid->coordinates);
    free(old_centroid);
}

double distance(vector *v1, vector *v2)
{
    double sum = 0;
    int i;
    for (i = 0; i < v1->dimension; i++)
    {
        sum += pow(v1->coordinates[i] - v2->coordinates[i], 2);
    }
    return sqrt(sum);
}

vector *sumVectors(vector *vectors, int num_of_vecs)
{
    int i;
    vector *sum_vec = (vector *)malloc(sizeof(vector));
    sum_vec->dimension = vectors[0].dimension;
    sum_vec->coordinates = (double *)calloc(sum_vec->dimension, sizeof(double));
    if (sum_vec == NULL)
    {
        errorHandling();
    }
    for (i = 0; i < sum_vec->dimension; i++)
    {
        int j;
        for (j = 0; j < num_of_vecs; j++)
        {
            sum_vec->coordinates[i] += vectors[j].coordinates[i];
        }
    }
    return sum_vec;
}

vector *mulByScalar(vector *v, double scalar)
{
    vector *mul_vec = (vector *)malloc(sizeof(vector));
    int i;
    mul_vec->dimension = v->dimension;
    mul_vec->coordinates = (double *)calloc(v->dimension, sizeof(double));
    if (mul_vec == NULL)
    {
        errorHandling();
    }
    for (i = 0; i < mul_vec->dimension; i++)
    {
        mul_vec->coordinates[i] = v->coordinates[i] * scalar;
    }
    return mul_vec;
}

void emptyCluster(cluster *clus)
{
    free(clus->members);
    clus->members = (vector *)malloc(sizeof(vector));
    clus->num_of_members = 0;
}

cluster *initiateClusters(all_vecs *all_vectors, int K)
{
    int i;
    cluster *cluster_array = (cluster *)malloc(K * sizeof(cluster));
    if (cluster_array == NULL)
    {
        errorHandling();
    }
    for (i = 0; i < K; i++)
    {
        int j;
        cluster_array[i].centroid = (vector *)malloc(sizeof(vector));
        cluster_array[i].centroid->dimension = all_vectors->all_vectors[i].dimension;
        cluster_array[i].centroid->coordinates = (double *)malloc(sizeof(double) * cluster_array[i].centroid->dimension);
        for (j = 0; j < cluster_array[i].centroid->dimension; j++)
        {
            cluster_array[i].centroid->coordinates[j] = all_vectors->all_vectors[i].coordinates[j];
        }
        cluster_array[i].num_of_members = 0;
        cluster_array[i].members = (vector *)malloc(sizeof(vector));
    }
    return cluster_array;
}
cluster *iterateAlgorithm(cluster *cluster_array, all_vecs *all_vectors, int K, int N, int iter)
{
    int i;
    for (i = 0; i < iter; i++)
    {
        int convergence_flag = 0;
        int j;
        for (j = 0; j < N; j++)
        {
            double min_dist = distance(&(all_vectors->all_vectors[j]), cluster_array[0].centroid);
            cluster *potencial_cluster = &(cluster_array[0]);
            int k;
            for (k = 1; k < K; k++)
            {
                double new_dist = distance(&(all_vectors->all_vectors[j]), cluster_array[k].centroid);
                if (new_dist < min_dist)
                {
                    min_dist = new_dist;
                    potencial_cluster = &(cluster_array[k]);
                }
            }
            assignVectorToCluster(&(all_vectors->all_vectors[j]), potencial_cluster);
        }
        for (j = 0; j < K; j++)
        {
            vector *old_centroid_copy = (vector *)malloc(sizeof(vector));
            int l;
            old_centroid_copy->dimension = cluster_array[j].centroid->dimension;
            old_centroid_copy->coordinates = (double *)malloc(sizeof(double) * old_centroid_copy->dimension);
            for (l = 0; l < old_centroid_copy->dimension; l++)
            {
                old_centroid_copy->coordinates[l] = cluster_array[j].centroid->coordinates[l];
            }
            updateCentroid(&(cluster_array[j]));
            convergence_flag += checkConvergence(old_centroid_copy, cluster_array[j].centroid);
            emptyCluster(&cluster_array[j]);
            free(old_centroid_copy->coordinates);
            free(old_centroid_copy);
        }
        if (convergence_flag == K)
            break;
    }
    return cluster_array;
}

all_vecs getInput()
{
    double n;
    char c;
    int i = 0, j = 0;
    all_vecs all_vectors;
    vector curr_vector;
    curr_vector.dimension = 0;
    curr_vector.coordinates = (double *)malloc(sizeof(double));
    if (curr_vector.coordinates == NULL)
    {
        errorHandling();
    }
    all_vectors.all_vectors = (vector *)malloc(sizeof(vector));
    if (all_vectors.all_vectors == NULL)
    {
        errorHandling();
    }
    while (scanf("%lf%c", &n, &c) == 2)
    {
        if (c == '\n')
        {
            vector new_vector;
            curr_vector.coordinates[j] = n;
            j++;
            if (i == 0)
                curr_vector.dimension++;
            all_vectors.all_vectors[i] = curr_vector;
            i++;
            all_vectors.all_vectors = (vector *)realloc(all_vectors.all_vectors, sizeof(vector) * (i + 1));
            if (all_vectors.all_vectors == NULL)
            {
                errorHandling();
            }
            new_vector.dimension = j;
            new_vector.coordinates = (double *)malloc(sizeof(double) * new_vector.dimension);
            if (new_vector.coordinates == NULL)
            {
                errorHandling();
            }
            curr_vector = new_vector;
            j = 0;
            continue;
        }
        curr_vector.coordinates[j] = n;
        j++;
        if (i == 0)
        {
            curr_vector.dimension++;
            curr_vector.coordinates = (double *)realloc(curr_vector.coordinates, sizeof(double) * (j + 1));
            if (curr_vector.coordinates == NULL)
            {
                errorHandling();
            }
        }
    }
    free(curr_vector.coordinates);
    all_vectors.num_vectors = i;
    return all_vectors;
}

void errorHandling()
{
    printf("An Error Has Occured\n");
}

void printOutput(cluster *clus, int K)
{
    int i;
    for (i = 0; i < K; i++)
    {
        printVector(clus[i].centroid);
    }
}

void freeMemory(cluster *cluster_array, all_vecs *all_vectors, int K, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        free(all_vectors->all_vectors[i].coordinates);
    }
    free(all_vectors->all_vectors);
    for (i = 0; i < K; i++)
    {
        free(cluster_array[i].members);
        free(cluster_array[i].centroid->coordinates);
        free(cluster_array[i].centroid);
    }
    free(cluster_array);
}

int main(int argc, char **argv)
{
    int K;
    double K_f;
    int iter = 400;
    double iter_f = 400;
    all_vecs all_vectors;
    cluster *cluster_array;
    int N;
    int iter_verification = 1;

    if (argc > 3)
    {
        errorHandling();
        return(1);
    }

    all_vectors = getInput();
    N = all_vectors.num_vectors;
    K = atoi(argv[1]);
    K_f = atof(argv[1]);
    if (K != K_f || !(K > 1 && K < N) || checkArg(argv[1])==0)
    {
        printf("Incorrect number of clusters!\n");
        return (1);
    }

    if (argc == 3)
    {
        iter = atoi(argv[2]);
        iter_f = atof(argv[2]);
        iter_verification=checkArg(argv[2]);
    }


    if (iter != iter_f || !(iter > 1 && iter < 1000) || iter_verification==0)
    {
        printf("Incorrect maximum iteration!\n");
        return (1);
    }
    cluster_array = initiateClusters(&all_vectors, K);
    cluster_array = iterateAlgorithm(cluster_array, &all_vectors, K, N, iter);
    printOutput(cluster_array, K);
    freeMemory(cluster_array, &all_vectors, K, N);
    return (0);
}
