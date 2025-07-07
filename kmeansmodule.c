# define PY_SSIZE_T_CLEAN
# include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
int checkConvergence(vector *v1, vector *v2, double eps);
void assignVectorToCluster(vector *v, cluster *clus);
void updateCentroid(cluster *clus);
double distance(vector *v1, vector *v2);
vector *sumVectors(vector *vectors, int num_of_vecs);
vector *mulByScalar(vector *v, double scalar);
void emptyCluster(cluster *clus);
cluster *initiateClusters(all_vecs *all_vectors, int num_of_clusters);
cluster *iterateAlgorithm(cluster *cluster_array, all_vecs *all_vectors, int K, int N, int iters, double eps);
//all_vecs getInput();
void errorHandling();
void printOutput(cluster *clus, int K);
void freeMemory(cluster *clus, all_vecs *all_vectors, all_vecs *all_centroids, int K, int N);
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

int checkConvergence(vector *v1, vector *v2, eps)
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
cluster *iterateAlgorithm(cluster *cluster_array, all_vecs *all_vectors, int K, int N, int iter, double eps)
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
            convergence_flag += checkConvergence(old_centroid_copy, cluster_array[j].centroid, eps);
            emptyCluster(&cluster_array[j]);
            free(old_centroid_copy->coordinates);
            free(old_centroid_copy);
        }
        if (convergence_flag == K)
            break;
    }
    return cluster_array;
}

/*
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
*/

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

void freeMemory(cluster *cluster_array, all_vecs *all_vectors, all_vecs *all_centroids, int K, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        free(all_vectors->all_vectors[i].coordinates);
    }
    for (i = 0; i < K; i++)
    {
        free(all_centroids->all_vectors[i].coordinates);
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

static PyObject* fit(PyObject* self, PyObject* args)
{
    PyObject *centroids; // צנטרואידים
    PyObject *points; // כלל הנקודות
    int iter; // מס׳ איטרציות
    double eps; // אפסילון
    all_vecs all_vectors; // רשימה של וקטורים - כל הנקודות - יש מצביע לרשימה של הטיפוס vector ויש אינט של כמה יש ברשימה
    all_vecs all_centroids; // רשימה של וקטורים - כל הצנטרואידים - יש מצביע לרשימה של הטיפוס vector ויש אינט של כמה יש ברשימה
    cluster *cluster_array; // מערך של קלאסטרים שלכל אחד מצביע לצנטרואיד ומצביע לרשימת נקודות של הצנטרואיד

    if (!PyArg_ParseTuple(args, "OOid", &centroids, &points, &iter, &eps)) { // שם את המערכים בתוך פייאובג׳ט ואת המספרים לפי טיפוסם
        errorHandling(); // תדפיס שגיאה אם ההקצאה נכשלה
        return NULL; // הפייתון כבר יצא עם אקסיט קוד 1
    }

    int K = PyArray_DIM((PyArrayObject*)centroids, 0); // כמה שורות יש במערך של הצנטרואידים = כמה צנטרודים/קלאסטרים יש = K
    int dim = PyArray_DIM((PyArrayObject*)centroids, 1); // כמה עמודות יש במערך של הצנטרואידים = מימד = dim
    int N = PyArray_DIM((PyArrayObject*)points, 0); // כמה נקודות יש בסה״כ = N

    // טיפול בנקודות
    all_vectors.num_vectors = N; // איתחול המס׳ של הווקטורים/נקודות בעצם שמייצג אותם
    all_vectors.all_vectors = (vector *)malloc(sizeof(vector) * N); // מקצים מקום ל-N וקטורים שזה בעצם נקודות
    if (all_vectors.all_vectors == NULL) { // אם ההקצאה נכשלה
        errorHandling(); // תדפיס שגיאה
        return NULL; // הפייתון כבר יצא עם אקסיט קוד 1
    }

    double *points_data = (double *)PyArray_DATA((PyArrayObject*)points); // שם את הנקודות אחת אחרי השניה במערך של דאבלים
    for (int i = 0; i < N; i++) { // שורה-שורה
        all_vectors.all_vectors[i].dimension = dim; // המימד של הנקודה ה-i נקבע להיות דים
        all_vectors.all_vectors[i].coordinates = (double *)malloc(sizeof(double) * dim); // מקצים מקום ל-dim דאבלים שיהיו הקורדינטות של הנקודה ה-i
        if (all_vectors.all_vectors[i].coordinates == NULL) { // אם ההקצאה נכשלה
            errorHandling(); // תדפיס שגיאה
            return NULL; // הפייתון כבר יצא עם אקסיט קוד 1
        }
        for (int j = 0; j < dim; j++) { // עמודה-עמודה
            all_vectors.all_vectors[i].coordinates[j] = points_data[i * dim + j]; // בוקטור האיי בקורדינטה הג׳יי נשים את הדאבל מספר שורה*מימד + עמודה במערך דאבלים
        }
    }

    // טיפול בצנטרואידים
    all_centroids.num_vectors = K; // איתחול של המס׳ של הצנטרואידים בעצם שמייצג אותם
    all_centroids.all_vectors = (vector *)malloc(sizeof(vector) * K); // מקצים מקום ל-K צנטרואידים
    if (all_centroids.all_vectors == NULL) { // אם ההקצאה נכשלת
        errorHandling(); // תדפיס שגיאה
        return NULL; // הפייתון כבר יצא עם אקסיט קוד 1
    }

    double *centroids_data = (double *)PyArray_DATA((PyArrayObject*)centroids); // שם את הצנטרואידים אחד אחרי השני במערך של דאבלים
    for (int i = 0; i < K; i++){ // שורה-שורה
        all_centroids.all_vectors[i].dimension = dim; // נקבע את המימד של הצנטרואיד האיי להיות דים
        all_centroids.all_vectors[i].coordinates = (double *)malloc(sizeof(double) * dim); מקצים מקום לדים דאבלים שיהיו הקורדינטות של הצנטרואיד האיי
        if (all_centroids.all_vectors[i].coordinates == NULL) { // אם ההקצאה נכשלה
            errorHandling(); // תדפיס שגיאה
            return NULL; // הפייתון כבר יצא עם אקסיט קוד 1
        }
        for (int j = 0; j < dim; j++) { // עמודה-עמודה
            all_centroids.all_vectors[i].coordinates[j] = centroids_data[i * dim + j]; // בצנטרואיד האיי בקורדינטה הג׳יי נשים את הדאבל מספר שורה*מימד + עמודה במערך דאבלים
        }
    }


    cluster_array = initiateClusters(&all_centroids, K);
    cluster_array = iterateAlgorithm(cluster_array, &all_vectors, K, N, iter, eps);


    PyObject *result = PyList_New(K); // יוצר pylist object שזה סוג של pyobject
    for (int i = 0; i < K; i++) { // לכל צנטרואיד
        PyObject *cent = PyList_New(dim); // יוצר רשימה מהמימד של הצנטרואיד
        for (int j = 0; j < dim; j++) { // לכל נקודה בצנטרואיד
            PyList_SetItem(cent, j, PyFloat_FromDouble(cluster_array[i].centroid->coordinates[j])); // מכניסה את הקורדינטה הג׳יי של הצנטרואיד האיי למקום הג׳יי בצנטרואיד בפייתון
        }
        PyList_SetItem(result, i, cent); // בתוצאה בשורה האיי אני שמה את הסנטרואיד שיצרתי
    }

    freeMemory(cluster_array, &all_vectors, &all_centroids, K, N);

    return result; // מחזיר לפייתון רשימה של רשימות
}



static PyMethodDef kmeansMethods[] = {
    {"fit",
     (PyCFunction) fit,
     METH_VARARGS,
     PyDoc_STR("run kmeans on centroids")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "kmeansmodule",
    NULL,
    -1,
    kmeansMethods
};

PyMODINIT_FUNC PyInit_kmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}

/* int main(int argc, char **argv)
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
} */
