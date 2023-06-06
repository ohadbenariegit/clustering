#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../common/helper_methods.h"

/*--------------------------------------------------------------------------------------------------------------------*/
#define MAX_ITER 300
#define EPS 0

/*--------------------------------------------------------------------------------------------------------------------*/
double get_squared_distance(int k, double *v1, double *v2) {
    double dist = 0;
    int i;

    for (i = 0; i < k; i++)
        dist += pow(v1[i] - v2[i], 2);

    return dist;
}

/* return index of centroid whose squared euclidean distance from datapoint is minimal */
int get_closest_centroid(int k, double *datapoint, double **centroids) {
    double dist, min_dist = HUGE_VAL;
    int cent, closest = -1;

    for (cent = 0; cent < k; cent++) {
        dist = get_squared_distance(k, datapoint, centroids[cent]);

        if (dist < min_dist) {
            min_dist = dist;
            closest = cent;
        }
    }

    return closest;
}


void set_centroids(int n, int k, double **datapoints, double **centroids) {
    double **new_centroids;
    int *sizes;
    double delta, max_delta = 0;
    int closest, count = 0;
    int i, j;

    new_centroids = create_matrix(k, k);

    sizes = (int *) create_array(k, sizeof(int));

    do {
        /* reset new centroids, cluster sizes, and deltas */
        for (i = 0; i < k; i++)
            memset(new_centroids[i], 0, k * sizeof(double));
        
        memset(sizes, 0, k * sizeof(int));

        /* calculate size and vector sum of clusters */
        for (i = 0; i < n; i++) {
            closest = get_closest_centroid(k, datapoints[i], centroids);

            for (j = 0; j < k; j++)
                new_centroids[closest][j] += datapoints[i][j];

            sizes[closest]++;
        }

        /* calculate cluster means, maximum delta, and set centroids to cluster means */
        for (i = 0; i < k; i++) {
            if (sizes[i] != 0) {
                for (j = 0; j < k; j++)
                    new_centroids[i][j] /= sizes[i];
            }

            delta = get_squared_distance(k, new_centroids[i], centroids[i]);

            if (delta > max_delta)
                max_delta = delta;

            memcpy(centroids[i], new_centroids[i], k * sizeof(double));
        }

        count++;
    } while (count < MAX_ITER && max_delta > EPS);

    free_matrix(k, new_centroids);
    free(sizes);
}