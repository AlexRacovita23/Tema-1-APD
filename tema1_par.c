// Author: APD team, except where source was noted
//#define _POSIX_C_SOURCE 200112L
//src

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

typedef struct {
    int ID; //thread ID
    int P; //number of threads
    ppm_image *image;
    ppm_image **contour_map;
    ppm_image *scaled_image;
    unsigned char **grid;
    pthread_barrier_t *barrier;
} thread_data;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

int min(int x, int y) {
    if (x < y) {
        return x;
    } else return y;
}

void update_image_parallel(ppm_image *image, ppm_image *contour, int x, int y) {

    // Nothing was changed here, also parallelising this function would make it so only some "stripes"
    // of the image would be updated, which would result in a mashup between the initial and the updated image
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

ppm_image *rescale_image_parallel(ppm_image *image, ppm_image *new_image, int ID, int P) {
    uint8_t sample[3];

    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        return new_image;
    }

    int start = ID * (double) new_image->x / P;
    int end = min((ID + 1) * (double)new_image->x / P, new_image->x);
    // use bicubic interpolation for scaling
    // We can parallelize the rescaling segment of the algorithm by dividing the image into horizontal strips
    // The image fragments are independent from each other 
    for (int i = start; i < end; i++) {
        for (int j = 0; j < new_image->y; j++) {
            float u = (float)i / (float)(new_image->x - 1);
            float v = (float)j / (float)(new_image->y - 1);
            sample_bicubic(image, u, v, sample);
            new_image->data[i * new_image->y + j].red = sample[0];
            new_image->data[i * new_image->y + j].green = sample[1];
            new_image->data[i * new_image->y + j].blue = sample[2];
        }
    }

    return new_image;
}

unsigned char **sample_grid_parallel(ppm_image *image, unsigned char **grid, int ID, int P) {
    int p = image->x / STEP;
    int q = image->y / STEP;
    
    int start_x = ID * (double) p / P;
    int end_x = min((ID + 1) * (double)p / P, p);

    // When sampling, the squares from which we sample are independent from each other so we can parallelize this process
    for (int i = start_x; i < end_x; i++) {
        for (int j = 0; j < q; j++) {
            
            ppm_pixel curr_pixel = image->data[i * STEP * image->y + j * STEP];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;
            if (curr_color > SIGMA) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    // Even though these "for" loops are not as intensive we can still divide the last line and column between the threads
    for (int i = start_x; i < end_x; i++) {
        ppm_pixel curr_pixel = image->data[i * STEP * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }
    
    int start_y = ID * (double) q / P;
    int end_y = min((ID + 1) * (double)q / P, q);

    for (int j = start_y; j < end_y; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * STEP];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;
        
        if (curr_color > SIGMA) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }
    return grid;
}

void march_parallel(ppm_image *image, unsigned char **grid, ppm_image **contour_map, int ID, int P) {
    int p = image->x / STEP;
    int q = image->y / STEP;
    int start = ID * (double) p / P;
    int end = min((ID + 1) * (double)p / P, p);

    // We can parallelize the marching segment of the algorithm by dividing the image into horizontal strips
    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image_parallel(image, contour_map[k], i * STEP, j * STEP);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

void *thread_function(void *arg) {
    thread_data* data = (thread_data*) arg;
    int ID = data->ID;
    int P = data->P;
    ppm_image **contour_map = data->contour_map;
    ppm_image *image = data->image;
    ppm_image *scaled_image = data->scaled_image;
    unsigned char **grid = data->grid;

    scaled_image = rescale_image_parallel(image, scaled_image, ID, P);
    pthread_barrier_wait(data->barrier); // Wait for rescaling to finish before sampling, otherwise unwanted behaviour can happen
    
    grid = sample_grid_parallel(scaled_image, grid, ID, P);
    pthread_barrier_wait(data->barrier);// Wait for sampling to finish before marching, otherwise unwanted behaviour can happen

    march_parallel(scaled_image, grid, contour_map, ID, P);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    int P = atoi((argv[3]));
    int r;
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, P); 

    pthread_t threads[P];
    thread_data *data = (thread_data *)malloc(sizeof(thread_data) * P);
    if (!data) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    ppm_image *image = read_ppm(argv[1]);
    ppm_image **contour_map = init_contour_map();

    // Rescaled image allocation
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // Allocate memory for new image data if rescaling happens
    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        new_image->x = image->x;
        new_image->y = image->y;
        new_image->data = image->data;
    } else {
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;
    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
        if (!new_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // Grid allocation
    unsigned char **grid = (unsigned char **)malloc((new_image->x + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    for (int i = 0; i <= new_image->x; i++) {
        grid[i] = (unsigned char *)malloc((new_image->y + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // Thread initialization
    for(int ID = 0; ID < P; ID++) {
        data[ID].ID = ID; //nothing to free
        data[ID].image = image; // freed later
        data[ID].contour_map = contour_map; //freed later
        data[ID].grid = grid; //freed later
        data[ID].scaled_image = new_image; //freed later
        data[ID].P = P; //nothing to free
        data[ID].barrier = &barrier;
        r = pthread_create(&threads[ID], NULL, thread_function, &data[ID]);
        if(r) {
            fprintf(stderr, "Error creating thread\n");
            exit(1);
        }
    }
    for (int ID = 0; ID < P; ID++) {
		int r = pthread_join(threads[ID], NULL);
		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", ID);
			exit(-1);
		}
	}
    // Write output
    write_ppm(data[0].scaled_image, argv[2]);
    for(int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);
    free(image->data);
    free(image);
    printf("segfault\n");
    if(new_image->data != image->data) {
       free(new_image->data);
    }
    
    for(int i = 0; i <= new_image->x; i++) {
        free(grid[i]);
    }
    free(new_image);
    free(grid);
    free(data);

    return 0;
}

