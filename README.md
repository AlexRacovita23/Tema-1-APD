# Tema-1-APD
Student: Racovita Alexandru Catalin
Class: 332 CB

For this project there were 3 main tasks that had to be done
in order to parallelise the Marching Squares Algorithm:
1. Create a data structure array which contained the information needed by the threads, allocate memory for the contour map, 
initial image, rescaled image when used and for the grid map.
2. Declare and initialize the pthreads, split the workload in all arrays fairly between the threads and wait for all threads
to finish before going to the next step of the algorithm. Workload split is explained in code comments.
3. Free all allocated memory, as each image contains at least 12.289KB and a maximum of 110.593 KB for the set of inputs used.
4K images for example can occupy much more memory so freeing it becomes more urgent the higher the size of the image.

Acceleration for the last 2 inputs, using 2 and respectively 4 threads:

Acceleration 1-2 Mappers: 1.78
Acceleration 1-4 Mappers: 3.33

Acceleration 1-2 Mappers: 1.91
Acceleration 1-4 Mappers: 3.72

