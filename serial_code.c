/*
Serial code
Authors: Achilleas Grammenos 1312
         Zelios Andreas 1326
*/
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include </mnt/nas/bmp/bmpfile.h>

#define IMAGE_MAX_X 3142
#define IMAGE_MAX_Y 2254

#define LAMBA_NAN 500
#define nrofbytes_per_stored_element 4

#define COORDINATES_MAX 100

//something negative to get out of bounds, in
//calculations we get abs of value to calculate
#define HEIGHT_BOUND -200

typedef struct coordinates
{
    float x;
    float y;
} coordinates;

//transform degrees to radians
float radian_from_degree(float degree)
{

	return degree*(M_PI/180.0);
}

//finds direction factor
float lambda_line(float degree)
{

	if( (int)degree == 90 )
	{

		return LAMBA_NAN;
	}

	return tan(radian_from_degree(degree));
}

// t / f if points are in the image returns TRUE
bool valid_point_inside_picture(int x ,int y)
{
	if( x <= 0 || x >= IMAGE_MAX_X || y <= 0 || y >= IMAGE_MAX_Y )
		return FALSE;

	else
		return TRUE;
}

//we get data matrix, x ,y and get height, if out of bounds return a large number
float getheight(void *data0 ,int x ,int y)
{
	float (*data)[IMAGE_MAX_X] = data0;
	float height;

	if( x < 0 || x > IMAGE_MAX_X || y < 0 || y > IMAGE_MAX_Y )
		height = (float) HEIGHT_BOUND;

	else
		height = data[y][x];

	return height;
}

//Find out the quarter
void compute_direction(float degrees ,int *x_direction ,int *y_direction)
{
	//reduction to 0-360 degrees, in case degrees > 360
	float degree = fmodf(degrees ,360);

	//first quarter
	if( degree > 0 && degree <= 90)
	{
		*x_direction = 1;
		*y_direction = 1;
	}
	else if( degree > 90 && degree <= 180) //second quarter
	{
		*x_direction = -1;
		*y_direction = 1;
	}
	else if( degree > 180 && degree <= 270)//third quarter
	{
		*x_direction = -1;
		*y_direction = -1;
	}
	else if( degree > 270 && degree <= 360)//fourth quarter
	{
		*x_direction = 1;
		*y_direction = -1;
	}
	else //in this case the function stops
	{
		*x_direction = 900;
		*y_direction = 900;
	}
}

//Initial point y0,x0 with degree
//old x,y values : x_old ,y_old
//calculate height y/x
//and ymax ,xmax of heighest point
float compute_height_at_x_y(void *data0 ,float currentheight ,float degree ,int x0 ,int y0 ,float x ,float y ,float x_old ,float y_old ,int *xmax ,int *ymax)
{
	int x_low ,x_high, i;
	int y_low ,y_high;
	int x_old_low ,y_old_low;
	float height = 0.0;
	float lambda;
	int x_direction ,y_direction;

	//Find the quarter
	compute_direction(degree ,&x_direction ,&y_direction);

	//lambda is the direction factor
	lambda = lambda_line(degree);

	//x_low, x_high, y_low, y_high are aproximate values for x, y
	x_low = (int)floor(x);
	x_high = (int)ceil(x);
	y_low = (int)floor(y);
	y_high = (int)ceil(y);

	x_old_low = (int)floor(x_old);
	y_old_low = (int)ceil(y_old);

	int levelsy;

	//only vertical movement of 2 neighboring points
	if(x_low == x_high && (y_old_low == (y_low + 1) || y_old_low == y_low))
	{
		height = (getheight(data0 ,y_low ,x_low) + getheight(data0 ,y_high ,x_high)) / 2; //interpolation
		*xmax=x_low;
		*ymax=y_low;
	}
	else
	{
		levelsy = abs(y_low - y_old_low);

		for (i = 0; i < levelsy; i++)
		{
			y_low = y_low + i*(y_direction); //move observer
			float h1 = getheight(data0 ,y_low ,x_old_low);
			float h2 = getheight(data0 ,y_low ,x_low);//our current point

			if (h1 == HEIGHT_BOUND || h2 == HEIGHT_BOUND)
			{
				height = HEIGHT_BOUND;

				//if one of the two pixels is in bound
				//choose one of them
				if (h1 == HEIGHT_BOUND && h2 != HEIGHT_BOUND)
					height = h2;

				if (h2 == HEIGHT_BOUND && h1 != HEIGHT_BOUND)
					height = h1;
			}

			else //interpolation
				height = (h1 + h2)/2;

 			y = (float)(h1 + h2)/2;


			//final x and it's low/heigh
			//lamda not valid, so we are vertical
			if (lambda == LAMBA_NAN)
				x = x0;

			else
				x = x0 + (y-y0)/lambda;

			x_low = (int) floor(x);
			x_high = (int) ceil(x);
			*xmax = x_low;
			*ymax = y_low;

			//if current height is greater than x0,y0 stop
			if (currentheight == (float) HEIGHT_BOUND)
				return currentheight;

			//if out of bounds
			if (valid_point_inside_picture(x_low,y_low) == FALSE)
				return currentheight;

			//if the height is lower stop
			//and return xmax/ymax/height
			if (currentheight < height)
				break;

			else
			{
				float h1 = getheight(data0 ,y_low ,x_high);
				float h2 = getheight(data0 ,y_high ,x_high);

				if (h1 == HEIGHT_BOUND && h2 != HEIGHT_BOUND)
					height = h2;

				else if (h2 == HEIGHT_BOUND && h1 != HEIGHT_BOUND)
					height = h1;

				else if (h1 == HEIGHT_BOUND && h2 == HEIGHT_BOUND)
					height = HEIGHT_BOUND;

				else
					height = (h1 + h2)/2;//interpolation
			}
		}
	}

	return height;
}

//Height calculation for each point x/x+1/x+2 and corresponding y
//also draws the Line of sight (LOS) in grayscale 0 - 255
void compute_height_degree(void *data0 ,float degree ,int x0 ,int y0 ,float observer_height ,int *xmax0 ,int *ymax0 ,float *pointheight0,int value,
							bmpfile_t *bmp)
{
	float (*data) [IMAGE_MAX_X] = data0;//transform 2-dimension array to 1-dimension array
	float currentheight;
	char we_are_higher = 1;
	float lambda;
	float x ,y;
	float x_old ,y_old;
	float pointheight;
	int xmax ,ymax;
	int y_direction ,x_direction;
	rgb_pixel_t pixel;
	rgb_pixel_t *pixel_old;

	currentheight = data[y0][x0] + observer_height;

	lambda = lambda_line(degree);

	x = (float)x0;
	y = (float)y0;

	//find the quarter of the next point
	compute_direction(degree ,&x_direction ,&y_direction);

	while (we_are_higher == 1)
	{
		if (valid_point_inside_picture(x ,y) == FALSE)
			break;

		x_old = x;
		y_old = y;

		//for tan90
		if ((int)degree == 90)
		{
			x = x0;
			y = y + y_direction;
		}
		else
		{
			x = x + x_direction;
			y = lambda * (x-x0) + y0;//line equation

		}

		if (valid_point_inside_picture(x ,y) == FALSE)
			break;

		//find height at (x,y)
		pointheight = compute_height_at_x_y(data0 ,currentheight ,degree ,x0 ,y0 ,x ,y ,x_old ,y_old ,&xmax ,&ymax);

		if (pointheight > currentheight){//if observer can't see

			we_are_higher = 0;//to stop while

			if (xmax < 0)
			{
				xmax = 0;
				pointheight = -10;//arbitary value
			}

			if (xmax >= IMAGE_MAX_X)
			{
				xmax = IMAGE_MAX_X - 1;
				pointheight = -10;
			}

			if(ymax < 0)
			{
				ymax = 0;
				pointheight = -10;
			}

			if (ymax >= IMAGE_MAX_Y)
			{
				ymax = IMAGE_MAX_Y - 1;
				pointheight = -10;
			}

			*ymax0 = ymax;
			*xmax0 = xmax;
			*pointheight0 = pointheight;
		}
		else{

			if (valid_point_inside_picture(x ,y) == FALSE)
                        	break;

			pixel_old = bmp_get_pixel(bmp , x, y);

			//in case we go above 255
			if((*pixel_old).red >= 255 && (*pixel_old).green >= 255 && (*pixel_old).blue >= 255)
			{
				(*pixel_old).red = 255;
                (*pixel_old).green = 255;
                (*pixel_old).blue = 255;
            }
			else{
				(*pixel_old).red += value;
				(*pixel_old).green += value;
				(*pixel_old).blue += value;

				pixel.red = (*pixel_old).red;
				pixel.green = (*pixel_old).green;
				pixel.blue = (*pixel_old).blue;
				bmp_set_pixel(bmp, x, y, pixel);
			}

		}
	}
}

int main(int argc ,char *argv[])
{

	int fd;
	float degree;
	float *data;
	int x0 ,y0;
	struct stat sbuf;
	FILE *fp;
	char *coordinate_file_name;
	int number_of_coordinates;
	bool all_pixels = FALSE;

	char filename[60];
	int ymax0 ,xmax0;
	float pointheight0;
	float degreemax;//max degrees of observers`s rotation
	float degreestep;//step of degree`s rise
	float observer_height;

	if (argc != 7)
		fprintf(stderr ,"Use ./a.out fname observer_height degreeStart degreeEnd degreeStep coordinate_file_name\n\n");

	strcpy(filename ,argv[1]);

	if ((fd = open(filename,O_RDONLY)) == -1)
	{
		perror("open");
		exit(1);
	}

	if (stat(filename ,&sbuf) == -1)
	{
		perror("Start");
		exit(1);
	}

	degree = atof(argv[3]);

	degreemax = atoi(argv[4]);

	if (degreemax > 360.0)
		degreemax = 360.0;

	degreestep = atof(argv[5]);

	observer_height = atof(argv[2]);

	coordinate_file_name = argv[6];

	fp = fopen (coordinate_file_name, "r");

	if (fp == NULL)
	{
		perror("Coordinates");
		exit(1);
	}

	coordinates *coordinates1 = malloc(sizeof(coordinates) * COORDINATES_MAX);

	int k;

	for (k=0; k<COORDINATES_MAX; k++)
	{
		coordinates1[k].x = -3;
		coordinates1[k].y = -3;
	}

	number_of_coordinates = 0;

	//read the coordinates that will be checked from the text file
	while(fscanf(fp, "%f %f", &coordinates1[number_of_coordinates].x, &coordinates1[number_of_coordinates].y) == 2)
    	number_of_coordinates++;

	for (k=0; k<number_of_coordinates; k++)
	{
		if (coordinates1[k].x == -1 && coordinates1[k].y == -1)
		{
			all_pixels = TRUE;
			x0 = 0;
			y0 = 0;
			printf("-- ALL pixels of the image will be viewshed/n/n");

			break;
		}
		else
		{
			x0 = coordinates1[0].x;
			y0 = coordinates1[0].y;
		}
	}

	fclose(fp);

	int value = 255/number_of_coordinates;// how much to change the rgb values

	int offset = y0 * IMAGE_MAX_X + x0;

	if (offset < 0 || offset > (sbuf.st_size - 1)/nrofbytes_per_stored_element)
	{
		fprintf(stderr, "mman :offset must be in range 0-%ld\n", (sbuf.st_size-1)/nrofbytes_per_stored_element);
		exit(1);
	}

	data = mmap ((caddr_t)0 ,sbuf.st_size ,PROT_READ ,MAP_SHARED ,fd ,0);

	if (data == (float*)(caddr_t)(-1))
	{
		perror("map");
		exit(1);
	}

	char filenamesave[120];
	bool rc;
	bmpfile_t *bmp;
	int depthbits = 8;

	//blue green red ,not rgb
	rgb_pixel_t pixelred = {0 ,0 ,225 ,0};
	rgb_pixel_t pixel = {0 ,0 ,0 ,0};

	struct timespec start ,end;

	char myprefix[] = "PostProcessBMP";

	sprintf(filenamesave ,"%dx%d@%4.2f_%s.bmp", y0 ,x0 ,observer_height ,myprefix);

	//create a bmp image
	if ((bmp = bmp_create(IMAGE_MAX_X, IMAGE_MAX_Y ,depthbits)) == NULL)
	{
		printf("FAIL in creating BMP\n");
		return 1;
	}

	int j,i;

	//firs the image is black
	for (j = 0; j < IMAGE_MAX_X; j++)
		for (i = 0; i < IMAGE_MAX_Y; i++)
			bmp_set_pixel(bmp ,j ,i ,pixel);
	//start time
	clock_gettime(CLOCK_MONOTONIC ,&start);

	printf("Observer height set to [%f]\n",observer_height);

	float idegree;

	if (all_pixels)
		number_of_coordinates = IMAGE_MAX_X * IMAGE_MAX_Y;

	for (i = 0; i < number_of_coordinates; i++)
	{
		if (!all_pixels)
		{
			x0 = coordinates1[i].x;
			y0 = coordinates1[i].y;
		}

		//in this case all pixels will be checked
		else
		{
			if (x0 < IMAGE_MAX_X)
				x0++;

			else
			{
				x0 = 0;
				y0++;
			}
		}

		printf ("Processing y0, x0 => (%d,%d)\n",y0 ,x0);


		//rotation
		for (idegree = degree; idegree < degreemax; idegree=idegree+degreestep){
			compute_height_degree(data ,idegree ,x0 ,y0 ,observer_height ,&xmax0 ,&ymax0 ,&pointheight0, value, bmp);

		}
	}

	clock_gettime(CLOCK_MONOTONIC ,&end);

    	const int DAS_NANO_SECONDS_IN_SEC = 1000000000;
    	long timeElapsed_s = end.tv_sec - start.tv_sec;
    	long timeElapsed_n = end.tv_nsec - start.tv_nsec;

    	//If we have a negative number in timeElapsed_n , borrow a carry from seconds
    	if ( timeElapsed_n < 0 )
	{
        	timeElapsed_n = DAS_NANO_SECONDS_IN_SEC + timeElapsed_n;
        	timeElapsed_s--;
    	}

	printf("Compute Time: %ld.%09ld secs \n",timeElapsed_s,timeElapsed_n);

	for (i = 0; i < number_of_coordinates; i++){
                if (!all_pixels){

                        x0 = coordinates1[i].x;
                        y0 = coordinates1[i].y;

			//draw a red cross in order to represent observer
        		for (j = -30; j < 30; j++){

                		bmp_set_pixel(bmp ,x0 + j ,y0 ,pixelred);
                		bmp_set_pixel(bmp ,x0 ,y0 + j ,pixelred);
        		}
		}
	}

	printf("Done... Saving to file %s: \n",filenamesave);
	rc = bmp_save(bmp, filenamesave);

	if (rc == FALSE)
		printf("FAIL\n");

	else
		printf("SUCCESS\n");

	bmp_destroy(bmp);

	return 0;
}
