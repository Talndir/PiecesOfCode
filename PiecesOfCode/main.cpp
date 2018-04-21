#include <stdint.h>
#include <vector>
#include <fstream>
#include <iostream>

#define SQRT2_RECIP 0.70710678118

const int WIDTH = 1024;	// 1280
const int HEIGHT = 682;	// 720
const int IMAGE_SIZE = WIDTH * HEIGHT * 3;
char* map;

void write_PPM(char image[IMAGE_SIZE], std::string name)
{
	std::ofstream ppm;
	ppm.open(name, std::ios::binary);
	ppm << "P6" << " " << WIDTH << " " << HEIGHT << " " << "255" << '\n';
	ppm.write(image, IMAGE_SIZE);
	ppm.close();
}

void read_PPM(char image[IMAGE_SIZE], std::string name)
{
	std::ifstream ppm;
	ppm.open(name, std::ios::binary);
	while (ppm.get() != '\n');
	ppm.read(image, IMAGE_SIZE);
	ppm.close();
}

void get_image(volatile uint32_t* hdmi, char image[IMAGE_SIZE])
{
	// Assuming accessing the HDMI input automatically scrolls it to the next pixel
	// Assuming there is no header/footer in each frame
	// Assuming scanning horizontally from top left to bottom right
	for (unsigned int i = 0; i < WIDTH * HEIGHT; ++i)
	{
		uint32_t pixel = *hdmi;
		// Assuming little endian (0BGR)
		image[i * 3 + 0] = pixel & 0x000000FF;
		image[i * 3 + 1] = (pixel & 0x0000FF00) >> 8;
		image[i * 3 + 2] = (pixel & 0x00FF0000) >> 16;
	}
}

void setPixel(char image[], int index, char value)
{
	image[index] = image[index + 1] = image[index + 2] = value;
}

void detect_fingers(volatile uint32_t* hdmi, std::vector<std::vector<int>*>& graph, int& exception)
{
	// Read current image
	char* image = new char[IMAGE_SIZE];
	//get_image(hdmi, image);
	read_PPM(image, "map_covered.ppm");

	// Subtract original image
	int sub = 0;
	int threshold = 8;
	for (unsigned int i = 0; i < IMAGE_SIZE; ++i)
	{
		sub = image[i] - map[i];
		sub < 0 ? sub = -sub : NULL;
		sub < threshold ? sub = 0 : NULL;
		image[i] = sub;
	}

	// De-noise
	int index;
	for (unsigned int i = 1; i < HEIGHT - 1; ++i)
	{
		for (unsigned int j = 1; j < WIDTH - 1; ++j)
		{
			index = (i * WIDTH + j) * 3;
			
			if (image[index] + image[index + 1] + image[index + 2] < 3 * threshold)
			{
				setPixel(image, index - WIDTH * 3 - 3, 0);
				setPixel(image, index - WIDTH * 3, 0);
				setPixel(image, index - WIDTH * 3 + 3, 0);
				setPixel(image, index - 3, 0);
				setPixel(image, index, 0);
			}
		}
	}

	for (unsigned int i = 0; i < IMAGE_SIZE; ++i)
	{
		if (image[i] + image[i + 1] + image[i + 2] > 0)
			setPixel(image, i, 255);
	}

	// Flood
	bool* hit = new bool[WIDTH * HEIGHT];
	for (unsigned int i = 0; i < WIDTH * HEIGHT; ++i)
		hit[i] = false;

	for (unsigned int i = 0; i < HEIGHT; ++i)
	{
		int hitIndex = i * WIDTH;
		int index = hitIndex * 3;
		int x = hitIndex + WIDTH;
		while (image[index] + image[index + 1] + image[index + 2] == 0 && hitIndex < x - 1)
		{
			hit[hitIndex] = true;
			++hitIndex;
			index += 3;
		}
		
		hitIndex = (i + 1) * WIDTH - 1;
		index = hitIndex * 3;
		x = hitIndex - WIDTH;
		while (image[index] + image[index + 1] + image[index + 2] == 0 && hitIndex > x + 1)
		{
			hit[hitIndex] = true;
			--hitIndex;
			index -= 3;
		}
	}
	for (unsigned int i = 0; i < WIDTH; ++i)
	{
		int hitIndex = i;
		int index = hitIndex * 3;
		int y = WIDTH * (HEIGHT - 1);
		while (image[index] + image[index + 1] + image[index + 2] == 0 && hitIndex < y)
		{
			hit[hitIndex] = true;
			hitIndex += WIDTH;
			index += WIDTH * 3;
		}

		index = IMAGE_SIZE - (3 * (WIDTH - i));
		hitIndex = index / 3;
		y = WIDTH;
		while (image[index] + image[index + 1] + image[index + 2] == 0 && hitIndex > y)
		{
			hit[hitIndex] = true;
			hitIndex -= WIDTH;
			index -= WIDTH * 3;
		}
	}
	
	for (unsigned int i = 0; i < WIDTH * HEIGHT; ++i)
	{
		if (!hit[i])
			setPixel(image, i * 3, 255);
	}

	// Create vector field
	int* vec = new int[WIDTH * HEIGHT * 2];
	int vecIndex, dist, diagDist;
	float diag;
	for (unsigned int i = 0; i < HEIGHT; ++i)
	{
		for (unsigned int j = 0; j < WIDTH; ++j)
		{
			index = (i * WIDTH + j) * 3;
			vecIndex = (i * WIDTH + j) * 2;

			if (!image[index])
			{
				vec[vecIndex] = 0;
				vec[vecIndex + 1] = 0;
				continue;
			}

			diag = SQRT2_RECIP;
			dist = diagDist = 1;
			bool flag = false;
			while (!flag)
			{
				// Right
				if (dist + j < WIDTH)
				{
					if (!image[index + 3 * dist])
					{
						vec[vecIndex] = dist;
						vec[vecIndex + 1] = 0;
						flag = true;
					}
				}

				if (diagDist + j < WIDTH)
				{
					// Right-Down
					if (diagDist + i < HEIGHT)
					{
						if (!image[index + 3 * diagDist * (WIDTH + 1)])
						{
							vec[vecIndex] = diagDist;
							vec[vecIndex + 1] = diagDist;
							flag = true;
						}
					}

					// Right-Up
					if (i - diagDist > 0)
					{
						if (!image[index - 3 * diagDist * (WIDTH - 1)])
						{
							vec[vecIndex] = diagDist;
							vec[vecIndex + 1] = -diagDist;
							flag = true;
						}
					}
				}

				// Left
				if (j - dist > 0)
				{
					if (!image[index - 3 * dist])
					{
						vec[vecIndex] = -dist;
						vec[vecIndex + 1] = 0;
						flag = true;
					}
				}

				if (j - diagDist > 0)
				{
					// Left-Down
					if (diagDist + i < HEIGHT)
					{
						if (!image[index + 3 * diagDist * (WIDTH - 1)])
						{
							vec[vecIndex] = -diagDist;
							vec[vecIndex + 1] = diagDist;
							flag = true;
						}
					}

					// Left-Up
					if (i - diagDist > 0)
					{
						if (!image[index - 3 * diagDist * (WIDTH + 1)])
						{
							vec[vecIndex] = -diagDist;
							vec[vecIndex + 1] = -diagDist;
							flag = true;
						}
					}
				}
				
				// Down
				if (dist + i < HEIGHT)
				{
					if (!image[index + 3 * dist * WIDTH])
					{
						vec[vecIndex] = 0;
						vec[vecIndex + 1] = dist;
						flag = true;
					}
				}

				// Up
				if (i - dist > 0)
				{
					if (!image[index - 3 * dist * WIDTH])
					{
						vec[vecIndex] = 0;
						vec[vecIndex + 1] = -dist;
						flag = true;
					}
				}
				
				++dist;
				diag += SQRT2_RECIP;
				diagDist = (int)diag;
			}
		}
	}

	// Calculate derivative of vector field
	int* grad = new int[WIDTH * HEIGHT * 2];
	for (unsigned int i = 1; i < HEIGHT - 1; ++i)
	{
		for (unsigned int j = 1; j < WIDTH - 1; ++j)
		{
			int dx = 0, dy = 0;
			int index = (i * WIDTH + j) * 2;
			
			dx += vec[index + 2] + SQRT2_RECIP * (vec[index - WIDTH * 2 + 2] + vec[index + WIDTH * 2 + 2]);
			dx -= vec[index - 2] + SQRT2_RECIP * (vec[index - WIDTH * 2 - 2] + vec[index + WIDTH * 2 - 2]);
			dy += vec[index + WIDTH * 2 + 1] + SQRT2_RECIP * (vec[index + WIDTH * 2 - 1] + vec[index + WIDTH * 2 + 3]);
			dy -= vec[index - WIDTH * 2 + 1] + SQRT2_RECIP * (vec[index - WIDTH * 2 - 1] + vec[index - WIDTH * 2 + 3]);

			grad[index] = dx / 2;
			grad[index + 1] = dy / 2;
		}
	}

	for (unsigned int i = 0; i < WIDTH; ++i)
	{
		grad[2 * i] = grad[2 * i + 1] = 0;
		grad[(WIDTH * HEIGHT - i - 1) * 2] = grad[(WIDTH * HEIGHT - i - 1) * 2 + 1] = 0;
	}

	for (unsigned int i = 0; i < HEIGHT; ++i)
	{
		grad[i * WIDTH * 2] = grad[i * WIDTH * 2 + 1] = 0;
		grad[(i + 1) * WIDTH * 2 - 2] = grad[(i + 1) * WIDTH * 2 - 1] = 0;
	}

	// Write divergence scalar field image
	for (unsigned int i = 0; i < WIDTH * HEIGHT; ++i)
	{
		
		int vx = vec[i * 2];
		int vy = vec[i * 2 + 1];
		image[i * 3] = (int)std::sqrt(vx * vx + vy * vy);
		
		int div = grad[i * 2] + grad[i * 2 + 1];
		image[i * 3] = div > 0 ? div : -div;
		image[i * 3 + 1] = image[i * 3 + 2] = 0;
	}

	write_PPM(image, "processed.ppm");
}

int main()
{
	// Read in original map
	map = new char[IMAGE_SIZE];
	read_PPM(map, "map.ppm");
	
	// Operate
	volatile uint32_t* hdmi = new volatile uint32_t;
	std::vector<std::vector<int>*> graph;
	int exception;
	detect_fingers(hdmi, graph, exception);

	return 0;
}