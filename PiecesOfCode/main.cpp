#include <stdint.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>

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

struct coord
{
	int x;
	int y;
	coord(int x_, int y_) : x(x_), y(y_) {};
};

void getFingertip(char image[IMAGE_SIZE], int& xpos, int& ypos, int borderSize)
{
	// Hands are white
	int x = 0, y = 0;
	xpos = ypos = 0;

	for (unsigned int i = borderSize; i < WIDTH - borderSize; ++i)
	{
		if (image[(borderSize * WIDTH + i) * 3])
		{
			x = i;
			y = borderSize;
			goto gotTip;
		}
		if (image[((WIDTH - borderSize) * HEIGHT - i - 1) * 3])
		{
			x = i;
			y = HEIGHT - borderSize;
			goto gotTip;
		}
		//std::cout << image[3 * i];
	}

	for (unsigned int i = borderSize; i < HEIGHT - borderSize; ++i)
	{
		if (image[(i * WIDTH + borderSize) * 3])
		{
			x = borderSize;
			y = i;
			goto gotTip;
		}
		if (image[((i + 1) * WIDTH - 1 - borderSize) * 3])
		{
			x = WIDTH - borderSize;
			y = i;
			goto gotTip;
		}
	}

	xpos = ypos = -1;
	return;

gotTip:

	// Flood hand
	image[(y * WIDTH + x) * 3] = 0;
	std::queue<coord> points;
	points.push(coord(x, y));
	int s = 0;

	while (!points.empty())
	{
		coord point = points.front();
		points.pop();
		x = point.x, y = point.y;
		int index = (y * WIDTH + x) * 3;

		// Right
		if (x < WIDTH - 1)
		{
			if (image[index + 3])
			{
				image[index + 3] = 0;
				points.push(coord(x + 1, y));
			}

			// Right-Down
			if (y < HEIGHT - 1)
			{
				if (image[index + 3 * (WIDTH + 1)])
				{
					image[index + 3 * (WIDTH + 1)] = 0;
					points.push(coord(x + 1, y + 1));
				}
			}

			// Right-Up
			if (y > 0)
			{
				if (image[index - 3 * (WIDTH - 1)])
				{
					image[index - 3 * (WIDTH - 1)] = 0;
					points.push(coord(x + 1, y - 1));
				}
			}
		}

		// Left
		if (x > 0)
		{
			if (image[index - 3])
			{
				image[index - 3] = 0;
				points.push(coord(x - 1, y));
			}

			// Left-Down
			if (y < HEIGHT - 1)
			{
				if (image[index + 3 * (WIDTH - 1)])
				{
					image[index + 3 * (WIDTH - 1)] = 0;
					points.push(coord(x - 1, y + 1));
				}
			}

			// Left-Up
			if (y > 0)
			{
				if (image[index - 3 * (WIDTH + 1)])
				{
					image[index - 3 * (WIDTH + 1)] = 0;
					points.push(coord(x - 1, y - 1));
				}
			}
		}

		// Down
		if (y < HEIGHT - 1)
		{
			if (image[index + 3 * WIDTH])
			{
				image[index + 3 * WIDTH] = 0;
				points.push(coord(x, y + 1));
			}
		}

		// Up
		if (y > 0)
		{
			if (image[index - 3 * WIDTH])
			{
				image[index - 3 * WIDTH] = 0;
				points.push(coord(x, y - 1));
			}
		}

		image[index + 1] = 0;
		++s;
		/*
		if (s % 100 == 0)
		{
			std::stringstream ss;
			ss << "test_hands/f" << s / 100 << ".ppm";
			write_PPM(image, ss.str());
			std::cout << "Wrote " << s / 100 << std::endl;
			std::cin.ignore();
		}
		*/
	}

	std::cout << s << std::endl;
	xpos = x, ypos = y;
	image[(y * WIDTH + x) * 3] = image[(y * WIDTH + x) * 3 + 1] = image[(y * WIDTH + x) * 3 + 2] = 255;
}

void detect_fingers(volatile uint32_t* hdmi, std::vector<std::vector<int>*>& graph, int& exception)
{
	// Read current image
	char* image = new char[IMAGE_SIZE];
	//get_image(hdmi, image);
	read_PPM(image, "test_hands/map_covered_hands.ppm");

	// Complex subtraction
	char* sub = new char[IMAGE_SIZE];
	int index;
	int w = 2, h = 2;
	int minr, ming, minb, min;
	for (unsigned int i = w; i < HEIGHT - w; ++i)
	{
		for (unsigned int j = h; j < WIDTH - h; ++j)
		{
			index = (i * WIDTH + j) * 3;
			min = 765;	// 255 * 3
			minr = ming = minb = 255;
			for (int k = -w; k <= w; ++k)
			{
				for (int l = -h; l <= h; ++l)
				{
					int temp = index + (k + l * WIDTH) * 3;
					int r = image[index] - map[temp];
					int g = image[index + 1] - map[temp + 1];
					int b = image[index + 2] - map[temp + 2];
					r = r > 0 ? r : -r;
					g = g > 0 ? g : -g;
					b = b > 0 ? b : -b;
					
					if (min > r + g + b)
					{
						minr = r;
						ming = g;
						minb = b;
						min = r + g + b;
					}
				}
			}

			sub[index] = minr;
			sub[index + 1] = ming;
			sub[index + 2] = minb;
		}
	}

	for (unsigned int i = 0; i < WIDTH; ++i)
	{
		for (unsigned int j = 0; j < h; ++j)
		{
			sub[3 * (i + j * WIDTH)] = sub[3 * (i + j * WIDTH) + 1] = sub[3 * (i + j * WIDTH) + 2] = 0;
			sub[(WIDTH * (HEIGHT - j) - i - 1) * 3] = sub[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 1] = sub[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 2] = 0;
		}
	}

	for (unsigned int i = 0; i < HEIGHT; ++i)
	{
		for (unsigned int j = 0; j < w; ++j)
		{
			sub[(i * WIDTH + j) * 3] = sub[(i * WIDTH + j) * 2 + 1] = sub[(i * WIDTH + j) * 2 + 2] = 0;
			sub[((i + 1) * WIDTH - j - 1 ) * 3] = sub[((i + 1) * WIDTH - j - 1) * 3 + 1] = sub[((i + 1) * WIDTH - j - 1) * 3 + 2] = 0;
		}
	}

	delete image;
	image = sub;

	write_PPM(image, "test_hands/subtract.ppm");

	// De-noise
	int threshold = 15;
	for (unsigned int i = 0; i < IMAGE_SIZE; i += 3)
	{
		if (image[i] + image[i + 1] + image[i + 2] > 3 * threshold)
			setPixel(image, i, 255);
		else
			setPixel(image, i, 0);
	}

	for (unsigned int i = 1; i < HEIGHT - 1; ++i)
	{
		for (unsigned int j = 1; j < WIDTH - 1; ++j)
		{
			index = (i * WIDTH + j) * 3;
			
			if (image[index] + image[index + 1] + image[index + 2] == 0)
			{
				setPixel(image, index - WIDTH * 3 - 3, 0);
				setPixel(image, index - WIDTH * 3, 0);
				setPixel(image, index - WIDTH * 3 + 3, 0);
				setPixel(image, index - 3, 0);
				setPixel(image, index, 0);
			}
		}
	}

	write_PPM(image, "test_hands/denoise.ppm");

	// Fill
	for (int p = 0; p < 3; ++p)
	{
		char* fill = new char[IMAGE_SIZE];
		for (unsigned int i = w; i < HEIGHT - w; ++i)
		{
			for (unsigned int j = h; j < WIDTH - h; ++j)
			{
				index = (i * WIDTH + j) * 3;
				setPixel(fill, index, 0);

				for (int k = -w; k <= w; ++k)
				{
					for (int l = -h; l <= h; ++l)
					{
						int temp = index + (k + l * WIDTH) * 3;
						if (image[temp])
						{
							setPixel(fill, index, 255);
							goto go;
						}
					}
				}

			go:
				continue;
			}
		}

		delete image;
		image = fill;

		for (unsigned int i = 0; i < WIDTH; ++i)
		{
			for (unsigned int j = 0; j < h; ++j)
			{
				image[3 * (i + j * WIDTH)] = image[3 * (i + j * WIDTH) + 1] = image[3 * (i + j * WIDTH) + 2] = 0;
				image[(WIDTH * (HEIGHT - j) - i - 1) * 3] = image[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 1] = image[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 2] = 0;
			}
		}

		for (unsigned int i = 0; i < HEIGHT; ++i)
		{
			for (unsigned int j = 0; j < w; ++j)
			{
				image[(i * WIDTH + j) * 3] = image[(i * WIDTH + j) * 3 + 1] = image[(i * WIDTH + j) * 3 + 2] = 0;
				image[((i + 1) * WIDTH - j - 1) * 3] = image[((i + 1) * WIDTH - j - 1) * 3 + 1] = image[((i + 1) * WIDTH - j - 1) * 3 + 2] = 0;
			}
		}
	}

	write_PPM(image, "test_hands/fill.ppm");

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

	delete hit;
	write_PPM(image, "test_hands/flood.ppm");

	// Get fingertips
	int x, y;
	getFingertip(image, x, y, w);
	std::cout << x << ", " << y << std::endl;
	write_PPM(image, "test_hands/f1.ppm");
	getFingertip(image, x, y, w);
	std::cout << x << ", " << y << std::endl;
	write_PPM(image, "test_hands/f2.ppm");
	std::cin.ignore();

	/*
	// Create vector field
	int* vec = new int[WIDTH * HEIGHT * 2];
	int vecIndex, dist, diagDist;
	float diag;
	for (int i = 0; i < HEIGHT; ++i)
	{
		for (int j = 0; j < WIDTH; ++j)
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
			while (!flag && dist < WIDTH + HEIGHT)
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

			//float mag = 1.0 / std::sqrt(vec[index] * vec[index] + vec[index + 1] * vec[index + 1]);
			float mag = 3;
			grad[index] = dx * mag;
			grad[index + 1] = dy * mag;
		}
	}

	delete vec;

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
		int div = grad[i * 2] + grad[i * 2 + 1];
		image[i * 3] = div > 0 ? div : -div;
		image[i * 3 + 1] = image[i * 3 + 2] = 0;
	}

	write_PPM(image, "test_hands/divergence.ppm");
	
	delete grad;
	*/

	delete image;
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

	delete map;

	return 0;
}