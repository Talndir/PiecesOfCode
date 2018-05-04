//#define DEBUG

#include <stdint.h>
#include <queue>
#include <map>
#include <algorithm>

#ifdef DEBUG
#include <fstream>
#include <iostream>
#include <sstream>
#endif

#define SQRT2_RECIP 0.70710678118

const unsigned int WIDTH = 1024;	// 1280
const unsigned int HEIGHT = 682;	// 720
const unsigned int IMAGE_SIZE = WIDTH * HEIGHT * 3;
char* map;

#ifdef DEBUG
// Write PPM image
void write_PPM(char image[IMAGE_SIZE], std::string name)
{
	std::ofstream ppm;
	ppm.open(name, std::ios::binary);
	ppm << "P6" << " " << WIDTH << " " << HEIGHT << " " << "255" << '\n';
	ppm.write(image, IMAGE_SIZE);
	ppm.close();
}

// Read PPM image
void read_PPM(char image[IMAGE_SIZE], std::string name)
{
	std::ifstream ppm;
	ppm.open(name, std::ios::binary);
	while (ppm.get() != '\n');
	ppm.read(image, IMAGE_SIZE);
	ppm.close();
}
#endif

// Get image frm HDMI
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

// Set all three channels to specified value
void setPixel(char image[], int index, char value)
{
	image[index] = image[index + 1] = image[index + 2] = value;
}

// 2D integer coordinate
struct coord
{
	int x;
	int y;
	coord() : x(0), y(0) {};
	coord(int x_, int y_) : x(x_), y(y_) {};
};

// Edge between two vertices in graph
struct edge
{
	coord a, b;			// x-y coord
	int p1, p2;			// Pixel index
	int v1, v2;	// Vertex index in graph
	int index;			// Index of b in a's adjacency list (+1 *3 to get actual index)
	float nx, ny;		// Normal vector (normalised)

	edge()
	{
		edge(0, 0, 0);
	};

	edge(int p1_, int p2_, int index_) : p1(p1_), p2(p2_), index(index_)
	{
		a = coord(p1 % WIDTH, p1 / WIDTH);
		b = coord(p2 % WIDTH, p2 / WIDTH);
		v1 = v2 = 0;

		calcNorm();
	}

	void calcNorm()
	{
		float ux = b.x - a.x;
		float uy = b.y - a.y;
		float mag = std::sqrt(ux * ux + uy * uy);
		if (mag != 0)
		{
			nx = uy / mag;
			ny = ux / mag;
		}
		else
		{
			nx = ny = 0;
		}
	}

	void setIndices(int v1_, int v2_)
	{
		v1 = v1_, v2 = v2_;
	}
};

// Subtract blank map with given image, match in an area
void subtract(char image[IMAGE_SIZE], int w, int h)
{
	char* sub = new char[IMAGE_SIZE];
	for (unsigned int i = w; i < HEIGHT - w; ++i)
	{
		for (unsigned int j = h; j < WIDTH - h; ++j)
		{
			int index = (i * WIDTH + j) * 3;
			int min = 765;	// 255 * 3
			int minr = 255, ming = 255, minb = 255;
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
		for (int j = 0; j < h; ++j)
		{
			sub[3 * (i + j * WIDTH)] = sub[3 * (i + j * WIDTH) + 1] = sub[3 * (i + j * WIDTH) + 2] = 0;
			sub[(WIDTH * (HEIGHT - j) - i - 1) * 3] = sub[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 1] = sub[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 2] = 0;
		}
	}

	for (unsigned int i = 0; i < HEIGHT; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			sub[(i * WIDTH + j) * 3] = sub[(i * WIDTH + j) * 2 + 1] = sub[(i * WIDTH + j) * 2 + 2] = 0;
			sub[((i + 1) * WIDTH - j - 1) * 3] = sub[((i + 1) * WIDTH - j - 1) * 3 + 1] = sub[((i + 1) * WIDTH - j - 1) * 3 + 2] = 0;
		}
	}

	delete image;
	image = sub;
}

// Remove noise
void denoise(char image[IMAGE_SIZE], int threshold)
{
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
			int index = (i * WIDTH + j) * 3;

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
}

// Slightly enlarge and fill some gaps in white regions
void fill(char image[IMAGE_SIZE], int w, int h)
{
	char* fill = new char[IMAGE_SIZE];
	for (unsigned int i = w; i < HEIGHT - w; ++i)
	{
		for (unsigned int j = h; j < WIDTH - h; ++j)
		{
			int index = (i * WIDTH + j) * 3;
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
		for (int j = 0; j < h; ++j)
		{
			image[3 * (i + j * WIDTH)] = image[3 * (i + j * WIDTH) + 1] = image[3 * (i + j * WIDTH) + 2] = 0;
			image[(WIDTH * (HEIGHT - j) - i - 1) * 3] = image[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 1] = image[(WIDTH * (HEIGHT - j) - i - 1) * 3 + 2] = 0;
		}
	}

	for (unsigned int i = 0; i < HEIGHT; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			image[(i * WIDTH + j) * 3] = image[(i * WIDTH + j) * 3 + 1] = image[(i * WIDTH + j) * 3 + 2] = 0;
			image[((i + 1) * WIDTH - j - 1) * 3] = image[((i + 1) * WIDTH - j - 1) * 3 + 1] = image[((i + 1) * WIDTH - j - 1) * 3 + 2] = 0;
		}
	}
}

// Fill unreachable areas such as holes
void flood(char image[IMAGE_SIZE])
{
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
}

// Return first found fingertip, each run produces new fingertip
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
	}

	xpos = x, ypos = y;
	image[(y * WIDTH + x) * 3] = image[(y * WIDTH + x) * 3 + 1] = image[(y * WIDTH + x) * 3 + 2] = 255;
}

// Create list of all edges in graph
void createEdges(int* graph[], edge edges[])
{
	int* vertex;
	int size = sizeof(graph) / sizeof(int);
	edge* edges2 = new edge[size * size];
	std::map<int, int> indices;
	int index = 0;

	for (unsigned int i = 0; i < size; ++i)
	{
		vertex = graph[i];
		int pix = vertex[0];
		indices.insert(std::make_pair(pix, (int)i));

		for (unsigned int j = 3; j < sizeof(vertex) / sizeof(int); j += 3)
			edges2[index++] = edge(pix, vertex[j], j / 3 - 1);
	}

	delete edges;
	edges = new edge[index];
	edges = edges2;
	for (unsigned int i = 0; i < index; ++i)
		edges[i].setIndices(indices.at(edges[i].p1), indices.at(edges[i].p2));
}

// Bisect nearest edge to point specified
void addVertex(int x, int y, int* graph[], edge edges[])
{
	float dist = WIDTH + HEIGHT;
	int index = -1;
	edge e;
	int size = sizeof(edges) / sizeof(edge);

	for (unsigned int i = 0; i < size; ++i)
	{
		e = edges[i];
		float d = (e.a.x - x) * e.nx + (e.a.y - y) * e.ny;

		if (d < dist)
		{
			dist = d;
			index = i;
		}
	}

	// Add new vertex and edges
	int gsize = sizeof(graph) / sizeof(int*);
	int *(*g2) = new int*[gsize + 1];
	std::copy(graph, graph + gsize, g2);
	delete graph;
	graph = g2;
	int* vertex = new int[9];
	graph[gsize] = vertex;
	e = edges[index];
	x += dist * e.nx;
	y += dist * e.ny;
	int pix = y * WIDTH + x;
	int weight = (graph[e.v1][1] + graph[e.v2][1]) / 2;
	vertex[0] = pix;
	vertex[1] = weight;	// Weight is average of connecting vertices
	vertex[2] = 1;		// Is starting point
	vertex[3] = e.p1;
	vertex[4] = graph[e.v1][1];
	vertex[5] = 0;
	vertex[6] = e.p2;
	vertex[7] = graph[e.v2][1];
	vertex[8] = 0;

	edge e2;
	int i = 0;
	for (; i < size; ++i)
	{
		if (edges[i].p2 == e.p1 && edges[i].p1 == e.p2)
		{
			e2 = edges[i];
			break;
		}
	}

	// Overwrite previous edges
	edge* edges2 = new edge[size + 2];
	std::copy(edges, edges + size, edges2);
	delete edges;
	edges = edges2;

	edges[index] = edge(pix, graph[e.v1][1], 0);
	edges[i] = edge(pix, graph[e.v2][1], 1);
	edges[size] = edge(e.p1, pix, e.index);
	edges[size + 1] = edge(e2.p1, pix, e2.index);

	int ind = (e.index + 1) * 3;
	vertex = graph[e.v1];
	vertex[ind] = pix;
	vertex[ind + 1] = weight;
	vertex[ind + 2] = 0;
	ind = (e2.index + 1) * 3;
	vertex[ind] = pix;
	vertex[ind + 1] = weight;
	vertex[ind + 2] = 0;
}

// Take image and add start/end vertices to graph where fingertips are detected
void detect_fingers(volatile uint32_t* hdmi, int* graph[], int& exception)
{
	// Read current image
	char* image = new char[IMAGE_SIZE];
	get_image(hdmi, image);
	
#ifdef DEBUG
	delete image;
	image = new char[IMAGE_SIZE];
	read_PPM(image, "test_hands/map_covered_hands.ppm");
#endif

	// Complex subtraction
	int w = 2;
	int h = w;
	subtract(image, w, h);
	//write_PPM(image, "test_hands/subtract.ppm");

	// De-noise
	int threshold = 15;
	denoise(image, threshold);

#ifdef DEBUG
	write_PPM(image, "test_hands/denoise.ppm");
#endif

	// Fill
	for (int p = 0; p < 3; ++p)
		fill(image, w, h);

#ifdef DEBUG
	write_PPM(image, "test_hands/fill.ppm");
#endif

	// Flood
	flood(image);
#ifdef DEBUG
	write_PPM(image, "test_hands/flood.ppm");
#endif

	// Get fingertips
	// Assuming immediate return if not exactly two fingers
	int x, y;
	getFingertip(image, x, y, w);

#ifdef DEBUG
	std::cout << x << ", " << y << std::endl;
	write_PPM(image, "test_hands/f1.ppm");
#endif

	int a, b;
	getFingertip(image, a, b, w);
	if (a == -1)
	{
		exception = -1;
		return;
	}

#ifdef DEBUG
	std::cout << a << ", " << b << std::endl;
	write_PPM(image, "test_hands/f2.ppm");
#endif

	int c, d;
	getFingertip(image, c, d, w);
	if (c != -1)
	{
		exception = 1;
		return;
	}

	// Add vertices to map
	edge* edges = nullptr;
	createEdges(graph, edges);
	addVertex(x, y, graph, edges);
	addVertex(a, b, graph, edges);

	delete image;
	exception = 0;
}

// Gets map from HDMI
void get_map(volatile uint32_t* hdmi)
{
	// Read in map
	map = new char[IMAGE_SIZE];
	get_image(hdmi, map);
}

int main()
{
	volatile uint32_t* hdmi = new volatile uint32_t;	// HDMI in

	// Read in original map
	get_map(hdmi);
#ifdef DEBUG
	delete map;
	map = new char[IMAGE_SIZE];
	read_PPM(map, "map.ppm");
#endif

	// Operate
	int*(*graph) = nullptr;						// Graph
	int exception;								// Exception
	detect_fingers(hdmi, graph, exception);

	delete map;

	return 0;
}