#include <stdint.h>
#include <vector>
#include <fstream>
#include <iostream>

const int WIDTH = 2;	// 1280
const int HEIGHT = 2;	// 720

void writePPM(char image[WIDTH * HEIGHT * 3], std::string name)
{
	std::ofstream ppm;
	ppm.open(name, std::ios::binary);
	ppm << "P6" << " " << WIDTH << " " << HEIGHT << " " << "255" << '\n';
	size_t size = WIDTH * HEIGHT * 3;
	ppm.write(image, size);
	ppm.close();
}

void readPPM(char image[WIDTH * HEIGHT * 3], std::string name)
{
	std::ifstream ppm;
	ppm.open(name, std::ios::binary);
	while (ppm.get() != '\n');
	size_t size = WIDTH * HEIGHT * 3;
	ppm.read(image, size);
	ppm.close();
}

void get_image(volatile uint32_t* hdmi, char image[WIDTH * HEIGHT * 3])
{
	// Assuming accessing the HDMI input automatically scrolls it to the next pixel
	// Assuming there is no header/footer in each frame
	// Assuming scanning horizontally from top left to bottom right
	for (unsigned int i = 0; i < WIDTH; ++i)
	{
		for (unsigned int j = 0; j < HEIGHT; ++j)
		{
			uint32_t pixel = *hdmi;
			// Assuming little endian (0GBR)
			image[j * WIDTH + i * HEIGHT + 0] = pixel & 0x000000FF;
			image[j * WIDTH + i * HEIGHT + 1] = (pixel & 0x0000FF00) >> 8;
			image[j * WIDTH + i * HEIGHT + 2] = (pixel & 0x00FF0000) >> 16;
		}
	}
}

void detect_fingers(volatile uint32_t* hdmi, std::vector<std::vector<int>*> graph, int& exception)
{

}

int main()
{
	char* image = new char[WIDTH * HEIGHT * 3];
	readPPM(image, "image.ppm");
	writePPM(image, "copy.ppm");

	return 0;
}