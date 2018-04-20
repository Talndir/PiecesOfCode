#include <stdint.h>
#include <vector>
#include <fstream>

const int WIDTH = 2;
const int HEIGHT = 2;

void imageToPPM(char image[WIDTH * HEIGHT * 3])
{
	std::ofstream ppm;
	ppm.open("image.ppm", std::ios::binary);
	ppm << "P6" << " " << WIDTH << " " << HEIGHT << " " << "255" << "\n";
	size_t size = WIDTH * HEIGHT * 3;
	ppm.write(image, size);
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
	image[0] = 100;
	image[1] = 200;
	image[2] = 255;
	image[3] = 100;
	image[4] = 200;
	image[5] = 255;
	image[6] = 100;
	image[7] = 200;
	image[8] = 255;
	image[9] = 10;
	image[10] = 20;
	image[11] = 255;
	imageToPPM(image);

	return 0;
}