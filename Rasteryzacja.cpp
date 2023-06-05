#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include <chrono>
#include <thread>
#include <iomanip>
#include <sstream>
#include <mutex>
#include <random>

#pragma pack(push, 1)
struct BMPHeader {
	char signature[2];
	uint32_t file_size;
	uint32_t reserved;
	uint32_t data_offset;
	uint32_t header_size;
	uint32_t width;
	uint32_t height;
	uint16_t color_planes;
	uint16_t bits_per_pixel;
	uint32_t compression;
	uint32_t image_size;
	uint32_t x_pixels_per_meter;
	uint32_t y_pixels_per_meter;
	uint32_t total_colors;
	uint32_t important_colors;
};
#pragma pack(pop)

struct Point {
	int x;
	int y;

	Point(int _x, int _y) : x(_x), y(_y) {}
};

std::mutex m;

std::vector<uint8_t> read_bmp(const std::string& file_path, uint32_t& width, uint32_t& height) {
	std::ifstream file(file_path, std::ios::binary);
	if (!file) {
		std::cerr << "Unable to open file: " << file_path << std::endl;
		return {};
	}

	BMPHeader header;
	file.read(reinterpret_cast<char*>(&header), sizeof(header));

	if (header.signature[0] != 'B' || header.signature[1] != 'M') {
		std::cerr << "Invalid BMP file: " << file_path << std::endl;
		return {};
	}

	if (header.bits_per_pixel != 24) {
		std::cerr << "Only 24-bit BMP files are supported." << std::endl;
		return {};
	}

	std::vector<uint8_t> image_data(header.image_size);
	file.read(reinterpret_cast<char*>(image_data.data()), header.image_size);

	if (!file) {
		std::cerr << "Error reading image data." << std::endl;
		return {};
	}

	width = header.width;
	height = header.height;

	return image_data;
}

void write_bmp(const std::string& file_path, const std::vector<uint8_t>& image_data, int width, int height) {
	BMPHeader header;
	std::fill_n(header.signature, 2, 0);
	header.signature[0] = 'B';
	header.signature[1] = 'M';
	header.file_size = sizeof(header) + image_data.size();
	header.reserved = 0;
	header.data_offset = sizeof(header);
	header.header_size = sizeof(header) - 14;
	header.width = width;
	header.height = height;
	header.color_planes = 1;
	header.bits_per_pixel = 24;
	header.compression = 0;
	header.image_size = image_data.size();
	header.x_pixels_per_meter = 0;
	header.y_pixels_per_meter = 0;
	header.total_colors = 0;
	header.important_colors = 0;

	std::ofstream file(file_path, std::ios::binary);
	if (!file) {
		std::cerr << "Unable to create file: " << file_path << std::endl;
		return;
	}

	file.write(reinterpret_cast<const char*>(&header), sizeof(header));
	file.write(reinterpret_cast<const char*>(image_data.data()), image_data.size());

	if (!file) {
		std::cerr << "Error writing image data." << std::endl;
	}
}

void drawLine(std::vector<uint8_t>& image_data, int width, int xA, int yA, int xB, int yB, int thickness) {
	int rest = (width * 3) % 4;
	int padding = 4 - rest;
	if (rest == 0) padding = 0;

	int dx = std::abs(xB - xA);
	int dy = std::abs(yB - yA);
	int sx = (xA < xB) ? 1 : -1;
	int sy = (yA < yB) ? 1 : -1;

	int err = dx - dy;
	int err2;

	int dirX = (xA < xB) ? 1 : -1;
	int dirY = (yA < yB) ? 1 : -1;
	int thicknessHalf = thickness / 2;

	while (true) {
		for (int i = -thicknessHalf; i <= thicknessHalf; i++) {
			for (int j = -thicknessHalf; j <= thicknessHalf; j++) {
				int xFinal = xA + i;
				int yFinal = yA + j;
				if (xFinal >= 0 && xFinal < width && yFinal >= 0 && yFinal < width) {
					int index = ((yFinal * width + xFinal) * 3) + (padding * yFinal);
					image_data[index] = 255;
					image_data[index + 1] = 255;
					image_data[index + 2] = 255;
				}
			}
		}

		if (xA == xB && yA == yB) {
			break;
		}

		err2 = 2 * err;

		if (err2 > -dy) {
			err -= dy;
			xA += sx;
		}

		if (err2 < dx) {
			err += dx;
			yA += sy;
		}
	}
}

void OMPDrawLine(std::vector<uint8_t>& image_data, int width, int xA, int yA, int xB, int yB, int thickness, int numThreads) {
	int rangeX = (xB - xA) / numThreads; // Podział odcinka na zakresy dla wątków
	int rangeY = (yB - yA) / numThreads;

	int X1, X2, Y1, Y2 = 0;

	omp_set_num_threads(numThreads);

#pragma omp parallel for
	for (int t = 0; t < numThreads; t++) {
		X1 = xA + t * rangeX;
		X2 = (t == numThreads - 1) ? xB : X1 + rangeX;
		Y1 = yA + t * rangeY;
		Y2 = (t == numThreads - 1) ? yB : Y1 + rangeY;
		drawLine(image_data, width, X1, Y1, X2, Y2, thickness);
	}
}

void ThreadDrawLine(std::vector<uint8_t>& image_data, int width, int xA, int yA, int xB, int yB, int thickness, int numThreads) {
	std::vector<std::thread> threads(numThreads);

	int rangeX = (xB - xA) / numThreads;
	int rangeY = (yB - yA) / numThreads;

	int X1, X2, Y1, Y2 = 0;

	for (int t = 0; t < numThreads; t++) {
		X1 = xA + t * rangeX;
		X2 = (t == numThreads - 1) ? xB : X1 + rangeX;
		Y1 = yA + t * rangeY;
		Y2 = (t == numThreads - 1) ? yB : Y1 + rangeY;

		threads[t] = std::thread(drawLine, std::ref(image_data), width, X1, Y1, X2, Y2, thickness);
	}

	for (auto& thread : threads) {
		thread.join();
	}
}

std::vector<Point> readPointsFromCSV(const std::string& filename) {
	std::vector<Point> points;

	std::ifstream file(filename);
	if (!file) {
		std::cout << "Error opening file: " << filename << std::endl;
		return points;
	}

	std::string line;
	while (std::getline(file, line)) {
		std::istringstream iss(line);
		std::string xStr, yStr;

		if (std::getline(iss, xStr, ',') && std::getline(iss, yStr, ',')) {
			int x = std::stoi(xStr);
			int y = std::stoi(yStr);
			Point point(x, y);
			points.push_back(point);
		}
		else {
			std::cout << "Invalid line format: " << line << std::endl;
		}
	}

	return points;
}

std::vector<Point> generateRandomPoints(int numPoints, int imageWidth) {
	std::vector<Point> points;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> distribution(0, imageWidth - 1);

	for (int i = 0; i < numPoints; ++i) {
		int x = distribution(gen);
		int y = distribution(gen);
		Point point(x, y);
		points.push_back(point);
	}

	return points;
}

int main() {
#pragma region Variables init and input reading
	int numbOfThreads;
	int lineThicknes;
	int numberOfPoints;
	double millisecondsSeq = 0.0;
	double millisecondsThread = 0.0;
	double millisecondsOMP = 0.0;
	double millisecondsSeq_Curve = 0.0;
	double millisecondsThread_Curve = 0.0;
	double millisecondsOMP_Curve = 0.0;
	double millisecondsSeq_Curve_Rand = 0.0;
	double millisecondsThread_Curve_Rand = 0.0;
	double millisecondsOMP_Curve_Rand = 0.0;

	std::string input_file = "xd.bmp";
	std::string output_file_Seq = "Line_Seq.bmp";
	std::string output_file_OMP = "Line_OMP.bmp";
	std::string output_file_Thread = "Line_Thread.bmp";
	std::string output_file_Seq_Curve = "Curve_Seq.bmp";
	std::string output_file_Thread_Curve = "Curve_Thread.bmp";
	std::string output_file_OMP_Curve = "Curve_OMP.bmp";
	std::string output_file_Seq_Curve_Rand = "Curve_Seq_Rand.bmp";
	std::string output_file_Thread_Curve_Rand = "Curve_Thread_Rand.bmp";
	std::string output_file_OMP_Curve_Rand = "Curve_OMP_Rand.bmp";

	std::vector<Point> points = readPointsFromCSV("points.csv");
	std::cout << "Input the number of threads to parallel processing \n";
	std::cin >> numbOfThreads;
	std::cout << std::endl;
	std::cout << "Input the line thickness \n";
	std::cin >> lineThicknes;
	std::cout << std::endl;
	std::cout << "Input the number of points to generate \n";
	std::cin >> numberOfPoints;
	std::cout << std::endl;
	// Read BMP file

	uint32_t width, height;
	std::vector<uint8_t> image_data_Seq = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_OMP = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_Thread = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_Seq_Curve = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_OMP_Curve = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_Thread_Curve = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_Seq_Curve_Rand = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_OMP_Curve_Rand = read_bmp(input_file, width, height);
	std::vector<uint8_t> image_data_Thread_Curve_Rand = read_bmp(input_file, width, height);

#pragma endregion

	if (!image_data_Seq.empty()) {
#pragma region Linear algorithm

#pragma region Sequential

		// Define points A and B and line thickness
		Point A(0, 0);
		Point B(5000, 5000);
		int thickness = 5;
		//***************************************Sequential***************************************
		// Start the timer
		auto beginSeq = std::chrono::high_resolution_clock::now();
		// Draw the line on the image data
		drawLine(image_data_Seq, width, A.x, A.y, B.x, B.y, thickness);
		auto endSeq = std::chrono::high_resolution_clock::now();
		auto elapsedSeq = std::chrono::duration_cast<std::chrono::nanoseconds>(endSeq - beginSeq);
		millisecondsSeq = elapsedSeq.count() * 1e-6;
		std::cout << "\nExecution time of SEQUENTIAL algorithm is: " << std::fixed << std::setprecision(6) << millisecondsSeq
			<< " [ms]\n";

#pragma endregion

#pragma region OpenMP

		//***************************************OpenMP***************************************
		// Start the timer
		auto beginOMP = std::chrono::high_resolution_clock::now();
		// Draw the line on the image data
		OMPDrawLine(image_data_OMP, width, A.x, A.y, B.x, B.y, thickness, numbOfThreads);
		auto endOMP = std::chrono::high_resolution_clock::now();
		auto elapsedOMP = std::chrono::duration_cast<std::chrono::nanoseconds>(endOMP - beginOMP);
		millisecondsOMP = elapsedOMP.count() * 1e-6;
		std::cout << "\nExecution time of OpenMP algorithm is: " << std::fixed << std::setprecision(6) << millisecondsOMP
			<< " [ms]\n";

#pragma endregion

#pragma region Thread

		//***************************************Thread***************************************
		// Start the timer
		auto beginThread = std::chrono::high_resolution_clock::now();
		// Draw the line on the image data
		ThreadDrawLine(image_data_Thread, width, A.x, A.y, B.x, B.y, thickness, numbOfThreads);
		auto endThread = std::chrono::high_resolution_clock::now();
		auto elapsedThread = std::chrono::duration_cast<std::chrono::nanoseconds>(endThread - beginThread);
		millisecondsThread = elapsedThread.count() * 1e-6;
		std::cout << "\nExecution time of Thread algorithm is: " << std::fixed << std::setprecision(6) << millisecondsThread
			<< " [ms]\n";

#pragma endregion

#pragma endregion
		//*************************************************CSV FILE*************************************************
#pragma region CSV file algorithm

#pragma region Sequential

		//***************************************Sequential***************************************
	  // Start the timer
		auto beginSeqCurve = std::chrono::high_resolution_clock::now();
		// Draw the line on the image data
		for (int i = 0; i < points.size() - 1; ++i) {
			const Point& point = points[i];
			const Point& nextPoint = points[i + 1];
			drawLine(image_data_Seq_Curve, width, point.x, point.y, nextPoint.x, nextPoint.y, thickness);
		}
		auto endSeqCurve = std::chrono::high_resolution_clock::now();
		auto elapsedSeqCurve = std::chrono::duration_cast<std::chrono::nanoseconds>(endSeqCurve - beginSeqCurve);
		millisecondsSeq_Curve = elapsedSeqCurve.count() * 1e-6;
		std::cout << "\nExecution time of SEQUENTIAL Curve algorithm is: " << std::fixed << std::setprecision(6) << millisecondsSeq_Curve
			<< " [ms]\n";

#pragma endregion

#pragma region Open MP

		//***************************************Open MP***************************************
		omp_set_num_threads(points.size()-1);
		auto beginOMPCurve = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
		for (int i = 0; i < points.size() - 1; ++i) {
			const Point& point = points[i];
			const Point& nextPoint = points[i + 1];
			drawLine(image_data_OMP_Curve, width, point.x, point.y, nextPoint.x, nextPoint.y, thickness);
		}

		auto endOMPCurve = std::chrono::high_resolution_clock::now();
		auto elapsedOMPCurve = std::chrono::duration_cast<std::chrono::nanoseconds>(endOMPCurve - beginOMPCurve);
		millisecondsOMP_Curve = elapsedOMPCurve.count() * 1e-6;

		std::cout << "\nExecution time of OprnMP Curve algorithm is: " << std::fixed << std::setprecision(6) << millisecondsOMP_Curve
			<< " [ms]\n";

#pragma endregion

#pragma region Thread

		auto beginThreadCurve = std::chrono::high_resolution_clock::now();

		// Define the number of threads to use
		int numThreads = points.size() - 1;
		std::vector<std::thread> threads(numThreads);

		for (int i = 0; i < numThreads; ++i) {
			const Point& point = points[i];
			const Point& nextPoint = points[i + 1];
			threads[i] = std::thread(drawLine, std::ref(image_data_Thread_Curve), width, point.x, point.y, nextPoint.x, nextPoint.y, thickness);
		}

		// Wait for all threads to finish
		for (auto& thread : threads) {
			thread.join();
		}

		auto endThreadCurve = std::chrono::high_resolution_clock::now();
		auto elapsedThreadCurve = std::chrono::duration_cast<std::chrono::nanoseconds>(endThreadCurve - beginThreadCurve);
		millisecondsThread_Curve = elapsedThreadCurve.count() * 1e-6;

		std::cout << "\nExecution time of Thread Curve algorithm is: " << std::fixed << std::setprecision(6) << millisecondsThread_Curve
			<< " [ms]\n";

#pragma endregion

#pragma endregion

		//************************************************* Random points *************************************************
#pragma region Random points algorithm

		std::vector<Point> randomPoints = generateRandomPoints(numberOfPoints, width);

#pragma region Sequential

		//***************************************Sequential***************************************
	  // Start the timer
		auto beginSeqCurve_Rand = std::chrono::high_resolution_clock::now();
		// Draw the line on the image data
		for (int i = 0; i < randomPoints.size() - 1; ++i) {
			const Point& point = randomPoints[i];
			const Point& nextPoint = randomPoints[i + 1];
			drawLine(image_data_Seq_Curve_Rand, width, point.x, point.y, nextPoint.x, nextPoint.y, thickness);
		}
		auto endSeqCurve_Rand = std::chrono::high_resolution_clock::now();
		auto elapsedSeqCurve_Rand = std::chrono::duration_cast<std::chrono::nanoseconds>(endSeqCurve_Rand - beginSeqCurve_Rand);
		millisecondsSeq_Curve_Rand = elapsedSeqCurve_Rand.count() * 1e-6;
		std::cout << "\nExecution time of SEQUENTIAL Curve_Rand algorithm is: " << std::fixed << std::setprecision(6) << millisecondsSeq_Curve_Rand
			<< " [ms]\n";

#pragma endregion

#pragma region Open MPb

		//***************************************Open MP***************************************
		omp_set_num_threads(randomPoints.size());
		auto beginOMPCurve_Rand = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
		for (int i = 0; i < randomPoints.size() - 1; ++i) {
			const Point& point = randomPoints[i];
			const Point& nextPoint = randomPoints[i + 1];
			drawLine(image_data_OMP_Curve_Rand, width, point.x, point.y, nextPoint.x, nextPoint.y, thickness);
		}

		auto endOMPCurve_Rand = std::chrono::high_resolution_clock::now();
		auto elapsedOMPCurve_Rand = std::chrono::duration_cast<std::chrono::nanoseconds>(endOMPCurve_Rand - beginOMPCurve_Rand);
		millisecondsOMP_Curve_Rand = elapsedOMPCurve_Rand.count() * 1e-6;

		std::cout << "\nExecution time of OprnMP Curve_Rand algorithm is: " << std::fixed << std::setprecision(6) << millisecondsOMP_Curve_Rand
			<< " [ms]\n";

#pragma endregion

#pragma region Thread

		auto beginThreadCurve_Rand = std::chrono::high_resolution_clock::now();

		// Define the number of threads to use
		numThreads = randomPoints.size() - 1;
		threads.clear();
		std::vector<std::thread> threadsRand(numThreads);

		for (int i = 0; i < numThreads; ++i) {
			const Point& point = randomPoints[i];
			const Point& nextPoint = randomPoints[i + 1];
			threadsRand[i] = std::thread(drawLine, std::ref(image_data_Thread_Curve_Rand), width, point.x, point.y, nextPoint.x, nextPoint.y, thickness);
		}

		// Wait for all threads to finish
		for (auto& thread : threadsRand) {
			thread.join();
		}

		auto endThreadCurve_Rand = std::chrono::high_resolution_clock::now();
		auto elapsedThreadCurve_Rand = std::chrono::duration_cast<std::chrono::nanoseconds>(endThreadCurve_Rand - beginThreadCurve_Rand);
		millisecondsThread_Curve_Rand = elapsedThreadCurve_Rand.count() * 1e-6;

		std::cout << "\nExecution time of Thread Curve_Rand algorithm is: " << std::fixed << std::setprecision(6) << millisecondsThread_Curve_Rand
			<< " [ms]\n";

#pragma endregion

#pragma endregion
		// Write BMP file
		write_bmp(output_file_Seq, image_data_Seq, width, height);
		write_bmp(output_file_OMP, image_data_OMP, width, height);
		write_bmp(output_file_Thread, image_data_Thread, width, height);
		write_bmp(output_file_Seq_Curve, image_data_Seq_Curve, width, height);
		write_bmp(output_file_OMP_Curve, image_data_OMP_Curve, width, height);
		write_bmp(output_file_Thread_Curve, image_data_Thread_Curve, width, height);
		write_bmp(output_file_Seq_Curve_Rand, image_data_Seq_Curve_Rand, width, height);
		write_bmp(output_file_OMP_Curve_Rand, image_data_OMP_Curve_Rand, width, height);
		write_bmp(output_file_Thread_Curve_Rand, image_data_Thread_Curve_Rand, width, height);
	}

	return 0;
}