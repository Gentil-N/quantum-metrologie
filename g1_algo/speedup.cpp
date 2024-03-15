#include <cassert>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <cmath>

void write_file(const char *filename, std::vector<float> &array)
{
    std::cout << "Start writing file...";
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        std::cerr << "Error opening the file " << filename << std::endl;
        exit(1);
    }
    for (size_t i = 0; i < array.size(); ++i)
    {
        fprintf(fp, "%f\n", array[i]);
    }
    fclose(fp);
    std::cout << "done" << std::endl;
}

void write_file_darray(const char *filename, std::vector<std::vector<float>> &array)
{
    std::ofstream outfile(filename, std::ofstream::binary);

    size_t approx_size = sizeof(char) * 10 * array[0].size();
    std::string num_list;
    num_list.reserve(approx_size);
    //std::cout << array[0].size() << " " << array.size() << std::endl;
    //exit(1);
    for (size_t i = 0; i < array.size(); ++i)
    {
        for (size_t j = 0; j < array[i].size(); ++j)
        {
            num_list += std::to_string(array[i][j]) + "\n";
        }
        std::cout << "Row written: " << i << std::endl;
        outfile.write(num_list.c_str(), num_list.size());
        num_list.clear();
    }
    outfile.close();
}

std::vector<float> read_file(const char *filename)
{
    std::cout << "Start reading file...";
    std::string line;
    std::ifstream file(filename);
    std::vector<float> array;
    if (file.is_open())
    {
        while (getline(file, line))
        {
            array.push_back(std::stof(line));
        }
        file.close();
    }
    else
    {
        std::cerr << "Error opening the file " << filename << std::endl;
    }
    std::cout << "done" << std::endl;
    return array;
}

std::vector<std::vector<float>> compute_g1_norm(std::vector<std::vector<float>> &f_list, std::vector<std::vector<float>> &g_list)
{
    std::cout << "Start computing g1 norm... preparing data" << std::endl;
    size_t sample_count = f_list[0].size();
    std::vector<std::vector<float>> g1_norm;
    g1_norm.resize(sample_count);
    size_t num_signal = f_list.size();
    assert(num_signal == g_list.size());

    std::cout << "Computing h" << std::endl;
    std::vector<float> h;
    h.reserve(sample_count);
    for (size_t t = 0; t < sample_count; ++t)
    {
        float denominator = 0.0;
        for (size_t i = 0; i < num_signal; ++i)
        {
            denominator += (f_list[i][t] * f_list[i][t] + g_list[i][t] * g_list[i][t]) / num_signal;
        }
        h[t] = denominator;
    }
    std::cout << "Entering main loop..." << std::endl;
    for (size_t t1 = 0; t1 < sample_count; ++t1)
    {
        std::vector<float> temporary_g1_t1;
        temporary_g1_t1.resize(sample_count);
        for (size_t t2 = 0; t2 < sample_count; ++t2)
        {
            float numerator_part_one = 0.0;
            float numerator_part_two = 0.0;
            for (size_t i = 0; i < num_signal; ++i)
            {
                numerator_part_one += (f_list[i][t1] * f_list[i][t2] + g_list[i][t1] * g_list[i][t2]) / num_signal;
                numerator_part_two += (g_list[i][t1] * f_list[i][t2] - f_list[i][t1] * g_list[i][t2]) / num_signal;
            }
            temporary_g1_t1[t2] = sqrtf((numerator_part_one * numerator_part_one + numerator_part_two * numerator_part_two) / (h[t1] * h[t2]));
        }
        g1_norm[t1] = temporary_g1_t1;
        std::cout << "t1 = " << t1 << "/" << sample_count << " done!" << std::endl;
    }
    return g1_norm;
}

int main()
{
    /*__m256 a = _mm256_set_ps(0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f);
    __m256 b = _mm256_set_ps(1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
    __m256 c = _mm256_add_ps(a, b);
    float *d = (float*)&c;
    for (size_t i = 0; i < 8; ++i)
    {
        printf("%.f\n", d[i]);
    }*/
    /*std::vector<float> array = {24739.642751, 0.7, 0.1};
    std::string filename("./test");
    write_file(filename.c_str(), array);
    std::vector<float> new_array = read_file(filename.c_str());*/

    const size_t NUM_FG_SIGNAL = 40;
    std::vector<std::vector<float>> f_list, g_list;
    for (size_t i = 0; i < NUM_FG_SIGNAL; ++i)
    {
        f_list.push_back(read_file((std::string("./f") + std::to_string(i)).c_str()));
    }
    for (size_t i = 0; i < NUM_FG_SIGNAL; ++i)
    {
        g_list.push_back(read_file((std::string("./g") + std::to_string(i)).c_str()));
    }
    std::vector<std::vector<float>> g1_norm = compute_g1_norm(f_list, g_list);
    std::string output_filename("./g1_norm_speedup");
    write_file_darray(output_filename.c_str(), g1_norm);
    return 0;
}
