#include <algorithm>
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
#include <thread>
#include <utility>

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

/**
 * Set cutoff to zero in order to remove the cutoff behavior
 */
std::vector<float> read_file(const char *filename, size_t cutoff)
{
    std::cout << "Start reading file...";
    std::string line;
    std::ifstream file(filename);
    std::vector<float> array;
    if (file.is_open())
    {
        size_t i = 0;
        while (getline(file, line))
        {
            if (i >= cutoff && cutoff != 0)
                break;
            array.push_back(std::stof(line));
            ++i;
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

/**
 * t1_end is not reached !
*/
void _compute_g1_norm_t1_block(const std::vector<std::vector<float>> &f_list, const std::vector<std::vector<float>> &g_list, const std::vector<float> &h, size_t t1_start, size_t t1_end, std::vector<std::vector<float>> &g1_norm_block, size_t &count)
{
    size_t sample_count = f_list[0].size();
    size_t num_signal = f_list.size();
    for (size_t t1 = t1_start; t1 < t1_end; ++t1)
    {
        g1_norm_block[t1 - t1_start].resize(sample_count);
        for (size_t t2 = 0; t2 < sample_count; ++t2)
        {
            float numerator_part_one = 0.0;
            float numerator_part_two = 0.0;
            for (size_t i = 0; i < num_signal; ++i)
            {
                numerator_part_one += (f_list[i][t1] * f_list[i][t2] + g_list[i][t1] * g_list[i][t2]) / num_signal;
                numerator_part_two += (g_list[i][t1] * f_list[i][t2] - f_list[i][t1] * g_list[i][t2]) / num_signal;
            }
            g1_norm_block[t1 - t1_start][t2] = sqrtf((numerator_part_one * numerator_part_one + numerator_part_two * numerator_part_two) / (h[t1] * h[t2]));
        }
        ++count;
    }
}

std::vector<std::vector<float>> compute_g1_norm(const std::vector<std::vector<float>> &f_list, const std::vector<std::vector<float>> &g_list)
{
    std::cout << "Start computing g1 norm... preparing data" << std::endl;
    size_t sample_count = f_list[0].size();
    size_t num_signal = f_list.size();
    assert(num_signal == g_list.size());

    std::cout << "Computing h" << std::endl;
    std::vector<float> h;
    h.resize(sample_count);
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
    size_t num_block = std::thread::hardware_concurrency();
    if (num_block == 0)
    {
        std::cerr << "Unable to determine core count!" << std::endl;
        exit(1);
    }
    std::vector<std::pair<float, float>> t1_boundaries_list;
    t1_boundaries_list.resize(num_block);
    {
        size_t sample_count_block = (size_t)(sample_count / num_block);
        size_t sample_count_remaining = sample_count - sample_count_block * num_block;
        size_t t1_start = 0;
        for (size_t i = 0; i < t1_boundaries_list.size(); ++i)
        {
            size_t t1_shift = i == t1_boundaries_list.size() - 1 ? sample_count_block  + sample_count_remaining : sample_count_block;
            t1_boundaries_list[i] = std::pair<float, float>(t1_start, t1_start + t1_shift);
            t1_start += sample_count_block;
        }
    }

    std::vector<std::vector<std::vector<float>>> g1_norm_block_list;
    g1_norm_block_list.resize(num_block);
    std::vector<size_t> thread_count_list;
    thread_count_list.resize(num_block);
    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_block; ++i)
    {
        std::cout << "Launching thread " << i << " out of " << num_block << " from " << t1_boundaries_list[i].first << " to " << t1_boundaries_list[i].second << std::endl;
        g1_norm_block_list[i].resize(t1_boundaries_list[i].second - t1_boundaries_list[i].first);
        threads.push_back(std::move(std::thread(_compute_g1_norm_t1_block, std::ref(f_list), std::ref(g_list), std::ref(h), t1_boundaries_list[i].first, t1_boundaries_list[i].second, std::ref(g1_norm_block_list[i]), std::ref(thread_count_list[i]))));
        std::cerr << "haha !" << std::endl;
    }
    bool are_all_terminated = false;
    while(!are_all_terminated)
    {
        size_t terminated = 0;
        for (size_t i = 0; i < thread_count_list.size(); ++i)
        {
            size_t max_sample = t1_boundaries_list[i].second - t1_boundaries_list[i].first;
            std::cout << "t(" << i << "): " << (float)((size_t)(10000.0 * (float)thread_count_list[i] / (float)(max_sample))) / 100.0 << "% - ";
            if (thread_count_list[i] >= max_sample -  1)
            {
                ++terminated;
            }
        }
        std::cout << std::endl;
        are_all_terminated = (terminated == thread_count_list.size());
    }
    for (size_t i = 0; i < threads.size(); ++i)
    {
        std::cout << "Waiting for thread " << i << std::endl;
        threads[i].join();
    }
    std::cout << "Replacing data..." << std::endl;
    std::vector<std::vector<float>> g1_norm;
    g1_norm.reserve(sample_count);
    for (size_t i = 0; i < g1_norm_block_list.size(); ++i)
    {
        std::cout << "Inserting block " << i << std::endl;
        g1_norm.insert(g1_norm.end(), g1_norm_block_list[i].begin(),g1_norm_block_list[i].end());
    }
    /*for (size_t t1 = 0; t1 < sample_count; ++t1)
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
    }*/

    return g1_norm;
}

float gaussian(float x, float sigma)
{
    return 1.0 / (std::sqrt(2.0 * M_PI) * sigma) * std::exp(-x*x / (2.0 * sigma*sigma));
}

std::vector<std::vector<float>> generate_gauss_kernel(size_t mx, size_t my)
{
    float sigma_x = 2.0 * (float)mx / 6.0;
    float sigma_y = 2.0 * (float)my / 6.0;
    std::vector<std::vector<float>> kernel;
    kernel.resize(mx * 2);
    float total = 0.0;
    for (size_t i = 0; i < kernel.size(); ++i)
    {
        kernel[i].resize(my * 2);
        for (size_t j = 0; j < kernel[i].size(); ++j)
        {
            float fx = gaussian((float)(i - mx + 0.5), sigma_x);
            float fy = gaussian((float)(j - my + 0.5), sigma_y);
            float res = fx * fy;
            total += res;
            kernel[i][j] = res;
        }
    }
    for (size_t i = 0; i < kernel.size(); ++i) // rescale to have one as the sum of the matrix
    {
        for (size_t j = 0; j < kernel.size(); ++j)
        {
            kernel[i][j] /= total;
        }
    }
    return kernel;
}

std::vector<std::vector<float>> downscale(const std::vector<std::vector<float>> &g1_norm, size_t div_width, size_t div_height)
{
    assert(g1_norm.size() % div_width == 0 && g1_norm[0].size() % div_height == 0);
    float mx = (float)div_width;
    float my = (float)div_height;
    std::vector<std::vector<float>> kernel = generate_gauss_kernel(mx, my);
    float width = (float)g1_norm.size() / (float)div_width;
    float height = (float)g1_norm[0].size() / (float)div_height;
    std::vector<std::vector<float>> downscaled;
    downscaled.resize(width);
    for (size_t i = 0; i < width ; ++i)
    {
        downscaled[i].resize(height);
        for (size_t j = 0; j < height; ++j)
        {
            //std::cout << "processing " << i << ", " << j << std::endl;
            float res = 0.0;
            for (size_t x = 0; x < (mx * 2); ++x)
            {
                for (size_t y = 0; y < (my * 2); ++y)
                {
                    float g1_norm_index_i = (i * mx) + (float)x - mx + 1;
                    float g1_norm_index_j = (j * my) + (float)y - my + 1;
                    if (g1_norm_index_i < 0.0 || g1_norm_index_i >= g1_norm.size() || g1_norm_index_j < 0.0 || g1_norm_index_j >= g1_norm[0].size())
                        continue;
                    //std::cout << g1_norm_index_i << " " << g1_norm_index_j << std::endl;
                    res += g1_norm[(size_t)g1_norm_index_i][(size_t)g1_norm_index_j] * kernel[x][y];
                }
            }
            downscaled[i][j] = res;
        }
    }
    return downscaled;
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

    const size_t START_SIGNAL = 1, STOP_SIGNAL = 10;
    const size_t CUTOFF = 10000;
    std::vector<std::vector<float>> f_list, g_list;
    for (size_t i = START_SIGNAL; i < STOP_SIGNAL + 1; ++i)
    {
        f_list.push_back(read_file((std::string("./f") + std::to_string(i)).c_str(), CUTOFF));
    }
    for (size_t i = START_SIGNAL; i < STOP_SIGNAL + 1; ++i)
    {
        g_list.push_back(read_file((std::string("./g") + std::to_string(i)).c_str(), CUTOFF));
    }
    std::vector<std::vector<float>> g1_norm = compute_g1_norm(f_list, g_list);
    std::cout << "Downscaling..." << std::endl;
    std::vector<std::vector<float>> g1_norm_downscaled = downscale(g1_norm, 10, 10);
    std::string output_filename("./g1_norm_speedup");
    write_file_darray(output_filename.c_str(), g1_norm_downscaled);
    return 0;
}
