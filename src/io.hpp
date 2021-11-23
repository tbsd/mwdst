#pragma once

#include <filesystem>

#include "mvdst.hpp"
#include "types.hpp"

namespace graph::io {

// Чтение из входного файла
std::unique_ptr<RawVertices> read(const std::filesystem::path& inPath);

// Вывод в файл
void write(const Mvdst& problem, const std::filesystem::path& outPath);

}  // namespace graph::io
